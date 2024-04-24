/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

// MOOSE includes
#include "DualProjectionAMR.h"
#include "MooseTypes.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "MultiApp.h"
#include "MooseMesh.h"

#include "MooseVariable.h"

// libMesh includes
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/system.h"
#include "libmesh/radial_basis_interpolation.h"

#include "libmesh/mesh_tools.h"

template<>
InputParameters validParams<DualProjectionAMR>()
{
    InputParameters params = validParams<MultiAppTransfer>();
    
    params.addRequiredParam<std::vector<VariableName> >("variable", "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<std::vector<VariableName> >("source_variable", "The variable to transfer from.");
    
    params.addRequiredParam<bool>("dual", "dual");
    params.addRequiredParam<bool>("add", "add");
    params.addRequiredParam<bool>("has_hanging_nodes", "has_hanging_nodes");
    
    return params;
}

DualProjectionAMR::DualProjectionAMR(const InputParameters & parameters) :
MultiAppTransfer(parameters),
_to_var_name(getParam<std::vector<VariableName> >("variable")),
_from_var_name(getParam<std::vector<VariableName> >("source_variable")),
_add(getParam<bool>("add")),
_dual(getParam<bool>("dual")),
_has_hanging_nodes(getParam<bool>("has_hanging_nodes"))
{
    _console<<"Constructor of DualProjectionAMR "<<name() << std::endl;
    if ( _from_var_name.size() != _to_var_name.size() )
    {
        _console<<"The number of from and to varaibles is different\n";
        exit(1);
    }
    
    for (int i=0; i<_from_var_name.size(); ++i)
        _console<<name()<<" is transferring "<<_from_var_name.at(i)<<" into "<<_to_var_name.at(i)<<std::endl;

    _console<<std::endl<<std::endl;
}

void
DualProjectionAMR::initialSetup()
{
    _console<<"initialSetup of DualProjectionAMR "<<name() << std::endl;
    
    if (_direction == TO_MULTIAPP)
    {
        _from_problem = & _multi_app->problemBase();
        _to_problem   = & _multi_app->appProblemBase(0);
    }
    else
    {
        _from_problem = &_multi_app->appProblemBase(0); // There is only one sub-app
        _to_problem   = & _multi_app->problemBase();
    }
    
    //_from_problem[0].setConstJacobian(true);
    //_to_problem[0].setConstJacobian(true);
    
    _from_mesh   = &_from_problem[0].mesh().getMesh();
    _to_mesh = &_to_problem[0].mesh().getMesh(); // here 0 in the index of the multiapp
    
    _local_from_nodes   =_from_mesh[0].n_local_nodes();
    _local_to_nodes     =  _to_mesh[0].n_local_nodes();
    
    _total_from_nodes   = _from_mesh[0].n_nodes();
    _total_to_nodes     = _to_mesh[0].n_nodes();
    
    _my_communicator = new Parallel::Communicator(   _communicator.get());
    _from_communicator = &_from_mesh[0].comm();
    _to_communicator = &_to_mesh[0].comm();
    
    _proj_matrix = new PetscMatrix<Number>(_my_communicator[0]);
    
    _from_vector = new PetscVector<Number>(_my_communicator[0],_total_from_nodes,_local_from_nodes);
    _to_vector   = new PetscVector<Number>(_my_communicator[0],_total_to_nodes  ,_local_to_nodes  );
    
    if (_dual)
    {
        _fine_mesh  =_from_mesh;
        _coarse_mesh=_to_mesh;
    }
    else
    {
        // if it is not dual, the fine mesh is the to
        _fine_mesh=_to_mesh;
        // if it is not dual, the coarse mesh is the from
        _coarse_mesh=_from_mesh;
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    constructInterpolationMatrix();
    auto t2 = std::chrono::high_resolution_clock::now();
    _console<<"constructProjectionMatrix() took "<< (std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count())/1000.0<< " seconds\n";
    
    if (_dual && _has_hanging_nodes)
    {
        _reorganization_matrix = new PetscMatrix<Number>(_my_communicator[0]);
        constructReorganizationMatrix();
    }
    
    if (0)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        _proj_matrix[0].print_matlab("projectionMatrix.m");
        auto t2 = std::chrono::high_resolution_clock::now();
        _console<<"print_matlab() took "<< (std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count())/1000.0<< " seconds\n";
    }
    
    //  if (_direction == TO_MULTIAPP)
    //    variableIntegrityCheck(_to_var_name);
    //  else
    //    variableIntegrityCheck(_from_var_name);

}

void DualProjectionAMR::execute()
{
    _console << "Beginning DualProjection " << name() << std::endl;
    
    // Construct references to projection matrix and meshes to have an easier access to names
    PetscMatrix<Number> & proj_matrix=_proj_matrix[0];
    
    // start
    for (int variable_iterator=0; variable_iterator<_from_var_name.size(); ++variable_iterator)
    {
    
    
        MooseVariable & from_var = _from_problem[0].getVariable(0, _from_var_name.at(variable_iterator));         ///////////////////// ATTENTO QUI
        SystemBase & from_system_base = from_var.sys();
        System & from_sys = from_system_base.system();
        unsigned int from_var_num = from_sys.variable_number(from_var.name());
        //const DofMap& from_dof_map = from_sys.get_dof_map();
        unsigned int from_sys_num = from_sys.number();
        NumericVector<Number> & from_solution =  *(from_sys.solution.get()); // *(fine_sys.current_local_solution.get()); //
        
        
        MooseVariable & to_var = _to_problem[0].getVariable(0, _to_var_name.at(variable_iterator));
        SystemBase & to_system_base = to_var.sys();
        System & to_sys = to_system_base.system();
        unsigned int to_var_num = to_sys.variable_number(to_var.name());
        //const DofMap& to_dof_map = to_sys.get_dof_map();
        unsigned int to_sys_num = to_sys.number();
        NumericVector<Number> & to_solution = *(to_sys.solution.get());
    
    
    
    MeshBase::const_node_iterator from_it = _from_mesh[0].local_nodes_begin();
    const MeshBase::const_node_iterator from_it_end = _from_mesh[0].local_nodes_end();
    for ( ; from_it != from_it_end; ++from_it)
    {
        const Node * node = *from_it;
        std::vector<dof_id_type> from_dof;
        std::vector<Number> value;
        from_dof.push_back( node[0].dof_number(from_sys_num, from_var_num, 0) );
        from_solution.get(from_dof,value);
        
        _from_vector[0].set(node->id(),value.at(0));
        
    }
    
    _from_vector[0].close();
    proj_matrix.vector_mult(_to_vector[0],_from_vector[0]);
    _to_vector[0].close();
    
    if (_dual && _has_hanging_nodes)
    {
        PetscVector<Number> _to_vector2(_my_communicator[0],_total_to_nodes  ,_local_to_nodes  );
        _reorganization_matrix[0].vector_mult(_to_vector2,_to_vector[0]);
        _to_vector2.close();
        
        _to_vector[0]=_to_vector2;
        _to_vector[0].close();
    }
    
    for (int iiii=_to_vector[0].first_local_index(); iiii<_to_vector[0].last_local_index(); ++iiii)
    {
        Node const & nodo=_to_mesh[0].node_ref(iiii);
        Number value=_to_vector[0](iiii);
        
        dof_id_type da_push=nodo.dof_number(to_sys_num, to_var_num , 0);
        if (_add)
        {
            to_solution.add(da_push,value);
        }
        else
        {
            to_solution.set(da_push,value);
        }
    }

    to_solution.close();
    }
    
    _console << "Finished InterpolationTransfer " << name() << std::endl;
}

bool DualProjectionAMR::is_inside()
{
    Real xc=_centroid(0);
    Real yc=_centroid(1);
    Real zc=_centroid(2);
    
    if (_xmin<xc && xc<_xmax)
        if (_ymin<yc && yc<_ymax)
            return true;
    
    return false;
    
}

void DualProjectionAMR::constructInterpolationMatrix()
{
    
    std::cout<<"\nStarted construction of interpolation matrix\n\n";
    // Construct references to projection matrix and meshes to have an easier access to names
    
    MeshBase const & coarse_mesh=_coarse_mesh[0];
    MeshBase const &   fine_mesh=  _fine_mesh[0];
    
    int _local_fine_nodes       =    _fine_mesh[0].n_local_nodes();
    int _local_coarse_nodes     =  _coarse_mesh[0].n_local_nodes();
    
    int _total_fine_nodes   =   _fine_mesh[0].n_nodes();
    int _total_coarse_nodes = _coarse_mesh[0].n_nodes();
    
    //std::cout<<"_total_fine_nodes "<<_total_fine_nodes<<std::endl;
    //std::cout<<"_total_coarse_nodes "<<_total_coarse_nodes<<std::endl;

    
    PetscMatrix<Number> & proj_matrix=_proj_matrix[0];
    
    proj_matrix.init(_total_fine_nodes,_total_coarse_nodes,
                     _local_fine_nodes,_local_coarse_nodes);
    
    PetscMatrix<Number> * remove_hanging_nodes;
    
    if (_has_hanging_nodes)
    {
        remove_hanging_nodes=new PetscMatrix<Number>( _my_communicator[0] );
        remove_hanging_nodes[0].init(_total_fine_nodes,_total_coarse_nodes,
                                     _local_fine_nodes,_local_coarse_nodes);
    }
    
    std::map<int,Point> coarse_point_map;
    std::map<int,Point> fine_point_map;
    
    int const levels_fine=libMesh::MeshTools::n_levels(fine_mesh);
    int const levels_coarse=libMesh::MeshTools::n_levels(coarse_mesh);
    
    std::cout<<"Coarse mesh : levels="<<levels_coarse<<" elements="<<coarse_mesh.n_active_elem()<<" nodes="<<_total_coarse_nodes<<std::endl;
    std::cout<<"  Fine mesh : levels="<<  levels_fine<<" elements="<<  fine_mesh.n_active_elem()<<" nodes="<<  _total_fine_nodes<<std::endl<<std::endl;
    
    MeshBase::const_element_iterator coarse_el_it=coarse_mesh.active_local_elements_begin();
    MeshBase::const_element_iterator const coarse_el_it_end=coarse_mesh.active_local_elements_end();
    
    std::map<int,int> coarse2fine;
    int elem_counter=0;
    
    std::cout<<"\n Construction of the coarse2fine map\n\n";
    
    for ( ; coarse_el_it != coarse_el_it_end; ++coarse_el_it)
    {
        ++elem_counter;
        if (elem_counter%100==0)
        {
            std::cout<<"\r"<<1.0*elem_counter/coarse_mesh.n_active_local_elem()*100.0<<"%"<< std::flush;
        }
        Elem *coarse_el=*coarse_el_it;
        RealVectorValue coarseCentroid=coarse_el[0].centroid();
        int level=coarse_el[0].level();
        
        bool _has_been_insered=0;
        MeshBase::const_element_iterator fine_el_it=fine_mesh.level_elements_begin(level);
        MeshBase::const_element_iterator const fine_el_it_end=fine_mesh.level_elements_end(level);
        
        for ( ; fine_el_it != fine_el_it_end; ++fine_el_it)
        {
            Elem *fine_el=*fine_el_it;
            RealVectorValue fineCentroid=fine_el[0].centroid();
            
            RealVectorValue diff=fineCentroid-coarseCentroid;
            Real distance=std::sqrt(diff*diff);
            if (distance<1e-12)
            {
                std::map<int,int>::iterator it;
                it = coarse2fine.find( coarse_el[0].id() );
                if (it != coarse2fine.end())
                {
                    std::cout<<"you are trying to insert an element that is already exiting\n ... exiting \n";
                    exit(1);
                }
                coarse2fine.insert ( std::pair<int,int>(coarse_el[0].id(),fine_el[0].id() ) );
                _has_been_insered=1;
                break;
            }
            
        }
        // here we just verify that everything went well
        if (_has_been_insered==0)
        {
            std::cout<<"an element on the coarse has not been found in the fine\n ... exiting \n";
            exit(1);
        }
    }

    if (coarse_mesh.n_active_local_elem() != coarse2fine.size())
    {
        std::cout<<"the size of the map is different from the number of elements \n";
        exit(1);
    }
    
    std::cout<<"\n\n done!\n\n";
    
    
//    std::cout<<coarse2fine.size()<<std::endl;
//    
//    for (std::map<int,int>::iterator my_it=coarse2fine.begin(); my_it!=coarse2fine.end(); ++my_it)
//    {
//        std::cout<<my_it->first<<" "<<my_it->second<<std::endl;
//        
//    }
    
    std::cout<<"\n Assembly of projection matrix\n\n";
    elem_counter=0;
    for (std::map<int,int>::iterator my_it=coarse2fine.begin(); my_it!=coarse2fine.end(); ++my_it)
    {
        ++elem_counter;
        if (elem_counter%100==0)
        {
            std::cout<<"\r"<<1.0*elem_counter/coarse2fine.size()*100.0<<"%"<< std::flush;
        }
        
        Elem const & coarse_el=coarse_mesh.elem_ref(my_it->first);
        Elem const &   fine_el=  fine_mesh.elem_ref(my_it->second);
        
        
        // JUST A CHECK
        if (coarse_el.active() != true)
        {
            std::cout<<"coarse elements have to be active \n";
            exit(1);
        }
        
        // We fill the coarse map
        for (int coarseNodeLocalId=0; coarseNodeLocalId<coarse_el.n_nodes(); ++coarseNodeLocalId)
        {
            Node const & coarseNode=coarse_el.node_ref(coarseNodeLocalId);
            Point coarsePoint=static_cast<const Point>(coarseNode);
            coarse_point_map.insert ( std::pair<int,Point>(coarseNode.id(),coarsePoint) );
        }
        
        // We fill the fine map, first we check if it has children
        if ( fine_el.has_children() )
        {
            if (fine_el.active()==true)
            {
                std::cout<<"fine elements have not to be active \n";
                exit(1);
            }
            
            for (int i=0; i<fine_el.n_children(); ++i)
            {
                Elem const * childOfFine=fine_el.child_ptr(i);
                for (int fine_node=0; fine_node<childOfFine[0].n_nodes(); ++fine_node)
                {
                    Node const & FineNode=childOfFine[0].node_ref(fine_node);
                    Point finePoint =static_cast<const Point>(FineNode);
                    fine_point_map.insert ( std::pair<int,Point>(FineNode.id(),finePoint) );
                }
                
            }

        }
        else
        {
            if (fine_el.active()!=true)
            {
                std::cout<<"fine elements have not be active \n";
                exit(1);
            }
            for (int fineNodeLocalId=0; fineNodeLocalId<fine_el.n_nodes(); ++fineNodeLocalId)
            {
                Node const & fineNode=fine_el.node_ref(fineNodeLocalId);
                Point finePoint=static_cast<const Point>(fineNode);
                fine_point_map.insert ( std::pair<int,Point>(fineNode.id(),finePoint) );
            }


        }
        // We have filled the fine map
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////
        ////////////////////////////////          HERE WE START THE LOCAL TO GLOBAL
        ////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        for (std::map<int,Point>::iterator fine_it=fine_point_map.begin(); fine_it!=fine_point_map.end(); ++fine_it)
        {
            
            for (std::map<int,Point>::iterator coarse_it=coarse_point_map.begin(); coarse_it!=coarse_point_map.end(); ++coarse_it)
            {
                RealVectorValue A=fine_it->second-coarse_it->second;
                Real norm=std::sqrt(A*A);
                
                if (norm < 1e-12)
                {
                    proj_matrix.set(fine_it->first,coarse_it->first,1.0);
                }
            }


            _centroid=coarse_el.centroid();
            RealVectorValue A=_centroid-fine_it->second;
            Real norm=std::sqrt(A*A);
            if (norm < 1e-12)
            {
                for (std::map<int,Point>::iterator local_it=coarse_point_map.begin(); local_it!=coarse_point_map.end(); ++local_it)
                {
                    proj_matrix.set(fine_it->first,local_it->first,1.0/4.0);
                }
                //break;
            }

            for (int i=0; i<4; ++i)    ///////////////////////////////  TOGLI IL 4
            {
                std::unique_ptr<Elem const> _side=coarse_el.side_ptr(i);
                Elem const * side = _side.get();
                RealVectorValue A=side[0].centroid()-fine_it->second;
                Real norm=std::sqrt(A*A);
                if (norm < 1e-12)
                {
                    for (int l=0; l<side[0].n_nodes(); ++l)
                    {
                        Node const & l_nodo=side[0].node_ref(l);
                        proj_matrix.set(fine_it->first,l_nodo.id(),1.0/2.0);
                        if (_has_hanging_nodes)
                        {
                            remove_hanging_nodes[0].set(fine_it->first,l_nodo.id(),-1.0/2.0);
                        }
                    }
                }
            }
        }// end for fine map

        
        // we just clear the maps
        
        fine_point_map.clear();
        coarse_point_map.clear();
        
    }

    std::cout<<"\n done!\n\n Closing matrix"<<std::flush;
    
    proj_matrix.close();
    
    std::cout<<"\n done!\n\n";
    
    if (_has_hanging_nodes)
    {
        remove_hanging_nodes[0].close();
    }
    
    std::cout<<"\n done!\n\n";
    
    if (_has_hanging_nodes)
    {
        for (int i=_from_vector[0].first_local_index(); i<_from_vector[0].last_local_index(); ++i)
        {
            _from_vector[0].set(i,1.0);
        }
        proj_matrix.vector_mult(_to_vector[0],_from_vector[0]);
        //std::cout<<_to_vector[0]<<std::endl;
       
        if ( _to_vector[0].first_local_index() != remove_hanging_nodes[0].row_start() )
        {
            std::cout<<"error\n";
            exit(1);
        }

        if ( _to_vector[0].last_local_index() != remove_hanging_nodes[0].row_stop() )
        {
            std::cout<<"error\n";
            exit(1);
        }
        
        std::vector<numeric_index_type> rows;
        for (int i=_to_vector[0].first_local_index(); i<_to_vector[0].last_local_index(); ++i)
        {
            if (_to_vector[0](i).real()<1.5)
            {
                rows.push_back(i);
            }
        }

        remove_hanging_nodes[0].zero_rows (rows, 0.0);
        remove_hanging_nodes[0].close();

        proj_matrix.add (1.0, remove_hanging_nodes[0]);
        
        proj_matrix.close();
        
        delete remove_hanging_nodes;
    }
    
//    MeshBase::const_element_iterator fine_it=fine_mesh.active_local_elements_begin();
//    for (; fine_it != fine_mesh.active_local_elements_end(); ++fine_it)
//    {
//        Elem *fine_el=*fine_it;
//        
//        for (int fine_node=0; fine_node<fine_el[0].n_nodes(); ++fine_node)
//        {
//            Node const & fine_nodo=fine_el[0].node_ref(fine_node);
//            Point fine_punto =static_cast<const Point>(fine_nodo);
//            
//            fine_point_map.insert ( std::pair<int,Point>(fine_nodo.id(),fine_punto) );
//            
//        }
//        
//        // we have to coarsen
//        if ( fine_el[0].level()==levels_fine-1 )
//        {
//            // we take the pointer to the coarse element
//            Elem *coarse_el=fine_el[0].parent();
//            
//            // We fill the map of the coarse element
//            for (int coarse_node=0; coarse_node<coarse_el[0].n_nodes(); ++coarse_node)
//            {
//                Node const & coarse_nodo=coarse_el[0].node_ref(coarse_node);
//                Point coarse_punto =static_cast<const Point>(coarse_nodo);
//                
//                coarse_point_map.insert ( std::pair<int,Point>(coarse_nodo.id(),coarse_punto) );
//                
//            }
//            
//            // in this for loop, we try to put the write entries in the projection matrix
//            
//            for (std::map<int,Point>::iterator fine_it=fine_point_map.begin(); fine_it!=fine_point_map.end(); ++fine_it)
//            {
//                
//                for (std::map<int,Point>::iterator coarse_it=coarse_point_map.begin(); coarse_it!=coarse_point_map.end(); ++coarse_it)
//                {
//                    RealVectorValue A=fine_it->second-coarse_it->second;
//                    Real norm=std::sqrt(A*A);
//                    
//                    if (norm < 1e-12)
//                    {
//                        if (coarse_it->first<_total_coarse_nodes)
//                            proj_matrix.set(fine_it->first,coarse_it->first,1.0);
//                        //break;
//                    }
//                }
//                
//                // centroide
//                
//                _centroid=coarse_el[0].centroid();
//                RealVectorValue A=_centroid-fine_it->second;
//                Real norm=std::sqrt(A*A);
//                if (norm < 1e-12)
//                {
//                    for (std::map<int,Point>::iterator local_it=coarse_point_map.begin(); local_it!=coarse_point_map.end(); ++local_it)
//                    {
//                        if (local_it->first<_total_coarse_nodes)
//                            proj_matrix.set(fine_it->first,local_it->first,1.0/4.0);
//                    }
//                    //break;
//                }
//                
//                for (int i=0; i<4; ++i)    ///////////////////////////////  TOGLI IL 4
//                {
//                    std::unique_ptr<Elem const> _side=coarse_el[0].side_ptr(i);
//                    Elem const * side = _side.get();
//                    RealVectorValue A=side[0].centroid()-fine_it->second;
//                    Real norm=std::sqrt(A*A);
//                    if (norm < 1e-12)
//                    {
//                        for (int l=0; l<side[0].n_nodes(); ++l)
//                        {
//                            Node const & l_nodo=side[0].node_ref(l);
//                            if (l_nodo.id()<_total_coarse_nodes)
//                                proj_matrix.set(fine_it->first,l_nodo.id(),1.0/2.0);
//                        }
//                    }
//                }
//            }// endfor over fine
//            
//            
//            
//        }
//        else
//        {
//            // we leave as it is, easy
//            for (std::map<int,Point>::iterator fine_it=fine_point_map.begin(); fine_it!=fine_point_map.end(); ++fine_it)
//            {
//                proj_matrix.set(fine_it->first,fine_it->first,1.0);
//            }
//            
//        }
//        
//        fine_point_map.clear();
//        coarse_point_map.clear();
//        
//    }
//
    if (_dual)
    {
        PetscMatrix<Number> _proj_matrix2(_my_communicator[0]);
        proj_matrix.get_transpose(_proj_matrix2);
        proj_matrix.swap(_proj_matrix2);
        proj_matrix.close();
    }
    
}

void DualProjectionAMR::constructReorganizationMatrix()
{
    MeshBase const & coarse_mesh=_coarse_mesh[0];
    
    int _local_coarse_nodes     =  _coarse_mesh[0].n_local_nodes();
    
    int _total_coarse_nodes = _coarse_mesh[0].n_nodes();

    PetscMatrix<Number> & reorganization_matrix=_reorganization_matrix[0];
    
    reorganization_matrix.init(_total_coarse_nodes,_total_coarse_nodes,
                               _local_coarse_nodes,_local_coarse_nodes);
    
    std::map<dof_id_type, std::vector<dof_id_type> > hanging_nodes;
    libMesh::MeshTools::find_hanging_nodes_and_parents(coarse_mesh,hanging_nodes);
    
    std::cout<<"Number of hanging nodes on the coarse level "<<hanging_nodes.size()<<std::endl;

    MeshBase::const_node_iterator it = coarse_mesh.local_nodes_begin();
    //MeshBase::const_node_iterator it = coarse_mesh.active_nodes_begin();
    const MeshBase::const_node_iterator it_end = coarse_mesh.local_nodes_end();
    //const MeshBase::const_node_iterator it_end = coarse_mesh.active_nodes_end();
    for ( ; it != it_end; ++it)
    {
        const Node * node = *it;
        dof_id_type const id=node[0].id();
        
        if ( hanging_nodes.find(id)==hanging_nodes.end() )
        {
            reorganization_matrix.set(id,id,1.0);
        }
        else
        {
            std::map<dof_id_type, std::vector<dof_id_type> >::iterator hang_it=hanging_nodes.find(id);
//            std::cout<<hang_it->second.size()<<std::endl;
//            std::cout<<id<<" "<<hang_it->first<<std::endl;
            
            for (int i=0;i<hang_it->second.size(); ++i)
            {
                dof_id_type hang_id=hang_it->second.at(i);
                reorganization_matrix.set(hang_id,id,0.5);
            }
        }
        
    }

    reorganization_matrix.close();
//    reorganization_matrix.print_matlab("reorganizationMatrix.m");


}


DualProjectionAMR::~DualProjectionAMR()
{
    if (_dual && _has_hanging_nodes)
    {
        delete _reorganization_matrix;
    }
    
    delete _my_communicator;
    delete _proj_matrix;
    delete _from_vector;
    delete _to_vector;
}
