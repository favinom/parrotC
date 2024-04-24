//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleFreqPoroelasticityInclusion.h"

//registerMooseObject("parrotApp", AssembleFreqPoroelasticityInclusion);

template <>
InputParameters
	validParams<AssembleFreqPoroelasticityInclusion>()
{
	InputParameters params = validParams<AssemblyUO>();
	params.addRequiredParam< std::string >("material_name","the user object where assembly is performed.");
	params.addParam< bool >("center_x","the user object where assembly is performed.");
	params.addParam< bool >("center_y","the user object where assembly is performed.");
	params.addParam< bool >("center_z","the user object where assembly is performed.");
	return params;
}

AssembleFreqPoroelasticityInclusion::AssembleFreqPoroelasticityInclusion(const InputParameters & parameters) :
AssemblyUO(parameters),
_materialName( getParam< std::string >("material_name") ),
_linearSystemAllocated(false),
_linearSystemAssembled(false),
_center_x( parameters.isParamValid("center_x") ? getParam< bool >("center_x") : false ),
_center_y( parameters.isParamValid("center_y") ? getParam< bool >("center_y") : false ),
_center_z( parameters.isParamValid("center_z") ? getParam< bool >("center_z") : false ),
_pi( 4.0*std::atan(1.0) ),
_imagUnit(0.0,1.0)
{
	std::cout<<"AssembleFreqPoroelasticityInclusion::AssembleFreqPoroelasticityInclusion() start\n";
		
	if ( _dim==2 && _nvars!=3 )
	{
		mooseError("The _dimension is 2 but the number of variables is not 3!");
	}
	if ( _dim==3 && _nvars!=4 )
	{
		mooseError("The _dimension is 3 but the number of variables is not 4!");
	}
	
	_d_var.resize(_dim);
	
	if ( !_tnis.has_variable("disp_x") )
	{
		mooseError("The system has no variable disp_x!");
	}
	else
	{
		_d_var.at(0)= _tnis.variable_number ("disp_x");
	}
	
	if ( !_tnis.has_variable("disp_y") )
	{
		mooseError("The system has no variable disp_y!");
	}
	else
	{
		_d_var.at(1)= _tnis.variable_number ("disp_y");
	}
	
	if ( !_tnis.has_variable("pressure") )
	{
		mooseError("The system has no variable pressure!");
	}
	else
	{
		_p_var    = _tnis.variable_number ("pressure");
	}
	
	if (_dim==3)
	{
		if ( !_tnis.has_variable("disp_z") )
		{
			mooseError("The system has no variable disp_z!");
		}
		else
		{
			_d_var.at(2)= _tnis.variable_number ("disp_z");
		}
	}

	_identity=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for (int d=0; d<_dim; ++d)
	{
		_identity(d,d)=1.0;
	}
	std::cout<<"AssembleFreqPoroelasticityInclusion::AssembleFreqPoroelasticityInclusion() stop\n";
}

void AssembleFreqPoroelasticityInclusion::execute()
{
	std::cout<<"AssembleFreqPoroelasticityInclusion::execute() start\n";
	if (!_linearSystemAllocated)
	{
		mallocLinearSystem();
		_linearSystemAllocated=true;
	}
	if (!_linearSystemAssembled)
	{
		assembleLinearSystem();
		_linearSystemAssembled=true;
	}

	std::cout<<"AssembleFreqPoroelasticityInclusion::execute() stop\n";
}

AssembleFreqPoroelasticityInclusion::~AssembleFreqPoroelasticityInclusion()
{
	std::cout<<"AssembleFreqPoroelasticityInclusion::~AssembleFreqPoroelasticityInclusion() start\n";
	if (_linearSystemAllocated)
	{
		deleteLinearSystem();
		_linearSystemAllocated=false;
	}
	std::cout<<"AssembleFreqPoroelasticityInclusion::~AssembleFreqPoroelasticityInclusion() stop\n";
}

void AssembleFreqPoroelasticityInclusion::solve( Real t_in )
{
	
	std::cout<<"AssembleFreqPoroelasticityInclusion::solve() start\n";
 	auto start = high_resolution_clock::now(); 

 	_mat_SM_ptr->zero();
 	_mat_SM_ptr->add(1.0,_static_stiffness[0]);
	Real scale=1.0/(2.0*_pi*std::pow(10.0,t_in));
 	_mat_SM_ptr->add(scale,_freqdep_stiffness[0]);

 	auto stop = high_resolution_clock::now(); 
 	auto duration = duration_cast<microseconds>(stop - start);

 	myout << "Summing matrices took "<<duration.count()/1000.0/1000.0 <<" s"<< std::endl;

 	_parrotSolver.setMatrixAndVectors( _mat_SM_ptr , _rhs_NV_ptr , _curr_NV_ptr );
 	_parrotSolver.solve();

 	//std::cout<<_curr_NV_ptr->type()<<std::endl;
 	//std::cout<<_dirichletFlag->type()<<std::endl;
	_dirichletFlag[0]=_curr_NV_ptr[0];
	
	_interp->vector_mult(_curr_NV_ptr[0],_dirichletFlag[0]);

	if (_center_x)
		centerVariable(0);
	if (_center_y)
		centerVariable(1);
	if (_center_z)
		centerVariable(2);

 	std::cout<<"AssembleFreqPoroelasticityInclusion::solve() stop\n";
}

void AssembleFreqPoroelasticityInclusion::centerVariable(int _variable)
{
	PetscVector<Number> max_real_vec(_comm);
	PetscVector<Number> min_real_vec(_comm);
	PetscVector<Number> max_imag_vec(_comm);
	PetscVector<Number> min_imag_vec(_comm);
	
	max_real_vec.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
	min_real_vec.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
	max_imag_vec.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
	min_imag_vec.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
	
	max_real_vec=-1e9;
	max_imag_vec=-1e9;
	min_real_vec=1e9;
	min_imag_vec=1e9;

	std::vector<unsigned int> dof_indices;
	MeshBase::const_node_iterator           nd = _mesh.local_nodes_begin ();
	MeshBase::const_node_iterator const end_nd = _mesh.local_nodes_end   ();

	for ( ; nd != end_nd; ++nd)
	{
		Node const * node = *nd;
		_dof_map.dof_indices (node, dof_indices, _d_var.at(_variable) );
		if (dof_indices.size() != 1)
		{
			mooseError("dof_indices.size() != 1");
		}
		Number temp=_curr_NV_ptr[0]( dof_indices.at(0) );
		max_real_vec.set(dof_indices.at(0) , temp.real() );
		min_real_vec.set(dof_indices.at(0) , temp.real() );
		max_imag_vec.set(dof_indices.at(0) , temp.imag() );
		min_imag_vec.set(dof_indices.at(0) , temp.imag() );
	}

 	max_real_vec.close();
 	min_real_vec.close();
 	max_imag_vec.close();
 	min_imag_vec.close();
	
 	Real maxReal=max_real_vec.max();
 	Real minReal=min_real_vec.min();
 	Real maxImag=max_imag_vec.max();
 	Real minImag=min_imag_vec.min();
		
 	Number shift( -0.5*(maxReal+minReal) , -0.5*(maxImag+minImag) );
	
 	max_real_vec.zero();
	
 	nd = _mesh.local_nodes_begin (); 	
 	for ( ; nd != end_nd; ++nd)
 	{
 		Node const * node = *nd;
 		_dof_map.dof_indices (node, dof_indices, _d_var.at(_variable) );
 		if (dof_indices.size() != 1)
 		{
 			mooseError("dof_indices.size() != 1");
 		}
 		max_real_vec.set( dof_indices.at(0), shift );
 	}
 	max_real_vec.close();
	
 	_curr_NV_ptr[0].add(max_real_vec);
	//sol_ptr[0].close();
}

void AssembleFreqPoroelasticityInclusion::mallocLinearSystem()
{
	std::cout<<"AssembleFreqPoroelasticityInclusion::mallocLinearSystem() start\n";
	_static_stiffness   = new PetscMatrix<Number>(_comm);
	_freqdep_stiffness  = new PetscMatrix<Number>(_comm);
	_interp             = new PetscMatrix<Number>(_comm);

	_static_stiffness  ->attach_dof_map(_dof_map);
	_freqdep_stiffness ->attach_dof_map(_dof_map);
	_interp            ->attach_dof_map(_dof_map);

	_static_stiffness  ->init();
	_freqdep_stiffness ->init();
	_interp            ->init();

	_dirichletFlag   = new PetscVector<Number>(_comm);
	_liftingFunction = new PetscVector<Number>(_comm);
	_dirichletFlag   ->init( _dof_map.n_dofs() , _dof_map.n_local_dofs() , false, _curr_NV_ptr->type());
	_liftingFunction ->init( _dof_map.n_dofs() , _dof_map.n_local_dofs() , false, _curr_NV_ptr->type());

	std::cout<<"AssembleFreqPoroelasticityInclusion::mallocLinearSystem() stop\n";	
}

void AssembleFreqPoroelasticityInclusion::deleteLinearSystem()
{
	std::cout<<"AssembleFreqPoroelasticityInclusion::deleteLinearSystem() start\n";
	
	delete _static_stiffness;
	delete _freqdep_stiffness;
	delete _interp;

	delete _dirichletFlag;
	delete _liftingFunction;

	std::cout<<"AssembleFreqPoroelasticityInclusion::deleteLinearSystem() stop\n";
}
