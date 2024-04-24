//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExportStressStrain.h"

#include "FEProblem.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "Material.h"

#include "ParrotSolver.h"

#include "FreqPoroelasticInclusion.h"

#include "libmesh/dof_map.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"

//typedef TransientExplicitSystem SystemClass;
typedef ExplicitSystem SystemClass;

//registerMooseObject("parrotApp", ExportStressStrain);

template <>
InputParameters
	validParams<ExportStressStrain>()
{
	InputParameters params = validParams<GeneralUserObject>();
	params.addRequiredParam< std::string >("output_filename","The name to be given to the output file.");
	params.addRequiredParam< std::string >("material_name","the user object where assembly is performed.");

	return params;
}

ExportStressStrain::ExportStressStrain(const InputParameters & parameters) :
GeneralUserObject(parameters),
_equationSystems(_fe_problem.es()),
_tnis( _equationSystems.get_system<TransientNonlinearImplicitSystem>("nl0") ),
_mesh( _equationSystems.get_mesh() ),
_comm ( _mesh.comm() ),
_dof_map( _tnis.get_dof_map() ),
_dim( _mesh.mesh_dimension() ),
_nvars( _tnis.n_vars() ),
_nl( _fe_problem.getNonlinearSystemBase() ),
_allocated(false),
_assembled(false),
_t_step(0),
_output_filename( getParam<std::string >("output_filename") ),
_pi( 4.0*std::atan(1.0) ),
_imagUnit(0.0,1.0),
_materialName( getParam< std::string >("material_name") )
{
	std::cout<<"ExportStressStrain::ExportStressStrain() start\n";
	
	// This error should never happen as the system nl0 has been retrieved above
	if( !_equationSystems.has_system("nl0") )
	{
		mooseError("The system nl0 does not exist!");
	}
	
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
	
	strain_name.resize(3);
	stress_name.resize(3);
	for (int i=0; i<3; ++i)
	{
		strain_name.at(i).resize(3);
		stress_name.at(i).resize(3);
	}
	
	strain_name.at(0).at(0)="eps_xx";
	strain_name.at(0).at(1)="eps_xy";
	strain_name.at(0).at(2)="eps_xz";
	strain_name.at(1).at(0)="eps_yx";
	strain_name.at(1).at(1)="eps_yy";
	strain_name.at(1).at(2)="eps_yz";
	strain_name.at(2).at(0)="eps_zx";
	strain_name.at(2).at(1)="eps_zy";
	strain_name.at(2).at(2)="eps_zz";
	
	stress_name.at(0).at(0)="sigma_xx";
	stress_name.at(0).at(1)="sigma_xy";
	stress_name.at(0).at(2)="sigma_xz";
	stress_name.at(1).at(0)="sigma_yx";
	stress_name.at(1).at(1)="sigma_yy";
	stress_name.at(1).at(2)="sigma_yz";
	stress_name.at(2).at(0)="sigma_zx";
	stress_name.at(2).at(1)="sigma_zy";
	stress_name.at(2).at(2)="sigma_zz";

	std::cout<<"ExportStressStrain::ExportStressStrain() stop\n";
}

void ExportStressStrain::initialize()
{
	std::cout<<"ExportStressStrain::initialize() start\n";
	
	UnstructuredMesh const & mesh_ref=dynamic_cast<UnstructuredMesh const &>(_mesh);
		
	if (!_allocated)
	{
		std::cout<<"inside !_allocated"<<std::endl;
		_mesh2=new MeshClass(mesh_ref);
		_equationSystems2=new EquationSystems(_mesh2[0]);
		_mat=new PetscMatrix<Number> (_comm);
		outputClass=new OutputClass(_mesh2[0]);
		_allocated=true;
	}

	if (!_assembled)
	{
		std::cout<<"inside !assembled"<<std::endl;
				
		auto const & mat = getMaterialByName(_materialName,true);
		FreqPoroelasticInclusion const & freqPoroelasticInclusion( dynamic_cast<FreqPoroelasticInclusion const &>(mat) );
	
		SystemClass & explicitSystem2 = _equationSystems2[0].add_system<SystemClass> ("Export");
		
		std::vector< std::vector<unsigned int> > eps_var;
		std::vector< std::vector<unsigned int> > sigma_var;
		eps_var.resize(_dim);
		sigma_var.resize(_dim);
		for (int i=0; i<_dim; ++i)
		{
			eps_var.at(i).resize(_dim);
			sigma_var.at(i).resize(_dim);
		}
	
		for (int i=0; i<_dim; ++i)
		{
			for (int j=0; j<_dim; ++j)
			{
				eps_var.at(i).at(j)=explicitSystem2.add_variable(strain_name.at(i).at(j).c_str(),CONSTANT,MONOMIAL);
			}
		}
	
		for (int i=0; i<_dim; ++i)
		{
			for (int j=0; j<_dim; ++j)
			{
				sigma_var.at(i).at(j)=explicitSystem2.add_variable(stress_name.at(i).at(j).c_str(),CONSTANT,MONOMIAL);
			}
		}
	
		_equationSystems2->init();
		
		outputClass->write_equation_systems ( _output_filename.c_str() , _equationSystems2[0]);
		outputClass->append(true);
		
		QBase const * const & qbase(_assembly.qRule());
		std::unique_ptr<QBase> qrule( QBase::build (qbase->type(),_dim,qbase->get_order()));
	
		FEType const fe_disp_type = _tnis.variable_type(_d_var.at(0));
		FEType const fe_pres_type = _tnis.variable_type(_p_var);
	
		UniquePtr<FEBase> fe_disp (FEBase::build(_dim, fe_disp_type));
		UniquePtr<FEBase> fe_pres (FEBase::build(_dim, fe_pres_type));
		fe_disp->attach_quadrature_rule (qrule.get());
		fe_pres->attach_quadrature_rule (qrule.get());
	
		const std::vector<Real>                       &  JxW     = fe_disp->get_JxW();
		const std::vector<std::vector<RealGradient> > & dphi     = fe_disp->get_dphi();
		const std::vector<std::vector<Real> >         &  psi     = fe_pres->get_phi();
		const std::vector<Point>                      &  q_point = fe_disp->get_xyz();
	
		std::vector<dof_id_type>                dof_indices_i;
		std::vector<dof_id_type>                dof_indices_j;
		std::vector<dof_id_type>                dof_indices_eps;
		std::vector<dof_id_type>                dof_indices_sigma;
		std::vector< std::vector<dof_id_type> > dof_indices_disp;
		std::vector<dof_id_type> 				dof_indices_pres;

		dof_indices_disp.resize(_dim);
	
	
		DofMap & dof_map2=explicitSystem2.get_dof_map();
	
		MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
		MeshBase::const_element_iterator const end_el = _mesh.active_local_elements_end();
	
		MeshBase::const_element_iterator       el2     = _mesh2->active_local_elements_begin();
		MeshBase::const_element_iterator const end_el2 = _mesh2->active_local_elements_end();
	
		_mat[0].init( dof_map2.n_dofs(), _dof_map.n_dofs(), dof_map2.n_local_dofs() , _dof_map.n_local_dofs() );
		//MatSetOption(mat.mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	
		Real area;
	
		for ( ; el != end_el; ++el, ++el2 )
		{
			Elem const * elem = *el;
			Elem const * elem2 = *el2;
			
			fe_disp->reinit (elem);
			fe_pres->reinit (elem);
		
			std::vector<Real> mu, la, al;
			mu.resize(q_point.size());
			la.resize(q_point.size());
			al.resize(q_point.size());
			
			area=0.0;

			for (int qp=0; qp<q_point.size();++qp)
 			{
 				area+=JxW[qp];
 				Point const & temp=q_point.at(qp);
 				Real dummy1, dummy2;
 				freqPoroelasticInclusion.computeQpProperties(temp,mu.at(qp),la.at(qp),al.at(qp),dummy1,dummy2);
 			}
					
			for (int i=0; i<_dim; ++i)
			{
				for (int j=0; j<_dim; ++j)
				{
					_dof_map.dof_indices ( elem , dof_indices_i  ,  _d_var.at(i) );
					_dof_map.dof_indices ( elem , dof_indices_j  ,  _d_var.at(j) );
					dof_map2.dof_indices ( elem2, dof_indices_eps, eps_var.at(i).at(j) );
								
					DenseMatrix<Number> Ke_i;
					DenseMatrix<Number> Ke_j;
					Ke_i.resize(dof_indices_eps.size(),dof_indices_i.size());
					Ke_j.resize(dof_indices_eps.size(),dof_indices_j.size());
				
					Ke_i.zero();
					Ke_j.zero();
				
					for (unsigned int ip=0; ip<dphi.size(); ip++)
					{
						for (unsigned int qp=0; qp<qrule->n_points(); qp++)
						{
							Ke_i(0,ip)+=0.5*JxW[qp] * dphi[ip][qp](j) /area;
							Ke_j(0,ip)+=0.5*JxW[qp] * dphi[ip][qp](i) /area;
						}
					}
				
					_mat[0].add_matrix(Ke_i, dof_indices_eps, dof_indices_i );
					_mat[0].add_matrix(Ke_j, dof_indices_eps, dof_indices_j );
				
				}
			}

			for (int i=0; i<_dim; ++i)
			{
				for (int j=0; j<_dim; ++j)
				{					
					_dof_map.dof_indices ( elem , dof_indices_i     ,    _d_var.at(i) );
					_dof_map.dof_indices ( elem , dof_indices_j     ,    _d_var.at(j) );
					dof_map2.dof_indices ( elem2, dof_indices_sigma , sigma_var.at(i).at(j) );
								
					DenseMatrix<Number> Ke_i;
					DenseMatrix<Number> Ke_j;
					Ke_i.resize(dof_indices_sigma.size(),dof_indices_i.size());
					Ke_j.resize(dof_indices_sigma.size(),dof_indices_j.size());
				
					Ke_i.zero();
					Ke_j.zero();				
									
					for (unsigned int ip=0; ip<dphi.size(); ip++)
					{
						for (unsigned int qp=0; qp<qrule->n_points(); qp++)
						{
							Ke_i(0,ip)+=mu.at(qp) * JxW[qp] * dphi[ip][qp](j) /area;
							Ke_j(0,ip)+=mu.at(qp) * JxW[qp] * dphi[ip][qp](i) /area;
						}
					}				
					
					_mat[0].add_matrix(Ke_i, dof_indices_sigma, dof_indices_i );
					_mat[0].add_matrix(Ke_j, dof_indices_sigma, dof_indices_j );
					
					if (i==j)
					{
						
						std::vector< DenseMatrix<Number> > bulk_re;
						std::vector< std::vector < unsigned int > > dof_indices_bulk;
						bulk_re.resize(_dim);
						dof_indices_bulk.resize(_dim);
						
						for (int d=0; d<_dim; ++d)
						{
							_dof_map.dof_indices ( elem , dof_indices_bulk.at(d) , _d_var.at(d) );
							
							bulk_re.at(d).resize( dof_indices_sigma.size() , dof_indices_bulk.at(d).size() );
							bulk_re.at(d).zero();
							
							for (unsigned int ip=0; ip<dphi.size(); ip++)
							{
								for (unsigned int qp=0; qp<qrule->n_points(); qp++)
								{
									bulk_re.at(d)(0,ip)+=la.at(qp)*JxW[qp]*dphi[ip][qp](d) /area;
								}
							}							
							_mat[0].add_matrix(bulk_re.at(d), dof_indices_sigma, dof_indices_bulk.at(d) );	
						}
						
						
						
						DenseMatrix<Number> pres_re;
						std::vector < unsigned int > dof_indices_pres;
						
						_dof_map.dof_indices ( elem , dof_indices_pres , _p_var );
						
						pres_re.resize( dof_indices_sigma.size() , dof_indices_pres.size()  );
						pres_re.zero();
						
						
						for (unsigned int ip=0; ip<psi.size(); ip++)
						{
							for (unsigned int qp=0; qp<qrule->n_points(); qp++)
							{
								pres_re(0,ip)-=al.at(qp) * JxW[qp]*psi[ip][qp] /area;
							}
						}
												
						_mat[0].add_matrix(pres_re, dof_indices_sigma, dof_indices_pres );	
						
					}// end if (i==j)
					
				}// end loop j
			}// end loop i
			
		}// end loop element
		_mat[0].close();
		_assembled=true;
	}//end if _assembled
		std::cout<<"ExportStressStrain::initialize() stop\n";
}

void ExportStressStrain::execute()
{	
	std::cout<<"ExportStressStrain::execute() start\n";
	_t_step++;
	if (_allocated && _assembled)
	{
		
		std::cout<<"inside writing()\n";
		
		// Here we read the solution of nl0
		std::unique_ptr<NumericVector<Number> > & sol1_uni_ptr(_tnis.solution);
		NumericVector<Number> * sol1_ptr = sol1_uni_ptr.get();
	
		// Here we read the solution of export
		SystemClass & explicitSystem2 = _equationSystems2[0].get_system<TransientExplicitSystem> ("Export");
		std::unique_ptr<NumericVector<Number> > & sol2_uni_ptr(explicitSystem2.solution);
		NumericVector<Number> * sol2_ptr = sol2_uni_ptr.get();		
	
		// we perform multiplication
		_mat[0].vector_mult(sol2_ptr[0],sol1_ptr[0]);
	
		// we write the solution
		outputClass->write_timestep ( _output_filename.c_str() , _equationSystems2[0], _t_step, _fe_problem.time());
		}
	std::cout<<"ExportStressStrain::execute() start\n";
}

ExportStressStrain::~ExportStressStrain()
{
	std::cout<<"ExportStressStrain::~ExportStressStrain() start\n";
	if (_allocated)
	{
		delete outputClass;
		delete _mat;
		delete _equationSystems2;
		delete _mesh2;
		_allocated=false;
	}
	std::cout<<"ExportStressStrain::~ExportStressStrain() stop\n";
}


//for (int ii=0; ii<dof_indices_eps.size(); ++ii)
//{
//	for (int jj=0; jj<dof_indices_i.size(); ++jj)
//	{
//		std::cout<<i<<" "<<j<<" "<<ii<<" "<<jj<<std::endl;
//		mat.set(dof_indices_eps.at(ii),dof_indices_i.at(jj),Ke_i(ii,jj));
//	}
//}
				
//std::cout<<dof_indices_eps.size()<<std::endl;
//std::cout<<dof_indices_i.size()<<std::endl;
//std::cout<<dof_indices_j.size()<<std::endl;
