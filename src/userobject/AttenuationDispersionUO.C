//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AttenuationDispersionUO.h"

#include <chrono>
using namespace std::chrono;

#include "FEProblem.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "Material.h"

#include "ParrotSolver.h"

#include "FreqPoroelasticInclusion.h"

#include "libmesh/dof_map.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/dirichlet_boundaries.h"

//registerMooseObject("parrotApp", AttenuationDispersionUO);

template <>
InputParameters
	validParams<AttenuationDispersionUO>()
{
	InputParameters params = validParams<GeneralUserObject>();
	params.addRequiredParam< std::string >("material_name","the user object where assembly is performed.");
	return params;
}

AttenuationDispersionUO::AttenuationDispersionUO(const InputParameters & parameters) :
GeneralUserObject(parameters),
_equationSystems(_fe_problem.es()),
_tnis( _equationSystems.get_system<TransientNonlinearImplicitSystem>("nl0") ),
_mesh( _equationSystems.get_mesh() ),
_comm ( _mesh.comm() ),
_dof_map( _tnis.get_dof_map() ),
_dim( _mesh.mesh_dimension() ),
_nvars( _tnis.n_vars() ),
_nl( _fe_problem.getNonlinearSystemBase() ),
_mat_SM(*_tnis.matrix),
//_mat_PM( dynamic_cast<PetscMatrix<Number> &>(_mat_SM) ),
_rhs_NV(*_tnis.rhs),
_sol_NV( _nl.solution() ),
_postprocessorVectorsAllocated(false),
_postprocessorVectorsAssembled(false),
_pi( 4.0*std::atan(1.0) ),
_imagUnit(0.0,1.0),
_materialName( getParam< std::string >("material_name") )
{
	std::cout<<"AttenuationDispersionUO::AttenuationDispersionUO start\n";
	
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

	std::cout<<"AttenuationDispersionUO::AttenuationDispersionUO stop\n";
}

void AttenuationDispersionUO::initialize()
{
	std::cout<<"AttenuationDispersionUO::initialize() start\n";

	if (!_postprocessorVectorsAllocated)
	{
		mallocPostprocessorVectors();
		_postprocessorVectorsAllocated=true;
	}
	if (!_postprocessorVectorsAssembled)
	{
		assemblePostprocessorVectors();
		_postprocessorVectorsAssembled=true;
	}

	std::cout<<"AttenuationDispersionUO::initialize() stop\n";
}

void AttenuationDispersionUO::execute()
{
	std::cout<<"AttenuationDispersionUO::execute() start\n";		
	NumericVector<Number> const * const & solution=_nl.currentSolution();

	// Be careful, here we could be computing useless values :-/
	for (int i=0; i<_dim; ++i)
	{
		for (int j=0; j<_dim; ++j)
		{
			_strainIntegral(i,j)=solution[0].dot(_strainVector[i][j][0]);
			_stressIntegral(i,j)=solution[0].dot(_stressVector[i][j][0]);
		}
	}

	// std::cout<<_strainVector[0][0][0].max()<<std::endl;
	// std::cout<<_strainVector[0][0][0].min()<<std::endl;
	// std::cout<<_strainIntegral<<std::endl;

	std::cout<<"AttenuationDispersionUO::execute() stop\n";
}

AttenuationDispersionUO::~AttenuationDispersionUO()
{
	std::cout<<"AttenuationDispersionUO::~AttenuationDispersionUO() start\n";
	deletePostprocessorVectors();
	std::cout<<"AttenuationDispersionUO::~AttenuationDispersionUO() stop\n";
}

void AttenuationDispersionUO::assemblePostprocessorVectors()
{
	std::cout<<"AttenuationDispersionUO::assemblePostprocessorVectors() start\n";
	auto start = high_resolution_clock::now(); 

	for (int i=0; i<_dim; ++i)
	{
		for (int j=0 ;j<_dim; ++j)
		{
			_strainVector[i][j][0]=0.0;
			_stressVector[i][j][0]=0.0;
		}
	}
	
	// We need this to check if we are into an inclusion/fracture
	Material const & mat = getMaterialByName(_materialName.c_str(),true);
	FreqPoroelasticInclusion const & freqPoroelasticInclusion( dynamic_cast<FreqPoroelasticInclusion const &>(mat) );
	
	QBase const * const & qbase(_assembly.qRule());
	std::unique_ptr<QBase> qrule( QBase::build (qbase->type(),_dim,qbase->get_order()));
	
	FEType const fe_disp_type = _tnis.variable_type(_d_var.at(0));
	FEType const fe_pres_type = _tnis.variable_type(_p_var);
		
	{
		FEType const temp = _tnis.variable_type(_d_var.at(1));
		if (fe_disp_type != temp)
		{
			mooseError("FEType of disp_x is different from disp_y!");
		}
	}
	if (_dim==3)
	{
		FEType const temp = _tnis.variable_type(_d_var.at(2));
		if (fe_disp_type != temp)
		{
			mooseError("FEType of disp_x is different from disp_z!");
		}
	}
	
	UniquePtr<FEBase> fe_disp (FEBase::build(_dim, fe_disp_type));
	UniquePtr<FEBase> fe_pres (FEBase::build(_dim, fe_pres_type));
	fe_disp->attach_quadrature_rule (qrule.get());
	fe_pres->attach_quadrature_rule (qrule.get());
	
	const std::vector<Real>& JxW = fe_disp->get_JxW();
	const std::vector<std::vector<RealGradient> > & dphi     = fe_disp->get_dphi();
	const std::vector<std::vector<Real> >         &  psi     = fe_pres->get_phi();
	const std::vector<Point>                      &  q_point = fe_disp->get_xyz();
	
	std::vector<dof_id_type>                dof_indices;
	std::vector< std::vector<dof_id_type> > dof_indices_disp;
	std::vector<dof_id_type> 				dof_indices_pres;

	dof_indices_disp.resize(_dim);
	
	MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
	MeshBase::const_element_iterator const end_el = _mesh.active_local_elements_end();
	
	for ( ; el != end_el; ++el)
	{
		Elem const * elem = *el;
		_dof_map.dof_indices (elem, dof_indices);
		for (int i=0;i<_dim; ++i)
		{
			_dof_map.dof_indices (elem, dof_indices_disp.at(i),_d_var.at(i));
		}
		_dof_map.dof_indices (elem, dof_indices_pres,_p_var);
				
		fe_disp->reinit (elem);
		fe_pres->reinit (elem);
		
		std::vector<Real> mu, la, al;
		mu.resize(q_point.size());
		la.resize(q_point.size());
		al.resize(q_point.size());
		
 		for (int qp=0; qp<q_point.size();++qp)
 		{
 			Point const & temp=q_point.at(qp);
 			Real dummy1, dummy2;
 			freqPoroelasticInclusion.computeQpProperties(temp,mu.at(qp),la.at(qp),al.at(qp),dummy1,dummy2);
 		}

		DenseVector<Number> Re_i;
		DenseVector<Number> Re_j;

		DenseVector<Number> shear_re_i;
		DenseVector<Number> shear_re_j;
		
		std::vector< DenseVector<Number> > bulk_re;
		bulk_re.resize(_dim);
		
		DenseVector<Number> pres_re;

		for (int i=0; i<_dim; ++i)
		{
			for (int j=0; j<_dim; ++j)
			{
				// These two numbers should be the same with our hypotheses
				int const loc_n_i=dof_indices_disp.at(i).size();
				int const loc_n_j=dof_indices_disp.at(j).size();
				
				int const loc_p=dof_indices_pres.size();
		
				Re_i.resize(loc_n_i);
				Re_j.resize(loc_n_j);
				Re_i.zero();
				Re_j.zero();
					
				shear_re_i.resize(loc_n_i);
				shear_re_j.resize(loc_n_j);
				shear_re_i.zero();
				shear_re_j.zero();
				
				pres_re.resize(loc_p);
				pres_re.zero();
				
				for (int d=0; d<_dim; ++d)
				{
					bulk_re.at(d).resize(loc_n_i);
					bulk_re.at(d).zero();
				}
				
				for (unsigned int ip=0; ip<dphi.size(); ip++)
				{
					for (unsigned int qp=0; qp<qrule->n_points(); qp++)
					{
						Real const v_i=0.5*JxW[qp] * dphi[ip][qp](j);
						Real const v_j=0.5*JxW[qp] * dphi[ip][qp](i);
						Re_i(ip) += v_i;
						Re_j(ip) += v_j;
						shear_re_i(ip)+=2.0*mu.at(qp)*v_i;
						shear_re_j(ip)+=2.0*mu.at(qp)*v_j;
						for (int d=0; d<_dim; ++d)
						{
							bulk_re.at(d)(ip)+=la.at(qp)*JxW[qp]*dphi[ip][qp](d);
						}
					}
				}
				
				for (unsigned int ip=0; ip<psi.size(); ip++)
				{
					for (unsigned int qp=0; qp<qrule->n_points(); qp++)
					{
						pres_re(ip)-=al.at(qp) * JxW[qp]*psi[ip][qp];
					}
				}
						
				_strainVector[i][j][0].add_vector(Re_i, dof_indices_disp.at(i));
				_strainVector[i][j][0].add_vector(Re_j, dof_indices_disp.at(j));
				_stressVector[i][j][0].add_vector(shear_re_i, dof_indices_disp.at(i));
				_stressVector[i][j][0].add_vector(shear_re_j, dof_indices_disp.at(j));
				if (i==j)
				{
					for (int d=0; d<_dim; ++d)
					{
						_stressVector[i][j][0].add_vector(bulk_re.at(d), dof_indices_disp.at(d));
					}
					_stressVector[i][j][0].add_vector( pres_re, dof_indices_pres );
				}
			}
		}
	}
	
	for (int i=0; i<_dim; ++i)
	{
		for (int j=0 ;j<_dim; ++j)
		{
			_strainVector[i][j][0].close();
			_stressVector[i][j][0].close();
		}
	}

auto stop = high_resolution_clock::now(); 
auto duration = duration_cast<microseconds>(stop - start);
std::cout << "Assembly AttenuationDispersionUO arrays took "<<duration.count() << std::endl;
std::cout<<"AttenuationDispersionUO::assemblePostprocessorVectors() stop\n";

	
}

Number const AttenuationDispersionUO::getStrainComponent(int i, int j) const
{
	return _strainIntegral(i,j);
}

Number const AttenuationDispersionUO::getStressComponent(int i, int j) const
{
	return _stressIntegral(i,j);
}

void AttenuationDispersionUO::mallocPostprocessorVectors()
{
	std::cout<<"AttenuationDispersionUO::mallocPostprocessorVectors() start\n";
	if ( !_postprocessorVectorsAllocated )
	{
		_strainVector.resize(_dim);
		_stressVector.resize(_dim);
		for (int i=0; i<_dim; ++i)
		{
			_strainVector.at(i).resize(_dim);
			_stressVector.at(i).resize(_dim);
		}
		for (int i=0; i<_dim; ++i)
		{
			for (int j=0; j<_dim; ++j)
			{
				_strainVector[i][j]=new PetscVector<Number>(_comm);
				_stressVector[i][j]=new PetscVector<Number>(_comm);
				_strainVector[i][j][0].init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
				_stressVector[i][j][0].init( _dof_map.n_dofs() , _dof_map.n_local_dofs() );
			}
		}
	}
	std::cout<<"AttenuationDispersionUO::mallocPostprocessorVectors() stop\n";
}

void AttenuationDispersionUO::deletePostprocessorVectors()
{
	std::cout<<"AttenuationDispersionUO::deletePostprocessorVectors() start\n";
	if ( _postprocessorVectorsAllocated )
	{
		for (int i=0; i<_dim; ++i)
		{
			for (int j=0; j<_dim; ++j)
			{
				delete _strainVector[i][j];
				delete _stressVector[i][j];
			}
		}
		_postprocessorVectorsAllocated=false;
	}
	std::cout<<"AttenuationDispersionUO::deletePostprocessorVectors() stop\n";
}
