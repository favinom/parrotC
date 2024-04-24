//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleFreqPoroelasticityInclusionPeriodicBC.h"

//registerMooseObject("parrotApp", AssembleFreqPoroelasticityInclusionPeriodicBC);

template <>
InputParameters
	validParams<AssembleFreqPoroelasticityInclusionPeriodicBC>()
{
	InputParameters params = validParams<AssembleFreqPoroelasticityInclusion>();
	params.addRequiredParam< std::vector<Real> >("bc_value_real","bc_value_real.");
	params.addRequiredParam< std::vector<Real> >("bc_value_imag","bc_value_imag.");

	MooseEnum direction("x y z");
	params.addRequiredParam<MooseEnum>("direction",direction,"the direction you want to compute");

	return params;
}

AssembleFreqPoroelasticityInclusionPeriodicBC::AssembleFreqPoroelasticityInclusionPeriodicBC(const InputParameters & parameters) :
AssembleFreqPoroelasticityInclusion(parameters),
_bc_value_real( getParam< std::vector<Real> >("bc_value_real") ),
_bc_value_imag( getParam< std::vector<Real> >("bc_value_imag") ),
_direction(getParam< MooseEnum >("direction") )
{
	if (_bc_value_real.size() != _dim)
	{
		mooseError("_bc_value_real.size() != _dim");
	}
	if (_bc_value_imag.size() != _dim)
	{
		mooseError("_bc_value_imag.size() != _dim");
	}
}

void AssembleFreqPoroelasticityInclusionPeriodicBC::assembleLinearSystem()
{
	std::cout<<"AssembleFreqPoroelasticityInclusionPeriodicBC::assembleLinearSystem() start\n";
	auto start = high_resolution_clock::now();

	RealVectorValue maxV;
	RealVectorValue minV;

	{
		PetscVector<Number> maxCoordinate(_comm);
		PetscVector<Number> minCoordinate(_comm);

		maxCoordinate.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() , false, _curr_NV_ptr->type());
		minCoordinate.init( _dof_map.n_dofs() , _dof_map.n_local_dofs() , false, _curr_NV_ptr->type());

		for (int d=0; d<_dim; ++d)
		{
			maxCoordinate=-1e9 ;
	    	minCoordinate= 1e9 ;

		    MeshBase::const_node_iterator       node_it     = _mesh.local_nodes_begin();
    	 	MeshBase::const_node_iterator const node_it_end = _mesh.local_nodes_end();
     		for ( ; node_it != node_it_end; ++node_it)
     		{
        		const Node * node = *node_it;
        		std::vector<dof_id_type> dof;
        		// is the variable 0 a P1 variable?
        		_dof_map.dof_indices (node, dof,0);
        		if (dof.size()!=1)
        		{
        			mooseError("dof.size()!=1");
        		}
        		maxCoordinate.set(dof.at(0) , node[0](d));
        		minCoordinate.set(dof.at(0) , node[0](d));
    		}
    		maxCoordinate.close();
    		minCoordinate.close();
    		maxV(d)=maxCoordinate.max();
    		minV(d)=minCoordinate.min();
		}

	}

	//std::cout<<maxV<<std::endl;
	//std::cout<<minV<<std::endl;

	_liftingFunction[0] = 0.0;

	std::vector<Number> value_bc(_dim);
	for (int i=0; i<_dim; ++i)
	{
		value_bc.at(i)=Number(_bc_value_real.at(i),_bc_value_imag.at(i));
	}

	{

		MeshBase::const_node_iterator       node_it     = _mesh.local_nodes_begin();
		MeshBase::const_node_iterator const node_it_end = _mesh.local_nodes_end();
		for ( ; node_it != node_it_end; ++node_it)
		{
			const Node * node = *node_it;
			std::vector<dof_id_type> dof;

			for (int v=0; v<_dim; ++v)
			{
				_dof_map.dof_indices (node, dof,_d_var.at(v));
				if (dof.size()!=1)
				{
					mooseError("dof.size()!=1");
				}
				Real temp=(node[0](_direction)-minV(_direction)) / (maxV(_direction)-minV(_direction));
				_liftingFunction[0].set(dof.at(0) ,  value_bc.at(v)*temp );
			}
		}
		_liftingFunction[0].close();
	}

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
	
	Material const & mat = getMaterialByName(_materialName.c_str(),true);
	FreqPoroelasticInclusion const & freqPoroelasticInclusion(dynamic_cast<FreqPoroelasticInclusion const &>(mat));

 	UniquePtr<FEBase> fe_disp (FEBase::build(_dim, fe_disp_type));
 	UniquePtr<FEBase> fe_pres (FEBase::build(_dim, fe_pres_type));
 	fe_disp->attach_quadrature_rule (qrule.get());
 	fe_pres->attach_quadrature_rule (qrule.get());
	
 	const std::vector<Real>& JxW = fe_disp->get_JxW();
 	const std::vector<std::vector<RealGradient> > & dphi     = fe_disp->get_dphi();
 	const std::vector<std::vector<Real> >         &  psi     = fe_pres->get_phi();
 	const std::vector<std::vector<RealGradient> > & dpsi     = fe_pres->get_dphi();
 	const std::vector<Point>                      &  q_point = fe_disp->get_xyz();
	
 	std::vector<dof_id_type>                dof_indices;
	std::vector<dof_id_type> 				dof_indices_interp;
 	std::vector< std::vector<dof_id_type> > dof_indices_disp;
 	std::vector<dof_id_type> 				dof_indices_pres;
 	
 	dof_indices_disp.resize(_dim);
	
 	DenseMatrix<Number> Ke_nonfreq;
 	DenseMatrix<Number> Ke_freq;
 	DenseMatrix<Number> Ke_interp;
	
 	std::vector< std::vector< DenseSubMatrix<Number> * > > A;
 	std::vector< DenseSubMatrix<Number> * >  B;
 	std::vector< DenseSubMatrix<Number> * >  BT;
	
 	DenseSubMatrix<Number> * M;
	
 	A.resize(_dim);
 	B.resize(_dim);
 	BT.resize(_dim);
	
 	for (int i=0; i<_dim; ++i)
 	{
 		A.at(i).resize(_dim);
 	}

 	for (int i=0; i<_dim; ++i)
 	{
 		for (int j=0; j<_dim; ++j)
 		{
 			A.at(i).at(j)=new DenseSubMatrix<Number>(Ke_nonfreq);
 		}
 		B.at(i)=new DenseSubMatrix<Number>(Ke_nonfreq);
 		BT.at(i)=new DenseSubMatrix<Number>(Ke_nonfreq);
 	}
 	M=new DenseSubMatrix<Number>(Ke_nonfreq);

	std::vector<Number> locVec;
        //_liftingFunction[0].get(dof_indices, values);
        _liftingFunction[0].localize(locVec);


 	MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
 	MeshBase::const_element_iterator const end_el = _mesh.active_local_elements_end();

 	for ( ; el != end_el; ++el)
 	{
 		Elem const * elem = *el;
 		fe_disp->reinit (elem);
 		fe_pres->reinit (elem);
		
 		_dof_map.dof_indices (elem, dof_indices);
 		dof_indices_interp=dof_indices;
 		for (int i=0;i<_dim; ++i)
 		{
 			int temp=_d_var.at(i);
 			_dof_map.dof_indices (elem, dof_indices_disp.at(i),temp);
 		}
 		_dof_map.dof_indices (elem, dof_indices_pres,_p_var);		

 		const unsigned int n_dofs   = dof_indices.size();
 		std::vector< unsigned int > n_d_dofs;
 		n_d_dofs.resize(_dim);
 		for (int i=0; i<_dim; ++i)
 		{
 			n_d_dofs.at(i)=dof_indices_disp.at(i).size();
 		}
 		unsigned int const n_p_dofs=dof_indices_pres.size();

 		Ke_nonfreq.resize (n_dofs, n_dofs);
 		Ke_nonfreq.zero();
 		Ke_freq.resize(n_p_dofs,n_p_dofs);
 		Ke_freq.zero();
		
 		for (int i=0; i<_dim; ++i)
 			for (int j=0; j<_dim; ++j)
 				A.at(i).at(j)[0].reposition (_d_var.at(i)* n_d_dofs.at(i) , _d_var.at(j)* n_d_dofs.at(j), n_d_dofs.at(i), n_d_dofs.at(j));
		
 		for (int i=0; i<_dim; ++i)
 		{
 			BT.at(i)[0].reposition ( _d_var.at(i)* n_d_dofs.at(i), _p_var       * n_p_dofs      , n_d_dofs.at(i), n_p_dofs       );
 			B .at(i)[0].reposition ( _p_var      * n_p_dofs      , _d_var.at(i) * n_d_dofs.at(i), n_p_dofs      , n_d_dofs.at(i) );
 		}
		
 		M[0].reposition ( _p_var* n_p_dofs, _p_var* n_p_dofs, n_p_dofs, n_p_dofs );
		
 		std::vector<Real> mu, la, al, di, om;
 		mu.resize(q_point.size());
 		la.resize(q_point.size());
 		al.resize(q_point.size());
 		di.resize(q_point.size());
 		om.resize(q_point.size());
		
 		for (int qp=0; qp<q_point.size();++qp)
 		{
 			Point const & temp=q_point.at(qp);
 			freqPoroelasticInclusion.computeQpProperties(temp,mu.at(qp),la.at(qp),al.at(qp),di.at(qp),om.at(qp));
 		}
 				
 		std::vector< std::vector< std::vector< RealTensorValue > > > E;
 		E.resize(_dim);
 		for (int d=0; d<_dim;++d)
 		{
 			E.at(d).resize( dphi.size() );
 		}
 		for (int d=0; d<_dim;++d)
 		{
 			for (int i=0; i<dphi.size(); ++i)
 			{
 				E.at(d).at(i).resize( qrule->n_points() );				
 			}

 		}
 		for (int d=0; d<_dim;++d)
 		{
 			for (int i=0; i<dphi.size(); ++i)
 			{
 				for (int qp=0; qp<qrule->n_points(); ++qp)
 				{
 					E.at(d).at(i).at(qp)=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
 					for (int c=0; c<_dim; ++c)
 					{
 						E.at(d).at(i).at(qp)(d,c)=dphi.at(i).at(qp)(c);
 					}
 					RealTensorValue & Et=E.at(d).at(i).at(qp);
 					Et=0.5*(Et+Et.transpose());
 				}				
 			}
 		}

		for (int id=0; id<_dim; ++id)
			for (int jd=0; jd<_dim; ++jd)
				for (int i=0; i<dphi.size(); ++i)
					for (int j=0; j<dphi.size(); ++j)
						for (int qp=0; qp<qrule->n_points(); ++qp)
		{
			RealTensorValue const & H=E.at(jd).at(j).at(qp);
			RealTensorValue const & V=E.at(id).at(i).at(qp);

			RealTensorValue sigma=2.0 * mu[qp]*H + la[qp]*H.tr()*_identity;

			A.at(id).at(jd)[0](i,j)+=JxW[qp]*sigma.contract(V);
		}

		
		for (int id=0; id<_dim; ++id)
				for (int i=0; i<dphi.size(); ++i)
					for (int j=0; j<psi.size(); ++j)
						for (int qp=0; qp<qrule->n_points(); ++qp)
		{
			RealTensorValue const & V=E.at(id).at(i).at(qp);
			BT.at(id)[0](i,j)+= -      1.0*JxW[qp]*( al[qp]*psi[j][qp]*V.tr() );
		}
		

		for (int id=0; id<_dim; ++id)
		{
			for (int i=0; i<BT.at(id)[0].m(); ++i)
			{
				for (int j=0; j<BT.at(id)[0].n(); ++j)
				{
					B.at(id)[0](j,i) = _imagUnit*BT.at(id)[0](i,j);
				}
			}
		}
		


		for (int i=0;i<psi.size();++i)
		{
			for (int j=0;j<psi.size();++j)
			{
				for (int qp=0; qp<qrule->n_points(); ++qp)
				{
					M[0](i,j)+=-_imagUnit*JxW[qp]*om[qp]* psi[i][qp]* psi[j][qp];
					Ke_freq(i,j)+=   -      1.0*JxW[qp]*di[qp]*dpsi[i][qp]*dpsi[j][qp];
				}
			}
		}
		
		DenseMatrix<Number> Ke_constr;
 		Ke_constr.resize(Ke_nonfreq.m(),Ke_nonfreq.n());
 		Ke_constr.zero();

		///////////////////////////////////////////////////////
		//std::vector<Number> values;
		Number * values;
		values=new Number[dof_indices.size()];
		
		for (unsigned ii=0; ii<dof_indices.size(); ++ii)
                {
			unsigned idl=dof_indices.at(ii);
			values[ii]=locVec.at(idl);
                }

		DenseVector<Number> lrhs;
		lrhs.resize(dof_indices.size());
		lrhs.zero();
		
		for (int r=0; r<Ke_nonfreq.m();++r)
		{
			for (int c=0; c<Ke_nonfreq.n();++c)
			{
				lrhs(r)+=Ke_nonfreq(r,c)*values[c];
			}
			lrhs(r)=-1.0*lrhs(r);
		}
		
		delete values;
		///////////////////////////////////////////////////////

 		_dof_map.constrain_element_matrix_and_vector (Ke_nonfreq , lrhs , dof_indices         , false   );
 		_dof_map.constrain_element_matrix            (Ke_freq           , dof_indices_pres    , false   );
 		_dof_map.constrain_element_matrix            (Ke_constr         , dof_indices_interp            );

 		Ke_interp.resize(dof_indices_interp.size(),dof_indices_interp.size());
		Ke_interp.zero();
				
		for (int i=0; i<Ke_constr.m(); ++i)
		{
			if (Ke_constr(i,i).real()>0.5)
			{
				for (int j=0; j<Ke_constr.n(); ++j)
				{
					if (i!=j)
					{
						Ke_interp(i,j)=-1.0*Ke_constr(i,j);
					}
				}
			}
			else
			{
				Ke_interp(i,i)=1.0;
			}
		}
		
		_static_stiffness  ->add_matrix (Ke_nonfreq , dof_indices      );
		_mat_SM_ptr		   ->add_matrix (Ke_nonfreq , dof_indices      );
 		_freqdep_stiffness ->add_matrix (Ke_freq    , dof_indices_pres );
 		_rhs_NV_ptr        ->add_vector (lrhs, dof_indices      );
		
 		for (int i=0; i<dof_indices.size(); ++i)
 		{
 			for (int j=0; j<dof_indices.size(); ++j)
 			{
 				_interp     ->set( dof_indices.at(i), dof_indices.at(j), Ke_interp(i,j) );
 			}			
 		}
 	} // end element loop

 	_static_stiffness  ->close();
 	_freqdep_stiffness ->close();
 	_mat_SM_ptr        ->close();
 	_interp            ->close();
 	_rhs_NV_ptr        ->close();

	std::vector<unsigned int> toDelete;
	for ( int i=_dirichletFlag->first_local_index(); i<_dirichletFlag->last_local_index(); ++i )
	{
		if ( _dirichletFlag[0](i).real()>0.5 )
		{
			toDelete.push_back(i);
		}
	}

 	_static_stiffness ->zero_rows(toDelete,1.0);
 	// we can skip the substitution of the diagonal element
 	_freqdep_stiffness->zero_rows(toDelete);
 
 	auto stop = high_resolution_clock::now();
 	auto duration = duration_cast<microseconds>(stop - start);

 	myout << "Assembly matrices took "<<duration.count()/1000.0/1000.0 <<" s"<<std::endl;
	std::cout<<"AssembleFreqPoroelasticityInclusionPeriodicBC::assembleLinearSystem() stop\n";
}

void AssembleFreqPoroelasticityInclusionPeriodicBC::solve( Real _t_in )
{
	std::cout<<"AssembleFreqPoroelasticityInclusionPeriodicBC::solve() start\n";
	 	auto start = high_resolution_clock::now(); 

 	_mat_SM_ptr->zero();
 	_mat_SM_ptr->add(1.0,_static_stiffness[0]);
	Real scale=1.0/(2.0*_pi*std::pow(10.0,_t_in));
 	_mat_SM_ptr->add(scale,_freqdep_stiffness[0]);

 	auto stop = high_resolution_clock::now(); 
 	auto duration = duration_cast<microseconds>(stop - start);

 	myout << "Summing matrices took "<<duration.count()/1000.0/1000.0 <<" s"<< std::endl;

 	_parrotSolver.setMatrixAndVectors( _mat_SM_ptr , _rhs_NV_ptr , _curr_NV_ptr );
 	_parrotSolver.solve();

	_dirichletFlag[0]=_curr_NV_ptr[0];
	
	_interp->vector_mult(_curr_NV_ptr[0],_dirichletFlag[0]);

	_curr_NV_ptr[0].add(_liftingFunction[0]);

	if (_center_x)
		centerVariable(0);
	if (_center_y)
		centerVariable(1);
	if (_center_z)
		centerVariable(2);

	std::cout<<"AssembleFreqPoroelasticityInclusionPeriodicBC::solve() stop\n";
		//_curr_NV_ptr[0]=_liftingFunction[0];
}
