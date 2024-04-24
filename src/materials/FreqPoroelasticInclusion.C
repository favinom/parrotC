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

#include "FreqPoroelasticInclusion.h"

//#include <sstream>
//#include "MooseMesh.h"

#include "InclusionsMeshModifier.h"

template<>
InputParameters validParams<FreqPoroelasticInclusion>()
{
	InputParameters params = validParams<Material>();
	params.addRequiredParam<Real>("mu_block", "second Lamé parameter");
	params.addRequiredParam<Real>("lambda_block", "first Lamé parameter");
	params.addRequiredParam<Real>("alpha_block", "another material parameter related to pressure");
	params.addRequiredParam<Real>("kappa_block", "permeability");
	params.addRequiredParam<Real>("eta_block", "viscosity");
	params.addRequiredParam<Real>("porosity_block", "porosity");
	params.addRequiredParam<Real>("kf_block", "bulk_modulus_fluid");
	params.addRequiredParam<Real>("ks_block", "bulk_modulus_solid");
	params.addRequiredParam< std::vector<Real> >("mu_fracture", "second Lamé parameter");
	params.addRequiredParam< std::vector<Real> >("lambda_fracture", "first Lamé parameter");
	params.addRequiredParam< std::vector<Real> >("alpha_fracture", "another material parameter related to pressure");
	params.addRequiredParam< std::vector<Real> >("kappa_fracture", "permeability");
	params.addRequiredParam< std::vector<Real> >("eta_fracture", "viscosity");
	params.addRequiredParam< std::vector<Real> >("porosity_fracture", "porosity");
	params.addRequiredParam< std::vector<Real> >("kf_fracture", "bulk_modulus_fluid");
	params.addRequiredParam< std::vector<Real> >("ks_fracture", "bulk_modulus_solid");

	params.addCoupledVar("disp_x"  , "complex first  disp component");
	params.addCoupledVar("disp_y"  , "complex second disp component");
	params.addCoupledVar("disp_z"  , "complex third  disp component");
	params.addCoupledVar("pressure", "complex pressure"             );
	
	params.addParam<std::string>("inclusion_meshmodifier","the name of the inclusion meshmodifier");

	return params;
}

FreqPoroelasticInclusion::FreqPoroelasticInclusion(const InputParameters & parameters) :
Material(parameters),
_mu_block(getParam<Real>("mu_block")),
_lambda_block(getParam<Real>("lambda_block")),
_alpha_block(getParam<Real>("alpha_block")),
_kappa_block(getParam<Real>("kappa_block")),
_eta_block(getParam<Real>("eta_block")),
_porosity_block(getParam<Real>("porosity_block")),
_kf_block(getParam<Real>("kf_block")),
_ks_block(getParam<Real>("ks_block")),
_mu_fracture(getParam< std::vector<Real> >("mu_fracture")),
_lambda_fracture(getParam< std::vector<Real> >("lambda_fracture")),
_alpha_fracture(getParam< std::vector<Real> >("alpha_fracture")),
_kappa_fracture(getParam< std::vector<Real> >("kappa_fracture")),
_eta_fracture(getParam< std::vector<Real> >("eta_fracture")),
_porosity_fracture(getParam< std::vector<Real> >("porosity_fracture")),
_kf_fracture(getParam< std::vector<Real> >("kf_fracture")),
_ks_fracture(getParam< std::vector<Real> >("ks_fracture")),
//_eps_real(declareProperty<RealTensorValue>("strain_real_inclusion")),
//_eps_imag(declareProperty<RealTensorValue>("strain_imag")),
//_tr_eps(declareProperty<Number>("trace_strain")),
//_sigma_real(declareProperty<RealTensorValue>("stress_real")),
//_sigma_imag(declareProperty<RealTensorValue>("stress_imag")),
// _omega(declareProperty<Real>("omega_property")),
// _mu(declareProperty<Real>("mu_property")),
// _lambda(declareProperty<Real>("lambda_property")),
// _alpha(declareProperty<Real>("alpha_property")),
// _kappa(declareProperty<Real>("kappa_property")),
// _eta(declareProperty<Real>("eta_property")),
// _porosity(declareProperty<Real>("porosity_property")),
// _kf(declareProperty<Real>("kf_property")),
// _ks(declareProperty<Real>("ks_property")),
// _inv_m(declareProperty<Real>("inverse_of_m")),
// _diffusion(declareProperty<Real>("diffusion_property")),
//_grad_disp_x(coupledGradient("disp_x")),
//_grad_disp_y(coupledGradient("disp_y")),
//_grad_disp_z(coupledGradient("disp_z")),
_pressure(  isParamValid("pressure") ? coupledValue("pressure") : _zero),
_imagUnit(0.0,1.0),
_pi( std::acos(-1.0) ),
_hasMeshModifier(parameters.isParamValid("inclusion_meshmodifier"))
{
	int dim=_mesh.dimension();
	_identity=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

	for (int i=0; i<dim; ++i)
	{
		_identity(i,i)=1.0;
	}

	_diffusion_block = _kappa_block   /_eta_block  ;
	_om_block        = _porosity_block   /_kf_block    + (_alpha_block   -_porosity_block   )/_ks_block;

	if (_hasMeshModifier)
	{
		_meshModifierName=getParam<std::string>("inclusion_meshmodifier");
		MeshModifier const & _myMeshModifier ( _app.getMeshModifier( _meshModifierName.c_str() )  );
		InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(_myMeshModifier) );

		int input_fn=inclusionsMeshModifier.get_fn();
		int total_fn=inclusionsMeshModifier.get_total_fn();

		if (_mu_fracture.size()<input_fn && _mu_fracture.size()!=1)
		{
			mooseError("_mu_fracture.size()");
		}
		else
		{
			_mu_fracture.resize(total_fn);	
		}
		if (_lambda_fracture.size()<input_fn && _lambda_fracture.size()!=1)
		{
			mooseError("_lambda_fracture.size()");
		}
		else
		{
			_lambda_fracture.resize(total_fn);	
		}
		if (_alpha_fracture.size()<input_fn && _alpha_fracture.size()!=1)
		{
			mooseError("_alpha_fracture.size()");
		}
		else
		{
			_alpha_fracture.resize(total_fn);	
		}
		if (_kappa_fracture.size()<input_fn && _kappa_fracture.size()!=1)
		{
			mooseError("_kappa_fracture.size()");
		}
		else
		{
			_kappa_fracture.resize(total_fn);	
		}
		if (_eta_fracture.size()<input_fn && _eta_fracture.size()!=1)
		{
			mooseError("_eta_fracture.size()");
		}
		else
		{
			_eta_fracture.resize(total_fn);	
		}
		if (_porosity_fracture.size()<input_fn && _porosity_fracture.size()!=1)
		{
			mooseError("_porosity_fracture.size()");
		}
		else
		{
			_porosity_fracture.resize(total_fn);	
		}
		if (_kf_fracture.size()<input_fn && _kf_fracture.size()!=1)
		{
			mooseError("_kf_fracture.size()");
		}
		else
		{
			_kf_fracture.resize(total_fn);	
		}
		if (_ks_fracture.size()<input_fn && _ks_fracture.size()!=1)
		{
			mooseError("_ks_fracture.size()");
		}
		else
		{
			_ks_fracture.resize(total_fn);	
		}

		_diffusion_fracture.resize(total_fn);
		_om_fracture.resize(total_fn);

		for (int i=0; i<total_fn; ++i)
		{
			_diffusion_fracture.at(i) = _kappa_fracture.at(i)/_eta_fracture.at(i);
			_om_fracture       .at(i) = _porosity_fracture.at(i)/_kf_fracture.at(i) + (_alpha_fracture.at(i)-_porosity_fracture.at(i))/_ks_fracture.at(i);
		}


	}



}

void
	FreqPoroelasticInclusion::computeQpProperties()
{

	// Real x_coord = _q_point[_qp](0);
	// Real y_coord = _q_point[_qp](1);


	// _omega[_qp] = 2.0*_pi*std::pow(10,_t);
	// // We set parameters to block parameters
	// _mu[_qp] = _mu_block;
	// _lambda[_qp] = _lambda_block;
	// _alpha[_qp] = _alpha_block;
	// _kappa[_qp] = _kappa_block;
	// _eta[_qp] = _eta_block;
	// _porosity[_qp] = _porosity_block;
	// _kf[_qp] = _kf_block;
	// _ks[_qp] = _ks_block;
    
    
    
    
    
	// for (int i = 0; i < _fn; i++)
	// {
	// 	Real temp =_a[i]*x_coord+_b[i]*y_coord + _c[i];
	// 	if ( std::fabs(temp) <=  _ft[i]/2.0 )
	// 	{
	// 		Real tempo=_ao[i]*x_coord+_bo[i]*y_coord + _co[i];
	// 		if (std::fabs(tempo) <=  _fl[i]/2.0 )
	// 		{
	// 			// std::cout<<"dentro "<<std::flush;
	// 			_mu[_qp] =  _mu_fracture;
	// 			_lambda[_qp] = _lambda_fracture;
	// 			_alpha[_qp] = _alpha_fracture;
	// 			_kappa[_qp] = _kappa_fracture;
	// 			//_kappa[_qp] = _ft[i]*_ft[i]/9.6e4;
	// 			//std::cout << _kappa[_qp] << std::endl;
	// 			_eta[_qp] = _eta_fracture;
	// 			_porosity[_qp] = _porosity_fracture;
	// 			_kf[_qp] = _kf_fracture;
	// 			_ks[_qp] = _ks_fracture;
	// 			break;

	// 		}
	// 	}
	// }
	// // The sign in the paper of Quintal et al. 2011 is wrong
	// _inv_m[_qp] = _porosity[_qp]/_kf[_qp] + (_alpha[_qp]-_porosity[_qp])/_ks[_qp];
	// _diffusion[_qp] = _kappa[_qp]/_eta[_qp]/_omega[_qp];
    
	// // Kinematics
	// Number temp00=_grad_disp_x[_qp](0);
	// _U_real(0,0) = temp00.real();
	// _U_imag(0,0) = temp00.imag();
	// Number temp01=_grad_disp_x[_qp](1);
	// _U_real(0,1) = temp01.real();
	// _U_imag(0,1) = temp01.imag();
	// // Third component
	// _U_real(0,2) = 0.0;
	// _U_imag(0,2) = 0.0;
    
	// Number temp10=_grad_disp_y[_qp](0);
	// _U_real(1,0) = temp10.real();
	// _U_imag(1,0) = temp10.imag();
	// Number temp11=_grad_disp_y[_qp](1);
	// _U_real(1,1) = temp11.real();
	// _U_imag(1,1) = temp11.imag();
	// // Third component
	// _U_real(1,2) = 0.0;
	// _U_imag(1,2) = 0.0;
    
	// _U_real(2,0) = 0.0;
	// _U_imag(2,0) = 0.0;
	// _U_real(2,1) = 0.0;
	// _U_imag(2,1) = 0.0;
	// _U_real(2,2) = 0.0;
	// _U_imag(2,2) = 0.0;
    
	// _eps_real[_qp] = 0.5*(_U_real+_U_real.transpose());
	// _eps_imag[_qp] = 0.5*(_U_imag+_U_imag.transpose());
    
	// _tr_eps[_qp] = _eps_real[_qp].tr()+_imagUnit*_eps_imag[_qp].tr();
    
	// // Stress
	// _sigma_real[_qp] = 2.0*_mu[_qp]*_eps_real[_qp]+_lambda[_qp]*_eps_real[_qp].tr()*_identity;
	// _sigma_imag[_qp] = 2.0*_mu[_qp]*_eps_imag[_qp]+_lambda[_qp]*_eps_imag[_qp].tr()*_identity;
    
    
}

void FreqPoroelasticInclusion::computeQpProperties(	Point const & point,
													Real & mu,
													Real & lambda,
													Real & alpha,
													Real & diffusion,
													Real & oneOverM) const
{
	mu=_mu_block;
	lambda=_lambda_block;
	alpha=_alpha_block;
	diffusion=_diffusion_block;
	oneOverM=_om_block;

	if (_hasMeshModifier)
	{
		MeshModifier const & _myMeshModifier ( _app.getMeshModifier( _meshModifierName.c_str() )  );
		InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(_myMeshModifier) );
		if ( inclusionsMeshModifier.isInside(point) )
		{
			mu       = _mu_fracture.at(0);
			lambda   = _lambda_fracture.at(0);
			alpha    = _alpha_fracture.at(0);
			diffusion= _diffusion_fracture.at(0);
			oneOverM = _om_fracture.at(0);
		}
	}
}

