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

#include "FreqPoroelasticSphere.h"

#include "MooseMesh.h"

template<>
InputParameters validParams<FreqPoroelasticSphere>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("x_center", "x coordinate of center");
  params.addRequiredParam<Real>("y_center", "y coordinate of center");
  params.addParam<Real>("z_center", "z coordinate of center");
  params.addRequiredParam<Real>("radius", "radius");
  params.addRequiredCoupledVar("disp_x", "real part of the first coupled component");
  params.addRequiredCoupledVar("disp_y", "real part of the second coupled component");
  params.addCoupledVar("disp_z", "real part of the third coupled component");

  return params;
}

FreqPoroelasticSphere::FreqPoroelasticSphere(const InputParameters & parameters) :
    Material(parameters),
    _dim(_mesh.dimension()),
    _x_center(getParam<Real>("x_center")),
    _y_center(getParam<Real>("y_center")),
    _z_center(_mesh.dimension() == 3 ? getParam<Real>("z_center") : 0.0),
    _radius(getParam<Real>("radius")),
    _mu(declareProperty<Real>("mu_property")),
    _lambda(declareProperty<Real>("lambda_property")),
    _alpha(declareProperty<Real>("alpha_property")),
    _sigma_real(declareProperty<RealTensorValue>("stress_real")),
    _sigma_imag(declareProperty<RealTensorValue>("stress_imag")),
    _eps_real(declareProperty<RealTensorValue>("strain_real")),
    _eps_imag(declareProperty<RealTensorValue>("strain_imag")),
    _diffusion(declareProperty<Real>("diffusion_property")),
    _inv_m(declareProperty<Real>("inverse_of_m")),
    _tr_eps(declareProperty<Number>("trace_strain")),
    _omega_property(declareProperty<Real>("omega_property")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y")),
    _grad_disp_z(_mesh.dimension() == 3 ? coupledGradient("disp_z") : _grad_zero),
    _imagUnit(0.0,1.0)
{
    
    _pi = std::acos(-1.0);
    
    if (_dim == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    
    if (_dim == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    
    // outside
    
    _mu_block = 3.0e9;         // ok
    _kd_block = 4.0e9;        // ok
    _lambda_block = _kd_block-2.0/3.0*_mu_block;
    _kappa_block = 1.0e-12;   // ok permeability
    _eta_block = 2.0e-5;      // ok
    _porosity_block = 0.25;   // ok also known as phi
    _kf_block = 0.04e9;       // ok
    _ks_block = 40.0e9;       // ok
    _alpha_block = 1.0-_kd_block/_ks_block;
    /// FROM here on there are the parameters of frac
    
//    Real E_local=_mu_block*(2.0*_mu_block+3.0*_lambda_block)/(_mu_block+_lambda_block);
//    Real diff_local=_kappa_block/_eta_block;
//    std::cout<<E_local<<std::endl<<diff_local<<std::endl;
//    int i;
//    std::cin>>i;
    
    // inside
    
  _mu_fracture = 3.0e9;       // ok
  _kd_fracture = 4.0e9;       // ok
  _lambda_fracture = _kd_fracture-2.0/3.0*_mu_fracture;
  _kappa_fracture = 1.0e-12;  // ok
  _eta_fracture = 0.003;      // ok
  _porosity_fracture = 0.25;  // ok
  _kf_fracture = 2.4e9;       // ok
  _ks_fracture = 40.0e9;      // ok
  _alpha_fracture = 1.0-_kd_fracture/_ks_fracture;
    
//    E_local=_mu_fracture*(2.0*_mu_fracture+3.0*_lambda_fracture)/(_mu_fracture+_lambda_fracture);
//    diff_local=_kappa_fracture/_eta_fracture;
//    std::cout<<E_local<<std::endl<<diff_local<<std::endl;
//    std::cin>>i;

}

void
FreqPoroelasticSphere::computeQpProperties()
{
    
    // BE CAREFUL...
    Real sf=1.0e9;
    
    _omega = 2.0*_pi*std::pow(10.0,_t);
    _omega_property[_qp]=_omega;
    
    Real x_coord = _q_point[_qp](0);
    Real y_coord = _q_point[_qp](1);
    Real z_coord = _q_point[_qp](2);

    Real dist_x=x_coord-_x_center;
    Real dist_y=y_coord-_y_center;
    Real dist_z=z_coord-_z_center;
    
    if (dist_x*dist_x+dist_y*dist_y+dist_z*dist_z<=_radius*_radius)
    // We are inside
    {
        // JUERG can you please ckeck?
        _mu[_qp] = _mu_fracture/sf;
        _lambda[_qp] = _lambda_fracture/sf;
        _alpha[_qp] = _alpha_fracture;
        _inv_m[_qp] = (_porosity_fracture/_kf_fracture + (_alpha_fracture-_porosity_fracture)/_ks_fracture)*sf;
        _diffusion[_qp] = _kappa_fracture/_eta_fracture/_omega*sf;
    }
    else
    // We are outside
    {
        // JUERG can you please ckeck?
        _mu[_qp] = _mu_block/sf;
        _lambda[_qp] = _lambda_block/sf;
        _alpha[_qp] = _alpha_block;
        _inv_m[_qp] = (_porosity_block/_kf_block + (_alpha_block-_porosity_block)/_ks_block)*sf;
        _diffusion[_qp] = _kappa_block/_eta_block/_omega*sf;
    }
    
    Number temp00=_grad_disp_x[_qp](0);
    _U_real(0,0) = temp00.real();
    _U_imag(0,0) = temp00.imag();
    Number temp01=_grad_disp_x[_qp](1);
    _U_real(0,1) = temp01.real();
    _U_imag(0,1) = temp01.imag();
    Number temp02=_grad_disp_x[_qp](2);
    _U_real(0,2) = temp02.real();
    _U_imag(0,2) = temp02.imag();
    
    Number temp10=_grad_disp_y[_qp](0);
    _U_real(1,0) = temp10.real();
    _U_imag(1,0) = temp10.imag();
    Number temp11=_grad_disp_y[_qp](1);
    _U_real(1,1) = temp11.real();
    _U_imag(1,1) = temp11.imag();
    Number temp12=_grad_disp_y[_qp](2);
    _U_real(1,2) = temp12.real();
    _U_imag(1,2) = temp12.imag();
    
    Number temp20=_grad_disp_z[_qp](0);
    _U_real(2,0) = temp20.real();
    _U_imag(2,0) = temp20.imag();
    Number temp21=_grad_disp_z[_qp](1);
    _U_real(2,1) = temp21.real();
    _U_imag(2,1) = temp21.imag();
    Number temp22=_grad_disp_z[_qp](2);
    _U_real(2,2) = temp22.real();
    _U_imag(2,2) = temp22.imag();
    
    _eps_real[_qp] = 0.5*(_U_real+_U_real.transpose());
    _eps_imag[_qp] = 0.5*(_U_imag+_U_imag.transpose());
    
    _tr_eps[_qp] = _eps_real[_qp].tr()+_imagUnit*_eps_imag[_qp].tr();
    _sigma_real[_qp] = 2.0*_mu[_qp]*_eps_real[_qp]+_lambda[_qp]*_eps_real[_qp].tr()*_identity;
    _sigma_imag[_qp] = 2.0*_mu[_qp]*_eps_imag[_qp]+_lambda[_qp]*_eps_imag[_qp].tr()*_identity;

}
