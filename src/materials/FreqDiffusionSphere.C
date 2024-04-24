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

#include "FreqDiffusionSphere.h"

#include "MooseMesh.h"

template<>
InputParameters validParams<FreqDiffusionSphere>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("x_center", "x coordinate of center");
  params.addRequiredParam<Real>("y_center", "y coordinate of center");
  params.addParam<Real>("z_center", "z coordinate of center");
  params.addRequiredParam<Real>("radius", "radius");
  //params.addRequiredCoupledVar("disp_x", "real part of the first coupled component");
  //params.addRequiredCoupledVar("disp_y", "real part of the second coupled component");
  //params.addCoupledVar("disp_z", "real part of the third coupled component");

  return params;
}

FreqDiffusionSphere::FreqDiffusionSphere(const InputParameters & parameters) :
    Material(parameters),
    _dim(_mesh.dimension()),
    _x_center(getParam<Real>("x_center")),
    _y_center(getParam<Real>("y_center")),
    _z_center(_mesh.dimension() == 3 ? getParam<Real>("z_center") : 0.0),
    _radius(getParam<Real>("radius")),
    //_mu(declareProperty<Real>("mu_property")),
    //_lambda(declareProperty<Real>("lambda_property")),
    //_alpha(declareProperty<Real>("alpha_property")),
    //_sigma_real(declareProperty<RealTensorValue>("stress_real")),
    //_sigma_imag(declareProperty<RealTensorValue>("stress_imag")),
    //_eps_real(declareProperty<RealTensorValue>("strain_real")),
    //_eps_imag(declareProperty<RealTensorValue>("strain_imag")),
    _diffusion(declareProperty<Real>("diffusion_property")),
    _inv_m(declareProperty<Real>("inverse_of_m")),
    //_tr_eps(declareProperty<Number>("trace_strain")),
    _omega_property(declareProperty<Real>("omega_property")),
    //_grad_disp_x(coupledGradient("disp_x")),
    //_grad_disp_y(coupledGradient("disp_y")),
    //_grad_disp_z(_mesh.dimension() == 3 ? coupledGradient("disp_z") : _grad_zero),
    _imagUnit(0.0,1.0)
{
    
    _pi = std::acos(-1.0);
    
    if (_dim == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    
    if (_dim == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    Real sf=1e0;
    
    
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

}

void
FreqDiffusionSphere::computeQpProperties()
{
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
        //_mu[_qp] = _mu_fracture;
        //_lambda[_qp] = _lambda_fracture;
        //_alpha[_qp] = _alpha_fracture;
        _inv_m[_qp] = _porosity_fracture/_kf_fracture + (_alpha_fracture-_porosity_fracture)/_ks_fracture;
        _diffusion[_qp] = _kappa_fracture/_eta_fracture/_omega;
    }
    else
    // We are outside
    {
        // JUERG can you please ckeck?
        //_mu[_qp] = _mu_block;
        //_lambda[_qp] = _lambda_block;
        //_alpha[_qp] = _alpha_block;
        _inv_m[_qp] = _porosity_block/_kf_block + (_alpha_block-_porosity_block)/_ks_block;
        _diffusion[_qp] = _kappa_block/_eta_block/_omega;
    }

}
