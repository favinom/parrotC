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

#include "FreqDiffusionCube.h"

#include "MooseMesh.h"

template<>
InputParameters validParams<FreqDiffusionCube>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("x_min", "x_min");
  params.addRequiredParam<Real>("x_max", "x_max");
  params.addRequiredParam<Real>("y_min", "y_min");
  params.addRequiredParam<Real>("y_max", "y_max");
//  params.addRequiredParam<Real>("z_min", "z_min");
//  params.addRequiredParam<Real>("z_max", "z_max");
  return params;
}

FreqDiffusionCube::FreqDiffusionCube(const InputParameters & parameters) :
    Material(parameters),
    _dim(_mesh.dimension()),
    _x_min(getParam<Real>("x_min")),
    _x_max(getParam<Real>("x_max")),
    _y_min(getParam<Real>("y_min")),
    _y_max(getParam<Real>("y_max")),
//    _z_center(_mesh.dimension() == 3 ? getParam<Real>("z_center") : 0.0),
    _diffusion(declareProperty<Real>("diffusion_property")),
    _inv_m(declareProperty<Real>("inverse_of_m")),
    _omega_property(declareProperty<Real>("omega_property")),
    _imagUnit(0.0,1.0)
{
    
    _pi = std::acos(-1.0);
    
    //Real sf=1e0;
    
    
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
FreqDiffusionCube::computeQpProperties()
{
    _omega = 2.0*_pi*std::pow(10.0,_t);
    _omega_property[_qp]=_omega;
    
    Real x_coord = _q_point[_qp](0);
    Real y_coord = _q_point[_qp](1);
    Real z_coord = _q_point[_qp](2);

    int c1=0,c2=0;
    
    if (1)// _x_min<x_coord &&  x_coord<_x_max && _y_min<y_coord && y_coord<_y_max )
    // We are inside
    {
        _inv_m[_qp] = _porosity_fracture/_kf_fracture + (_alpha_fracture-_porosity_fracture)/_ks_fracture;
        _diffusion[_qp] = _kappa_fracture/_eta_fracture/_omega;
        ++c1;
    }
    else
    // We are outside
    {
        _inv_m[_qp] = _porosity_block/_kf_block + (_alpha_block-_porosity_block)/_ks_block;
        _diffusion[_qp] = _kappa_block/_eta_block/_omega;
        ++c2;
    }
    
    if (c1!=0 && c2!=0)
    {
        exit(1);
    }

}
