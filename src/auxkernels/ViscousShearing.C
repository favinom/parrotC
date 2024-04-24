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
#include "ViscousShearing.h"

template <>
InputParameters
validParams<ViscousShearing>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("disp_x", "real part of the first  coupled component");
  params.addRequiredCoupledVar("disp_y", "real part of the second coupled component");
  return params;
}

ViscousShearing::ViscousShearing(const InputParameters & parameters)
  : AuxKernel(parameters),
    _gradient_disp_x(coupledGradient("disp_x")),
    _gradient_disp_y(coupledGradient("disp_y")),
    _eta(getMaterialProperty<Real>("eta_property")),
    _omega(getMaterialProperty<Real>("omega_property")),
    _imagUnit(0.0,1.0)
{
}

Number
ViscousShearing::computeValue()
{
  return _eta[_qp]*_imagUnit*_omega[_qp]*(_gradient_disp_x[_qp](1)+_gradient_disp_y[_qp](0));
}
