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

#include "DirichletBCC.h"

template <>
InputParameters
validParams<DirichletBCC>()
{
  InputParameters p = validParams<NodalBC>();
  p.addParam<Real>("value_r", 0.0, "Value of the BC");
  p.addParam<Real>("value_i", 0.0, "Value of the BC");
//  p.declareControllable("value");
  return p;
}

DirichletBCC::DirichletBCC(const InputParameters & parameters) :
NodalBC(parameters),
_value(getParam<Real>("value_r"), getParam<Real>("value_i"))
{
}

Number
DirichletBCC::computeQpResidual()
{
	  return _u[_qp] - _value;
}
