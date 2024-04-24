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
#include "MatOut.h"

template <>
InputParameters
validParams<MatOut>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

MatOut::MatOut(const InputParameters & parameters)
  : AuxKernel(parameters),
    _eta(getMaterialProperty<Real>("eta_property"))
{
}

Number
MatOut::computeValue()
{
  return _eta[_qp];
}
