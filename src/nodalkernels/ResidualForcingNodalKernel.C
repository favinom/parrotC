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

#include "ResidualForcingNodalKernel.h"


template<>
InputParameters validParams<ResidualForcingNodalKernel>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addRequiredCoupledVar("residual", "residual");
  return params;
}

ResidualForcingNodalKernel::ResidualForcingNodalKernel(const InputParameters & parameters) :
    NodalKernel(parameters),
    _residual(coupledValue("residual"))
{
}

Number
ResidualForcingNodalKernel::computeQpResidual()
{
  return _residual[_qp];
}
