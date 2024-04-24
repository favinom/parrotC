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

#ifndef RESIDUALFORCINGNODALKERNEL_H
#define RESIDUALFORCINGNODALKERNEL_H

#include "NodalKernel.h"

//Forward Declarations
class ResidualForcingNodalKernel;

template<>
InputParameters validParams<ResidualForcingNodalKernel>();

/**
 * Represents the rate in a simple ODE of du/dt = f
 */
class ResidualForcingNodalKernel : public NodalKernel
{
public:
  /**
   * Constructor grabs the Function
   */
  ResidualForcingNodalKernel(const InputParameters & parameters);

protected:
  virtual Number computeQpResidual() override;

  VariableValue const & _residual;
};

#endif
