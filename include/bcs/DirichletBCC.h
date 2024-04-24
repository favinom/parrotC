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

#ifndef DIRICHLETBCC_H
#define DIRICHLETBCC_H

#include "NodalBC.h"

class DirichletBCC;

template <>
InputParameters validParams<DirichletBCC>();

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 */
class DirichletBCC : public NodalBC
{
public:
  DirichletBCC(const InputParameters & parameters);

protected:
  virtual Number computeQpResidual() override;

  /// The value for this BC
  const Number _value;
};

#endif /* DIRICHLETBC_H */
