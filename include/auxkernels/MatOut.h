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

#ifndef MATOUT_H
#define MATOUT_H

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class MatOut;

template <>
InputParameters validParams<MatOut>();

class MatOut : public AuxKernel
{
public:
  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  MatOut(const InputParameters & parameters);

protected:
  virtual Number computeValue() override;

private:
  /// Reference to the gradient of the coupled variable
  const MaterialProperty<Real> & _eta;

};

#endif // MATOUT_H
