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

#ifndef VISCOUSSHEARING_H
#define VISCOUSSHEARING_H

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class ViscousShearing;

template <>
InputParameters validParams<ViscousShearing>();

class ViscousShearing : public AuxKernel
{
public:
  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  ViscousShearing(const InputParameters & parameters);

protected:
  virtual Number computeValue() override;

private:
  /// Reference to the gradient of the coupled variable
  const VariableNumberGradient & _gradient_disp_x;
  const VariableNumberGradient & _gradient_disp_y;
  const MaterialProperty<Real> & _eta;
  const MaterialProperty<Real> & _omega;
  Number const _imagUnit;

};

#endif // VISCOUSSHEARING_H
