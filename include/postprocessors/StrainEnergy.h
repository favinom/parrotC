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

#ifndef STRAIN_ENERGY_H
#define STRAIN_ENERGY_H

#include "ElementIntegralPostprocessor.h"

class StrainEnergy;

template<>
InputParameters validParams<StrainEnergy>();

class StrainEnergy : public ElementIntegralPostprocessor
{
public:
  StrainEnergy(const InputParameters & parameters);

protected:
  virtual Number computeQpIntegral() override;

  const MaterialProperty<Number> & _tr_eps;
  const MaterialProperty<RealTensorValue> & _eps_real;
  const MaterialProperty<RealTensorValue> & _eps_imag;
  const MaterialProperty<RealTensorValue> & _sigma_real;
  const MaterialProperty<RealTensorValue> & _sigma_imag;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _inv_m;
  const VariableValue & _pres;

};

#endif
