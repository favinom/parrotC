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

#ifndef MEANMATERIALPROPERTY_H
#define MEANMATERIALPROPERTY_H

#include "ElementIntegralPostprocessor.h"

class MeanMaterialProperty;

template<>
InputParameters validParams<MeanMaterialProperty>();

class MeanMaterialProperty : public ElementIntegralPostprocessor
{
public:
  MeanMaterialProperty(const InputParameters & parameters);

protected:
  virtual Number computeQpIntegral();

  const MaterialProperty<RealTensorValue> & _eps_real;
  const MaterialProperty<RealTensorValue> & _eps_imag;
  const MaterialProperty<RealTensorValue> & _sigma_real;
  const MaterialProperty<RealTensorValue> & _sigma_imag;
  // bool _has_alpha;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _omega;
  const VariableValue & _p;
  const unsigned int _i;
  const unsigned int _j;
  const unsigned int _dostress;
  const unsigned int _doreal;
  const unsigned int _dofreq;
};

#endif
