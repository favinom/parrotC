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

#ifndef ATTENUATIONDISPERSION_H
#define ATTENUATIONDISPERSION_H

#include "GeneralPostprocessor.h"

class AttenuationDispersion;

template<>
InputParameters validParams<AttenuationDispersion>();

class AttenuationDispersion : public GeneralPostprocessor
{
public:
  AttenuationDispersion(const InputParameters & parameters);
  virtual void initialize() {}
  virtual void execute() {}

protected:
  virtual Number getValue();

  const PostprocessorValue & _pp_val_stress_real;
  const PostprocessorValue & _pp_val_stress_imag;
  const PostprocessorValue & _pp_val_strain_real;
  const PostprocessorValue & _pp_val_strain_imag;
  const unsigned int _doatt;
  const unsigned int _doshear;
};

#endif
