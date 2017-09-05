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

#include "AttenuationDispersion.h"

template<>
InputParameters validParams<AttenuationDispersion>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("pp_name_real_stress", "Get the name of the postprocessor that computed the real part of the stress");
  params.addRequiredParam<PostprocessorName>("pp_name_imag_stress", "Get the name of the postprocessor that computed the imaginary part of the stress");
  params.addRequiredParam<PostprocessorName>("pp_name_real_strain", "Get the name of the postprocessor that computed the real part of the strain");
  params.addRequiredParam<PostprocessorName>("pp_name_imag_strain", "Get the name of the postprocessor that computed the imaginary part of the strain");
  params.addRequiredParam<unsigned int>("doatt", "flag for computing the attenuation");
  params.addRequiredParam<unsigned int>("doshear", "flag for computing the shear-wave dispersion");
  return params;
}

AttenuationDispersion::AttenuationDispersion(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _pp_val_stress_real(getPostprocessorValue("pp_name_real_stress")),
    _pp_val_stress_imag(getPostprocessorValue("pp_name_imag_stress")),
    _pp_val_strain_real(getPostprocessorValue("pp_name_real_strain")),
    _pp_val_strain_imag(getPostprocessorValue("pp_name_imag_strain")),
    _doatt(getParam<unsigned int>("doatt")),
    _doshear(getParam<unsigned int>("doshear"))
{}

Number
AttenuationDispersion::getValue()
{
  if (_doatt == 1) {
    Number nom = _pp_val_strain_real*_pp_val_stress_imag - _pp_val_stress_real*_pp_val_strain_imag;
    nom = nom/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);
    Number denom = _pp_val_stress_real*_pp_val_strain_real + _pp_val_stress_imag*_pp_val_strain_imag;
    denom = denom/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);
    return nom/denom;
  } else {
    if (_doshear == 1) {
      Number temp = _pp_val_stress_real*_pp_val_strain_real + _pp_val_stress_imag*_pp_val_strain_imag;
      return temp/(2.0*(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag));
    } else {
      Number temp = _pp_val_stress_real*_pp_val_strain_real + _pp_val_stress_imag*_pp_val_strain_imag;
      return temp/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);
    }
  }
}

