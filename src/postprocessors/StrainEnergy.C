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

#include "StrainEnergy.h"
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<StrainEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("pres", "pressure");
  return params;
}

StrainEnergy::StrainEnergy(const InputParameters & parameters) :
    ElementIntegralPostprocessor(parameters),
    _tr_eps(getMaterialProperty<Number>("trace_strain")),
    _eps_real(getMaterialProperty<RealTensorValue>("strain_real")),
    _eps_imag(getMaterialProperty<RealTensorValue>("strain_imag")),
    _sigma_real(getMaterialProperty<RealTensorValue>("stress_real")),
    _sigma_imag(getMaterialProperty<RealTensorValue>("stress_imag")),
    _alpha(getMaterialProperty<Real>("alpha_property")),
    _inv_m(getMaterialProperty<Real>("inverse_of_m")),
    _pres(coupledValue("pres"))
{}

Number
StrainEnergy::computeQpIntegral()
{
    Real pres_real=_pres[_qp].real();
    Real pres_imag=_pres[_qp].imag();
    Real tr_eps_real=_tr_eps[_qp].real();
    Real tr_eps_imag=_tr_eps[_qp].imag();

    Real z_real=pres_real*_inv_m[_qp]+_alpha[_qp]*tr_eps_real;
    Real z_imag=pres_imag*_inv_m[_qp]+_alpha[_qp]*tr_eps_imag;
    
    Real a=_sigma_real[_qp].contract(_eps_real[_qp]);
    Real b=_sigma_imag[_qp].contract(_eps_imag[_qp]);
    Real c=pres_real*z_real;
    Real d=pres_imag*z_imag;
    
    
    return (a+b+c+d)/4.0;
}

