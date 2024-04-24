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

#include "MeanMaterialProperty.h"
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<MeanMaterialProperty>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<unsigned int>("i", "first component");
  params.addRequiredParam<unsigned int>("j", "second component");
  params.addRequiredParam<unsigned int>("dostress", "flag for stress");
  params.addRequiredParam<unsigned int>("doreal", "flag for real part");
  params.addRequiredParam<unsigned int>("dofreq", "flag for printing the frequency");
  params.addCoupledVar("pressure", "real or imaginary pressure");
  return params;
}

MeanMaterialProperty::MeanMaterialProperty(const InputParameters & parameters) :
    ElementIntegralPostprocessor(parameters),
    _eps_real(getMaterialProperty<RealTensorValue>("strain_real")),
    _eps_imag(getMaterialProperty<RealTensorValue>("strain_imag")),
    _sigma_real(getMaterialProperty<RealTensorValue>("stress_real")),
    _sigma_imag(getMaterialProperty<RealTensorValue>("stress_imag")),
    // _has_alpha(hasMaterialProperty<Real>("alpha_property")),
    // _alpha(_has_alpha ? getMaterialProperty<Real>("alpha_property"):_zero),
    // _alpha(_zero),
    _alpha(getMaterialProperty<Real>("alpha_property")),
    _omega(getMaterialProperty<Real>("omega_property")),
    _p(parameters.isParamValid("pressure") ? coupledValue("pressure"):_zero),
    _i(getParam<unsigned int>("i")),
    _j(getParam<unsigned int>("j")),
    _dostress(getParam<unsigned int>("dostress")),
    _doreal(getParam<unsigned int>("doreal")),
    _dofreq(getParam<unsigned int>("dofreq"))
{}

Number
MeanMaterialProperty::computeQpIntegral()
{
  if (_dofreq == 1)
  {
    //std::cout << _qrule->n_points() << std::endl;
    return _omega[_qp];
  }
  else
  {
      if (_dostress == 1)
      {
          if (_doreal == 1)
          {
              if (_i == _j)
              {
                  Number const & local_p = _p[_qp];
                  return _sigma_real[_qp](_i,_j) - _alpha[_qp]*local_p.real();
                  // return _sigma_real[_qp](_i,_j);
              }
              else
              {
                  return _sigma_real[_qp](_i,_j);
              }
          }
          else
          {
              if (_i == _j)
              {
                  Number const & local_p = _p[_qp];
                  return _sigma_imag[_qp](_i,_j) - _alpha[_qp]*local_p.imag();
                  // return _sigma_imag[_qp](_i,_j);
              }
              else
              {
                  return _sigma_imag[_qp](_i,_j);
              }
          }
          return 0.0;
      }
      else
      {
          if (_doreal == 1)
          {
              return _eps_real[_qp](_i,_j);
          }
          else
          {
              return _eps_imag[_qp](_i,_j);
          }
      }
  }
}

