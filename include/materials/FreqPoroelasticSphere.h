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

#ifndef FREQPOROELASTICSPHERE_H
#define FREQPOROELASTICSPHERE_H

#include "Material.h"

//Forward Declarations
class FreqPoroelasticSphere;

template<>
InputParameters validParams<FreqPoroelasticSphere>();

/**
 * Example material class that defines a few properties.
 */
class FreqPoroelasticSphere : public Material
{
public:
  FreqPoroelasticSphere(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  int _dim;
  Real _x_center,_y_center,_z_center,_radius;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _lambda;
  MaterialProperty<Real> & _alpha;
  MaterialProperty<RealTensorValue> & _sigma_real;
  MaterialProperty<RealTensorValue> & _sigma_imag;
  // observe the next two properties are not needed in the kernels but in the postprocessors
  MaterialProperty<RealTensorValue> & _eps_real;
  MaterialProperty<RealTensorValue> & _eps_imag;
  MaterialProperty<Real> & _diffusion;
  MaterialProperty<Real> & _inv_m;
  MaterialProperty<Number> &_tr_eps;
  // observe the next two properties are not needed in the kernels but in the postprocessors
  MaterialProperty<Real> &_omega_property;
  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  VariableGradient const & _grad_disp_z;

  Real _pi;
  Number const _imagUnit;
  RealTensorValue _identity;

  Real _mu_block;
  Real _lambda_block;
  Real _alpha_block;
  Real _kappa_block;
  Real _eta_block;
  Real _porosity_block;
  Real _kf_block;
  Real _ks_block;
  Real _kd_block;
    
  Real _mu_fracture;
  Real _lambda_fracture;
  Real _alpha_fracture;
  Real _kappa_fracture;
  Real _eta_fracture;
  Real _porosity_fracture;
  Real _kf_fracture;
  Real _ks_fracture;
  Real _kd_fracture;
  
  Real _omega;
    
  RealTensorValue _U_real;
  RealTensorValue _U_imag;

};

#endif //FREQPOROELASTICSPHERE_H
