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

#ifndef FREQDIFFUSIONCUBE_H
#define FREQDIFFUSIONCUBE_H

#include "Material.h"

//Forward Declarations
class FreqDiffusionCube;

template<>
InputParameters validParams<FreqDiffusionCube>();

/**
 * Example material class that defines a few properties.
 */
class FreqDiffusionCube : public Material
{
public:
  FreqDiffusionCube(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  int _dim;
    Real _x_min,_x_max,_y_min,_y_max;//,_z_min,_zmin;
  MaterialProperty<Real> & _diffusion;
  MaterialProperty<Real> & _inv_m;
  // observe the next two properties are not needed in the kernels but in the postprocessors
  MaterialProperty<Real> &_omega_property;

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
};

#endif //FREQPOROELASTICSPHERE_H
