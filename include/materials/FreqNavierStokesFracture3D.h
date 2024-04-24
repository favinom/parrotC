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

#ifndef FreqNavierStokesFracture3D_H
#define FreqNavierStokesFracture3D_H

#include "Material.h"

//Forward Declarations
class FreqNavierStokesFracture3D;

template<>
InputParameters validParams<FreqNavierStokesFracture3D>();

/**
 * Example material class that defines a few properties.
 */
class FreqNavierStokesFracture3D : public Material
{
public:
  FreqNavierStokesFracture3D(const InputParameters & parameters);
  ~FreqNavierStokesFracture3D();

protected:
  virtual void computeQpProperties();
//  virtual void initQpStatefulProperties();

  void ComputeNormalsFromAngles(RealVectorValue const & angles,
                                RealVectorValue & n1,
                                RealVectorValue & n2,
                                RealVectorValue & n3);

  bool is_inside(RealVectorValue const & point);

private:
  int _fn; 
  std::string _fx_string; 
  std::string _fy_string; 
  std::string _fz_string; 
  std::string _fa1_string;
  std::string _fa2_string;
  std::string _fa3_string;
  std::string _fd1_string;
  std::string _fd2_string;
  std::string _fd3_string;

  RealVectorValue * _center;
  RealVectorValue * _rotation;
  RealVectorValue * _dimension;

  RealVectorValue ** _n;

  RealVectorValue *_d;

  Real _pi;
    
  RealTensorValue _U_real;
  RealTensorValue _U_imag;
  MaterialProperty<RealTensorValue> & _eps_real;
  MaterialProperty<RealTensorValue> & _eps_imag;
  MaterialProperty<Number> &_tr_eps;
  MaterialProperty<RealTensorValue> & _sigma_real;
  MaterialProperty<RealTensorValue> & _sigma_imag;
  Real _mu_block;
  Real _K_block;
  Real _eta_block;
  Real _mu_fracture;
  Real _K_fracture;
  Real _eta_fracture;
  MaterialProperty<Real> & _omega;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _K;
  MaterialProperty<Real> & _eta;
  MaterialProperty<Number> & _alpha;
  MaterialProperty<Number> & _beta;
  MaterialProperty<Real> & _dummy_alpha;
  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  VariableGradient const & _grad_disp_z;
  RealTensorValue _identity;
  Number const _imagUnit;
};

#endif //FreqNavierStokesFracture3D_H
