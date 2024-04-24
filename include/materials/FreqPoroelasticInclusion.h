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

#ifndef FreqPoroelasticInclusion_H
#define FreqPoroelasticInclusion_H

#include "Material.h"

//Forward Declarations
class FreqPoroelasticInclusion;

template<>
InputParameters validParams<FreqPoroelasticInclusion>();

/**
 * Example material class that defines a few properties.
 */
class FreqPoroelasticInclusion : public Material
{
public:
  FreqPoroelasticInclusion(const InputParameters & parameters);
  ~FreqPoroelasticInclusion() { _console<<"~FreqPoroelasticInclusion()\n"; };

  void computeQpProperties( Point const &, Real &,Real &,Real &,Real &,Real &) const;
			
protected:
  virtual void computeQpProperties();
          
//  virtual void initQpStatefulProperties();

private:
    
  Real const _mu_block;
  Real const _lambda_block;
  Real const _alpha_block;
  Real const _kappa_block;
  Real const _eta_block;
  Real const _porosity_block;
  Real const _kf_block;
  Real const _ks_block;
  Real _diffusion_block;
  Real _om_block;
   
  std::vector<Real> _mu_fracture;
  std::vector<Real> _lambda_fracture;
  std::vector<Real> _alpha_fracture;
  std::vector<Real> _kappa_fracture;
  std::vector<Real> _eta_fracture;
  std::vector<Real> _porosity_fracture;
  std::vector<Real> _kf_fracture;
  std::vector<Real> _ks_fracture;
  std::vector<Real> _diffusion_fracture;
  std::vector<Real> _om_fracture;

//  RealTensorValue _U_real;
//  RealTensorValue _U_imag;
//  MaterialProperty<RealTensorValue> & _eps_real;
//  MaterialProperty<RealTensorValue> & _eps_imag;
//  MaterialProperty<Number> &_tr_eps;
//  MaterialProperty<RealTensorValue> & _sigma_real;
//  MaterialProperty<RealTensorValue> & _sigma_imag;
//  MaterialProperty<Real> & _omega;
//  MaterialProperty<Real> & _mu;
//  MaterialProperty<Real> & _lambda;
//  MaterialProperty<Real> & _alpha;
//  MaterialProperty<Real> & _kappa;
//  MaterialProperty<Real> & _eta;
//  MaterialProperty<Real> & _porosity;
//  MaterialProperty<Real> & _kf;
//  MaterialProperty<Real> & _ks;
//  MaterialProperty<Real> & _inv_m;
//  MaterialProperty<Real> & _diffusion;
//    VariableGradient const & _grad_disp_x;
//    VariableGradient const & _grad_disp_y;
//    VariableGradient const & _grad_disp_z;
    VariableValue    const & _pressure;
    
    Number const _imagUnit;
    Real const _pi;
    RealTensorValue _identity;

    std::string _meshModifierName;
    bool const _hasMeshModifier;

};

#endif //FreqPoroelasticInclusion_H
