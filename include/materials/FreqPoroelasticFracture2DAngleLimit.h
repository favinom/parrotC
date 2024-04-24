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

#ifndef FreqPoroelasticFracture2DAngleLimit_H
#define FreqPoroelasticFracture2DAngleLimit_H

#include "Material.h"

//Forward Declarations
class FreqPoroelasticFracture2DAngleLimit;

template<>
InputParameters validParams<FreqPoroelasticFracture2DAngleLimit>();

/**
 * Example material class that defines a few properties.
 */
class FreqPoroelasticFracture2DAngleLimit : public Material
{
public:
  FreqPoroelasticFracture2DAngleLimit(const InputParameters & parameters);
  ~FreqPoroelasticFracture2DAngleLimit()
    {
        delete [] _fx;
        delete [] _fy;
        delete [] _ft;
        delete [] _fl;
        delete [] _fa;
        delete [] _a;
        delete [] _b;
        delete [] _c;
        delete [] _ao;
        delete [] _bo;
        delete [] _co;
        
        std::cout<<"Called desctructor\n";
    };

protected:
  virtual void computeQpProperties();
//  virtual void initQpStatefulProperties();

private:
  int _fn; 
  std::string _fx_string; 
  std::string _fy_string; 
  std::string _fl_string; 
  std::string _ft_string; 
  std::string _fa_string;

  Real * _fx;
  Real * _fy;
  Real * _fl;
  Real * _ft;
  Real * _fa;

    
    Real * _a;
    Real * _b;
    Real * _c;

    Real * _ao;
    Real * _bo;
    Real * _co;

    
    Real _pi;
    
  RealTensorValue _U_real;
  RealTensorValue _U_imag;
  MaterialProperty<RealTensorValue> & _eps_real;
  MaterialProperty<RealTensorValue> & _eps_imag;
  MaterialProperty<Number> &_tr_eps;
  MaterialProperty<RealTensorValue> & _sigma_real;
  MaterialProperty<RealTensorValue> & _sigma_imag;
  Real _mean_angle_set2;
  Real _halfspread_angle_set2;
  Real _mu_block;
  Real _lambda_block;
  Real _alpha_block;
  Real _kappa_block;
  Real _eta_block;
  Real _porosity_block;
  Real _kf_block;
  Real _ks_block;
  Real _mu_fracture1;
  Real _lambda_fracture1;
  Real _alpha_fracture1;
  Real _kappa_fracture1;
  Real _eta_fracture1;
  Real _porosity_fracture1;
  Real _kf_fracture1;
  Real _ks_fracture1;
  Real _mu_fracture2;
  Real _lambda_fracture2;
  Real _alpha_fracture2;
  Real _kappa_fracture2;
  Real _eta_fracture2;
  Real _porosity_fracture2;
  Real _kf_fracture2;
  Real _ks_fracture2;
  MaterialProperty<Real> & _omega;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _lambda;
  MaterialProperty<Real> & _alpha;
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _eta;
  MaterialProperty<Real> & _porosity;
  MaterialProperty<Real> & _kf;
  MaterialProperty<Real> & _ks;
  MaterialProperty<Real> & _inv_m;
  MaterialProperty<Real> & _diffusion;
  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  RealTensorValue _identity;
  std::string _f_center_x_string;
  std::string _f_center_z_string;
  std::string _f_thick_string;
  std::string _f_length_string;
  std::string _f_dip_string;
  Number const _imagUnit;
};

#endif //FreqPoroelasticFracture2DAngleLimit_H
