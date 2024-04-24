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

#ifndef FreqPoroelasticFracture2DInvasion_H
#define FreqPoroelasticFracture2DInvasion_H

#include "Material.h"

//Forward Declarations
class FreqPoroelasticFracture2DInvasion;

template<>
InputParameters validParams<FreqPoroelasticFracture2DInvasion>();

/**
 * Example material class that defines a few properties.
 */
class FreqPoroelasticFracture2DInvasion : public Material
{
public:
  FreqPoroelasticFracture2DInvasion(const InputParameters & parameters);
  ~FreqPoroelasticFracture2DInvasion()
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
        delete [] _invasion_data;        

        std::cout<<"Called destructor\n";
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
  
  std::string _filename_invasion;
  int _nx_invasion; 
  Real _dx_invasion; 
  short * _invasion_data;

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
  Real _mu_block;
  Real _lambda_block;
  Real _alpha_block;
  Real _kappa_block;
  Real _eta_block;
  Real _porosity_block;
  Real _kf_block;
  Real _ks_block;
  Real _mu_fracture;
  Real _lambda_fracture;
  Real _alpha_fracture;
  Real _kappa_fracture;
  Real _eta_fracture_invaded;
  Real _eta_fracture_brine;
  Real _porosity_fracture;
  Real _kf_fracture_invaded;
  Real _kf_fracture_brine;
  Real _ks_fracture;
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

#endif //FreqPoroelasticFracture2DInvasion_H
