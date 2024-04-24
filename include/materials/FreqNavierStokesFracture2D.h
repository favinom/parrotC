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

#ifndef FreqNavierStokesFracture2D_H
#define FreqNavierStokesFracture2D_H

#include "Material.h"

//Forward Declarations
class FreqNavierStokesFracture2D;

template<>
InputParameters validParams<FreqNavierStokesFracture2D>();

/**
 * Example material class that defines a few properties.
 */
class FreqNavierStokesFracture2D : public Material
{
public:
  FreqNavierStokesFracture2D(const InputParameters & parameters);
  ~FreqNavierStokesFracture2D()
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
  RealTensorValue _identity;
  std::string _f_center_x_string;
  std::string _f_center_z_string;
  std::string _f_thick_string;
  std::string _f_length_string;
  std::string _f_dip_string;
  Number const _imagUnit;
};

#endif //FreqNavierStokesFracture2D_H
