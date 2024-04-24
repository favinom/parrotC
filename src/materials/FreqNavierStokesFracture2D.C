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

#include "FreqNavierStokesFracture2D.h"

#include <sstream>
#include "MooseMesh.h"

template<>
InputParameters validParams<FreqNavierStokesFracture2D>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<int>("fn", "number of fractures");
  params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fl_string", "length of fractures");
  params.addRequiredParam<std::string>("ft_string", "thickness of fractures");
  params.addRequiredParam<std::string>("fa_string", "angle of fractures");
  params.addRequiredParam<Real>("mu_block", "shear modulus");
  params.addRequiredParam<Real>("K_block", "bulk modulus");
  params.addRequiredParam<Real>("eta_block", "shear viscocity");
  params.addRequiredParam<Real>("mu_fracture", "shear modulus");
  params.addRequiredParam<Real>("K_fracture", "bulk modulus");
  params.addRequiredParam<Real>("eta_fracture", "shear viscocity");
  params.addRequiredCoupledVar("disp_x", "complex first  coupled component");
  params.addRequiredCoupledVar("disp_y", "complex second coupled component");

  return params;
}

FreqNavierStokesFracture2D::FreqNavierStokesFracture2D(const InputParameters & parameters) :
    Material(parameters),
    _fn(getParam<int>("fn")), 
    _fx_string(getParam<std::string>("fx_string")), 
    _fy_string(getParam<std::string>("fy_string")), 
    _fl_string(getParam<std::string>("fl_string")), 
    _ft_string(getParam<std::string>("ft_string")), 
    _fa_string(getParam<std::string>("fa_string")), 
    _eps_real(declareProperty<RealTensorValue>("strain_real")),
    _eps_imag(declareProperty<RealTensorValue>("strain_imag")),
    _tr_eps(declareProperty<Number>("trace_strain")),
    _sigma_real(declareProperty<RealTensorValue>("stress_real")),
    _sigma_imag(declareProperty<RealTensorValue>("stress_imag")),
    _mu_block(getParam<Real>("mu_block")),
    _K_block(getParam<Real>("K_block")),
    _eta_block(getParam<Real>("eta_block")),
    _mu_fracture(getParam<Real>("mu_fracture")),
    _K_fracture(getParam<Real>("K_fracture")),
    _eta_fracture(getParam<Real>("eta_fracture")),
    _omega(declareProperty<Real>("omega_property")),
    _mu(declareProperty<Real>("mu_property")),
    _K(declareProperty<Real>("K_property")),
    _eta(declareProperty<Real>("eta_property")),
    _alpha(declareProperty<Number>("alpha_coefficient")),
    _beta(declareProperty<Number>("beta_coefficient")),
    _dummy_alpha(declareProperty<Real>("alpha_property")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y")),
    _imagUnit(0.0,1.0)
{
    
    std::cout<<"Constructor called\n";
    
    if (_mesh.dimension() != 2)
    {
        std::cout<<"You cannot use this material in dimension different from 2\n";
        exit(1);
    }

    _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    
    _pi = acos(-1.0);
    
    _fx = new Real [_fn];
    _fy = new Real [_fn];
    _fl = new Real [_fn];
    _ft = new Real [_fn];
    _fa = new Real [_fn];
    
    _a  = new Real [_fn];
    _b  = new Real [_fn];
    _c  = new Real [_fn];

    _ao = new Real [_fn];
    _bo = new Real [_fn];
    _co = new Real [_fn];
 
    std::istringstream fx_ss(_fx_string);
    std::istringstream fy_ss(_fy_string);
    std::istringstream fl_ss(_fl_string);
    std::istringstream ft_ss(_ft_string);
    std::istringstream fa_ss(_fa_string);
    
    std::string token;
    
    for (int i=0; i<_fn; ++i)
    {
        std::getline(fx_ss, token, ',');
        _fx[i]=atof(token.c_str());
        std::getline(fy_ss, token, ',');
        _fy[i]=atof(token.c_str());
        std::getline(fl_ss, token, ',');
        _fl[i]=atof(token.c_str());
        std::getline(ft_ss, token, ',');
        _ft[i]=atof(token.c_str());
        std::getline(fa_ss, token, ',');
        _fa[i]=atof(token.c_str());
        _fa[i]=(90.0-_fa[i])/180.0*pi;
        
        _a[i] =  std::sin(_fa[i]);
        _b[i] = -std::cos(_fa[i]);
        _c[i] = -std::sin(_fa[i])*_fx[i]+std::cos(_fa[i])*_fy[i];

        _ao[i] =  std::cos(_fa[i]);
        _bo[i] =  std::sin(_fa[i]);
        _co[i] = -std::cos(_fa[i])*_fx[i]-std::sin(_fa[i])*_fy[i];
        
        
    }
    
    
    // To be sure, we should put some check on the parsing of string streams
    // std::cout<<_fa[_fn-1]<<std::endl;
    
}

void
FreqNavierStokesFracture2D::computeQpProperties()
{
    _dummy_alpha[_qp] = 0.0;

    Real x_coord = _q_point[_qp](0);
    Real y_coord = _q_point[_qp](1);


    _omega[_qp] = 2.0*_pi*std::pow(10,_t);
    // We set parameters to block parameters
    _mu[_qp] = _mu_block;
    _K[_qp] = _K_block;
    _eta[_qp] = _eta_block;
    
    
    for (int i = 0; i < _fn; i++)
    {
        Real temp =_a[i]*x_coord+_b[i]*y_coord + _c[i];
        if ( std::fabs(temp) <=  _ft[i]/2.0 )
        {
            Real tempo=_ao[i]*x_coord+_bo[i]*y_coord + _co[i];
            if (std::fabs(tempo) <=  _fl[i]/2.0 )
            {
                //std::cout<<"dentro\n";
                _mu[_qp] =  _mu_fracture;
                _K[_qp] = _K_fracture;
                _eta[_qp] = _eta_fracture;
                break;

            }
        }
    }
    
    // Kinematics
    Number temp00=_grad_disp_x[_qp](0);
    _U_real(0,0) = temp00.real();
    _U_imag(0,0) = temp00.imag();
    Number temp01=_grad_disp_x[_qp](1);
    _U_real(0,1) = temp01.real();
    _U_imag(0,1) = temp01.imag();
    // Third component
    _U_real(0,2) = 0.0;
    _U_imag(0,2) = 0.0;
    
    Number temp10=_grad_disp_y[_qp](0);
    _U_real(1,0) = temp10.real();
    _U_imag(1,0) = temp10.imag();
    Number temp11=_grad_disp_y[_qp](1);
    _U_real(1,1) = temp11.real();
    _U_imag(1,1) = temp11.imag();
    // Third component
    _U_real(1,2) = 0.0;
    _U_imag(1,2) = 0.0;
    
    _U_real(2,0) = 0.0;
    _U_imag(2,0) = 0.0;
    _U_real(2,1) = 0.0;
    _U_imag(2,1) = 0.0;
    _U_real(2,2) = 0.0;
    _U_imag(2,2) = 0.0;
    
    _eps_real[_qp] = 0.5*(_U_real+_U_real.transpose());
    _eps_imag[_qp] = 0.5*(_U_imag+_U_imag.transpose());
    
    _tr_eps[_qp] = _eps_real[_qp].tr()+_imagUnit*_eps_imag[_qp].tr();

    // Complex material properties Alpha and Beta
    _alpha[_qp] = _mu[_qp]+_imagUnit*_omega[_qp]*_eta[_qp];
    _beta[_qp]  = _K[_qp]-2.0/3.0*_mu[_qp]-_imagUnit*_omega[_qp]*2.0/3.0*_eta[_qp];

    // Stress
    _sigma_real[_qp] = 2.0*(_alpha[_qp].real()*_eps_real[_qp]-_alpha[_qp].imag()*_eps_imag[_qp]);
    _sigma_real[_qp]+= (_beta[_qp].real()*_eps_real[_qp].tr()-_beta[_qp].imag()*_eps_imag[_qp].tr())*_identity;
    _sigma_imag[_qp] = 2.0*(_alpha[_qp].real()*_eps_imag[_qp]+_alpha[_qp].imag()*_eps_real[_qp]);
    _sigma_imag[_qp]+= (_beta[_qp].real()*_eps_imag[_qp].tr()+_beta[_qp].imag()*_eps_real[_qp].tr())*_identity;
    
}

