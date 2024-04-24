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

#include "FreqNavierStokesFracture3D.h"

#include <sstream>
#include "MooseMesh.h"

template<>
InputParameters validParams<FreqNavierStokesFracture3D>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<int>("fn", "number of fractures");
  params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fz_string", "z-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fa1_string", "rotation along z axis");
  params.addRequiredParam<std::string>("fa2_string", "rotation along y axis");
  params.addRequiredParam<std::string>("fa3_string", "rotation along x axis");
  params.addRequiredParam<std::string>("fd1_string", "fracture dimension 1");
  params.addRequiredParam<std::string>("fd2_string", "fracture dimension 2");
  params.addRequiredParam<std::string>("fd3_string", "fracture dimension 3");
  params.addRequiredParam<Real>("mu_block", "shear modulus");
  params.addRequiredParam<Real>("K_block", "bulk modulus");
  params.addRequiredParam<Real>("eta_block", "shear viscocity");
  params.addRequiredParam<Real>("mu_fracture", "shear modulus");
  params.addRequiredParam<Real>("K_fracture", "bulk modulus");
  params.addRequiredParam<Real>("eta_fracture", "shear viscocity");
  params.addRequiredCoupledVar("disp_x", "complex first  coupled component");
  params.addRequiredCoupledVar("disp_y", "complex second coupled component");
  params.addRequiredCoupledVar("disp_z", "complex second coupled component");

  return params;
}

FreqNavierStokesFracture3D::FreqNavierStokesFracture3D(const InputParameters & parameters) :
    Material(parameters),
    _fn(getParam<int>("fn")), 
    _fx_string(getParam<std::string>("fx_string")), 
    _fy_string(getParam<std::string>("fy_string")), 
    _fz_string(getParam<std::string>("fz_string")), 
    _fa1_string(getParam<std::string>("fa1_string")),
    _fa2_string(getParam<std::string>("fa2_string")),
    _fa3_string(getParam<std::string>("fa3_string")),
    _fd1_string(getParam<std::string>("fd1_string")),
    _fd2_string(getParam<std::string>("fd2_string")),
    _fd3_string(getParam<std::string>("fd3_string")),
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
    _grad_disp_z(coupledGradient("disp_z")),
    _imagUnit(0.0,1.0)
{
    
    
    if (_mesh.dimension() != 3)
    {
        std::cout<<"You cannot use this material in dimension different from 3\n";
        exit(1);
    }

    _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    
    
    _pi = acos(-1.0);
    
    _center   =new RealVectorValue [_fn];
    _rotation =new RealVectorValue [_fn];
    _dimension=new RealVectorValue [_fn];
    _d        =new RealVectorValue [_fn];

    _n = new RealVectorValue * [_fn];
    for (int i=0; i<_fn; ++i)
    {
        _n[i]=new RealVectorValue [3];
    }

    std::istringstream  fx_ss( _fx_string);
    std::istringstream  fy_ss( _fy_string);
    std::istringstream  fz_ss( _fz_string);
    std::istringstream fa1_ss(_fa1_string);
    std::istringstream fa2_ss(_fa2_string);
    std::istringstream fa3_ss(_fa3_string);
    std::istringstream fd1_ss(_fd1_string);
    std::istringstream fd2_ss(_fd2_string);
    std::istringstream fd3_ss(_fd3_string);

    std::string token;

    for (int i=0; i<_fn; ++i)
    {
        std::getline(fx_ss, token, ',');
        _center[i](0)=std::atof(token.c_str());
        std::getline(fy_ss, token, ',');
        _center[i](1)=std::atof(token.c_str());
        std::getline(fz_ss, token, ',');
        _center[i](2)=std::atof(token.c_str());

        std::getline(fa1_ss, token, ',');
        _rotation[i](0)=atof(token.c_str());
        std::getline(fa2_ss, token, ',');
        _rotation[i](1)=atof(token.c_str());
        std::getline(fa3_ss, token, ',');
        _rotation[i](2)=atof(token.c_str());
        _rotation[i]=_rotation[i]/180.0*pi;

        std::getline(fd1_ss, token, ',');
        _dimension[i](0)=atof(token.c_str());
        std::getline(fd2_ss, token, ',');
        _dimension[i](1)=atof(token.c_str());
        std::getline(fd3_ss, token, ',');
        _dimension[i](2)=atof(token.c_str());

        ComputeNormalsFromAngles(_rotation[i],_n[i][0],_n[i][1],_n[i][2]);

        for (int j=0; j<3; ++j)
        {
            _d[i](j)=_n[i][j]*_center[i];
        }
     }

    _fx_string.clear();
    _fy_string.clear();
    _fz_string.clear();
    _fa1_string.clear();
    _fa2_string.clear();
    _fa3_string.clear();
    _fd1_string.clear();
    _fd2_string.clear();
    _fd3_string.clear();

    delete [] _center;
    delete [] _rotation;
}

void
FreqNavierStokesFracture3D::computeQpProperties()
{
    _dummy_alpha[_qp] = 0.0;

    RealVectorValue point;
    point(0) = _q_point[_qp](0);
    point(1) = _q_point[_qp](1);
    point(2) = _q_point[_qp](2);

    _omega[_qp] = 2.0*_pi*std::pow(10,_t);
    // We set parameters to block parameters
    _mu[_qp] = _mu_block;
    _K[_qp] = _K_block;
    _eta[_qp] = _eta_block;
    
    if ( is_inside(point))
    {
        _mu[_qp] =  _mu_fracture;
        _K[_qp] = _K_fracture;
        _eta[_qp] = _eta_fracture;
    }
    
    // Kinematics
    Number temp00=_grad_disp_x[_qp](0);
    _U_real(0,0) = temp00.real();
    _U_imag(0,0) = temp00.imag();
    Number temp01=_grad_disp_x[_qp](1);
    _U_real(0,1) = temp01.real();
    _U_imag(0,1) = temp01.imag();
    Number temp02=_grad_disp_x[_qp](2);
    _U_real(0,2) = temp02.real();
    _U_imag(0,2) = temp02.imag();
    
    Number temp10=_grad_disp_y[_qp](0);
    _U_real(1,0) = temp10.real();
    _U_imag(1,0) = temp10.imag();
    Number temp11=_grad_disp_y[_qp](1);
    _U_real(1,1) = temp11.real();
    _U_imag(1,1) = temp11.imag();
    Number temp12=_grad_disp_y[_qp](2);
    _U_real(1,2) = temp12.real();
    _U_imag(1,2) = temp12.imag();
    
    Number temp20=_grad_disp_z[_qp](0);
    _U_real(2,0) = temp20.real();
    _U_imag(2,0) = temp20.imag();
    Number temp21=_grad_disp_z[_qp](1);
    _U_real(2,1) = temp21.real();
    _U_imag(2,1) = temp21.imag();
    Number temp22=_grad_disp_z[_qp](2);
    _U_real(2,2) = temp22.real();
    _U_imag(2,2) = temp22.imag();
    
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

void FreqNavierStokesFracture3D::ComputeNormalsFromAngles(RealVectorValue const & angles,
                                              RealVectorValue & n1,
                                              RealVectorValue & n2,
                                              RealVectorValue & n3)
{
    RealTensorValue R1;
    RealTensorValue R2;
    RealTensorValue R3;

    R1(0,0)=std::cos(angles(0));
    R1(0,1)=-std::sin(angles(0));
    R1(0,2)=0.0;
    R1(1,0)=std::sin(angles(0));
    R1(1,1)=std::cos(angles(0));
    R1(1,2)=0.0;
    R1(2,0)=0.0;
    R1(2,1)=0.0;
    R1(2,2)=1.0;

    R2(0,0)=std::cos(angles(1));
    R2(0,1)=0.0;
    R2(0,2)=-std::sin(angles(1));
    R2(1,0)=0.0;
    R2(1,1)=1.0;
    R2(1,2)=0.0;
    R2(2,0)=std::sin(angles(1));
    R2(2,1)=0.0;
    R2(2,2)=std::cos(angles(1));

    R3(0,0)=1.0;
    R3(0,1)=0.0;
    R3(0,2)=0.0;
    R3(1,0)=0.0;
    R3(1,1)=std::cos(angles(2));
    R3(1,2)=-std::sin(angles(2));
    R3(2,0)=0.0;
    R3(2,1)=std::sin(angles(2));
    R3(2,2)=std::cos(angles(2));

    RealTensorValue R=R1*R2*R3;
    std::cout<<R<<std::endl;

    for (int i=0; i<3; ++i)
    {
        n1(i)=R(i,0);
        n2(i)=R(i,1);
        n3(i)=R(i,2);
    }
}

bool FreqNavierStokesFracture3D::is_inside(RealVectorValue const & point)
{
    for (int i=0; i<_fn; ++i)
    {
        Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
        if (temp1<_dimension[i](0)/2.0)
        {
            Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
            if (temp2<_dimension[i](1)/2.0)
            {
                Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
                if (temp3<_dimension[i](2)/2.0)
                {
                    return true;
                }

            }
        }
    }
    return false;
}

FreqNavierStokesFracture3D::~FreqNavierStokesFracture3D()
{
    delete [] _dimension;
    for (int i=0; i<_fn; ++i)
    {
        delete [] _n[i];
    }
    delete [] _n;
    delete [] _d;
}

