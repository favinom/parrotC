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

#include "FreqPoroelasticFracture2DInvasion.h"

#include <sstream>
#include <fstream>
#include "MooseMesh.h"

template<>
InputParameters validParams<FreqPoroelasticFracture2DInvasion>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<int>("fn", "number of fractures");
  params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fl_string", "length of fractures");
  params.addRequiredParam<std::string>("ft_string", "thickness of fractures");
  params.addRequiredParam<std::string>("fa_string", "angle of fractures");
  params.addRequiredParam<std::string>("filename_invasion", "filename of the invasion file");
  params.addRequiredParam<int>("nx_invasion", "number of elements in invasion matrix");
  params.addRequiredParam<Real>("dx_invasion", "spacing of invasion matrix");
  params.addRequiredParam<Real>("mu_block", "second Lamé parameter");
  params.addRequiredParam<Real>("lambda_block", "first Lamé parameter");
  params.addRequiredParam<Real>("alpha_block", "another material parameter related to pressure");
  params.addRequiredParam<Real>("kappa_block", "permeability");
  params.addRequiredParam<Real>("eta_block", "viscosity");
  params.addRequiredParam<Real>("porosity_block", "porosity");
  params.addRequiredParam<Real>("kf_block", "bulk_modulus_fluid");
  params.addRequiredParam<Real>("ks_block", "bulk_modulus_solid");
  params.addRequiredParam<Real>("mu_fracture", "second Lamé parameter");
  params.addRequiredParam<Real>("lambda_fracture", "first Lamé parameter");
  params.addRequiredParam<Real>("alpha_fracture", "another material parameter related to pressure");
  params.addRequiredParam<Real>("kappa_fracture", "permeability");
  params.addRequiredParam<Real>("eta_fracture_invaded", "viscosity of invaded fracture");
  params.addRequiredParam<Real>("eta_fracture_brine", "viscosity of not invaded fracture");
  params.addRequiredParam<Real>("porosity_fracture", "porosity");
  params.addRequiredParam<Real>("kf_fracture_invaded", "bulk_modulus_fluid of invaded fracture");
  params.addRequiredParam<Real>("kf_fracture_brine", "bulk_modulus_fluid of not invaded fracture");
  params.addRequiredParam<Real>("ks_fracture", "bulk_modulus_solid");
  params.addRequiredCoupledVar("disp_x", "complex first  coupled component");
  params.addRequiredCoupledVar("disp_y", "complex second coupled component");

  return params;
}

FreqPoroelasticFracture2DInvasion::FreqPoroelasticFracture2DInvasion(const InputParameters & parameters) :
    Material(parameters),
    _fn(getParam<int>("fn")), 
    _fx_string(getParam<std::string>("fx_string")), 
    _fy_string(getParam<std::string>("fy_string")), 
    _fl_string(getParam<std::string>("fl_string")), 
    _ft_string(getParam<std::string>("ft_string")), 
    _fa_string(getParam<std::string>("fa_string")), 
    _filename_invasion(getParam<std::string>("filename_invasion")), 
    _nx_invasion(getParam<int>("nx_invasion")), 
    _dx_invasion(getParam<Real>("dx_invasion")),
    _eps_real(declareProperty<RealTensorValue>("strain_real")),
    _eps_imag(declareProperty<RealTensorValue>("strain_imag")),
    _tr_eps(declareProperty<Number>("trace_strain")),
    _sigma_real(declareProperty<RealTensorValue>("stress_real")),
    _sigma_imag(declareProperty<RealTensorValue>("stress_imag")),
    _mu_block(getParam<Real>("mu_block")),
    _lambda_block(getParam<Real>("lambda_block")),
    _alpha_block(getParam<Real>("alpha_block")),
    _kappa_block(getParam<Real>("kappa_block")),
    _eta_block(getParam<Real>("eta_block")),
    _porosity_block(getParam<Real>("porosity_block")),
    _kf_block(getParam<Real>("kf_block")),
    _ks_block(getParam<Real>("ks_block")),
    _mu_fracture(getParam<Real>("mu_fracture")),
    _lambda_fracture(getParam<Real>("lambda_fracture")),
    _alpha_fracture(getParam<Real>("alpha_fracture")),
    _kappa_fracture(getParam<Real>("kappa_fracture")),
    _eta_fracture_invaded(getParam<Real>("eta_fracture_invaded")),
    _eta_fracture_brine(getParam<Real>("eta_fracture_brine")),
    _porosity_fracture(getParam<Real>("porosity_fracture")),
    _kf_fracture_invaded(getParam<Real>("kf_fracture_invaded")),
    _kf_fracture_brine(getParam<Real>("kf_fracture_brine")),
    _ks_fracture(getParam<Real>("ks_fracture")),
    _omega(declareProperty<Real>("omega_property")),
    _mu(declareProperty<Real>("mu_property")),
    _lambda(declareProperty<Real>("lambda_property")),
    _alpha(declareProperty<Real>("alpha_property")),
    _kappa(declareProperty<Real>("kappa_property")),
    _eta(declareProperty<Real>("eta_property")),
    _porosity(declareProperty<Real>("porosity_property")),
    _kf(declareProperty<Real>("kf_property")),
    _ks(declareProperty<Real>("ks_property")),
    _inv_m(declareProperty<Real>("inverse_of_m")),
    _diffusion(declareProperty<Real>("diffusion_property")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y")),
    _imagUnit(0.0,1.0)
{
    
    // std::cout<<"Constructor called\n";
    
    if (_mesh.dimension() != 2)
    {
        std::cout<<"You cannot use this material in dimension different from 2\n";
        exit(1);
    }

    _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
  
    std::cout << "This is cpu " << _communicator.rank() << std::endl;

    std::streampos filesize;
	short temp=0;
    

    std::ifstream myfile(_filename_invasion, std::ifstream::binary); 
    if(myfile.fail())
    {
        std::cerr<<"Cannot open the bin file, exiting...\n";
        exit(1);
    }
    myfile.seekg(0, std::ios::end);
    filesize = myfile.tellg();
    _invasion_data = new short [filesize/2];
    myfile.seekg(0, std::ios::beg);
    for (int i=0; i<filesize/2; i++) {
        myfile.read((char*)&temp,sizeof(temp)); // the char-pointer is used, so that exactly 1 byte is read
        _invasion_data[i] = temp;
    }
    myfile.close();
    // std::cout << "Showing data" << std::endl;
    // std::cout << _invasion_data[3999278] << std::endl;

    // std::string filename_out = "temp_out";
    // filename_out=filename_out+std::to_string(_communicator.rank())+".txt";
    // std::ofstream outfile (filename_out);
    // if (outfile.is_open())
    // {
    //     for (int i=0; i<filesize/2; i++) {
    //         outfile << _invasion_data[i] << "\n";
    //     }
    //     outfile.close();
    // }
 
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
}

void
FreqPoroelasticFracture2DInvasion::computeQpProperties()
{

    Real x_coord = _q_point[_qp](0);
    Real y_coord = _q_point[_qp](1);


    _omega[_qp] = 2.0*_pi*std::pow(10,_t);
    // We set parameters to block parameters
    _mu[_qp] = _mu_block;
    _lambda[_qp] = _lambda_block;
    _alpha[_qp] = _alpha_block;
    _kappa[_qp] = _kappa_block;
    _eta[_qp] = _eta_block;
    _porosity[_qp] = _porosity_block;
    _kf[_qp] = _kf_block;
    _ks[_qp] = _ks_block;
    
    
    int temp_xcor = 0;
    int temp_ycor = 0;
    Real temp_shift = _nx_invasion*_dx_invasion/2.0;   

     
    for (int i = 0; i < _fn; i++)
    {
        Real temp =_a[i]*x_coord+_b[i]*y_coord + _c[i];
        if ( std::fabs(temp) <=  _ft[i]/2.0 )
        {
            Real tempo=_ao[i]*x_coord+_bo[i]*y_coord + _co[i];
            if (std::fabs(tempo) <=  _fl[i]/2.0 )
            {
                _mu[_qp] =  _mu_fracture;
                _lambda[_qp] = _lambda_fracture;
                _alpha[_qp] = _alpha_fracture;
                _kappa[_qp] = _kappa_fracture;
                //_kappa[_qp] = _ft[i]*_ft[i]/9.6e4;
                _porosity[_qp] = _porosity_fracture;
                _ks[_qp] = _ks_fracture;
                
                _eta[_qp] = _eta_fracture_brine;
                _kf[_qp] = _kf_fracture_brine;
                temp_ycor = round((temp_shift+y_coord)/_dx_invasion);
                temp_xcor = round((temp_shift+x_coord)/_dx_invasion);
                if (_invasion_data[(temp_xcor+1)*_nx_invasion-temp_ycor-1] == 1) 
                { // saturated fracture (CO2)
                    _eta[_qp] = _eta_fracture_invaded;
                    _kf[_qp] = _kf_fracture_invaded;
                }
                break;
            }
        }
    }
    // The sign in the paper of Quintal et al. 2011 is wrong
    _inv_m[_qp] = _porosity[_qp]/_kf[_qp] + (_alpha[_qp]-_porosity[_qp])/_ks[_qp];
    _diffusion[_qp] = _kappa[_qp]/_eta[_qp]/_omega[_qp];
    
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
    
    // Stress
    _sigma_real[_qp] = 2.0*_mu[_qp]*_eps_real[_qp]+_lambda[_qp]*_eps_real[_qp].tr()*_identity;
    _sigma_imag[_qp] = 2.0*_mu[_qp]*_eps_imag[_qp]+_lambda[_qp]*_eps_imag[_qp].tr()*_identity;
    
    
}

