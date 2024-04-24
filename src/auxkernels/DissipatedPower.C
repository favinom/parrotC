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
/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Auxiliary Kernel to visualize the fibers.                    */
/****************************************************************/

#include "DissipatedPower.h"


template <>
InputParameters
validParams<DissipatedPower>()
{
  // inherit the parameters of AuxKernel:
  InputParameters params = validParams<AuxKernel>();

  //specify for which component we consider:
    params.addRequiredCoupledVar("pres", "pressure");
    params.addRequiredParam<PostprocessorName>("strain_energy_pp", "");
  
  return params;
}

DissipatedPower::DissipatedPower(const InputParameters &parameters):
    AuxKernel(parameters),
    _kappa(getMaterialProperty<Real>("kappa_property")),
    _eta(getMaterialProperty<Real>("eta_property")),
    _omega(getMaterialProperty<Real>("omega_property")),
    _grad_pres(coupledGradient("pres")),
    _pp_value(getPostprocessorValue("strain_energy_pp"))
{
}

Number
DissipatedPower::computeValue()
{
    
    RealVectorValue grad_real;
    RealVectorValue grad_imag;
    for (int ii=0; ii<3; ++ii)
    {
        Number comp=_grad_pres[_qp](ii);
        grad_real(ii)=comp.real();
        grad_imag(ii)=comp.imag();
    }

    Number a= 1.0/2.0*_kappa[_qp]/_eta[_qp]*(grad_real*grad_real + grad_imag*grad_imag);
    
//    std::string filename_out = "temp_out";
//    filename_out=filename_out+std::to_string(_communicator.rank())+".txt";
//    std::ofstream outfile; 
//    outfile.open(filename_out,std::ios_base::app);
//    if (outfile.is_open())
//    {
//        outfile << a/_pp_value/2.0/_omega[_qp] << "\n";
//        outfile.close();
//    }
    
    return a/_pp_value/2.0/_omega[_qp];
//   return 1.0; 
}
