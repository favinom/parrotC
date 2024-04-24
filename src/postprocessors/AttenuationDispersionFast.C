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

#include "AttenuationDispersionFast.h"

#include "NonlinearSystem.h"

#include "AttenuationDispersionUO.h"

template<>
InputParameters validParams<AttenuationDispersionFast>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  
  params.addRequiredParam<UserObjectName>("UserObjectName", "UserObjectName");
  params.addRequiredParam<unsigned int>("i", "i");
  params.addRequiredParam<unsigned int>("j", "j");

  MooseEnum operation_type("ATTENUATION DISPERSION STRAIN STRESS");
  params.addRequiredParam<MooseEnum>("operation_type",operation_type,"what is the postprocessor computing");

  MooseEnum component_type("REAL IMAG ALL","ALL");
  params.addRequiredParam<MooseEnum>("component",component_type,"the component you want to compute");

  return params;
}

AttenuationDispersionFast::AttenuationDispersionFast(const InputParameters & parameters) :
GeneralPostprocessor(parameters),
_userObjectName(getParam<UserObjectName>("UserObjectName")),
_operation_type(getParam<MooseEnum     >("operation_type") ),
_component_type(getParam<MooseEnum     >("component") ),
_i             (getParam<unsigned int  >("i")),
_j             (getParam<unsigned int  >("j"))
{}

Number
AttenuationDispersionFast::getValue()
{
	AttenuationDispersionUO const & attenuationDispersionUO=_fe_problem.getUserObject<AttenuationDispersionUO>(_userObjectName);
	
	Number strain=attenuationDispersionUO.getStrainComponent(_i,_j);
	Number stress=attenuationDispersionUO.getStressComponent(_i,_j);

	Number H=stress/strain; 
	if (_i != _j)
		H=H/2.0;

	Number temp;

	switch (_operation_type)
	{
		case 0: //ATTENUATION

		temp = H.imag()/H.real();
		
		break;
		case 1: // DISPERSION
		{
		Real _pp_val_stress_real=stress.real();
		Real _pp_val_stress_imag=stress.imag();
		Real _pp_val_strain_real=strain.real();
		Real _pp_val_strain_imag=strain.imag();
		
		temp = _pp_val_stress_real*_pp_val_strain_real + _pp_val_stress_imag*_pp_val_strain_imag;

		temp=temp/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);	
		}
		break;
		case 2: // STRAIN
		temp = strain;
		break;
		case 3: // STRESS
		temp = stress;
		break;
		// we should put a default with an error
	}

	switch (_component_type)
	{
		case 0: // REAL

		temp = temp.real();
		
		break;
		case 1: // IMAG
		temp = temp.imag();
		break;

	}
	return temp;
}

// Real _pp_val_strain_real=strain.real();
// Real _pp_val_strain_imag=strain.imag();
// Real _pp_val_stress_real=stress.real();
// Real _pp_val_stress_imag=stress.imag();
// std::cout<<"_pp_val_strain_real "<<_pp_val_strain_real<<std::endl;
// std::cout<<"_pp_val_strain_imag "<<_pp_val_strain_imag<<std::endl;
// std::cout<<"_pp_val_stress_real "<<_pp_val_stress_real<<std::endl;
// std::cout<<"_pp_val_stress_imag "<<_pp_val_stress_imag<<std::endl;
// Real nom = _pp_val_strain_real*_pp_val_stress_imag - _pp_val_stress_real*_pp_val_strain_imag;
// nom = nom/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);
// Real denom = _pp_val_stress_real*_pp_val_strain_real + _pp_val_stress_imag*_pp_val_strain_imag;
// denom = denom/(_pp_val_strain_real*_pp_val_strain_real + _pp_val_strain_imag*_pp_val_strain_imag);
// //return nom/denom;
