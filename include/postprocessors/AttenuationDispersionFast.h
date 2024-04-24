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

#pragma once

#include "GeneralPostprocessor.h"

class AttenuationDispersionFast;

template<>
InputParameters validParams<AttenuationDispersionFast>();

class AttenuationDispersionFast : public GeneralPostprocessor
{
public:
  AttenuationDispersionFast(const InputParameters & parameters);
  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}

protected:
  virtual Number getValue();
  
  UserObjectName _userObjectName;

  MooseEnum const _operation_type;
  MooseEnum const _component_type;
  
  unsigned int const _i;
  unsigned int const _j;
};


//   virtual Real getValue() override;
//   virtual void initialize() override {}
//   virtual void execute() override {}
//   
