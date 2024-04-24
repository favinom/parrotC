//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AssembleFreqPoroelasticityInclusion.h"

// Forward declarations
class AssembleFreqPoroelasticityInclusionPeriodicBC;

template <>
InputParameters validParams<AssembleFreqPoroelasticityInclusionPeriodicBC>();

/**
 * Base class for user-specific data
 */
class AssembleFreqPoroelasticityInclusionPeriodicBC : public AssembleFreqPoroelasticityInclusion
{
protected:


	virtual void assembleLinearSystem();
	virtual void solve( Real _t_in );
	
	std::vector<Real> _bc_value_real;
	std::vector<Real> _bc_value_imag;
	MooseEnum const   _direction;
	    
public:
  AssembleFreqPoroelasticityInclusionPeriodicBC(const InputParameters & params);
    
};
