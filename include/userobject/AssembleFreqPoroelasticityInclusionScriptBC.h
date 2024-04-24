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
class AssembleFreqPoroelasticityInclusionScriptBC;

template <>
InputParameters validParams<AssembleFreqPoroelasticityInclusionScriptBC>();

/**
 * Base class for user-specific data
 */
class AssembleFreqPoroelasticityInclusionScriptBC : public AssembleFreqPoroelasticityInclusion
{
protected:


	virtual void assembleLinearSystem();
	    
public:
  AssembleFreqPoroelasticityInclusionScriptBC(const InputParameters & params);
    
};
