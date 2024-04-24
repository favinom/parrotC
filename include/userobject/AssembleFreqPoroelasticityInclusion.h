//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <chrono>
using namespace std::chrono;

#include "AssemblyUO.h"

#include "FreqPoroelasticInclusion.h"


#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/dirichlet_boundaries.h"

#include "FEProblem.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "Material.h"


// Forward declarations
class AssembleFreqPoroelasticityInclusion;

template <>
InputParameters validParams<AssembleFreqPoroelasticityInclusion>();

/**
 * Base class for user-specific data
 */
class AssembleFreqPoroelasticityInclusion : public AssemblyUO
{
protected:

	std::vector<int> _d_var;
	int _p_var;

	std::string const _materialName;

	bool _linearSystemAllocated;
	bool _linearSystemAssembled;

	bool const _center_x,_center_y,_center_z;

	Real _pi;
	Number _imagUnit;
	RealTensorValue _identity;

	PetscMatrix<Number> * _static_stiffness;
	PetscMatrix<Number> * _freqdep_stiffness;
	PetscMatrix<Number> * _interp;
	
	PetscVector<Number> * _dirichletFlag;
	PetscVector<Number> * _liftingFunction;
		
	void mallocLinearSystem();
	void deleteLinearSystem();

	void centerVariable(int _variable);

	virtual void assembleLinearSystem()=0;
	    
public:
  AssembleFreqPoroelasticityInclusion(const InputParameters & params);
  
  ~AssembleFreqPoroelasticityInclusion();

  virtual void execute() override;
  virtual void solve( Real _t_in );
  
};

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
//    virtual void threadJoin(const UserObject & uo) override;
