//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblem.h"

#include "AssembleFreqPoroelasticityInclusion.h"

class UserObjectAssemblyProblem;

template <>
InputParameters validParams<UserObjectAssemblyProblem>();

class UserObjectAssemblyProblem : public FEProblem
{
public:
    UserObjectAssemblyProblem(const InputParameters & parameters);    
    
    virtual void computeJacobian(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & ,
                            SparseMatrix<Number> & );
    
    virtual void computeResidual(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & ,
                            NumericVector<Number> & );

	virtual void solve();
	
	bool converged();
    
	UserObjectName const userObjectName;
	
	//AssembleFreqPoroelasticityInclusion * _assembleVectorsAndMatrices;

};

