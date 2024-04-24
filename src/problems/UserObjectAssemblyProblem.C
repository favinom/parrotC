//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UserObjectAssemblyProblem.h"

#include "NonlinearSystemBase.h"

//#include "petscmat.h"

//#include "chrono"

//registerMooseObject("parrotApp", UserObjectAssemblyProblem);

template <>
InputParameters
validParams<UserObjectAssemblyProblem>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<UserObjectName>("assembly_userobject","the user object where assembly is performed.");
    return params;
}


UserObjectAssemblyProblem::UserObjectAssemblyProblem(const InputParameters & parameters) :
FEProblem(parameters),
userObjectName ( getParam<UserObjectName>("assembly_userobject") )
{}


void
UserObjectAssemblyProblem::computeResidual(NonlinearImplicitSystem & /*sys*/,
                                   const NumericVector<Number> & ,
                                   NumericVector<Number> & )
{
	_console<<"UserObjectAssemblyProblem::computeResidualSys\n";
}

void
UserObjectAssemblyProblem::computeJacobian(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & ,
                                  SparseMatrix<Number> & )
{
		_console<<"UserObjectAssemblyProblem::computeJacobianSys\n";
}

void UserObjectAssemblyProblem::solve(){
	_console<<"UserObjectAssemblyProblem::solve()\n";
	AssembleFreqPoroelasticityInclusion const & _assembleVectorsAndMatrices=getUserObject<AssembleFreqPoroelasticityInclusion>(userObjectName);
	AssembleFreqPoroelasticityInclusion & pj = const_cast<AssembleFreqPoroelasticityInclusion &>(_assembleVectorsAndMatrices);
	pj.solve(_time);
	_nl->update();
}


bool UserObjectAssemblyProblem::converged()
{
	_console<<"UserObjectAssemblyProblem::converged()\n";
	return true;
}
