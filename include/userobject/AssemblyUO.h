//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once


#include "GeneralUserObject.h"

// MOOSE includes

#include "NonlinearSystemBase.h"

// libMesh includes

#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"

// libMesh includes linear algebra

#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

// include my solver
#include "ParrotSolver.h"



class AssemblyUO;

template <>
InputParameters validParams<AssemblyUO>();

/**
 * Base class for user-specific data
 */
class AssemblyUO : public GeneralUserObject
{
protected:

	// MOOSE objects
	NonlinearSystemBase & _nl;

	// libmesh objects
	EquationSystems & _equationSystems;
	TransientNonlinearImplicitSystem & _tnis;
	MeshBase const & _mesh;
	Parallel::Communicator const & _comm;
	DofMap const & _dof_map;
	int const _dim;
	int const _nvars;

	// linear albegra objects
	SparseMatrix<Number>  * & _mat_SM_ptr;
	NumericVector<Number> * & _rhs_NV_ptr;

	NumericVector<Number> & _sol_NV;
	NumericVector<Number> const * & _currentSolution_ptr;
	std::unique_ptr<NumericVector<Number> > & _curr_NV_uni_ptr;
	NumericVector<Number> * _curr_NV_ptr;
	NumericVector<Number> & _curr_NV;

	ParrotSolver _parrotSolver;

	std::ostream myout;

public:
  AssemblyUO(const InputParameters & params);

  ~AssemblyUO();

  virtual void initialize() {}  ;  
  virtual void execute   () = 0 ;
  virtual void finalize  () {}  ;

  virtual void assembleLinearSystem()=0;
  virtual void solve(Real _t) = 0;

};
