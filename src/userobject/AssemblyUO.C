//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssemblyUO.h"

//registerMooseObject("parrotApp", AssemblyUO);

template <>
InputParameters
	validParams<AssemblyUO>()
{
	InputParameters params = validParams<GeneralUserObject>();
	return params;
}

AssemblyUO::AssemblyUO(const InputParameters & parameters) :
GeneralUserObject(parameters),
_nl( _fe_problem.getNonlinearSystemBase() ),
_equationSystems(_fe_problem.es()),
_tnis( _equationSystems.get_system<TransientNonlinearImplicitSystem>("nl0") ),
_mesh( _equationSystems.get_mesh() ),
_comm ( _mesh.comm() ),
_dof_map( _tnis.get_dof_map() ),
_dim( _mesh.mesh_dimension() ),
_nvars( _tnis.n_vars() ),
_mat_SM_ptr(_tnis.matrix),
_rhs_NV_ptr(_tnis.rhs),
//_mat_SM(*_tnis.matrix),
//_rhs_NV(*_tnis.rhs),
_sol_NV( _nl.solution() ),
_currentSolution_ptr( _nl.currentSolution() ),
_curr_NV_uni_ptr(_tnis.solution),
_curr_NV_ptr( _curr_NV_uni_ptr.get() ),
_curr_NV(_curr_NV_ptr[0]),
_parrotSolver(1, _comm ),
myout(std::cout.rdbuf())
{
	std::cout<<"AssemblyUO::AssemblyUO() start\n";
	// This error should never happen as the system nl0 has been retrieved above
	if( !_equationSystems.has_system("nl0") )
	{
		mooseError("The system nl0 does not exist!");
	}

	myout<<std::fixed;
	myout<<std::setprecision(1);

	std::cout<<"AssemblyUO::AssemblyUO() stop\n";
}
 
AssemblyUO::~AssemblyUO(){}
