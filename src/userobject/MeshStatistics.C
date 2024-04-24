//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshStatistics.h"

#include "FEProblem.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"


template <>
InputParameters
validParams<MeshStatistics>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

MeshStatistics::MeshStatistics(const InputParameters & parameters) :
GeneralUserObject(parameters),
_fe_problem(parameters.get<FEProblem *>("_fe_problem"))
{
    std::cout<<"constructor called from MeshStatistics\n";
}

void MeshStatistics::execute()
{
    std::cout<<"execute\n";
    _mesh = &_fe_problem[0].mesh().getMesh();
    int nnodes = _mesh[0].n_nodes();
    int nelem = _mesh[0].n_active_elem();
    MeshBase const & mesh = _mesh[0];
    std::map<dof_id_type, std::vector<dof_id_type> > hanging_nodes;
    libMesh::MeshTools::find_hanging_nodes_and_parents(mesh,hanging_nodes);
    std::cout<<"Total nodes: "<<nnodes<<"\n";
    std::cout<<"Hanging nodes: "<<hanging_nodes.size()<<"\n";
    std::cout<<"Total elements: "<<nelem<<"\n";
    

    EquationSystems & _equationSystems(_fe_problem[0].es()) ;
    int nES=_equationSystems.n_systems();
//    std::cout<<"The EQ has "<<nES<<std::endl;
    LinearImplicitSystem & lis0=_equationSystems.get_system<LinearImplicitSystem> (0);
//    LinearImplicitSystem & lis1=_equationSystems.get_system<LinearImplicitSystem> (1);
//    std::cout<<"System lis0 has "<<lis0.n_matrices ()<<" matrices\n";
    auto &_m_sys = _fe_problem[0].es().get_system<TransientNonlinearImplicitSystem>("nl0");
//    std::cout<<_fe_problem->getCurrentExecuteOnFlag()<<std::endl;
    PetscMatrix<Number> * petsc_mat_m = dynamic_cast<PetscMatrix<Number> * >(_m_sys.matrix);
//    std::cout<<"before writing\n";
    petsc_mat_m[0].print_matlab("daUO.m");
//    std::cout<<"after writing\n";
//    exit(1);
}

void MeshStatistics::initialize()
{
    std::cout<<"initialize\n";
}

void MeshStatistics::finalize()
{
    std::cout<<"finalize\n";
}

