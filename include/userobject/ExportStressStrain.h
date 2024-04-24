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

//#include "Material.h"
#include "FreqPoroelasticFracture2D.h"

#include "libmesh/distributed_mesh.h"

#include "libmesh/transient_system.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"

#include "libmesh/serial_mesh.h"

//typedef DistributedMesh MeshClass;
//typedef SerialMesh MeshClass;
typedef ReplicatedMesh MeshClass;
typedef Nemesis_IO OutputClass;
//typedef ExodusII_IO OutputClass;

// Forward declarations
class ExportStressStrain;

template <>
InputParameters validParams<ExportStressStrain>();

/**
 * Base class for user-specific data
 */
class ExportStressStrain : public GeneralUserObject
{
protected:
		
    EquationSystems & _equationSystems;
	TransientNonlinearImplicitSystem & _tnis;
	MeshBase & _mesh;
	Parallel::Communicator const & _comm;
	DofMap const & _dof_map;
	int const _dim;
	int const _nvars;
	
	NonlinearSystemBase & _nl;
		
	std::vector<int> _d_var;
	int _p_var;
	
	std::vector< std::vector < std::string > > strain_name;
	std::vector< std::vector < std::string > > stress_name;
	
	MeshClass * _mesh2;
	EquationSystems *_equationSystems2;

	bool _allocated;
	bool _assembled;
	
	int _t_step;
	
	std::string _output_filename;
	
	std::vector< std::vector< PetscVector<Number> * > > _strainVector;
    std::vector< std::vector< PetscVector<Number> * > > _stressVector;

	PetscMatrix<Number> * _mat;
	OutputClass * outputClass;
	
	Real _pi;
	Number _imagUnit;
	
	std::string _materialName;

    void assemblePostprocessorVectors();

	    
public:
  ExportStressStrain(const InputParameters & params);
  
  ~ExportStressStrain();
  virtual void initialize() override;
  virtual void execute() override ;
  virtual void finalize() override {};
  
};
