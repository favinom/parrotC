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
//#include "FreqPoroelasticFracture2D.h"

#include "libmesh/transient_system.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

// Forward declarations
class AttenuationDispersionUO;

template <>
InputParameters validParams<AttenuationDispersionUO>();

/**
 * Base class for user-specific data
 */
class AttenuationDispersionUO : public GeneralUserObject
{
protected:
    
	// Reference to the material (parent class)
	//Material const & _mat;
	// Reference to the specific (freqPoroelasticFracture2D)
	//FreqPoroelasticFracture2D const & freqPoroelasticFracture2D;
	
    EquationSystems & _equationSystems;
	TransientNonlinearImplicitSystem & _tnis;
	MeshBase const & _mesh;
	Parallel::Communicator const & _comm;
	DofMap const & _dof_map;
	int const _dim;
	int const _nvars;
	
	NonlinearSystemBase & _nl;
	
	SparseMatrix<Number>  & _mat_SM;
//	PetscMatrix<Number>   & _mat_PM;
	NumericVector<Number> & _rhs_NV;
//	PetscVector<Number>   & _rhs_PV;
	NumericVector<Number> & _sol_NV;
//	PetscVector<Number>   & _sol_PV;
	
	std::vector<int> _d_var;
	int _p_var;
	
	
	std::vector< std::vector< PetscVector<Number> * > > _strainVector;
    std::vector< std::vector< PetscVector<Number> * > > _stressVector;

	bool _postprocessorVectorsAllocated;
	bool _postprocessorVectorsAssembled;

	NumberTensorValue _strainIntegral;
	NumberTensorValue _stressIntegral;
	
	Real _pi;
	Number _imagUnit;

	std::string const _materialName;
	
	void mallocPostprocessorVectors();
	void deletePostprocessorVectors();
	void mallocLinearSystem();
	void deleteLinearSystem();

    void assemblePostprocessorVectors();
	    
public:
  AttenuationDispersionUO(const InputParameters & params);
  
  ~AttenuationDispersionUO();
  virtual void initialize() override;
  virtual void execute() override ;
  virtual void finalize() override {};
  
  Number const getStrainComponent(int i, int j) const;
  Number const getStressComponent(int i, int j) const;

};

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
//    virtual void threadJoin(const UserObject & uo) override;
