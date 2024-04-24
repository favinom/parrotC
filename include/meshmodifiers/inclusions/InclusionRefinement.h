//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MeshModifier.h"
#include "FEProblem.h"
//#include "MooseEnum.h"

#include "libmesh/distributed_mesh.h"

#include "InclusionsMeshModifier.h"


// Forward declerations
class InclusionRefinement;

template <>
InputParameters validParams<InclusionRefinement>();

/**
 * MeshModifier for defining a Subdomain inside or outside of a bounding box
 */
class InclusionRefinement : public MeshModifier
{
public:
    InclusionRefinement(const InputParameters & parameters);
    
    virtual void modify() override;

    void doAMR();
    void doUMR(int i);
    void doRefine(std::vector<int> const & refinements);
    
protected:
    
    std::string const _meshModifierName;
   	MeshModifier const & _myMeshModifier;
	InclusionsMeshModifier const & _inclusionsMeshModifier;
	MeshBase * _meshBase;
	//DistributedMesh * distributedMesh;
	UnstructuredMesh * distributedMesh;
	MeshRefinement * meshRefinement;

	bool _hasRefinementVector;
	std::vector<int> _refinements;

	bool _hasOutputFileName;
	std::string _outputFileName;
	
	bool _doBoundary;
};
