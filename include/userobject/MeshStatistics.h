//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MESHSTATISTICS
#define MESHSTATISTICS

#include "GeneralUserObject.h"
#include "FEProblem.h"
#include "MooseMesh.h"

// Forward declarations
class MeshStatistics;

template <>
InputParameters validParams<MeshStatistics>();

class MeshStatistics : public GeneralUserObject
{
protected:
    FEProblem * _fe_problem;
//    EquationSystems & _equationSystems;
    MeshBase const * _mesh;
    
public:
  MeshStatistics(const InputParameters & params);

  virtual void execute() override ;


  virtual void initialize() override ;


  virtual void finalize() override ;

};

#endif
