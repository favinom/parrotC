#include "parrotcApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "DiffusionC.h"
#include "DiffusionMaterialsC.h"
#include "FreqPoroelasticityFast.h"
#include "FreqStaticLinearizedNavierStokesFast.h"

#include "ViscousShearing.h"
#include "MatOut.h"
#include "DissipatedPower.h"

#include "ResidualForcingNodalKernel.h"

#include "FreqPoroelasticSphere.h"
#include "FreqPoroelasticFracture2D.h"
#include "FreqPoroelasticFracture2DAngleLimit.h"
#include "FreqPoroelasticFracture2DInvasion.h"
#include "FreqNavierStokesFracture2D.h"
#include "FreqNavierStokesFracture3D.h"
#include "FreqDiffusionSphere.h"
#include "FreqDiffusionCube.h"
#include "FreqPoroelasticInclusion.h"

#include "DirichletBCC.h"

#include "SphereMarker.h"
#include "CubeMarker.h"
#include "SimpleMarker.h"
#include "SimpleMarker3D.h"

#include "DualProjectionAMR.h"
#include "DualProjectionAMR3D.h"

#include "AttenuationDispersion.h"
#include "AttenuationDispersionFast.h"
#include "MeanMaterialProperty.h"
#include "StrainEnergy.h"

#include "ZeroPredictor.h"

#include "MeshStatistics.h"
#include "AssembleFreqPoroelasticityInclusionScriptBC.h"
#include "AssembleFreqPoroelasticityInclusionPeriodicBC.h"
#include "ExportStressStrain.h"

#include "UserObjectAssemblyProblem.h"
#include "AttenuationDispersionUO.h"

#include "InclusionsMeshModifier.h"
#include "InclusionRefinement.h"

template <>
InputParameters
validParams<parrotcApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

parrotcApp::parrotcApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  parrotcApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  parrotcApp::associateSyntax(_syntax, _action_factory);
}

parrotcApp::~parrotcApp() {}

// External entry point for dynamic application loading
extern "C" void
parrotcApp__registerApps()
{
  parrotcApp::registerApps();
}
void
parrotcApp::registerApps()
{
  registerApp(parrotcApp);
}

// External entry point for dynamic object registration
extern "C" void
parrotcApp__registerObjects(Factory & factory)
{
  parrotcApp::registerObjects(factory);
}
void
parrotcApp::registerObjects(Factory & factory)
{
    registerKernel(DiffusionC);
    registerKernel(DiffusionMaterialsC);
    registerKernel(FreqPoroelasticityFast);
    registerKernel(FreqStaticLinearizedNavierStokesFast);

    registerAux(ViscousShearing);
    registerAux(MatOut);
    registerAux(DissipatedPower);
    
    registerNodalKernel(ResidualForcingNodalKernel);
    
    registerMaterial(FreqPoroelasticSphere);
    registerMaterial(FreqPoroelasticFracture2D);
    registerMaterial(FreqPoroelasticFracture2DAngleLimit);
    registerMaterial(FreqPoroelasticFracture2DInvasion);
    registerMaterial(FreqNavierStokesFracture2D);
    registerMaterial(FreqNavierStokesFracture3D);
    registerMaterial(FreqDiffusionSphere);
    registerMaterial(FreqDiffusionCube);
    registerMaterial(FreqPoroelasticInclusion);
    
    registerBoundaryCondition(DirichletBCC);
    
    registerMarker(SphereMarker);
    registerMarker(CubeMarker);
    registerMarker(SimpleMarker);
    registerMarker(SimpleMarker3D);
    
    registerTransfer(DualProjectionAMR);
    registerTransfer(DualProjectionAMR3D);
    
    registerPostprocessor(AttenuationDispersion);
	registerPostprocessor(AttenuationDispersionFast);
    registerPostprocessor(MeanMaterialProperty);
    registerPostprocessor(StrainEnergy);
    
    registerPredictor(ZeroPredictor);
    
    registerUserObject(MeshStatistics);
    registerUserObject(AssembleFreqPoroelasticityInclusionScriptBC);
    registerUserObject(AssembleFreqPoroelasticityInclusionPeriodicBC);
    registerUserObject(ExportStressStrain);
    registerUserObject(AttenuationDispersionUO);
	
	registerProblem(UserObjectAssemblyProblem);
	
	registerMeshModifier(InclusionsMeshModifier);
	registerMeshModifier(InclusionRefinement);
}

// External entry point for dynamic syntax association
extern "C" void
parrotcApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  parrotcApp::associateSyntax(syntax, action_factory);
}
void
parrotcApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
