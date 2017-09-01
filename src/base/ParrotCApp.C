#include "ParrotCApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "FreqPoroelasticityFast.h"
#include "FreqPoroelasticSphere.h"
#include "DirichletBCC.h"


template<>
InputParameters validParams<ParrotCApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ParrotCApp::ParrotCApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ParrotCApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ParrotCApp::associateSyntax(_syntax, _action_factory);
}

ParrotCApp::~ParrotCApp()
{
}

// External entry point for dynamic application loading
extern "C" void ParrotCApp__registerApps() { ParrotCApp::registerApps(); }
void
ParrotCApp::registerApps()
{
  registerApp(ParrotCApp);
}

// External entry point for dynamic object registration
extern "C" void ParrotCApp__registerObjects(Factory & factory) { ParrotCApp::registerObjects(factory); }
void
ParrotCApp::registerObjects(Factory & factory)
{
    registerKernel(FreqPoroelasticityFast);
    registerMaterial(FreqPoroelasticSphere);
    registerBoundaryCondition(DirichletBCC);

}

// External entry point for dynamic syntax association
extern "C" void ParrotCApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ParrotCApp::associateSyntax(syntax, action_factory); }
void
ParrotCApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
