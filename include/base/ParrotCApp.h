#ifndef PARROTCAPP_H
#define PARROTCAPP_H

#include "MooseApp.h"

class ParrotCApp;

template<>
InputParameters validParams<ParrotCApp>();

class ParrotCApp : public MooseApp
{
public:
  ParrotCApp(InputParameters parameters);
  virtual ~ParrotCApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* PARROTCAPP_H */
