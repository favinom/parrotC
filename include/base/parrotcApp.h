#ifndef PARROTCAPP_H
#define PARROTCAPP_H

#include "MooseApp.h"

class parrotcApp;

template <>
InputParameters validParams<parrotcApp>();

class parrotcApp : public MooseApp
{
public:
  parrotcApp(InputParameters parameters);
  virtual ~parrotcApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* PARROTCAPP_H */
