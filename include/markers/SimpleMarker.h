/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef SIMPLEMARKER_H
#define SIMPLEMARKER_H

#include"Marker.h"

class SimpleMarker;

template<>
InputParameters validParams<SimpleMarker>();

class SimpleMarker : public Marker
{
public:
  SimpleMarker(const InputParameters & parameters);
  virtual ~SimpleMarker(){};

  virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

private:
  int _fn;
  std::string _fx_string;
  std::string _fy_string;
  std::string _fl_string;
  std::string _ft_string;
  std::string _fa_string;

};

#endif /* SIMPLEMARKER_H */
