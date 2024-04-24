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

#ifndef SPHEREMARKER_H
#define SPHEREMARKER_H

#include"Marker.h"

class SphereMarker;

template<>
InputParameters validParams<SphereMarker>();

class SphereMarker : public Marker
{
public:
  SphereMarker(const InputParameters & parameters);
  //virtual ~SphereMarker(){};

  //virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

private:
    Real _x_center,_y_center,_z_center,_radius;
};

#endif /* SIMPLEMARKER_H */
