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

#ifndef CUBEMARKER_H
#define CUBEMARKER_H

#include"Marker.h"

class CubeMarker;

template<>
InputParameters validParams<CubeMarker>();

class CubeMarker : public Marker
{
public:
  CubeMarker(const InputParameters & parameters);
  //virtual ~SphereMarker(){};

  //virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

private:
    Real _x_min,_x_max,_y_min,_y_max,_z_min,_z_max;
    
};

#endif /* SIMPLEMARKER_H */
