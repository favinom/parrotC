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

#include "CubeMarker.h"

// libMesh includes
#include "libmesh/error_vector.h"

template<>
InputParameters validParams<CubeMarker>()
{
  InputParameters params = validParams<Marker>();
  params.addRequiredParam<Real>("x_min", "x coordinate of the center");
  params.addRequiredParam<Real>("x_max", "x coordinate of the center");
  params.addRequiredParam<Real>("y_min", "x coordinate of the center");
  params.addRequiredParam<Real>("y_max", "x coordinate of the center");

  return params;
}


CubeMarker::CubeMarker(const InputParameters & parameters) :
Marker(parameters),
_x_min(getParam<Real>("x_min")),
_x_max(getParam<Real>("x_max")),
_y_min(getParam<Real>("y_min")),
_y_max(getParam<Real>("y_max"))
{}


Marker::MarkerValue
CubeMarker::computeElementMarker()
{
    Point _centroid = (*_current_elem).centroid();
    Real x_center=_centroid(0);
    Real y_center=_centroid(1);

    if (_x_min<x_center && x_center<_x_max)
    {
        Real dist1=std::fabs(y_center-_y_min);
        Real dist2=std::fabs(y_center-_y_max);
        Real dist=std::min(dist1,dist2);
        
        if (dist<1.25*(*_current_elem).hmax())
            return REFINE;
    }

    if (_y_min<y_center && y_center<_y_max)
    {
        Real dist1=std::fabs(x_center-_x_min);
        Real dist2=std::fabs(x_center-_x_max);
        Real dist=std::min(dist1,dist2);
        
        if (dist<1.25*(*_current_elem).hmax())
            return REFINE;
    }

    return DO_NOTHING;
  
}

