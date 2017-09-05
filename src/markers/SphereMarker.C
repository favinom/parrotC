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

#include "SphereMarker.h"

// libMesh includes
#include "libmesh/error_vector.h"

template<>
InputParameters validParams<SphereMarker>()
{
  InputParameters params = validParams<Marker>();
  params.addRequiredParam<Real>("x_center", "x coordinate of the center");
  params.addRequiredParam<Real>("y_center", "y coordinate of the center");
  params.addRequiredParam<Real>("z_center", "z coordinate of the center");
  params.addRequiredParam<Real>("radius", "radius of the sphere");
  return params;
}


SphereMarker::SphereMarker(const InputParameters & parameters) :
    Marker(parameters),
    _x_center(getParam<Real>("x_center")),
    _y_center(getParam<Real>("y_center")),
    _z_center(getParam<Real>("z_center")),
    _radius(getParam<Real>("radius"))
{
}


Marker::MarkerValue
SphereMarker::computeElementMarker()
{
    Real min_distance= 1e15;
    Real max_distance=-1e15;
    Real x_node,y_node,z_node;
    Real x_dist,y_dist,z_dist,distance;
    
    for (int i = 0; i < (*_current_elem).n_nodes(); ++i)
    {
        x_node=(*_current_elem).point(i)(0);
        y_node=(*_current_elem).point(i)(1);
        z_node=(*_current_elem).point(i)(2);
        x_dist=x_node-_x_center;
        y_dist=y_node-_y_center;
        z_dist=z_node-_z_center;
        
        distance=std::sqrt( x_dist*x_dist+ y_dist*y_dist + z_dist*z_dist );
        
        min_distance=std::min(distance,min_distance);
        max_distance=std::max(distance,max_distance);
    }
    
    if (min_distance<=_radius && _radius <= max_distance)
    {
//        std::cout<<"ciao"<<std::endl;
          return REFINE;
    }
    return DO_NOTHING;
  
}

