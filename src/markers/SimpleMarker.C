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

#include "SimpleMarker.h"

// libMesh includes
#include "libmesh/error_vector.h"

template<>
InputParameters validParams<SimpleMarker>()
{
  InputParameters params = validParams<Marker>();
  params.addRequiredParam<int>("fn", "number of fractures");
  params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
  params.addRequiredParam<std::string>("fl_string", "length of fractures");
  params.addRequiredParam<std::string>("ft_string", "thickness of fractures");
  params.addRequiredParam<std::string>("fa_string", "angle of fractures");

  return params;
}


SimpleMarker::SimpleMarker(const InputParameters & parameters) :
    Marker(parameters),
    _fn(getParam<int>("fn")),
    _fx_string(getParam<std::string>("fx_string")),
    _fy_string(getParam<std::string>("fy_string")),
    _fl_string(getParam<std::string>("fl_string")),
    _ft_string(getParam<std::string>("ft_string")),
    _fa_string(getParam<std::string>("fa_string"))

{
}

void
SimpleMarker::markerSetup()
{
}

Marker::MarkerValue
SimpleMarker::computeElementMarker()
{
  Real lbx, lby, ubx, uby;
  Real brcx, brcy, blcx, blcy;
  Real trcx, trcy, tlcx, tlcy;
  Real w1, w2, test_lb_y_right, test_lb_y_left, test_ub_y_right, test_ub_y_left;
  Real maxnum, minnum;
  // Amount of nodes of element
  //Real max_z = (*_current_elem).point(0)(2);
  //Real min_z = (*_current_elem).point(0)(2);
  Real max_y = (*_current_elem).point(0)(1);
  Real min_y = (*_current_elem).point(0)(1);
  Real max_x = (*_current_elem).point(0)(0);
  Real min_x = (*_current_elem).point(0)(0);
  for (int i = 1; i < (*_current_elem).n_nodes(); i++) 
  {
    //if ((*_current_elem).point(i)(2)>max_z)
    //{
    //  max_z = (*_current_elem).point(i)(2);
    //}
    //if ((*_current_elem).point(i)(2)<min_z)
    //{
    //  min_z = (*_current_elem).point(i)(2);
    //}
    if ((*_current_elem).point(i)(1)>max_y)
    {
      max_y = (*_current_elem).point(i)(1);
    }
    if ((*_current_elem).point(i)(1)<min_y)
    {
      min_y = (*_current_elem).point(i)(1);
    }
    if ((*_current_elem).point(i)(0)>max_x)
    {
      max_x = (*_current_elem).point(i)(0);
    }
    if ((*_current_elem).point(i)(0)<min_x)
    {
      min_x = (*_current_elem).point(i)(0);
    }
  }
  //int npoints = 2;
  //Real fx[npoints] = {0.0,0.0}; // x-coordinate of fracture center point
  //Real fy[npoints] = {0.0,0.0}; // y-coordinate of fracture center point
  //Real fl[npoints] = {200.0,200.0}; // fracture length
  //Real ft[npoints] = {0.1,0.1}; // fracture thickness
  //Real fa[npoints] = {45.0,135.0}; // fracture dip angle
  Real fx[_fn];
  Real fy[_fn];
  Real fl[_fn];
  Real ft[_fn];
  Real fa[_fn];
  int commavec[_fn+1];
  // Convert the fx-string to double
  // First, determine where the commas are located
  int i = 1;
  int si = 1;
  int stopcon = 0;
  commavec[0] = -1;
  while (stopcon == 0) {
    if (_fx_string[si] == ',') {
      commavec[i] = si;
      i+=1;
    } else if (_fx_string[si] == '\0') {
      commavec[i] = si;
      stopcon = 1;
    }
    si+=1;
  }
  // Second, do the actual conversion
  for (int i = 0; i < _fn; i++) {
    fx[i] = atof(_fx_string.substr(commavec[i]+1,commavec[i+1]-commavec[i]-1).c_str());
  }  
  // Convert the fy-string to double
  // First, determine where the commas are located
  i = 1;
  si = 1;
  stopcon = 0;
  commavec[0] = -1;
  while (stopcon == 0) {
    if (_fy_string[si] == ',') {
      commavec[i] = si;
      i+=1;
    } else if (_fy_string[si] == '\0') {
      commavec[i] = si;
      stopcon = 1;
    }
    si+=1;
  }
  // Second, do the actual conversion
  for (int i = 0; i < _fn; i++) {
    fy[i] = atof(_fy_string.substr(commavec[i]+1,commavec[i+1]-commavec[i]-1).c_str());
  }  
  // Convert the fl-string to double
  // First, determine where the commas are located
  i = 1;
  si = 1;
  stopcon = 0;
  commavec[0] = -1;
  while (stopcon == 0) {
    if (_fl_string[si] == ',') {
      commavec[i] = si;
      i+=1;
    } else if (_fl_string[si] == '\0') {
      commavec[i] = si;
      stopcon = 1;
    }
    si+=1;
  }
  // Second, do the actual conversion
  for (int i = 0; i < _fn; i++) {
    fl[i] = atof(_fl_string.substr(commavec[i]+1,commavec[i+1]-commavec[i]-1).c_str());
  }  
  // Convert the ft-string to double
  // First, determine where the commas are located
  i = 1;
  si = 1;
  stopcon = 0;
  commavec[0] = -1;
  while (stopcon == 0) {
    if (_ft_string[si] == ',') {
      commavec[i] = si;
      i+=1;
    } else if (_ft_string[si] == '\0') {
      commavec[i] = si;
      stopcon = 1;
    }
    si+=1;
  }
  // Second, do the actual conversion
  for (int i = 0; i < _fn; i++) {
    ft[i] = atof(_ft_string.substr(commavec[i]+1,commavec[i+1]-commavec[i]-1).c_str());
  }  
  // Convert the fa-string to double
  // First, determine where the commas are located
  i = 1;
  si = 1;
  stopcon = 0;
  commavec[0] = -1;
  while (stopcon == 0) {
    if (_fa_string[si] == ',') {
      commavec[i] = si;
      i+=1;
    } else if (_fa_string[si] == '\0') {
      commavec[i] = si;
      stopcon = 1;
    }
    si+=1;
  }
  // Second, do the actual conversion
  for (int i = 0; i < _fn; i++) {
    fa[i] = atof(_fa_string.substr(commavec[i]+1,commavec[i+1]-commavec[i]-1).c_str());
  }  
  //std::cout << "markers: fa[0] = " << fa[0] << "; fa[1] = " << fa[1] << std::endl;
  for (int i = 0; i < _fn; i++) {
      // Determine the corner points of the fracture
      // Determine center point of lower boundary (lb)
      if (fa[i]<=90.0) {
          lbx = fx[i]+sin((90.0-fa[i])/180.0*pi)*ft[i]/2.0;
          lby = fy[i]-cos((90.0-fa[i])/180.0*pi)*ft[i]/2.0;
      } else {
          lbx = fx[i]-sin((fa[i]-90.0)/180.0*pi)*ft[i]/2.0;
          lby = fy[i]-cos((fa[i]-90.0)/180.0*pi)*ft[i]/2.0;
      }
      // Determine center point of upper boundary (ub)
      if (fa[i]<=90.0) {
          ubx = fx[i]-sin((90.0-fa[i])/180.0*pi)*ft[i]/2.0;
          uby = fy[i]+cos((90.0-fa[i])/180.0*pi)*ft[i]/2.0;
      } else {
          ubx = fx[i]+sin((fa[i]-90.0)/180.0*pi)*ft[i]/2.0;
          uby = fy[i]+cos((fa[i]-90.0)/180.0*pi)*ft[i]/2.0;
      }
      // Determine bottom right corner (brc)
      if (fa[i]<=90.0) {
          brcx = lbx+sin(fa[i]/180.0*pi)*fl[i]/2.0;
          brcy = lby+cos(fa[i]/180.0*pi)*fl[i]/2.0;
      } else {
          brcx = lbx+sin((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
          brcy = lby-cos((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
      }
      // Determine bottom left corner (blc)
      if (fa[i]<=90.0) {
          blcx = lbx-sin(fa[i]/180.0*pi)*fl[i]/2.0;
          blcy = lby-cos(fa[i]/180.0*pi)*fl[i]/2.0;
      } else {
          blcx = lbx-sin((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
          blcy = lby+cos((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
      }
      // Determine top right corner (trc)
      if (fa[i]<=90.0) {
          trcx = ubx+sin(fa[i]/180.0*pi)*fl[i]/2.0;
          trcy = uby+cos(fa[i]/180.0*pi)*fl[i]/2.0;
      } else {
          trcx = ubx+sin((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
          trcy = uby-cos((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
      }
      // Determine top left corner (tlc)
      if (fa[i]<=90.0) {
          tlcx = ubx-sin(fa[i]/180.0*pi)*fl[i]/2.0;
          tlcy = uby-cos(fa[i]/180.0*pi)*fl[i]/2.0;
      } else {
          tlcx = ubx-sin((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
          tlcy = uby+cos((180.0-fa[i])/180.0*pi)*fl[i]/2.0;
      }
      // Determine the largest x-value of the fracture
      maxnum = brcx;
      if (maxnum < blcx) maxnum = blcx;
      if (maxnum < trcx) maxnum = trcx;
      if (maxnum < tlcx) maxnum = tlcx;
      // Determine the smallest x-value of the fracture 
      minnum = brcx;
      if (minnum > blcx) minnum = blcx;
      if (minnum > trcx) minnum = trcx;
      if (minnum > tlcx) minnum = tlcx;
      // Check if the fracture is in the x-dimension in the range of the current element 
      if (maxnum>min_x && minnum<max_x){
        // Now, check also the y-dimension
        if (fa[i]<0.1 || fa[i]>189.9) {
          // Determine the largest y-value of the fracture
          maxnum = brcy;
          if (maxnum < blcy) maxnum = blcy;
          if (maxnum < trcy) maxnum = trcy;
          if (maxnum < tlcy) maxnum = tlcy;
          // Determine the smallest y-value of the fracture 
          minnum = brcy;
          if (minnum > blcy) minnum = blcy;
          if (minnum > trcy) minnum = trcy;
          if (minnum > tlcy) minnum = tlcy;
          if (maxnum>min_y && minnum<max_y){
              return REFINE;
          }
        } else {
          // Check if the fracture is at the right side of the element inside the element
          w1 = (max_x-blcx)/(brcx-blcx);
          w2 = (brcx-max_x)/(brcx-blcx);
          test_lb_y_right = w2*blcy+w1*brcy;
          w1 = (max_x-tlcx)/(trcx-tlcx);
          w2 = (trcx-max_x)/(trcx-tlcx);
          test_ub_y_right = w2*tlcy+w1*trcy;
          if (test_lb_y_right<max_y && test_ub_y_right>min_y) {
              return REFINE;
          }
          // Check if the fracture is at the left side of the element inside the element
          w1 = (min_x-blcx)/(brcx-blcx);
          w2 = (brcx-min_x)/(brcx-blcx);
          test_lb_y_left = w2*blcy+w1*brcy;
          w1 = (min_x-tlcx)/(trcx-tlcx);
          w2 = (trcx-min_x)/(trcx-tlcx);
          test_ub_y_left = w2*tlcy+w1*trcy;
          if (test_lb_y_left<max_y && test_ub_y_left>min_y) {
              return REFINE;
          }
          // Check if the fracture is so steep that it crosses the element on the top and the bottom
          // for an increasing fracture
          if (test_ub_y_left<min_y && test_lb_y_right>max_y) {
              return REFINE;
          }
          // for an decreasing fracture
          if (test_ub_y_right<max_y && test_lb_y_left>min_y) {
              return REFINE;
          }
        }
      }
  }
  return DO_NOTHING;
  
}

