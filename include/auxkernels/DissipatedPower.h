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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Auxiliary Kernel to visualize the fibers.                    */
/****************************************************************/


#ifndef DISSIPATEDPOWER_H
#define DISSIPATEDPOWER_H

#include "AuxKernel.h"

class DissipatedPower;

template <>
InputParameters validParams<DissipatedPower>();

class DissipatedPower : public AuxKernel
{
public:
  DissipatedPower(const InputParameters &parameters);

protected:
  virtual Number computeValue() override;

  MaterialProperty<Real> const & _kappa;
  MaterialProperty<Real> const & _eta;
  MaterialProperty<Real> const & _omega;
  VariableNumberGradient const & _grad_pres;
  PostprocessorValue     const & _pp_value;

};
#endif
