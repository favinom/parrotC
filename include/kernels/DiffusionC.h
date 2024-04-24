#ifndef DIFFUSIONC_H
#define DIFFUSIONC_H

#include "Kernel.h"

class DiffusionC;

template<>
InputParameters validParams<DiffusionC>();


class DiffusionC : public Kernel
{
public:
    DiffusionC( InputParameters const & params);
protected:
    virtual Number computeQpResidual();
    virtual Number computeQpJacobian();
    //virtual Number computeQpOffDiagJacobian(unsigned int jvar){};

    Number _imagUnit;
    Real _omega;
    
};

#endif 
