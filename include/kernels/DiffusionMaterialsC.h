#ifndef DIFFUSIONMATERIALSC_H
#define DIFFUSIONMATERIALSC_H

#include "Kernel.h"

class DiffusionMaterialsC;

template<>
InputParameters validParams<DiffusionMaterialsC>();


class DiffusionMaterialsC : public Kernel
{
public:
    DiffusionMaterialsC( InputParameters const & params);
protected:
    virtual Number computeQpResidual();
    virtual Number computeQpJacobian();
    //virtual Number computeQpOffDiagJacobian(unsigned int jvar){};

    const MaterialProperty<Real> & _diffusion;
    const MaterialProperty<Real> & _inv_m;
    
    Number _imagUnit;
    Real _omega;
    
};

#endif 
