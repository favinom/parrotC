#include "DiffusionC.h"

//#include "SystemBase.h"
//#include "MooseMesh.h"
//#include "MooseVariable.h"
//#include "Assembly.h"
//#include "libmesh/quadrature.h"

template <>
InputParameters validParams<DiffusionC>() {
  InputParameters params = validParams<Kernel>();
  return params;
}

DiffusionC::DiffusionC(const InputParameters & params) :
Kernel(params),
_imagUnit(0.0,1.0)
{}

Number DiffusionC::computeQpResidual()
{
    _omega=std::pow(10.0,_t);
    return 1.0/_omega*_grad_u[_qp]*_grad_test[_i][_qp]+_imagUnit*_u[_qp]*_test[_i][_qp];
}

Number DiffusionC::computeQpJacobian()
{
    _omega=std::pow(10.0,_t);
    return 1.0/_omega*_grad_phi[_j][_qp]*_grad_test[_i][_qp]+_imagUnit*_phi[_j][_qp]*_test[_i][_qp];
}
