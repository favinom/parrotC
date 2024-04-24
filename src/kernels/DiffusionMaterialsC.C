#include "DiffusionMaterialsC.h"

//#include "SystemBase.h"
//#include "MooseMesh.h"
//#include "MooseVariable.h"
//#include "Assembly.h"
#include "libmesh/quadrature.h"
//#include "libmesh/mesh_tools.h"

template <>
InputParameters validParams<DiffusionMaterialsC>() {
  InputParameters params = validParams<Kernel>();
  return params;
}

DiffusionMaterialsC::DiffusionMaterialsC(const InputParameters & params) :
Kernel(params),
_diffusion(getMaterialProperty<Real>("diffusion_property")),
_inv_m(getMaterialProperty<Real>("inverse_of_m")),
_imagUnit(0.0,1.0)
{}

Number DiffusionMaterialsC::computeQpResidual()
{
//    for (int i=0; i<4; ++i)
//    {
//        Node const & node=_current_elem->node_ref(i);
//        std::cout<<node<<std::endl;
//    }
//    for (int myqp=0; myqp<_qrule->n_points(); ++myqp)
//    {
//        std::cout<<_q_point[myqp]<<" "<<_JxW[myqp]<<std::endl;
//    }
//    exit(1);
    
    _omega=std::pow(10.0,_t);
    return _diffusion[_qp]*_grad_u[_qp]*_grad_test[_i][_qp]+_imagUnit*_inv_m[_qp]*_u[_qp]*_test[_i][_qp];
}

Number DiffusionMaterialsC::computeQpJacobian()
{
    _omega=std::pow(10.0,_t);
    return _diffusion[_qp]*_grad_phi[_j][_qp]*_grad_test[_i][_qp]+_imagUnit*_inv_m[_qp]*_phi[_j][_qp]*_test[_i][_qp];
}
