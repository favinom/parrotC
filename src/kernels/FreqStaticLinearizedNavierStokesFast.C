#include "FreqStaticLinearizedNavierStokesFast.h"

#include "SystemBase.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<FreqStaticLinearizedNavierStokesFast>() {
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "real part of the first  coupled component");
  params.addRequiredCoupledVar("disp_y", "real part of the second coupled component");
  params.addCoupledVar        ("disp_z", "real part of the third  coupled component");

  return params;
}

FreqStaticLinearizedNavierStokesFast::FreqStaticLinearizedNavierStokesFast(const InputParameters & params) :
    Kernel(params),
    _var_disp_x(_sys.getVariable(_tid, params.get<NonlinearVariableName>("variable"))),
    _P2(     _var_disp_x.phi()     ),
    _grad_P2(_var_disp_x.gradPhi() ),
    _dim(_mesh.dimension()),
    _disp_x_var(coupled("disp_x")),
    _disp_y_var(coupled("disp_y")),
    _disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
    _alpha(getMaterialProperty<Number>("alpha_coefficient")),
    _beta(getMaterialProperty<Number>("beta_coefficient")),
    _sigma_real(getMaterialProperty<RealTensorValue>("stress_real")),
    _sigma_imag(getMaterialProperty<RealTensorValue>("stress_imag")),
    _tr_eps(getMaterialProperty<Number>("trace_strain")),
    _imagUnit(0.0,1.0)
{

    if (_dim == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);

    if (_dim == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    
    _V               = new RealTensorValue **[_dim];
    _eps_lin         = new RealTensorValue **[_dim];
    _sigma_lin_real  = new RealTensorValue **[_dim];
    _sigma_lin_imag  = new RealTensorValue **[_dim];
    
    for (int d=0; d<_dim; ++d)
    {
        _V              [d] = new RealTensorValue *[64];
        _eps_lin        [d] = new RealTensorValue *[64];
        _sigma_lin_real [d] = new RealTensorValue *[64];
        _sigma_lin_imag [d] = new RealTensorValue *[64];
    }
    
    for (int d=0; d<_dim; ++d)
        for (int j = 0; j < 64 ; ++j)
        {
            _V              [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
            _eps_lin        [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
            _sigma_lin_real [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
            _sigma_lin_imag [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
        }
    
    // we allocate the local vectors and matrices
    _f_local= new DenseVector<Number> [_dim];
    _Elasticity = new DenseMatrix<Number> *[_dim];
    
    for (int d=0; d<_dim; ++d)
    {
        _Elasticity[d]=new DenseMatrix<Number>[_dim];
    }

}

void FreqStaticLinearizedNavierStokesFast::computeResidual()
{
    if (_dim==2)
        computeResidual2D();
    if (_dim==3)
        computeResidual3D();

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////
//////////////////////                 RESIDUAL 2D
//////////////////////
//////////////////////////////////////////////////////////////////////////////////


void FreqStaticLinearizedNavierStokesFast::computeResidual2D()
{
    // Here we take the reference to the residual entries
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    
    // Resize and initialization of local vectors
    // NOTE we assume the displacement variables are of the same order
    for (int i=0; i<_dim; ++i)
    {
        _f_local[i].resize(f_x.size());
        _f_local[i].zero();
    }
    
    // Assembling of local RHSs
    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _P2.size(); _i++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                NumberTensorValue Sigma =_sigma_real[_qp];
                Sigma+=_imagUnit*_sigma_imag[_qp];
                
                _f_local[d](_i) += _JxW[_qp] * _coord[_qp] *
                (
                 Sigma(d , 0) * _grad_P2[_i][_qp] (0)+
                 Sigma(d , 1) * _grad_P2[_i][_qp] (1)
                 );
            }
    
    // Some of local contributions of RHS to global RHSs
    f_x += _f_local[0];
    f_y += _f_local[1];
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////
//////////////////////                 RESIDUAL 3D
//////////////////////
//////////////////////////////////////////////////////////////////////////////////

void FreqStaticLinearizedNavierStokesFast::computeResidual3D()
{
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    DenseVector<Number> & f_z = _assembly.residualBlock(_disp_z_var);
    
    for (int i=0; i<_dim; ++i)
    {
        _f_local[i].resize(f_x.size());
        _f_local[i].zero();
    }

    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _P2.size(); _i++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                NumberTensorValue Sigma=_sigma_real[_qp];
                Sigma+=_imagUnit*_sigma_imag[_qp];

                _f_local[d](_i) += _JxW[_qp] * _coord[_qp] *
                (
                 Sigma(d , 0) * _grad_P2[_i][_qp] (0)+
                 Sigma(d , 1) * _grad_P2[_i][_qp] (1)+
                 Sigma(d , 2) * _grad_P2[_i][_qp] (2)
                 );
                
            }
    
    f_x += _f_local[0];
    f_y += _f_local[1];
    f_z += _f_local[2];
    
}

void FreqStaticLinearizedNavierStokesFast::computeOffDiagJacobian(unsigned int jvar)
{
    if ( jvar==_var.number() )
    {
        if (_dim==2)
        {
            computeJacobian2D();
        }
        else
        {
            if (_dim==3)
            {
                computeJacobian3D();
            }
            else
            {
                mooseError("FreqStaticLinearizedNavierStokesFast::computeOffDiagJacobian");
            }
        }
    }
}

void FreqStaticLinearizedNavierStokesFast::computeJacobian2D()
{
 
    // References to the right components of the stiffness matrix
    DenseMatrix<Number> & A_xx = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
    DenseMatrix<Number> & A_xy = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
    
    DenseMatrix<Number> & A_yx = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
    DenseMatrix<Number> & A_yy = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
    
    for (int i=0; i<_dim; ++i)
    {
        for (int j=0; j<_dim; ++j)
        {
            _Elasticity[i][j].resize(A_xx.m(),A_xx.n());
            _Elasticity[i][j].zero();
        }
    }
    
    initTensorVariables();
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
            for (_i = 0; _i < _P2.size(); ++_i)
                for (_j = 0; _j < _P2.size(); ++_j)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                    {
                        RealTensorValue const & Sigma_real=_sigma_lin_real[jdim][_j][_qp];
                        RealTensorValue const & Sigma_imag=_sigma_lin_imag[jdim][_j][_qp];
                        RealTensorValue & V=_V[idim][_i][_qp];
                        _Elasticity[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*(
			Sigma_real.contract(V)+_imagUnit*Sigma_imag.contract(V));
                    }
    
    A_xx  += _Elasticity[0][0];
    A_xy  += _Elasticity[0][1];
    A_yx  += _Elasticity[1][0];
    A_yy  += _Elasticity[1][1];
    
}

void FreqStaticLinearizedNavierStokesFast::computeJacobian3D()
{
    // References to the right components of the stiffness matrix
    DenseMatrix<Number> &  A_xx = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
    DenseMatrix<Number> &  A_xy = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
    DenseMatrix<Number> &  A_xz = _assembly.jacobianBlock(_disp_x_var, _disp_z_var);
    
    DenseMatrix<Number> &  A_yx = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
    DenseMatrix<Number> &  A_yy = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
    DenseMatrix<Number> &  A_yz = _assembly.jacobianBlock(_disp_y_var, _disp_z_var);

    DenseMatrix<Number> &  A_zx = _assembly.jacobianBlock(_disp_z_var, _disp_x_var);
    DenseMatrix<Number> &  A_zy = _assembly.jacobianBlock(_disp_z_var, _disp_y_var);
    DenseMatrix<Number> &  A_zz = _assembly.jacobianBlock(_disp_z_var, _disp_z_var);
    
    for (int i=0; i<_dim; ++i)
    {
        for (int j=0; j<_dim; ++j)
        {
            _Elasticity[i][j].resize(A_xx.m(),A_xx.n());
            _Elasticity[i][j].zero();
        }
    }
    
    initTensorVariables();
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
            for (_i = 0; _i < _P2.size(); ++_i)
                for (_j = 0; _j < _P2.size(); ++_j)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                    {
                        RealTensorValue const & Sigma_real=_sigma_lin_real[jdim][_j][_qp];
                        RealTensorValue const & Sigma_imag=_sigma_lin_imag[jdim][_j][_qp];
                        RealTensorValue & V=_V[idim][_i][_qp];
                        _Elasticity[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*(
			Sigma_real.contract(V)+_imagUnit*Sigma_imag.contract(V));
                    }
    
    A_xx  += _Elasticity[0][0];
    A_xy  += _Elasticity[0][1];
    A_xz  += _Elasticity[0][2];
    
    A_yx  += _Elasticity[1][0];
    A_yy  += _Elasticity[1][1];
    A_yz  += _Elasticity[1][2];

    A_zx  += _Elasticity[2][0];
    A_zy  += _Elasticity[2][1];
    A_zz  += _Elasticity[2][2];
    
}


void FreqStaticLinearizedNavierStokesFast::initTensorVariables()
{
        for (int d=0; d<_mesh.dimension(); ++d)
            for ( _j = 0; _j < _P2.size() ; ++_j)
                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                {
                    _V[d][_j][_qp]=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                    for (int j=0; j<_mesh.dimension(); ++j)
                        _V[d][_j][_qp](d,j)=_grad_P2[_j][_qp](j);
                    
                    RealTensorValue & V=_V[d][_j][_qp];
                    _eps_lin[d][_j][_qp]=0.5*( V+V.transpose() );
                    
                    RealTensorValue & eps=_eps_lin[d][_j][_qp];
                    _sigma_lin_real[d][_j][_qp]=2.0 * _alpha[_qp].real()*eps + _beta[_qp].real()*eps.tr()*_identity;
                    _sigma_lin_imag[d][_j][_qp]=2.0 * _alpha[_qp].imag()*eps + _beta[_qp].imag()*eps.tr()*_identity;
                }
    
}
