#include "FreqPoroelasticityFast.h"

#include "SystemBase.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<FreqPoroelasticityFast>() {
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "real part of the first  coupled component");
  params.addRequiredCoupledVar("disp_y", "real part of the second coupled component");
  params.addCoupledVar        ("disp_z", "real part of the third  coupled component");

  params.addRequiredCoupledVar("pres"  , "real part of pressure");

  params.addRequiredParam<NonlinearVariableName>("variable_p", "THE pressure variable");

  return params;
}

FreqPoroelasticityFast::FreqPoroelasticityFast(const InputParameters & params) :
    Kernel(params),
    _var_disp_x(_sys.getVariable(_tid, params.get<NonlinearVariableName>("variable"))),
    _var_pres  (_sys.getVariable(_tid, params.get<NonlinearVariableName>("variable_p"))),
    _P2(     _var_disp_x.phi()     ),
    _grad_P2(_var_disp_x.gradPhi() ),
    _P1(     _var_pres.phi()       ),
    _grad_P1(_var_pres.gradPhi()   ),
    _dim(_mesh.dimension()),
    _disp_x_var(coupled("disp_x")),
    _disp_y_var(coupled("disp_y")),
    _disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
    _p_var(coupled("pres")),
    _pres(coupledValue("pres")),
    _grad_pres(coupledGradient("pres")),
    _mu(getMaterialProperty<Real>("mu_property")),
    _lambda(getMaterialProperty<Real>("lambda_property")),
    _alpha(getMaterialProperty<Real>("alpha_property")),
    _sigma_real(getMaterialProperty<RealTensorValue>("stress_real")),
    _sigma_imag(getMaterialProperty<RealTensorValue>("stress_imag")),
    _diffusion(getMaterialProperty<Real>("diffusion_property")),
    _inv_m(getMaterialProperty<Real>("inverse_of_m")),
    _tr_eps(getMaterialProperty<Number>("trace_strain")),
    _imagUnit(0.0,1.0)
{

    if (_dim == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);

    if (_dim == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    
    _V          = new RealTensorValue **[_dim];
    _eps_lin    = new RealTensorValue **[_dim];
    _sigma_lin  = new RealTensorValue **[_dim];
    
    for (int d=0; d<_dim; ++d)
    {
        _V         [d] = new RealTensorValue *[64];
        _eps_lin   [d] = new RealTensorValue *[64];
        _sigma_lin [d] = new RealTensorValue *[64];
    }
    
    for (int d=0; d<_dim; ++d)
        for (int j = 0; j < 64 ; ++j)
        {
            _V         [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
            _eps_lin   [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
            _sigma_lin [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
        }
    
    // we allocate the local vecotrs and matrices
    _f_local= new DenseVector<Number> [_dim];
    _Elasticity = new DenseMatrix<Number> *[_dim];
    
    for (int d=0; d<_dim; ++d)
    {
        _Elasticity[d]=new DenseMatrix<Number>[_dim];
    }

    _Gradient= new DenseMatrix<Number> [_dim];
    _Divergence=new DenseMatrix<Number> [_dim];

    
    
}

void FreqPoroelasticityFast::computeResidual()
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


void FreqPoroelasticityFast::computeResidual2D()
{
    // Here we take the reference to the residual entries
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    DenseVector<Number> & f_p = _assembly.residualBlock(_p_var);
    
    // Resize and initialization of local vectors
    // NOTE we assume the displacement variables are of the same order
    for (int i=0; i<_dim; ++i)
    {
        _f_local[i].resize(f_x.size());
        _f_local[i].zero();
    }
    _f_p_local.resize(f_p.size());
    _f_p_local.zero();
    
    // Assembling of local RHSs
    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _P2.size(); _i++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                NumberTensorValue P=_sigma_real[_qp];
                P+=_imagUnit*_sigma_imag[_qp];
                P-=_alpha[_qp]*_pres[_qp]*_identity;
                
                _f_local[d](_i) += _JxW[_qp] * _coord[_qp] *
                (
                 P(d , 0) * _grad_P2[_i][_qp] (0)+
                 P(d , 1) * _grad_P2[_i][_qp] (1)
                 );
            }
    
    for (_i=0; _i<_P1.size(); ++_i)
        for (_qp=0; _qp<_qrule->n_points(); ++_qp)
        {
            _f_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_imagUnit*(_alpha[_qp]*_tr_eps[_qp]+_inv_m[_qp]*_pres[_qp])*_P1[_i][_qp];
            _f_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_diffusion[_qp]*_grad_pres[_qp]*_grad_P1[_i][_qp];
        }
    
    // Some of local contributions of RHS to global RHSs
    f_x += _f_local[0];
    f_y += _f_local[1];
    f_p += _f_p_local;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////
//////////////////////                 RESIDUAL 3D
//////////////////////
//////////////////////////////////////////////////////////////////////////////////

void FreqPoroelasticityFast::computeResidual3D()
{
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    DenseVector<Number> & f_z = _assembly.residualBlock(_disp_z_var);
    
    DenseVector<Number> & f_p = _assembly.residualBlock(_p_var);

    for (int i=0; i<_dim; ++i)
    {
        _f_local[i].resize(f_x.size());
        _f_local[i].zero();
    }
    _f_p_local.resize(f_p.size());
    _f_p_local.zero();

    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _P2.size(); _i++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                NumberTensorValue P=_sigma_real[_qp];
                P+=_imagUnit*_sigma_imag[_qp];
                P-=_alpha[_qp]*_pres[_qp]*_identity;

                _f_local[d](_i) += _JxW[_qp] * _coord[_qp] *
                (
                 P(d , 0) * _grad_P2[_i][_qp] (0)+
                 P(d , 1) * _grad_P2[_i][_qp] (1)+
                 P(d , 2) * _grad_P2[_i][_qp] (2)
                 );
                
            }
    
    for (_i=0; _i<_P1.size(); ++_i)
        for (_qp=0; _qp<_qrule->n_points(); ++_qp)
        {
            _f_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_imagUnit*(_alpha[_qp]*_tr_eps[_qp]+_inv_m[_qp]*_pres[_qp])*_P1[_i][_qp];
            _f_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_diffusion[_qp]*_grad_pres[_qp]*_grad_P1[_i][_qp];
        }
    
    f_x += _f_local[0];
    f_y += _f_local[1];
    f_z += _f_local[2];
    f_p += _f_p_local;
    
}

void FreqPoroelasticityFast::computeOffDiagJacobian(unsigned int jvar)
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
                mooseError("FreqPoroelasticityFast::computeOffDiagJacobian");
            }
        }
    }
}

void FreqPoroelasticityFast::computeJacobian2D()
{
 
    // References to the right components of the stiffness matrix
    DenseMatrix<Number> & A_xx = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
    DenseMatrix<Number> & A_xy = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
    DenseMatrix<Number> & BT_xp = _assembly.jacobianBlock(_disp_x_var, _p_var);
    
    DenseMatrix<Number> & A_yx = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
    DenseMatrix<Number> & A_yy = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
    DenseMatrix<Number> & BT_yp = _assembly.jacobianBlock(_disp_y_var, _p_var);

    DenseMatrix<Number> & B_px = _assembly.jacobianBlock(_p_var,_disp_x_var);
    DenseMatrix<Number> & B_py = _assembly.jacobianBlock(_p_var,_disp_y_var);
    DenseMatrix<Number> & C_pp = _assembly.jacobianBlock(_p_var,_p_var);
    
    for (int i=0; i<_dim; ++i)
    {
        for (int j=0; j<_dim; ++j)
        {
            _Elasticity[i][j].resize(A_xx.m(),A_xx.n());
            _Elasticity[i][j].zero();
        }
        _Gradient[i].resize(BT_xp.m(),BT_xp.n());
        _Gradient[i].zero();
    }
    _Diffusion.resize(C_pp.m(),C_pp.n());
    _Diffusion.zero();
    _Mass.resize(C_pp.m(),C_pp.n());
    _Mass.zero();
    
    initTensorVariables();
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=idim; jdim<_dim; ++jdim)
            for (_i = 0; _i < _P2.size(); ++_i)
                for (_j = 0; _j < _P2.size(); ++_j)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                    {
                        RealTensorValue & sigma=_sigma_lin[jdim][_j][_qp];
                        RealTensorValue & V=_V[idim][_i][_qp];
                        _Elasticity[idim][jdim](_i,_j)+=_JxW[_qp]*sigma.contract(V); // * _coord[_qp]
                    }
    
    for (int idim=0; idim<_dim; ++idim)
        for (_i = 0; _i < _P2.size(); ++_i)
            for (_j = 0; _j < _P1.size(); ++_j)
                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                {
                    RealTensorValue & V=_V[idim][_i][_qp];
                    _Gradient[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*_alpha[_qp]*_P1[_j][_qp]*V.tr();
                }
    
    _Gradient[0].get_transpose(_Divergence[0]);
    _Gradient[1].get_transpose(_Divergence[1]);
    _Divergence[0]*=_imagUnit;
    _Divergence[1]*=_imagUnit;
    
    for (int i=0; i<_Elasticity[1][0].m(); ++i)
        for (int j=0; j<_Elasticity[1][0].m(); ++j)
            _Elasticity[1][0](i,j)=_Elasticity[0][1](j,i);
    
    for (_i=0;_i<_P1.size(); ++_i)
        for (_j=0;_j<_P1.size(); ++_j)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                _Diffusion(_i,_j)+=_JxW[_qp]*_diffusion[_qp]*_grad_P1[_i][_qp]*_grad_P1[_j][_qp]; // * _coord[_qp]

    for (_i=0;_i<_P1.size(); ++_i)
        for (_j=0;_j<_P1.size(); ++_j)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                _Mass(_i,_j)+=_JxW[_qp] * _coord[_qp]*_inv_m[_qp]*_P1[_i][_qp]*_P1[_j][_qp];

  
    A_xx  += _Elasticity[0][0];
    A_xy  += _Elasticity[0][1];
    A_yx  += _Elasticity[1][0];
    A_yy  += _Elasticity[1][1];
    
    BT_xp -= _Gradient[0];
    BT_yp -= _Gradient[1];

    B_px  -= _Divergence[0];
    B_py  -= _Divergence[1];
    
    C_pp  -= _Diffusion;
    _Mass *= -1.0*_imagUnit;
    C_pp  +=_Mass;
    
}

void FreqPoroelasticityFast::computeJacobian3D()
{
    // References to the right components of the stiffness matrix
    DenseMatrix<Number> &  A_xx = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
    DenseMatrix<Number> &  A_xy = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
    DenseMatrix<Number> &  A_xz = _assembly.jacobianBlock(_disp_x_var, _disp_z_var);
    DenseMatrix<Number> & BT_xp = _assembly.jacobianBlock(_disp_x_var,      _p_var);
    
    DenseMatrix<Number> &  A_yx = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
    DenseMatrix<Number> &  A_yy = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
    DenseMatrix<Number> &  A_yz = _assembly.jacobianBlock(_disp_y_var, _disp_z_var);
    DenseMatrix<Number> & BT_yp = _assembly.jacobianBlock(_disp_y_var,      _p_var);

    DenseMatrix<Number> &  A_zx = _assembly.jacobianBlock(_disp_z_var, _disp_x_var);
    DenseMatrix<Number> &  A_zy = _assembly.jacobianBlock(_disp_z_var, _disp_y_var);
    DenseMatrix<Number> &  A_zz = _assembly.jacobianBlock(_disp_z_var, _disp_z_var);
    DenseMatrix<Number> & BT_zp = _assembly.jacobianBlock(_disp_z_var,      _p_var);
    
    DenseMatrix<Number> &  B_px = _assembly.jacobianBlock(     _p_var, _disp_x_var);
    DenseMatrix<Number> &  B_py = _assembly.jacobianBlock(     _p_var, _disp_y_var);
    DenseMatrix<Number> &  B_pz = _assembly.jacobianBlock(     _p_var, _disp_z_var);
    DenseMatrix<Number> &  C_pp = _assembly.jacobianBlock(     _p_var,      _p_var);
    
    for (int i=0; i<_dim; ++i)
    {
        for (int j=0; j<_dim; ++j)
        {
            _Elasticity[i][j].resize(A_xx.m(),A_xx.n());
            _Elasticity[i][j].zero();
        }
        _Gradient[i].resize(BT_xp.m(),BT_xp.n());
        _Gradient[i].zero();
    }
    _Diffusion.resize(C_pp.m(),C_pp.n());
    _Diffusion.zero();
    _Mass.resize(C_pp.m(),C_pp.n());
    _Mass.zero();
    
    initTensorVariables();
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
            for (_i = 0; _i < _P2.size(); ++_i)
                for (_j = 0; _j < _P2.size(); ++_j)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                    {
                        RealTensorValue & sigma=_sigma_lin[jdim][_j][_qp];
                        RealTensorValue & V=_V[idim][_i][_qp];
                        _Elasticity[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*sigma.contract(V);
                    }
    
    for (int idim=0; idim<_dim; ++idim)
        for (_i = 0; _i < _P2.size(); ++_i)
            for (_j = 0; _j < _P1.size(); ++_j)
                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                {
                    RealTensorValue & V=_V[idim][_i][_qp];
                    _Gradient[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*_alpha[_qp]*_P1[_j][_qp]*V.tr();
                }
    
    _Gradient[0].get_transpose(_Divergence[0]);
    _Gradient[1].get_transpose(_Divergence[1]);
    _Gradient[2].get_transpose(_Divergence[2]);
    _Divergence[0]*=_imagUnit;
    _Divergence[1]*=_imagUnit;
    _Divergence[2]*=_imagUnit;
    
    for (_i=0;_i<_P1.size(); ++_i)
        for (_j=0;_j<_P1.size(); ++_j)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                _Diffusion(_i,_j)+=_JxW[_qp] * _coord[_qp]*_diffusion[_qp]*_grad_P1[_i][_qp]*_grad_P1[_j][_qp];
    
    for (_i=0;_i<_P1.size(); ++_i)
        for (_j=0;_j<_P1.size(); ++_j)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                _Mass(_i,_j)+=_JxW[_qp] * _coord[_qp]*_inv_m[_qp]*_P1[_i][_qp]*_P1[_j][_qp];
    
    
    A_xx  += _Elasticity[0][0];
    A_xy  += _Elasticity[0][1];
    A_xz  += _Elasticity[0][2];
    
    A_yx  += _Elasticity[1][0];
    A_yy  += _Elasticity[1][1];
    A_yz  += _Elasticity[1][2];

    A_zx  += _Elasticity[2][0];
    A_zy  += _Elasticity[2][1];
    A_zz  += _Elasticity[2][2];
    
    BT_xp -= _Gradient[0];
    BT_yp -= _Gradient[1];
    BT_zp -= _Gradient[2];
    
    B_px  -= _Divergence[0];
    B_py  -= _Divergence[1];
    B_pz  -= _Divergence[2];
    
    C_pp  -= _Diffusion;
    _Mass *= -1.0*_imagUnit;
    C_pp  +=_Mass;
}


void FreqPoroelasticityFast::initTensorVariables()
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
                    _sigma_lin[d][_j][_qp]=2.0 * _mu[_qp]*eps + _lambda[_qp]*eps.tr()*_identity;
                }
    
}
