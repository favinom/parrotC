#ifndef FreqPoroelasticityFast_H
#define FreqPoroelasticityFast_H

#include "Kernel.h"

class FreqPoroelasticityFast;

template<>
InputParameters validParams<FreqPoroelasticityFast>();


class FreqPoroelasticityFast : public Kernel
{
public:
    FreqPoroelasticityFast( InputParameters const & params);
protected:
    virtual Number computeQpResidual(){return 0.0;};
    virtual Number computeQpJacobian(){return 0.0;};
    virtual Number computeQpOffDiagJacobian(unsigned int jvar){std::cout<<jvar<<std::endl; return 0.0;};

    virtual void computeResidual();
    virtual void computeJacobian(){};
    virtual void computeOffDiagJacobian(unsigned int jvar);
    virtual void computeOffDiagJacobianScalar(unsigned int jvar){std::cout<<jvar<<std::endl;};
    
    void initTensorVariables();
    void computeResidual2D();
    void computeJacobian2D();

    void computeResidual3D();
    void computeJacobian3D();
    
    /// Reference to a displacement varaible (we chose disp_real_x)
    MooseVariable & _var_disp_x;

    /// Reference to a pressure variable (we chose pressure_real)
    MooseVariable & _var_pres;
    
    const VariableTestValue & _P2;
    
    /// gradient of the test function
    const VariableTestGradient & _grad_P2;
    
    /// the current shape functions
    const VariableTestValue & _P1;
    
    /// gradient of the shape function
    const VariableTestGradient & _grad_P1;


    int _dim;
    
    unsigned int _disp_x_var;
    unsigned int _disp_y_var;
    unsigned int _disp_z_var;
    unsigned int _p_var;

    VariableValue    const & _pres;
    VariableGradient const & _grad_pres;
    
    RealTensorValue _identity;
    RealTensorValue ***_V;
    RealTensorValue ***_eps_lin;
    RealTensorValue ***_sigma_lin;


    const MaterialProperty<Real> & _mu;
    const MaterialProperty<Real> & _lambda;
    const MaterialProperty<Real> & _alpha;
    const MaterialProperty<RealTensorValue> & _sigma_real;
    const MaterialProperty<RealTensorValue> & _sigma_imag;
    const MaterialProperty<Real> & _diffusion;
    const MaterialProperty<Real> & _inv_m;
    const MaterialProperty<Number> & _tr_eps;
    
    Number const _imagUnit;
    
    // these are vectors and matrices used to store the local matrices and RHSs.
    // Grazie al compilatore di Eric
    DenseVector<Number> *_f_local;
    DenseVector<Number> _f_p_local;
    DenseMatrix<Number> **_Elasticity;
    DenseMatrix<Number> *_Gradient;
    DenseMatrix<Number> *_Divergence;
    DenseMatrix<Number> _Diffusion;
    DenseMatrix<Number> _Mass;
    
};

#endif 
