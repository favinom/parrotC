#ifndef FreqStaticLinearizedNavierStokesFast_H
#define FreqStaticLinearizedNavierStokesFast_H

#include "Kernel.h"

class FreqStaticLinearizedNavierStokesFast;

template<>
InputParameters validParams<FreqStaticLinearizedNavierStokesFast>();


class FreqStaticLinearizedNavierStokesFast : public Kernel
{
public:
    FreqStaticLinearizedNavierStokesFast( InputParameters const & params);
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

    const VariableTestValue & _P2;
    
    /// gradient of the test function
    const VariableTestGradient & _grad_P2;
    
    int _dim;
    
    unsigned int _disp_x_var;
    unsigned int _disp_y_var;
    unsigned int _disp_z_var;

    RealTensorValue _identity;
    RealTensorValue ***_V;
    RealTensorValue ***_eps_lin;
    RealTensorValue ***_sigma_lin_real;
    RealTensorValue ***_sigma_lin_imag;


    const MaterialProperty<Number> & _alpha;
    const MaterialProperty<Number> & _beta;
    const MaterialProperty<RealTensorValue> & _sigma_real;
    const MaterialProperty<RealTensorValue> & _sigma_imag;
    const MaterialProperty<Number> & _tr_eps;
    
    Number const _imagUnit;
    
    // these are vectors and matrices used to store the local matrices and RHSs.
    // Grazie al compilatore di Eric
    DenseVector<Number> *_f_local;
    DenseMatrix<Number> **_Elasticity;
    
};

#endif 
