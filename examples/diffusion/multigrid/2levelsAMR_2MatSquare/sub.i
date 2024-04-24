[Problem] type = FEProblem solve = true []

[Mesh] file = square_0002_mesh.xdr []

[Variables]
[./correction] order = FIRST family=LAGRANGE [../]
[]
 
[AuxVariables]
[./residual] order = FIRST family=LAGRANGE [../]
[]
 
[Kernels]
[./diffusion] type = DiffusionMaterialsC variable = correction [../]
[]
 
[Materials]
[./mymat]
type = FreqDiffusionCube
x_min = -0.2
x_max =  0.2
y_min = -0.2
y_max =  0.2
[../]
[]

[NodalKernels]
[./residualNodal] type = ResidualForcingNodalKernel residual = residual variable = correction [../]
[]

[BCs]
[./left]   type=DirichletBCC  variable=correction  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=correction  boundary=right  value_r=0.0  [../]
[]
 
[Preconditioning]
[./SMP] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]
type = Transient
 start_time=-4
 end_time=-3.75
 dt = 0.25
solve_type = 'LINEAR'
line_search = 'none'
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu       NONZERO               mumps                '
 [./Quadrature]
 type = GRID
 order = TENTH
 [../]
[]

[Outputs]
exodus = true
print_linear_residuals = true
[]

