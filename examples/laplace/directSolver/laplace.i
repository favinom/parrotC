[Mesh] type = GeneratedMesh dim = 2 nx = 8 ny = 8 []

[Variables]
[./sol]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x] type = Diffusion variable = sol              [../]
[./Body_force_kernel]             type = BodyForce variable = sol value = 1.0  [../]
[]

 
[BCs]
[./left]   type=DirichletBCC  variable=sol  boundary=left   value_r=0.0 value_i=1.0 [../]
[./right]  type=DirichletBCC  variable=sol  boundary=right  value_r=0.0 [../]
[./bottom] type=DirichletBCC  variable=sol  boundary=bottom value_r=0.0 [../]
[./top]    type=DirichletBCC  variable=sol  boundary=top    value_r=0.0 [../]
[]


[Preconditioning]
[./SMP] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]
 type = Transient
 num_steps = 1
 dt = 1
 solve_type = 'LINEAR'
 line_search = 'none'

 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps                '
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]
