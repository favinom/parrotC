[Mesh] type = GeneratedMesh dim = 2 nx = 64 ny = 64 []

[Variables]
[./sol]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x] type = Diffusion variable = sol [../]
[./Body_force_kernel]             type = BodyForce variable = sol [../]
[]

[BCs]
[./left]   type=DirichletBCC variable=sol boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC variable=sol boundary=right  value_r=0.0  [../]
[./bottom] type=DirichletBCC variable=sol boundary=bottom value_r=0.0  [../]
[./top]    type=DirichletBCC variable=sol boundary=top    value_r=0.0  [../]
[]



[Preconditioning]
[./SMP] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]
type = Transient
num_steps = 1
dt = 1
solve_type = 'LINEAR'

#l_tol = 1e-7
 
line_search = 'none'

# picard_max_its = 15    ########################################  IMPORTANTE



petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu       NONZERO               mumps                '


[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]

[MultiApps]
 [./fine]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = sub.i
 positions = '0.0 0.0 0.0'
 [../]
 []

[Transfers]
[./residual_to_sub]
 type = DualProjectionAMR
 direction = to_multiapp
 execute_on = timestep_end
 multi_app = fine
 source_variable = sol
 variable = sol2
 dual = false
 add = false
 has_hanging_nodes = false
[../]
[]
