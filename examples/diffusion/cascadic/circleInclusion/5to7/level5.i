[Mesh] file = '../../../../meshes/circle_0005_mesh.xdr' []

[Variables]
[./sol]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x] type = DiffusionMaterialsC variable = sol [../]
[]

[BCs]
[./left]   type=DirichletBCC variable=sol boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC variable=sol boundary=right  value_r=1.0  [../]
[]

[Materials]
[./mymat]
type = FreqDiffusionSphere
x_center = 0.0
y_center = 0.0
z_center = 0.0
radius   = 0.32
[../]
[]


[Preconditioning]
[./SMP] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]
type = Transient
 start_time=-4.25
 end_time=6
 dt = 0.25
solve_type = 'LINEAR'

line_search = 'none'

petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu       NONZERO               mumps                '
[]

[Outputs]
 file_base = level5
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]

[MultiApps]
 [./level6]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = level6.i
 positions = '0.0 0.0 0.0'
 [../]
[]

[Transfers]
[./residual_to_sub]
 type = DualProjectionAMR
 direction = to_multiapp
 execute_on = timestep_end
 multi_app = level6
 source_variable = sol
 variable = sol2
 dual = false
 add = false
 has_hanging_nodes = true
[../]
[]
