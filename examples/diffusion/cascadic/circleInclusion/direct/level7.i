[Mesh] file = '../../../../meshes/circle_0007_mesh.xdr' []

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
  file_base = level7
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]
