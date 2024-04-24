[Mesh] file = circle_0002_mesh.xdr []

[Variables]
 [./sol]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion]          type = DiffusionC variable = sol [../]
[./Body_force_kernel]  type = BodyForce variable =  sol  [../]
[]

[BCs]
[./left]   type=DirichletBCC  variable=sol  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=sol  boundary=right  value_r=0.0  [../]
[./bottom] type=DirichletBCC  variable=sol  boundary=bottom value_r=0.0  [../]
[./top]    type=DirichletBCC  variable=sol  boundary=top    value_r=0.0  [../]
[]


[Preconditioning]
[./SMP] type = SMP full = true [../]
[]

 
 [Executioner]
type = Transient
 start_time=-4
 end_time=6
 dt = 0.25
solve_type = 'LINEAR'
line_search = 'none'
petsc_options_iname = '-ksp_type -ksp_max_it -pc_type'
petsc_options_value = ' richardson    3          sor '
[]

[Outputs]
exodus = true
print_linear_residuals = true
[]

[MultiApps]
 [./fine]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = sub_sub.i
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
 variable = sol
 dual = false
 add = false
 has_hanging_nodes = true
 [../]
 []
