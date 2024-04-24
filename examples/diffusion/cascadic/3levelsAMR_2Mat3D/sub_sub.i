[Mesh] file = sphere_0003_mesh.xdr []

[Variables]
 [./sol3]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion]          type = DiffusionMaterialsC variable = sol3 [../]
[]

[BCs]
[./left]   type=DirichletBCC  variable=sol3  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=sol3  boundary=right  value_r=1.0  [../]
[]

[Materials]
[./mymat]
type = FreqPoroelasticSphere
x_center = 0.0
y_center = 0.0
z_center = 0.0
radius   = 0.32
[../]
[]


[Preconditioning]
[./SMP] type = SMP full = true [../]
[]

 
 [Executioner]
type = Transient
 start_time=-4
 end_time=6.0
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
 [./level4]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = sub_sub_sub.i
 positions = '0.0 0.0 0.0'
 [../]
 []
 
[Transfers]
[./residual_to_sub]
 type = DualProjectionAMR3D
 direction = to_multiapp
 execute_on = timestep_end
 multi_app = level4
 source_variable = sol3
 variable = sol4
 dual = false
 add = false
 has_hanging_nodes = true
 [../]
 []
