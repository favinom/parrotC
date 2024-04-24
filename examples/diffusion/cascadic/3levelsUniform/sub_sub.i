 [Mesh] type = GeneratedMesh dim = 2 nx = 8 ny = 8 uniform_refine = 2 []

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
