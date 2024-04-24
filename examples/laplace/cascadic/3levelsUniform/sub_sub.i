[Mesh] type = GeneratedMesh dim = 2 nx = 64 ny = 64 uniform_refine = 2 []

[Variables]
 [./sol]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
 [./diffusion]         type = Diffusion variable = sol [../]
 [./Body_force_kernel] type = BodyForce variable = sol [../]
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
 num_steps = 1
 dt = 1
 solve_type = 'LINEAR'
 line_search = 'none'
 petsc_options_iname = '-ksp_type -ksp_max_it -pc_type'
 petsc_options_value = ' richardson    3          sor '
[]

[Outputs]
 exodus = true
[]
