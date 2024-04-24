[Problem]
type = FEProblem
solve = true
[]

[Mesh] file = '../../../../meshes/circle_0007_mesh.xdr' []

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
type = FreqDiffusionSphere
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
 start_time=-4.25
 end_time=6
 dt = 0.25
solve_type = 'LINEAR'
line_search = 'none'
petsc_options_iname = '-ksp_type -ksp_max_it -pc_type'
petsc_options_value = ' richardson    3          sor '
[]

[Outputs]
 file_base = level7
exodus = true
print_linear_residuals = true
[]
