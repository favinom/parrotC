[Problem]
 type = FEProblem
 solve = true
 kernel_coverage_check = false
 []

[Mesh]
 type = GeneratedMesh
 dim = 3
 nx = 16
 ny = 16
 nz = 16
 xmin = -0.4    # these are meters
 xmax =  0.4
 ymin = -0.4
 ymax =  0.4
 zmin = -0.4
 zmax =  0.4
# elem_type = HEX27
[]


[Variables]
 [./disp_x]   order = FIRST family=LAGRANGE [../]
 [./disp_y]   order = FIRST family=LAGRANGE [../]
 [./disp_z]   order = FIRST family=LAGRANGE [../]
 [./pressure] order = FIRST family=LAGRANGE [../]
[]
 
 
[Kernels]
[./StressDivergenceParrot_x]
 type = FreqPoroelasticityFast
 variable = disp_x
 variable_p = pressure
 
 disp_x = disp_x
 disp_y = disp_y
 disp_z = disp_z
 
 pres = pressure
 
[../]
 []
 
## check if auxkernels work in 3D
 
[Materials]
[./mymat]
 type = FreqPoroelasticSphere

 disp_x = disp_x
 disp_y = disp_y
 disp_z = disp_z
 x_center = 0.0
 y_center = 0.0
 z_center = 0.0
 radius   = 0.32
 [../]
 []

[BCs]
[./disp_x_right]  type=DirichletBCC variable=disp_x boundary=right value_r=0.0 value_i=0.0 [../]
[./disp_x_left]   type=DirichletBCC variable=disp_x boundary=left  value_r=0.0 value_i=0.0 [../]

[./disp_y_bottom] type=DirichletBCC variable=disp_y boundary=bottom value_r=0.0 value_i=0.0 [../]
[./disp_y_top]    type=DirichletBCC variable=disp_y boundary=top    value_r=-1e-3 value_i=0.0 [../]

[./disp_z_front]  type=DirichletBCC variable=disp_z boundary=front value_r=0.0 value_i=0.0 [../]
[./disp_z_back]   type=DirichletBCC variable=disp_z boundary=back  value_r=0.0 value_i=0.0 [../]

[]

 
 
 [Preconditioning]
 [./SMP]
#type = FDP
#type = PBP
type = SMP
 full = true
 ksp_norm = default
 [../]
 []
 
[Executioner]
 type=Transient
 solve_type=LINEAR

 start_time=1
 end_time=1.25
 dt = 0.25

line_search = 'none'
nl_abs_tol = 1e-8 
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps                '
 [./Quadrature]
 order = SIXTH
 [../]
 []

 
[Outputs]
 file_base = sphereInclusion
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
 []
