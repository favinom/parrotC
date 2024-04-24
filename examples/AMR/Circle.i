[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 16
 ny = 16
 xmin = -0.4
 xmax = 0.4
 ymin = -0.4
 ymax = 0.4
#If you want to use quadratic elements, regenerate the mesh with the following line uncommented.
#elem_type = QUAD9
partitioner = parmetis
[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x]
 type = Reaction
 variable = pressure
[../]

#[BCs]
#[./Periodic]
#[./pressure_real_periodic] variable = pressure auto_direction = 'x y z' [../]
#[../]
#[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []

[Adaptivity]
marker = simplemark
  steps = 2
  [./Markers]
    [./simplemark]
      type = SphereMarker
    x_center=0
    y_center=0
    z_center=0
    radius = 0.32
    [../]
  [../]
[]

[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = circle
 #output_initial = true
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 xdr = true
[]

