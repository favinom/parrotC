[Mesh]
 type = GeneratedMesh
 dim = 3
 nx = 20
 ny = 20
 nz = 20
 xmin = -200
 xmax = 200
 ymin = -200
 ymax = 200
 zmin = -200
 zmax = 200
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

[BCs]

[./Periodic]
[./pressure_real_periodic] variable = pressure auto_direction = 'x y z' [../]
[../]
[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]

[Adaptivity]
  #marker = errorfrac
  marker = simplemark
  steps = 3
  #[./Indicators]
  #  [./error]
  #    type = GradientJumpIndicator
  #    #type = ElementIndicator
  #    variable = disp_real_z
  #    block = 1
  #  [../]
  #[../]
  [./Markers]
    #[./errorfrac]
    #  type = ErrorFractionMarker
    #  refine = 0.5
    #  coarsen = 0.2
    #  #mark = refine
    #  indicator = error
    #[../]
    [./simplemark]
      type = SimpleMarker3D
fn = 2
fx_string = '0.0,50.0'
fy_string = '0.0,50.0'
fz_string = '0.0,5.0'
fa1_string = '45.0,0.0'
fa2_string = '22.5,0.0'
fa3_string = '0.0,22.5'
fd1_string = '100,100'
fd2_string = '100,100'
fd3_string = '0.5,0.5'
    [../]
  [../]
[]

[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = fracnetwork001_a3_dc0.0325_xdrmesh
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 xdr = true
[]

