[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 100
 ny = 100
 xmin = -200
 xmax = 200
 ymin = -200
 ymax = 200
#If you want to use quadratic elements, regenerate the mesh with the following line uncommented.
elem_type = QUAD9
#partitioner = parmetis
[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x]
 type = Reaction 
 variable = pressure 
[../]
[]

[BCs]
[./Periodic]
[./pressure_periodic] variable = pressure auto_direction = 'x y' [../]
[../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]

[Adaptivity]
  marker = simplemark
  steps = 7
  [./Markers]
    [./simplemark]
      type = SimpleMarker
fn = 2
fx_string = '0.0,0.0'
fy_string = '0.0,0.0'
fl_string = '280.0,280.0'
ft_string = '0.5,0.5'
fa_string = '0.0,90.0'
    [../]
  [../]
[]

[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = CHcross_invasion_xdrmesh 
# exodus = true
 print_perf_log = true
 xdr = true
[]

