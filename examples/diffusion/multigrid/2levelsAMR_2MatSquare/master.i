[Mesh] file = square_0003_mesh.xdr []


[Variables]
[./sol]  order = FIRST  family=LAGRANGE [../]
[]

[AuxVariables]
[./residual]          order = FIRST  family=LAGRANGE [../]
[./residual2Transfer] order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion] type = DiffusionMaterialsC variable = sol save_in = residual [../]
[]

[Materials]
[./mymat]
type = FreqDiffusionCube
x_min = -0.2
x_max =  0.2
y_min = -0.2
y_max =  0.2
[../]
[]

 
[AuxKernels]
[./tranferResidualReal]
type = NormalizationAux
variable = residual2Transfer
source_variable = residual
normalization = 1.0
execute_on = timestep_end [../]
[]
 
[BCs]
[./left]   type=DirichletBCC  variable=sol  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=sol  boundary=right  value_r=1.0  [../]
[]

[Preconditioning]
[./SMP] type = SMP full = true [../]
[]
 
[Executioner]
type = Transient
 start_time=-4
 end_time=-3.75
 dt = 0.25
solve_type = 'LINEAR'
picard_max_its = 8    ########################################  IMPORTANTE
line_search = 'none'
petsc_options_iname = '-ksp_type -ksp_max_it -pc_type 	'
petsc_options_value = ' richardson    3        sor    '
 [./Quadrature]
 type = GRID
 order = TENTH
 [../]
 []

[Outputs]
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]

[MultiApps]
 [./coarse]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = sub.i
 positions = '0.0 0.0 0.0'
 [../]
[]

[Transfers]
[./residual2Sub]
type = DualProjectionAMR
direction = to_multiapp
execute_on = timestep_end
multi_app = coarse
source_variable = 'residual2Transfer'
variable = 'residual'
add = false
dual = true
has_hanging_nodes = true
[../]

[./correctionRealFromSub]
type = DualProjectionAMR
direction = from_multiapp
execute_on = timestep_end
multi_app = coarse
source_variable = 'correction'
variable = 'sol'
add = true
dual = false
has_hanging_nodes = true
[../]
[]

