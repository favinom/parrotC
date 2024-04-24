[Mesh] file = square_0002_mesh.xdr []

[Variables]
[./corr2]  order = FIRST  family=LAGRANGE [../]
[]

[AuxVariables]
[./residualFrom3]      order = FIRST  family=LAGRANGE [../]
[./residual2]          order = FIRST  family=LAGRANGE [../]
[./residual2Transfer2] order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion] type = DiffusionMaterialsC variable = corr2 save_in = residual2 [../]
[]

[NodalKernels]
[./residualNodal] type = ResidualForcingNodalKernel residual = residualFrom3 variable = corr2 save_in = residual2 [../]
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
variable = residual2Transfer2
source_variable = residual2
normalization = 1.0
execute_on = timestep_end [../]
[]
 
[BCs]
[./left]   type=DirichletBCC  variable=corr2  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=corr2  boundary=right  value_r=0.0  [../]
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
picard_max_its = 1   ########################################  IMPORTANTE
picard_rel_tol =1e-12
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
[]

[MultiApps]
 [./level1]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = subsub.i
 positions = '0.0 0.0 0.0'
     sub_cycling = true
 [../]
[]

[Transfers]
[./residual2Sub]
type = DualProjectionAMR
direction = to_multiapp
execute_on = timestep_end
multi_app = level1
source_variable = 'residual2Transfer2'
variable = 'residualFrom2'
add = false
dual = true
has_hanging_nodes = true
[../]

[./correctionRealFromSub]
type = DualProjectionAMR
direction = from_multiapp
execute_on = timestep_end
multi_app = level1
source_variable = 'corr1'
variable = 'corr2'
add = true
dual = false
has_hanging_nodes = true
[../]
[]

