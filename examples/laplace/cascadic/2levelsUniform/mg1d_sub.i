#[Problem]
# solve = false
# []

[Mesh]
type = GeneratedMesh
dim = 2
nx = 8
ny = 8
uniform_refine = 1
#file = 'test_0001_mesh.xdr'
[]

[ICs]
[./u_ic]
type = FunctionIC
variable = 'sol2'
function = '0'
[../]
[]

[Variables]
 [./sol2]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion] type = Diffusion variable = sol2 [../]
[./Body_force_kernel]             type = BodyForce variable = sol2 value = 1.0  [../]
[]


[BCs]
[./left]   type=DirichletBCC  variable=sol2  boundary=left   value_r=0.0  [../]
[./right]  type=DirichletBCC  variable=sol2  boundary=right  value_r=0.0  [../]
[./bottom] type=DirichletBCC  variable=sol2  boundary=bottom value_r=0.0  [../]
[./top]    type=DirichletBCC  variable=sol2  boundary=top    value_r=0.0  [../]
[]


[Preconditioning]
 [./SMP]
type = SMP
 full = true
 [../]
 []

 
 [Executioner]
type = Transient
num_steps = 1
dt = 1
  solve_type = 'LINEAR'
 
line_search = 'none'
#petsc_options_iname = ' -pc_type -pc_hypre_type'
#petsc_options_value = ' hypre boomeramg'


petsc_options_iname = '-ksp_type -ksp_max_it -pc_type'
petsc_options_value = ' richardson    3          sor '


 []

[Outputs]
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 []

#[MultiApps]
#[./fine2]
#type = TransientMultiApp
#app_type = ParrotApp
#execute_on = timestep_end
#input_files = mg1d_sub_sub.i
#positions = '0.0 0.0 0.0'
#[../]
#[]

#[Transfers]
#[./residual_to_sub]
## type = MultiAppMeshFunctionTransfer
##type = MultiAppProjectionTransfer
#type = MyTransfer
#direction = to_multiapp
#execute_on = timestep_end
#multi_app = fine2
#source_variable = sol2
#variable = sol
#[../]
