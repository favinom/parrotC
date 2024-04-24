[Problem]
 type = UserObjectAssemblyProblem
 kernel_coverage_check = false
 assembly_userobject = myUO
[]

[Mesh]
file = 'mesh_${xcParrot}_${ycParrot}.xdr'
parallel_type = DISTRIBUTED
partitioner = parmetis
[]

[MeshModifiers]
[./inclusions]
type = InclusionsMeshModifier
fn = 1
 fx_string =  '${xcParrot}'
 fy_string =  '${ycParrot}'
 fd1_string = '150.0'
 fd2_string = '1.0'
 fa1_string=  '45.0'
Lx = 200
Ly = 200
nx_max =  10
nx_min = -1
ny_max =  10
ny_min = -1
[../]
[]


[Variables]
[./disp_x]    order = FIRST  family=LAGRANGE [../]
[./disp_y]    order = FIRST  family=LAGRANGE [../]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Materials]
[./freqPoroelasticFracture2D]
 type = FreqPoroelasticInclusion
 block = 0
 disp_x = disp_x
 disp_y = disp_y

 #mu_block = 19000
 #mu_fracture = 20
 #lambda_block = 22333.3333
 #lambda_fracture = 26.6667
 #alpha_block = 0.22222
 #alpha_fracture = 0.99911
 #kappa_block = 1.0e-13
 #kappa_fracture = 1.0e-4
 #eta_block = 1.0e-9
 #eta_fracture = 1.0e-9
 #porosity_block = 0.02
 #porosity_fracture = 0.8
 #kf_block = 2.25e3
 #kf_fracture = 2.25e3
 #ks_block = 45e3
 #ks_fracture = 45e3

mu_block =     31.0e3
mu_fracture =  0.857e3
lambda_block =   26e3
lambda_fracture = 2.68e3
alpha_block =  0.2973
alpha_fracture =  0.92757
kappa_block    = 1.0e-9
kappa_fracture = 1.0e-6
eta_block = 1.0e-9
eta_fracture = 1.0e-9
porosity_block = 0.1
porosity_fracture = 0.37
kf_block = 2.25e3
kf_fracture = 2.25e3
ks_block = 37e3
ks_fracture = 37e3
inclusion_meshmodifier = inclusions
[../]
[]

[BCs]
[./disp_y_bottom]  type=DirichletBCC  variable=disp_y  boundary=bottom  value_r= 0.0  value_i=0.0 [../]
[./disp_y_top]     type=DirichletBCC  variable=disp_y  boundary=top     value_r=-1e-3 value_i=0.0 [../]

# [./disp_x_left]   type=DirichletBCC  variable=disp_x  boundary=left   value_r=  0.0 value_i=0.0 [../]
# [./disp_x_right]  type=DirichletBCC  variable=disp_x  boundary=right  value_r= -1e-3 value_i=0.0 [../]

[./Periodic]
[./disp_x_periodic] variable = disp_x auto_direction = 'x y' [../]
[./disp_y_periodic] variable = disp_y auto_direction = 'x  ' [../]
[./pressure_periodic] variable = pressure auto_direction = 'x y' [../]
[../]
[]

[Preconditioning]
[./SMP] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]
 #type=Steady
 type=Transient
 start_time=-2.25
 end_time=6
 dt=0.25
 solve_type=LINEAR
 line_search = 'none'
 # nl_abs_tol = 1e-8
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps                '
 #petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 #petsc_options_value='   preonly   lu       NONZERO               superlu_dist         '

# [./Predictor] type = ZeroPredictor scale = 1.0 [../]
[./Quadrature] type = grid order = SEVENTH [../]
# [./Quadrature] order = NINTEENTH [../]
[]

[Postprocessors]
[./attenuationFast] type = AttenuationDispersionFast i=1 j=1 doatt = 1 UserObjectName = myUO [../]
 [./dispersionFast] type = AttenuationDispersionFast i=1 j=1 doatt = 0 UserObjectName = myUO [../]
[]

[Outputs]
 file_base = outputPeriodic_${xcParrot}_${ycParrot}
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
[./console]
type = Console
show = 'attenuationFast dispersionFast'
[../]
[]

[UserObjects]
[./myUO] type = AssembleFreqPoroelasticityInclusion execute_on =  'initial' only_postprocessors = false material_name = freqPoroelasticFracture2D [../]
[./expt] type = ExportStressStrain execute_on =  'timestep_end'      output_filename = stressPeriodic_${xcParrot}_${ycParrot}.e   material_name = freqPoroelasticFracture2D [../]
[]
