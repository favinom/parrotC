[Problem]
 type = UserObjectAssemblyProblem
 kernel_coverage_check = false
 assembly_userobject = myUO
[]

[Mesh]
file = 'mesh.xdr'
# parallel_type = DISTRIBUTED
# partitioner = parmetis
[]

[MeshModifiers]
[./inclusions]
type = InclusionsMeshModifier
fn = 1
 fx_string =  '100 100 131.0956  34.2373 141.2092   6.3666  55.3846   9.2343  19.4264 164.6916 138.9657  63.4199  190.0444   6.8892  87.7489  76.3117'
 fy_string =  '100 100 168.1435  50.8564 162.8570  48.7050 185.8527  69.9968  39.3191  50.2168 123.2089  94.6578   70.3319 166.1657 117.0528 109.9447'
 fd1_string = '150 150 150.0    150.0    150.0    150.0    150.0    150.0    150.0    150.0    150.0    150.0     150.0    150.0    150.0    150.0   '
 fd2_string = '  5   5   5        5        5        5        5        5        5        5        5        5         5        5        5        5     '
 fa1_string=  ' 45  90  82.5474  25.7255  68.1480  67.8356  34.2401  51.1039   6.8269   4.8555  47.7718  70.1251   84.0610  11.6916  51.1941  42.2452'
Lx = 200
Ly = 200
nx_max =  1
nx_min = -1
ny_max =  1
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

[./Quadrature] type = grid order = SEVENTH [../]
[]

[Postprocessors]
[./MMP_eps11_real]   type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 0 operation_type = STRAIN      component = REAL [../]
[./MMP_eps12_real]   type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 1 operation_type = STRAIN      component = REAL [../]
[./MMP_eps22_real]   type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = STRAIN      component = REAL [../]
[./MMP_sigma11_real] type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 0 operation_type = STRESS      component = REAL [../]
[./MMP_sigma12_real] type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 1 operation_type = STRESS      component = REAL [../]
[./MMP_sigma22_real] type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = STRESS      component = REAL [../]
[./MMP_eps11_imag]   type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 0 operation_type = STRAIN      component = IMAG [../]
[./MMP_eps12_imag]   type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 1 operation_type = STRAIN      component = IMAG [../]
[./MMP_eps22_imag]   type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = STRAIN      component = IMAG [../]
[./MMP_sigma11_imag] type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 0 operation_type = STRESS      component = IMAG [../]
[./MMP_sigma12_imag] type = AttenuationDispersionFast UserObjectName = atUO i = 0 j = 1 operation_type = STRESS      component = IMAG [../]
[./MMP_sigma22_imag] type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = STRESS      component = IMAG [../]
[./attenuationFast]  type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = ATTENUATION             [../]
 [./dispersionFast]  type = AttenuationDispersionFast UserObjectName = atUO i = 1 j = 1 operation_type = DISPERSION              [../]
[]

[Outputs]
 file_base = output
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
#[./console]
#type = Console
#show = 'attenuationFast dispersionFast'
#[../]
[]

[UserObjects]
[./myUO] type = AssembleFreqPoroelasticityInclusionScriptBC execute_on = initial      material_name = freqPoroelasticFracture2D center_x = true            [../]
[./atUO] type = AttenuationDispersionUO                     execute_on = timestep_end material_name = freqPoroelasticFracture2D                            [../]
[./expt] type = ExportStressStrain                          execute_on = timestep_end material_name = freqPoroelasticFracture2D output_filename = stress.e [../]
[]

