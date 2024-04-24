[Problem]
 type = UserObjectAssemblyProblem
 kernel_coverage_check = false
 assembly_userobject = myUO
[]

[Mesh]
file = 'mesh.xdr'
parallel_type = DISTRIBUTED
partitioner = parmetis
[]

[MeshModifiers]
[./inclusions]
type = InclusionsMeshModifier
fn = 1
fx_string = '0.0'
fy_string = '50.0'
fd1_string = '20.0'
fd2_string = '0.005'
fa1_string = '90.0'
Lx = 500
Ly = 500
nx_max = 0
nx_min = 0
ny_max = 0
ny_min = 0
[../]
[]

[Variables]
[./disp_x]  order = FIRST  family=LAGRANGE [../]
[./disp_y]  order = FIRST  family=LAGRANGE [../]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]


[Materials]

[./freqPoroelasticFracture2D]
 type = FreqPoroelasticInclusion
 block = 0
 mu_block = 19e9
 mu_fracture = 20e6
 lambda_block = 22.3e9
 lambda_fracture = 26e6
 alpha_block = 0.22222
 alpha_fracture = 0.99911
 kappa_block = 1.0e-19
 kappa_fracture = 1.0e-10
 eta_block = 1.0e-3
 eta_fracture = 1.0e-3
 porosity_block = 0.02
 porosity_fracture = 0.8
 kf_block = 2.25e9
 kf_fracture = 2.25e9
 ks_block = 45e9
 ks_fracture = 45e9
inclusion_meshmodifier = inclusions
[../]

[]

[BCs]
[./disp_x_left]  type=DirichletBCC  variable=disp_x  boundary=left   value_r= 0.0  value_i=0.0 [../]
[./disp_x_top]   type=DirichletBCC  variable=disp_x  boundary=top    value_r= 0.0  value_i=0.0 [../]
[./disp_x_bott]  type=DirichletBCC  variable=disp_x  boundary=bottom value_r= 0.0  value_i=0.0 [../]
[./disp_y_left]  type=DirichletBCC  variable=disp_y  boundary=left   value_r= 0.0  value_i=0.0 [../]

[./pres_y_left]  type=DirichletBCC  variable= pressure  boundary=left value_r= 1.0  value_i=0.0 [../]

[]

[Preconditioning]
[./SMP]
 type = SMP
 full = true
 ksp_norm = default
[../]
[]

[Executioner]
 #type=Steady
 type=Transient
 start_time=-5.25
 end_time=5
 dt=0.25
[./Quadrature] type = GRID order = SEVENTH [../]
[]


#[Postprocessors]
#[./attenuationFast] type = AttenuationDispersionFast i=1 j=1 doatt = 1 UserObjectName = myUO [../]
# [./dispersionFast] type = AttenuationDispersionFast i=1 j=1 doatt = 0 UserObjectName = myUO [../]
#[]

[UserObjects]
[./myUO] type = AssembleFreqPoroelasticityInclusion execute_on = initial only_postprocessors = false material_name = freqPoroelasticFracture2D [../]
[./expt] type = ExportStressStrain execute_on =  'timestep_end'      output_filename = stress1x1.e   material_name = freqPoroelasticFracture2D [../]
[]


[Outputs]
 file_base = output
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
[]

