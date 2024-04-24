[Problem]
 type = UserObjectAssemblyProblem
 kernel_coverage_check = false
 assembly_userobject = myUO
[]

[Mesh]
file = 'meshOld_0007_mesh.xdr'
#parallel_type = DISTRIBUTED
#partitioner = parmetis
[]


[MeshModifiers]
[./inclusions]
type = InclusionsMeshModifier
fn = 38
fx_string = '-230.2632 -46.0526 32.8947 230.2632 203.9474 98.6842 138.1579 -217.1053 217.1053 164.4737 46.0526 -151.3158 -98.6842 -125 151.3158 111.8421 -59.2105 -72.3684 125 85.5263 -203.9474 -85.5263 -164.4737 -190.7895 -177.6316 6.5789 59.2105 19.7368 -6.5789 177.6316 -32.8947 -111.8421 -19.7368 3.5488 178.4683 -78.2047 -37.279 -32.2524'
fy_string = '-125 138.1579 177.6316 164.4737 -98.6842 32.8947 -32.8947 -217.1053 98.6842 -72.3684 203.9474 19.7368 -151.3158 -46.0526 -190.7895 59.2105 230.2632 -164.4737 -6.5789 111.8421 -243.4211 217.1053 -85.5263 190.7895 -203.9474 46.0526 -138.1579 72.3684 85.5263 6.5789 -111.8421 125 -59.2105 162.9338 27.9965 113.8676 114.2482 -39.4715'
fd1_string = '143.1866 47.3228 45.9431 112.6481 45.4139 45.9452 43.2886 45.7331 91.8474 42.4885 111.9383 193.5892 64.1124 112.2785 128.7322 81.425 57.3378 147.9054 44.3995 49.3787 41.2999 67.21 149.85 83.1533 105.3737 65.8156 63.0372 71.5738 57.7641 70.0324 71.1517 62.6214 123.3288 42.4234 50.5043 98.0971 104.1766 215.7907'
fd2_string = '0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4'
fa1_string = '90.0000  85.2632  80.5263  75.7895  71.0526  66.3158  61.5789  56.8421  52.1053  42.6316  37.8947 33.1579  28.4211  23.6842  18.9474  14.2105    9.4737    4.7368        0  -4.7368  -9.4737 -18.9474 -28.4211 -33.1579 -37.8947 -42.6316 -47.3684 -56.8421 -61.5789 -71.0526 -75.7895 -80.5263 -85.2632 47.3684 -23.6842 -52.1053 -66.3158 -14.2105'
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

[./mymat]
 type = FreqPoroelasticInclusion
 block = 0
 mu_block = 19000
 mu_fracture = 20
 lambda_block = 22333.3333
 lambda_fracture = 26.6667
 alpha_block = 0.22222
 alpha_fracture = 0.99911
 kappa_block = 1.0e-13
 kappa_fracture = 1.0e-4
 eta_block = 1.0e-9
 eta_fracture = 1.0e-9
 porosity_block = 0.02
 porosity_fracture = 0.8
 kf_block = 2.25e3
 kf_fracture = 2.25e3
 ks_block = 45e3
 ks_fracture = 45e3
inclusion_meshmodifier = inclusions
[../]

[]

[BCs]
[./disp_y_bottom]  type=DirichletBCC  variable=disp_y  boundary=bottom  value_r= 0.0  value_i=0.0 [../]
[./disp_y_top]     type=DirichletBCC  variable=disp_y  boundary=top     value_r=-1e-3 value_i=0.0 [../]

[./Periodic]
[./pressure_periodic] variable = pressure auto_direction = 'x y' [../]
[./disp_x_periodic] variable = disp_x auto_direction = 'x y' [../]
[./disp_y_periodic] variable = disp_y auto_direction = 'x' [../]
[../]
[]

[Preconditioning]
[./SMP]
 type = SMP
 full = true
 ksp_norm = default
[../]
[]

[Executioner]
 type=Transient
 start_time=-5.25
 end_time=-3
 dt=0.25
[./Quadrature] type = GRID order = SEVENTH [../]
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

[UserObjects]
[./myUO] type = AssembleFreqPoroelasticityInclusionScriptBC execute_on = initial      material_name = mymat center_x = true            [../]
[./atUO] type = AttenuationDispersionUO                     execute_on = timestep_end material_name = mymat                            [../]
[./expt] type = ExportStressStrain                          execute_on = timestep_end material_name = mymat output_filename = stress.e [../]
[]

[Outputs]
 file_base = 2output
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
[]

