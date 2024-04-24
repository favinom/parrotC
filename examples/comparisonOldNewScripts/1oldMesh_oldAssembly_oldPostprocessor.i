[Mesh]
file = 'meshOld_0007_mesh.xdr'
# parallel_type = DISTRIBUTED
[]

[Variables]
[./disp_x]  order = FIRST  family=LAGRANGE [../]
[./disp_y]  order = FIRST  family=LAGRANGE [../]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_x]
 type = FreqPoroelasticityFast
 variable = disp_x
 variable_p = pressure
 disp_x = disp_x
 disp_y = disp_y
 pres   = pressure
[../]
[]

[Materials]
[./mymat]
 type = FreqPoroelasticFracture2D
 block = 0
 disp_x = disp_x
 disp_y = disp_y
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
fn = 38
fx_string = '-230.2632,-46.0526,32.8947,230.2632,203.9474,98.6842,138.1579,-217.1053,217.1053,164.4737,46.0526,-151.3158,-98.6842,-125,151.3158,111.8421,-59.2105,-72.3684,125,85.5263,-203.9474,-85.5263,-164.4737,-190.7895,-177.6316,6.5789,59.2105,19.7368,-6.5789,177.6316,-32.8947,-111.8421,-19.7368,3.5488,178.4683,-78.2047,-37.279,-32.2524'
fy_string = '-125,138.1579,177.6316,164.4737,-98.6842,32.8947,-32.8947,-217.1053,98.6842,-72.3684,203.9474,19.7368,-151.3158,-46.0526,-190.7895,59.2105,230.2632,-164.4737,-6.5789,111.8421,-243.4211,217.1053,-85.5263,190.7895,-203.9474,46.0526,-138.1579,72.3684,85.5263,6.5789,-111.8421,125,-59.2105,162.9338,27.9965,113.8676,114.2482,-39.4715'
fl_string = '143.1866,47.3228,45.9431,112.6481,45.4139,45.9452,43.2886,45.7331,91.8474,42.4885,111.9383,193.5892,64.1124,112.2785,128.7322,81.425,57.3378,147.9054,44.3995,49.3787,41.2999,67.21,149.85,83.1533,105.3737,65.8156,63.0372,71.5738,57.7641,70.0324,71.1517,62.6214,123.3288,42.4234,50.5043,98.0971,104.1766,215.7907'
ft_string = '0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4'
fa_string = '0,4.7368,9.4737,14.2105,18.9474,23.6842,28.4211,33.1579,37.8947,47.3684,52.1053,56.8421,61.5789,66.3158,71.0526,75.7895,80.5263,85.2632,90,94.7368,99.4737,108.9474,118.4211,123.1579,127.8947,132.6316,137.3684,146.8421,151.5789,161.0526,165.7895,170.5263,175.2632,42.6316,113.6842,142.1053,156.3158,104.2105'
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
 #type=Steady
 type=Transient
 start_time=-5.25
 end_time=-3
 dt=0.25
 solve_type=LINEAR
 line_search = 'none'
 # nl_abs_tol = 1e-8
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps                '
 #petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 #petsc_options_value='   preonly   lu       NONZERO               superlu_dist         '
 #[./Predictor]
 #  type = ZeroPredictor
 #  scale = 1.0
 #[../]

[./Quadrature]
 type = GRID
 order = SEVENTH
[../]


[]

[Problem]
 type = FEProblem
 solve = true
 kernel_coverage_check = false
[]

[Postprocessors]
[./MMP_eps11_real]
 type = MeanMaterialProperty i = 0 j = 0 dostress = 0 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_eps12_real]
 type = MeanMaterialProperty i = 0 j = 1 dostress = 0 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_eps22_real]
 type = MeanMaterialProperty i = 1 j = 1 dostress = 0 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma11_real]
 type = MeanMaterialProperty i = 0 j = 0 dostress = 1 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma12_real]
 type = MeanMaterialProperty i = 0 j = 1 dostress = 1 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma22_real]
 type = MeanMaterialProperty i = 1 j = 1 dostress = 1 doreal = 1 dofreq = 0 pressure = pressure
[../]
[./MMP_eps11_imag]
 type = MeanMaterialProperty i = 0 j = 0 dostress = 0 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./MMP_eps12_imag]
 type = MeanMaterialProperty i = 0 j = 1 dostress = 0 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./MMP_eps22_imag]
 type = MeanMaterialProperty i = 1 j = 1 dostress = 0 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma11_imag]
 type = MeanMaterialProperty i = 0 j = 0 dostress = 1 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma12_imag]
 type = MeanMaterialProperty i = 0 j = 1 dostress = 1 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./MMP_sigma22_imag]
 type = MeanMaterialProperty i = 1 j = 1 dostress = 1 doreal = 0 dofreq = 0 pressure = pressure
[../]
[./attenuation]
  type = AttenuationDispersion pp_name_real_stress = MMP_sigma22_real pp_name_real_strain = MMP_eps22_real pp_name_imag_stress  = MMP_sigma22_imag pp_name_imag_strain = MMP_eps22_imag  doatt = 1 doshear = 0
[../]
[./dispersion]
  type = AttenuationDispersion pp_name_real_stress = MMP_sigma22_real pp_name_real_strain = MMP_eps22_real pp_name_imag_stress  = MMP_sigma22_imag pp_name_imag_strain = MMP_eps22_imag  doatt = 0 doshear = 0
[../]
[]

[Outputs]
 file_base = oldMesh_oldAssembly_oldPostprocessor
 #exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
[]

