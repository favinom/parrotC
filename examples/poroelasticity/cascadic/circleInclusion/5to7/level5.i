[Problem]
 type = FEProblem
 solve = true
 kernel_coverage_check = false
 []

[Mesh] file = '../../../../meshes/circle_0005_mesh.xdr' []

[Variables]
 [./disp_x]   order = FIRST family=LAGRANGE [../]
 [./disp_y]   order = FIRST family=LAGRANGE [../]
 [./pressure] order = FIRST family=LAGRANGE [../]
[]
 
 
[Kernels]
[./StressDivergenceParrot_x]
 type = FreqPoroelasticityFast
 variable = disp_x
 variable_p = pressure
 
 disp_x = disp_x
 disp_y = disp_y
 
 pres = pressure
 
[../]
 []
 
## check if auxkernels work in 3D
 
[Materials]
[./mymat]
 type = FreqPoroelasticSphere
 disp_x = disp_x
 disp_y = disp_y
 x_center = 0.0
 y_center = 0.0
 z_center = 0.0
 radius   = 0.32
 [../]
 []

[BCs]
[./disp_x_right]  type=DirichletBCC variable=disp_x boundary=right value_r=0.0 value_i=0.0 [../]
[./disp_x_left]   type=DirichletBCC variable=disp_x boundary=left  value_r=0.0 value_i=0.0 [../]

[./disp_y_bottom] type=DirichletBCC variable=disp_y boundary=bottom value_r=0.0 value_i=0.0 [../]
[./disp_y_top]    type=DirichletBCC variable=disp_y boundary=top    value_r=-1e-3 value_i=0.0 [../]

[]

 
 
 [Preconditioning]
 [./SMP]
#type = FDP
#type = PBP
type = SMP
 full = true
ksp_norm = default
 [../]
 []
 
[Executioner]
 type=Transient
 solve_type=LINEAR

 start_time=-4.25
 end_time=6.0
 dt = 0.25


line_search = 'none'
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps                '
#[./Quadrature]
# type = GRID
# order = SIXTH
# [../]
[]

 
[Outputs]
 file_base = level5
 exodus = true
 csv = true
 print_linear_residuals = true
 print_perf_log = true
 []

 [MultiApps]
 [./level6]
 type = TransientMultiApp
 app_type = parrotcApp
 execute_on = timestep_end
 input_files = level6.i
 positions = '0.0 0.0 0.0'
 [../]
 []
 
[Transfers]
[./residual_to_sub]
 type = DualProjectionAMR
 direction = to_multiapp
 execute_on = timestep_end
 multi_app = level6
 source_variable = 'disp_x disp_y pressure'
 variable = 'disp_x disp_y pressure'
 dual = false
 add = false
 has_hanging_nodes = true
 [../]
 []

 
[Postprocessors]
[./omega_times_volume] type = MeanMaterialProperty i = 1 j = 1 dostress = 0 doreal = 0 dofreq = 1 pressure = pressure [../]
[./MMP_eps_real] type = MeanMaterialProperty i = 1 j = 1 dostress = 0 doreal = 1 dofreq = 0 pressure = pressure [../]
[./MMP_sigma_real] type = MeanMaterialProperty i = 1 j = 1 dostress = 1 doreal = 1 dofreq = 0 pressure = pressure [../]
[./MMP_eps_imag] type = MeanMaterialProperty i = 1 j = 1 dostress = 0 doreal = 0 dofreq = 0 pressure = pressure [../]
[./MMP_sigma_imag] type = MeanMaterialProperty i = 1 j = 1 dostress = 1 doreal = 0 dofreq = 0 pressure = pressure [../]
[./attenuation]
 type = AttenuationDispersion pp_name_real_stress = MMP_sigma_real pp_name_real_strain = MMP_eps_real pp_name_imag_stress  = MMP_sigma_imag pp_name_imag_strain = MMP_eps_imag  doatt = 1 doshear = 0
[../]
[./dispersion]
type = AttenuationDispersion pp_name_real_stress = MMP_sigma_real pp_name_real_strain = MMP_eps_real pp_name_imag_stress  = MMP_sigma_imag pp_name_imag_strain = MMP_eps_imag  doatt = 0 doshear = 0
[../]
[]
