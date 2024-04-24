[Problem]
 kernel_coverage_check = false
 solve = false
[]

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 40 # 100
 ny = 40 # 100
 xmin = 0
 xmax = 200
 ymin = 0
 ymax = 200
# If you want to use quadratic elements, regenerate the mesh with the following line uncommented.
# elem_type = QUAD9
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

[./inclusionRefinement]
type = InclusionRefinement
fractureMeshModifier = inclusions
doBoundaryRefinement = true
refinements = '6 0'
outputFileName = mesh.xdr
[../]

[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]
