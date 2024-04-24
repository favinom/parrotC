[Problem]
 kernel_coverage_check = false
 solve = false
[]

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 100
 ny = 100
 xmin = 0
 xmax = 100
 ymin = 0
 ymax = 100
#If you want to use quadratic elements, regenerate the mesh with the following line uncommented.
#elem_type = QUAD9
#partitioner = parmetis
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

[./inclusionRefinement]
type = InclusionRefinement
fractureMeshModifier = inclusions
doBoundaryRefinement = true
refinements = '11 0'
outputFileName = mesh.xdr
[../]

[]


[Executioner]
 type=Steady
 solve_type=LINEAR
[]
