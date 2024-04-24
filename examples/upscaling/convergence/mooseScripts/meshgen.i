[Problem]
 kernel_coverage_check = false
 solve = false
[]

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = ${nel} # 100
 ny = ${nel} # 100
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
 fx_string =  '100.0'
 fy_string =  '100.0'
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

[./inclusionRefinement]
type = InclusionRefinement
fractureMeshModifier = inclusions
doBoundaryRefinement = ${dbr}
refinements = '${amr} 0'
outputFileName = mesh_${nel}_${amr}_${dbrString}.xdr
[../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]
