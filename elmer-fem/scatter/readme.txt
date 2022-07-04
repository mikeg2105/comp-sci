
Helmholtz scattering problem
box with barrier and two holes


create mesh with blender 
split components into separate parts and import each part as .stl

read into gmsh and save as gmsh file

use elmer grid to convert each gmsh file

 ElmerGrid 14 2 wavebox1a -autoclean
 ElmerGrid 14 2 wavebox1b -autoclean


Assign each part a material type

 ElmerGrid 2 2 wavebox1a -boundorder -bulktype 1 17 1 -out wb1a
 ElmerGrid 2 2 wavebox1b -boundorder -bulktype 1 13 2 -out wb1b


unite the parts together (now that each has been given a different label)

 ElmerGrid 2 2 wb1a -in wb1b -out finalmesh -unite -merge 1.0e-10
