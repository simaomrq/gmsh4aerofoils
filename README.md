# gmsh4aerofoils
Aerofoil meshing using GMSH

Python script to generate meshes for SU2 CFD solver

usage: pyGMSH_foil.py [-h] [--plot PLOT] foilPath

to run code:
	pyGMSH_foil.py crm.eta65.unswept31.5deg.sharp.te.txt 
	
or to plot the aerofoil before generating the mesh:
	pyGMSH_foil.py --plot crm.eta65.unswept31.5deg.sharp.te.txt

output files generated:
		aerofoil.geo_unrolled	-> gmsh script
		mesh_out_s.su2 		-> if sharp t.e.
	or
		mesh_out_b.su2 		-> if blunt t.e.
