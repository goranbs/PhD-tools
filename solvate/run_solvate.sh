#!/bin/bash
####################################################################################
#
# script: run_solvate.sh
#
# dependency: vmd, GROMACS 4.6.5, solvate.sh, spc216.gro
#
# Description:
#
# Solvate system(s) with SPC water using GROMACS
# Guess angles in water molecules using TopoTools
# Dump LAMMPS data file
####################################################################################
#
# Change file path and final dimention of output file:
# NB: x,y,z must be larger or equal to the system dimenstion of the input file!

FILES=./calcites/*.pdb     # path/to/files/*.gro  (.gro or .pdb)
x=8                        # final box dimension x [nm]
y=3.9827                   # final box dimension y [nm]  
z=18                       # final box dimension z [nm]

####################################################################################

for file in $FILES ; do
	./solvate.sh $file $x $y $z
done
##-- EOF
