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

#FILES=./aragonites/*.pdb     # path/to/files/*.gro  (.gro or .pdb)
#FILES=cal_rot0_sep27.pdb
#FILES=big_cal_rot0_sep36.pdb
#FILES=./calcites/*.pdb

#FILES=./gold-nocharge/*.pdb
#FILES=./gold-q22x/*.pdb
FILES=./gold-confinement-charge/*.pdb

#x=8                        # final box dimension x [nm]
#y=3.9853                   # final box dimension y [nm]  
#z=18                       # final box dimension z [nm]
#y=4.0721
#z=20
x=5.00000 # [nm]
y=2.44692 # [nm]
z=15.0000 # [nm]
####################################################################################

for file in $FILES ; do
	./solvate.sh $file $x $y $z
done
##-- EOF
