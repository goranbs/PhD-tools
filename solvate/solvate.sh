#!/bin/bash
############################################################################
#                               SOLVATE
# 
# Solvate a .pdb or .gro file using GROMACS' genbox command.
# Use GROMACS' editconf command to create a vacuum. (that is, the provided
# x,y,z dimentions must be larger or equal to the original system size)
# Output a solvated system in .pdb format, readable for VMD
# Create temporary input file for VMD.
# Call VMD with this input file in order to output a .data file for LAMMPS.
############################################################################

args=("$@")               # input arguments

##-- input arguments:
ifile=${args[0]}          # input file (.pdb or .gro)
x=${args[1]}              # size of simulation box in x dir (nm)
y=${args[2]}              # -//- y
z=${args[3]}              # -//- z

##-- default values:
vdwd=0.11                        # van der waals distance between solute and solvate
solvent="spc216.gro"             # solvent
vmdfile=${ifile}'-Solvated.pdb'  # input file for VMD
lammpsfile=${ifile}'.data'       # lammps data file: final output
outfile='tmp.vmd'                # temporary input file for VMD
tmpfile1="tmp.gro"               # temporary output: .gro file
tmpfile2="tmp.gro"               # temporary output: solvated system

##-- make sure it is a .gro file:
editconf -f ${ifile} -o ${tmpfile1}
##-- solvate the system:
genbox -vdwd ${vdwd} -cp ${tmpfile1} -cs ${solvent} -o ${tmpfile2}
##-- change system box size:
editconf -box $x $y $z -f ${tmpfile2} -o ${vmdfile}
rm ${tmpfile1} ${tmpfile2}

# -----------------------------------------------------------

echo "mol new ${vmdfile}" >> ${outfile}
echo "topo guessangles -sel waters" >> ${outfile}
#echo 'set sel [atomselect top "resname SRF"]'
#echo "atomselect0 set charge 0.0509259259259259"
echo "topo writelammpsdata ${lammpsfile} full" >> ${outfile}
echo "quit" >> ${outfile}

vmd -e ${outfile}  # run VMD with the generated temporary file
rm ${outfile}

echo "input file           : " ${ifile}
echo "pdb solvated system  : " ${vmdfile}
echo "lammps data file     : " ${lammpsfile}
#echo "tmpfile1  : " ${tmpfile1}
#echo "tmpfile2  : " ${tmpfile2}
##-- EOF
