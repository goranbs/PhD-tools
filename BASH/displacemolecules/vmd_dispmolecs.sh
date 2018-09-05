#!/bin/bash
############################################################################
#
# bash script vmd_displace_molecules.sh
#
# Description:
#
# Read two vmd compatible molecules, find geometrical center and displace
# into new periodic box with dimention x,y,z, and align in z-direction. 
# Displace according to thickness of lower molecule and the distance 
# between molecule1 and molecule2.
# Write output file with degree of rotation and estimated COM separation.
#
#      _______________
#      |   ______    |
#      |   |    |    |-> COM upper molecule
#      |   ------    |--
#      |             |  sep
#      |   ______    |__
#      |   |    |    |-> COM lower molecule
#      |   ------    |
#      |             |
#      ---------------
#
#                                  Author: Gøran Brekke Svaland 2016-06-22
###########################################################################

args=("$@")          # input arguments

# ---------------- Set values from input args --------------------------- #

molecule1=${args[0]} # upper molecule m1
molecule2=${args[1]} # lower molecule m2

x=${args[2]}         # final pbc x
y=${args[3]}         # final pbc y
z=${args[4]}         # final pbc z

sep=${args[5]}       # interfacial separation
deg=${args[6]}       # angle of horizontal (z) rotation for upper molecule
name=${args[7]}      # prefix of output filenames
t1=${args[8]}        # thickness of m1
t2=${args[9]}       # thickness of m2

#-- midpoints of pbc system
xm=$(echo "scale=5; $x/2" | bc)  # x centering of slabs
ym=$(echo "scale=5; $y/2" | bc)  # y centering of slabs
zm=$(echo "scale=5; $z/2" | bc)  # z centering of slabs


#-- output name of new pdb configuration
output="${name}_rot${deg}_sep${sep}" # acc. to rotation and separation

#-- print input info
echo "molecule1: ${molecule1}"      # m2
echo "molecule2: ${molecule2}"      # m1
echo "x   =$x"                      # final pbc x
echo "y   =$y"                      # final pbc y
echo "z   =$z"                      # final pbc z
echo "COMd=${sep}"                  # interfacial sep
echo "deg =${deg}"                  # angle rot around y
echo "t1  =${t1}"                   # thickness of molecule1
echo "t2  =${t2}"                   # thickness of molecule2

# ------------- create readable file for vmd  ------------- #

#-- temporary file
outfile="tmp.vmd"  # temporary input file for vmd

# --------------------------------------------------------- #
#-- proc for computing the geometrical center of a selection
#   commented out because I found that this method were
#   already implemented in VMD; measure center $selection
#echo -e 'proc geom_center {selection} {
#        set gc [veczero]
#        foreach coord [$selection get {x y z}] {
#           set gc [vecadd $gc $coord]
#        }
#        return [vecscale [expr 1.0 /[$selection num]] $gc]
#}' >> ${outfile}
# You may still use this metod with the call:
# set Rg [geom_center atomselect0]
# lassign $C a b c # assign elements in C to variables a,b,c
# --------------------------------------------------------- #

#-- load molecule, rotate an angle deg around z axis, and move
#-- to center of pbc box:
echo "mol new ${molecule1}" >> ${outfile}
echo "set sel [atomselect 0 all]" >> ${outfile}
echo "set C [measure center atomselect0]" >> ${outfile}
echo "lassign $C a b c" >> ${outfile}
echo "set A [pbc get]" >> ${outfile} # could alternatively be used? If I had managed...
echo 'lassign $A x y z j k l' >> ${outfile}
echo 'atomselect0 moveby [vecinvert $C]' >> ${outfile}
echo "atomselect0 move [transaxis z ${deg}]" >> ${outfile}

#-- load initial molecule to origin, and move into pbc box:
echo "mol new ${molecule2}" >> ${outfile}
echo "set B [pbc get]" >> ${outfile}
echo 'lassign $B xx yy zz jj kk ll' ${outfile}
echo "pbc set {${x} ${y} ${z}}" >> ${outfile}
echo "set sel [atomselect 1 all]" >> ${outfile}
echo "set D [measure center atomselect2]" >> ${outfile}
echo 'lassign $D d e f' >> ${outfile}
echo 'atomselect2 moveby [vecinvert $D]' >> ${outfile}
echo "atomselect2 moveby {0 0 ${sep}}" >> ${outfile}
echo "atomselect2 moveby {0 0 ${t1}}" >> ${outfile}

#-- compute final geometrical positions of COM
echo "measure center atomselect0" >> ${outfile}
echo "measure center atomselect2" >> ${outfile}

#-- create new atomselection containing the two molecules,
#-- and move to center of simulation box:
echo "set kk {}" >> ${outfile}
echo "lappend kk atomselect0 atomselect2" >> ${outfile}
echo 'set mol [::TopoTools::selections2mol $kk]' >> ${outfile}
echo "set sel [atomselect 2 all]" >> ${outfile}
half_t1=$(echo "scale=5; ${t1}/2.0 " | bc)
half_t2=$(echo "scale=5; ${t2}/2.0 " | bc)
z0=$(echo "scale=5; (${zm}-${t1}-${half_t2})" | bc)
echo ${half_t1} ${z0}
echo "atomselect10 moveby {${xm} ${ym} ${half_t1}}" >> ${outfile}
echo "atomselect10 moveby {0 0 ${z0}}" >> ${outfile}
echo "mol delete 0" >> ${outfile}
echo "mol delete 1" >> ${outfile}
echo "pbc set {${x} ${y} ${z}}" >> ${outfile}

#-- write output file:
echo "animate write pdb ${output}.pdb 2" >> ${outfile}
#echo "topo writelammpsdata ${output}.data full" >> ${outfile}
echo exit >> ${outfile}

vmd -e ${outfile} # run vmd with commands in ${outfile}
rm ${outfile}     # delete temporary file: ${outfile}


#-- EOF
