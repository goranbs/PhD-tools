#!/bin/bash
##############################################################################
#
# script: displacemolecules_vmd.sh
#
# dependency: vmd, vmd_dispmolecs.sh
#
# Description:
#
# input molecule files; m1,m2; files must be readable by vmd
# new system size pbc: x,y,z; displace molecules to center of pbc system.
# displace molecule m1, as upper molecule, and rotate by angle; deg.
# Lower molecule m2, is displaced z0 in z-direction.
# Upper molecule m1, is displaced with an geometrical interfacial distance
# z0 from from molecule m2, according to the geometrical thickness of the 
# molecules.
#
# General description:
#
# Read two vmd compatible molecules, find geometrical center and displace
# into new periodic box with dimention x,y,z, and align in z-direction. 
# Displace according to thickness of lower molecule and the distance 
# between molecule1 and molecule2.
# Write output file with degree of rotation and estimated interfacial 
# separtation.
#
#       _______________
#   |   |   ______    |
#   |   |   |    |    |-> COM upper molecule
#   |   |   ------    |--
#   z   |             |  dist
#   |   |   ______    |__
#   |   |   |    |    |-> COM lower molecule
#   |   |   ------    |
#   |   |             |
#       ---------------
#
# TODO/bugs:
#
# 1) Using set sel with the "atomselect" command, returns a selection,
#    this selection has been hard coded into the bash script, and might
#    cause problems if not the same selection number is returned.
#
# 2) Only pdb files has been tested.
#
#                                  Author: Gøran Brekke Svaland 2016-06-22

###############################################################################

#-- Input parameters
m1='calcite_slab.pdb'           # upper molecule
m2='calcite_slab.pdb'           # lower molecule

name='cal'                      # prefix of output filenames

x=80.0                          # final pbc x
y=39.7                          # final pbc y
z=83.1                          # final pbc z

t1=10.39 #17.218 #10.3979	# thickness of m1
t2=10.39 #17.218 #10.3979	# thickness of m2 

z0=25.0                         # displacement of lower molecule, molec2
d=17.0                          # initial interfacial distance between m1, m2
deg=45.0                        # angle of rotation for upper molecule

nsep=3				# number of individual surface-surface separations
res=0.5				# resolution [angstrom]

# ------------------- Perform backups of molecules -------------------------- #
cp ${m1} backup-${m1} # backup of molecule1
cp ${m2} backup-${m2} # backup of molecule2
# --------------------------------------------------------------------------- #

#-- Run vmd_dispmolecs.sh to create systems of z displaced molecules m1 and m2.
#   Decreasing molecule distance separations by i, from initial separation d, 
#   create systems:
#   for i in 0 0.1 0.2 0.3; do

# --------------------------------------------------------------------------- #
sep=${d}                                           # initial separation
fsep=$(echo "scale=5; ${d}+${nsep}*${res} " | bc)  # final separation
i=1                                                # counter 
# --------------------------------------------------------------------------- #


until [ $i -gt ${nsep} ] ; do                
	echo $i $sep
	zm1=$(echo "scale=5; ${z0}+$i" | bc)
        
	./vmd_dispmolecs.sh ${m1} ${m2} $x $y $z ${sep} ${deg} ${name} ${t1} ${t2}

        sep=$(echo "scale=5; ${sep}+${res} " | bc)
        let i+=1

done

#-- EOF
