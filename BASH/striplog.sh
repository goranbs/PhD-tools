#!/bin/bash

# -------------------------------------------------------------------- #
#
# Give LAMMPS log file as argument, return .out file containing only
# section between "Step" and "Loop time" in the thermo output.
#
# -------------------------------------------------------------------- #

args=("$@")

file=${args[0]}

# find first line matching "Step" then print all subsequent lines;
sed -i.bak -n '/Step/,$p' $file

# print all lines from the first until "Loop time" is found;
sed -i -n '1,/Loop time/ p' $file

# print all lines exept the last;
head -n -1 $file > tmp.dat

# print all starting from line nr 2;
tail -n +2 tmp.dat > ${file}.out

rm tmp.dat

# EOF
