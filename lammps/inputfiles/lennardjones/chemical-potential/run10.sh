#!/bin/bash

# -- Create directories and equilibrate LJ-fluid at different densities using the LAMMPS input script: in.eq

T=1.3265	# temperature
rho_i=0.25	# lowest density
rho_f=0.35	# highest density
N=10		# number of points
i=0

# ------------- #

while [ $i -le $N ] ; do
	d=$(bc <<< "scale=2;($i*(${rho_f}-${rho_i})/($N))")
	s=$(bc <<< "scale=2;(${rho_i}+$d)")
	echo $s

	cd rho${s}
	cp ../in.run10 .

	sed -i -e '18c\''variable T equal '${T}'' in.run01
	sed -i -e '23c\''read_data eq.data extra/atom/types 1' in.run10

	mpirun -n 12 lammps -in in.run10

	# --------------------------------------------------

	cd ..

	let i+=1
done
