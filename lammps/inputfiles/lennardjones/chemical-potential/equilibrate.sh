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

	mkdir rho${s}
	cd rho${s}
	cp ../in.eq .

	sed -i -e '15c\''variable rho equal '${s}'' in.eq
        sed -i -e '16c\''variable T equal '${T}'' in.eq

	mpirun -n 12 lammps -in in.eq

	# --------------------------------------------------

	cd ..

	let i+=1
done
