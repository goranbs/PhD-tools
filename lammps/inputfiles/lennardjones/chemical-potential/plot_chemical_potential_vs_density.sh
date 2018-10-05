#!/bin/bash

# -- Create directories and equilibrate LJ-fluid at different densities using the LAMMPS input script: in.eq

T=1.3265	# temperature
rho_i=0.25	# lowest density
rho_f=0.35	# highest density
N=10		# number of points
i=0

# ------------- #
filename="mu_vs_rho.dat"
# ---------------------- #

echo "# chemical potential (mu) vs density (rho) at T=${T}" > ${filename}
echo "# rho mu_up mu_down" >> ${filename}

while [ $i -le $N ] ; do
	d=$(bc <<< "scale=2;($i*(${rho_f}-${rho_i})/($N))")
	s=$(bc <<< "scale=2;(${rho_i}+$d)")
	echo $s

	u=$(fep ${T} < rho${s}/fep01-1.lmp)
#        v=$(fep ${T} < rho${s}/fep10-1.lmp)

	
	echo "${s} ${u}" >> ${filename}
	
	# --------------------------------------------------


	let i+=1
done

echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'mu_vs_rho.eps'

set encoding iso_8859_1

#set title 'T*=${T}'
set xlabel 'rho*'
set ylabel '{/Symbol m}*'

plot 'mu_vs_rho.dat' u 1:2 w lp t 'T*=${T}'

" > mu_vs_rho.gnu

#gnuplot mu_vs_rho.gnu

