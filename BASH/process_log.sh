#!/bin/bash
# --------------------------------------------------------------------------------------- #
# -------------------------------- process_log.sh --------------------------------------- #
#
# USAGE:
#		>> bash process_log.sh filename [log.lammps] terminal [x11,eps]
#
#
# DESCRIPTION:
#
# 1) Strip LAMMPS log file to make it useful for awk
# 2) Perform rolling averages of data fields
# 3) Generate gnuplot input files
# 	a) If terminal is set to "postscript"; create EPS files of graphs
#	b) If terminal is set to "x11"       ; do nothing 
#	c) Default terminal: x11
#
# DEPENDENCIES:
#		kk.awk
#		striplog.sh
#
# --------------------------------------------------------------------------------------- #

args=("$@")

file=${args[0]}		# input LAMMPS log file
term=${args[1]}		# optional terminal type: 1) x11 2) ps,eps,postscript

bash striplog.sh $file

awk -f kk.awk ${file}.out > ${file}.ravg

rm ${file}.out		# delete tmp file created by striplog.sh

# --------------------------------------------------------------------------------------- #
# ---------------------- WRITE GNUPLOT FILES -------------------------------------------- #

# -- SET TERMINAL:

PS="PS EPS ps eps postscript"
X11="X11 x11 x-11 X-11"

terminal="set terminal x11" 	# default
for i in $PS ; do
	#echo $i
	if [ "$i" == "$term" ] ; then
		term="eps"
		terminal="set terminal postscript eps enhanced color font 'Helvetica,24' linewidth 3"
		echo "terminal = postscript"
	fi
done

for i in $X11 ; do
	#echo $i
	if [ "$i" == "$term" ] ; then
		terminal="set terminal x11"
		echo "terminal = x11"
	fi
done
	
# --------------------------------------------------------------------------------------- #

#-- energy-fe-fse;
echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'Energy-fe-fse.eps'

set title 'Energy of slab-slab and slab-water interactions'
set xlabel 'time [ns]'
set ylabel 'energy [10^3 kJ/mol]'

dt=2.0			# time step
ns=dt*1.0e-6		# fsec to ns
MkJ=4.184*1.0e-3	# kcal to 10^6 kJoule

f(x) = mean_a
g(x) = mean_b
fit [0.5:*] f(x) '${file}.ravg' u (\$2*ns):(\$7*MkJ) via mean_a
fit [0.5:*] g(x) '${file}.ravg' u (\$2*ns):(\$9*MkJ) via mean_b
set key outside center above

plot 	'${file}.ravg' u (\$2*ns):(\$7*MkJ) w l t gprintf(\"f_e=%.3f 10^3 kJ/mol\",mean_a), \
	'${file}.ravg' u (\$2*ns):(\$9*MkJ) w l t gprintf(\"fs_e=%.3f 10^3 kJ/mol\",mean_b)

" > energy-fe-fse.gnu

#-- forces;
echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'Forces.eps'

set title 'Forces'
set xlabel 'time [ns]'
set ylabel 'f/A [MPa]'

Acal = 1710.282084      	# surface area of calcite slab [angstrom^2]
MPa = 6947.695295538502/Acal	# kcal/mol/angstrom to MPa
dt=2.0				# time step
ns=dt*1.0e-6			# fsec to ns

f(x) = mean_a
g(x) = mean_b
fit [0.5:*] f(x) '${file}.ravg' u (\$2*ns):(\$8*MPa) via mean_a
fit [0.5:*] g(x) '${file}.ravg' u (\$2*ns):(\$10*MPa) via mean_b
set key outside center above

plot 	'${file}.ravg' u (\$2*ns):(\$8*MPa) w l t gprintf(\"f_z=%.3f MPa\",mean_a), \
	'${file}.ravg' u (\$2*ns):(\$10*MPa) w l t gprintf(\"fs_z=%.3f MPa\",mean_b)
" > forces.gnu

#-- nconf;
echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'Nconf.eps'

set title '# confined water molecules'
set xlabel 'time [ns]'
set ylabel '# molecules'

dt=2.0				# time step
ns=dt*1.0e-6			# fsec to ns

f(x) = mean_a
fit [0.5:*] f(x) '${file}.ravg' u (\$2*ns):4 via mean_a

plot 	'${file}.ravg' u (\$2*ns):4 w l t gprintf(\"N=%.3f\",mean_a)
" > nconf.gnu

#-- pe;
echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'Pe.eps'

set title 'Potential energy of system'
set xlabel 'time [ns]'
set ylabel 'energy [10^6 kJ/mol]'

dt=2.0				# time step
ns=dt*1.0e-6			# fsec to ns
MkJ=4.184*1.0e-6		# kcal to 10^6 kJoule

plot 	'${file}.ravg' u (\$2*ns):(\$3*MkJ) w l t 'pe'
" > pe.gnu

#-- upper-lower-COM;
echo "
#!/usr/bin/gnuplot -persist

${terminal}
set output 'Upper-lower-COM.eps'
set encoding iso_8859_1

set title 'Center of mass'
set xlabel 'time [ns]'
set ylabel 'z [{\305}]'

dt=2.0				# time step
ns=dt*1.0e-6			# fsec to ns

f(x) = mean_a
g(x) = mean_b
fit [0.5:*] f(x) '${file}.ravg' u (\$2*ns):5 via mean_a
fit [0.5:*] g(x) '${file}.ravg' u (\$2*ns):6 via mean_b
set key outside center above

plot 	'${file}.ravg' u (\$2*ns):5 w l t gprintf(\"lower slab z=%.3f {\305}\",mean_a), \
	'${file}.ravg' u (\$2*ns):6 w l t gprintf(\"upper slab z=%.3f {\305}\",mean_b)
" > upper-lower-COM.gnu

# --------------------------------------------------------------------------------------- #
# ---------------------- CREATE EPS of graphs ------------------------------------------- #

if [ "$term" == "eps" ] ; then
	gnuplot energy-fe-fse.gnu
	gnuplot forces.gnu
	gnuplot nconf.gnu
	gnuplot pe.gnu
	gnuplot upper-lower-COM.gnu
fi

# EOF

