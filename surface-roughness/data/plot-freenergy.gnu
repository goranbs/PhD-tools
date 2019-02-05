
set terminal x11

set xlabel "h [{\305}]"
set ylabel "F [mJ/m^2]"

#########################################################################
Joule=4.1868		# 1 calorie = 4.1868 Joule
Acal=1710.282084	# surface area calcite slab [angstrom^2]
Aara=1673.826		# surface area aragonite slab [angstrom^2]
Na = 6.022		# avogadros constant [*10^23]

# convert from kcal/mol to mJ/m^2
# kJ/mol * 1/angstrom^2 * 1/Na = 10^3 * 1/Na * 1/A [J/m^2]

fcal = 1000*Joule/Na/Acal	# conversion factor calcite

# MPa*m = 10^6Joule
# MPa*angstrom = 10^6 * N/m^2 * 10^(-10) m = 0.1 mJ/m^2

fgrace = 0.1
#########################################################################

oCal  = 15.26      # offset from COM value at shortest sep.
shift = 169.79975  # offset from zero at large separations

plot 	'calcite_inregistry.dat' u ($1-oCal):(fcal*$6) w lp t 'calcite_{in}', \
	'grace-out-freenergy.dat' u 1:(fgrace*(-$2-shift)) w lp t 'grace'

#unset terminal
#set table 'freenergy.txt'
#plot 'grace-out-freenergy.dat' u 1:(fgrace*(-$2-shift)) w lp t 'grace'

#-- EOF
