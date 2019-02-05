
set terminal x11

set ylabel 'disjoining pressure [MPa]'
set xlabel 'd [angstrom]'

set pointsize 2

plot 	'calcite_inregistry.dat' u ($1-15.3):4 w l lc rgb 'red' t 'original data', \
	'disjoiningpressure.txt' u 1:2 w p pt 4 lc rgb 'blue' t 'grabit points'
