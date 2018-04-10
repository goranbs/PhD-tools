#!/usr/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 6    last modified September 2014
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2014
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')

#set terminal x11 
set terminal postscript eps enhanced color font 'Helvetica,18' linewidth 4
set output 'Dimensions.eps'

set key top left
set xlabel 'log(1/r)'
set ylabel 'log(N)'

f(x,a) = a*x + b
a=1

fit [0:2.5] f(x,a) '1D_0000.dat' u 6:9 via a,b
D1 = a
fit [0:2.5] f(x,a) '2D_0000.dat' u 6:9 via a,b
D2 = a
fit [0:2.5] f(x,a) '3D_0000.dat' u 6:9 via a,b
D3 = a

plot	'1D_0000.dat' u 6:9 w lp t sprintf('1D Df=%.2f', D1), \
	'2D_0000.dat' u 6:9 w lp t sprintf('2D Df=%.2f', D2), \
	'3D_0000.dat' u 6:9 w lp t sprintf('3D Df=%.2f', D3), \
	f(x,D1) lt 0 not, \
	f(x,D2) lt 0 not, \
	f(x,D3) lt 0 not


#    EOF
