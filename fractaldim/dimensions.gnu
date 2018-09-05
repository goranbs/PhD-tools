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

set terminal x11 
#set terminal postscript eps enhanced color font 'Helvetica,18' linewidth 4
#set output 'Dimensions.eps'

set key top left spacing 1.2
set xlabel 'log(1/r)'
set ylabel 'log(N)'

f(x,a,b) = a*x + b
a=1
b=0

shiftM = 5.4
shifthM = 6.3

fit [0:2.5] f(x,a,b) '1D_0000.dat' u 6:9 via a,b
D1 = a
b1 = b
stdD1 = 0 #sqrt(FIT_WSSR / (FIT_NDF + 1 ))

fit [0:2.5] f(x,a,b) '2D_0000.dat' u 6:9 via a,b
D2 = a
b2=b
stdD2 = 0 #sqrt(FIT_WSSR / (FIT_NDF + 1 ))

fit [0:2.5] f(x,a,b) '3D_0000.dat' u 6:9 via a,b
D3 = a
b3=b
stdD3 = 0 #sqrt(FIT_WSSR / (FIT_NDF + 1 ))

a=2
b=1
fit [2:3] f(x,a,b) 'menger_sponge_0000.dat' u ($6+shiftM):9 via a,b
Ms1 = a
bs1 = b
stdMs1 = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

fit [2:4] f(x,a,b) 'huge-menger-1-n100_0000.dat' u ($6+shifthM):9 via a,b
Ms2 = a
bs2 = b
stdMs2 = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

fit [2:4] f(x,a,b) 'huge-menger-2-n100_0000.dat' u ($6+shifthM):9 via a,b
Ms3 = a
bs3 = b
stdMs3 = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

plot	'1D_0000.dat' u 6:9 w lp t sprintf('1D Df=%.2f', D1), \
	'2D_0000.dat' u 6:9 w lp t sprintf('2D Df=%.2f', D2), \
	'3D_0000.dat' u 6:9 w lp t sprintf('3D Df=%.2f', D3), \
	'menger_sponge_0000.dat' u ($6+shiftM):9 w lp t sprintf('Ms_1 Df=%.2f', Ms1), \
	'huge-menger-1-n100_0000.dat' u ($6+shifthM):9 w lp t sprintf('Ms_2 Df=%.2f', Ms2), \
	'huge-menger-2-n100_0000.dat' u ($6+shifthM):9 w lp t sprintf('Ms_3 Df=%.2f', Ms3), \
	f(x,D1,b1) lt 0 not, \
	f(x,D2,b2) lt 0 not, \
	f(x,D3,b3) lt 0 not, \
        f(x,Ms1,bs1) lt 0 not, \
        f(x,Ms2,bs2) lt 0 not, \
        f(x,Ms3,bs3) lt 0 not


set print 'Dimensions.dat' append
print "# Data_name alpha std_alpha beta"
print "D1 ", D1, stdD1, b1
print "D2 ", D2, stdD2, b2
print "D3 ", D3, stdD3, b3
print "Ms1 ", Ms1, stdMs1, bs1
print "Ms2 ", Ms2, stdMs2, bs2
print "Ms3 ", Ms3, stdMs3, bs3


#    EOF
