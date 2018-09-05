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


# -- Plotting output results from fractaldim.py

set terminal x11 

list = system('ls *D*.dat')

set xlabel 'log(1/r)'
set ylabel 'log(N)'

plot for [file in list] file u 6:9 w lp t file

replot 'menger_sponge_0000.dat' u ($6+5.38):9 w lp t 'Menger sponge'
#    EOF
