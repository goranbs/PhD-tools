#!/bin/bash

mol="Ca"
sep=25

dF_4="dF-res4-${mol}-${sep}.dat"
dF_7="dF-res7-${mol}-${sep}.dat"
dF_8="dF-res8-${mol}-${sep}.dat"

rho_file="../data/spce-res4-${mol}-${sep}.profile"


shift_dF_4=$(awk 'NR==2 {print $5}' ${dF_4})
shift_dF_7=$(awk 'NR==2 {print $5}' ${dF_7})
shift_dF_8=$(awk 'NR==2 {print $5}' ${dF_8})


shift=$(echo "scale=5; (${shift_dF_4}+${shift_dF_7}+${shift_dF_8})/3.0" | bc) 
#echo $shift

echo "
#set terminal x11
#set terminal pngcairo enhanced color font 'Helvetica,18' linewidth 4 size 1000,1000
#set output 'Plot-${mol}-${sep}.png'
set terminal postscript eps enhanced color font 'Helvetica,18' linewidth 4
set output 'Plot-${mol}-${sep}.eps'
set encoding iso_8859_1

shift=${shift}

set ylabel 'dF [k_BT]'
set y2label '{/Symbol r} [g/cm^3]'
set xlabel 'z [{\305}]'

set xrange [0:11]
set yrange [-5:5]
set y2range [0:3.2]

set xtics 0,1,11
set ytics -5,1,5 nomirror
set y2tics 0,1,3 nomirror
set key at 2.5,4.7

PT1=4
PT2=8
PT3=10

plot 	'${dF_4}' u 2:3 w lp lc rgb 'red' pt PT1 t '1' axis x1y1, \\
	'${dF_7}' u 2:3 w lp lc rgb 'coral' pt PT2 t '2' axis x1y1, \\
        '${dF_8}' u 2:3 w lp lc rgb '#DC143C' pt PT3 t '3' axis x1y1, \\
	'${rho_file}' u (\$2-shift):(\$4*1.9) w l lc rgb 'blue' t '{/Symbol r}_w' axis x1y2


" > plot-${mol}-${sep}.gnu

gnuplot plot-${mol}-${sep}.gnu

