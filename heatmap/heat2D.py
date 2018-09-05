#!/usr/bin/python
"""
Created on Mon Apr 24 11:46:17 2017

@author: goranbs

This script is perform 2D binning of input data

"""

import numpy as np                         # numpy library
import argparse                            # handle command line arguments
import matplotlib.pyplot as plt

wa_parser = argparse.ArgumentParser(description='Perform 2D histogramming')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='input file with columnar data')
wa_parser.add_argument('-o', '--output', metavar=('prefix'), type=str, nargs=1,
                       help='output prefix')                       
wa_parser.add_argument('--val', metavar=('x','y'), nargs=2,
                       help='columns in data file for 2d histogramming. Columns start counting on column 1.')                       
wa_parser.add_argument('--xlim', metavar=('xlo','xhi'), nargs=2,
                       help='upper and lower limits for values in column x. Default: evaluate limits from columnar data.')
wa_parser.add_argument('--ylim', metavar=('ylo','yhi'), nargs=2,
                       help='upper and lower limits for values in column y. Default: evaluate limits from columnar data.')                       
wa_parser.add_argument('-n', '--nbins', metavar=('n'), nargs=1,
                       help='number of histogram bins used in x,y. Default: 1000')                       
wa_parser.add_argument('-s', '--scale', metavar=('s'), nargs=1,
                       help='scale box counts by s. Default: 1')
wa_parser.add_argument('-c', '--cmax', metavar=('c_hi'), nargs=1,
                       help='set upper limit for colour box. Default: max_count_per_box*scale')                       
wa_parser.add_argument('-p', '--plot', action='store_true',
                       help='show plot. Default: do not show')                  
wa_parser.add_argument('--xlabel', metavar=('string'), nargs=1,
                       help='print this xlabel')
wa_parser.add_argument('--ylabel', metavar=('string'), nargs=1,
                       help='print this ylabel')
wa_parser.add_argument('--clabel', metavar=('string'), nargs=1,
                       help='print this colour label')
wa_parser.add_argument('--title', metavar=('string'), nargs=1,
                       help='print this title')                       
args = wa_parser.parse_args()

# ----------------------------------------------------------------------------

# -- Evaluate which columns to use from input data file
if args.val is None:
    print "Warning! No columns provided through --val flag!"
    print "...trying columns 1 2 3"
    x=0
    y=1
else:
    x=int(args.val[0])-1 # x (start indexing at 0 in python)
    y=int(args.val[1])-1 # y

# -- Number of bins to use for histogram
if args.nbins is None:
    nbins=1000
else:
    nbins=int(args.nbins[0])

# -- Read data file columns and evaluate histogram limits
data = np.loadtxt(args.inputfile[0],comments='#')

if (args.xlim) is None:
    xmin=float(min(data[:,x]))
    xmax=float(max(data[:,x]))
else:    
    xmin=float(args.xlim[0])
    xmax=float(args.xlim[1])
    
if (args.ylim) is None:
    ymin=float(min(data[:,y]))
    ymax=float(max(data[:,y]))
else:
    ymin=float(args.ylim[0])
    ymax=float(args.ylim[1])

# -- Scale number of counts by value
if args.scale is None:
    scale=1.0
else:
    scale=float(args.scale[0])
    

# -- Labels
if args.xlabel is None:
    xlabel='x'
else:
    xlabel=args.xlabel[0]

if args.ylabel is None:
    ylabel='y'
else:
    ylabel=args.ylabel[0]

if args.clabel is None:
    clabel='counts'
else:
    clabel=args.clabel[0]

# -- histogram bins
xbins=np.linspace(xmin,xmax,nbins)
ybins=np.linspace(ymin,ymax,nbins)
xcenter = (xbins[0:-1]+xbins[1:])/2.0
ycenter = (ybins[0:-1]+ybins[1:])/2.0
#aspectratio = 1.0*(xmax - 0)/(1.0*ymax - 0)

H, xedges, yedges = np.histogram2d(data[:,x], data[:,y], bins=(xbins,ybins))  
  
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
Sum=np.sum(H)   # number_of_data_points
Max=np.amax(H)  # max_count_per_box
#print 'Number of data points : ', Sum
#print 'Max count in one box  : ', Max
#print 'Scale value           : ', scale

if args.cmax is None:
    cmax=Max*scale
else:
    cmax=float(args.cmax[0])
    
#print 'Colour bar limits     :  0,', cmax
#print 'X-axis limits         : ', xmin, xmax
#print 'Y-axis limits         : ', ymin, ymax
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

# -- Plotting the
add_background_img=True
if add_background_img:
    background="c1-surf-x20-60-y0-40.png"
    im = plt.imread(background)
    plt.imshow(im, extent=[xmin,xmax,ymin,ymax])
    
plt.imshow(H*scale, extent=[xmin,xmax,ymin,ymax],
       interpolation='nearest', origin='lower',aspect='equal',alpha=0.8)

cbox = plt.colorbar()
cbox.set_label(clabel)
plt.clim(vmin=0,vmax=cmax)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
if args.title is not None:
    plt.title(args.title[0])

if (args.plot):
    plt.show()

if args.output is None:
    fname='out-{}.eps'.format(args.inputfile[0])
else:
    fname='{}-{}.eps'.format(args.output[0],args.inputfile[0])
    
plt.savefig(fname,format='eps')
      
# ----------------------------------------------------------- EOF    
    