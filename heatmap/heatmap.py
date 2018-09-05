#!/usr/bin/python
"""
Created on Mon Apr 24 11:46:17 2017

@author: goranbs

This script is perform 2D binning of input data

"""

import numpy as np                         # numpy library
import argparse                            # handle command line arguments

wa_parser = argparse.ArgumentParser(description='Perform 2D binning of dump1d output')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs='+',
                       help='dump1d output file from eldip')
wa_parser.add_argument('-p', '--plot', action='store_true',
                       help='plot result')
wa_parser.add_argument('--val', metavar=('x','y','z'), nargs=3,
                       help='columns in data file for 2d histogramming. Columns start counting on column 1.')                       
wa_parser.add_argument('--nbins', metavar=('n'), nargs=1,
                       help='number of histogram bins used in x,y,z. Default 200')
wa_parser.add_argument('--xlim', metavar=('xlo','xhi'), nargs=2,
                       help='upper and lower limits for values in column x.')
wa_parser.add_argument('--ylim', metavar=('ylo','yhi'), nargs=2,
                       help='upper and lower limits for values in column y.')
wa_parser.add_argument('--zlim', metavar=('zlo','zhi'), nargs=2,
                       help='upper and lower limits for values in column z.')                       
                       
args = wa_parser.parse_args()

# ----------------------------------------------------------------------------

nbins=200

ninputs=np.size(args.inputfile)

if args.val is None:
    print "Warning! No columns provided through --val flag!"
    print "...trying columns 1 2 3"
    a=0
    b=1
    c=2
else:
    a=int(args.val[0])-1 # x (start indexing at 0 in python)
    b=int(args.val[1])-1 # y
    c=int(args.val[2])-1 # z

if args.nbins is None:
    nbins=200
else:
    nbins=int(args.nbins[0])

xlo=float(args.xlim[0])
xhi=float(args.xlim[1])
ylo=float(args.ylim[0])
yhi=float(args.ylim[1])
zlo=float(args.zlim[0])
zhi=float(args.zlim[1])

for i in range(ninputs):

    dump1d = np.loadtxt(args.inputfile[i],comments='#')
    H1, xedges1, yedges1 = np.histogram2d(dump1d[:,a], dump1d[:,b], bins=nbins,range=[[xlo,xhi],[ylo,yhi]])  # xy
    H2, xedges2, yedges2 = np.histogram2d(dump1d[:,a], dump1d[:,c], bins=nbins,range=[[xlo,xhi],[zlo,zhi]])  # xz
    # xz
    out1= 'out_xy_' + args.inputfile[i]
    out2= 'out_xz_' + args.inputfile[i]
    np.savetxt(out1,np.transpose(H1))
    np.savetxt(out2,np.transpose(H2))
    print np.shape(H1), np.shape(xedges1), np.shape(yedges1)
    
    #print xedges1, yedges1
    #print xedges2, yedges2
    
    if (args.plot):
        import matplotlib.pyplot as plt
        #fig, (ax1, ax2) = plt.subplots(2,sharey=True)
        fig = plt.figure()
        
        plt.imshow(H1, interpolation='nearest', origin='low', extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]]) # xy
        plt.xlim((xedges1[0],xedges1[-1]))
        plt.ylim((yedges1[0],yedges1[-1]))
        plt.xlabel('x')
        plt.ylabel('y')
                        
        plt.draw()
        
if (args.plot):
    plt.show()    
