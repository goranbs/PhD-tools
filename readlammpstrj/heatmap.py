#!/usr/bin/python
"""
Created on Mon Apr 24 11:46:17 2017

@author: goranbs

Read dump1d output from eldip (neweldip) and write heatmap matrix

dump1d performs accumulation of dipole moments.
This simple script is supposed to perform simple 2d binning of this data
and output a resulting 2d heatmap.

"""

import numpy as np                         # numpy library
import argparse                            # handle command line arguments
wa_parser = argparse.ArgumentParser(description='Perform 2D binning of dump1d output')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs='+',
                       help='dump1d output file from eldip')
wa_parser.add_argument('-p', '--plot', action='store_true',
                       help='plot result')
wa_parser.add_argument('--val', metavar=('x','y','z'), nargs=3,
                       help='Indices in data file for 2d histogramming.')                       

args = wa_parser.parse_args()

# ----------------------------------------------------------------------------

nbins=200

ninputs=np.size(args.inputfile)
a=int(args.val[0]) # x
b=int(args.val[1]) # y
c=int(args.val[2]) # z

for i in range(ninputs):

    dump1d = np.loadtxt(args.inputfile[i],comments='#')
    H1, xedges1, yedges1 = np.histogram2d(dump1d[:,a], dump1d[:,b], bins=nbins)  # xy
    H2, xedges2, yedges2 = np.histogram2d(dump1d[:,a], dump1d[:,c], bins=nbins)  # xz
    out1= 'out_xy_' + args.inputfile[i]
    out2= 'out_xz_' + args.inputfile[i]
    np.savetxt(out1,np.transpose(H1))
    np.savetxt(out2,np.transpose(H2))
    
    if (args.plot):
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2,sharey=True)
        ax1.imshow(H1, interpolation='nearest', origin='low', extent=[0, nbins, 0, nbins]) # xy
        ax2.imshow(H2, interpolation='nearest', origin='low', extent=[0, nbins, 0, nbins]) # xz
        ax1.set_ylabel('x')
        ax1.set_xlabel('y')
        ax2.set_ylabel('x') 
        ax2.set_xlabel('z')
        plt.draw()
        
if (args.plot):
    plt.show()    
    
# ---------------------------------------------- Scatter plotting:
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax = fig.add_subplot(111, projection='3d')
#for i in range(ninputs):
#   dump1d = np.loadtxt(args.inputfile[i],comments='#')   
    #plt.xlim(-0.5,0.5)
    #plt.ylim(-0.5,0.5)
    #ax.scatter(dump1d[:,a], dump1d[:,c], alpha=0.1)
    #circle = plt.Circle((0,0), 0.5, color='r', fill=False)
    #ax.add_artist(circle)
    #ax.scatter(dump1d[:,a], dump1d[:,b], dump1d[:,c], alpha=0.1)
# EOF
