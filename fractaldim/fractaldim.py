#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 15:24:28 2018

@author: goranbs (Gøran Brekke Svaland)

fractaldim.py
    computes the number of boxes N, occupied by particles, based on their coordinates.
    
Output:
    prefix_XXXX.dat [chunkNr,rx,ry,rz,N,log(1/rx),log(1/rx),log(1/rx),log(N)]
        this file is output for each frame in the trajectory, and can be used to 
        compute the fractal dimension of the system: Df*log(N) ~ log(1/r), 
        where Df is the fractal dimension of the system, by evaluating the slope
        of log(N) vs. log(1/r).
        
        If the coordinates of the particles in the LAMMPS trajectory is in 
        fractional coordinates [xs,ys,zs], fractional coordinates will be used.
        If the coords are unscaled [x,y,z], unscaled coordinates will be used.
        Unwrapped coordinates [xsu, ysu, zsu], are not supported.

Dependencies:
    * readlammpstrj.py to read particle positions
    * numpy
    * argparse
    * sys

* Triclinic box shapes are not supported!
* Unwrapped coordinates are not supported!
* Trajectories up to 9999 frames are supported

"""

import numpy as np                         # numpy library
from readlammpstrj import LAMMPStrj as trj # read lmp trj
import argparse                            # handle command line arguments
import sys

# -- Input argument parser

wa_parser = argparse.ArgumentParser(description='The box counting algorithim')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='Required! LAMMPS trajectory file containing atom information: id type x y z xs ys zs')
wa_parser.add_argument('-t','--types', metavar='[1,2,...]', type=str, nargs=1,
                       help='Required! List of atom types to use e.g. [1,2,4] in lammpstrj file. White spaces = error!')
wa_parser.add_argument('-n','--nbox', metavar='N', type=int, nargs=1,
                       help='Partition system into N boxes in all dimensions. Default: N=50')
wa_parser.add_argument('-o', '--outputprefix', metavar=('filename'), type=str, nargs=1,
                       help='Output file name(s) prefix')     
wa_parser.add_argument('-s', '--skipframes', metavar=('Nf'), type=int, nargs=1,
                       help='Skip the first Nf frames of the input trajectory. Default: Nf=0')
wa_parser.add_argument('-e', '--every', metavar=('Ne'), type=int, nargs=1,
                       help='Use every Ne frame in trajectory. Default: Ne=1')                       

args = wa_parser.parse_args()

# --------------------------------------------------------------------------- #
# -- Evaluate input arguments from argparser

if not (args.inputfile):
    wa_parser.error("No Input file provided! Provide input file through --inputfile")

if not (args.types):
    wa_parser.error("No atom types provided! Provide types through --types")
        
if not (args.nbox):
    nbox = 50
    nboxcubed = nbox**3
else:
    nbox = args.nbox[0]
    nboxcubed = nbox**3

if not (args.every):
    every = 1
else:
    every = args.every[0]
    
if (args.types):
    types = args.types[0].replace(']', '')
    types = types.replace('[', '')
    types = types.split(',')
    types = [int(atype) for atype in types]

if not (args.outputprefix):
    prefix = 'output'
else:
    prefix = args.outputprefix[0]
    
# --------------------------------------------------------------------------- #


def get_system_boundaries(lmp_obj=None):
    
    if (lmp_obj == None):
        print "Warning! function takes readlammpstrj object as input"
        return 0
    else:
        system_boundaries = np.zeros((3,2))   # xlo,xhi;ylo,yhi;zlo,zhi  # initial system boundaries 
        system_size = np.zeros(3)

        xl = obj.system_size[0,0]
        xh = obj.system_size[0,1]
        yl = obj.system_size[1,0]
        yh = obj.system_size[1,1]
        zl = obj.system_size[2,0]
        zh = obj.system_size[2,1]

        system_boundaries[0,0] = xl
        system_boundaries[0,1] = xh
        system_boundaries[1,0] = yl
        system_boundaries[1,1] = yh
        system_boundaries[2,0] = zl
        system_boundaries[2,1] = zh
        
        system_size[0] = xh-xl
        system_size[1] = yh-yl
        system_size[2] = zh-zl
    
        return system_boundaries, system_size
                
  
def write_to_file(output, data, datafields):
    """ Write output file 
        len(headers) == len(datafields)
        headers = names of data fields
    """
    if (len(data) != len(datafields)):
        print "Warning! number of data fields != number of headers!"
        print 'len:   ', len(data), len(datafields)
        print 'shape: ', np.shape(data), np.shape(datafields)

    ofile = open(output,'w')
    header = "# chunk "
    for element in datafields:
        header += element + " "

    header = header + '\n'
    ofile.write(header)
    
    it = 0
    for i in xrange(len(data[0])):
        line = str(it) + " "
        it += 1
        for j in xrange(len(data)):
            line += str(float(data[j][i])) + " "
        line += "\n"
        ofile.write(line)
            
    ofile.close()
    print "Finished writing file: ", output
      
if args.inputfile[0]:
    """ Use the box counting algorithm to compute N and r.
        Df*log(N) ~ log(1/r)
    """
    
    obj = trj(args.inputfile[0]) # create LAMMPStrj object 
    nframes = obj.nframes        # number of time frames in trajectory    
    ## ----------------------------------------------------- ##
    ID='id'                      # atom id
    TYPE='type'                  # atom type
    X='x'                        # unscaled atom position x
    Y='y'                        # unscaled atom position y
    Z='z'                        # unscaled atom position z
    Xs='xs'                      # scaled atom position x 
    Ys='ys'                      # scaled atom position y
    Zs='zs'                      # scaled atom position z
    ## ----------------------------------------------------- ##
    # check if we are using fractional coordinates or unscaled:
    useCoords = None
    if obj.data.has_key(Xs):
        if obj.data.has_key(Ys):
            if obj.data.has_key(Zs):
                useCoords = 'fractional'
            
    if obj.data.has_key(X):
        if obj.data.has_key(Y):
            if obj.data.has_key(Z):
                useCoords = 'unscaled'
                Xs = X
                Ys = Y
                Zs = Z

    if useCoords is None:
        print "!!!Error!!! Coordinates missing!"
        print "Is the system not three dimensional?\nCoordinates in file:"
        for i in [Xs,Ys,Zs,X,Y,Z]:
            print i, obj.data.has_key(i)
        print "\nExiting..."
        sys.exit()
    ## ----------------------------------------------------- ##
    if args.skipframes is None:
        skipframes = 0
    else:
        skipframes = args.skipframes[0]

    if (nframes <= skipframes):
        wa_parser.error("--skipframes N < nframes")    
    ## ----------------------------------------------------- ##

    # -- Perform the box counting algorithm        
    print "\n   Computing...\n"        
    for i in range(nframes):
        # -- Skip frames?
        if (i >= skipframes and i%every == 0):
            data = obj.get_data()
            natoms = obj.natoms[-1]
            system_boundaries, system_size = get_system_boundaries(obj)
            
            #boxes = np.zeros((nbox-1, nboxcubed)) # memory inefficient!
            N = np.zeros((nbox-1,)) # more memory efficient
            R = np.zeros((nbox-1,3))

            # parallellizable part?        
            for j in xrange(1,nbox):
                    
                    n = np.zeros((j*j*j,))
    
                    if useCoords == 'fractional':        
                        # fractional system size
                        dr = 1.0/float(j)
                        r = np.array([dr,dr,dr])
                        R[j-1] = r
                        system_boundaries = np.zeros(np.shape(system_boundaries))
                        system_boundaries[:,1] = 1
                        system_size = np.ones(np.shape(system_size))
                    if useCoords == 'unscaled':
                        # unscaled coordinates
                        r = system_size/float(j)
                        R[j-1] = r
                    
                    for atom in range(natoms):
                        t = data[TYPE][atom]
                        if (t in types):
                    
                            xs = data[Xs][atom]
                            ys = data[Ys][atom]
                            zs = data[Zs][atom]
                            
                            # -- wrap particle into simulation box if outside box boundaries
                            if xs < system_boundaries[0,0]:
                                xs = xs + system_size[0]
                            if xs > system_boundaries[0,1]:
                                xs = xs - system_size[0]
                                
                            if ys < system_boundaries[1,0]:
                                ys = ys + system_size[1]
                            if ys > system_boundaries[1,1]:
                                ys = ys - system_size[1]
                                
                            if zs < system_boundaries[2,0]:
                                zs = zs + system_size[2]
                            if zs > system_boundaries[2,1]:
                                zs = zs - system_size[2]
                            
                            xb = int(xs/r[0])
                            yb = int(ys/r[1])
                            zb = int(zs/r[2])
                            
                            index = xb*j**2 + yb*j + zb
                            
                            n[index] = 1 # box is occupied

                    N[j-1] = np.sum(n,axis=0)

            rx = R[:,0]
            ry = R[:,1]
            rz = R[:,2]
            
            columns = ['Rx','Ry','Rz','N','logRx','logRy','logRz','logN']
            outdata = [rx, ry, rz, N, np.log(1.0/rx), np.log(1.0/ry), np.log(1.0/rz), np.log(N)]
            
            suffix = '_{:04d}.dat'.format(i)
            filename = prefix + suffix
            
            write_to_file(filename,outdata, columns)
    
    obj.close_trj
        

# ------------------------------------------------------------------- EOF