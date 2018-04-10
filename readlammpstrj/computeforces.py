#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:29:25 2016

@author: goranbs
--------------------------------------

Based on:

eldip.py

-------------------------------------

Compute forces on group of atoms.
The group of atoms needs to be well defined. HOW?
- atomID
- location

1) Define group of atoms in same manner as LAMMPS:
    group A box xmin xmax ymin ymax zmin zmax
    group B type a b c ..
    group C intersect A B

display help:
>> computeforces.py -h

This will be a brure force script for computing the forces on specific groups
of atoms. The CaCO3 slabs.

-------------------------------------

TODO:

-------------------------------------

"""

import numpy as np                         # numpy library
from readlammpstrj import LAMMPStrj as trj # read lmp trj
import argparse                            # handle command line arguments

wa_parser = argparse.ArgumentParser(description='Compute average forces from LAMMPS trajectory')
wa_parser.add_argument('-f','--forces', action='store_true', 
                    help='Method: Compute total force on group defined by types given in -t and within x,y,z bounds.')
wa_parser.add_argument('-t','--types', metavar='[1,2,...]', type=str, nargs=1,
                       help='list of types in lammpstrj file. White spaces = error!')
wa_parser.add_argument('-x', '--xrange', metavar=('xmin','xmax'), type=float, nargs=2,
                       help='box boundary x-direction, xmin > xlo, xmax < xhi')
wa_parser.add_argument('-y', '--yrange', metavar=('ymin','ymax'), type=float, nargs=2,
                       help='box boundary y-direction, ymin > ylo, ymax < yhi')
wa_parser.add_argument('-z', '--zrange', metavar=('zmin','zmax'), type=float, nargs=2,
                       help='box boundary z-direction, zmin > zlo, zmax < zhi')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='LAMMPS trajectory file containing atom information: id mol type q x y z')

args = wa_parser.parse_args()

types = args.types[0].replace(']', '')
types = types.replace('[', '')
types = types.split(',')
types = [int(atype) for atype in types]

if args.xrange is None:
    args.xrange = [None,None]
if args.yrange is None:
    args.yrange = [None,None]
if args.zrange is None:
    args.zrange = [None,None]

print "compute forces       : ", args.forces
print "xrange               : ", args.xrange
print "yrange               : ", args.yrange
print "zrange               : ", args.zrange
print "types                : ", types
print "input file name      : ", args.inputfile

# -------------------------------------------------------

#########################################################
#
# Computes
#
#########################################################
    

def get_system_boundaries(lmp_obj=None):
    
    if (lmp_obj == None):
        print "Error! function takes readlammpstrj object as input"
        return 0
    else:
        system_boundaries = np.zeros((3,2))   # xlo,xhi;ylo,yhi;zlo,zhi  # initial system boundaries 
        user_def_boundaries = np.zeros((3,2)) # xlo,xhi;ylo,yhi;zlo,zhi  # user defined system boundaries

        xl = obj.isystem_size[0,0]          # initial system sizes
        xh = obj.isystem_size[0,1]
        yl = obj.isystem_size[1,0]
        yh = obj.isystem_size[1,1]
        zl = obj.isystem_size[2,0]
        zh = obj.isystem_size[2,1]
        system_boundaries[0,0] = xl         # initial system boundaries
        system_boundaries[0,1] = xh
        system_boundaries[1,0] = yl
        system_boundaries[1,1] = yh
        system_boundaries[2,0] = zl
        system_boundaries[2,1] = zh

        xmin = args.xrange[0]               # given ranges
        ymin = args.yrange[0]
        zmin = args.zrange[0]
        xmax = args.xrange[1]
        ymax = args.yrange[1]
        zmax = args.zrange[1]
        
        if xmin is None:
            xmin = xl
            xmax = xh
        else:
            if xmin < xl:
                xmin = xl
            if xmax > xh:
                xmax = xh
            
        user_def_boundaries[0,0] = xmin
        user_def_boundaries[0,1] = xmax
        
        if ymin is None:
            ymin = yl
            ymax = yh
        else:
            if ymin < yl:
                ymin = yl
            if ymax > yh:
                ymax = yh
            
        user_def_boundaries[1,0] = ymin
        user_def_boundaries[1,1] = ymax
        
        if zmin is None:
            zmin = zl
            zmax = zh
        else:
            if zmin < zl:
                zmin = zl
            if zmax > zh:
                zmax = zh
            
        user_def_boundaries[2,0] = zmin
        user_def_boundaries[2,1] = zmax
        
        #print user_def_boundaries
        #print system_boundaries
        
        return system_boundaries, user_def_boundaries
        
def get_number_dim_from_dim(dim='x'):
    """ return 0 for x, 1 for y and 2 for z """    
    d=0
    if (dim == 'x'):
        d=0
    if (dim == 'y'):
        d=1
    if (dim == 'z'):
        d=2
    if (dim not in ['x','y','z']):
        print "Error! given dim is not valid! dim = ", dim
        print "Using default: dim = x = 0"
        d = 0
    
    return d
    

if args.forces:
    """ Compute dipole moment D, of molecules in system.
        Compute orientational order parameter; Tp = 1/2 < 3cos^2(theta) - 1>
        Compute angle theta; |D||n|cos(theta) = D dot n
        
        Restrictions:
        1) Molecule ID (mol) and atom ID's (id) must be present in lammpstrj file
        2) Positions of atoms must be present in lammpstrj file (x,y,z)
        3) Atom types (type) must be present in lammpstrj file
        4) The given atom types must belong to a molecule ID (mol)
        5) Geometrical center of molecules is used in stead of center of mass.
            Method is therefore approximate for molecules with charge =! 0.
    """
        
    obj = trj(args.inputfile[0]) # create LAMMPStrj object 
    ## ----------------------------------------------------- ##
    ID='id'                      # atom id
    MOL='mol'                    # molecule id
    TYPE='type'                  # atom type
    Q='q'                        # charge
    X='x'                        # unscaled atom position x
    Y='y'                        # unscaled atom position y
    Z='z'                        # unscaled atom position z
    FX='fx'                      # force component x
    FY='fy'                      # force component y
    FZ='fz'                      # force component z
    ## ----------------------------------------------------- ##    
    
    nframes = obj.nframes        # number of time frames in trajectory    

    F_tot = [0,0,0]              # total force on group
    
    for i in range(nframes):
        #print "## ------ TIMEFRAME ", (i+1), "/", nframes, "------ ##"
        
        data = obj.get_data()
        natoms = obj.natoms[-1]
        system_boundaries, ud_boundaries = get_system_boundaries(obj)
        ## ----------------------------------------------------- ##    
        j = 0
        for atom in range(natoms):
            t = data[TYPE][atom]
            if (t in types):
                x = data[X][atom]
                y = data[Y][atom]
                z = data[Z][atom]
                if ( ud_boundaries[0,0] <= x <= ud_boundaries[0,1] and ud_boundaries[1,0] <= y <= ud_boundaries[1,1] and ud_boundaries[2,0] <= z <= ud_boundaries[2,1]):                    

                    j += 1

                    fx = data[FX][atom]
                    fy = data[FY][atom]
                    fz = data[FZ][atom]
                    
                    F_tot[0] += fx
                    F_tot[1] += fy
                    F_tot[2] += fz
                    #print data[FZ][atom]
        #print j
        
    F_tot = np.array(F_tot)/nframes
    print " # ---------------------------------------------------------- #"
    print "Force ", F_tot, " kcal/mol/angstrom"
    nN = 0.06952
    MPa = 6952
    F_tot *= nN
    print "Force ", F_tot, " nN"
    print " # ---------------------------------------------------------- #"
        
    obj.close_trj()
    ##-- End of method

                
        

                
            
            
        
        
        
        
        
        


        
        