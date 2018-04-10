#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:29:25 2016

@author: goranbs
--------------------------------------

Based on:

eldip.py

-------------------------------------

Compute electric dipolemoment of molecules across a molecular system
- 1D binning
- 2D binning
compute orientational order parameter Tp = 1/2 <3cos(theta)**2 -1>, Tp E [-0.5,1.0]
Tp is computed the orientational order parameter relative to a vector n, 
provided through the commandline. 
Tp = -0.5 : electric dipole moment otrhogonal to n
Tp =  1.0 : electric dipole parallel with n

display help:
>> eldip.py -h

-------------------------------------

neweldip: Add functionality to dump all dipole moments in each bin
? How should this information be written to file?
- In chunks?
- One file for every chunk?
    - Makes it easy to display with gnuplot
    - Could be many files...

-------------------------------------

eldip.py reads a LAMMPS trajectory using the class: readlammpsdata.py
- charge (q), coordinates (x,y,z), atom id (id) and molecule id (mol) must 
  be given in the LAMMPS trajectory.

"""

import numpy as np                         # numpy library
from readlammpstrj import LAMMPStrj as trj # read lmp trj
import argparse                            # handle command line arguments

wa_parser = argparse.ArgumentParser(description='Compute average xyz behaviour from LAMMPS trajectory')
wa_parser.add_argument('-d','--dipolemoment', action='store_true', 
                    help='Method: Compute dipole moment of molecules (same molID) defined by group of types given in -t')
wa_parser.add_argument('-f', '--dump1d', action='store_true',
                       help='Dump one file for each bin containing dipole moments of all molecules. Cartesian and spherical coordinates.')
wa_parser.add_argument('-t','--types', metavar='[1,2,...]', type=str, nargs=1,
                       help='list of types in lammpstrj file. White spaces = error!')
wa_parser.add_argument('-n', '--normal', metavar=('x','y','z'), type=float, nargs=3,
                       help='tuple: x y z. Default: 0 0 1')
wa_parser.add_argument('-x', '--xrange', metavar=('xmin','xmax'), type=float, nargs=2,
                       help='box boundary x-direction, xmin > xlo, xmax < xhi')
wa_parser.add_argument('-y', '--yrange', metavar=('ymin','ymax'), type=float, nargs=2,
                       help='box boundary y-direction, ymin > ylo, ymax < yhi')
wa_parser.add_argument('-z', '--zrange', metavar=('zmin','zmax'), type=float, nargs=2,
                       help='box boundary z-direction, zmin > zlo, zmax < zhi')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='LAMMPS trajectory file containing atom information: id mol type q x y z')
wa_parser.add_argument('-o', '--outputprefix', metavar=('filename'), type=str, nargs=1,
                       help='output file name(s) prefix')     
wa_parser.add_argument('--bin1d', metavar=('dim', 'origin', 'delta'), type=str, nargs=3,
                       help='1D binning. Default: z lower 1. dim=x,y,z origin=lower,upper delta=thickness of spatial bins.')

args = wa_parser.parse_args()

types = args.types[0].replace(']', '')
types = types.replace('[', '')
types = types.split(',')
types = [int(atype) for atype in types]
prefix = args.outputprefix

if args.xrange is None:
    args.xrange = [None,None]
if args.yrange is None:
    args.yrange = [None,None]
if args.zrange is None:
    args.zrange = [None,None]
if args.normal is None:
    normal = (0.0, 0.0, 1.0)
if args.normal is not None:
    normal = (float(args.normal[0]), float(args.normal[1]), float(args.normal[2]))

print "compute dipole moment: ", args.dipolemoment
print "xrange               : ", args.xrange
print "yrange               : ", args.yrange
print "zrange               : ", args.zrange
print "normal               : ", normal
print "types                : ", types
print "input file name      : ", args.inputfile
print "output prefix        : ", prefix

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
        
def get_1d_nbins(system_boundaries=None, dim='x', origin='lower', delta=1.0):
    """ --bin1d dim origin delta
        dim = x or y or z
        origin = lower or upper
        delta = thickness of spatial bins
    """

    d = get_number_dim_from_dim(dim)
    
    l = system_boundaries[d,1] - system_boundaries[d,0]
    delta = float(delta)
    
    if (delta > l):
        print "Warning! Given delta > system size! Resetting delta = l = ", l
        delta = l
    
    nbins = int(l/float(delta) + 1) # number of bins
    bins = np.zeros((nbins,10))      # array of bins; [px,py,pz,P**2,theta,nm,phi]
    baskets = np.empty((nbins,),dtype=object)
    for i in xrange(nbins): baskets[i] = []  # fill with empty lists
    
    # center position of spatial bins from origin:
    lloc = system_boundaries[d,0] + delta
    hloc = system_boundaries[d,1] - delta
    if (origin == 'upper'):
        pos = np.linspace(hloc, lloc, nbins)
    else:
        pos = np.linspace(lloc, hloc, nbins) # Default: origin = 'lower'
        
    return pos, bins, nbins, baskets
    

def compute_dipolemoment(atoms,normal,system_boundaries):
    """ takes a list of tuples [(id1,q1,x1,y1,z1), (id2,q2,x2,y2,z2), ...]
        1) computes geometrical center of molecule
        2) computes electric dipole moment of molecule (collection of atoms
        3) computes cos(theta) = (p dot n )/ (p dot p)
        4) If molecule compose of 3 atoms; compute H-O-H angle (hardcoded for water)
    """
    dX = (system_boundaries[0][1] - system_boundaries[0][0])/2.
    dY = (system_boundaries[1][1] - system_boundaries[1][0])/2.
    dZ = (system_boundaries[2][1] - system_boundaries[2][0])/2.

    molecule = np.zeros(np.shape(atoms)) # whole molecule
    
    p = np.zeros((3,), dtype=float)   # electric dipole moment
    r = np.zeros((3,), dtype=float)   # initialize geometric center of molecule
    pp = np.zeros((9,), dtype=float)  # Dx Dy Dz Tp Tpx Tpy Tpz theta phi
    
    natoms = len(atoms)
    for i in xrange(natoms):
        molecule[i] = atoms[i] # initialize molecules atom coordinates
    
    for i in xrange(1,natoms,1):
        dx = molecule[i-1][2] - molecule[i][2]
        dy = molecule[i-1][3] - molecule[i][3]
        dz = molecule[i-1][4] - molecule[i][4]
        

        if (dx > dX):
            dx -= dX*2
            #print "dx > dX: ", molecule[i-1][0], molecule[i][0], molecule[i-1][2], molecule[i][2]
            molecule[i][2] = molecule[i-1][2] - dx
        if (dx <= -dX):
            dx += dX*2
            #print "dx < dX: ", molecule[i-1][0], molecule[i][0], molecule[i-1][2], molecule[i][2]
            molecule[i][2] = molecule[i-1][2] - dx
            
        if (dy > dY):
            dy -= dY*2
            #print "dy > dY: ", molecule[i-1][0], molecule[i][0], molecule[i-1][3], molecule[i][3]
            molecule[i][3] = molecule[i-1][3] - dy
        if (dy <= -dY):
            dy += dY*2
            #print "dy < dY: ", molecule[i-1][0], molecule[i][0], molecule[i-1][3], molecule[i][3]
            molecule[i][3] = molecule[i-1][3] - dy
            
        if (dz > dZ):
            dz -= dZ*2
            #print "dz > dZ: ", molecule[i-1][0], molecule[i][0], molecule[i-1][4], molecule[i][4]
            molecule[i][4] = molecule[i-1][4] - dz
        if (dz <= -dZ):
            dz += dZ*2
            #print "dz < dZ: ", molecule[i-1][0], molecule[i][0], molecule[i-1][4], molecule[i][4]
            molecule[i][4] = molecule[i-1][4] - dz
                 
    it = 0
    rp = np.zeros((3,natoms), dtype=float)
    for atom in molecule:
        x = atom[2]
        y = atom[3]
        z = atom[4]
        #print np.size(rp), np.shape(rp)
        rp[it][0] = x
        rp[it][1] = y
        rp[it][2] = z
        r[0] += x
        r[1] += y
        r[2] += z
        it += 1
                
    r = r/natoms                     # geometric center of molecule
    
    for atom in molecule:
        q = atom[1]                  # charge of atom
        p[0] += q*(atom[2] - r[0])   # p = q_i*(r_i - r_c)
        p[1] += q*(atom[3] - r[1])
        p[2] += q*(atom[4] - r[2])
        
    pdotn = np.dot(p,normal)
    lp = np.linalg.norm(p)
    ln = np.linalg.norm(normal)
    costheta = pdotn/(lp*ln)
    theta = np.arccos(costheta)*180/np.pi
    
    pdotx = np.dot(p,(1,0,0))       # x 
    pdoty = np.dot(p,(0,1,0))       # y
    pdotz = np.dot(p,(0,0,1))       # z
    costhetax = pdotx/lp            # cos(theta) relative to x-component
    costhetay = pdoty/lp            # cos(theta) relative to y-component
    costhetaz = pdotz/lp            # cos(theta) relative to z-component
    
    pp[0] = p[0]                    # el. dipole mom. x
    pp[1] = p[1]                    # el. dipole mom. y
    pp[2] = p[2]                    # el. dipole mom. z
    pp[3] = costheta*costheta       # cos^2(t) of ang. betw. el.dip.mom. and normal vec.
    pp[4] = costhetax*costhetax     # -//- relative to x vector
    pp[5] = costhetay*costhetay     # -//- relative to y vector
    pp[6] = costhetaz*costhetaz     # -//- relative to z vector
    pp[7] = theta                   # angle betw. el. dip.mom. and normal (n) vec.

    # Only for computing Phi
    if (natoms == 3):
        # Find oxygen atom from its charge:
        if (atoms[0][1]<0):
            iO = 0
            iH1 = 1
            iH2 = 2
        if (atoms[1][1]<0):
            iO = 1
            iH1 = 0
            iH2 = 2
        else:
            iO = 2
            iH1 = 1
            iH2 = 0
        
        rOH1 = rp[iH1] - rp[iO]
        rOH2 = rp[iH2] - rp[iO]
        
        OH1dotOH2 = np.dot(rOH1,rOH2)
        ROH1 = np.linalg.norm(rOH1)
        ROH2 = np.linalg.norm(rOH2)
        cosphi = OH1dotOH2/(ROH1*ROH2)
        phi = np.arccos(cosphi)*180/np.pi
    else:
        phi = 666
    
    pp[8] = phi  # H-O-H angle of water molecule
    
    return pp, r
    
    
def add_to_bin(p, r, bins, delta, boundaries, dim, baskets):
    """ Check if position is within given boundaries
        add dipolemoment to correct bin in bins
    """
    # 1) check if r is within boundaries
    # 2) compute what bin to put it in

    d = get_number_dim_from_dim(dim)
    
    xmin = boundaries[0,0]
    xmax = boundaries[0,1]
    ymin = boundaries[1,0]
    ymax = boundaries[1,1]
    zmin = boundaries[2,0]
    zmax = boundaries[2,1]
    
    if (xmin < r[0] and r[0] < xmax ):
        if (ymin < r[1] and r[1] < ymax ):
            if (zmin < r[2] and r[2] < zmax ):
                # within the limits, so we'll include the molecule
                rd = r[d] - boundaries[d,0]
                index = int(rd/delta)  # compute bin index of the molecule
                #print index, rd, (lmin-lmax), lmin, lmax, r[d]
                bins[index][0] += p[0]  # Dx
                bins[index][1] += p[1]  # Dy
                bins[index][2] += p[2]  # Dz
                bins[index][3] += p[3]  # costhetasquared (relatice to n)
                bins[index][4] += p[4]  # costhetasquared (relatice to x (100))
                bins[index][5] += p[5]  # costhetasquared (relatice to y (010))
                bins[index][6] += p[6]  # costhetasquared (relatice to z (001))
                bins[index][7] += p[7]  # theta
                bins[index][8] += p[8]  # phi
                bins[index][9] += 1     # number of molecules in bin
                baskets[index].append((p[0],p[1],p[2]))


def write_to_file(output, data, datafields):
    """ Write output file 
        len(headers) == len(datafields)
        headers = names of data fields
    """
    if (len(data) != len(datafields)):
        print "Error! number of data fields != number of headers!"
        print 'len:   ', len(data), len(datafields)
        print 'shape: ', np.shape(data), np.shape(datafields)

    ofile = open(output,'w')
    header = "chunk "
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
    

if args.dipolemoment:
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
    Xs='xs'                      # scaled atom position x # not supported
    Ys='ys'                      # scaled atom position y # not supported
    Zs='zs'                      # scaled atom position z # not supported
    ## ----------------------------------------------------- ##    
    
    nframes = obj.nframes        # number of time frames in trajectory    

    allpos = []                  # bin1d bins
    allbins = []
    allnbins = []
    allnatoms = []
    allbaskets = []
    
    aposX = []                   # bin2d bins
    aposY = []
    a2Dbins = [] 
    
    for i in range(nframes):
        print "## ------ TIMEFRAME ", i, "/", nframes, "------ ##"
        
        data = obj.get_data()
        natoms = obj.natoms[-1]
        system_boundaries, ud_boundaries = get_system_boundaries(obj)
        
        ## ----------------------------------------------------- ##
        molecules = dict()
        
        for atom in range(natoms):
            t = data[TYPE][atom]
            if (t in types):
                mol = data[MOL][atom]
                if mol in molecules:
                    id_ = data[ID][atom]
                    q = data[Q][atom]
                    x = data[X][atom]
                    y = data[Y][atom]
                    z = data[Z][atom]
                    t = data[TYPE][atom]
                    
                    molecules[mol].append((id_,q,x,y,z,t))
                    
                else:
                    id_ = data[ID][atom]
                    q = data[Q][atom]
                    x = data[X][atom]
                    y = data[Y][atom]
                    z = data[Z][atom]
                    t = data[TYPE][atom]
                    
                    molecules.update({mol:[(id_,q,x,y,z,t)]})
        
        ##-- 1D binning
        if (args.bin1d is not None):
            dim = args.bin1d[0]            # dimention to perform binning
            origin = args.bin1d[1]         # origin of bins (lower, upper)
            delta = float(args.bin1d[2])   # bin spacing
            pos, bins, nbins, baskets = get_1d_nbins(ud_boundaries,dim,origin,delta)
            
            for molecule in molecules:
                # compute electric dipole moment of molecules and add to 1D bin
                p,r = compute_dipolemoment(molecules[molecule], normal, system_boundaries)
                add_to_bin(p, r, bins, delta, ud_boundaries, dim, baskets)
            
            allpos.append(pos)             # append to bucket of bins:
            for j in xrange(nbins):
                divisor = bins[j][9]       # number of molecules in bin: N
                if (divisor == 0):         # if no particles in bin
                    bins[j][3] = 1/3.      # this assures Tp(z) = -0.5
                    bins[j][4] = 1/3.      # this assures Tp(z) = -0.5
                    bins[j][5] = 1/3.      # this assures Tp(z) = -0.5
                    bins[j][6] = 1/3.      # this assures Tp(z) = -0.5
                    bins[j][8] = 90        # phi if no samples
                else:
                    bins[j][3] /= divisor  # divide costheta**2 by number of molecules
                    bins[j][4] /= divisor  # average angle theta of molecules in bin
                    bins[j][5] /= divisor  # average angle theta of molecules in bin
                    bins[j][6] /= divisor  # average angle theta of molecules in bin
                    bins[j][7] /= divisor  # average theta
                    bins[j][8] /= divisor  # average phi
                
            allbins.append(bins)
            allnbins.append(nbins)
            allbaskets.append(baskets)
        
            
    if (args.bin1d is not None):        
        # averaging bin1d bins:
    
    
        nbins = allnbins[0]
        px = np.zeros( (nbins, 1) )     # electric dipolemoment x-dir
        py = np.zeros( (nbins, 1) )     # y-dir
        pz = np.zeros( (nbins, 1) )     # z-dir
        Tp = np.zeros( (nbins, 1) )     # orientational order parameter relative to input; dim
        Tpx = np.zeros( (nbins, 1) )    # orientational order parameter relative to (100)
        Tpy = np.zeros( (nbins, 1) )    # orientational order parameter relative to (010)
        Tpz = np.zeros( (nbins, 1) )    # orientational order parameter relative to (001)
        theta = np.zeros( (nbins, 1) )  # theta
        phi = np.zeros( (nbins, 1) )    # phi
        Nmolecs = np.zeros( (nbins,1) ) # Average number of molecules in bin
        d = get_number_dim_from_dim(dim)# dim
        
        baskets = np.empty((nbins,),dtype=object)
        for i in xrange(nbins): baskets[i] = []  # fill with empty lists
        
        for i in xrange(len(allbins)):
            abin = allbins[i]
            baskets += allbaskets[i]
            for j in xrange(len(abin)):
                px[j] += abin[j][0]
                py[j] += abin[j][1]
                pz[j] += abin[j][2]
                Tp[j] += (3*abin[j][3] - 1) # ref DOI: 10.1103/PhysRevLett.101.056102
                Tpx[j] += (3*abin[j][4] - 1) # ref DOI: 10.1103/PhysRevLett.101.056102
                Tpy[j] += (3*abin[j][5] - 1) # ref DOI: 10.1103/PhysRevLett.101.056102
                Tpz[j] += (3*abin[j][6] - 1) # ref DOI: 10.1103/PhysRevLett.101.056102
                theta[j] += abin[j][7]
                phi[j] += abin[j][8]
                Nmolecs[j] += abin[j][9]
           
        
        px = px/nframes
        py = py/nframes
        pz = pz/nframes
        Tp = Tp/(2*nframes)
        Tpx = Tpx/(2*nframes)
        Tpy = Tpy/(2*nframes)
        Tpz = Tpz/(2*nframes)
        theta = theta/nframes
        phi = phi/nframes
        nmol = Nmolecs/nframes
        
        ## write to file:
        coord = 'coord_' + dim
        Tpdim = 'Tp_' + str(int(normal[0])) + str(int(normal[1])) + str(int(normal[2]))
        thetadim = 'theta_' + str(int(normal[0])) + str(int(normal[1])) + str(int(normal[2]))
        
        datafields = [coord, 'px', 'py', 'pz', 'Tpx', 'Tpy', 'Tpz', Tpdim, thetadim, 'phi', 'Nmolecules']
        data = [allpos[0],px,py,pz,Tpx,Tpy,Tpz,Tp,theta,phi,nmol]
        
        if prefix == None:
            outfile = "ed_bin1d.dat"
        else:
            outfile = prefix[0] + "_ed_bin1d.dat"
            
        write_to_file(outfile, data, datafields)
        
        if (args.dump1d):
            # dump electric dipole moments. One file for each bin...            
            IT=0
            for basket in baskets:
                IT+=1
                outf = 'dump1d_{:04d}.dat'.format(IT)
                data = ['px','py','pz']
                
                ofile = open(outf,'w')
                header = '#id px py pz phi theta r\n'
                ofile.write(header)    
                #print np.shape(basket), len(basket)
                it=0
                for px, py, pz in basket:
                    it+=1
                    r = np.sqrt(px**2 + py**2 + px**2)      # radius
                    if (px == 0):
                        th = 1.5707963267948966             # lim(atan(x)) x -> inf
                    else:
                        th = np.arctan(py/px)               # polar (theta)
                    
                    try:
                        pzor = pz/r
                    except:
                        print 'r = {}'.format(r)
                        pzor = pz
                        
                    if (pzor < -1):
                        az = np.pi                          # acos(-1) = pi
                    if (pzor > 1):
                        az = np.pi                          # acos(1) = pi
                    else:
                        az = np.arccos(pz/r)                # azimuthal (phi)
                    
                    ofile.write("{} {} {} {} {} {} {}\n".format(it, px, py, pz, az, th, r))

            
        #-- End of bin1d compute ----------------------------------------------        
        
    obj.close_trj()
    ##-- End of method

                
        

                
            
            
        
        
        
        
        
        


        
        