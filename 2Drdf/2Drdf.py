#!/usr/bin/python
"""
Created on Mon Jul  9 16:46:39 2018

@author: goranbs

compute 2D radial distribution functions from LAMMPS trajectory
compute the rdf in the xy-plane

g_{A}{B}(R)

A,B: give list A of atom types, and list B of atom types.
dz : give dz to determine the thickness of the plane.
p  : give (x0,y0,z0) for a point on the plane.

TODO How to account for periodic boundaries in the g(r)?
    - We could replicate the system internally. One to the left and one to the right in x,y,z (26 replicas).

TODO How can we make the script more general? Define a plane by a surface normal,
    and compute the rdf in that plane.
    - The challenge is to compute the reference volume of the slice.
    this slice is the intersection between the surface and the simulation box.
    - however, this is just a problem of correct normalization.
  
n = (a,b,c)                       # plane normal vector
p = (x0,y0,z0)                    # point on the plane
r = (x,y,z)                       # test vector (coordinate of an atom)

n * (r - p) = 0                   # a point r is in the plane if this holds
P(r) = aX + bY + cZ + d = 0       # plane equation
D = abs(P(r))/sqrt(a^2+b^2+c^2)   # distance D of point r from plane
n cross (r cross n)               # components of r projected onto plane

"""

import numpy as np                          # numpy library
import argparse                             # handle command line arguments
import sys, os                              # sys and os
PATH = os.path.dirname(os.path.realpath(__file__))
trjPATH = PATH + '/../readlammps' 
sys.path.insert(0, trjPATH)        # append relative path to readlammpstrj
from readlammpstrj import LAMMPStrj as trj  # import module for reading LAMMPS trajectories
import matplotlib.pyplot as plt             # plotting


wa_parser = argparse.ArgumentParser(description='Compute the 2D radial distribution function in the xy-plane')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='LAMMPS trajectory file containing atom information: id mol type q x y z')
wa_parser.add_argument('-p', '--point', metavar='[x0,y0,z0]', type=str, nargs=1,
                       help='Vector point on surface with surface normal n=(0,0,1).')     
wa_parser.add_argument('-a','--A', metavar='[1,2,...]', type=str, nargs=1,
                       help='List A of atom types e.g. [1,2,4] in lammpstrj file. To compute g_{A,B}(r). Default: A = [1]')
wa_parser.add_argument('-b','--B', metavar='[1,2,...]', type=str, nargs=1,
                       help='List B of atom types e.g. [1,2,4] in lammpstrj file. To compute g_{A,B}(r). Default: B = [1]')                       
#wa_parser.add_argument('-n', '--normal', metavar='[x,y,z]', type=str, nargs=1,
 #                      help='Surface normal [x,y,z]. Default: (0,0,1)')
wa_parser.add_argument('-d', '--dz', metavar='dz', type=float, nargs=1,
                       help='Sampling region away from xy surface plane. Thickness of sampling region is 2*dz. Default: dz = 1')
wa_parser.add_argument('-o', '--outputprefix', metavar=('filename'), type=str, nargs=1,
                       help='Output file name prefix. Default: input filename')
wa_parser.add_argument('-e', '--every', metavar=('Ne'), type=int, nargs=1,
                       help='Use every Ne frame. Default: Ne = 1')                       
wa_parser.add_argument('-s', '--skipframes', metavar=('Ns'), type=int, nargs=1,
                       help='Skip the first Ns frames. Default: Ns = 0')
wa_parser.add_argument('--nbins', metavar=('Nb'), type=int, nargs=1,
                       help='Define the number of bins used for the histogram. Default: Nb = 1000')
wa_parser.add_argument('--rmax', metavar=('R'), type=int, nargs=1,
                       help='Define the maximum radius used for the histogram. Default: R = max(dr)')                       
wa_parser.add_argument('--plot', action='store_true', 
                    help='Show plot of the result.')      

args = wa_parser.parse_args()

# --------------------------------------------------------------------------- #

# ---------------------------------- CHECK for input file
if not (args.inputfile):
    wa_parser.error("No Input file provided! Provide input file through --inputfile")

# ---------------------------------- CHECK for reference point
if not (args.point):
    wa_parser.error("No reference point for the surface plane provided! Example: --point [3,2,1]")
else:
    P = args.point[0].replace(']','')
    P = P.replace('[','')
    P = P.split(',')
    Pvec = np.array([float(p) for p in P])
    if (len(P)!=3):
        wa_parser.error("Error! Reference point p = {}".format(args.point[0]))

# ---------------------------------- CHECK for atom types in lists A and B
if (args.A):
    A = args.A[0].replace(']', '')
    A = A.replace('[', '')
    A = A.split(',')
    Atypes = [int(atype) for atype in A]
else:
    Atypes = 1    

if (args.B):
    B = args.B[0].replace(']', '')
    B = B.replace('[', '')
    B = B.split(',')
    Btypes = [int(atype) for atype in B]
else:
    Btypes = 1
    
# ---------------------------------- CHECK for surface normal
# A surface normal could be used to make the code more general,
    # however, it turned out to be difficult to compute the volume of the 
    # slice intersecting the simulation box.
#if (args.normal):
#    N = args.normal[0].replace(']','')
#    N = N.replace('[','')
#    N = N.split(',')
#    Nvec = np.array([float(n) for n in N])
#    if (len(N)!=3):
#        wa_parser.error("Error! Surface normal N = {}".format(args.normal[0]))
#else:
#    Nvec = np.array([0,0,1])
    
# ---------------------------------- CHECK for thickness of plane
if (args.dz):
    dz = float(args.dz[0])
    if (dz <= 0):
            wa_parser.error("Error! dz must be greater than 0")
else:
    dz = 1.0
    
# ---------------------------------- CHECK for output prefix
if (args.outputprefix):
    prefix = args.outputprefix[0]
else:
    prefix = args.inputfile[0]
    prefix = prefix.split('.')
    prefix = prefix[0]

# ---------------------------------- CHECK for Ne and skipframes
if (args.every):
    Ne = args.every[0]
else:
    Ne = 1
if (args.skipframes):
    Ns = args.skipframes[0]
else:
    Ns = 0

# ---------------------------------- CHECK for radius range (R) and number of bins (Nb)
if (args.nbins):
    nbins = int(args.nbins[0])
    if (nbins <= 0):
        wa_parser.error("Error! nbins = {}".format(args.nbins[0]))
else:
    nbins = 1000
if (args.rmax):
    Rmax = float(args.rmax[0])
    if (Rmax <= 0):
        wa_parser.error("Error! rmax = {}".format(args.rmax[0]))
else:
    Rmax = None
# --------------------------------------------------------------------------- #

print "\n# -- -- -- -- -- -- -- -- Input settings -- -- -- -- -- -- -- -- #"
print "reference point P                 :", Pvec
print "list of atom types A              :", Atypes
print "list of atom types B              :", Btypes
#print "surface normal N                  :", Nvec
print "thickness of sampling region 2*dz :", 2*dz
print "use every Ne data frame           :", Ne
print "skip first Ns data frames         :", Ns
print "output data file name             : {}.dat".format(prefix)
print "# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #"

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

def printframe(frame,endframe):
    """
    Print time frame to terminal.
    """
    line = "\r    timeframe: {:d} / {:d}".format(frame, endframe)
    #print(line),
    sys.stdout.write(line)
    sys.stdout.flush()

def distance_from_plane(n,p,r,nnorm=None):
    """
    Compute the distance D = abs(P(r))/sqrt(a^2+b^2+c^2)
    n = surface normal
    p = point on surface
    r = point to be measured its distance from the surface.
    """
    #return np.abs(np.dot(n,(p-r)))/np.linalg.norm(n)
    #return np.abs(np.dot(n,(p-r)))/nnorm
    # the normal vector is already a unit vector!
    return np.abs(np.dot(n,(p-r)))

def distance_from_xy_plane(p,r):
    """
    Compute the distance in the z-component.
    
    """
    return np.abs(p[2]-r[2])
    
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
    ofile.write("# g(r) in the xy-plane from 2Drdf.py\n")
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


# --------- Compute 2D radial distribution function

obj = trj(args.inputfile[0]) # create LAMMPStrj object 
nframes = obj.nframes        # number of time frames in trajectory    
## ----------------------------------------------------- ##
ID='id'                      # atom id
MOL='mol'                    # molecule id
TYPE='type'                  # atom type
Q='q'                        # charge
X='x'                        # unscaled atom position x
Y='y'                        # unscaled atom position y
Z='z'                        # unscaled atom position z
## ----------------------------------------------------- ##    


D = []
D_plane = []
#nnorm = np.linalg.norm(Nvec)
#Nvec = Nvec/nnorm            # make sure Nvec is a unit vector
Nvec = np.array([0,0,1])
Na=0
Nb=0
V_slice=0
count=0

print "\n   Computing...\n"
for i in range(nframes):
    data = obj.get_data()
    if (i >= Ns and i%Ne == 0):
        printframe(i+1,nframes)
        natoms = obj.natoms[-1]
        system_boundaries, system_size = get_system_boundaries(obj)
        V_slice += system_size[0]*system_size[1]*2*dz
        count += 1
        for i in range(natoms-1):
            IDi = data[ID][i]
            ti = data[TYPE][i]
            if (ti in Atypes):
                ri = np.array([data[X][i],data[Y][i],data[Z][i]])
                if (distance_from_xy_plane(Pvec,ri) < dz):
                    Na += 1
                    for j in range(i+1,natoms):
                        IDj = data[ID][j]
                        tj = data[TYPE][j]
                        if (tj in Btypes):
                            if (IDi != IDj):
                                rj = np.array([data[X][j],data[Y][j],data[Z][j]])
                                if (distance_from_xy_plane(Pvec,rj) < dz):                                
                                    Nb += 1
                                    dij = (ri - rj)                                    
                                    dij_n = np.cross(Nvec,np.cross(dij,Nvec)) # components of vector projected onto plane
                                    D.append(np.linalg.norm(dij))
                                    D_plane.append(np.linalg.norm(dij_n))

obj.close_trj()
# --------------------------------------------------------------------------- #

print "\n"
V_slice = V_slice/count
Na = float(Na)/count
Nb = 2*float(Nb)/Na/count
D = np.array(D)
D_plane = np.array(D_plane)
print "Na: ", Na, "Nb: ", Nb, "V_slice: ", V_slice

if (args.rmax):
    hD, binsD = np.histogram(D,nbins,range=(0.0,Rmax),density=True)
    hD_plane, binsD_plane = np.histogram(D_plane,nbins,range=(0.0,Rmax),density=True)
else:
    hD, binsD = np.histogram(D,nbins,density=True)
    hD_plane, binsD_plane = np.histogram(D,nbins,density=True)
    

#weights = (Nb/V_slice)*np.pi*2*dz*np.array(binsD[1:]**2 - binsD[:-1]**2)
weights = np.array(binsD[1:]**2 - binsD[:-1]**2)

normed_histD = hD/weights
normed_histD_plane = hD_plane/weights

filename = prefix + ".dat"
data = [binsD[:-1], normed_histD, normed_histD_plane]
datafields = ["r", "gAB(r)", "gAB_inplane(r)"]
write_to_file(filename, data, datafields)

if (args.plot):
    plt.figure()
    plt.plot(binsD[:-1],normed_histD, label='spatial distances')
    plt.hold(True)
    plt.plot(binsD[:-1], normed_histD_plane, label='projected onto plane')
    plt.xlabel('r')
    plt.ylabel('$g^{A,B}(r)$')
    plt.title("2D rdf")
    plt.legend()
    plt.hold(False)
    plt.show()

# --------------------------------------------------------------------------- #