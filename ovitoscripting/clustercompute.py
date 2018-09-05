
"""
@author: goranbs

Ovito python script for analysing a high-Mg-constentration cluster in a
calcite crystal. This script utilizes the DXA algorithm provided in ovito
to characterize particles still belonging to the fcc-lattice structure of
the Mg-Ca ions, and the 'amorphous' layer enclosing the high-Mg-density cluster.
The radius of the inner sphere (still fcc-lattice) is computed
The radius of the outer sphere (non fcc-lattice) is computed

useage:
    >> ./ovitos clustercompute.py ifile.lammpstrj ofile.dat 19.0 2.2 [19,31.25,32.24] [3,4] 
    
input:
    1) lammps trajectory: ifile.lammpstrj
    2) output data: ofile.dat
    3) initial radius of the high-Mg-density cluster: sphere_radius
    4) approximate width of the amorphous layer around the Mg cluster: radius_threshold
    5) center of the Mg cluster: sphere_center=[x,y,z]
    6) particle types to exclude in the calculation: excludetypes=[a,b,c,...]

output:
    ofile.dat
        output data containing columns:
        1)  chunk       : frame number
        2)  vol_inner   : volume of inner sphere
        3)  vol_shell   : volume of shell around sphere
        4)  area_inner  : surface area of inner sphere
        5)  area_outer  : surface area of outer sphere
        6)  r1_inner    : radius of inner sphere (from volume)
        7)  r2_inner    : radius of inner sphere (from surface area)
        8)  r1_outer    : -//- outer sphere
        9)  r2_inner    : -//- outer sphere
        10) nrho_inner  : number density of inner sphere
        11) nrho_outer  : number density of outer shell
        12) nrho_rest   : number density of system outside outer sphere
        13) nrho_system : number density of system
        

"""

import numpy as np
import sys
import ovito

# --------------------------------------------------------------------- #

# call:
# >> ./ovitos hello.py ifile.lammpstrj ofile.dat 19.0 2.2 [19,31.25,32.24] [3,4] 
# --------------------------------------------------------------------- #

ifile = str(sys.argv[1])
ofile = str(sys.argv[2])
sphere_radius = float(sys.argv[3])
radius_threshold = float(sys.argv[4])
sphere_center = np.array([float(i) for i in sys.argv[5].strip('[]').split(',')])
excludetypes = set([int(i) for i in sys.argv[6].strip('[]').split(',')])

print("#"+"-"*40+"#")
print("input file name    : {}".format(ifile))
print("output file name   : {}".format(ofile))
print("sphere radius      : {:.2f} angstrom".format(sphere_radius))
print("radius threshold   : {:.2f} angstrom".format(radius_threshold))
print("center of sphere   : [{:.2f}, {:.2f}, {:.2f}]".format(sphere_center[0],sphere_center[1],sphere_center[2]))
print("exclude atom types : {}".format(sys.argv[6]))
print("#"+"-"*40+"#")

# --------------------------------------------------------------------- #

#ifile = "test2.lammpstrj"
#ofile = "out.dat"

#sphere_radius = 18.0                            # radius of sphere
#radius_threshold = 2.0                          # approx width of amorphous layer around sphere
#sphere_center = np.array([19, 31.25, 32.24])    # center of sphere
#particletypes = set([3,4])                      # set of particle types to not consider.

# --------------------------------------------------------------------- #

def volume(a,b,c,cosa,cosb,cosg):
    """ Compute volume of trigonal simulation box"""
    V = a*b*c*((1 - cosa*cosa - cosb*cosb - cosg*cosg) + 2*(cosa*cosb*cosg))**(1/2.)
    return V

def sphere_vol2r(volume):
    return np.power( 3*volume/(4*np.pi) , 1/3. )
    
def sphere_area2r(area):
    return np.power( area/(4*np.pi) , 1/2. )
    
def write_to_file(output, data, datafields):
    """ Write output file 
        len(headers) == len(datafields)
        headers = names of data fields
    """
    if (len(data) != len(datafields)):
        print("Error! number of data fields != number of headers!")
        print("len:   {} {}".format(len(data), len(datafields)))
        print("shape: {} {}".format(np.shape(data), np.shape(datafields)))

    ofile = open(output,'w')
    ofile.write("# output from ovitos.py\n# Description: inner=inner sphere, outer=shell around inner sphere, rest=everything outside outer sphere, system=global\n")
    header = "# chunk "
    for element in datafields:
        header += element + " "

    header = header + '\n'
    ofile.write(header)
    
    it = 0
    for i in range(len(data[0])):
        line = str(it) + " "
        it += 1
        for j in range(len(data)):
            line += str(float(data[j][i])) + " "
        line += "\n"
        ofile.write(line)
            
    ofile.close()
    print("Finished writing file: {}".format(output))
    
# --------------------------------------------------------------------- #

# -- perform calculations
    
radius_center_amorph = sphere_radius + radius_threshold
expression_outer_shell = "StructureType==1 && (Position.X - {})^2 + (Position.Y - {})^2 + (Position.Z - {})^2 > {}^2".format(sphere_center[0],sphere_center[1],sphere_center[2],radius_center_amorph)
expression_inner_shell = "StructureType==1 && (Position.X - {})^2 + (Position.Y - {})^2 + (Position.Z - {})^2 < {}^2".format(sphere_center[0],sphere_center[1],sphere_center[2],radius_center_amorph)


# -- node for outer shell: node1
node1 = ovito.io.import_file(ifile, multiple_frames = True)
data = node1.source
nparticles = data.number_of_particles
data = node1.compute()
nexcluded = 0
ntypes = []
nexcluded_string = " of which: "
for atype in excludetypes:
    ni = np.count_nonzero( data.particle_properties.particle_type.array == atype )
    nexcluded += ni
    ntypes.append(ni)
    nexcluded_string += "{} of type {}, ".format(ni, atype)
    
    
nincluded = nparticles - nexcluded
print("nparticles   : {}".format(nparticles))
print("nincluded    : {}".format(nincluded))
print("nexcluded    : {} {}".format(nexcluded,nexcluded_string))
print("#"+"-"*40+"#")

node1.modifiers.append(ovito.modifiers.SelectParticleTypeModifier(property = "Particle Type", types = excludetypes))
node1.compute()
node1.modifiers.append(ovito.modifiers.DeleteSelectedParticlesModifier())
node1.compute()
node1.modifiers.append(ovito.modifiers.DislocationAnalysisModifier())
node1.compute()
node1.modifiers.append(ovito.modifiers.SelectExpressionModifier(expression=expression_outer_shell))
node1.compute()
node1.modifiers.append(ovito.modifiers.DeleteSelectedParticlesModifier())
node1.compute()
node1.modifiers.append(ovito.modifiers.ConstructSurfaceModifier(radius=7, smoothing_level=8))
node1.compute()

# -- node for inner shell: node2
node2 = ovito.io.import_file(ifile, multiple_frames = True)
node2.modifiers.append(ovito.modifiers.SelectParticleTypeModifier(property = "Particle Type", types = excludetypes))
node2.compute()
node2.modifiers.append(ovito.modifiers.DeleteSelectedParticlesModifier())
node2.compute()
node2.modifiers.append(ovito.modifiers.DislocationAnalysisModifier())
node2.compute()
node2.modifiers.append(ovito.modifiers.SelectExpressionModifier(expression=expression_inner_shell))
node2.compute()
node2.modifiers.append(ovito.modifiers.InvertSelectionModifier())
node2.compute()
node2.modifiers.append(ovito.modifiers.DeleteSelectedParticlesModifier())
node2.compute()
node2.modifiers.append(ovito.modifiers.ConstructSurfaceModifier(radius=7, smoothing_level=8))
node2.compute()

# -- inner=inner sphere, outer=shell around inner sphere
datafields = ["vol_inner", "vol_outer","area_inner","area_outer","r1_inner","r2_inner","r1_outer","r2_outer","nrho_inner","nrho_outer","nrho_rest","nrho_system"]
ndatafields = len(datafields)
data = np.empty((ndatafields,),dtype=object)
for i,v in enumerate(data): 
    data[i] = []         # array of empty lists


# trigonal unit cell parameters
# a = lx
# b = np.sqrt(ly**2 + xy**2)
# c = np.sqrt(lz**2 + xz**2 + yz**2)
#
# cosa = (xy*xz + ly*yz)/(b*c)
# cosb = xz/c
# cosg = xy/b
#
# alpha = np.arccos(cosa)
# beta = np.arccos(cosb)
# gamma = np.arccos(cosg)
# vol_cell = volume(a,b,c,cosa,cosb,cosg)

for frame in range(node1.source.num_frames):

    print("running frame: {}".format(frame))
    
    node1.compute(frame)
    node2.compute(frame)
    tot_volume = node1.output.cell.volume
#    cell = node1.output.cell.matrix
#    lx = cell[0][0]
#    ly = cell[1][1]
#    lz = cell[2][2]
#    xy = cell[0][1]
#    xz = cell[2][0]
#    yz = cell[2][1]
#    
#    a = lx
#    b = np.sqrt(ly**2 + xy**2)
#    c = np.sqrt(lz**2 + xz**2 + yz**2)
#
#    cosa = (xy*xz + ly*yz)/(b*c)
#    cosb = xz/c
#    cosg = xy/b
#
#    alpha = np.arccos(cosa)
#    beta = np.arccos(cosb)
#    gamma = np.arccos(cosg)
#    vol_cell = volume(a,b,c,cosa,cosb,cosg)
#    print(lx,ly,lz,xy,xz,yz)
#    print(tot_volume)
#    print(vol_cell)
    
    # -- outer sphere
    vol_outer = node1.output.attributes['ConstructSurfaceMesh.solid_volume']
    area_outer = node1.output.attributes['ConstructSurfaceMesh.surface_area']
    n_selected1 = node1.output.attributes['SelectExpression.num_selected']
    r1_outer = sphere_vol2r(vol_outer)
    r2_outer = sphere_area2r(area_outer)
    
    # -- inner sphere
    vol_inner = node2.output.attributes['ConstructSurfaceMesh.solid_volume']
    area_inner = node2.output.attributes['ConstructSurfaceMesh.surface_area'] 
    n_selected2 = node2.output.attributes['SelectExpression.num_selected']
    r1_inner = sphere_vol2r(vol_inner)
    r2_inner = sphere_area2r(area_inner)
    
    # -- number of particles in system
    n_inner = n_selected2                           # particles in inner sphere
    n_outer = nincluded - n_selected1 - n_inner     # particles in outer shell
    n_rest = nincluded - n_outer - n_inner          # particles in outside shell
    n_tot = n_rest + n_outer + n_inner              # total number of particles in system
    ntest = nincluded - n_tot                       # test if n_tot adds up
    if (ntest!=0):
        print("Number of particles does not add up! dN = {}".format(ntest))
        print("particles tot: {} a+b+c: {} a: {} b: {} c: {}".format(nincluded,n_tot,n_rest,n_outer,n_inner))
    
    # -- volumes
    vol_rest = tot_volume - vol_outer
    vol_shell = vol_outer - vol_inner
    v_tot = vol_rest + vol_shell + vol_inner
    vtest = tot_volume - v_tot
    if (np.floor(vtest)!=0):
        print("Volumes in system does not add up! dV = {}".format(vtest))
        print("volume       : {} a+b+c: {} a: {} b: {} c: {}".format(tot_volume,v_tot,vol_rest,vol_shell,vol_inner))
    
    nrho_outer = n_outer/vol_shell                  # number density of outer shell
    nrho_inner = n_inner/vol_inner                  # number density of inner sphere
    nrho_rest = n_rest/vol_rest                     # number density of particles not in nanoparticle inner or outer shell
    nrho_system = nincluded/tot_volume              # number density of system
    
    tmpdata = [vol_inner, vol_shell, area_inner, area_outer, r1_inner, r2_inner, r1_outer, r2_outer, nrho_inner, nrho_outer, nrho_rest, nrho_system]
    for i in range(ndatafields):
        data[i].append(tmpdata[i])


write_to_file(ofile,data,datafields)
    
    
# --------------------------------------------------------------------- #
