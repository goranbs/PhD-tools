# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:08:10 2018

@author: goranbs

Compute the diffusion coefficient from the average velocity components v = (vx,vy,vz)

## useage: >> python diffusion_vacf_py ifile.dat dt decorrelationtime nstep
## exampl: >> python diffusion_vacf_py ifile.dat 2.0 100 1000

## output: 
    1) diffusion-ifile.out  : Running average of diffusion coefficients
    2) vacf-ifile.out       : The average of the VACF's

The relation between the diffusion coefficient D, and the velocity autocorrelation function (VACF) is

    D = integral(VACF*dt)/d

where d is the dimensionality of the system (3D -> d=3) and the VACF is:

    VCAF = v0*vi

where v0 is the initial velocity, and vi is the velocity at time step i.

"""

import sys
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# useage: >> python diffusion_vacf_py ifile.dat 2.0 100 1000

# ------------------------------------- #
fname = str(sys.argv[1])                # input file name
dt = float(sys.argv[2])                 # time step
decorrelationtime = int(sys.argv[3])    # decrorrelation time of particle velocities
nstep = int(sys.argv[4])                # number of time steps to use in VACF
# ------------------------------------- #

#fname="vacf-c25-2.dat"      # input file name
#dt=2.0                      # time step used in simulation
#decorrelationtime = 100     # number of time steps before vcaf is decorrelated
#nstep=1000

x=1     # column of vx velocities in input file
y=2     # -//- vy
z=3     # -//- vz

data = np.genfromtxt(fname,comments="#")
ndatapoints = len(data[:,0])

nv0 = int(ndatapoints/decorrelationtime)
ni = nv0-nstep/decorrelationtime

time = np.linspace(0,nstep*dt,nstep)

vacfs = np.zeros((5,nstep))
Timestep = np.zeros(ni)

Itot = np.zeros(ni)     # Time evolution of total diffusion coefficient
Ixy = np.zeros(ni)      # Time evolution of diffusion coefficient in xy-plane
Iz = np.zeros(ni)       # Time evolution of diffusion coefficient in z

#plt.figure()
#plt.hold(True)
for i in range(ni):
    
    j = i*decorrelationtime
    k = j+nstep
    
    v0x = data[j,x]
    v0y = data[j,y]
    v0z = data[j,z]
    vx = np.array(data[j:k,x])
    vy = np.array(data[j:k,y])
    vz = np.array(data[j:k,z])
    vacf_vx = v0x*vx
    vacf_vy = v0y*vy
    vacf_vz = v0z*vz
    
    vacfs[0] += vacf_vx
    vacfs[1] += vacf_vy
    vacfs[2] += vacf_vz
    vacfs[3] += vacf_vx + vacf_vy + vacf_vz
    vacfs[4] += vacf_vx + vacf_vy
    
    if i > 0:
        vi = vacfs[3]/i
        vi_xy = vacfs[4]/i
        vi_z = vacfs[2]/i
        itoti = integrate.simps(vi,time)/3.0
        itoti_xy = integrate.simps(vi_xy,time)/2.0
        itoti_z = integrate.simps(vi_z,time)
        Itot[i] = itoti
        Ixy[i] = itoti_xy
        Iz[i] = itoti_z
        Timestep[i] = j
        
#    print np.shape(time)
#    print np.shape(vacf_vx)

#    plt.plot(time,vacf_vx)

#plt.hold(False)
#plt.show()

vacfx = vacfs[0]/ni
vacfy = vacfs[1]/ni
vacfz = vacfs[2]/ni
vacf = vacfs[3]/ni
vacf_xy = vacfs[4]/ni

ix = integrate.simps(vacfx,time)
iy = integrate.simps(vacfy,time)
iz = integrate.simps(vacfz,time)
itot = integrate.simps(vacf,time)/3.0
itot_xy = integrate.simps(vacf_xy,time)/2.0

header="step Dxyz Dxy Dz # final Dxyz={}, Dxy={}, Dx={}, Dy={}, Dz={}".format(itot,itot_xy,ix,iy,iz)
X = (Timestep,Itot,Ixy,Iz)
np.savetxt("diffusion-{}.out".format(fname), np.transpose(X), header=header)

header="step vacf vacf_x vacf_y vacf_z"
Y = (time, vacf, vacfx, vacfy, vacfz)
np.savetxt("vacf-{}.out".format(fname),np.transpose(Y),header=header)

plot_figure=False

if (plot_figure):
    plt.figure()
    plt.hold(True)
    plt.plot(time,vacfx, label='x')
    plt.plot(time,vacfy, label='y')
    plt.plot(time,vacfz, label='z')
    plt.plot(time,vacf, label='tot')
    plt.hold(False)
    plt.legend()
    plt.show()


# ---------------------- # How to produce the input for this script, using LAMMPS:

#group           one id 10 # just an example
#compute         1 one property/atom vx vy vz
#compute         2 one reduce sum c_1[1] c_1[2] c_1[3]
#thermo          1
#thermo_style    custom step c_2[1] c_2[2] c_2[3]
#
# or using ave/time:
#fix            AVE one ave/time 1 1 1 c_2[1] c_2[2] c_2[3] file ifile.dat

# ---------------------- #
    
