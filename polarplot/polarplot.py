# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 13:21:22 2018

@author: goranbs (Gøran Brekke-Svaland)

Description:

This script plots the center-to-surface interaction energy U(r,theta),  of
a nanoparticle with diameter "sigma" and interaction strength "epsilon".

A polar plot is produced for a number of surface domains (slices), and the
average U(r) is plotted in a separate plot.

The global roughness of the nanoparticle is "kappa", and the domain step size
is "stepsize", which is the descrete kink-heights of the nanoparticle.

"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# --------------------

slices = 200          # number of surface slices/domains
kappa=3               # global surface roughness [nm]
stepsize=0.3          # step size [nm]

# --------------------

epsilon=130.0       # interaction strength
sigma=60.0          # particle diameter [nm]
max_sigma=60        # max value of r, in polar plot
a=24                # mie exponent
# --------------------

resolution = 200    # r,theta resolution

# -------------------------------------------- #

rad = np.linspace(0, max_sigma, resolution)
azm = np.linspace(0, 2*np.pi, resolution)

r, th = np.meshgrid(rad, azm)

delta = np.zeros((resolution,resolution))
npoints = resolution/slices


sumkappasquared = 0
stepsizes = []
for i in xrange(slices):
    randomNr = np.random.normal(0,kappa)    # Gaussian normal distr. with std=kappa
    k = stepsize*np.round(randomNr/stepsize)
    sumkappasquared += k*k
    stepsizes.append(k)
    for j in xrange( npoints*(i-1), npoints*i ):
        delta[j,:] = k

roughness_rms = np.sqrt(sumkappasquared/slices)
roughness_std = np.std(delta[:,1])

# ---------------------------------------- #
print "Global roughness"
print "rms roughness:  ", roughness_rms
print "std roughness:  ", roughness_std
# ---------------------------------------- #

# -- mie potential, center-to-surface distance (r): U(r,theta)
z = 4*epsilon*((sigma/(r+delta+sigma/2.0))**(a) - (sigma/(r+delta+sigma/2.0))**(a/2.0))

def stepsizes_hist(stepsizes):
    nbins=10
    #np.histogram(stepsize,bins=nbins,normed=False)
    plt.hist(stepsizes,bins=nbins,normed=False)
    plt.xlim([-nbins,nbins])
    plt.ylabel('number of steps')
    plt.xlabel('step size [nm]')
    plt.show()

def average_potential_to_1d(z_rtheta,rad,epsilon,max_sigma,slices,kappa):
    z = np.average(z_rtheta,axis=0)
    plt.figure()
    plt.plot(rad,z)
    plt.axis([(sigma/6), max_sigma, -epsilon, epsilon])
    plt.xlabel(r'$r \, / \, [nm]$')
    plt.ylabel(r'$U_r(r,\theta ) \, / \, [k_BT]$')
    name = "avgInteraction_n{}_k{}.eps".format(slices,kappa)
    plt.savefig(name,format='eps')
    plt.show()

def plotting_stuff(th,r,z,azm,epsilon,slices,kappa):
    plt.figure()
    #ax = Axes3D(fig)
    plt.subplot(projection="polar")
    plt.pcolormesh(th, r, z, vmin=-epsilon, vmax=epsilon, cmap='jet')
    plt.plot(azm, r, color='k', ls='none')
    cbar = plt.colorbar()
    cbar.set_label(r'$U_r(r,\theta ) \, / \, [k_BT]$')
    plt.grid()
    name = "polarplot_n{}_k{}.eps".format(slices,kappa)
    plt.savefig(name,format='eps')
    plt.show()
    
stepsizes_hist(stepsizes)
plotting_stuff(th,r,z,azm,epsilon,slices,kappa)
average_potential_to_1d(z,rad,epsilon,max_sigma,slices,kappa)
# ------------------------------------------------------------- #