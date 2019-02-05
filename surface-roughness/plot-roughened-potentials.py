# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 17:47:17 2018

@author: goranbs

Load free energy profile (integral of disjoining pressure)

1) plot unitless sphere-sphere potential using the derjaguin approximation

2) use a gaussian convolution to produce a roughened potential

3) use a discrete gaussian convolution to produce more realistic potential
"""


import numpy as np
import scipy.integrate as integrate
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


####################  VARIABLES  ########################################
sigma       = 60        #  = 1 sigma = diameter of nanoparticle [nm]
T           = 300       #  temperature

rho_array   = [0.03, 0.06, 0.09, 0.12, 0.15, 0.18] # [nm]
rho         = max(rho_array) # [nm]
#########################################################################

####################  CONSTANTS  ########################################
k       = 1.38064852            # [10^(-23)k = kB]
fac     = 100*sigma**2/(k*T)    # conversion factor 10^(2) * [sigma^2 / kT]
reff    = sigma/4.0             # Effective radius [nm]
Reff    = reff/sigma            # Effective radius [sigma]
Ds      = 0.3/sigma            # descrete surface heights/steps [sigma]

Rho_array = np.array(rho_array)/sigma # distribution of surface roughness [sigma]
Rho     = max(Rho_array)        # max roughness
#########################################################################

####################  ABOUT THE CONVERSION FACTOR #######################
#
# input data in:
#   [mJ/m^2 = 10^(-3) * 10^(-18) J/nm^2 = 10^(-21) J/nm^2]
# input data to unitless: 
#   [J/nm^2] * [sigma^2 / kBT]
#   10^(-21) * 10^(23) [J/nm^2] * [sigma^2 / kT]
#   100 * [J/nm^2] * [sigma^2 / kT]

####################  PLOTTING GRAPHS  ##################################

plot_initial_data                   = True
plot_unitless_initial_data          = True
plot_gauss_smooth                   = True
plot_force                          = False
plot_potential                      = False
plot_interpl_and_rough_potential    = True
plot_interpl_and_rough_force        = True
#########################################################################


####################  INPUT DATA ########################################
#
# data is the free energy per unit area [mJ/m^2] = 10^(-3) * [J/m^2]
#
# the data was produced by taking the disjoining pressure, 
# 1) using MATLAB grabit, to add more points to the curve
# 2) taking the integral of the curve using grace

data = np.loadtxt('./data/freenergy.txt', comments='#', usecols=(0,1))

if (plot_initial_data):
    plt.figure()
    plt.plot(data[:,0], data[:,1], 'b--o')
    plt.ylabel('W(h) [mJ/m^2]')
    plt.xlabel('h [angstrom]')
    plt.show()
    
# 1) make unitless!
xi = data[:,0]
yi = data[:,1]
xshift = xi[0]
xu = (xi-xshift)/(sigma*10)
yu = (2*np.pi*Reff*fac)*yi

## -- 2) shift data set to zero
xu = xu-xu[0]
xdata_last = xu[-1]

## -- 4) fit closest data points to a repulsive linear function
def repulsive(r,A,B):
    return A*r + B

A = (yu[1] - yu[0])/(xu[1]-xu[0])
B = xu[0]

A2 = 1*A    # use a larger slope of the repulsive linear function
#A2 = -4.2*10**(9)
### -- Youngs modulus of calcite: 80 GPa = 4.2*10**(9) sigma (sigma = 60 nm)

# extra data points (x0,y0) and (xn,yn)
# these points are added 5*Rho to the left and to the right of the original data points, respectively
x0 = -5*Rho
y0 = repulsive(x0,A2,B)
xn = xu[-1] + 5*Rho
yn = 0

x = np.zeros(len(xu)+2)
y = np.zeros(len(xu)+2)
x[0] = x0
x[1:-1] = xu
x[-1] = xn 
y[0] = y0
y[1:-1] = yu
y[-1] = yn

## -- 5) need to make the same discretization of the data. Use interpolation.

# -- a) create data points in line:
nline=5    # number of datapoints in line fit
xline = np.linspace(-5*Rho,B,nline)
yline = repulsive(xline,A2,B)

#plt.figure()
#plt.plot(xline,yline, 'y-o')
#plt.show()

## -- interpolate full data curve
nh = 1000    # number of data points in interpolated data
h = np.linspace(-5*Rho,x[-1]+5*Rho,nh)
y_interp = np.interp(h,x,y)

if(plot_unitless_initial_data):
    plt.figure()
    plt.plot(x,y,'r-', label='W(h)')
    plt.plot(xline,yline, 'yo', label='linear repulsive')
    plt.plot(h,y_interp, 'b--', label='W(h) + linear repulsive')
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('$2\pi R_{eff}$ W*(h*) [$k_BT/\sigma$]')
    plt.legend()
    plt.xlim((h[0],h[-1]))
    plt.show()



#############################################################################
#############################################################################
# Gaussian distribution

def Gaussian(x, sigma=1.0, a=0.0):
    return np.exp(-((x-a)/sigma)**2/2)/(sigma*np.sqrt(2*np.pi))

if (plot_gauss_smooth):
    z = np.linspace(-Rho*5,Rho*5,len(h))
    dz = z[1]-z[0]
    rz = np.round(z/Ds)*Ds
    gauss = Gaussian(z, Rho, a=0.0)
    rgauss = Gaussian(rz, Rho, a=0.0)
    rgauss = rgauss/sum(rgauss)/dz
    Z = z*sigma
    dZ = dz*sigma
    rZ = rz*sigma
    Gauss = Gaussian(Z, rho, a=0.0)
    rGauss = Gaussian(rZ, rho, a=0.0)
    rGauss = rGauss/sum(rGauss)/dZ
    A1 = integrate.simps(gauss,z)
    A2 = integrate.simps(Gauss,Z)
    plt.figure()
    plt.plot(z,gauss, 'r--', label='smooth')
    plt.plot(z,rgauss, 'b--', label='rough')
    plt.title('Area under curve = {}'.format(A1))
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('$\phi$')
    plt.show()
    plt.figure()
    plt.plot(Z,Gauss, 'r--', label='smooth')
    plt.plot(Z,rGauss, 'b--', label='rough')
    nrandom=1000
    random_nr = np.random.normal(0,rho,nrandom)
    random_shift = np.round(random_nr/Ds)*Ds
    plt.hist(random_shift,bins=5,normed=True)
    plt.title('Area under curve = {}'.format(A2))
    plt.xlabel('h [nm]')
    plt.ylabel('$\phi$')
    plt.show()

#############################################################################    
#############################################################################
# plot the smooth force and smooth potential (sphere-sphere)
# W(h) : Free energy of separation of two flat surfaces (already scaled by 2*pi*Reff)
# F(h) : Force between two spherical particles: F(h) = 2*pi*Reff*W(h)
# U(h) : Potential energy between two sphereical particles of radius sigma

addshift = 1

dh=h[1]-h[0]

F = y_interp
U = -integrate.cumtrapz(F,dx=dh)
U = np.append(U,U[-1])
U = U - U[-1]

fmin=min(F)
ifmin = np.argmin(F)
delta = h[ifmin]
#print ifmin, delta, fmin, F[ifmin]
hshift = addshift - delta

h = h + hshift

if (plot_force):
    plt.figure()
    plt.plot(h,F, 'b--', label='$F(h) = 2\pi R_{eff}W(h)$')
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('F*(h*) [$k_BT/\sigma$]')
    plt.title('Force sphere-sphere')
    plt.legend()
    plt.ylim((fmin,-1.5*fmin))
    plt.xlim((hshift,xdata_last+hshift))
    plt.show()
    
if (plot_potential):
    plt.figure()
    plt.plot(h,U, 'b--', label='$U(h) = -\int \, F(h)\, dh + U_0$')
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('U*(h*) [$k_BT$]')
    plt.title('Potential sphere-sphere')
    plt.legend(loc='best')
    plt.xlim((h[0],h[-1]))
    plt.show()


#############################################################################    
#############################################################################
# compute the roughnened forces and the roughened potentials (sphere-sphere)

nrough  = len(Rho_array)

#-- Roughened potentials/forces from a smooth gaussian convolution:
convF   = np.zeros((nh,nrough))     # convolution on F    
convU   = np.zeros((nh,nrough))     # convolution on U

# -- Roughened potentials/forces from descrete surface heights
nrandom = 2000
rh      = np.zeros((nh,nrandom,nrough))
rFl     = np.zeros((nh,nrandom,nrough))
rUl     = np.zeros((nh,nrandom,nrough))

for j in xrange(nrough):
    roughness = Rho_array[j]
    
    # -- gaussian convolution
    for i in xrange(nh):
        xi = h[i]   
        g = Gaussian(h, roughness, a=xi)
        convF[i,j] = sum(F*g)
        convU[i,j] = sum(U*g)
    
    # -- Create nrandom shifts in h, drawn from a gaussian distribution
    # -- find the force and the potential corresponding to this shift (using interpolation)
    random_nr = np.random.normal(0,roughness,nrandom)
    random_shift = np.round(random_nr/Ds)*Ds
    
    for i in xrange(nrandom):
        rh[:,i,j] = h+random_shift[i]
        rFl[:,i,j] = np.interp(rh[:,i,j],h,F)
        rUl[:,i,j] = np.interp(rh[:,i,j],h,U)


convF   = convF*dh
convU   = convU*dh

rF      = np.zeros((nh,nrough))
rU      = np.zeros((nh,nrough))
for j in xrange(nrough):
    rF[:,j] = np.mean(rFl[:,:,j], axis=1)
    rU[:,j] = np.mean(rUl[:,:,j], axis=1)

#print np.shape(rF)
#print np.shape(convF)

#############################################################################    
#############################################################################
# -- plot the roughened potentials and forces

l1 = mlines.Line2D([],[], color='black', label='Continuous', linestyle='--', linewidth=3)
l2 = mlines.Line2D([],[], color='black', label='Discrete', linestyle=':', linewidth=3)

## -- color scheme
Red     = 1.0
Green   = 0.2
Blue    = 0.9
alpha   = 0.8

linewidth=1.5

# -- sphere-sphere Force plot:
lineList_forces = []
if (plot_interpl_and_rough_force):
    plt.figure()
    minF = min(F)
    plt.plot(h,F,'b-',linewidth=3)
    data_key = mlines.Line2D([],[],color='blue', label=r'$\rho =$ 0.0 nm',linewidth=3)
    lineList_forces.append(data_key)
    for j in xrange(nrough):
        red=Red/np.sqrt(nrough-j+1)
        blue=Blue/np.sqrt(j+1)
        green=Green
        plt.plot(h,convF[:,j], '--', color=(red,green,blue,alpha),linewidth=linewidth)
        plt.plot(h,rF[:,j], ':', color=(red,green,blue,alpha),linewidth=linewidth)
        data_key = mlines.Line2D([],[],color=(red,green,blue,alpha), label=r'$\rho = $ {} nm'.format(rho_array[j]),linewidth=linewidth)
        lineList_forces.append(data_key)
    
    plt.xlim( (addshift, xdata_last+hshift) )
    plt.ylim( (minF, -1.5*minF) )
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('F*(h*) [$k_BT/\sigma$]')
    plt.title('$\sigma $ = {} nm'.format(sigma))
    lineList_forces.append(l1)
    lineList_forces.append(l2)
    plt.legend(bbox_to_anchor=(1.0, 0.0, 0.5, 1.0), handles=lineList_forces, borderaxespad=0, mode='expand',loc='lower left', frameon=False)
    plt.tight_layout(rect=[0,0,0.78,1.0])    
    plt.savefig('force-sphere-sphere.pdf', format='pdf')
    plt.show()

############### This is not a valid shift, however, need to make the plot look nice:
h = h+3*delta
###############

# -- sphere-sphere Potential plot:
lineList_potential = []
if (plot_interpl_and_rough_potential):
    plt.figure()
    minU = min(U)
    plt.plot(h,U,'b-',linewidth=3)
    data_key = mlines.Line2D([],[],color='blue', label=r'$\rho =$ 0.0 nm',linewidth=3)
    lineList_potential.append(data_key)
    for j in xrange(nrough):
        red=Red/np.sqrt(nrough-j+1)
        blue=Blue/np.sqrt(j+1)
        green=Green
        plt.plot(h,convU[:,j], '--', color=(red,green,blue,alpha),linewidth=linewidth)
        plt.plot(h,rU[:,j], ':', color=(red,green,blue,alpha),linewidth=linewidth)
        data_key = mlines.Line2D([],[],color=(red,green,blue,alpha), label=r'$\rho = $ {} nm'.format(rho_array[j]),linewidth=linewidth)
        lineList_potential.append(data_key)        
    
    plt.xlim( (addshift, xdata_last+hshift) )
    plt.ylim( (minU, -1.5*minU) )    
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('U*(h*) [$k_BT$]')
    plt.title('$\sigma$ = {} nm'.format(sigma))
    lineList_potential.append(l1)
    lineList_potential.append(l2)
    plt.legend(bbox_to_anchor=(1.0, 0.0, 0.5, 1.0), handles=lineList_potential, borderaxespad=0, mode='expand',loc='lower left', frameon=False)
    plt.tight_layout(rect=[0,0,0.78,1.0])    
    plt.savefig('potential-sphere-sphere.pdf', format='pdf')
    plt.show()
    

##############################################################################