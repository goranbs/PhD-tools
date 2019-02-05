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
from matplotlib import pyplot as plt, rc
#import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'size'   : 16}

rc('font',**font)


####################  VARIABLES  ########################################
sigma   = 60        #  = 1 sigma = diameter of nanoparticle [nm]
rho     = 0.3        #  distribution of surface roughness [nm]
T       = 300       #  temperature
#########################################################################

####################  CONSTANTS  ########################################
k       = 1.38064852            # [10^(-23)k = kB]
fac     = 100*sigma**2/(k*T)    # conversion factor 10^(2) * [sigma^2 / kT]
reff    = sigma/4.0             # Effective radius [nm]
Reff    = reff/sigma            # Effective radius [sigma]
Ds      = 0.3/sigma            # descrete surface heights/steps [sigma]
Rho     = rho/sigma             # distribution of surface roughness [sigma]
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

calcite_data = np.loadtxt('./data/calcite_inregistry.dat', comments='#', usecols=(0,5), skiprows=1)
oCal=15.26
Acal=1710.282084
Joule=4.1868
Na=6.022
fCal=1000*Joule/Na/Acal
if (plot_initial_data):
    plt.figure()
    plt.plot(data[:,0], data[:,1], 'b--o')
    plt.plot(calcite_data[:,0]-oCal, calcite_data[:,1]*fCal, 'r--')
    plt.ylabel('$W(h)$ [$mJ/m^2$]')
    plt.xlabel('$h$ [$\AA$]')
    plt.xlim((data[0,0]-1, data[-1,0]))
    plt.ylim((min(data[:,1]),max(data[:,1])))
    plt.savefig('separation-free-energy.pdf',format='pdf',bbox_inches = "tight")
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
nh = 400    # number of data points in interpolated data
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
    #nrandom=1000
    #random_nr = np.random.normal(0,rho,nrandom)
    #random_shift = np.round(random_nr/Ds)*Ds
    #nbins = int(10*rho/(Ds*sigma))
    #print nbins
    #plt.hist(random_shift,bins=nbins,normed=True,rwidth=0.5)
    ds = Ds*sigma
    histy, index = np.unique(rGauss, return_index=True)
    histx = Z[index]
    histx = np.append(histx+0.5*ds,-histx-0.5*ds)
    histy = np.append(histy,histy)
    plt.bar(histx,histy,width=ds*0.4, align='center')
    plt.title(r'$\rho$ = {} nm'.format(rho))
    plt.xlabel('h [nm]')
    plt.ylabel('$\phi$')
    plt.savefig('histogram.pdf',format='pdf',bbox_inches = "tight")
    plt.show()
    
#############################################################################
#
# W(h) : Free energy of separation of two flat surfaces (already scaled by 2*pi*Reff)
# F(h) : Force between two spherical particles: F(h) = 2*pi*Reff*W(h)
# U(h) : Potential energy between two sphereical particles of radius sigma

addshift = 0
h = h+addshift

dh=h[1]-h[0]

F = y_interp
U = -integrate.cumtrapz(F,dx=dh)
U = np.append(U,U[-1])
U = U - U[-1]

if (plot_force):
    plt.figure()
    plt.plot(h,F, 'b--', label='$F(h) = 2\pi R_{eff}W(h)$')
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('F*(h*) [$k_BT/\sigma$]')
    plt.title('Force sphere-sphere')
    plt.legend()
    fmin=min(F)
    plt.ylim((fmin,-1.5*fmin))
    plt.xlim((addshift,xdata_last+addshift))
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

convF   = np.zeros(nh)     # convolution on F    
convU   = np.zeros(nh)     # convolution on U
for i in xrange(nh):
    j = i
    xj = h[j]   
    g = Gaussian(h, Rho, a=xj)
    convF[i] = sum(F*g)
    convU[i] = sum(U*g)


F_rough = convF*dh
U_rough = convU*dh



# -- Create nrandom shifts in h, drawn from a gaussian distribution
# -- find the force and the potential corresponding to this shift (using interpolation)

nrandom = 1000
rh      = np.zeros((nh,nrandom))
rFl      = np.zeros((nh,nrandom))
rUl      = np.zeros((nh,nrandom))

random_nr = np.random.normal(0,Rho,nrandom)
random_shift = np.round(random_nr/Ds)*Ds

for i in xrange(nrandom):
    rh[:,i] = h+random_shift[i]
    rFl[:,i] = np.interp(rh[:,i],h,F)
    rUl[:,i] = np.interp(rh[:,i],h,U)

rF = np.mean(rFl,axis=1)
rU = np.mean(rUl,axis=1)    

if (plot_interpl_and_rough_force):
    plt.figure()
    plt.plot(h,F, 'b-', label='$F(h) = 2\pi R_{eff}W(h)$')
    plt.plot(h,F_rough, 'r--', label=r'smooth $\rho$ = {:3f} $\sigma$'.format(Rho))
    plt.plot(h,rF, 'y--', label=r'rough $\rho$ = {:3f} $\sigma$'.format(Rho))
    #plt.plot(h,F_disc, 'y-.', label=r'discrete $\rho$ = {:3f} $\sigma$'.format(Rho))
    #plt.plot(rh_unique,rcF, 'g-', label=r'discrete $\rho$ = {:3f} $\sigma$'.format(Rho))
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('F*(h*) [$k_BT/\sigma$]')
    plt.title('Force sphere-sphere')
    plt.legend()
    fmin=min(F)
    plt.ylim((fmin,-1.5*fmin))
    plt.xlim((addshift,xdata_last+addshift))
    plt.show()

if (plot_interpl_and_rough_potential):
    plt.figure()
    plt.plot(h,U, 'b-', label='original')
    plt.plot(h,U_rough, 'r--', linewidth=3, label=r'smooth $\rho$ = {:3f} $\sigma$'.format(Rho))
    plt.plot(h,rU, 'y--', linewidth=3, label=r'rough $\rho$ = {:3f} $\sigma$'.format(Rho))
    #plt.plot(h,U_disc, 'y-.', linewidth=3, label=r'discrete $\rho$ = {:3f} $\sigma$'.format(Rho))
    #plt.plot(rh_unique,rcU, 'g-', linewidth=4, label=r'discrete $\rho$ = {:3f} $\sigma$'.format(Rho))
    #plt.plot(h,Gaussian(h,Rho,a=0.002),'g--')
    plt.legend()
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('U*(h*) [$k_BT$]')
    plt.title('Roughened and smooth potential')
    plt.xlim((h[0],h[-1]))
    #plt.xlim((-0.002,0.02))
    #plt.xlim((addshift,xdata_last+addshift))
    #plt.ylim((-130,130))
    plt.ylim((-20,20))
    plt.show()
    

##############################################################################