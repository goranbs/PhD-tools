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
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


####################  VARIABLES  ########################################
sigma_array = [40,50,60]    #  = 1 sigma = diameter of nanoparticle [nm]
sigma_max   = max(sigma_array) 

Rho     = 0.1/sigma_max     #
T       = 300               #  temperature
addshift = 1
#########################################################################

####################  PLOTTING GRAPHS  ##################################

plot_initial_data                   = False
plot_unitless                       = True

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


H = []
W = []
last_xdatas = []

for sigma in sigma_array:

    ####################  CONSTANTS  ########################################
    k       = 1.38064852            # [10^(-23)k = kB]
    fac     = 100*sigma**2/(k*T)    # conversion factor 10^(2) * [sigma^2 / kT]
    reff    = sigma/4.0             # Effective radius [nm]
    Reff    = reff/sigma            # Effective radius [sigma]
    #########################################################################
    
    ####################  ABOUT THE CONVERSION FACTOR #######################
    #
    # input data in:
    #   [mJ/m^2 = 10^(-3) * 10^(-18) J/nm^2 = 10^(-21) J/nm^2]
    # input data to unitless: 
    #   [J/nm^2] * [sigma^2 / kBT]
    #   10^(-21) * 10^(23) [J/nm^2] * [sigma^2 / kT]
    #   100 * [J/nm^2] * [sigma^2 / kT]
        
    # 1) make unitless!
    xi = data[:,0]
    yi = data[:,1]
    xshift = xi[0]
    xu = (xi-xshift)/(sigma*10)
    yu = (2*np.pi*Reff*fac)*yi
    
    ## -- 2) shift data set to zero
    xu = xu-xu[0]
    xdata_last = xu[-1]
    last_xdatas.append(xdata_last)
    
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

    dh=h[1]-h[0]

    U = -integrate.cumtrapz(y_interp,dx=dh)
    U = np.append(U,U[-1])
    U = U - U[-1]    
    
    W.append(U)
    H.append(h+addshift)

    
##############################################################################
#############################################################################
# -- plot the smooth potential for different particle sizes (sigma)

## -- color scheme
Red     = 1.0
Green   = 0.2
Blue    = 0.9
alpha   = 0.8

linewidth=1.5
nsigma = len(sigma_array)

# -- sphere-sphere Force plot:
lineList = []

minU = np.min(W)

if (plot_unitless):
    plt.figure()
    for j in xrange(nsigma):
        red=Red/np.sqrt(nsigma-j+1)
        blue=Blue/np.sqrt(j+1)
        green=Green
        plt.plot(H[j],W[j], '--', color=(red,green,blue,alpha),linewidth=linewidth)
        data_key = mlines.Line2D([],[],color=(red,green,blue,alpha), label=r'$\sigma = $ {}'.format(sigma_array[j]),linewidth=linewidth)
        lineList.append(data_key)
    
    plt.xlim( (addshift-0.003, xdata_last+addshift) )
    plt.ylim( (minU, -0.5*minU) )
    plt.xlabel('h* [$\sigma$]')
    plt.ylabel('U*(h*) [$k_BT$]')
    plt.legend(bbox_to_anchor=(1.0, 0.0, 0.5, 1.0), handles=lineList, borderaxespad=0, mode='expand',loc='lower left', frameon=False)
    plt.tight_layout(rect=[0,0,0.78,1.0])
    plt.savefig('sigma.pdf', format='pdf')
    plt.show() 


