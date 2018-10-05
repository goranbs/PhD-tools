#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 19:12:10 2018

@author: goranbs

quick integral

# specify number of lines to skip
# integrate over chosen columns

"""

import numpy as np                         # numpy library
import scipy.integrate as integrate        # scipy integration methods
import matplotlib.pyplot as plt            # plotting
import sys

plot_result=True

Pa=101325.0         # conversion factor to [Pa]
GPa=0.000101325     # conversion factor to [GPa]

filename = sys.argv[1]
skipNrows = int(sys.argv[2])
X = int(sys.argv[3])-1           # python starts counting at zero!
Y = int(sys.argv[4])-1           # python starts counting at zero!

#print skipNrows, X, Y

data = np.loadtxt(filename, comments='#', skiprows=skipNrows)

x = data[:,X]
x = x[x<0.2]
y = -data[:,Y]*GPa          # NB! flipping sign, and myltipy by GPa
IDmax = np.argmax(y)

#print IDmax

x_ = x[0:IDmax]
y_ = y[0:IDmax]

Itot_ = integrate.simps(x_,y_)
Icum_ = integrate.cumtrapz(x_,y_, initial=0)

print Itot_

if (plot_result):

    Itot = integrate.simps(x,y)
    Icum = integrate.cumtrapz(x,y, initial=0)

    plt.figure()
    #plt.plot(x, y, 'b-o', label='data')
    plt.plot(x, Icum, 'rx--', label='I(data)={}'.format(Itot))
    plt.hold(True)
    plt.plot(x_, Icum_, 'bx--', label='I(data)={}'.format(Itot_))
    plt.hold(False)
    plt.xlabel('x')
    plt.ylabel('GJ/m^3')
    plt.legend()
    plt.show()


