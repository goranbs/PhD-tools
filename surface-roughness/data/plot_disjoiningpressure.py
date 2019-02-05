# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 17:18:55 2018

@author: goranbs

"""

import numpy as np
#from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib import rc

font = {'size'   : 16}

rc('font',**font)

calcite_in_registry = np.loadtxt('calcite_inregistry.dat', comments='#', usecols=(0,3))
disjoiningpressure = np.loadtxt('disjoiningpressure.txt', comments='#', usecols=(0,1))

shift=15.3  # angstrom

x1 = calcite_in_registry[:,0]
x1 = x1 - shift
y1 = calcite_in_registry[:,1]
x2 = disjoiningpressure[:,0]
y2 = disjoiningpressure[:,1]


lw=1.5      # linewidth


plt.figure()
plt.plot(x1,y1,'r-', label='original dataset', linewidth=lw)
plt.plot(x2,y2, 'bo', label='extended dataset', linewidth=lw)
plt.xlabel('$h$ [$\AA$]')
plt.ylabel(r'$\Pi(h)$ [MPa]')
plt.legend()
plt.xlim((5,20))
plt.savefig('disjoiningpressure.pdf', format='pdf', bbox_inches = "tight")
plt.show()


