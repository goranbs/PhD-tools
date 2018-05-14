# -*- coding: utf-8 -*-
"""
Created on Mon May 14 10:43:45 2018

@author: goranbs

compute, plot and write the histogram of a dataset

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ------------------------------------------------ #

if len(sys.argv) < 3:
    print "\nRead columnar data set to compute histogram"
    print "Useage: hist.py filename colnr [nbins]\n"
    print "lines ignored in 'filename' starts with: #"
    print "filename : data file with columnar data"
    print "colnr    : column to use in data set from file 'filename'"
    print "nbins    : number of bins used in histogram [default nbins=100]"
    sys.exit()

# ------------------------------------------------ #

name = sys.argv[1]
colnr = int(sys.argv[2]) - 1
try:
    nbins = int(sys.argv[3])
except:
    nbins = 100

#print name, colnr, nbins
    
def write_data(x,y,gauss_fit,popt):
    pref = name[:name.rfind(".")]
    outfile = pref + '.out'
    fopen = open(outfile,'w')
    
    header = "# hist.py name={} colnr={} nbins={}\n".format(name,(colnr+1),nbins)
    subhead1 = "# Gaussian fit: a0={} mean={} sigma={}\n".format(popt[0],popt[1],popt[2])
    subhead2 = "# x count gauss_fit\n"
    
    fopen.write(header)
    fopen.write(subhead1)
    fopen.write(subhead2)
    for i in xrange(len(x)):
        fopen.write("{:f} {:f} {:f}\n".format(x[i],y[i],gauss_fit[i]))
        
    fopen.close()
    
def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

# -- do the magic
data = np.loadtxt(fname=name,comments='#')    
data = data[:,colnr]

hist, bin_edges = np.histogram(data,bins=nbins,normed=False)

# -- initial guess for the gaussian
a_0=1
mean_0=1
sigma_0=1

#print np.shape(hist), np.shape(bin_edges)
# -- curve fit

delta=bin_edges[1]-bin_edges[0]
bin_centers=(bin_edges[:-1]+delta)

popt, pcov = curve_fit(gauss,bin_centers,hist,p0=[a_0,mean_0,sigma_0])

print '--'*20
print 'gaussian fit:'
print ['a0','mean','sigma']
print popt
print '--'*20

plotting='no'

if plotting=='yes':
    plt.figure()
    plt.hist(data, bins=nbins, normed=False, label='hist')
    plt.plot(bin_centers, gauss(bin_centers,*popt), '--', label='gauss')
    plt.xlabel('x')
    plt.ylabel('count')
    plt.legend()
    plt.show()
    
write_data(bin_centers,hist,gauss(bin_centers,*popt),popt)
# ------------------------------------------------------- EOF