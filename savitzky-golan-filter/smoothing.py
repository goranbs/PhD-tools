#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Feb 15 18:24:30 2018

@author: goranbs (Gøran Brekke Svaland)

Description:

Apply Savitzky-Golay filter to noisy data.
Output the filtered data input.

* Input data file with columnar data: f(x) x
* Assuming comments in data file starts with '#'

Savitzky-Golay filter from: http://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------ #

if len(sys.argv) < 4:
    print "\nApply Savitzky-Golay smoothing filter to noisy data"
    print "Useage: smoothen.py filename colx coly [wsize] [polyn]\n"
    print "lines ignored in 'filename' starts with: #"
    print "filename : data file with columnar data"
    print "colx     : column along x-axis. First column in data set is 1"
    print "coly     : column along y-axis (f(x))."
    print "wsize    : window size of SG filter (optional). Positive odd number. wsize=1001 (default)"
    print "polyn    : polynomial order of SG filter (optional). polyn=3 (default)"
    sys.exit()

# ------------------------------------------------ #

name = sys.argv[1]
colx = int(sys.argv[2]) - 1
coly = int(sys.argv[3]) - 1
try:
    wsize = sys.argv[4]
except:
    wsize = 101
try:
    polyn = sys.argv[5]
except:
    polyn = 3

print name, colx, coly, wsize, polyn

data = np.loadtxt(fname=name,comments='#')
    
index = np.argsort(data[:,1])
sx = data[index,1]
sy = data[index,3]

# ------------------------------------------------ #

def write_data(x, y, name, colx, xoly, wsize, polyn):
    pref = name[:name.rfind(".")]
    outfile = pref + '.out'
    fopen = open(outfile,'w')
    
    header = "# smoothing.py name={} colx={} coly={} wsize={} polyn={}\n".format(name,(colx+1),(coly+1),wsize,polyn)
    subhead= "# colx smoothy\n"
    
    fopen.write(header)
    fopen.write(subhead)
    for i in xrange(len(x)):
        fopen.write("{:f} {:f}\n".format(x[i],y[i]))
        
    fopen.close()

def plot_data(x,y,yhat,colx,coly,wsize,polyn):
    plt.figure()
    plt.plot(x,y,'or', label='data')
    plt.plot(x,yhat,'-pb', label='smooth')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.title("colx={} coly={} wsize={} polyn={}".format((colx+1),(coly+1),wsize,polyn))
    plt.show()

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
# ------------------------------------------------ #

yhat = savitzky_golay(sy, wsize, polyn)     # run Savitzky-Golan filter
write_data(sx,sy,name,colx,coly,wsize,polyn)    # write output
plot_data(sx,sy,yhat,colx,coly,wsize,polyn)     # plot data and filtered data

# ------------------------------------------------ #
    
