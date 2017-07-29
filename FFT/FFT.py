# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 17:54:08 2017

@author: goranbs

FFT.py

Load signal from data file and perform FFT.
Write FFT of signal to "FFT" file.
Write found amplitudes and corresponding frequencies to "ampfreq" file.

Useage:

>> FFT.py -i ifile.dat -c [1,2] -p
>> FFT.py -i ifile.dat -c [1,2]

The input file is given through the -i flag.
The columns to use is given through the -c flag.
The indexing of columns starts at 1 being the first column of the input file.
The columns are interpreted as: [t, f(t)], i.e. f(t) is the amplitude at time t.

"""

import numpy as np                         # numpy library
#from scipy.optimize import curve_fit       # curve fitting
import scipy.signal as argrelextrema       # scipy integration methods
from datetime import datetime              # time and date
import argparse                            # handle command line arguments
import sys

wa_parser = argparse.ArgumentParser(description='Perform FFT of signal. Write results to two separate files: FFT and ampfreq.')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs=1,
                       help='Input columnar data file. Data must be space or tab separated.')
wa_parser.add_argument('-c','--columns', metavar='[1,2]', type=str, nargs=1,
                       help='[x, f(x)] columns in input file to use for FFT.')
wa_parser.add_argument('-p', '--plot', action='store_true',
                       help='Plot results? Add the -p flag')
                       
args = wa_parser.parse_args()

# -- What columns in the input data files should be used?
cols = args.columns[0].replace(']', '')
cols = cols.replace('[', '')
cols = cols.split(',')
cols = [(int(atype)-1) for atype in cols] # Since we define the first column as column nr 1. In python, first column is nr 0.

for col in cols:
    if (col < 0):
        sys.exit("Error! Cannot use columns in data file(s): -c [{:d} {:d}]".format(cols[0]+1,cols[1]+1))
        
# --------------------------------------------------------------- #

def write_FFT(freq, amp, ifile='FFT.out'):
    """
    Write FFT signal to file
    """
    
    nf = len(freq)
    na = len(amp)
    if (nf != na):
        sys.exit("Error! When trying to write FFT output file: n(freq) != n(amplitude)")
    
    
    tmp = ifile.split('.')
    prefix = tmp[0]
    suffix = tmp[-1]
    outfile = prefix + '-FFT.' + suffix
    ofile = open(outfile,'w')
    
    now = datetime.now()
    tt = now.timetuple()
    header1 = '# FFT signal from input file: {} using columns: {} {}\n'.format(ifile,cols[0],cols[1])
    header2 = '# Produced date: {:02d}/{:02d}-{} {:02d}:{:02d}:{:02d}\n'.format(tt[2],tt[1],tt[0],tt[3],tt[4],tt[5])
    header3 = '# Frequency Amplitude\n'
    
    ofile.write(header1)
    ofile.write(header2)
    ofile.write(header3)
    
    for i in xrange(nf):
        line = ' {:.9f} {:.9f}\n'.format(freq[i],amp[i])
        ofile.write(line)
        
    ofile.close()
    # ----------------------- DONE write_FFT()

def write_ampfreq(freq,amp, ifile='freqamp.out'):
    """
    Write found amplitudes with corresponding frequencies, from highest to lowest amplitude.
    """
    tmp = ifile.split('.')
    prefix = tmp[0]
    suffix = tmp[-1]
    outfile = prefix + '-ampfreq.' + suffix
    ofile = open(outfile,'w')
    
    now = datetime.now()
    tt = now.timetuple()
    header1 = '# Found amplitudes and frequencies from signal: {} using columns: {} {} \n'.format(ifile,cols[0]+1,cols[1]+1)
    header2 = '# Produced date: {:02d}/{:02d}-{} {:02d}:{:02d}:{:02d}\n'.format(tt[2],tt[1],tt[0],tt[3],tt[4],tt[5])
    header3 = '# Amplitude Frequency\n'
    
    ofile.write(header1)
    ofile.write(header2)
    ofile.write(header3)
    
    for i in xrange(len(amp)):
        line = ' {:.9f} {:.9f}\n'.format(amp[i],freq[i])
        ofile.write(line)
        
    ofile.close()    
    # ------------------------ DONE write_ampfreq()
        
# --------------------------------------------------------------- #
    
# -- Load data from file:
data = np.loadtxt(args.inputfile[0], usecols=(cols[0],cols[1]))
xdata = data[:,0]
ydata = data[:,1]
dt = xdata[1]-xdata[0]  # time step

N = xdata[-1]           # signal length
#T = 1.0/N              # signal period
L = len(xdata)          # number of data points

Y = np.fft.fft(ydata)   # FFT of signal
Y2 = Y[1:L/2]           # use half of the signal
Y2 = 2*np.abs(Y2/L)     # normalize

t = xdata[1:L/2]            # time
f = np.fft.fftfreq(L,d=dt)  # frequencies
f = f[1:L/2]                # len(f) = len(Y2)
#print len(f), f.size, f.shape
#print len(Y2), Y2.size, Y2.shape

# -- find distinct amplitudes:
max_amplitude = max(Y2)
# We will here only care about amplitudes larger than:
threshold=0.1*max_amplitude
idx = argrelextrema.argrelmax(Y2,order=5)   # locate extremas
indexes = []
extremas = Y2[idx]
# Store values within threshold:
for index in idx[0]:
    extrema = Y2[index]
    if (extrema > threshold):
        indexes.append(index)

## -- write output to file:
write_FFT(f,Y2,args.inputfile[0])
write_ampfreq(f[indexes],Y2[indexes],args.inputfile[0])
#-------------------------------------------------------- #
# -- plotting?
if (args.plot):
    import matplotlib.pyplot as plt            # plotting
    
    plt.figure(1)
    plt.plot(xdata,ydata,'o-')

    plt.figure(2)
    plt.plot(f,Y2,linewidth=3)
    plt.plot(f[indexes],Y2[indexes],'ro',markersize=10)
    plt.hold(False)
    plt.xlabel('frequency')
    plt.ylabel('amplitude')
    plt.show()
    
#-------------------------------------------------------- #

#amplitudes = Y2[indexes]
#frequencies = f[indexes]
#print amplitudes
#print frequencies
# ------------------------------------------------------------------ #

# ---------------------- TEST -------------------------------------- #
# this test gives the right amplitude of the signal

run_test = False

if (run_test):
    import matplotlib.pyplot as plt            # plotting
    
    # -- Create signal from amplitudes and frequencies:
    A1=3
    A2=10
    A3=0.5
    omega1=300
    omega2=200
    omega3=60
    
    threshold = 0.4 # in this case it is easy to set the threshold...
    
    N = 1000        # sampling freq
    T = 1.0/N       # sampling period  
    L = 1500        # length of signal
    
    t = np.linspace(0,L,L)*T  # time
    f = (float(N)/L)*np.linspace(1,L/2,L/2-1) # frequencies
    
    S = A1*np.sin(2*np.pi*omega1*t) + A2*np.sin(2*np.pi*omega2*t) + A3*np.sin(2*np.pi*omega3*t) + np.random.normal(0,0.1*A1,len(t))
    
    Y = np.fft.fft(S)
    Y2 = Y[1:L/2]
    Y2 = 2*np.abs(Y2/L)
    
    idx = argrelextrema.argrelmax(Y2,order=5)
    
    indexes = []
    extremas = Y2[idx]
    for index in idx[0]:
        extrema = Y2[index]
        if (extrema > threshold):
            indexes.append(index)
    
    plt.figure()
    plt.plot(f,Y2)
    plt.hold(True)
    plt.plot(f[indexes],Y2[indexes],'ro',markersize=10)
    plt.hold(False)
    plt.xlabel('frequency')
    plt.ylabel('amplitude')
    plt.show()

    write_FFT(f,Y2,args.inputfile[0])
    write_ampfreq(f[indexes],Y2[indexes],args.inputfile[0])
    
    print "Amplitudes : ", Y2[indexes]
    print "Frequencies: ", f[indexes]
    
# ----------------------------------------------------------------- # EOF