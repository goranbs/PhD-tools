#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:05:41 2017

@author: goranbs

Integrate LAMMPS profile outputs :)

integrateProfile.py -f table.profile -c [a,b]

    Returns the integral F = integral(f(x)dx), where x is in column a,
    and f(x) is in column b (read: flag -c [a,b])
    
"""


import numpy as np                         # numpy library
import scipy.integrate as integrate        # scipy integration methods
import argparse                            # handle command line arguments
import matplotlib.pyplot as plt            # plotting
from datetime import datetime              # date and time


wa_parser = argparse.ArgumentParser(description='Perform integration of tabular values')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs='+',
                       help='input data file(s). Files must be space or tab separated, columnar data')
wa_parser.add_argument('-c','--columns', metavar='[1,2]', type=str, nargs=1,
                       help='list of [x,f(x)] columns in input file. F = integral(f(x)dx)')
wa_parser.add_argument('-q','--cumfactors', metavar='[q1,q2,...]', type=str, nargs=1,
                       help='individual integration factors [q1,q2,...], one for each input file. Default values are q=1')
wa_parser.add_argument('-o', '--output', metavar=('filename'), type=str, nargs=1, required=False,
                       help='output file name')

args = wa_parser.parse_args()

Ni = len(args.inputfile) # number of input files

cols = args.columns[0].replace(']', '')
cols = cols.replace('[', '')
cols = cols.split(',')
cols = [int(atype) for atype in cols]

if args.cumfactors is not None:
    cumf = args.cumfactors[0].replace(']', '')
    cumf = cumf.replace('[', '')
    cumf = cumf.split(',')
    cumf = [float(atype) for atype in cumf]
else:
    cumf = np.ones(Ni)
    
print cumf

# TODO: Raise error of len(cols) > 2
# TODO: Raise error of len(cumf) != len(args.inputfile)

########################################################

Ecums = np.zeros((Ni,),dtype=list)

xes = []      # titles for x values in each data file
fofxes = []   # titles for f(x) values in each data file
nchunks = []  # number of rows (size of chunk) in each data file

for n in range(Ni):
    # open input file for reading:

    name = args.inputfile[n]
    ifile = open(args.inputfile[n],'r')    
    
    ifile.readline()
    header1 = ifile.readline()
    header2 = ifile.readline()
    header3 = ifile.readline()
    header1 = header1.split()
    header2 = header2.split()
    header3 = header3.split()
    header1.pop(0)
    header2.pop(0)

    print "---------------------------------------------------------------"
    print name
    print header1
    print header2
    print header3
    
    Nchunks = int(header3[1])
    Ncount = int(header3[2])
    Ncols = len(header2)
    
    matrixdata = np.zeros((Nchunks, Ncols), dtype=float)
    
    for i in range(Nchunks):
        line = ifile.readline()
        line = line.split()
        
        for j in range(Ncols):
            matrixdata[i,j] = float(line[j])
        
    print ""
    print "Performing integration using:"
    print "x    = ", header2[cols[0]]
    print "f(x) = ", header2[cols[1]]
    xes.append(header2[cols[0]])
    fofxes.append(header2[cols[1]])
    nchunks.append(Nchunks)
    
    x = matrixdata[:,cols[0]]
    f = matrixdata[:,cols[1]]
    f2 = f*cumf[n]              # multiplied by input factor
    
    Etot = integrate.simps(f2,x)
    Ecum = integrate.cumtrapz(f2,x, initial=0)
    
    print "Etot(x) = ", Etot
    print ""

    plt.plot(x, Ecum, label=str(name))
    Ecums[n] = Ecum

Emean = np.mean(Ecums,axis=0)
plt.plot(x, Emean, label='E(x)')
plt.legend()
#plt.show()

##############################################################
# -- Write output file:
if args.output is None:
    ofilename = "output.dat"
else:
    ofilename = args.output[0]
ofile = open(ofilename, 'w')

# TODO: Test if all files has the same column title for x!
# TODO: Test if all files has the same ranges

def checkEqual2(iterator):
   return len(set(iterator)) <= 1

if not (checkEqual2(xes)):
    print "Warning!"
    
if not(checkEqual2(nchunks)):
    print "Warning! Not all files has the same number of rows"

# -- write output file header:
header = 'File gererated by integrateProfile.py: {:s}\n# -- Files used:\n'.format(datetime.now().isoformat(' '))
ofile.write(header)
header = '{:s} '.format(xes[0])
for i in range(Ni):
    name = args.inputfile[i]
    name += ' q={:.5f}\n'.format(cumf[i])
    ofile.write(name)
    header += 'F({:s}) '.format(fofxes[i])
header += 'Ftot({:s}) \n'.format(xes[0])
entries = '# -- Entries:\n'
ofile.write(entries)
ofile.write(header)

for i in range(Nchunks):
    line = '{:.10f} '.format(x[i])
    for j in range(Ni):
        nr = '{:.10f} '.format(Ecums[j][i])
        line += nr
    line += '{:.10f}\n'.format(Emean[i])
    ofile.write(line)
    
ofile.close()

    



    





