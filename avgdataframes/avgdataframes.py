# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 11:48:27 2018

@author: goranbs


* columns with header 'std' are currently removed by hard code!! NB!!

"""

import numpy as np                         # numpy library
import argparse                            # handle command line arguments
import pandas as pd                        # python pandas
import matplotlib.pyplot as plt            # plotting

wa_parser = argparse.ArgumentParser(description='Read multiple input data files and perform averaging over data frame points. This assumes input files to have the same structure.')
wa_parser.add_argument('-i', '--inputfile', metavar=('filename'), type=str, nargs='+',
                       help='input files with columnar data')
wa_parser.add_argument('-o', '--output', metavar=('prefix'), type=str, nargs=1,
                       help='output prefix')                       
                       
args = wa_parser.parse_args()
                       
nfiles = len(args.inputfile[:])


# data file contains columns:
# nr Temp std v_Nconf std PotEng std c_Lcom[3] std c_Ucom[3] std c_Icom[1] std c_Icom[2] std c_Icom[3] std v_fe std v_fz std v_fse std v_fsz std v_fie std v_fiz std
# 0    1   2    3      4    5     6     7       8      9      10   11      12     13     14    15      16   17   18   19  20  21   22   23   24   25   26    27   28


dataframes = []
nrows = []
# -- read data files
for i in xrange(nfiles):

    with ( open(args.inputfile[i],'r') ) as ifile:
        
        header = ifile.readline().split()
        header = header[1:]
        ncols = len(header)
        dataframe = {}
        for j in xrange(ncols):
            dataframe.update({header[j] : [] })
        
        nlines = 0
        for line in ifile:
            nlines +=1
            line = line.split()
            for k in xrange(ncols):
                entry = header[k]
                dataframe[entry].append(float(line[k]))
        
        #dataframes.append(dataframe)
        dataframe.pop('std')    # remove std entry !! NB !!
        dataframes.append(pd.DataFrame(data=dataframe))
        nrows.append(nlines)

#print dataframe['nr']        
#print dataframe['Temp']

# -- perform averaging of data points accross dataframes

cat = pd.concat(dataframes)
dog = cat.groupby(cat.index)

mean = dog.mean()       # compute mean over data frames
std = dog.std()         # compute standard deviation over data frames

colnames = list(mean.columns.values)    # get all column names
stdnames = []
for colname in colnames:
    stdnames.append(colname + '_std')   # add 'std' to names of stdev columns

std.columns = stdnames                  # rename columns

fulltable = pd.concat([mean,std],axis=1)    # concatinate mean and stdev tables

# -- write data frames to file
if args.output is None:
    ofile = 'fulltable.txt'
else:
    ofile = args.output[0] + '.txt'    

#colnames = list(fulltable.columns.values)
#for name in colnames:
#    print name
fulltable.to_csv(ofile, sep='\t') 
#mean.to_csv('mean.txt', sep='\t')
#std.to_csv('std.txt', sep='\t')

# -- plotting:

# nr            1
#PotEng         2
#Temp           3
#c_Icom[1]      4 v
#c_Icom[2]      5 v
#c_Icom[3]      6 v
#c_Lcom[3]      7 v
#c_Ucom[3]      8 v
#nr             9
#v_Nconf        10 v
#v_fe           11
#v_fie          12
#v_fiz          13 v
#v_fse          14
#v_fsz          15 v
#v_fz           16 v
#PotEng_std     17
#Temp_std       18
#c_Icom[1]_std  19
#c_Icom[2]_std  20
#c_Icom[3]_std  21
#c_Lcom[3]_std  22
#c_Ucom[3]_std  23
#nr_std         24
#v_Nconf_std    25
#v_fe_std       26
#v_fie_std      27
#v_fiz_std      28 v
#v_fse_std      29
#v_fsz_std      30 v
#v_fz_std       31 v


#fulltable.plot(x='c_Lcom[3]',y='v_Nconf')
#plt.show()

#gnuplot> plot 'avg-Ca-2-3-4.txt' u ($8-$7):($16+$15+$13):30 w errorbars 
#gnuplot> replot 'avg-Ca-2-3-4.txt' u ($8-$7):13:28 w errorbars 

