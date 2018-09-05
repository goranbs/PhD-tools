# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 09:32:44 2018

@author: goranbs


Estimate the free energy barrier of a single ion moving from one side
of the confinement to the other.

beta*dF = -ln(p(z))

where p(z) is the probability density of finding an ion at a specific location

p(z) is found from histogramming the positions of the ion throughout the simulation.


a) How to find the center of the confinement?
    1) Identify the four tallest peaks in the density profile
    2) The center point between the two central peaks is the center of the confinement.

b) Compute histogram of ion positions and normalize it

c) Look at the symmetry of p(z) -> plots

d) compute total p(z) reflecting the particle positions around the center of the confinement

f) write data to file

g) BASH/gnuplotting/*makeplots-dF.sh are shell scripts for plotting output data

# -- data fields in log file:
# Step Temp v_Nconf PotEng c_Lcom[3] c_Ucom[3] c_Icom[1] c_Icom[2] c_Icom[3] v_fe v_fz v_fse v_fsz v_fie v_fiz

"""

import numpy as np
import matplotlib.pyplot as plt

resNr=8
molecule="Ca"
datafield=8

# -- plot stuff?
plot_water_densityprofile=True
plot_probability_density=True
plot_free_energy=True
show_plots=False
dpi=200

# ------------------------------------------------- #

for nr in reversed(range(13,26)):

    fname="../data/log.lammps-res{}-{}-{}.out".format(resNr,molecule,nr)
    data = np.loadtxt(fname,dtype=float)
    
    dname="../data/spce-res4-{}-{}.profile".format(molecule,nr)
    ddata = np.loadtxt(dname,dtype=float,skiprows=4)


    IDs = np.argpartition(ddata[:,3],-8)
    
    dmax = ddata[IDs[-8:],1]
    rmax = ddata[IDs[-8:],3]
    
    mask = np.argsort(dmax)
    Dmax = dmax[mask]
    Rmax = rmax[mask]
    
    delta_Dmax = Dmax[1:]-Dmax[:-1]
    print Dmax
    zpos = []
    zid = []
    for i in range(1,len(Dmax)):
        if ( Dmax[i] > (Dmax[i-1]+2.0)):
            if (len(zpos) < 2):
                zpos.append(Dmax[i])
                zid.append(i)
    
    nbins=100
    hnbins=nbins/2
    zmin=zpos[0]
    zmax=zpos[1]
    dz = (zmax-zmin)
    zmid = zmin+dz/2.0
    
    if (plot_water_densityprofile):
        plt.figure()
        plt.plot(ddata[:,1],ddata[:,3])
        plt.hold(True)
        plt.plot(Dmax[zid],Rmax[zid], 'o')
        plt.hold(False)
        plt.title("{} res{} sep{} zmid={}".format(molecule,resNr,nr,zmid))
        plt.savefig('water-densityprofile-res{}-{}-{}'.format(resNr,molecule,nr),dpi=dpi)
        if (show_plots):
            plt.show()
    
        print "zmin = ", zmin 
        print "zmax = ", zmax
        print "zmid = ", zmid
        print "dz = ", dz
        
    
    hist, bin_edges = np.histogram(data[:,datafield],bins=nbins,range=(zmin,zmax), density=True)
    bins = np.linspace(-dz/2.0, dz/2.0, nbins)
    
    ndatapoints = len(data[:,datafield])
    reflected_data = np.zeros(ndatapoints)
    for i in range(ndatapoints):
        if (data[i,8]<zmid):
            reflected_data[i] = 2*zmid - data[i,datafield] -zmid
        else:
            reflected_data[i] = data[i,datafield] -zmid
    
    # -- plot data points
    #plt.figure()
    #plt.plot(reflected_data, '.')
    #plt.hold(True)
    #plt.plot(data[:,datafield]-zmid, 'x')
    #plt.hold(False)
    #plt.show()
        
    nbins2 = nbins/2.
    h2, b2 = np.histogram(reflected_data,bins=nbins2,range=(0,dz/2.0), density=True)
    bins2 = np.linspace(0,dz/2.0,nbins2)
    
    # -- XKCD style plot
    #with plt.xkcd(scale=1,length=100,randomness=2):
    #    plt.plot(bins,hist)
    #    plt.show()
    
    if (plot_probability_density):
        plt.figure()
        plt.plot(bins[:hnbins],hist[:hnbins],'b-',label='left')
        plt.hold(True)
        plt.plot(bins[hnbins:],hist[hnbins:],'r-',label='right')
        plt.plot(-bins[:hnbins],hist[:hnbins],'m-',label='left reflected')
        plt.plot(-bins[hnbins:],hist[hnbins:],'c-',label='right reflected')
        plt.plot(bins2,h2,'bx--', label='total')
        plt.hold(False)
        plt.xlabel("z [$\AA{}$]")
        plt.ylabel("p(z)")
        #plt.legend()
        plt.title("{} res{} sep{}".format(molecule,resNr,nr))
        plt.savefig('pz-res{}-{}-{}'.format(resNr,molecule,nr),dpi=dpi)
        if (show_plots):
            plt.show()
    
    # ----------------- #
    
    T=300               # Temperature
    kB=1.3806e-23       # [kg*m^2/K/s^2]
    beta=1.0/(kB*T)     # 1/(kT)
    beta=1.0
    
    # -- 1) 
    dFl = -np.log(hist[:hnbins])/beta    # [mJ] // [kg*m^2/s^2]
    
    centerFl = dFl[hnbins-1]
    dFl = dFl-centerFl
    
    mdFl = np.min(dFl)
    idFl = np.where(dFl==mdFl)
    
    # -- 2)
    dFr = -np.log(hist[hnbins:])/beta    # [mJ] // [kg*m^2/s^2]
    
    centerFr = dFr[0]
    dFr = dFr-centerFr
    
    mdFr = np.min(dFr)
    idFr = np.where(dFr==mdFr)
    
    # -- 3) Total
    dF = -np.log(h2)/beta
    cF = dF[0]
    dF = dF - cF
    mF = np.min(dF)
    idF = np.where(dF==mF)
    
    if (plot_free_energy):
        plt.figure()
        plt.plot(bins[:hnbins],dFl,'b-', label='left')
        plt.hold(True)
        plt.plot(bins[hnbins:],dFr,'r-', label='right')
        plt.plot(-bins[:hnbins],dFl,'m-', label='left reflected')
        plt.plot(-bins[hnbins:],dFr,'c-', label='right reflected')
        plt.plot(bins[idFl[0]],mdFl,'bo', label='min left')
        plt.plot(bins[idFr[0]+hnbins],mdFr, 'ro', label='min right')
        plt.plot(bins2,dF,'bx--', label='total')
        plt.plot(bins2[idF[0]],mF,'rd', label='min total')
        plt.hold(False)
        plt.xlabel("z [$\AA{}$]")
        plt.ylabel("F [$k_BT$]")
        plt.title("{} res{} sep{}".format(molecule,resNr,nr))
        plt.savefig('dF-res{}-{}-{}'.format(resNr,molecule,nr),dpi=dpi)
        if (show_plots):
            plt.show()
    
        print "dF = ", mF
        print "dFr = ", mdFr
        print "dFl = ", mdFl
    
    # ------------------------------ # Write to file:
    
    def write_to_file(output, data, datafields, subheader):
        """ Write output file 
            len(headers) == len(datafields)
            headers = names of data fields
        """
        if (len(data) != len(datafields)):
            print "Error! number of data fields != number of headers!"
            print 'len:   ', len(data), len(datafields)
            print 'shape: ', np.shape(data), np.shape(datafields)
    
        ofile = open(output,'w')
        ofile.write("# Free energy profile of ion in nanopore\n")
        ofile.write(subheader)
        header = "# chunk "
        for element in datafields:
            header += element + " "
    
        header = header + '\n'
        ofile.write(header)
        
        it = 0
        for i in xrange(len(data[0])):
            line = str(it) + " "
            it += 1
            for j in xrange(len(data)):
                line += str(float(data[j][i])) + " "
            line += "\n"
            ofile.write(line)
                
        ofile.close()
        print "Finished writing file: ", output   
        
    output = "dF-res{}-{}-{}.dat".format(resNr,molecule,nr)
    write_to_file(output,[bins2,dF],['z','dF'],subheader="# min(dF)= {} zmid= {}\n".format(mF,zmid))
