#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:29:25 2016

@author: goranbs
--------------------------------------

Read LAMMPS trajectory: readlammpstrj.py

Description:

    Read LAMMPS trajectory and perform systematic analysis on each snapshot
    
Usage:

    >> readlammpstrj.py input.lammpstrj method keywords
    
-------------------------------------

"""

import numpy as np # numpy library
import os          # path and os system handeling
import copy

#-- Idea: Create script that uses readlammpstrj or perhaps even readlammpsdata
#         and performs some method on the read data.
#import argparse    # handle command line arguments
#wa_parser = argparse.ArgumentParser(description='Read LAMMPS trajectory')
#wa_parser.add_argument('-w','--waterangle', metavar='XYZ', type=str, nargs=1, 
#                    help='Compute angle of water molecules in X,Y or Z direction')                    
#wa_parser.add_argument('-m', '--minmaxpos', metavar='R', type=float, nargs=2,
#                           help='Define min and max position in system from where method should be applied. R=float')
#args = wa_parser.parse_args()
#print args.minmaxpos
#print args.waterangle[0]
#print type(args.waterangle[0])

# -------------------------------------------------------

class LAMMPStrj(object):
    """
        Read LAMMPS trajectory file and store header info.
        Add functionality for reading timestep by timestep.
        Method for fetching data at the given timestep.
        
        properties:
            itime = initial time step
            atoms = total number of atoms
            system_size = boundary conditions
            data = data values in colums of trajectory
            
        attributes:
            itime
                - return initial time step
            get_time()
                - return time step of snapshot
            natoms
                - return number of atoms
            system_size
                - return simulation box boundaries
            data
                - return dictionary of data fields
        
    """

    def __init__(self,LAMMPStrjfile):
        self.trj = LAMMPStrjfile
        self._errors = False  # Errors encountered?
        self.inatoms = 0      # number of atoms initially
        self.itime = 0        # time of initial time step
        self.data = {}        # data fields in main bodies
        self._header = 0      # number of lines in headers
        self.timesteps = []   # timesteps
        self.natoms = []      # number of atoms
        self.nframes = 0      # number of time frames in trajectory
        self.system_size = np.zeros((3,2))  # system size
        self.isystem_size = np.zeros((3,2)) # initial system size
        self._txt = self._extract_data()    # call to function for extracting header data
        self._nentries = len(self.data)     # number of body data fields
        
    def __str__(self):
        " Info that is printed if class name is printed."
        return "Class {}".format(self.__class__.__name__)
    
    def __repr__(self):
        my_repr = "{}(input={})".format(self.__class__.__name__,\
                                        " 'LAMMPS-trajectory-file' ")
        return my_repr
        
    def _seekfor(self, ofile, this="txt"):
        it = 0
        value = 0
        while True:
            it += 1
            line = ofile.readline()
            line = line.split()
            line = " ".join(line[1:])
            if (line == this):
                val = ofile.readline()
                value = int(val)
                break
            # break if more than 50 lines read or line = "ITEM: ATOMS"
            if (line == "ATOMS" or it > 50):
                self._errors = True # error occured
                print "Error: Did not find section: ", this, " in header"
                break
        ofile.seek(0,0) # go back to top of file
        return value
        
    def _seekfor_system_size(self, ofile):
        it = 0
        while True:
            it += 1
            line = ofile.readline()
            line = line.split()
            try:
                test = line[2]
            except:
                test = None
            if (test == "BOUNDS"):
                x = ofile.readline()
                y = ofile.readline()
                z = ofile.readline()
                x = x.split()
                y = y.split()
                z = z.split()
                try:
                    self.isystem_size[0,0] = float(x[0])
                    self.isystem_size[0,1] = float(x[1])
                    self.isystem_size[1,0] = float(y[0])
                    self.isystem_size[1,1] = float(y[1])
                    self.isystem_size[2,0] = float(z[0])
                    self.isystem_size[2,1] = float(z[1])
                except:
                    self._errors = True
                    print "Error: Could not convert system sizes to float:"
                    print x, y, z
                    break
                break
            if (it > 50):
                break
        ofile.seek(0,0)
        
    def _seekfor_data(self,ofile):
        it = 0
        self._header = 0
        while True:
            it += 1
            line = ofile.readline()
            line = line.split()
            try:
                test = line[1]
            except:
                test = line[0]
            
            if (test == "ATOMS"):
                # create dictionary with keys in first header:
                for i in range(len(line[2:])):
                    j = i+2
                    self.data.update( {line[j] : [] } )
                self._header = it
                break
            if (it > 50):
                self._noerrors = False
                print "Error: Could not find ITEM: ATOMS section. Final line read:"
                print line
                break
        ofile.seek(0,0)
        
    def _open_trj(self, ofile):
        try:
            txt = open(ofile,'r')
        except:
            print "\n I/O - Error. Could not read file: ", ofile, '\n'
            raise IOError
            
        basename, file_extension = os.path.splitext(ofile)

        if file_extension != ".lammpstrj":
            print "I/O - Error:\nFile suffix lammpstrj !=", file_extension
            raise IOError
        
        return txt
        
    def close_trj(self):
        """ Close file """
        self._txt.close()
        
    def _extract_data(self):
        """ Extract data from header section """
        
        txt = self._open_trj(self.trj)
    
        self.itime = self._seekfor(txt, "TIMESTEP")
        self.inatoms = self._seekfor(txt,'NUMBER OF ATOMS')
        self._seekfor_system_size(txt)
        self._seekfor_data(txt)
        
        if (self._errors):
            print "Error occured!"
            txt.close()
        
        ntimeframes = 0
        # we have to read through the file in order to find the number of time frames:
        for line in txt:
            if (line[0] == "I"):
                line = line.split()
                line = " ".join(line[1:])
                if (line == "TIMESTEP"):
                    ntimeframes += 1
        
        self.nframes = ntimeframes
        txt.seek(0,0)
        return txt
        
    def get_data(self):
        """
            Read next data block from trajectory file.
            Return data according to data fields in main body
            - what to do when reached EOF?
        """
        
        # create new temporary dictionary with self.data as template.
        tmpdata = copy.deepcopy(self.data)
        
        # 1) Read header of section
        key = "no"     # default key
        it = 0         # reset iterator
        entries = []   # reset entries
        while (it < self._header):
            line = self._txt.readline()
            it += 1
            line = line.split()
            try:
                key = str(line[1])
            except:
                key = "no"
                
            if ( key == "TIMESTEP"):
                time = self._txt.readline()
                it += 1
                self.timesteps.append(int(time))
            if (key == "NUMBER"):
                natoms = self._txt.readline()
                it += 1
                Natoms = int(natoms)
                self.natoms.append(Natoms)
            if (key == "BOX"):
                for i in xrange(3):
                    xx = self._txt.readline()
                    it += 1
                    xx = xx.split()
                    self.system_size[i,0] = float(xx[0])
                    self.system_size[i,1] = float(xx[1])    
            if ( key == "ATOMS"):
                entries = line[2:]

        # 2) Read body of section
        for i in range(Natoms):
            line = self._txt.readline()
            line = line.split()
            for i in xrange(self._nentries):
                entry = entries[i]
                tmpdata[entry].append(float(line[i]))

        
        return tmpdata
        
        


        
        