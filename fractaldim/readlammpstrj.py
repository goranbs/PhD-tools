#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:29:25 2016

@author: goranbs
------------------------------------------------------------------------------

Read LAMMPS trajectory: readlammpstrj.py

Python class for reading LAMMPS dump files.
The headers of the dump files is assumed to have a specific order, like this:
    
    ITEM: TIMESTEP
    X
    ITEM: NUMBER OF ATOMS
    N
    ITEM: BOX BOUNDS xy xz yz pp pp pp
    xlo xhi xy
    ylo yhi xz
    zlo zhi yz
    ITEM: ATOMS id type q x y z fx fy fz vx vy vz

where, X is the time step, N is the number of atoms.
The bounding box can be both triclinic (as shown in this example above), 
and orthorombic. The ATOMS section can have any kind of custom layout,
the ATOMS section parameters (here: id type q x y z fx fy fz vx vy vz),
will be identified and stored in a dictionary.
The values in the dictionary can be accessed at the different time steps.

Useage of readlammpstrj.py:

>> from readlammpstrj import LAMMPStrj as trj
>> obj = trj("lammpsfile.dump")
>> data = obj.get_data()
>> volume = obj.get_volume()
>> natoms = obj.natoms[-1]  # number of atoms at current time step
...
    
"""

# -------------------------------------------------------

import numpy as np # numpy library
import os          # path and os system handeling
import copy        # enable deepcopy

# -------------------------------------------------------

class LAMMPStrj(object):
    """
        Read LAMMPS trajectory file and store header info.
        Add functionality for reading timestep by timestep.
        Method for fetching data at the given timestep.
        
        methods:
        
            close_trj()         # close input file for reading
            get_data()          # read and return data in next time step section
            get_volume()        # return simulation box volume
        
        attributes:

            a                   # cell vector a
            b                   # cell vector b
            c                   # cell vector c
            alpha               # inclination angle alpha
            beta                # inclination angle beta
            gamma               # inclination angle gamma
            
            itime               # initial time step
            inatoms             # initial number of atoms
            isystem_size        # initial system size
    
            timesteps           # list of time steps timesteps[-1] is current time step
            natoms              # list of number of atoms natoms[-1] is current number
            system_size         # current system size
            
            data                # dictionary of data fields in ATOMS body section.

    """

    def __init__(self,LAMMPStrjfile):
        self.a = 0                              # cell vector a
        self.b = 0                              # cell vector b
        self.c = 0                              # cell vector c
        self.alpha = 0                          # inclination angle alpha
        self.beta = 0                           # inclination angle alpha
        self.gamma = 0                          # inclination angle alpha
        self._cosa = 0                          # cos(alpha)
        self._cosb = 0                          # cos(beta)
        self._cosg = 0                          # cos(gamma)
        self.trj = LAMMPStrjfile                # input file
        self._box = "orthorhombic"              # orthorhombic or trigonal
        self._errors = False                    # Errors encountered?
        self.inatoms = 0                        # number of atoms initially
        self.itime = 0                          # time of initial time step
        self.data = {}                          # data fields in main bodies
        self._header = 0                        # number of lines in headers
        self.timesteps = []                     # list of timesteps
        self.natoms = []                        # number of atoms
        self.nframes = 0                        # number of time frames in trajectory
        self.system_size = np.zeros((3,3))      # system size
        self.isystem_size = np.zeros((3,3))     # initial system size
        self._txt = self._extract_data()        # call to function for extracting header data
        self._nentries = len(self.data)         # number of body data fields
        
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
                if (len(x) == 2):
                    # Orthogonal box
                    #print "Orthogonal box"
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
                if (len(x) == 3):
                    # Trigonal box
                    #print " Trigonal box"
                    self._box = "trigonal"
                    try:
                        self.isystem_size[0,0] = float(x[0]) # xlo
                        self.isystem_size[0,1] = float(x[1]) # xhi
                        self.isystem_size[1,0] = float(y[0]) # ylo
                        self.isystem_size[1,1] = float(y[1]) # yhi
                        self.isystem_size[2,0] = float(z[0]) # zlo
                        self.isystem_size[2,1] = float(z[1]) # zhi
                        
                        self.isystem_size[0,2] = float(x[2]) # xy
                        self.isystem_size[1,2] = float(y[2]) # xz
                        self.isystem_size[2,2] = float(z[2]) # yz
                    except:
                        self._errors = True
                        print "Error: Could not convert system sizes to float:"
                        print x, y, z
                        break
                    break                
            if (it > 50):
                break
        try:
            self.system_size[0,0] = self.isystem_size[0,0]
            self.system_size[0,1] = self.isystem_size[0,1]
            self.system_size[1,0] = self.isystem_size[1,0]
            self.system_size[1,1] = self.isystem_size[1,1]
            self.system_size[2,0] = self.isystem_size[2,0]
            self.system_size[2,1] = self.isystem_size[2,1]
            self.system_size[0,2] = self.isystem_size[0,2]
            self.system_size[1,2] = self.isystem_size[1,2]
            self.system_size[2,2] = self.isystem_size[2,2]
        except:
            self._errors = True
            print "Error: Set set system size!"            
        
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
        self._set_simulationbox_parameters()
        self._volume = self._compute_volume()
        
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
        
    def _set_simulationbox_parameters(self):
        """ set cell vectors and inclination angles """
        
        xl = self.system_size[0,0]
        xh = self.system_size[0,1]
        yl = self.system_size[1,0]
        yh = self.system_size[1,1]
        zl = self.system_size[2,0]
        zh = self.system_size[2,1]
        
        xy = self.system_size[0,2]
        xz = self.system_size[1,2]
        yz = self.system_size[2,2]
        
        lx = (xh-xl)
        ly = (yh-yl)
        lz = (zh-zl)
                
        self.a = lx
        self.b = np.sqrt(ly**2 + xy**2)
        self.c = np.sqrt(lz**2 + xz**2 + yz**2)
        
        self._cosa = (xy*xz + ly*yz)/(self.b*self.c)
        self._cosb = xz/self.c
        self._cosg = xy/self.b
        
        self.alpha = np.arccos(self._cosa)
        self.beta = np.arccos(self._cosb)
        self.gamma = np.arccos(self._cosg)
        
    def _compute_volume(self):
        """ update the volume of the simulation box """
        self.volume = self.a*self.b*self.c*\
        (1 - self._cosa*self._cosa - self._cosb*self._cosb - self._cosg*self._cosg) + \
        2*(self._cosa*self._cosb*self._cosg)
        
    def get_volume(self):
        return self.volume
        
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
                if (self._box == "trigonal"):
                    for i in xrange(3):
                        xx = self._txt.readline()
                        it += 1
                        xx = xx.split()
                        self.system_size[i,0] = float(xx[0])
                        self.system_size[i,1] = float(xx[1])     
                        self.system_size[i,2] = float(xx[2])
                else:
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

        # 3) Update cell vectors, inclination angles and box volume:
        self._set_simulationbox_parameters()
        self._volume = self._compute_volume()
        
        return tmpdata
        
        


        
        