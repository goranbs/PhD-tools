#!/usr/bin/python
"""
Created on Mon Aug 22 14:04:11 2016

@author: goranbs

unit-testing of 'readlammpstrj.py'
----------------------------------
execute test folder:

>>nosetests test_readlammpstrj.py
"""

global AFILE
AFILE = "freerun-3.lammpstrj"

import sys
sys.path.insert(0, '../')

from readlammpstrj import LAMMPStrj

def test_numberofatoms():
    """ Find nr of atoms.
        Testing of initialization of lammps data file. 
    """
    filename = AFILE            # test file
    obj = LAMMPStrj(filename)
    
    Natoms = obj.inatoms

    one = (Natoms == 24474)
    
    assert one
    
def test_timestep():
    """ Find initial timestep """
    filename = AFILE            # test file
    obj = LAMMPStrj(filename)
    
    Itime = obj.itime

    itime = (Itime == 0)
    
    obj.close_trj()
    assert itime

    
def test_seekfor_system_size():
    """ Test if readlammps manage to read system size correctly. """

    filename = AFILE # testfile
    obj = LAMMPStrj(filename)
    
    # system size of testfile
    xloxhi = (-0.134995, 79.865)
    yloyhi = (-0.043499, 39.7835)
    zlozhi = (0.9, 180.9)
    
    xl = obj.isystem_size[0,0] == xloxhi[0]
    xh = obj.isystem_size[0,1] == xloxhi[1]
    yl = obj.isystem_size[1,0] == yloyhi[0]
    yh = obj.isystem_size[1,1] == yloyhi[1]
    zl = obj.isystem_size[2,0] == zlozhi[0]
    zh = obj.isystem_size[2,1] == zlozhi[1]
    
    testlist = []
    testlist.append(xl), testlist.append(xh), testlist.append(yl)
    testlist.append(yh), testlist.append(zl), testlist.append(zh)
    
    ok = True
    if False in testlist:
        ok = False

    obj.close_trj()
    assert ok
    
def test_seekfor():
    """ Test if _seekfor(ofile,this) finds correct locations in file"""
    
    filename = AFILE # testfile
    obj = LAMMPStrj(filename)
    
    txt = obj._open_trj(obj.trj)

    itime = obj._seekfor(txt, "TIMESTEP")
    inatoms = obj._seekfor(txt,'NUMBER OF ATOMS')
    
    ok=True
    if (itime != 0 or inatoms != 24474):
        ok=False
        
    obj.close_trj()
    assert ok
    
    
def test_get_data():
    """ Test if get_data() reads next data blocks in trajectory"""
    
    # this test checks only the two first data blocks...
    
    filename = AFILE # testfile
    obj = LAMMPStrj(filename)

    data1 = obj.get_data()
    data2 = obj.get_data()
    
    keys1 = data1.keys()
    keys2 = data2.keys()

    # does data2 contain the same keys as data1?
    keysIn_data2 = True
    for key in keys1:
        haskey = data2.has_key(key)
        if haskey == False:
            keysIn_data2 = False  # then key not in data2
            print "not in data2: ", key
    
    keysIn_data1 = True
    for key in keys2:
        haskey = data1.has_key(key)
        if haskey == False:
            keysIn_data2 = False  # then key not in data1
            print "not in data1: ", key

    ok=False                
    if (keysIn_data1 and keysIn_data2):
        ok=True
        
    obj.close_trj()
    
    assert ok
    

def test_readsAllAtoms():
    """ Test if the get_data() function reads as many lines as there are atoms
        in the system. If the number of atoms change, this should be taken
        into account.
    """
    # !!! Assumption: all atoms has got an atom id, and that this id is named "id"
    ID="id"
    # Test only the two first data sections:
    
    filename = AFILE # testfile
    obj = LAMMPStrj(filename)
    

    data1 = obj.get_data()
    data2 = obj.get_data()
    
    natoms = obj.natoms     # list of #atoms at time1 and time2
    
    ids1 = data1[ID]
    ids2 = data2[ID]

    idpresent1 = natoms[0] in ids1
    idpresent2 = natoms[1] in ids2
    
    ok=False
    if idpresent1 and idpresent2:
        ok=True

    obj.close_trj()    
    
    assert ok