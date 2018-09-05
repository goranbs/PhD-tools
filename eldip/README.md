# eldip.py

`eldip.py`

Compute electric dipole moment of system, or individual atoms in the system. Has been used to compute the electric dipole moment (EDP) of SPC/E water and binning the EDPs in the system.

* utilizes readlammpstrj.py

* lammpstrajectory should contain information about:
1. x,y,z
2. q
3. id
4. mol
3. fx,fy,fz

* computes the electric dipolemoment of atoms within the same molecule (only for three atoms at a time)
* 1D and 2D binning
* computes orientational order parameter Tp
* computes angle between dipole and provided normal vector
* computes total dipolemoment of system given atom types
* computes total forces on atoms given types
* computes volume, cell vectors and inclination angles

## Usage:
>> eldip -h
>> eldip -i file.lammpstrj -t [1,2,3] -x 10 20 -d -s -f -p --dump1d z lower 0.1

##Install:
* requires Python 2.7
* requires readlammpstrj.py
* requires python numpy
* requires python argparse
* requires python sys

You can install through terminal in the following manner:
>> sudo ln -s ~/path/to/eldip.py /usr/local/bin/eldip

Or you can create an alias in the .bashrc file:
>> alias 2Drdf='~/mypath/to/2Drdf.py' 

