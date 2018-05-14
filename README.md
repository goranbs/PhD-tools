# PhD-tools
Useful scripts and tweaks

1) displacemolecules:
Displace and align two vmd readable molecules in z-direction 
and apply some angle of rotation in the xy-plane for the
upper molecule.

2) solvate:
Solvate molecules with SPC water.
Possibility to create vacuum in some direction.
Output LAMMPS data files and pdb files.

3) readlammpstrj:
	-readlammpstrj.py
		python class for reading LAMMPS trajectory, frame by frame
	-eldip.py
		python script for computing electric dipolemoment and Tp, an
		orientational order parameter:

		Tp = 1/2<(3u^2 - 1)>, Tp E [-0.5,1.0]

		where u is the local dipole moment of some molecules.
		Tp indicates the tendency of the dipolemoment being aligned
		with some given vector n, that is provided through the 
		commandline.

		1D and 2D binning in x,y,z dimentions available
	
		display usage information:
		>> python eldip.py -h

4) hist-gauss
	- compute histogram of input data set and make a fit to a gaussian

