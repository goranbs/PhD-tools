## 3D Fractal calculation

fractaldim.py

	computes the number of boxes N, occupied by particles, based on their coordinates.

Usage:

	fractaldim.py -h
    
Output:

	prefix_XXXX.dat 
	
	output data columns:
	[chunkNr,rx,ry,rz,N,log(1/rx),log(1/rx),log(1/rx),log(N)]

	The output file is written for each frame in the trajectory, and can be used to compute the fractal dimension of the system: log(N) ~ Df*log(1/r), where Df is the fractal dimension of the system.
        
	If the coordinates of the particles in the LAMMPS trajectory is in fractional coordinates [xs,ys,zs], fractional coordinates will be used. 
	If the coords are unscaled [x,y,z], unscaled coordinates will be used. 
	Unwrapped coordinates [xsu, ysu, zsu], are not supported.

Dependencies:
* readlammpstrj.py to read particle positions
* numpy
* argparse
* sys

Files:
* fractaldim.py	: 3D fractal calculation from LAMMPS trajectory
* readlammpstrj.py : read LAMMPS trajectories
* in.generate_system : LAMMPS input file to generate test systems
	* fcc-1D,2D,3D.lammpstrj : test systems
	* 1D,2D,3D.dat : results generated with fractaldim.py
* dimensions.gnu : plot results from test systems


### Plot of log(N) vs. log(1/r) from the generated test systems

1D_0000.dat, 2D_0000.dat and 3D_0000.dat files used to produce the plot below. See dimensions.gnu.

![Plot of result from 1D, 2D and 3D test systems](Dimensions.png?raw=true "Dimensions")

