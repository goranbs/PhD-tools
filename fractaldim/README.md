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

1D_0000.dat, 2D_0000.dat and 3D_0000.dat, menger_sponge_0000.dat files used to produce the plot below. See dimensions.gnu. In dimensions.gnu, data from two huge Menger sponges are also included and which for clarity is not shown in the plot below.

![Plot of result from 1D, 2D, 3D and Menger sponge test systems](Dimensions.png?raw=true "Dimensions")

Image of a Menger sponge with theoretical fractal dimension Df = 2.726

![Plot of result from 1D, 2D and 3D test systems](Menger_sponge-3D-side-view.png?raw=true "Menger sponge")

The Menger sponge can be perodically replicated to produce a fractal dimension which is closer to the theoretical value. 
Below, an image of a 2x2x2 replica of the above Menger sponge is shown. This configuration, which consists of 2592000 particles, 
results in a fractal dimension Df = 2.72(9), where (9) is the standard deviation in the last digit. Out of curiosity, shifting the particle positions by a vector V, within a larger periodic simulation
box, we obtained  a fractal dimension Df = 2.71(7).

![Plot of result from 1D, 2D and 3D test systems](Menger_sponge-huge-3D-side-view.png?raw=true "Menger sponge")

The Menger sponge can be created from following a Moltemplate example found in the LAMMPS directory: lammps/tools/moltemplate/examples/misc_examples/menger_sponge/
Or see link: https://github.com/lammps/lammps/tree/master/tools/moltemplate/examples/misc_examples/menger_sponge
