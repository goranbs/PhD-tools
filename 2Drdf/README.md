# 2Drdf.py

`2Drdf.py`

Computes the two dimensional radial distribution function in the xy-plane (disks).

## Usage:
>> 2Drdf -h
>> 2Drdf -i file.lammpstrj -p [1,2,3] -a [1,2,3] -b [3,4,5,6] -d 0.2 -o output_file_name -e 10 -n 100 --nbins 100 --rmax 20 --plot

## Install:

* requires Python 2.7
* requires readlammpstrj.py
* requires python numpy
* requires python argparse
* requires python sys
* matplotlib

You can create a symbolic link to the program in the following manner:
>> sudo ln -s ~/path/to/2Drdf.py /usr/local/bin/2Drdf

Or you can create an alias in the .bashrc file:
>> alias 2Drdf='~/mypath/to/2Drdf.py' 

## TODO:

* How to account for periodic boundaries in the g(r)?
	* We could replicate the system in the loop. Looping over relative positions in the xy-plane, just adding the system box size in the x and y directions.

* How can we make the script more general? 
	* Define a plane by a surface normal, and compute the rdf in that plane.
	* The challenge is to compute the reference volume of the slice. This slice is the intersection between the surface and the simulation box. However, this is just a problem of correct normalization.
