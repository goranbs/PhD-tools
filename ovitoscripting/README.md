# clustercompute.py

`clustercompute.py`

Utilizes Ovito's scripting language and its DXA algorithm for computing dislocations.

Setting:
The point of the code is to characterize a calcium carbonate crystal system which contains
a high-Mg-concentration cluster (HMC). The HMC has been created by converting a fraction of
the Ca ions to Mg in a sphere of radius r0 with center [x0,y0,z0]. The HMC has in some cases
been rotated an angle theta inside the calcite crystal. In this case, a disrupted layer is
created around the HMC. We would like to characterize this layer and the inner HMC sphere.
Using the DXA algorithm we characterize what particles belong to a fcc-lattice, the ones that
does not is shown to belong to the disrupted layer. We can thus compute the volume of the sphere
containing the disrupted layer and the sphere within containing particles still belonging to a
fcc-lattice symmetry.

## Usage:
>> ./ovitos clustercompute.py ifile.lammpstrj ofile.dat r0 dr [x0,y0,z0] [a,b,c,...]
where
* r0 = radius of the initial sphere
* dr = approximal thickness of the amporphous layer around the initial sphere
* [x0,y0,z0] = center point of the initial sphere
* [a,b,c] = list of particles to exclude from the calculation

## Install:

* requires Ovito
* requires ovitos
* requires python 3.4
* requires python numpy
* requires python sys

You can create a symbolic link to the ovitos executable in the following manner:
>> sudo ln -s ~/path/to/ovito/ovitos /usr/local/bin/ovitos

Whereby you will be able to run the command shown under useage as for example:
>> ovitos clustercompute.py ifile.lammpstrj ofile.dat 19 2.0 [10,10,10] [3,4]


