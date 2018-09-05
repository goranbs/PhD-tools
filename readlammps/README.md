# readlammps.py

## python class for reading LAMMPS dump files, frame by frame

* returns a dictionary at a timestep with the entries in the trajectory
* access number of atoms at every time step
* access system volume, cell vectors and inclination angles at every time step
* access system boundaries at every time step

Import readlammpstrj in another python script (see eldip.py):
```
from readlammpstrj import LAMMPStrj as trj
```

Calling readlammpstrj with an input file:

```
obj = trj(args.inputfile[0]) # create LAMMPStrj object
```

