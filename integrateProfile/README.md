integrateProfile.py

integrateProfile.py is a python script that integrates LAMMPS .profile files.

Integrations can be done over several .profile files at a time, summing up the total integral in an additional column.
Columns can be weighted using the -q flag. This is intended to be used for number densities of particles, so that the electric field from these particles can be found.
Two integrations over number densities gives you the electric potential.

Example of usage: \n

integrateProfile.py -h \n
integrateProfile.py -i Na.profile Cl.profile -q [1,-1] -c [2,5] -o efield-NaCl.profile \n
integrateProfile.py -i efield-NaCl.profile -q [1] -c [1,3] -o potential-NaCl.profile \n

Restriction: Assumes one chunk, so averages in LAMMPS should be done over one simulation run. E.g: \n

variable Ne equal 100                     # Nevery; use every Ne in averages \n
variable Ns equal 100000                  # Nstep; number of time steps \n
variable Nr equal ${Ns}/${Ne}             # Number of values in average \n

compute DensityOfWater spce chunk/atom bin/1d z lower ${resolution} units box region CenterReg \n
fix DensityProfileWater spce ave/chunk ${Ne} ${Nr} ${Ns} Dwater density/mass density/number norm all file spce-${name}-${runNr}.profile \n


