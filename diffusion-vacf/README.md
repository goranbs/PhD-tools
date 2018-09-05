# diffusion_vacf.py

## python script for computing the diffusion coefficient from the velocity autocorrelation function (VACF)

```
python diffusion_vacf.py ifile.dat dt decorrelationtime nstep
```
where ifile.dat is the input file containing columnar data with velocity components: vx,vy,vz.
dt is the time step used in the simulation. The decorrelation time is the number of time steps
it takes before the velocities are decorrelated (the VACF are on average zero), this number can
in principle be low, it would just increase the sampling (number of VACF's being used in the average).
The nstep variable is the number of time steps to be used for each VACF. It is important that nstep is
large enough to capture the tail of the VACF so that the integral in this region becomes zero.

## About this script

This script computes the diffusion coefficient from the average velocity components v = (vx,vy,vz)

### useage
* useage: >> python diffusion_vacf_py ifile.dat dt decorrelationtime nstep
* exampl: >> python diffusion_vacf_py ifile.dat 2.0 100 1000

### output: 
* diffusion-ifile.out  : Running average of diffusion coefficients
* vacf-ifile.out       : The average of the VACF's

The relation between the diffusion coefficient D, and the velocity autocorrelation function (VACF) is

*    D = integral(VACF*dt)/d

where d is the dimensionality of the system (3D -> d=3) and the VACF is:

*    VCAF = v0*vi

where v0 is the initial velocity, and vi is the velocity at time step i.

