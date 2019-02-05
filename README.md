# PhD-tools

Useful codes, scripts and tweaks 
written by Goran Brekke Svaland

##  Repository content:

* `readlammpstrj.py`
	* python class for reading LAMMPS trajectory, frame by frame
* `eldip.py`
	* python script for computing electric dipolemoment and orientational order parameter Tp, using readlammpstrj.py
* `integrateProfile.py`
	* python script for integrating LAMMPS profiles generated from ave/time or ave/chunk
* `heatmap.py`
	* python script for making 2D heatmaps from scattered data
* percolation
	* `run_voro++_cl.sh` is a script to modify the input .gro file and run voro++ to find neighbors of each atom
	* `percolation_cl.py` is python script that takes the output of voro++ to cluster atoms based on whether their voronoi cells share a face, and dumps a vmd file that can be sourced by vmd to create atom representations of the clusters found by the code.
	* `run_cl_traj.sh` is a script to run the percolation code for a trajectory.
	* scripting: sub directory containing some random overlay scripts used to run voro++ and the percolation code.
* BASH-scripts
	* `displacemolecules` displace two molecules using VMD
	* `process_lammpslog` process LAMMPS log files
	* `GROMACSconfmerge` merge two gromacs configurations into one
	* `gnuplotting` use BASH to plot data via gnuplot
* savitzky-golan-filter
	* `smoothing.py` filter noisy data using a Savitzky-Golan filter
* avgdataframes
	* `avgdataframes.py` using pandas Data Frames to average over several data frames.
* 2Drdf
	* `2Drdf.py` computes the two dimensional radial distribution function in the xy-plane
* ovitoscripting
	* `clustercompute.py` utilizes the DXA algorithm in ovito for dislocation analysis
* diffusion-vacf
	* `diffusion_vacf.py` compute diffusion coefficient from average velocities. Based on the velocity autocorrelation function.
* fepstuff
	* `fepstuff.py` compute free energy profile of ion in confinement
* piechart
	* example script for plotting a pie chart
* solvate
	* solvate two mineral slabs with water using GROMACS
* fractaldim
	* `fractaldim.py` compute the fractal dimension of a system using the box counting algorithm
* hist-gauss
	* compute histogram of columnar dataset and make a fit to a gaussian bell curve
* surface-roughness
	* scripts for computing colloidal interaction potential
