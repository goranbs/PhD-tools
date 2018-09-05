# fepstuff

Estimate the free energy barrier of a single ion moving from one side
of the confinement to the other.

beta*dF = -ln(p(z))

where p(z) is the probability density of finding an ion at a specific location

p(z) is found from histogramming the positions of the ion throughout the simulation.


a) How to find the center of the confinement?
    1) Identify the four tallest peaks in the density profile
    2) The center point between the two central peaks is the center of the confinement.

b) Compute histogram of ion positions and normalize it

c) Look at the symmetry of p(z) -> plots

d) compute total p(z) reflecting the particle positions around the center of the confinement

f) write data to file

g) BASH/gnuplotting/*makeplots-dF.sh are shell scripts for plotting output data

