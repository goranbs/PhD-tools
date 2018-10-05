#############################################################################################
# system of LJ particles
#############################################################################################

# ---------------   System information   --------------------
units		lj
atom_style	atomic
boundary	p p p
neighbor	0.3 bin
neigh_modify	every 1 delay 1 check yes one 3000

# -----------------------------------------------------------
log		log.lammps-eq

variable	rho equal 0.3 			# density
variable	T equal 0.75			# temperature

# -----------------------------------------------------------

variable 	S equal 10			# system size
variable	V equal ${S}*${S}*${S}		# volume
variable	N equal round(${rho}*${V})	# number of particles

region		SimBox block 0 ${S} 0 ${S} 0 ${S} units box
create_box	1 SimBox

fix		1 all deposit ${N} 1 1 18127 attempt 1000 region SimBox near 1.0
group		liquid type 1

variable	cutoff equal 2.5
pair_style	lj/cut ${cutoff}
pair_coeff	1 1 1 1 ${cutoff} 
mass		1 1 

variable	dt equal 0.0005
variable	tc equal 200*${dt}
variable	tp equal 1000*${dt}
variable	Ns equal 100000
variable	Np equal ${Ns}/1000

timestep	${dt}    # default lj units
thermo		${Np}
thermo_style	custom step temp press pe ke density

#dump		Dump all atom ${Np} lj.lammpstrj

# -- run enough timesteps to insert all atoms via fix deposit
run		${N}

# -- then run equilibration in NVT

velocity         all create 1 10753 dist gaussian

#fix		NPT all npt temp ${T} ${T} ${tc} iso ${P} ${P} ${tp}
fix		NVT all nvt temp ${T} ${T} ${tc}

run		${Ns}

write_data	eq.data nocoeff

