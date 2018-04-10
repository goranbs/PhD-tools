#############################################################################################
#
# filename: "in.run"
#
#############################################################################################

# ---------------   System information   --------------------
units            real
atom_style       full
boundary         p p p
neighbor         0.3 bin
neigh_modify     every 1 delay 1 check yes one 3000
kspace_style     pppm 1.0e-4
variable         cutoff equal 12.5
pair_style       hybrid/overlay lj/cut ${cutoff} coul/long ${cutoff}
bond_style       harmonic
angle_style      harmonic

# ---------------------------  Gold and  spce  ----------------------------------------------

read_data "gold-773-8nm.data"
include "ff-Au-spce-Cl.in.settings"       # Force field
balance 1.0 shift z 10 1.0                  # balance cpus

# -------------------------------------------------------------------------------------------
variable name string 773-8nm
variable runNr equal 2
log log.lammps-${name}-${runNr}
# -------------------------------------------------------------------------------------------

# -- set execution variables:
variable Ne equal 100                     # Dump log info every
variable Nprint equal 1000                # Dump trajectory every
variable Ns equal 30000                    # Number of time steps

variable Nr equal ${Ns}/${Ne}             # Number of values in average

variable T equal 300                      # Temperature [K]
variable dt equal 2.0                     # time step [fsec]
variable tc equal 100*${dt}               # Temprature damping factor
variable Ang equal 2                      # tot Angstrom to move upper slab during sim

# ------- Define regions, groups and calculations -------------------------------------------
variable         zmid equal (zhi-zlo)/2
region           lower block EDGE EDGE EDGE EDGE EDGE ${zmid} side in  # lower region
region           upper block EDGE EDGE EDGE EDGE ${zmid} EDGE side in  # upper region
group            Lower region lower                                    # all particles lower
group            Upper region upper                                    # all particles upper
group            AuL intersect Lower Au                                # lower slab
group            AuU intersect Upper Au                                # upper slab
compute          Lcom AuL com                                          # center of mass lower
compute          Ucom AuU com                                          # center of mass upper
#group            sub id 127 306 559 738                                # four Au atoms
#fix              freeze sub move linear 0.0 0.0 0.0                    # constrain sub group
#group            otherAu subtract Au sub                               # all other Au, not in sub


##-- Vertical density profile
variable resolution equal 0.1   # resolution of bins in z-dir

compute Dwater spce chunk/atom bin/1d z lower ${resolution} units box
compute DAu Au chunk/atom bin/1d z lower ${resolution} units box
compute Dow Ow chunk/atom bin/1d z lower ${resolution} units box
compute Dhw Hw chunk/atom bin/1d z lower ${resolution} units box

fix DproW spce ave/chunk ${Ne} ${Nr} ${Ns} Dwater density/mass density/number norm all file spce-${name}-${runNr}.profile
fix DproAu Au ave/chunk ${Ne} ${Nr} ${Ns} DAu density/mass density/number norm all file Au-${name}-${runNr}.profile
fix DproWO Ow ave/chunk ${Ne} ${Nr} ${Ns} Dow density/mass density/number norm all file Ow-${name}-${runNr}.profile
fix DproWH Hw ave/chunk ${Ne} ${Nr} ${Ns} Dhw density/mass density/number norm all file Hw-${name}-${runNr}.profile

# ------- Thermodynamic output --------------------------------------------------------------

timestep         ${dt}
thermo           ${Ne}


# ------- Thermostatting -----------------------------------------------------------------

dump  Dump all custom ${Nr} ${name}-${runNr}.lammpstrj id mol type q x y z fx fy fz vx vy vz
velocity         all create ${T} 11753 dist gaussian rot yes

##################################################################################
# -- For displacing the upper slab:
#variable         vel equal ${Ang}/(${Ns}*${dt})  # velocity [Angstrom/fsec]
#fix              MV AuU move linear 0 ${vel} 0 units box
##################################################################################

#fix              NVTau otherAu nvt temp ${T} ${T} ${tc}
fix              NVTw spce nvt temp ${T} ${T} ${tc}

run              ${Ns}

# -----------------------------------------------------------------------------------------

write_data       ${name}-${runNr}.data

# ----------------------- EOF -------------------------------------------------------------
