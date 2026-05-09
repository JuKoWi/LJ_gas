source $env(VMDDIR)/plugins/noarch/tcl/pbctools2.8/pbctools.tcl
package require pbctools
# Load the trajectory file
mol new traj.lammpstrj type lammpstrj waitfor all

# Set drawing style
mol delrep 0 top                  ;# delete default representation
mol addrep top                    ;# add new one
mol modstyle 0 top VDW 0.1 12.0  ;# VDW spheres, radius 0.3, resolution 12

# Optional: color by atom index
#mol modcolor 0 top Index

# Set background to white
#color Display Background white

pbc box -center origin


# Play the trajectory
animate goto 0        ;# rewind to start
mol top [molinfo top] ;# make sure it's the active mol
animate forward       ;# hit play
