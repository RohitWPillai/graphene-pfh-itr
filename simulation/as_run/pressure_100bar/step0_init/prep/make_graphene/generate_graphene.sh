#!/bin/bash

# Input parameters
N=${1:-7}     # Number of sheets - should be odd a number
xL=${2:--12}   # Length in x-direction (lattice units)
xH=${3:-12}    # Length in x-direction (lattice units)
yL=${4:--6}    # Length in y-direction (lattice units)
yH=${5:-6}     # Length in y-direction (lattice units)
sSep=3.35      # Spacing between sheets (in Angstroms)
epsCC=0.0046   # Intralayer interaction between graphene sheets
sigCC=3.276
rcut=15
mCarb=12.0170 # g / mol
seedrndm=12423
temp=300.0

# Output LAMMPS input file
output="in.run"     

# Header
cat << EOF > $output
#############################################################################################
# LAMMPS input script to equilibrate a graphite-water system                                #
# Author : Thejas Hulikal Chakrapani                                                        #
# Institute for Multiscale Thermofluids                                                     #
# University of Edinburgh                                                                   # 
# Date : 20th February, 2025                                                                #    
#                                                                                           #
# References :                                                                              #
# [1] Units Conversion  :   https://physics.nist.gov/cuu/Constants/energy.html              #
#                           https://physics.nist.gov/cuu/Constants/factorlist.html          #
# [2] C-C interaction   :  Lindsay, L., and D. A. Broido.                                   #
#                          “Optimized Tersoff and Brenner Empirical Potential Parameters    #
#                          for Lattice Dynamics and Phonon Thermal Transport in Carbon      #
#                          Nanotubes and Graphene.”                                         #
#                          Physical Review B 81, no. 20 (May 27, 2010): 205441.             #
#                          https://doi.org/10.1103/PhysRevB.81.205441.                      #
#############################################################################################

#------------------
# General commands 
#------------------

log             log.graphene_generate
dimension       3            
units           metal        
boundary        p p p
atom_style      full
newton          on

#--------------------------------------------------------------
# Variables for converting real to metal units
#--------------------------------------------------------------
# Please refer to NIST [1] website for (basic) units conversion
# -------------------------------------------------------------

variable    mCarb           equal 12.0170 # g / mol

variable bond equal 1.42
variable latC equal 2.46                # Lattice constant in the basal plane graphite

variable third     equal   1./3.
variable twoThird  equal   2./3.
variable half      equal   1./2.
variable sixth     equal   1./6.

variable zL equal -${N}/2.
variable zH equal  ${N}/2.

variable a equal sqrt(3.0)/2.0          # Primitive lattice vectors in graphite
variable b equal 1./2.
variable c equal -sqrt(3.0)/2.0
variable d equal sqrt(3.0)*\${latC}

lattice custom  1                                               &
                a1      \${latC}      0               0          &
                a2      0            \${d}            0          &
                a3      0            0               ${sSep}    &
                basis   0            0               0          &
                basis   \${half}     \${sixth}         0          &
                basis   \${half}     \${half}          0          &
                basis   0           \${twoThird}      0

region simbox block ${xL} ${xH} ${yL} ${yH} \${zL} \${zH} units lattice
create_box $((N)) simbox  # Allocate types for all sheets

EOF

# Initial z-coordinate range
zLower=$(echo "(1./2.)*(${N}-1)-0.1" | bc -l)
zUpper=$(echo "(1./2.)*(${N}-1)+0.1" | bc -l)

# Displace sheets for stacking
for i in $(seq 1 $N); do
    zDisplacement=$(echo "-(${N}-${i})" | bc -l)
    if (( i % 2 == 0 )); then
        # B layer: displace in both z and y directions
        cat << EOF >> $output
region sheet$i block ${xL} ${xH} ${yL} ${yH} ${zLower} ${zUpper} units lattice
create_atoms $i region sheet$i
group graf$i region sheet$i
displace_atoms graf$i move 0 \${third} ${zDisplacement} units lattice

EOF
    else
        # A layer: displace only in z direction
        cat << EOF >> $output
region sheet$i block ${xL} ${xH} ${yL} ${yH} ${zLower} ${zUpper} units lattice
create_atoms $i region sheet$i
group graf$i region sheet$i
displace_atoms graf$i move 0 0 ${zDisplacement} units lattice

EOF
    fi
done

# Define hybrid force field and pair_coeffs
cat << EOF >> $output

EOF

echo "      " >> $output

# Footer
cat << EOF >> $output

# Set Masses 
mass *     ${mCarb}

# Set Charges 
set type *     charge   0.0000 # C

# Set Velocities 
velocity all  create ${temp} ${seedrndm} dist gaussian  mom yes

run 1

write_data out.graphene types numeric pair ij

EOF

echo "LAMMPS input file '$output' created successfully."
