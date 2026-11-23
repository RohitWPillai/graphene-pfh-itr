#!/bin/bash
# graphene_PFH_npt_workflow.sh
# This script automates the workflow for preparing graphene, generating a PFH slab,
# converting the packmol output to a LAMMPS-readable format, and setting up an NPT simulation.
#
# Folder structure assumptions:
#   runFiles/                  
#      SiC.tersoff
#      submit
#      graphene/generate_graphene.sh
#      h2o/calc.h2o_slab_len.py
#      h2o/h2o_gen_lmp_in.cpp
#      h2o/in.h2o_packmol
#      h2o/in.h2o_tip4p
#      in.forcefield        <-- original force field template
#
#   prep/                      --> will be created if not present
#       make_graphene/         --> will contain generate_graphene.sh and graphene output files
#       make_PFH/            --> will contain PFH-related files (copied from runFiles)
#       combine/               --> will contain the combined system input file and copied in.graphene/in.PFH/in.forcefield
#
#   npt/                       --> should exist and contain necessary simulation files (e.g. in.varDefs)
#
# User-defined inputs:
SHEETS=9
slab_width=60      # PFH slab thickness in Å (this controls pfhboxlen)
sheet_sep=3.35
LMP_MPI="/work/<project>/<project>/<user>/lammps-stable_29Aug2024/src/lmp_mpi"
molTemp="/work/<project>/<project>/<user>/moltemplate/moltemplate/scripts/moltemplate.sh"
cleanup="/work/<project>/<project>/<user>/moltemplate/moltemplate/scripts/cleanup_moltemplate.sh"
# LMP_MPI="/path/to/lammps-stable_29Aug2024/src/lmp_mpi"

# User-defined NPT variables:
temp_val=300
Ndump_val="1e4"
Nrst_val="2e6"
Ninpt_val="2e4"
Ninpt_long_val="10e6"
Nevery_val="1"
Nrep_val="5e4"
Nfreq_val="5e4"
Nt_val="1e4"

home=$(pwd)

# (Density is prescribed in calc.h2o_slab_len.py)

set -e  # Exit immediately if any command fails.

# Function to print an error message and exit.
error_exit() {
    echo "Error: $1"
    exit 1
}

#######################################
# Directory Setup: prep and make_graphene
#######################################
if [ ! -d "prep" ]; then
    echo "Directory 'prep' not found. Creating 'prep' and 'prep/make_graphene'..."
    mkdir -p prep/make_graphene || error_exit "Failed to create prep/make_graphene directory."
    cp runFiles/graphene/generate_graphene.sh prep/make_graphene/ || error_exit "Failed to copy generate_graphene.sh to prep/make_graphene."
else
    if [ ! -d "prep/make_graphene" ]; then
        echo "Directory 'prep/make_graphene' not found. Creating it and copying generate_graphene.sh..."
        mkdir -p prep/make_graphene || error_exit "Failed to create prep/make_graphene directory."
        cp runFiles/graphene/generate_graphene.sh prep/make_graphene/ || error_exit "Failed to copy generate_graphene.sh to prep/make_graphene."
    fi
fi

#######################################
# Directory Setup: prep/make_PFH (PFH files)
#######################################
if [ ! -d "prep/make_PFH" ]; then
    echo "Directory 'prep/make_PFH' not found. Creating it..."
    mkdir -p prep/make_PFH || error_exit "Failed to create prep/make_PFH directory."
    cp -r runFiles/PFH/*  prep/make_PFH/ || error_exit "Failed to copy PFH-related files to prep/make_PFH."  
fi

#############################
# Stage 1: Graphene Preparation
#############################

echo "Stage 1: Preparing graphene sheets..."
GRAPHENE_DIR="prep/make_graphene"
if [ ! -d "$GRAPHENE_DIR" ]; then
    error_exit "Directory $GRAPHENE_DIR does not exist."
fi

pushd "$GRAPHENE_DIR" > /dev/null

./generate_graphene.sh ${SHEETS} -10 10 -5 5 || error_exit "generate_graphene.sh failed."
echo "Graphene input file (in.run) generated."

if [ ! -x "$LMP_MPI" ]; then
    error_exit "LAMMPS executable not found or not executable at $LMP_MPI"
fi

${LMP_MPI} -in "in.run" > lammps_graphene.log 2>&1 || error_exit "LAMMPS simulation for graphene failed."
echo "LAMMPS simulation completed: out.graphene generated (log in lammps_graphene.log)."
popd > /dev/null

#############################
# Stage 1.5: PFH Slab Dimensions & Molecule Count Setup
#############################

echo "Stage 1.5: Extracting graphene dimensions for PFH slab setup..."
GRAPHENE_OUT="prep/make_graphene/out.graphene"
if [ ! -f "$GRAPHENE_OUT" ]; then
    error_exit "$GRAPHENE_OUT not found."
fi

read gxlo gxhi <<< $(grep "xlo xhi" "$GRAPHENE_OUT" | awk '{print $1, $2}')
read gylo gyhi <<< $(grep "ylo yhi" "$GRAPHENE_OUT" | awk '{print $1, $2}')

x_bound="$gxhi"
y_bound="$gyhi"

echo "Extracted lateral dimensions from out.graphene:"
echo "  xlo = $gxlo, xhi = $gxhi  -->  |x| = $x_bound Å"
echo "  ylo = $gylo, yhi = $gyhi  -->  |y| = $y_bound Å"

half_slab=$(echo "scale=4; $slab_width/2" | bc)
echo "PFH slab thickness (Lz) is set to ${slab_width} Å (half slab = ${half_slab} Å)."

Lx_total=$(echo "$x_bound * 2" | bc -l)
Ly_total=$(echo "$y_bound * 2" | bc -l)

pushd prep/make_PFH > /dev/null
N_required=$(python3 calc.PFH_slab_len.py "$Lx_total" "$Ly_total" "$slab_width")
echo "Number of PFH molecules required (rounded to nearest 100): $N_required"

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^number .*/number $N_required/" packmol_files/in.pfh_packmol       
    sed -i '' "s/xxx/$N_required/"               moltemplate_files/system.lt
    echo "Route 1 : Updated in.pfh_packmol with number of PFH molecules: $N_required"
else
    sed -i "s/^number .*/number $N_required/" packmol_files/in.pfh_packmol
    sed -i "s/xxx/$N_required/"               moltemplate_files/system.lt
    echo "Route 2 : Updated in.pfh_packmol with number of PFH molecules: $N_required"
fi

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^inside box .*/inside box -${x_bound} -${y_bound} -${half_slab} ${x_bound} ${y_bound} ${half_slab}/" packmol_files/in.pfh_packmol
    sed -i '' "s/XLOW/-${x_bound}/"  moltemplate_files/system.lt
    sed -i '' "s/XHIGH/${x_bound}/"  moltemplate_files/system.lt
    sed -i '' "s/YLOW/-${y_bound}/"  moltemplate_files/system.lt
    sed -i '' "s/YHIGH/${y_bound}/"  moltemplate_files/system.lt
    sed -i '' "s/ZLOW/-${half_slab}/"  moltemplate_files/system.lt
    sed -i '' "s/ZHIGH/${half_slab}/"  moltemplate_files/system.lt
else
    sed -i "s/^inside box .*/inside box -${x_bound} -${y_bound} -${half_slab} ${x_bound} ${y_bound} ${half_slab}/" packmol_files/in.pfh_packmol
    sed -i "s/XLOW/-${x_bound}/"    moltemplate_files/system.lt
    sed -i "s/XHIGH/${x_bound}/"    moltemplate_files/system.lt
    sed -i "s/YLOW/-${y_bound}/"    moltemplate_files/system.lt
    sed -i "s/YHIGH/${y_bound}/"    moltemplate_files/system.lt
    sed -i "s/ZLOW/-${half_slab}/"  moltemplate_files/system.lt
    sed -i "s/ZHIGH/${half_slab}/"  moltemplate_files/system.lt
fi
echo "Updated in.pfh_packmol with box dimensions: inside box -${x_bound} -${y_bound} -${half_slab} ${x_bound} ${y_bound} ${half_slab}"


cd packmol_files/
./packmol < in.pfh_packmol > packmol_PFH.log 2>&1 || error_exit "Packmol run failed."
echo "Packmol output generated: out.pfh_packmol."
cp out.pfh_packmol ../moltemplate_files/
cd ../moltemplate_files

echo "Running moltemplate related files..."

outFile="out.pfh_packmol"
${molTemp} -xyz ${outFile} system.lt 

echo "Lets cleanup moltemplate related files."    
${cleanup}
echo "Cleanup of moltemplate files done."    

echo "copying all relevant files."    
mv -f system.data system.in* ../
rm -rf output_ttree/
rm -f run.in.EXAMPLE  
# mv -f warning*.txt ../ 2> /dev/null
echo "done and moving on"

cd ${home}

# popd > /dev/null

#############################
# Stage 2: Combining Graphene and PFH
#############################

echo "Stage 2: Combining graphene sheets and PFH slab..."
mkdir -p prep/combine || error_exit "Failed to create prep/combine directory."

# Copy out.graphene as in.graphene and out.PFH as in.PFH into prep/combine.
cp prep/make_graphene/out.graphene prep/combine/in.graphene || error_exit "Failed to copy out.graphene as in.graphene"
cp prep/make_PFH/system.data prep/combine/in.PFH || error_exit "Failed to copy out.PFH as in.PFH"
echo "Copied in.graphene and in.PFH to prep/combine."

# Extract z-dimensions from out.graphene.
read gzlo gzhi <<< $(grep "zlo zhi" prep/make_graphene/out.graphene | awk '{print $1, $2}')
grafboxlen=$(echo "scale=4; 2*$gzhi" | bc -l)
pfhboxlen=$(echo "${slab_width}+2.5" | bc -l)

Lz_total=$(echo "scale=4; 2*$grafboxlen + $pfhboxlen" | bc -l)
zL=$(echo "scale=4; -0.5*$Lz_total" | bc -l)
zH=$(echo "scale=4; 0.5*$Lz_total" | bc -l)

echo "Graphene sheet spans ${grafboxlen} angstroms"
echo "PFH slab spans ${pfhboxlen} angstroms"
echo "zL equal ${zL} angstroms"
echo "zH equal ${zH} angstroms"

xL_val=$(echo "-1 * $x_bound" | bc -l)
xH_val="$x_bound"
yL_val=$(echo "-1 * $y_bound" | bc -l)
yH_val="$y_bound"

offset_first=3
offset_second=$(echo "$SHEETS + 3" | bc)

shift3=$(echo "scale=4; -0.5*$grafboxlen - 0.5*$pfhboxlen" | bc -l)
shift1=0.0

# Compute the total number of atom types: 2*SHEETS+3.
# Terminal carbon in PFH - Type 1
# Carbon in the backbone of PFH - Type 2 
# Fluorine in PFH - Type 3 
# 9 graphene sheets * 2 = 18 types 
# Total = 2*SHEETS + 3 = 21 
atom_types=$(( SHEETS*2 + 3 ))

cat << EOF > prep/combine/in.run
############################################################################################
# Combined LAMMPS input script for graphene-PFH system
# Generated by graphene_PFH_npt_workflow.sh
############################################################################################
log             log.combined
dimension       3
units           metal 
boundary        p p p
atom_style      full
echo            both
processors      * * 1

variable grafboxlen equal ${grafboxlen}
variable pfhboxlen  equal ${pfhboxlen}

variable shift1 equal ${shift1}
variable shift3 equal ${shift3}

variable xL equal ${xL_val}
variable xH equal ${xH_val}
variable yL equal ${yL_val}
variable yH equal ${yH_val}
variable zL equal ${zL}
variable zH equal ${zH}

region mybox block \${xL} \${xH} \${yL} \${yH} \${zL} \${zH} units box
create_box ${atom_types} mybox bond/types 2 angle/types 3 dihedral/types 3 extra/bond/per/atom 4 extra/angle/per/atom 6 extra/dihedral/per/atom 11

include in.varDefs
include in.forcefield

# Read data for the bottom graphene block.
read_data in.graphene add append offset ${offset_first} 0 0 0 0 shift 0 0 \${shift3}
# Read data for the PFH slab.
read_data in.PFH add append offset 0 0 0 0 0 shift 0 0 \${shift1}
# Read data for the top graphene block (offset = SHEETS+2).
read_data in.graphene add append offset ${offset_second} 0 0 0 0 shift 0 0 \${shift3}
EOF

start_type=$(echo "$SHEETS + 4" | bc)
end_type=$(echo "2*$SHEETS + 3" | bc)
rotate_types=""
for t in $(seq $start_type $end_type); do
    rotate_types="$rotate_types $t"
done

cat << EOF >> prep/combine/in.run
group rotateSheet type ${rotate_types}
displace_atoms rotateSheet rotate 0 0 0 1 0 0 180

# charges for atoms in PFH
set type 1 charge 0.36      # C1
set type 2 charge 0.24      # C2
set type 3 charge -0.120    # F 
set type 4* charge 0.0000   # C-Tersoff

# Energy minimization
minimize 1.0e-6 1.0e-8 1000 10000

write_data combined_system.dat

quit
EOF

echo "Combined input file generated in prep/combine/in.run."

cd ${home}
cp runFiles/npt/in.SiC.tersoff prep/combine/
cp runFiles/npt/in.forcefield  prep/combine/
cp runFiles/npt/in.varDefs     prep/combine/

#############################
# Stage 4: Running Final Simulation with lmp_mpi
#############################
echo "Stage 4: Running final simulation using lmp_mpi on 10 cores..."
pushd prep/combine > /dev/null
${LMP_MPI} -in in.run -l log.sim 2>&1  || error_exit "LAMMPS simulation for graphene failed."
popd > /dev/null
echo "Final simulation completed successfully."
