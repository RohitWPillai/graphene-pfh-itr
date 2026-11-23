#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=long
#SBATCH --account=<project>
#SBATCH --nodes=3
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=pfh_npt
#SBATCH --output=log.%x.therm
#SBATCH --error=log.%x.error

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export SLURM_CPU_FREQ_REQ=2250000

press=100

# Define LAMMPS executable
LMP_MPI="/work/<project>/<project>/<user>/lammps-stable_29Aug2024/src/lmp_mpi"

#------------ Step 1: NPT Equilibration ------------#
echo "Starting srun for step1_npt"
srun --hint=nomultithread --distribution=block:block ${LMP_MPI} -in in.run -v "targetPress" ${press} 
echo "npt finished"
if [ $? -ne 0 ] || [ ! -f "config.end_of_npt" ]; then
    echo "ERROR: NPT step failed. Check step1_npt/log.sim"
    exit 1
fi
