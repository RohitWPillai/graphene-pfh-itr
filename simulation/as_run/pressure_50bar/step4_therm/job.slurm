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

# Define LAMMPS executable
LMP_MPI="/work/<project>/<project>/<user>/lammps-stable_29Aug2024/src/lmp_mpi"

RND1=1
RND2=2
RND3=3
RND4=4

#------------ Step 3: NVE Simulation ------------#
echo "Starting srun for step4_therm"
srun --hint=nomultithread --distribution=block:block ${LMP_MPI} -in in.run -v "RND1" ${RND1} -v "RND2" ${RND2} -v "RND3" ${RND3} -v "RND4" ${RND4}> log.sim 2>&1
echo "therm finished"
cd ..
