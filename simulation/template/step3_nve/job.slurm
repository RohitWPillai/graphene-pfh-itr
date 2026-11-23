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

#------------ Step 3: NVE Simulation ------------#
echo "Starting srun for step3_nve"
xH=$(grep "xlo xhi :" ../step2_nvt/log.lammps | tail -1 | awk '{print $5}')
yH=$(grep "ylo yhi :" ../step2_nvt/log.lammps | tail -1 | awk '{print $5}')
zH=$(grep "zlo zhi :" ../step2_nvt/log.lammps | tail -1 | awk '{print $5}')
srun --hint=nomultithread --distribution=block:block ${LMP_MPI} -in in.run -v "halfLz" ${zH} -v "xH" ${xH} -v "yH" ${yH} > log.sim 2>&1
echo "nve finished"
cd ..
