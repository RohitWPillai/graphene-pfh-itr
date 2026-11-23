#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=<project>
#SBATCH --job-name=350bar
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=log.%x.pfh_sdhf
#SBATCH --error=err.%x.pfh_sdhf
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Print the start time
echo "Job started at: $(date)"

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Propagate the cpus-per-task setting from script to srun commands
#    By default, Slurm does not propagate this setting from the sbatch
#    options to srun commands in the job script. If this is not done,
#    process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load cray-fftw

# Compile once with OpenMP support
EXEC="exec_solid"
# CC -I/opt/cray/pe/fftw/3.3.8.8/x86_rome/include solid_sdhf.cpp utils.cpp -o ${EXEC} -O3 -std=c++11

# --- User-configurable roots & parameters ---
SRC="/scratch-nvme/<project>/<project>/<user>/graf_pfh/"         # base directory containing pressure folders
PRESSURES=("350bar" )      # process all pressure cases
SUFFIXES=(H.lwall H.rwall)                          # process left- and right-wall dumps
SEGLEN=2500                                         # segment length in frames
LOGDIR="logs/solid"
OUTDIR="outs/solid"
NDUMPS=80
FRAMES=50001
SEGSTART=1
SEGEND=3198
# NATOMS=7200
# SEGMAX equals 1598 (for SEGLEN=5000)
# SEGMAX equals 3198 (for SEGLEN=5000)

mkdir -p "${LOGDIR}"
mkdir -p "${OUTDIR}"

# Launch one multi-threaded task per (pressure, suffix)
for p in "${PRESSURES[@]}"; do
  for suf in "${SUFFIXES[@]}"; do
      for i in $(seq ${SEGSTART} ${SEGEND}); do

        echo "Running sdhf analysis for ${p} (${suf})..."

        LOG="${LOGDIR}/${p}_${i}_${suf}.log"
        # IN="${SRC}/${p}/rep_1/dump${i}.${suf}"
        IN="${SRC}/${p}/rep_1/"
        # OUT="${OUTDIR}/${p}_${i}_${suf}.dat"

       srun --nodes=1 --ntasks=1 --exact --cpus-per-task=$SLURM_CPUS_PER_TASK --mem-per-cpu=1500M \
         ./${EXEC} "${IN}" "${suf}" "outs/solid/${p}/${suf}" "${SEGLEN}" "50" "1" "${i}" > "$LOG" &
     done
  done
done

wait
echo "All SDHF analyses completed."

# mv *pfh ${OUTDIR}/

# Print the end time
echo "Job ended at: $(date)"
