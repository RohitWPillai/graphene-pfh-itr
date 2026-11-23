#!/bin/bash
#SBATCH --job-name=computeTemp
#SBATCH --time=01:00:00
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Replace [budget code] below with your budget code (e.g. <project>)
#SBATCH --account=<project>
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set the number of threads to 1
# This prevents any threaded system libraries from automatically
# using threading.
export OMP_NUM_THREADS=1
# Propagate the cpus-per-task setting from script to srun commands
# By default, Slurm does not propagate this setting from the sbatch
# options to srun commands in the job script. If this is not done,
# process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

start=$(date +%s)

execFile="exec.sheet_temp"
CC -std=c++17 -I/path/to/fftw/include compute_sheetTemp.cpp -o ${execFile} -O3
zlo=-500
zhi=500
blocksize=50000
nDumps=80
press=( "50" "100" "200" "250" "300" "350" "400" )

for p in ${press[@]}
do
    mkdir -p outs/${p}bar
    mkdir -p logs/${p}bar
done

# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do
    for i in $(seq 1 ${nDumps})
    do
        inFile="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.C.lwall"
        outFile="outs/${p}bar/out${i}.C.lwall"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
              --exact --mem=1500M ${execFile} ${inFile} ${outFile} ${zlo} ${zhi} ${blocksize} > logs/${p}bar/log.sim${i} &
    done
done

#wait
#echo "All cold left walls done"
#
# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do
    for i in $(seq 1 ${nDumps})
    do
        inFile="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.C.rwall"
        outFile="outs/${p}bar/out${i}.C.rwall"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
              --exact --mem=1500M ${execFile} ${inFile} ${outFile} ${zlo} ${zhi} ${blocksize} > logs/${p}bar/log.sim${i} &
    done
done

#wait
#echo "All cold right walls done"
#
# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do
    for i in $(seq 1 ${nDumps})
    do
        inFile="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.H.lwall"
        outFile="outs/${p}bar/out${i}.H.lwall"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
              --exact --mem=1500M ${execFile} ${inFile} ${outFile} ${zlo} ${zhi} ${blocksize} > logs/${p}bar/log.sim${i} &
    done
done

#wait
#echo "All hot left walls done"
# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do
    for i in $(seq 1 ${nDumps})
    do
        inFile="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.H.rwall"
        outFile="outs/${p}bar/out${i}.H.rwall"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
              --exact --mem=1500M ${execFile} ${inFile} ${outFile} ${zlo} ${zhi} ${blocksize} > logs/${p}bar/log.sim${i} &
    done
done
wait
end=$(date +%s)
runtime=$((end - start))
echo "Total runtime: $runtime seconds"
echo "All jobs done"
