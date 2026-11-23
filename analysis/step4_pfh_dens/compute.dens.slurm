#!/bin/bash 

#SBATCH --job-name=computeDens
#SBATCH --time=3:00:00
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Replace [budget code] below with your budget code (e.g. <project>)
#SBATCH --account=<project>
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Propagate the cpus-per-task setting from script to srun commands
#    By default, Slurm does not propagate this setting from the sbatch
#    options to srun commands in the job script. If this is not done,
#    process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

start=$(date +%s)

execFile="exec.dens"

if [ ! -f ${execFile} ];then
    echo "${execFile} not found, so compiling"
    CC computeDens.24Jun25.cpp -o ${execFile} -O3
fi

# zHot=$(python3 find_wall_lims.py | grep 'left wall hot'  | awk '{print $5}' )      # hot wall on the left domain
# zCold=$(python3 find_wall_lims.py | grep 'left wall cold' | awk '{print $5}' )      # cold wall on the right domain
nbins=5000
blocksize=50000
nDumps=80
press=( "50" "100" "200" "250" "300" "350" "400" )
posFile="../step2_sheet_pos_and_temp/positions_means.dat"   # use only absolute paths 

mkdir -p outs
mkdir -p logs 

# create directories 
for p in ${press[@]}
do
    mkdir -p outs/${p}bar
    mkdir -p logs/${p}bar
done

echo "----Left wall-----"
domainSide="L"
# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do

    zCold=$(awk -v p="$p" '$1 == p {print $2}' "${posFile}")    # right cold wall position
     zHot=$(awk -v p="$p" '$1 == p {print $3}' "${posFile}")    # left hot wall position (zLo is always hot)

    echo "------- Pressure=${p} bar -------"
    echo "zHot=${zHot} and zCold=${zCold} for the left domain"

    for i in $(seq 1 ${nDumps}) 
    do

        in="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.lpfh"
        out="outs/${p}bar/out${i}.lpfh.dens"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
             --exact --mem=1500M ${execFile} ${in} ${out} ${zHot} ${zCold} ${nbins} ${blocksize} ${domainSide} > logs/${p}bar/log.lpfh.sim${i} &
    done
    echo "----------------------------------"
    echo "                                  "
done

echo "----Right wall-----"
domainSide="R"
# Loop over 80 subjobs pinning each to a different core
for p in ${press[@]}
do

    zHot=$(awk -v p="$p" '$1 == p {print $4}' "${posFile}")    # left cold wall position
    zCold=$(awk -v p="$p" '$1 == p {print $5}' "${posFile}")    # right hot wall position

    echo "------- Pressure=${p} bar -------"
    echo "zCold=${zCold} and zHot=${zHot} for the left domain"

    for i in $(seq 1 ${nDumps}) 
    do

        in="/scratch-nvme/<project>/<project>/<user>/PFH/${p}bar/rep_1/dump${i}.rpfh"
        out="outs/${p}bar/out${i}.rpfh.dens"
        srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
             --exact --mem=1500M ${execFile} ${in} ${out} ${zHot} ${zCold} ${nbins} ${blocksize} ${domainSide} > logs/${p}bar/log.rpfh.sim${i} &
    done
    echo "----------------------------------"
    echo "                                  "
done

# Wait for all subjobs to finish
wait
end=$(date +%s)
runtime=$((end - start))
echo "Total runtime: $runtime seconds"
echo "All jobs done"
