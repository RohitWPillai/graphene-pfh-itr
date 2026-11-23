#!/bin/bash 

module load cray-fftw

EXEC="./exec_solid"
#CC -I/opt/cray/pe/fftw/3.3.8.8/x86_rome/include solid_sdhf.cpp utils.cpp -o ${EXEC} -O3 -std=c++11

${EXEC} "/scratch-nvme/<project>/<project>/<user>/graf_pfh/50bar/rep_1/" "H.lwall" "H.lwall" 2500 50 1 0
