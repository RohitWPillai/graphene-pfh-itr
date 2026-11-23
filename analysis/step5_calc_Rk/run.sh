#!/bin/bash 

press=( "50" "100" "200" "250" "300" "350" "400" )
# press=( "50" )    # run single case 
nSets=80

dest="./outs"

mkdir -p ${dest}
mkdir -p plots

for p in ${press[@]}
do
    echo "--- p = $p ----"
    #---- sheet temperature ----#
    ./comp.sheet_temp.py ${p} ${nSets}

    #---- pfh temperature ----#
    ./comp.pfh_temp.py ${p} ${nSets}

    #---- pfh densities ----#
    ./comp.pfh_dens.py ${p} ${nSets}

    #---- Rk calculation ----#
    fluxFile="../step1_heat_flux/output.dat"
    flux=$(awk -v p="$p" '$1 == p {print $4}' "${fluxFile}")    # left cold wall position
    echo "${flux}"
   
    ./comp.Rk.block.py ${flux} ${p}     # block stats 
    ./comp.Rk.py ${flux} ${p}           # one average value 

done

# plot Rk as a function of pressure 
