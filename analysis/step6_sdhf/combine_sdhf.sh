#!/bin/bash 

EXEC="combine"
CC combine_sdhfs.cpp -o ${EXEC} -O3 -std=c++11

startSeg=1
endSeg=3198
solid=1
fluid=0

#----------- Solid -----------------#
if [ "${solid}" -eq 1 ]
then
    press=( "50bar" "100bar" "200bar" "250bar" "300bar" "350bar" "400bar" )
    src="outs/solid"
    exts="H.rwall,H.lwall"
    outsuf="solid"

    for p in ${press[@]}
    do
        ./${EXEC} ${src}/${p}/ ${exts} 2500 1 3198 sdhf_consolidated/${p}_${outsuf}.dat
    done
fi

#----------- Fluid -----------------#
if [ "${fluid}" -eq 1 ]
then
    press=( "50bar" "100bar" "200bar" "250bar" "300bar" "350bar" "400bar" )
    src="outs/fluid"
    exts="lpfh,rpfh"
    outsuf="fluid"

    for p in ${press[@]}
    do
        ./${EXEC} ${src}/${p}/ ${exts} 2500 1 3198 sdhf_consolidated/${p}_${outsuf}.dat
    done
fi
