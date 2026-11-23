#!/bin/bash 

press=( "50" "100" "200" "250" "300" "350" "400" )

mkdir -p qInOut

for p in ${press[@]}
do
    scp /work/<project>/<project>/<user>/graf_pfh_workflow/pressure_${p}/step5_dump/iter4/out.qInOut ./qInOut/out.${p}_bar
done
