#!/bin/bash 

EXEC="exec.vdos"
#CC -I/path/to/fftw/include calcVDOSPFH.11Aug25.cpp -o ${EXEC} -O3 &

#rm nohup.out
IN="/scratch-nvme/<project>/<project>/<user>/graf_pfh/50bar/rep_1/dump1.rpfh"
ZLO=29
ZHI=39
SEGLEN=2000
dt="1e-15"

SEGNUM=1
OUT="out.${SEGNUM}.vdos"
./${EXEC} ${IN} ${OUT} ${ZLO} ${ZHI} ${dt} ${SEGLEN} ${SEGNUM} &

# SEGNUM=2
# OUT="out.${SEGNUM}.vdos"
# ${EXEC} ${IN} ${OUT} ${ZLO} ${ZHI} ${dt} ${SEGLEN} ${SEGNUM} &

# SEGNUM=3
# OUT="out.${SEGNUM}.vdos"
# ${EXEC} ${IN} ${OUT} ${ZLO} ${ZHI} ${dt} ${SEGLEN} ${SEGNUM} &

# SEGNUM=4
# OUT="out.${SEGNUM}.vdos"
# ${EXEC} ${IN} ${OUT} ${ZLO} ${ZHI} ${dt} ${SEGLEN} ${SEGNUM} &

# SEGNUM=5
# OUT="out.${SEGNUM}.vdos"
# ${EXEC} ${IN} ${OUT} ${ZLO} ${ZHI} ${dt} ${SEGLEN} ${SEGNUM} &
