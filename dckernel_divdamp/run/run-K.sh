#! /bin/bash -x
#
# for K computer
#
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2
#export fu07bf=1

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_divdamp/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_divdamp.exe .
ln -svf ${HMDIR}/dckernel_divdamp/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_divdamp/data/snapshot.dc_divdamp3d.pe000000 .

rm -rf ./prof*

fapp -C -Ihwm -Hevent=Statistics -d prof -L 10 ./dckernel_divdamp.exe
