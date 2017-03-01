#! /bin/bash -x
#
# for K(micro)
#
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:30:00"
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

OUTDIR=${HMDIR}/physicskernel_radiation/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/physicskernel_radiation.exe .
ln -svf ${HMDIR}/physicskernel_radiation/data/vgrid94.dat .
ln -svf ${HMDIR}/physicskernel_radiation/data/PARA.bnd29ch111sp .

ln -svf ${HMDIR}/physicskernel_radiation/data/snapshot.radiation.pe000003 .
ln -svf ${HMDIR}/physicskernel_radiation/data/check.radiation.pe000003 .

rm -rf ./prof*

#fapp -C -Ihwm -Hevent=Statistics -d prof -L 10 ./physicskernel_radiation.exe || exit 1
./physicskernel_radiation.exe || exit 1
