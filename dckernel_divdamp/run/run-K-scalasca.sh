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
. /work/aics_apps/scalasca/Env_scalasca
/opt/FJSVXosPA/bin/xospastop
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2
export SCAN_ANALYZE_OPTS="-i -s"
metrics="L1_MISS:L1_I_MISS:L1_D_MISS:L2_MISS:TLB_MISS:TLB_I_MISS:TLB_D_MISS:FLOATING_POINT"

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_divdamp/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_divdamp.exe .
ln -svf ${HMDIR}/dckernel_divdamp/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_divdamp/data/snapshot.dc_divdamp3d.pe000000 .

rm -rf ./epik_trace

scan -t -m ${metrics} -e epik_trace ./dckernel_divdamp.exe
