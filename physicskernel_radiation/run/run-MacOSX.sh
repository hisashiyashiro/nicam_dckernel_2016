#! /bin/bash -x
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/physicskernel_radiation/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/physicskernel_radiation.exe .
ln -svf ${HMDIR}/physicskernel_radiation/data/vgrid94.dat .
ln -svf ${HMDIR}/physicskernel_radiation/data/PARA.bnd29ch111sp .

ln -svf ${HMDIR}/physicskernel_radiation/data/snapshot.radiation.pe000003 .
ln -svf ${HMDIR}/physicskernel_radiation/data/check.radiation.pe000003 .

./physicskernel_radiation.exe
