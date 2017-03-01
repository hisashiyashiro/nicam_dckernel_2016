#! /bin/bash -x
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/physicskernel_microphysics/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/physicskernel_microphysics.exe .
ln -svf ${HMDIR}/physicskernel_microphysics/data/vgrid94.dat .

ln -svf ${HMDIR}/physicskernel_microphysics/data/snapshot.microphysics.pe000003 .
ln -svf ${HMDIR}/physicskernel_microphysics/data/check.microphysics.pe000003 .

./physicskernel_microphysics.exe
