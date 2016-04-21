#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_diffusion/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_diffusion.exe .
ln -svf ${HMDIR}/dckernel_diffusion/data/snapshot.dc_diffusion.pe000000 .

./dckernel_diffusion.exe
