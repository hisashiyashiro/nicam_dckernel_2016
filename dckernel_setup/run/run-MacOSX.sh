#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_setup/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_setup.exe .
ln -svf ${HMDIR}/dckernel_setup/data/snapshot.dc_setup.pe000000 .

./dckernel_setup.exe
