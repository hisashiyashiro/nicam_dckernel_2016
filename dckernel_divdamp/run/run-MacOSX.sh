#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_divdamp/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_divdamp.exe .
ln -svf ${HMDIR}/dckernel_divdamp/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_divdamp/data/snapshot.dc_divdamp3d.pe000000 .

./dckernel_divdamp.exe
