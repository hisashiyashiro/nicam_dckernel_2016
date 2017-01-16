#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_communication/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_communication.exe .
ln -svf ${HMDIR}/dckernel_communication/conf/communication.cnf .

mpirun -np 2 ./dckernel_communication.exe
