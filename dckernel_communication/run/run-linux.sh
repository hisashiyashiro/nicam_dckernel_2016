#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=2
#PBS -N dckernel_communication
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_communication/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_communication.exe .
ln -svf ${HMDIR}/dckernel_communication/conf/communication.cnf .

mpirun -np 2 ./dckernel_communication.exe
