#! /bin/bash -x
#
# for DKRZ Mistral
#
#SBATCH --partition=compute
#SBATCH --job-name=NICAMknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:10:00
set -e

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_communication/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_communication.exe .
ln -svf ${HMDIR}/dckernel_communication/conf/communication.cnf .

mpirun -np 2 ./dckernel_communication.exe
