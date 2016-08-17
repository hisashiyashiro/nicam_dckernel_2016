#! /bin/bash -x
#
# for DKRZ Mistral with SLURM
#
#SBATCH --partition=compute
#SBATCH --job-name=NICAMknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
set -e

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_setup/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_setup.exe .
ln -svf ${HMDIR}/dckernel_setup/data/snapshot.dc_setup.pe000000 .

./dckernel_setup.exe
