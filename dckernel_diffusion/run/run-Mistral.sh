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

OUTDIR=${HMDIR}/dckernel_diffusion/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_diffusion.exe .
ln -svf ${HMDIR}/dckernel_diffusion/data/snapshot.dc_diffusion.pe000000 .

./dckernel_diffusion.exe
