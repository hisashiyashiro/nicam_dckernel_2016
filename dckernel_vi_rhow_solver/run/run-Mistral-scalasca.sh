#! /bin/bash -x
#
# for DKRZ Mistral
#
#SBATCH --partition=compute
#SBATCH --job-name=NICAMknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
set -e

export OMP_NUM_THREADS=4
#export SCOREP_METRIC_PAPI=PAPI_FP_OPS

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_vi_rhow_solver/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_vi_rhow_solver.exe .
ln -svf ${HMDIR}/dckernel_vi_rhow_solver/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_vi_rhow_solver/data/snapshot.dc_vi_rhow_solver.pe000000 .

rm -rf ./epik_trace

scan -t -e epik_trace ./dckernel_vi_rhow_solver.exe
