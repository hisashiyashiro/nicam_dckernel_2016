#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N dckernel_vi_rhow_solver
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_vi_rhow_solver/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_vi_rhow_solver.exe .
ln -svf ${HMDIR}/dckernel_vi_rhow_solver/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_vi_rhow_solver/data/snapshot.dc_vi_rhow_solver.pe000000 .

./dckernel_vi_rhow_solver.exe
