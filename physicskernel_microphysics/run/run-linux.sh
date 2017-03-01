#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N physicskernel_microphysics
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/physicskernel_microphysics/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/physicskernel_microphysics.exe .
ln -svf ${HMDIR}/physicskernel_microphysics/data/vgrid94.dat .

ln -svf ${HMDIR}/physicskernel_microphysics/data/snapshot.microphysics.pe000003 .
ln -svf ${HMDIR}/physicskernel_microphysics/data/check.microphysics.pe000003 .

./physicskernel_microphysics.exe
