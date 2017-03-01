#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N physicskernel_radiation
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/physicskernel_radiation/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/physicskernel_radiation.exe .
ln -svf ${HMDIR}/physicskernel_radiation/data/vgrid94.dat .
ln -svf ${HMDIR}/physicskernel_radiation/data/PARA.bnd29ch111sp .

ln -svf ${HMDIR}/physicskernel_radiation/data/snapshot.radiation.pe000003 .
ln -svf ${HMDIR}/physicskernel_radiation/data/check.radiation.pe000003 .

./physicskernel_radiation.exe
