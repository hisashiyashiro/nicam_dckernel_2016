#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N dckernel_divdamp
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_divdamp/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_divdamp.exe .
ln -svf ${HMDIR}/dckernel_divdamp/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/dckernel_divdamp/data/snapshot.dc_divdamp3d.pe000000 .

./dckernel_divdamp.exe
