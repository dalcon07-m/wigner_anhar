#!/bin/bash
#PBS -u macias
#PBS -N vib6_-0.10
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
#PBS -m be
#PBS -r n

FILE="vib6_-0.10"
ScrDir="/scr/macias/${PBS_JOBID}_$FILE"
Wdir="/home/macias/hoso/wigner_anh/vib6/-0.10"
. /soft/g09.c01/g09/bsd/g09.profile

mkdir -p $ScrDir
cd $ScrDir
g09 < $Wdir/$FILE.com > $ScrDir/$FILE.log
cp *.log $Wdir
cp *.o $Wdir
cp *.e $Wdir
exit
