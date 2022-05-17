#!/bin/csh
#PBS -S /bin/csh
#PBS -N 750R200bE9
#PBS -l nodes=moore10:ppn=1
#PBS -o mainRDME.out
#PBS -e mainRDME.err


cd $PBS_O_WORKDIR

./rdme
