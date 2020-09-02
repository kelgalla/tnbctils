#!/bin/sh

#PBS -M kelgalla@iupui.edu
#PBS -m abe
#PBS -l vmem=10gb,walltime=01:00:00
#PBS -N wrapper
#PBS -o wrapper.out
#PBS -e wrapper.err
#PBS -S /bin/sh
#PBS -V

cd $PBS_O_WORKDIR

perl findDataType.pl