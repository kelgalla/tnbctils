#!/bin/sh

#PBS -M kelgalla@iupui.edu
#PBS -m abe
#PBS -l vmem=10gb,walltime=24:00:00
#PBS -N gdc-data
#PBS -o ../gdc-data.out
#PBS -e ../gdc-data.err
#PBS -S /bin/sh
#PBS -V

cd $PBS_O_WORKDIR

gdc-client download -m ../gdc_manifest.2017-05-10T00_30_13.024848.tsv -t ../gdc-user-token.2017-05-15T18-35-24-04-00.txt