#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --time-min=08:00:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --job-name=tempo_geo_vlidort
#SBATCH -A s1412
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv TEMPOBIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/TEMPO 
setenv OUTDIR /discover/nobackup/pcastell/TEMPO/DATA
setenv RUN_CMD  "mpirun -np 28"

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $TEMPOBIN

##################################################################
######
######         Perform single iteration of VLIDORT Run
###### Make sure to clean up shared memory segements before and after the job
######
##################################################################
./clean_mem.sh
$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc
./clean_mem.sh 

rm Aod_EOS.rc
rm Chem_MieRegistry.rc
rm ExtData
rm geo_vlidort.x
#rm geo_vlidort.rc
mv *nc4 $OUTDIR