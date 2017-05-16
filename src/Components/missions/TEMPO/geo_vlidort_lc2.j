#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=00:30:00
#SBATCH --time-min=00:10:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=tempo_geo_vlidort
#SBATCH -A s1412
#SBATCH --array=1-3
#SBATCH --output=slurm_%A_%a.out
#SBATCH --error=slurm_%A_%a.err
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv GEOBIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/TEMPO 
setenv OUTDIR /discover/nobackup/pcastell/TEMPO/DATA
setenv ADDOUTDIR /discover/nobackup/pcastell/TEMPO/DATA
setenv RUN_CMD  "mpirun -np 28"

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $GEOBIN

##################################################################
######
######         Perform single iteration of VLIDORT Run
###### Make sure to clean up shared memory segements before and after the job
######
##################################################################
./clean_mem.sh
$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc ${SLURM_ARRAY_TASK_ID}
./clean_mem.sh 

#rm Aod_EOS.rc
#rm Chem_MieRegistry.rc
#rm ExtData
#rm geo_vlidort.x
#rm clean_mem.sh
#rm geo_vlidort.rc
#set size=`stat -c %s slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err`
#set DONT_MOVE=0
#foreach size (`stat -c %s slurm_${SLURM_ARRAY_JOB_ID}*.err`)
#	if ($size > 0) then
#		set DONT_MOVE=1
#	endif
#end


#if ($DONT_MOVE == 0) then
#	foreach fn (`ls tempo-g5nr.lc2.*${SLURM_ARRAY_TASK_ID}.nc4`)
#	  mv $fn $OUTDIR/
#	end
#
#	foreach fn (`ls tempo-g5nr.lc.*${SLURM_ARRAY_TASK_ID}.nc4`)
#	  mv $fn $ADDOUTDIR/
#	end
#endif
