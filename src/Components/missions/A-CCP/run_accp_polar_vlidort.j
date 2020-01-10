#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=01:00:00
#SBATCH --ntasks=8 --cpus-per-task=1 --ntasks-per-node=8
#SBATCH --job-name=accp_vlidort
#SBATCH -A s1180
#SBATCH --qos=pace
#SBATCH --partition=preops
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv BIN $PWD 
setenv TMPDIR $LOCAL_TMPDIR 

mkdir -p $LOCAL_TMPDIR/ExtData
cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtData/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtData ./ExtData

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtDataOsku/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtDataOsku ./ExtDataOsku

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/Chem_MieRegistry.rc $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/Chem_MieRegistry.rc ./Chem_MieRegistry.rc


source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $BIN

##################################################################
######
######         Do Sampling
######
##################################################################
python -u mp_run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T00:00 2006-01-02T00:00 lidar_files.pcf gpm.pcf polar07.pcf 550 > slurm_${SLURM_JOBID}_py.out &

rm ExtData
rm ExtDataOsku
rm Chem_MieRegistry.rc
