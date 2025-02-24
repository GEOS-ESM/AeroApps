#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=120 --cpus-per-task=1 --ntasks-per-node=120
#SBATCH --constraint=mil
#SBATCH --job-name=sbg_vlidort
#SBATCH -A s2190
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
source $HOME/.cshrc
source ./setup_env

setenv TMPDIR $LOCAL_TMPDIR
setenv x /discover/nobackup/pcastell/workspace/aero_work/aist/sbg

mkdir -p $LOCAL_TMPDIR/ExtData_oldRC
cp -r ${x}/ExtData_oldRC/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtData_oldRC ./ExtData_oldRC

cp -r ${x}/ExtDataOsku/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtDataOsku ./ExtDataOsku

cp -r ${AERODIR}/install/bin/missions/ACCP/Chem_MieRegistry.rc $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/Chem_MieRegistry.rc ./Chem_MieRegistry.rc

##################################################################
######
######         Do Sampling
######
##################################################################
python -u sbg_vlidort.py -D 1 2006-01-16T17:32 2006-01-16T17:33 imager_files.pcf ssd650.pcf sbg.pcf 200 > slurm_${SLURM_JOBID}_py.out
rm -rf $LOCAL_TMPDIR/*
rm ExtData_oldRC
rm ExtDataOsku
rm Chem_MieRegistry.rc
rm -rf rc
