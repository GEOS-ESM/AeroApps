#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --job-name=accp_vlidort
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
source setup_env
setenv TMPDIR $LOCAL_TMPDIR 

mkdir -p $LOCAL_TMPDIR/ExtData
cp -r ${x}/ExtData/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtData ./ExtData

cp -r ${x}/ExtDataOsku/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtDataOsku ./ExtDataOsku

cp -r ${AERODIR}/install/bin/missions/ACCP/Chem_MieRegistry.rc $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/Chem_MieRegistry.rc ./Chem_MieRegistry.rc

##################################################################
######
######         Do Sampling
######
##################################################################
python -u ./accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T00:00 2006-01-02T00:00 lidar_files.pcf gpm.pcf polar07.pcf 550 > slurm_${SLURM_JOBID}_py.out &
rm -rf $LOCAL_TMPDIR/*
rm ExtData
rm ExtDataOsku
rm Chem_MieRegistry.rc
