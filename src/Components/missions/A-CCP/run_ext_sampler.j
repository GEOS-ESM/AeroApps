#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=00:15:00
#SBATCH --ntasks=12 --cpus-per-task=1 --ntasks-per-node=12
#SBATCH --job-name=lidar_sampler
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
setenv G5BIN /discover/nobackup/pcastell/workspace/GAAS/Linux/bin
setenv BIN $PWD
setenv TMPDIR $LOCAL_TMPDIR 

mkdir -p $LOCAL_TMPDIR/ExtDataOsku
cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtDataOsku/* $LOCAL_TMPDIR/ExtDataOsku
ln -s $LOCAL_TMPDIR/ExtDataOsku ./ExtDataOsku


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
python -u ${G5BIN}/ext_sampler.py --input=INFILE --output=OUTFILE --format=NETCDF4_CLASSIC --channel=CHANNEL --rc=Aod_EOS.rc --intensive > slurm_${SLURM_JOBID}_py.out
rm ExtDataOsku
