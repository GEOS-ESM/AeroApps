#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=accp_sdist
#SBATCH --account=s2190
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

#######################################################################
#                    USER INPUTS
#######################################################################
setenv AEROBIN $PWD 


#######################################################################
#                  System Environment Variables
#######################################################################
umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src

mkdir -p $LOCAL_TMPDIR/ExtData
cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtData/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtData ./ExtData

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtDataOsku/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/ExtDataOsku ./ExtDataOsku

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/Chem_MieRegistry.rc $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/Chem_MieRegistry.rc ./Chem_MieRegistry.rc

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/rc/Aod_EOS.Osku.rc $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/Aod_EOS.Osku.rc  ./Aod_EOS.rc


source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $AEROBIN

##################################################################
######
######         Run the size distribution code
######
##################################################################
./accp_sdist.py 2006-02-01T00 2006-02-01T04 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-01T04 2006-02-01T08 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-01T08 2006-02-01T12 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-01T12 2006-02-01T16 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-01T16 2006-02-01T20 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-01T20 2006-02-02T00 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T00 2006-02-02T04 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T04 2006-02-02T08 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T08 2006-02-02T12 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T12 2006-02-02T16 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T16 2006-02-02T20 lidar_files.pcf gpm.pcf &
./accp_sdist.py 2006-02-02T20 2006-02-03T00 lidar_files.pcf gpm.pcf 


rm -rf $LOCAL_TMPDIR/*
rm ExtData
rm ExtDataOsku
rm Chem_MieRegistry.rc
rm Aod_EOS.rc
