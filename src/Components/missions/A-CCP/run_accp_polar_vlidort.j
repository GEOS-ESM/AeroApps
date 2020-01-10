#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=01:00:00
#SBATCH --ntasks=11 --ntasks-per-node=1
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

mkdir -p $LOCAL_TMPDIR/rc
cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/rc/* $LOCAL_TMPDIR/rc
ln -s $LOCAL_TMPDIR/rc ./rc

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
srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_550.pcf > slurm_${SLURM_JOBID}_py_550.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_670.pcf > slurm_${SLURM_JOBID}_py_670.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_360.pcf > slurm_${SLURM_JOBID}_py_360.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_380.pcf > slurm_${SLURM_JOBID}_py_380.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_410.pcf > slurm_${SLURM_JOBID}_py_410.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_870.pcf > slurm_${SLURM_JOBID}_py_870.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_940.pcf > slurm_${SLURM_JOBID}_py_940.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_1230.pcf > slurm_${SLURM_JOBID}_py_1230.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_1380.pcf > slurm_${SLURM_JOBID}_py_1380.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_1550.pcf > slurm_${SLURM_JOBID}_py_1550.out &

srun -N 1 -n 1 python -u run_accp_polar_vlidort.py -v --DT_hours 1 2006-01-01T01 2006-01-01T02 accp_polar_vlidort_1650.pcf > slurm_${SLURM_JOBID}_py_1650.out &

wait

rm rc
rm ExtData
rm ExtDataOsku
rm Chem_MieRegistry.rc
