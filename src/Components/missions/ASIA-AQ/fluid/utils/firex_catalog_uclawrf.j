#!/bin/csh -fx

#SBATCH --job-name=firex_catalog_uclawrf
#SBATCH --account=s1321
#SBATCH --time=1:00:00
#SBATCH --qos=daohi
#SBATCH --ntasks=28
#SBATCH --export=NONE
#SBATCH --constraint=hasw
#SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_catalog_UCLAWRF_20190827_0.log

limit stacksize unlimited
source /home/dao_ops/jardizzo/FLUID/firex-aq/utils/pyg_modules
cd /home/dao_ops/jardizzo/FLUID/firex-aq/utils

firex_catalog.sh 20190827 0 UCLAWRFchem 3 12 3
