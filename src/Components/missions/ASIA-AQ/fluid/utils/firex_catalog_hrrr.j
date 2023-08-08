#!/bin/csh -fx

#SBATCH --job-name=firex_catalog_hrrr
#SBATCH --account=s1321
#SBATCH --time=1:00:00
#SBATCH --qos=daohi
#SBATCH --ntasks=28
#SBATCH --export=NONE
#SBATCH --constraint=hasw
#SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_catalog_HRRR_20190827_0.log

limit stacksize unlimited
source /home/dao_ops/jardizzo/FLUID/firex-aq/utils/pyg_modules
cd /home/dao_ops/jardizzo/FLUID/firex-aq/utils

firex_catalog.sh 20190827 0 HRRR 2 0 3
