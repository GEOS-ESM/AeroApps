#!/bin/csh -fx

#SBATCH --job-name=firex_catalog_geoscf
#SBATCH --account=s1321
#SBATCH --time=2:00:00
#SBATCH --qos=daohi
#SBATCH --ntasks=28
#SBATCH --export=NONE
#SBATCH --constraint=hasw
#SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_catalog_GEOSCF_20190826_12.log

limit stacksize unlimited
source /home/dao_ops/jardizzo/FLUID/firex-aq/utils/pyg_modules
cd /home/dao_ops/jardizzo/FLUID/firex-aq/utils

firex_catalog.sh 20190826 120000 GEOSCF 5 0 3
