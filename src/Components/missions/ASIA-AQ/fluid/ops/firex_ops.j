#!/bin/csh -fx

#SBATCH --job-name=firex_ops
#SBATCH --account=s1321
#SBATCH --time=1:00:00
#SBATCH --qos=daohi
#SBATCH --ntasks=28
#SBATCH --export=NONE
#SBATCH --constraint=hasw
#SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_ops_20190713_0.log

limit stacksize unlimited
source /home/dao_ops/jardizzo/FLUID/firex-aq/ops/pyg_modules
cd /home/dao_ops/jardizzo/FLUID/firex-aq/ops

pyfile_iter.py 20190713 0 models.rc models/*.yml
