#!/bin/csh -fx

#SBATCH --job-name=firex_catalog_$model
#SBATCH --account=s1321
#SBATCH --time=2:00:00
#SBATCH --qos=daohi
#SBATCH --ntasks=28
#SBATCH --export=NONE
#SBATCH --constraint=hasw
#SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_catalog_$model_%iy4%im2%id2_%ih2.log

limit stacksize unlimited
source /home/dao_ops/jardizzo/FLUID/firex-aq/utils/pyg_modules
cd /home/dao_ops/jardizzo/FLUID/firex-aq/utils

firex_catalog.sh %iy4%im2%id2 %ih20000 $model %iy4%im2%id2%ih2 %y4%m2%d2%h2 $tau_inc
