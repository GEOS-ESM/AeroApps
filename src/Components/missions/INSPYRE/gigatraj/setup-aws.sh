#!/bin/bash
module purge
module load netcdf/4.9.3-gcc11.4.0 netcdf-fortran/4.6.2-gcc11.4.0-ompi4.1.7 openmpi/4.1.7  miniconda3/25.5.1-w2aqjdr
conda activate inspyre
export BASEDIR=/dev/null
export ROOT=/shared/users-local/adasilva/local
