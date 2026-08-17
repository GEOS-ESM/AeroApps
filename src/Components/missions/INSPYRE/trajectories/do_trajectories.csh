#!/bin/tcsh
module purge
module load comp/gcc/14.2.0
module load mpi/openmpi/4.1.6/gcc-14.2.0
module load python/GEOSpyD/Min23.5.2-0_py3.11
setenv BASEDIR /discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-8.24.0/x86_64-pc-linux-gnu/gfortran_14.2.0-openmpi_4.1.6/Linux
setenv ROOT "/home/pcolarco/INSPYRE_25/pcolarco/gigatraj/installed/"

# Define the forecast cycle and the base output path
set FcstCycle = `date +"%Y%m%d"`_00
set BaseOutDir = "/discover/nobackup/projects/gmao/iesa/pub/aerosol/inspyre/trajectories"
set FinalOutDir = "${BaseOutDir}/${FcstCycle}"
mkdir -p ${FinalOutDir}

./proc_fire_table.py
./initialize_cases.csh
./run_cases.csh
./plot_cases.csh

chgrp s3339 *png
chmod g+w *png
#Copy the pngs into the specific forecast cycle folder
\mv -f fcst${FcstCycle}*png ${FinalOutDir}/

#set group permissions on the new directory itself just in case
chgrp s3339 ${FinalOutDir}
chmod g+ws ${FinalOutDir}
