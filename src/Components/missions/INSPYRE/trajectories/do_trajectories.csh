#!/bin/tcsh
module purge
module load comp/gcc/14.2.0
module load mpi/openmpi/4.1.6/gcc-14.2.0
module load python/GEOSpyD/Min23.5.2-0_py3.11
setenv BASEDIR /discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-8.24.0/x86_64-pc-linux-gnu/gfortran_14.2.0-openmpi_4.1.6/Linux
setenv ROOT "/home/pcolarco/INSPYRE_25/pcolarco/gigatraj/installed/"

./proc_fire_table.py
./initialize_cases.csh
./run_cases.csh
./plot_cases.csh
