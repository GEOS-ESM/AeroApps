#!/bin/tcsh
#SBATCH --time=5:00:00
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --job-name=cpexcv_sample_regrid
#SBATCH --account=s3339

cd /home/pcolarco/geos_aerosols/pcolarco/AeroApps/src/Components/missions/CPEX-CV/

# setup environment
  setenv GEOSMODEL /home/pcolarco/geos_aerosols/pcolarco/AeroApps/
  source $GEOSMODEL/env@/g5_modules
  set path = ( . $GEOSMODEL/install/bin $path )

  set REGRID_EXEC = $GEOSMODEL/install/bin/Regrid_Util.x
  set GRID_SPEC   = "PC1152x721-DC"

#  set coll = "aer_inst_3hr_glo_C360x360x6_v72"
  set coll = "asm_inst_1hr_glo_C360x360x6_v72"
  
#  set m21cdir = "/home/pcolarco/MERRA21C/e5303_m21c_jan18/chem/Y2022/M09/"
  set m21cdir = "/home/pcolarco/MERRA21C/e5303_m21c_jan18/diag/Y2022/M09/"
  foreach file (`\ls -1 $m21cdir/e5303_m21c_jan18.${coll}.2022*T0[0,3,6,9]* $m21cdir/e5303_m21c_jan18.${coll}.2022*T1[2,5,8]* $m21cdir/e5303_m21c_jan18.${coll}.2022*T21*`)
  set file_ = `echo $file:t`
  set ofile = `echo "$file_:gs/C360x360x6/L1152x721/"`
  echo $ofile
  
            mpirun -np 6 ${REGRID_EXEC} \
                -i "${file}" \
                -o "${ofile}" \
                -ogrid "${GRID_SPEC}" \
                -file_weights 
  end

#  set m21cdir = "/home/pcolarco/MERRA21C/e5303_m21c_jan18/chem/Y2022/M10/"
  set m21cdir = "/home/pcolarco/MERRA21C/e5303_m21c_jan18/diag/Y2022/M10/"
  foreach file (`\ls -1 $m21cdir/e5303_m21c_jan18.${coll}.2022-10-0*T0[0,3,6,9]* $m21cdir/e5303_m21c_jan18.${coll}.2022-10-0*T1[2,5,8]* $m21cdir/e5303_m21c_jan18.${coll}.2022-10-0*T21*`)

  set file_ = `echo $file:t`
  set ofile = `echo "$file_:gs/C360x360x6/L1152x721/"`
  echo $ofile
  
            mpirun -np 6 ${REGRID_EXEC} \
                -i "${file}" \
                -o "${ofile}" \
                -ogrid "${GRID_SPEC}" \
                -file_weights 
  end
