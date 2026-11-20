#!/bin/csh

# Set base location
  cd /home/pcolarco/ARCSIX/magpie

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

set model = "c180R_arcsix.inst3d_aer_v"
set model = "fp.inst3_3d_aer_Nv"
set model = "MERRA2.inst3_3d_aer_Nv"
set model = "MERRA2.tavg1_2d_aer_Nx"
set dir = "stn_samples"
  
# 
  foreach MM ("05" "06" "07" "08")
  foreach DD ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" \
              "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" \
              "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")

# Fix time to be record dimension
  set file = $dir/$model.stations.2024$MM$DD.nc
  if( -f $file) then
    echo $file
    ncks -O --mk_rec_dmn time $file -o $file >> /dev/null
    set attstr = "minutes since 2024-"$MM"-"$DD" 00:00:00"
    set unitstr = \"`echo $attstr`\"
    echo ncatted -O -a units,time,o,c,$unitstr $file > att.txt
    chmod 755 att.txt
    ./att.txt >> /dev/null
    rm -f att.txt
  endif

  end
  end

# And finally ncrcat
  ncrcat $dir/$model.stations.2024????.nc $dir/$model.stations.2024.nc
  
