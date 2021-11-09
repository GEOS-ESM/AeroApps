#!/bin/csh
# Script to package vlidort OMI simulator output into a copy
# of OMAERUV Level2 files so that the same readers can be
# used. Adds some additional variables output from the 
# simulator

# Set up environment
  set EXPID = c180R_v202_aura_schill
  cd $NOBACKUP/aura/${EXPID}

# OMAERUV Level 2 file location
  set inp = '/home/pcolarco/piesa/data/OMI/Level2/OMAERUV/'

# Create local output directory
  mkdir -p OMAERUV

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# Copy OMI L2 file locally and stuff it with model produced radiances
  foreach file (`\ls -1 ai.*Full.npz`)
   set datestr = `echo $file:r | rev | cut -c5-18 | rev`
   set yyyy = `echo $datestr | cut -c1-4`
   set mm   = `echo $datestr | cut -c6-7`
   set dd   = `echo $datestr | cut -c8-9`
   set inpfil = `\ls -1 $inp/$yyyy/$mm/$dd/OMI-Aura_L2-OMAERUV_${datestr}*he5`
   set file_ = $EXPID.$inpfil:t:r.vl_rad.geos5_pressure.he5
   /bin/cp -f $inpfil OMAERUV/$file_
   set rfile = radiance.`echo $file | cut -c1-3 --complement`

# notes for python processing (ipython --pylab)
# This gets the "official" OMI file, opens it and reads the "NormRadiance" and then fills it
# with the values fromt the VLIDORT simulation
# PRC: add on 5/28/15 the "UVAerosolIndex"
cat > tmp.py <<EOF
import h5py
import numpy as np
# Open the file to modify and get the radiances
f = h5py.File('OMAERUV/$file_',mode='r+')
g = f.get('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields')
v = g.get('NormRadiance')
# Open the model return and stuff the radiances
data = np.load('$rfile')
x = data['rad_VL'].shape
x = x[0]/60
v[:,:,:] = np.reshape(data['rad_VL'],(x,60,3))

# Now get and stuff the "uncorrected" AI (called Residue)
data = np.load('$file')
g = f.get('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields')
v = g.get('Residue')
v[:,:] = np.reshape(data['residue'],(x,60))

# Now get and stuff the "uncorrected" LER at 388 nm only
g = f.get('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields')
v = g.get('Reflectivity')
v[:,:,1] = np.reshape(data['refl_VL'],(x,60))

# Now add some new variables
sp = np.reshape(data['spher'],(x,60))
g = f.create_dataset('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SphericalAlbedo',data=sp)
tr = np.reshape(data['trans'],(x,60))
g = f.create_dataset('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Transmissivity',data=tr)
ic = np.reshape(data['i388calc'],(x,60))
g = f.create_dataset('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Radiance388Calc',data=ic)
aod = np.reshape(data['AOT'],(x,60,3))
g = f.create_dataset('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MERRAero_AOT',data=aod)
ssa = np.reshape(data['SSA'],(x,60,3))
g = f.create_dataset('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MERRAero_SSA',data=ssa)

# Now get the "corrected" AI and stuff into UVAerosolIndex
g = f.get('HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields')
v = g.get('UVAerosolIndex')
v[:,:] = np.reshape(data['AI_VL'],(x,60))

## Now stuff the surface pressure into the TerrainPressure
#data3 = np.load('./inst3d_aer_v/${EXPID}_aer_L2-OMAERUV_${datestr}Full.npz')
#g = f.get('HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields')
#v = g.get('TerrainPressure')
#v[:,:] = np.reshape(data3['ps']/100,(x,60))

# Now close all
f.flush()
f.close()
EOF

    echo $file_
    python ./tmp.py
    if ($status) rm -f OMAERUV/$file_
    echo $status

end


