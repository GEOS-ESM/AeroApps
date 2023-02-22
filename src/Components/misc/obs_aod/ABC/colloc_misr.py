"""

Collocate MISR and Aeronet.

"""

import os
from datetime import datetime, timedelta

from glob  import glob
from numpy import sort

from pyobs.mapss import *

if __name__ == "__main__":

    # Root directories
    # ----------------
    t_dir = '/nobackup/MAPSS'
    a_dir = t_dir+'/AERONET/AerosolOpticalDepth'
    g_dir = t_dir+'/MISR/Geometry'
    r_dir = t_dir+'/MISR/RegEqRefl'
    m_dir = t_dir+'/MISR/RegBestEstimateSpectralOptDepth'


    one_day = timedelta(hours=24)

    for year in range(2000,2015):
        for month in range(1,13):

            print(" ")
            print("Processing ", year, month)

            if os.path.isfile('mapss.misr_maod.%04d%02d.npz'%(year,month)):
                print("  ... already processed, skipping")
                continue


            # Collect list of matching files for this month
            # ---------------------------------------------
            a_files, g_files, r_files, m_files = ([],[],[],[]) # empty lists
            t0 = datetime(year,month,1)
            for day in range(1,32):
                t = t0 + (day-1) * one_day
                if t.month != month:
                    continue

                a_file = glob(a_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                g_file = glob(g_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                r_file = glob(r_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                m_file = glob(m_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))


                if  (len(a_file)!=len(g_file)) or \
                    (len(a_file)!=len(r_file)) or \
                    (len(a_file)!=len(m_file)):
                    print("- No matching 4-plets for this day ", year, month, day)
                    continue

                if len(a_file)==0:
                    print("- No data for this day ", year, month, day)
                    continue

                a_files += a_file
                g_files += g_file
                r_files += r_file
                m_files += m_file

            # Read files for this month
            # -------------------------
            if len(a_files)==0:
                 print("- No data for this month ", year, month)
                 continue
            else:
                 a, g, r, m = (ANET(a_files),MISR_GEOM(g_files),
                               MISR_MREF(r_files),MISR_AOD(m_files),
                               )

            # Make sure d and r files are collocated
            # --------------------------------------
            if any(g.lon!=m.lon):
                raise ValueError('>>> Inconsistent D and M files for %d %d'%(year,month))
            if any(g.lon!=r.lon):
                raise ValueError('>>> Inconsistent G and R files for %d %d'%(year,month))

            # Collocate r and a files
            # -----------------------
            print('- Collocating %d reflectances with %d AERONET obs (%d files)'\
                  %(r.nobs,a.nobs,len(a_files)))
            Cndx = r.collocate(a)

            # Colapse objects to collocated points only
            # -----------------------------------------
            print('- Colapsing with %d collocations'%Cndx[Cndx>0].size)
            I = Cndx>=0
            J = Cndx[I]
            a.colapse(J)
            g.colapse(I)
            r.colapse(I)
            m.colapse(I)

            # Store collocation indices
            # -------------------------
            a.Cndx = J
            g.Cndx = I
            r.Cndx = I
            m.Cndx = I
            
            # Save collocated files
            # ---------------------
            print('- Saving files with %d collocations'%Cndx[Cndx>0].size)
            a.savez('mapss.anet_misr.%04d%02d.npz'%(year,month))
            g.savez('mapss.misr_geom.%04d%02d.npz'%(year,month))
            r.savez('mapss.misr_mref.%04d%02d.npz'%(year,month))
            m.savez('mapss.misr_maod.%04d%02d.npz'%(year,month))


