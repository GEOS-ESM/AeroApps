"""

Collocate Deep Blue and Aeronet.

"""

from datetime import datetime, timedelta

from glob  import glob
from numpy import sort

from pyobs.mapss import *

if __name__ == "__main__":

    # Root directories
    # ----------------
    t_dir = '/nobackup/MAPSS'
    t_dir = '/Users/adasilva/workspace/Data_Analysis/DeepBlue'
    a_dir = t_dir+'/Aeronet_AOD'
    d_dir = t_dir+'/DeepBlue_AOD_Land'
    r_dir = t_dir+'/DeepBlue_Reflectance_Land'
    s_dir = t_dir+'/DeepBlue_Surface_Reflectance_Land'

    one_day = timedelta(hours=24)

    for year in range(2002,2011):
        for month in range(1,13):

            print " "
            print "Processing ", year, month

            # Collect list of matching files for this month
            # ---------------------------------------------
            a_files, d_files, r_files, s_files = ([],[],[],[]) # empty lists
            t0 = datetime(year,month,1)
            for day in range(1,32):
                t = t0 + (day-1) * one_day
                if t.month != month:
                    continue

                a_file = glob(a_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                d_file = glob(d_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                r_file = glob(r_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))
                s_file = glob(s_dir+'/%d/mapss.*.%d.%02d.%02d.*.dat'%(year,year,month,day))

                if (len(a_file)!=len(d_file)) or \
                    (len(a_file)!=len(r_file)) or \
                    (len(a_file)!=len(s_file)):
                    print "- No matching 4-plets for this day ", year, month, day
                    continue

                if len(a_file)==0:
                    print "- No data for this day ", year, month, day
                    continue

                a_files += a_file
                d_files += d_file
                r_files += r_file
                s_files += s_file

            # Read files for this month
            # -------------------------
            if len(a_files)==0:
                 print "- No data for this month ", year, month
                 continue
            else:
                 a, d, r, s = (ANET(a_files),DEEP_AOD(d_files),
                               DEEP_MREF(r_files),DEEP_SREF(s_files),
                               )

            # Make sure d and r files are collocated
            # --------------------------------------
            if any(d.lon!=r.lon):
                raise ValueError, '>>> Inconsistent R and D files for %d %d'%(year,month)
            if any(s.lon!=r.lon):
                raise ValueError, '>>> Inconsistent S and D files for %d %d'%(year,month)

            # Collocate r and a files
            # -----------------------
            print '- Collocating %d reflectances with %d AERONET obs (%d files)'\
                  %(d.nobs,a.nobs,len(a_files))
            Cndx = r.collocate(a)

            # Colapse objects to collocated points only
            # -----------------------------------------
            print '- Colapsing with %d collocations'%Cndx[Cndx>0].size
            I = Cndx>=0
            J = Cndx[I]
            a.colapse(J)
            d.colapse(I)
            r.colapse(I)
            s.colapse(I)

            # Store collocation indices
            # -------------------------
            a.Cndx = J
            d.Cndx = I
            r.Cndx = I
            s.Cndx = I
            
            # Save collocated files
            # ---------------------
            print '- Saving files with %d collocations'%Cndx[Cndx>0].size
            a.savez('mapss.anet.%04d%02d.npz'%(year,month))
            d.savez('mapss.dtau.%04d%02d.npz'%(year,month))
            r.savez('mapss.dref.%04d%02d.npz'%(year,month))
            s.savez('mapss.sref.%04d%02d.npz'%(year,month))


