#!/usr/bin/env python3
"""
Splits AVHRR into synoptic chunks.
"""

import os
import sys
from pyobs.avhrr import AVHRR_L2B, getSyn, granules

from glob      import glob
from datetime  import date, datetime, timedelta
from numpy     import savez, array

#---------------------------------------------------------------------
def makethis_dir(filename):
    """Creates the relevant directory if necessary."""
    path, filen = os.path.split(filename)
    if path != '':
        rc = os.system('mkdir -p '+path)
        if rc:
            raise IOError("could not create directory "+path)
        
#---------------------------------------------------------------------
if __name__ == "__main__":

    RootDir = '/nobackup/AVHRR/Level2/PATMOSX'
    NpzDir  = '/nobackup/AVHRR/Level2/Synoptic'

    ident = 'patmosx_v05r02'
    orb = 'des'
    
    oneday = timedelta(seconds=24*60*60)
    Verb = False

    if len(sys.argv)<2:
        print("Usage:")
        print("        synoptic_avhrr.py  year1 [year2 year3 ...]")
        sys.exit(1)
    else:
        Years = array(sys.argv[1:]).astype(int)

    # Loop over years
    # ---------------
    for year in Years:
        tyme0 = datetime(year,1,1)
        toy = datetime(year,12,31) - tyme0
        ndays = 1 + int(toy.total_seconds()/(24*60*60))
        for doy in range(1,ndays+1):

            tyme = tyme0 + (doy-1) * oneday

            nymd = tyme.year*10000 + tyme.month*100  + tyme.day

            fname = NpzDir+'/Y%d/M%02d/D%02d/%s.%s.%d%02d%02d_21z.npz'%\
                           (tyme.year,tyme.month,tyme.day,ident,orb,tyme.year,tyme.month,tyme.day)

            if os.path.exists(fname):
                print("<> Skipping ", tyme.date())
                continue
            else:
                print("<> Processing ", tyme.date(), doy)

            # Read the bracketing data for this day
            # -------------------------------------
            Files = granules('n??',orb,tyme,RootDir=RootDir,bracket=True)
            if len(Files) > 0:
                a = AVHRR_L2B(Files,Verb=Verb)
            else:
                for syn_hour in range(0,24,3):
                    t = tyme + timedelta(seconds=syn_hour*60*60)
                    print("   x Writing placeholder file on", t)
                    fname = NpzDir+'/Y%d/M%02d/D%02d/%s.%s.%d%02d%02d_%02dz.npz'%\
                            (t.year,t.month,t.day,ident,orb,t.year,t.month,t.day,t.hour)
                    makethis_dir(fname)
                    savez(fname,nobs=0)
                continue
                
            # Sanity check
            # ------------
            iValid = (a.tau_630  > -0.01) &\
                     (a.ref_630  >  0)    &\
                     (a.ref_860  >  0)    
#                    (a.tau_860  > -0.01) &\
#                    (a.ref_1600 >  0) 

            dt_secs = timedelta(seconds=3*60*60) # 3 hours

            # Loop over synoptic times
            # ------------------------
            Nobs = 0
            for syn_hour in range(0,24,3):

                t = tyme + timedelta(seconds=syn_hour*60*60)
                t1 = t - dt_secs/2
                t2 = t + dt_secs/2

                # Down selection
                # --------------
                I = iValid & (a.tyme>=t1) & (a.tyme<t2)

                nobs  = len(a.lon[I])
                Nobs += nobs

                # Save NPZ file
                # -------------
                fname = NpzDir+'/Y%d/M%02d/D%02d/%s.%s.%d%02d%02d_%02dz.npz'%\
                               (t.year,t.month,t.day,ident,orb,t.year,t.month,t.day,t.hour)

                if os.path.exists(fname):
                    print("   - Skipping writing %d obs on "%nobs,t)
                    continue

                makethis_dir(fname)
                if any(I):
                    print("   + Writing %d obs on "%nobs,t)
                    a.writeNPZ(fname,I=I)
                else:
                    print("   o Writing empty file on",t)
                    savez(fname,nobs=0)
