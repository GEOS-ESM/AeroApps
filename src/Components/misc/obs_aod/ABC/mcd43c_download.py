#!/usr/bin/env python3
"""
Download MCD43C data for a date range
"""

import os, sys, subprocess
from   datetime               import datetime, timedelta
from   dateutil.parser        import parse as isoparse
from   glob                   import glob
import argparse



def downloadFile(args,tyme):

    if args.coll == '006':
        HTTP = 'http://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD43C1.006/'
    elif args.coll == '061':
        HTTP = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD43C1.061/'

    command = 'wget -q -r -nH -nd -np -R "index.html*" -R "*.xml" -P '    

    YY = tyme.strftime('%Y')
    MM = tyme.strftime('%m')
    DD = tyme.strftime('%d')
    doy  = tyme.strftime('%j')

    inFileList = glob("{}/Y{}/M{}/*A{}{}*.hdf".format(args.rootDir,YY,MM,YY,doy))

    if len(inFileList) != 1:
        Outdir = "{}/Y{}/M{}/".format(args.rootDir,YY,MM)
        dd = '{}.{}.{}'.format(YY,MM,DD)
        print('Downloading '+dd)
        cmd = command + Outdir + ' ' + HTTP + dd + '/'
        subprocess.call(cmd,shell=True)        
        inFileList = glob("{}/Y{}/M{}/*A{}{}*.hdf".format(args.rootDir,YY,MM,YY,doy))
        if len(inFileList) != 1:
            print(cmd)
            print('problem downloading '+ dd)
#            raise Exception('problem downloading '+ dd)
            
        
if __name__ == '__main__':
    #   Parse command line options
    #   --------------------------
    rootDir = '/nobackup/3/pcastell/MODIS/MCD43C1/061'
    coll   = '061'

    parser = argparse.ArgumentParser()
    parser.add_argument("startDate",help="Start date of download")

    parser.add_argument("endDate",help="End date of download")

    parser.add_argument('-c','--coll',dest='coll',default=coll,
                        help="MODIS collection version (default=%s)"%coll)

    parser.add_argument("--root",dest="rootDir",default=rootDir,
                        help="where to write data to (default=%s)"%rootDir)

    parser.add_argument("-o","--overwrite",default=False,
                        action="store_true",
                        help="overwrite existing data (default=False)")



    args = parser.parse_args()   

    if not os.path.exists(args.rootDir):
        os.makedirs(args.rootDir)

    sdate = isoparse(args.startDate)
    edate = isoparse(args.endDate)

    while sdate <= edate:
        downloadFile(args,sdate)

        sdate += timedelta(days=1)
