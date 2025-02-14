#! /usr/bin/env python

import sys
import ftplib
from multiprocessing import Pool

def upload(x): 

    ftp = ftplib.FTP('dissemination.ecmwf.int')
    ftp.login("xinxin.ye","e3fp20JX")
    ftp.cwd('DATA/CAMS_GLOBAL/2019061000')
    f = open(x,'wb')
    ftp.retrbinary('RETR ' + x, f.write)
    f.close()
    ftp.quit()

if __name__ == '__main__':
    pool = Pool(4)
    pool.map(upload,sys.argv[1:])
