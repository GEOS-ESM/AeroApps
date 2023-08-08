#! /usr/bin/env python

import sys
import ftplib
from multiprocessing import Pool

def upload(file): 

    ftp = ftplib.FTP('dissemination.ecmwf.int')
    ftp.login("xinxin.ye","e3fp20JX")
    f = open(file,'wb')
    ftp.retrbinary('RETR ' + file, f)
    f.close()
    ftp.quit()

print sys.argv[1:]
sys.exit(0)

pool = Pool(3)
pool.map(upload,sys.argv[1:])
