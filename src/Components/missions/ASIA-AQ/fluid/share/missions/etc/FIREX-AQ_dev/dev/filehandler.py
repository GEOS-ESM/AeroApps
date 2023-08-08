import os
import ftplib
import urllib
import datetime as dt
from string import *
from multiprocessing import Pool

def get_handler(protocol):

    handlers = { 'ftp': FTPHandler,
                 'wget': WGETHandler
               }

    return handlers.get(protocol, FileHandler)()

def ftp_upload(request):

    machine     = request['machine']
    login       = request['login']
    password    = request['password']
    remote_file = os.path.basename(request['remote_file'])
    remote_dir  = os.path.dirname(request['remote_file'])
    local_file  = os.path.basename(request['local_file'])
    local_dir   = os.path.dirname(request['local_file'])

    print machine
    print login
    print password
    print remote_file
    print remote_dir
    print local_file
    print local_dir
    print ''

    return

    try:
        os.makedirs(local_dir, 0755)
    except:
        pass

    ftp = ftplib.FTP(machine)
    ftp.login(login,password)
    ftp.cwd(remote_dir)
    f = open(local_file,'wb')
    ftp.retrbinary('RETR ' + remote_file, f.write)
    f.close()
    ftp.quit()

def wget_upload(request):

    url         = request['url']
    remote_file = os.path.join(url,request['remote_file'])
    local_file  = request['local_file']
    local_dir   = os.path.dirname(local_file)

    try:
        os.makedirs(local_dir, 0755)
    except:
        pass

    print 

    urllib.request.urlretrieve(remote_file, local_file)

class FileHandler(object):

    def __init__(self): pass

    def list(self, request):

        req        = dict(request)

        idate      = str(request['idate'])
        itime      = '%06d'%(request['itime'],)
        fcst_dt    = dt.datetime.strptime(idate+itime,'%Y%m%d%H%M%S')
        tau_start  = request['tau_start']
        tau_end    = request['tau_end']
        tau_inc    = request['tau_inc']

        assert tau_inc > 0, "Bad time increment: tau_inc"

        local_dir  = request.get('local_dir',  '')
        remote_dir = request.get('remote_dir', '')

        tau = tau_start
        while tau <= tau_end:

            defs    = {'tau':"%03d"%(tau,)}
            time_dt = fcst_dt + dt.timedelta(hours=tau)

            for name in request['remote_files']:

                local  = os.path.join(local_dir, name)
                remote = os.path.join(remote_dir, name)

                files = remote + ' ' + local

                files = time_dt.strftime(files)
                files = fcst_dt.strftime(files)

                files = Template(files).safe_substitute(defs)
                files = tuple(files.split())

                req['remote_file'] = files[0]
                req['local_file']  = files[1]

                if not os.path.isfile(files[1]): yield dict(req)

            tau += tau_inc

class FTPHandler(FileHandler):

    def get(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(ftp_upload, self.list(request))

class WGETHandler(FileHandler):

    def get(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(wget_upload, self.list(request))
