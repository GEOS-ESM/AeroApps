import os
import ftplib
import subprocess
import datetime as dt
from string import *
from multiprocessing import Pool

from taskmanager import *

def get_handler(protocol):

    handlers = { 'ftp': FTPHandler,
                 'wget': WGETHandler,
                 'cams': CAMSHandler
               }

    return handlers.get(protocol, FileHandler)()

def ftp_upload(request):

    machine     = request['machine']
    login       = request['login']
    password    = request['password']
    remote_file = os.path.basename(request['remote_file'])
    remote_dir  = os.path.dirname(request['remote_file'])
    local_file  = request['local_file']
    local_dir   = os.path.dirname(request['local_file'])

    try:
        os.makedirs(local_dir, 0755)
    except:
        pass

    try:
        ftp = ftplib.FTP(machine)
        ftp.login(login,password)
        ftp.cwd(remote_dir)
        f = open(local_file,'wb')
        ftp.retrbinary('RETR ' + remote_file, f.write)
        f.close()
        ftp.quit()
    except:
        pass

def cams_process(request):

    files       = request['files']
    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    merged_file = os.path.join(output_dir, request['merged_file'])
    prep_file   = os.path.join(output_dir, request['prep_file'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    options = '-d latitude,15.0,70.0 -d longitude,190.0,310.0'.split()

    file = os.path.join(input_dir, files[0])
    subprocess.call(['ncks'] + ['-O'] + options + [file, merged_file])

    for file in files[1:]:
        file = os.path.join(input_dir, file)
        subprocess.call(['ncks', '-h', '-A'] + options + [file, merged_file])

    cmd = 'cams.py -v -o ' + prep_file + ' ' + merged_file
    subprocess.call(cmd.split())

class FileHandler(object):

    def __init__(self):

        self.fcst_dt = None
        self.time_dt = None
        self.defs    = None    

    def iter(self, request):

        req = dict(request)

        idate        = str(request['idate'])
        itime        = '%06d'%(request['itime'],)
        self.fcst_dt = dt.datetime.strptime(idate+itime,'%Y%m%d%H%M%S')

        local_dir  = request.get('local_dir',  '')
        remote_dir = request.get('remote_dir', '')

        tau_start  = request['tau_start']
        tau_end    = request['tau_end']
        tau_inc    = request['tau_inc']

        assert tau_inc > 0, "Bad time increment: tau_inc"

        tau = tau_start
        while tau <= tau_end:

            self.defs    = { 'tau':"%03d"%(tau,), 'tau2':"%02d"%(tau,) }
            self.time_dt = self.fcst_dt + dt.timedelta(hours=tau)

            for name in request['remote_files']:

                local  = self.resolve(os.path.join(local_dir, name))
                remote = self.resolve(os.path.join(remote_dir, name))

                req['remote_file'] = remote
                req['local_file']  = local

                if not os.path.isfile(local): yield dict(req)

            tau += tau_inc

    def resolve(self, s):

        s_out = self.time_dt.strftime(s)
        s_out = self.fcst_dt.strftime(s_out)
        s_out = Template(s_out).safe_substitute(self.defs)
        
        return s_out

class FTPHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(ftp_upload, self.iter(request))

class WGETHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        task = TaskManager(sockets)

        for r in self.iter(request):

            url         = r['url']
            remote_file = os.path.join(url,r['remote_file'])
            local_file  = r['local_file']
            local_dir   = os.path.dirname(local_file)

            try:
                os.makedirs(local_dir, 0755)
            except:
                pass

            task.spawn('wget -O ' + local_file + ' ' + remote_file)

        task.wait()

class CAMSHandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        idate      = str(request['idate'])
        itime      = '%06d'%(request['itime'],)
        fcst_dt    = dt.datetime.strptime(idate+itime,'%Y%m%d%H%M%S')
        tau_start  = request['tau_start']
        tau_end    = request['tau_end']
        tau_inc    = request['tau_inc']

        assert tau_inc > 0, "Bad time increment: tau_inc"

        tau = tau_start
        while tau <= tau_end:

            time_dt       = fcst_dt + dt.timedelta(hours=tau)
            self.fcst_dt  = fcst_dt
            self.time_dt  = time_dt
            self.defs     = { 'tau':"%03d"%(tau,), 'tau2':"%02d"%(tau,) }

            req['merged_file'] = self.resolve(request['merged_file'])
            req['prep_file']   = self.resolve(request['prep_file'])
            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])

            prep_file = os.path.join(req['output_dir'],req['prep_file'])

            if os.path.isfile(prep_file): continue

            req['files'] = []
            for file in request['files']:
                req['files'].append(self.resolve(file))

            yield dict(req)

            tau += tau_inc

    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(cams_process, self.iter(request))
