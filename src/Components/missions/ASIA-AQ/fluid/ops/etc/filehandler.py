import os
import sys
import ftplib
import subprocess
from string import *
from shutil import copyfile
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

import gradstime as gt
import gradsdataservice as dataservice

from handlers import *
from taskmanager import *

def get_handler(protocol, config):

    handlers = { 'ftp': FTPHandler,
                 'wget': WGETHandler,
                 'copy': CopyHandler,
                 'cams': CAMSHandler,
                 'cam_chem': CAMChemHandler,
                 'arqi': ARQIhandler,
                 'hrrr-smoke': HRRRhandler,
                 'raqms': RAQMShandler,
                 'wrf-chem': WRFChemHandler,
                 'geosfp': GEOSFPhandler,
                 'geoscf': GEOSCFhandler
               }

    return handlers.get(protocol, NoHandler)(config)

class FileHandler(object):

    def __init__(self, config=None):

        self.gtime = None
        self.tau   = 0

        self.defs  = {}
        if config: self.defs = dict(config)

    def iter(self, request):

        req = dict(request)

        local_dir  = request.get('local_dir',  '')
        remote_dir = request.get('remote_dir', '')

        for tau in self.itertime(request):

            for name in request['remote_files']:

                local  = self.resolve(os.path.join(local_dir, name))
                remote = self.resolve(os.path.join(remote_dir, name))

                req['remote_file'] = remote
                req['local_file']  = local

                if not os.path.isfile(local):
                    yield dict(req)
                elif os.stat(local).st_size == 0:
                    yield dict(req)

    def status(self, request):

        status = 0
        for x in self.iter(request): status += 1

        return status

    def resolve(self, s, time_handler=None):

        if not time_handler: time_handler = self.gtime.strftime

        s_out = Template(s).safe_substitute(self.defs)
        s_out = Template(s_out).safe_substitute(self.defs)
        s_out = time_handler(s_out, self.tau)
        
        return s_out

    def itertime(self, request):

        self.defs.update(request)

        idate      = request['idate']
        itime      = request['itime']
        shift_dt   = int(request.get('shift_dt', 0))
        tau_start  = request['tau_start']
        tau_end    = request['tau_end']
        tau_inc    = request['tau_inc']
        self.gtime = gt.GradsTime(idate,itime,shift_dt)

        assert tau_inc > 0, "Bad time increment: tau_inc"

        tau = tau_start
        while tau <= tau_end:
            self.tau = tau
            yield tau
            tau += tau_inc

    def grads_write(self, request):

        grads_dir     = self.resolve(request['grads_dir'])
        grads_file    = self.resolve(request['grads_file'])
        grads_file    = os.path.join(grads_dir, grads_file)
        dummy_file    = self.resolve(request['dummy_file'])
        grads_records = request.get('grads_records', [])

        try:
            os.makedirs(grads_dir, 0755)
        except:
            pass

        with open(grads_file, 'w') as f:

            for rec in grads_records:
                f.write(self.resolve(rec, self.gtime.stritime) + '\n')

#       ds = dataservice.Service()
#       try:
#           ret = ds.open(grads_file)
#       except:
#           copyfile(dummy_file, grads_file)

class NoHandler(FileHandler):

    def handler(self, request): pass

    def status(self, request): return 1

class FTPHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(ftp_upload, self.iter(request))

        return self.status(request)

class CopyHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(copy_upload, self.iter(request))

        return self.status(request)

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

        return self.status(request)

class CAMSHandler(FileHandler):

    def iter(self, request, **kwargs):

        req = dict(request)
        req.update(kwargs)

        force = req.get('force', False)

        for tau in self.itertime(request):

            req['merged_file'] = self.resolve(request['merged_file'])
            req['prep_file']   = self.resolve(request['prep_file'])
            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])

            prep_file = os.path.join(req['output_dir'],req['prep_file'])

            if not force and os.path.isfile(prep_file): continue
            print prep_file

            req['files'] = []
            for file in request['files']:
                req['files'].append(self.resolve(file))

            yield dict(req)

    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(cams_process, self.iter(request))

        self.grads_write(request)

        return 0

class CAMChemHandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        for tau in self.itertime(request):

            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])
            req['output_file'] = self.resolve(request['output_file'])

            output_file = os.path.join(req['output_dir'],req['output_file'])

            if os.path.isfile(output_file): continue

            req['file'] = self.resolve(request['files'][0])

            yield dict(req)


    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(cam_chem_process, self.iter(request))

        self.grads_write(request)

        return 0

class ARQIhandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        for tau in self.itertime(request):

            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])
            req['output_file'] = self.resolve(request['output_file'])

            output_file = os.path.join(req['output_dir'],req['output_file'])

            if os.path.isfile(output_file): continue

            req['file'] = self.resolve(request['files'][0])

            yield dict(req)


    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(arqi_process, self.iter(request))

        self.grads_write(request)

        return 0

class HRRRhandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        for tau in self.itertime(request):

            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])
            req['output_file'] = self.resolve(request['output_file'])

            output_file = os.path.join(req['output_dir'],req['output_file'])

            if os.path.isfile(output_file): continue

            req['file'] = self.resolve(request['files'][0])

            yield dict(req)


    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(hrrr_process, self.iter(request))

        self.grads_write(request)

        return 0

    def grads_write(self, request):

        file          = self.resolve(request['files'][0])
        file_dir      = self.resolve(request['output_dir'])
        input_file    = os.path.join(file_dir, file)
        grads_dir     = self.resolve(request['grads_dir'])
        grads_file    = self.resolve(request['grads_file'])
        grads_tmp     = os.path.join(grads_dir, grads_file) + '.tmp'
        grads_ctl     = os.path.join(grads_dir, grads_file) + '.ctl'
        grads_index   = os.path.join(grads_dir, grads_file) + '.idx'
        grads_records = request.get('grads_records', [])

        try:
            os.makedirs(grads_dir, 0755)
        except:
            pass

        f = open(grads_tmp, "w")
        cmd = 'alt_g2ctl ' + input_file
        subprocess.check_call(cmd.split(), stdout=f)
        f.close()

        lines = []
        with open(grads_tmp, 'r') as f: lines = f.readlines()

        with open(grads_ctl, 'w') as f_out:

            for rec in lines:

                key = rec.split()[0]

                if key in grads_records:
                    rec = key + ' ' + grads_records[key] 
                    f_out.write(self.resolve(rec, self.gtime.stritime) + '\n')
                    if key == 'dset': f_out.write('options template\n')
                else:
                    f_out.write(rec)

        os.remove(grads_tmp)

        cmd = 'alt_gmp ' + grads_ctl
        subprocess.call(cmd.split())

class RAQMShandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        for tau in self.itertime(request):

            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])
            req['output_file'] = self.resolve(request['output_file'])

            output_file = os.path.join(req['output_dir'],req['output_file'])

            if os.path.isfile(output_file): continue

            req['file'] = self.resolve(request['files'][0])

            yield dict(req)

    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        pool = Pool(tasks)
        pool.map(raqms_process, self.iter(request))

        self.grads_write(request)

        return 0

class WRFChemHandler(FileHandler):

    def iter(self, request):

        req = dict(request)

        for tau in self.itertime(request):

            req['input_dir']   = self.resolve(request['input_dir'])
            req['output_dir']  = self.resolve(request['output_dir'])
            req['output_file'] = self.resolve(request['output_file'])

            output_file = os.path.join(req['output_dir'],req['output_file'])

            if os.path.isfile(output_file): continue

            req['file'] = self.resolve(request['files'][0])

            print 'iterating'
            yield dict(req)

    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

#       tp = ThreadPool(tasks)
#       for r in self.iter(request):
#           tp.apply_async(wrfchem_process, (r,))

#       tp.close()
#       tp.join()

        pool = Pool(tasks)
        pool.map(wrfchem_process, self.iter(request))

        self.grads_write(request)

        return 0

class GEOSFPhandler(FileHandler):

    def handler(self, request):

        tau = next(self.itertime(request))

        for collection in request['collections']:

            data_dir    = self.resolve(request['data_dir'])
            input_dir   = self.resolve(request['input_dir'])
            output_dir  = self.resolve(request['output_dir'])
            dummy_file  = self.resolve(request['dummy_file'])
            file        = self.resolve(collection + '.%iy4%im2%id2_%ih2')

            input_file  = os.path.join(input_dir, collection, file)
            output_file = os.path.join(output_dir, file)

            try:
                os.makedirs(output_dir, 0755)
            except:
                pass

            if not os.path.isfile(input_file):
                copyfile(dummy_file, output_file)
                return

            lines = []
            with open(input_file, 'r') as f: lines = f.readlines()

            with open(output_file, 'w') as f_out:

                for rec in lines:

                    key = rec.split()[0]

                    if key.lower() == 'dset':
                        name = os.path.basename(rec.split()[-1])
                        pathname = os.path.join(data_dir, name)
                        f_out.write('dset ' + pathname + '\n')
                    else:
                        f_out.write(rec)

        return 0

class GEOSCFhandler(FileHandler):

    def handler(self, request):

        tau = next(self.itertime(request))

        for collection in request['collections']:

            data_dir    = self.resolve(request['data_dir'])
            input_dir   = self.resolve(request['input_dir'])
            output_dir  = self.resolve(request['output_dir'])
            dummy_file  = self.resolve(request['dummy_file'])
            file        = self.resolve(collection + '.%iy4%im2%id2_%ih2z')

            input_file  = os.path.join(input_dir, collection, file)
            output_file = os.path.join(output_dir, file)

            try:
                os.makedirs(output_dir, 0755)
            except:
                pass

            if not os.path.isfile(input_file):
                copyfile(dummy_file, output_file)
                print 'Missing collection: ',collection
                return 0

            lines = []
            with open(input_file, 'r') as f: lines = f.readlines()

            with open(output_file, 'w') as f_out:

                for rec in lines:

                    key = rec.split()[0]

                    if key.lower() == 'dset':
                        name = os.path.basename(rec.split()[-1])
                        pathname = os.path.join(data_dir, name)
                        f_out.write('dset ' + pathname + '\n')
                    elif key.lower() == 'tdef':
                        tinc = rec.split()[-1]
                        f_out.write(' '.join(rec.split()[0:-1]) + ' 1hr' + '\n')
                    else:
                        f_out.write(rec)

        return 0
