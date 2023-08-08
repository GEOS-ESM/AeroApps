import os
import tempfile
import sys
import ftplib
import subprocess
from string import *
from shutil import copyfile
from multiprocessing import Pool

import gradstime as gt
import gradsdataservice as dataservice

from handlers import *
from taskmanager import *

def get_handler(protocol, config, status=0):

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
                 'geoscf': GEOSCFhandler,
                 'uclawrf': UCLAWRFhandler,
                 'uiowawrf': UIOWAWRFChemHandler,
                 'firework': FireWorkHandler,
                 'ncarwrfchem': NCARWRFChemHandler,
               }

    return handlers.get(protocol, NoHandler)(config, status)

class FileHandler(object):

    def __init__(self, config=None, status=0):

        self.gtime  = None
        self.tau    = 0
        self.status = status

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

    def get_status(self, request):

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

        idate       = request['idate']
        itime       = request['itime']
        shift_dt    = int(request.get('shift_dt', 0))
        tau_start   = request['tau_start']
        tau_end     = request['tau_end']
        tau_inc     = request['tau_inc']
        self.gtime  = gt.GradsTime(idate,itime,shift_dt)
        self.gstart = gt.GradsTime(idate,itime,shift_dt+tau_start)

        self.defs['gstart'] = self.gstart.stritime('%ih2:00z%id2%imc%iy4')

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
        grads_records = request.get('grads_records', [])

        try:
            os.makedirs(grads_dir, 0755)
        except:
            pass

        if self.status != 0:
            bname = os.path.basename(grads_file)
            dummy_file = bname.split('.')[0]
            dummy_file = os.path.join(grads_dir, dummy_file)
            copyfile(dummy_file, grads_file)
            return

        with open(grads_file, 'w') as f:

            for rec in grads_records:
                f.write(self.resolve(rec, self.gtime.stritime) + '\n')

class PlotHandler(FileHandler):

    def handler(self, request):

        idate      = request['idate']
        itime      = request['itime']
        shift_dt   = int(request.get('shift_dt', 0))
        tau_start  = request['tau_start']
        tau_end    = request['tau_end']
        tau_inc    = request['tau_inc']
        gtime      = gt.GradsTime(idate,itime,shift_dt)

        batch_job  = request['BATCH']

        ftmp = tempfile.NamedTemporaryFile()
        with open(ftmp.name, 'w') as f:
            for line in batch_job:

                line = Template(line).safe_substitute(self.defs)
                line = Template(line).safe_substitute(request)
                line  = gtime.strftime(line, tau_end)
                f.write(line+'\n')

        os.system('cp ' + ftmp.name + ' temp.j')
        os.system('qsub ' + ftmp.name)

class NoHandler(FileHandler):

    def handler(self, request): pass

    def get_status(self, request): return 1

class FTPHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(ftp_upload, self.iter(request))

        return self.get_status(request)

class CopyHandler(FileHandler):

    def handler(self, request):

        sockets = request.get('sockets', 4)
        assert sockets >= 1 and sockets <=10, "Bad socket request"

        pool = Pool(sockets)
        pool.map(copy_upload, self.iter(request))

        return self.get_status(request)

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

        return self.get_status(request)

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(cams_process, self.iter(request))
        else:
            for tau in self.iter(request): break

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(cam_chem_process, self.iter(request))
        else:
            for tau in self.iter(request): break

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(arqi_process, self.iter(request))
        else:
            for r in self.iter(request): break

        self.grads_write(request)

        return 0

class FireWorkHandler(FileHandler):

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(firework_process, self.iter(request))
        else:
            for r in self.iter(request): break

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(hrrr_process, self.iter(request))
        else:
            for tau in self.iter(request): break

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

        if self.status != 0:
            bname = os.path.basename(grads_ctl)
            dummy_file = bname.split('.')[0]
            dummy_file = os.path.join(grads_dir, dummy_file)
            copyfile(dummy_file, grads_ctl)
            return

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(raqms_process, self.iter(request))
        else:
            for tau in self.iter(request): break

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

            yield dict(req)

    def handler(self, request):

        tasks = request.get('tasks', 4)
        assert tasks >= 1 and tasks <=10, "Maximum tasks exceeded"

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(wrfchem_process, self.iter(request))
        else:
            for tau in self.iter(request): break

        self.grads_write(request)

        return 0

class UIOWAWRFChemHandler(FileHandler):

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(uiowawrfchem_process, self.iter(request))
        else:
            for tau in self.iter(request): break

        self.grads_write(request)

        return 0

class NCARWRFChemHandler(FileHandler):

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

        if self.status == 0:
            pool = Pool(tasks)
            pool.map(ncarwrfchem_process, self.iter(request))
        else:
            for tau in self.iter(request): break

        self.grads_write(request)

        return 0

class GEOSFPhandler(FileHandler):

    def handler(self, request):

        tau = next(self.itertime(request))

        for collection in request['collections']:

            data_dir    = self.resolve(request['data_dir'])
            input_dir   = self.resolve(request['input_dir'])
            output_dir  = self.resolve(request['output_dir'])
            file        = self.resolve(collection + '.%iy4%im2%id2_%ih2')

            input_file  = os.path.join(input_dir, collection, file)
            output_file = os.path.join(output_dir, file)

            try:
                os.makedirs(output_dir, 0755)
            except:
                pass

            dummy_file = os.path.join(output_dir, collection)

            if not os.path.isfile(input_file):
                copyfile(dummy_file, output_file)
                continue

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
            file        = self.resolve(collection + '.%iy4%im2%id2_%ih2z')

            input_file  = os.path.join(input_dir, collection, file)
            output_file = os.path.join(output_dir, file)

            try:
                os.makedirs(output_dir, 0755)
            except:
                pass

            dummy_file = os.path.join(output_dir, collection)

            if not os.path.isfile(input_file):
                copyfile(dummy_file, output_file)
                continue

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

class UCLAWRFhandler(FileHandler):

    def handler(self, request):

        tau = next(self.itertime(request))

        self.grads_write(request)

        return 0

