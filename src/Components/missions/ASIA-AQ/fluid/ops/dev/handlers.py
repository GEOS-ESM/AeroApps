import os
import shlex
import ftplib
import subprocess
from shutil import copyfile
from multiprocessing import Pool

from taskmanager import *

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

    print remote_file, local_file

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

def copy_upload(request):

    remote_file = request['remote_file']
    local_file  = request['local_file']
    local_dir   = os.path.dirname(request['local_file'])

    try:
        os.makedirs(local_dir, 0755)
    except:
        pass

    print remote_file, local_file

    if os.path.isfile(remote_file):
        try:
            copyfile(remote_file, local_file)
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

def arl_process(request):

    files       = request['files']
    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    merged_file = os.path.join(output_dir, request['merged_file'])
    prep_file   = os.path.join(output_dir, request['prep_file'])
    vars        = ','.join(request['vars'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    file = os.path.join(input_dir, files[0])
    options = '-h -O -v PBL2 -a LAY'.split()
    subprocess.call(['ncwa'] + options + [file, merged_file])

    file = os.path.join(input_dir, files[1])
    options = ('-h -A -v ' + vars).split()
    subprocess.call(['ncks'] + options + [file, merged_file])

    file = os.path.join(input_dir, files[2])
    options = '-h -A -d LAY,0,19 -v PRES,QV,QC,ZF,TA'.split()
    subprocess.call(['ncks'] + options + [file, merged_file])

def cam_chem_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])
    tmp_file    = output_file + '.tmp.nc'

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    options = '-d lat,15.0,70.0 -d lon,190.0,310.0'.split()
    print ' '.join(['ncks'] + ['-O'] + options + [input_file, tmp_file])
    subprocess.call(['ncks'] + ['-O'] + options + [input_file, tmp_file])

    cmd = 'cam_waccm2p.py -v -o ' + output_file + ' ' + tmp_file
    print cmd
    subprocess.call(cmd.split())

    os.remove(tmp_file)

def arqi_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    ncl_script  = request['ncl_script']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    tfile1  = output_file + '.tmp1.nc'
    tfile2  = output_file + '.tmp2.nc4'

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'ncl' + ' \'file_in="' + input_file + \
              '"\' \'file_out="' + tfile1 + '"\' ' + ncl_script

    print cmd
    subprocess.call(shlex.split(cmd))

    cmd = 'arqi_lv1_2p.py -v -o ' + output_file + ' ' + tfile1
    print cmd
    subprocess.call(cmd.split())

    cmd = 'arqi_lv3_2p.py -v -o ' + tfile2 + ' ' + tfile1
    print cmd
    subprocess.call(cmd.split())

    cmd = 'ncks -h -A ' + tfile2 + ' ' + output_file
    print cmd
    subprocess.call(cmd.split())

    if os.path.isfile(tfile1): os.remove(tfile1)
    if os.path.isfile(tfile2): os.remove(tfile2)

def firework_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    ncl_script  = request['ncl_script']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    tfile1  = output_file + '.tmp1.nc'
    tfile2  = output_file + '.tmp2.nc4'

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'ncl' + ' \'file_in="' + input_file + \
              '"\' \'file_out="' + tfile1 + '"\' ' + ncl_script

    print cmd
    subprocess.call(shlex.split(cmd))

    cmd = 'firework_lv1_2p.py -v -o ' + output_file + ' ' + tfile1
    print cmd
    subprocess.call(cmd.split())

    cmd = 'firework_lv2_2p.py -v -o ' + tfile2 + ' ' + tfile1
    print cmd
    subprocess.call(cmd.split())

    cmd = 'ncks -h -A ' + tfile2 + ' ' + output_file
    print cmd
    subprocess.call(cmd.split())

    if os.path.isfile(tfile1): os.remove(tfile1)
    if os.path.isfile(tfile2): os.remove(tfile2)

def ncarwrfchem_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    ncl_script  = request['ncl_script']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    tfile1  = output_file + '.tmp1.nc'

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'ncl' + ' \'file_in="' + input_file + \
              '"\' \'file_out="' + tfile1 + '"\' ' + ncl_script

    print cmd
    subprocess.call(shlex.split(cmd))

    cmd = 'ncarwrfchem2llp.py -v -o ' + output_file + ' ' + tfile1
    print cmd
    subprocess.call(cmd.split())

    if os.path.isfile(tfile1): os.remove(tfile1)

def hrrr_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    copyfile(input_file, output_file)

def raqms_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'raqms2p.py -v -o ' + output_file + ' ' + input_file
    print cmd
    subprocess.call(cmd.split())

def wrfchem_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'wrf2llp.py -v -o ' + output_file + ' ' + input_file
    print cmd
    subprocess.call(cmd.split())

def uiowawrfchem_process(request):

    input_dir   = request['input_dir']
    output_dir  = request['output_dir']
    input_file  = os.path.join(input_dir,  request['file'])
    output_file = os.path.join(output_dir, request['output_file'])

    try:
        os.makedirs(output_dir, 0755)
    except:
        pass

    cmd = 'uiowawrfchem2llp.py -v -o ' + output_file + ' ' + input_file
    print cmd
    subprocess.call(cmd.split())
