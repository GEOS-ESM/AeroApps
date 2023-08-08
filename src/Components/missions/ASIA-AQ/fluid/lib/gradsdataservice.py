import gdsvml
from dataservice import *
#from grads.ganum import GaNum
from mygrads.ganum import GaNum

class Service(GaNum,DataService):

    def __init__(self, config=None, window=False):

        self.files   = []
        self.fileID  = 1
        self.cmd_stack = []
        self.grid = {1:0}
        self.dfile = 1
        self.lonvals = None
        self.mpvals = None

        DataService.__init__(self, config)

        try:
            GaNum.__init__(self, Bin='/opt/opengrads/grads --with-version 2.1.0.oga.1',Window=window, Echo=False)
        except:
            try:
                GaNum.__init__(self, Bin='/opt/opengrads-2.0.2/grads --with-version 2.0.2.oga.2',Window=window, Echo=False)
            except:
                GaNum.__init__(self, Bin='grads',Window=window, Echo=False)

    def open(self, request, **kwargs):

        if isinstance(request,basestring):
            return super(Service,self).open(request)

        hash = dict(request)
        hash.update(kwargs)
        
        stream     = hash['stream']
        collection = hash['collection']
        file       = self.inquire(hash)

        for fh in self.files:

            if fh.file != file: continue
            self.fileID = fh.fileinfo.fid
            return fh

        #print "**OPENING** ",file
        fh   = super(Service,self).open(file)

        self.fileID = fh.fid
        self.grid[fh.fid] = self.config(['stream', stream, 'grid'], 0)
        self.cmd("set dfile %d"%fh.fid)
        self.cmd("set t 1")
        self.cmd("set z 1")

        qf = self.query("file")
        qd = self.query("dims")
        qt = self.query("time")
        qc = self.query("ctlinfo")

        if qc.ztype == 'linear':
            qc.zlevs = [ qc.z0 + i * qc.dz for i in range(0,qc.nz) ]

        fh = FileHandle(file=file,
                        fileinfo=qf,
                        diminfo=qd,
                        timeinfo=qt,
                        ctlinfo=qc,
                        stream=stream,
                        collection=collection)

        self.files.append(fh)

        return fh

    def handle(self, id=None):

        if id is None:
            index = self.fileID - 1
        else:
            index = id - 1

        if len(self.files) >= index and index >= 0:
            return self.files[index]

    def cmd2(self, command, **kwargs):

        for cmd in command.split('\n'):

            cmd = self.check_cmd(cmd)
            self.cmd_stack.append(cmd)
            super(Service, self).cmd(cmd, **kwargs)

    def check_cmd(self, cmd):

        grid = self.grid[self.dfile]

        if 'set lon' in cmd:

            self.lonvals = [float(v) for v in cmd.split()[2:]]
            lons = self.check_lon(self.lonvals)
            cmd  = 'set lon'
            for lon in lons: cmd += ' ' + str(lon)
            return cmd

        if 'set mpvals' in cmd:

            self.mpvals = [float(v) for v in cmd.split()[2:]]
            lons   = self.check_lon(self.mpvals[0:2])
            cmd    = 'set mpvals %f %f %f %f'%tuple(lons+self.mpvals[2:])
            return cmd

        if 'set dfile' in cmd:

            self.dfile = int(cmd.split()[2])

            cmd = [cmd]

            if self.lonvals:
                lons = self.check_lon(self.lonvals)
                c = 'set lon'
                for lon in lons: c += ' ' + str(lon)
                cmd += [c]

            if self.mpvals:
                lons = self.check_lon(self.mpvals[0:2])
                cmd += ['set mpvals %f %f %f %f'%tuple(lons+self.mpvals[2:])]

            return '\n'.join(cmd)

        return cmd

    def check_lon(self, lons):

        grid = self.grid[self.dfile]
        if grid == 0: return lons

        lons = list(lons)

        if grid == -1:
            for i,lon in enumerate(lons):
                if lon > 180.0: lons[i] -= 360.0
            if len(lons) == 2 and lons[0] == lons[1]: lons = [-180, 180]

        elif grid == 1:
            for i,lon in enumerate(lons):
                if lon < 0.0: lons[i] += 360.0
            if len(lons) == 2 and lons[0] == lons[1]: lons = [0, 360]

        return lons
        
    def qstack(self):

        return list(self.cmd_stack)

    __call__ = cmd2
