#--------------------------------------------------------------------------
#
#    Copyright (C) 2006-2008 by Arlindo da Silva <dasilva@opengrads.org>
#    All Rights Reserved.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation# using version 2 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY# without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program# if not, please consult  
#              
#              http://www.gnu.org/licenses/licenses.html
#
#    or write to the Free Software Foundation, Inc., 59 Temple Place,
#    Suite 330, Boston, MA 02111-1307 USA
#
#------------------------------------------------------------------------

"""
This module implements a Python interface to GrADS by means of
bi-directional pipes. It requires a small patch to GrADS to unbuffer
stdin/stdout and to allow capture of error codes. This modules
provides the basic GrADS class which implements the basic functionality to
start grads, send commands to it and to retrieve the text output
produced by *grads* in response to such command.  
"""

__version__ = '1.2.1'

import sys
import os

from time     import sleep
from datetime import datetime
from calendar import timegm
from string   import *
from types    import *
from math     import floor, ceil
from time     import time
from array    import array as array
from gahandle import *

# If possible, uses the subprocess module (Python 2.4 & new)
try:
    from subprocess import Popen, PIPE
    if os.name == 'nt':
        Kwds = dict( shell=True, bufsize=0,
                     stdin=PIPE, stdout=PIPE )
    elif os.name == 'java':
        Kwds = dict( shell=True, bufsize=4,
                     stdin=PIPE, stdout=PIPE )
    else:
        Kwds = dict( shell=True, bufsize=0,
                     stdin=PIPE, stdout=PIPE,
                     close_fds=True)
        
# If Python 2.3 or older, use the deprecated popen2()
except:
    from popen2     import popen2 as Popen
    if os.name == 'nt':
        Kwds = dict( mode='b' )
    else:
        Kwds = dict( bufsize=0 )

try:
    import jarray
    from java import util as jutil
    from org.opengrads.interfaces import GrADSObject
except:
    GrADSObject = object


class GrADSError(Exception):
    """
    Defines GrADS general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class GaCore(GrADSObject):
    """
    This is the core GrADS client class. It provides the basic functionality 
    to start the GrADS application, send commands to it and parse its
    output.
 
    _Methods provided:
        cmd    - sends a generic command to GrADS
        open   - opens a file, returning metadata in a dictionary
        query  - queries the GrADS state, returning properties in
                 a dictionary
        rline  - returns a given line from the GrADS stdout
        rword  - returns a given word from the GrADS stdout
        setdim - set dimension from dictionary returned by query()
        flush  - fluxes the communication pipes

    Methods for exchanging data arrays with GrADS as well as an interface
    to Python based graphics are provided in derived classes such as
    *GanNum* and *GaLab*.
    """
    
    def __init__ (self, 
                  Bin='grads', Echo=True, Opts='', Port=False, 
                  Strict=False, Verb=0, Window=None):
        """
        Starts the GrADS process using Popen function. Optional input
        parameters are:
        Bin     the name of the GrADS binary; you can enter a full path
                here, say, "/opt/mygrads/bin/gradshdf", or simply the
                name of other executables, e.g., "gradsdods" if this is
                on your PATH.
        Echo    if True output from GrADS will not be echoed on the screen
        Port    If set to True GrADS will start in portrait mode; default is
                to start GrADS in landscape mode.
        Opts    Use that to send additional options to GrADS, such as -c
        Strict  By default, it behaves a bit like perl: not raising
                exceptions for index errors, etc. If Strict=1, then
                exceptions will be raised.
        Verb    provides information about each GrADS command executed. If
                set to zero, it will be real quiet. When Verb=1 the rc code
                of failed GrADS commands will be printed. If verb=2 then
                the rc for all GrADS commands will be printed.
        Window  If True GrADS will start with a graphics window; the
                default is to start GrADS in batch mode.
        """

#       Default foe graphical window
#       ----------------------------
        if Window==None:
            try:
                if sys.stdin.isatty:  Window=True
                else:                 Window=False
            except:
                Window=True

#       Build GrADS command line
#       ------------------------
        cmdline = Bin + ' -u'
        if Window != True: cmdline = cmdline + ' -b'
        if Port   == True: cmdline = cmdline + ' -p'
        else:              cmdline = cmdline + ' -l'
        cmdline = cmdline + ' ' + Opts

#       Spawn GrADS process with bi-directional pipes
#       ---------------------------------------------
        try:
            # http://bit.ly/qvkUQK
            self.p = Popen(cmdline, **Kwds )                  # could be popen2() or Popen)

            if type(self.p) is TupleType:
                Reader, Writer = self.p                       # older popen2() 
            else:
                Reader, Writer = self.p.stdout, self.p.stdin  # newer Popen

        except Exception, e:
            raise e
            raise GrADSError, "Could not start the GrADS process with <"\
                           + cmdline + ">"

#       Record state
#       ------------
        self.Reader = Reader
        self.Writer = Writer
        self.Echo = Echo
        self.Verb = Verb
        self.Strict = Strict
        self.rc = 0

#       Parse out inital splash screen
#       -------------------------------
        rc = self._parseReader()

#       Record GrADS version
#       --------------------
        self.cmd('q config',Quiet=True)
        self.Version = self.rword(1,2)
        self.byteorder = self.rword(1,4)
        if ( self.byteorder!='big-endian' and self.byteorder!='little-endian'):
            self.byteorder = self.rword(1,3) # grads v2

#       Check availability of extensions
#       --------------------------------
        self.HAS_UDXT = False
        self.HAS_UDCT = False # old v1.9.0-rc1 style (deprecated)
        self.cmd('q udxt',Quiet=True)
        if self.rword(1,1) is not 'Invalid':
            self.HAS_UDXT = True
        else:
            self.cmd('q udct',Quiet=True) # as in v1.9.0-rc1
            if self.rword(1,1) is not 'Invalid':
                self.HAS_UDXT = True
                self.HAS_UDCT = True

#       Do we have the IPC extensions?
#       ------------------------------
        if self.HAS_UDXT is True:
            try:
                self.cmd('ipc_close',Quiet=True)
                self.HAS_IPC = True
            except GrADSError:
                self.HAS_IPC = False

#       OK, ready to go
#       ---------------
        if Verb: print "Started <"+cmdline+">, rc = ",rc


#........................................................................

    def __del__ ( self ):
        """ Sends GrADS the 'quit' command and close the pipes."""

        if self.Verb > 0: print "Stopping the GrADS process..."

        try:
#           self.cmd('quit')
            self.Reader.close()
            self.Writer.close()
        except:
            pass

#........................................................................

    def cmd ( self, gacmd, Quiet=False, Block=True, **kwopt ):
        """
        Sends a command to GrADS. When Block=True, the output is captured 
        and can be retrieved by methods rline() and rword(). On input,
        
        gacmd    can be a single GrADS command or several
                 commands separated by a new line ('\n').

        Quiet    if true does not echo anything to screen;
                 otherwise it honors the Echo preference
                 when the object was created.

        Block    whether the command blocks; if Block=False it will
                 return even if GrADS has not finished executing
                 the command; use it with caution and if you know
                 what you are doing.

        kwopt    used as variables for string interpolation. For example,

                      ga.cmd('set time $m$y', m='jan', y=1998)

                 is equivalent to:

                      ga.cmd('set time jan1998')

        IMPORTANT: Notice the names "$Quiet" and "$Block" cannot be used
                   for string interpolation.
                   
        """
        
        if len(kwopt)>0:
            Cmds = Template(gacmd).substitute(kwopt).split('\n')
        else:
            Cmds = gacmd.split('\n')
            
        Verb = self.Verb
        for cmd_ in Cmds:
            cmd = cmd_ + '\n'
            self.Writer.write(cmd)
            self.Writer.flush()
            if Block:
                rc = self._parseReader(Quiet)
                if rc != 0: 
                    if Verb==1:   print "rc = ", rc, ' for ' + cmd_
                    raise GrADSError, 'GrADS returned rc=%d for <%s>'%(rc,cmd_)
                else:
                    if Verb>1:    print "rc = ", rc, ' for ' + cmd_
        return
        
    __call__ = cmd

#........................................................................

    def flush ( self ):
        """ Flushes the GrADS pipes. This is primarily an internal 
        method that sometimes can be useful to the end-user. If for
        some reason the echo you see on the screen seems out of sync
        with the command you just issued, a call to *flush()* may be your
        ticket. (This situation is usually a symptom of an internal
        error.) 
        """

#       Issue a unique GrADS command
#       ----------------------------
        cmd = 'query config <TIME> %s </TIME>'%time()
        self.Writer.write(cmd+'\n')
        self.Writer.flush()

#       Position stream pointer after <IPC> marker
#       ------------------------------------------
        s = '<IPC> ' + cmd
        n = len(s)
        got = ''
        while got[:n] != s:
            got = self.Reader.readline()
            if got == '':
                raise GrADSError, "GrADS terminated while waiting for response"


#       OK, just issue another regular command and let the
#       parser do its job
#       --------------------------------------------------
        self.cmd('query config',Quiet=True)

#........................................................................

    def open ( self, fname, ftype='default', Quiet=False ):
        """
        Opens a GrADS file, returning the relevant metadata. On input,
        ftype can be used to specify which GrADS open command to open
        the file:
        ctl      the classic GrADS control file; uses the "open" command
        sdf      NetCDF, HDF or a DODS URL; uses the "sdfopen" command
        xdf      NetCDF, HDF or a DODS URL; uses the "xdfopen" command
        default  Open command is Guessed from the file name; the heuristic
                 algorithm works pretty well but is not perfect, therefore
                 the options above for when it fails.
        """

#       Determine the file opener
#       -------------------------
        FNAME = upper(fname) 
        FTYPE = upper(ftype) 
        if FTYPE == 'DEFAULT':
            if    FNAME[:7]  == 'HTTP://' or\
                  FNAME[-4:] == '.HDF'    or\
                  FNAME[-4:] == '.NC4'    or\
                  FNAME[-3:] == '.NC':          opener = 'sdfopen'
            elif  FNAME[-4:] == '.DDF':         opener = 'xdfopen'
            elif  FNAME[-4:] == '.XTL':         opener = 'xdfopen'
            else:                               opener = 'open'
        elif FTYPE == 'SDF':                    opener = 'sdfopen'
        elif FTYPE == 'XDF':                    opener = 'xdfopen'
        else:                                   opener = 'open'

#       Issue the GrADS command
#       -----------------------
        gacmd = opener + ' ' + fname
        try:
            self.cmd ( gacmd, Quiet=Quiet )
        except GrADSError:
            raise GrADSError, 'GrADS cannot open file <'+fname+'>'

#       Next, fill in the file handle
#       -----------------------------
        fid = -1
        for i in range(self.nLines, 0, -1):
            line = self.rline(i)
            if line[:9] == 'Data file' or line[:8] == 'SDF file':
                fid = int(self.rword(i,8))
                break
            
        fh = self.query('file %d'%fid, Quiet=True)
        qh = self.query('ctlinfo %d'%fid, Quiet=True)
        fh.undef = qh.undef # will need that for fwrite later

        return fh
                      
#........................................................................

    def jopen ( self, fname, ftype='default', Quiet=False ):
        """
        Similar to open() but returns a Java HashMap instead.
        This method only makes sense under Jython.
        """
        qh = self.open(fname,ftype,Quiet)
        return _toHashMap(qh)

#........................................................................

    def query ( self, what, Quiet=False ):
        """
        Queries GrADS internal state and returns a GaHandle object
        with the results of the query:

           qh = self.query(what)

        where *what* is the name of the properpty being queried.
        When Quiet=True echo of GrADS commands is temporarily 
        suppressed. The following properties are implemented:

        qh = self.query ( "file #" )
 
           Returns information about the file number "#". If the file number
           is ommitted the default file is used. The following *qh* attributes
           are defined:
 
              fid        - the file Id number
              title      - the title of the dataset
              desc       - the description of the dataset
              bin        - the binary file name
              type       - file type
              nx         - number of longitudinal points
              ny         - number of latitudinal points
              nz         - number of vertical levels
              nt         - number of times on file
              ne         - number of ensember members in file
              nvars      - number of variables on file
              vars       - list of variable names in the file
              var_levs   - list of number of levels for each variable
              var_titles - list of long names for each variable
              var_info   - list of tuplet triplets (var,levs,title) for each 
                           variable to allow easy iterating:
                             qh = ga.query('file')
                             for varname, nlevels, vartitle in qh.var_info:
                                 # ...

        qh = self.query ( "ctlinfo" )
                                                                                
           Returns information about the file coordinate variables. The following
           *qh* attributes are defined:
             
           dset   - binary dataset name
           title  - dataset title
           dtype  - dataset type
           undef  - undefined value

           nx     - number of longitudes
           xtype  - "linear" or "levels"
           x0     - longitude origin (only if xtype = "linear")
           dx     - longitude increment (only if xtype = "linear")
           xlevs  - longitude levels (only if xtype = "levels")

           ny     - number of latitudes
           ytype  - "linear" or "levels"
           y0     - latitude origin (only if ytype = "linear")
           dy     - latitude increment (only if ytype = "linear")
           ylevs  - latitude levels (only if ytype = "levels")

           nz     - number of vertical levels
           ztype  - "linear" or "levels"
           z0     - vertical level origin (only if ztype = "linear")
           dz     - vertical level increment (only if ztype = "linear")
           zlevs  - vertical level levels (only if ztype = "levels")

           nt     - number of times
           ttype  - "linear" or "levels"
           t0     - time origin 
           dt     - time increment (in minutes)

           ne     - number of ensembles
           etype  - "linear" or "levels"
           e0     - ensemble origin (only if etype = "linear")
           de     - ensemble increment (only if etype = "linear")
           elevs  - ensemble levels (only if etype = "levels")
              
        qh = self.query ( "dims" )
                                                                                
           Returns information about the dimension environment. The following
           *qh* attributes are defined:
                                                                                
              dfile   - default file number

              x_state - x-coordinate state: "fixed" or "varying"
              lon     - longitudinal range
              x       - x-index range (may not be an integer)
              xi      - bracketing x-index range (always an integer)
              nx      - number of longitudes (xi[1]-xi[0]+1)

              y_state - y-coordinate state: "fixed" or "varying"
              lat     - latitudinal range
              y       - y-index range (may not be an integer)
              yi      - bracketing y-index range (always an integer)
              ny      - number of latitudes (yi[1]-yi[0]+1)

              z_state - z-coordinate state: "fixed" or "varying"
              lev     - level range
              z       - z-index range (may not be an integer)
              zi      - bracketing z-index range (always an integer)
              nz      - number of levels (zi[1]-zi[0]+1)

              t_state - time-coordinate state: "fixed" or "varying"
              time    - time range
              t       - t-index range
              ti      - bracketing t-index range (always an integer)
              nt      - number of times (ti[1]-ti[0]+1)

              e_state - ens-coordinate state: "fixed" or "varying"
              ens     - ens range
              e       - e-index range
              ei      - bracketing e-index range (always an integer)
              ne      - the number of ensemble members (ei[1]-ei[0]+1)

            The ensemble dimension information will be reported for both 
            GrADS v1.x and GrADS v2.0. Since GrADS v1.x does not support
            ensemble dimensions, the number of ensemble members *ne* is
            set to 1 in this case.

        qh = self.query ( "udc" )
                                                                                
           Returns information about User Defined Commands. The following
           *qh* attributes are defined:

           udcs - list of UDC names

        qh = self.query ( "udf" )
                                                                                
           Returns information about User Defined Functions. The following
           *qh* attributes are defined:

           udfs - list of UDF names

        qh = self.query ( "time" )
                                                                                
           Returns information about current range of time dimension, including
		   unix timestamp values. 

        """

        qh = GaHandle(what) # handle to hold results

#       Query GrADS
#       -----------
        try:
            self.cmd('query '+what,Quiet)
            qh.rc = self.rc
        except GrADSError: 
            raise GrADSError, 'Cannot query GrADS about <'+what+'>'

#       Parse output
#       ------------
        tokens = split(what)
        what = tokens[0]

#       query dims
#       ----------
        if what == 'dims':
            
            qh.dfile = self.rword(1,5)
            qh.rank = 0
            qh.x_state = self.rword(2,3)
            if qh.x_state == 'fixed':
                qh.lon = (float(self.rword(2,6)), float(self.rword(2,6)))
                qh.x   = (float(self.rword(2,9)), float(self.rword(2,9)))
            else:
                qh.rank = qh.rank+1
                qh.lon = (float(self.rword(2,6)), float(self.rword(2,8)))
                qh.x   = (float(self.rword(2,11)), float(self.rword(2,13)))
            qh.y_state = self.rword(3,3)
            if qh.y_state == 'fixed':
                qh.lat = (float(self.rword(3,6)), float(self.rword(3,6)))
                qh.y   = (float(self.rword(3,9)), float(self.rword(3,9)))
            else:
                qh.rank = qh.rank+1
                qh.lat = (float(self.rword(3,6)), float(self.rword(3,8)))
                qh.y   = (float(self.rword(3,11)), float(self.rword(3,13)))
            qh.z_state = self.rword(4,3)
            if qh.z_state == 'fixed':
                qh.lev = (float(self.rword(4,6)), float(self.rword(4,6)))
                qh.z   = (int(float(self.rword(4,9))), int(float(self.rword(4,9))))
            else:
                qh.rank = qh.rank+1
                qh.lev = (float(self.rword(4,6)), float(self.rword(4,8)))
                qh.z   = (int(float(self.rword(4,11))), int(float(self.rword(4,13))))

            qh.t_state = self.rword(5,3)
            if qh.t_state == 'fixed':
                qh.time = (self.rword(5,6), self.rword(5,6))
                qh.t    = (int(float(self.rword(5,9))), int(float(self.rword(5,9))))
            else:
                qh.rank = qh.rank+1
                qh.time = (self.rword(5,6),  self.rword(5,8))
                qh.t    = (int(float(self.rword(5,11))), int(float(self.rword(5,13))))

#           Ensemble dimension only supported in GrADS v2.0
#           -----------------------------------------------
            if self.Version[1] is '2':
                qh.e_state = self.rword(6,3)
                if qh.e_state == 'fixed':
                    qh.ens = (self.rword(6,6), self.rword(6,6))
                    qh.e    = (int(self.rword(6,9)), int(self.rword(6,9)))
                else:
                    qh.rank = qh.rank+1
                    qh.ens  = (self.rword(6,6),  self.rword(6,8))
                    qh.e    = (int(self.rword(6,11)), int(self.rword(6,13)))

#           Fake a dimension environment in GrADS v1.x
#           ------------------------------------------
            else:
                qh.e_state = 'fixed'
                qh.ens = (1,1)
                qh.e   = (1,1)

#           Index space: 
#             x  ... indices as reported by GrADS (not always a whole number)
#             xi ... bracketing indices (always a whole number)
#           -----------------------------------------------------------------
            qh.xi = (int(floor(qh.x[0])),int(ceil(qh.x[1])))
            qh.yi = (int(floor(qh.y[0])),int(ceil(qh.y[1])))
            qh.zi = (int(floor(qh.z[0])),int(ceil(qh.z[1])))
            qh.ti = (int(floor(qh.t[0])),int(ceil(qh.t[1])))
            qh.ei = (int(floor(qh.e[0])),int(ceil(qh.e[1])))

#           Number is based on bracketing integer indices
#           ---------------------------------------------
            qh.nx = qh.xi[1] - qh.xi[0] + 1
            qh.ny = qh.yi[1] - qh.yi[0] + 1
            qh.nz = qh.zi[1] - qh.zi[0] + 1
            qh.nt = qh.ti[1] - qh.ti[0] + 1
            qh.ne = qh.ei[1] - qh.ei[0] + 1

            qh.tyme = [gat2dt(t) for t in qh.time]

#       query file [fid]
#       ----------------
        elif what == 'file':

            tokens = split(self.rline(1))
            qh.title = join(tokens[3:])

            qh.fid   = int(self.rword(1,2))
            qh.desc  = self.rword(2,2)
            qh.bin   = self.rword(3,2)
            qh.type  = self.rword(4,3)
            qh.nx    = int(self.rword(5,3))
            qh.ny    = int(self.rword(5,6))
            qh.nz    = int(self.rword(5,9))
            qh.nt    = int(self.rword(5,12))

            if self.Version[1] is '2':
                qh.ne = int(self.rword(5,15))
            else:
                qh.ne = 1

            qh.nvars = int(self.rword(6,5))

            vars       = [ ] 
            var_levs   = [ ] 
            var_titles = [ ]
            var_dims = []
            
            for i in range(7, 7+qh.nvars):
                line  = self.rline(i)
                words = line.strip().split(' ')
                nw=len(words)
                currvar   = str(words[0])
                nc=1
                if(words[nc]==''): nc=nc+1
                currlev   = int(words[nc])
                nc=nc+1
                if(words[nc]==''): nc=nc+1
                currdims = words[nc] 
                nc=nc+1
                currtitle = ' '.join( words[nc:] ).lstrip()
                
                vars.append(currvar)
                var_levs.append(currlev)
                var_dims.append(currdims)
                var_titles.append(currtitle)

            qh.vars       = list(vars)        # deprecated
            qh.var_levs   = list(var_levs)    # deprecated
            qh.var_titles = list(var_titles)  # deprecated

            qh.var_info   = zip(vars, var_levs, var_dims, var_titles)

            qh.Vars = {}
            for v,levs,dims,title in qh.var_info:
                qh.Vars[v] = GaHandle('v')
                qh.Vars[v].name = v
                qh.Vars[v].levs = levs
                qh.Vars[v].dims = dims
                qh.Vars[v].long_name = title
            
            
#       query UDX table
#       ---------------
        elif what == 'udc':
            qh.udcs = [ ]
            if self.nLines>6:
                for i in range(self.nLines-7):
                    qh.udcs.append(self.rword(i+6,1))

        elif what == 'udf':
            qh.udfs = [ ]
            if self.nLines>6:
                for i in range(self.nLines-7):
                    qh.udfs.append(self.rword(i+6,1))

        elif what == 'time':

            qh.t1 = self.rword(1,3)
            qh.t2 = self.rword(1,5)

            qh.tyme1 = gat2dt(qh.t1)
            qh.t1_unix=timegm(qh.tyme1.timetuple())

            qh.tyme2 = gat2dt(qh.t2)
            qh.t2_unix=timegm(qh.tyme2.timetuple())
			

#       query dims
#       ----------
        elif what == 'ctlinfo':

            for i in range(self.nLines):
                if self.rword(i,1) == "dset":  i_dset = i
                if self.rword(i,1) == "title": i_title = i
                if self.rword(i,1) == "undef": i_undef = i
                if self.rword(i,1) == "dtype": i_dtype = i
                if self.rword(i,1) == "xdef":  i_xdef = i
                if self.rword(i,1) == "ydef":  i_ydef = i
                if self.rword(i,1) == "zdef":  i_zdef = i
                if self.rword(i,1) == "tdef":  i_tdef = i
                if self.rword(i,1) == "edef":  i_edef = i

            qh.dset = self.rword(i_dset,2)
            qh.title = join(self.Words[i_title][1:])
            qh.undef = float(self.rword(i_undef,2))
            try:
                qh.dtype = self.rword(i_dtype,2)
            except:
                qh.dtype = "ctl"
            qh.nx = int(self.rword(i_xdef,2))
            qh.xtype = self.rword(i_xdef,3)
            if qh.xtype == 'linear':
                qh.x0 = float(self.rword(i_xdef,4))
                qh.dx = float(self.rword(i_xdef,5))
            elif qh.xtype == 'levels':
                qh.xlevs = [float(x_) for x_ in self.Words[i_xdef][3:]]
                cnt = 0
                while ( len(qh.xlevs) != qh.nx ):
                    cnt+=1
                    qh.xlevs += [float(x_) for x_ in self.Words[i_xdef+cnt][:]]
                 
            qh.ny = int(self.rword(i_ydef,2))
            qh.ytype = self.rword(i_ydef,3)
            if qh.ytype == 'linear':
                qh.y0 = float(self.rword(i_ydef,4))
                qh.dy = float(self.rword(i_ydef,5))
            elif qh.ytype == 'levels':
                qh.ylevs = [float(y_) for y_ in self.Words[i_ydef][3:]]
                cnt = 0
                while ( len(qh.ylevs) != qh.ny ):
                    cnt+=1
                    qh.ylevs += [float(y_) for y_ in self.Words[i_ydef+cnt][:]]
                 
            qh.nz = int(self.rword(i_zdef,2))
            qh.ztype = self.rword(i_zdef,3)
            if qh.ztype == 'linear':
                qh.z0 = float(self.rword(i_zdef,4))
                qh.dz = float(self.rword(i_zdef,5))
            elif qh.ztype == 'levels':
                qh.zlevs = [float(z_) for z_ in self.Words[i_zdef][3:]]
                cnt = 0
                while ( len(qh.zlevs) != qh.nz ):
                    cnt+=1
                    qh.zlevs += [float(z_) for z_ in self.Words[i_zdef+cnt][:]]

            qh.nt = int(self.rword(i_tdef,2))
            qh.ttype = self.rword(i_tdef,3)
            qh.t0 = self.rword(i_tdef,4)
            qh.dt = float(self.rword(i_tdef,5)[:-2])
                 
            try:
                qh.ne = int(self.rword(i_edef,2))
                qh.etype = self.rword(i_edef,3)
                if qh.etype == 'linear':
                    qh.e0 = float(self.rword(i_edef,4))
                    qh.de = float(self.rword(i_edef,5))
                elif qh.etype == 'levels':
                    qh.elevs = [float(e_) for e_ in self.Words[i_edef][3:]]
                    cnt = 0
                    while ( len(qh.elevs) != qh.ne ):
                        cnt+=1
                        qh.elevs += [float(e_) for e_ in self.Words[i_edef+cnt][:]]
            except:
                qh.ne = 1
                qh.etype = "linear"
                qh.e0 = 1
                qh.de = 1

#       Not implemented yet
#       -------------------
        else:
            qh.rc = -1
            print "query('%s') method not fully implemented yet"%what

        return qh

#........................................................................

    def jquery ( self, what, Quiet=False ):
        """
        Queries GrADS internal state and returns a Java Hashtable with 
        the results. This method only makes sense under Jython.
        """
        qh = self.query(what,Quiet)
        return _toHashMap(qh)

#........................................................................

    def coords ( self ):
        """
        Returns a GaHandle object with attributes defining the
        following coordinate variables associated with the current
        dimension environmnet:

        name  ---  set to "coords"
        undef ---  based on output undef value (v2.0.a2 or newer)
                   based on default file, otherwise
        dims  ---  list with names of varying dimensions
        denv  ---  dimension environment: result of query("dims")
        ens   ---  list with ensemble, e.g., ['e1', 'e2' ]
        time  ---  list with time, e.g., ['00Z01JAN1987', '12Z01JAN1987' ]
        lev   ---  float array with vertical levels
        lat   ---  float array with latitudes
        lon   ---  float array with longitudes

        Notice that this method defines many of the attributes of
        a GaGrid object.

        """

        dh = self.query("dims",Quiet=True)

#       Construct handle
#       ----------------
        ch = GaHandle("coords")
        ch.denv = dh

#       Varying dimensions
#       ------------------
        ch.dims = [ 'ens','time', 'lev', 'lat', 'lon' ]
        if ch.denv.ne==1: ch.dims.remove('ens')
        if ch.denv.nt==1: ch.dims.remove('time')
        if ch.denv.nz==1: ch.dims.remove('lev')
        if ch.denv.ny==1: ch.dims.remove('lat')
        if ch.denv.nx==1: ch.dims.remove('lon')
        ch.shape = [ ch.denv.ne, ch.denv.nt, ch.denv.nz, ch.denv.ny, ch.denv.nx]
        if ch.denv.ne==1: ch.shape.remove(1)
        if ch.denv.nt==1: ch.shape.remove(1)
        if ch.denv.nz==1: ch.shape.remove(1)
        if ch.denv.ny==1: ch.shape.remove(1)
        if ch.denv.nx==1: ch.shape.remove(1)

        self.cmd("set x 1",Quiet=True)
        self.cmd("set y 1",Quiet=True)
        self.cmd("set z 1",Quiet=True)
        if self.Version[1] is '2':
            self.cmd("set e 1",Quiet=True)

#       ensemble coordinates
#       --------------------
        if self.Version[1] is '2':
            ch.ens = []
            for n in range(dh.ne):
                e = dh.ei[0] + n
                self.cmd("set e %d"%e,Quiet=True)
                self.cmd("q ens",Quiet=True)
                ch.ens.append(self.rword(1,3))
                self.cmd("set e 1",Quiet=True)
                
#       Time coordinates
#       ----------------
        ch.time = []
        for n in range(dh.nt):
            t = dh.ti[0] + n
            self.cmd("set t %d"%t,Quiet=True)
            self.cmd("q time",Quiet=True)
            ch.time.append(self.rword(1,3))
        self.cmd("set t 1",Quiet=True)

#       Level coordinates
#       -----------------
        self.cmd("set z %d %d"%dh.zi,Quiet=True)
        ch.lev  = self.eval('lev')
        self.cmd("set z 1",Quiet=True)

#       Latitude coordinates
#       --------------------
        self.cmd("set y %d %d"%dh.yi,Quiet=True)
        ch.lat  = self.eval('lat')
        self.cmd("set y 1",Quiet=True)

#       Longitude coordinates
#       ---------------------
        self.cmd("set x %d %d"%dh.xi,Quiet=True)
        ch.lon = self.eval('lon')
        self.cmd("set x 1",Quiet=True)

#       Retore dimension environment
#       ----------------------------
        if self.Version[1] is '2':
            self.cmd("set e %d %d"%dh.ei,Quiet=True)
        self.cmd("set t %d %d"%dh.ti,Quiet=True)
        self.cmd("set z %d %d"%dh.zi,Quiet=True)
        self.cmd("set y %d %d"%dh.yi,Quiet=True)
        self.cmd("set x %d %d"%dh.xi,Quiet=True)

#       Undef
#       -----
        try:
            self.cmd('q undef',Quiet=True)
            ch.undef = float(self.rword(1,7))
        except: # for legacy reasons, before q undef was available
            self.cmd("q ctlinfo",Quiet=True)
            for i in range(self.nLines):
                if self.rword(i+1,1) == 'undef':
                    ch.undef = float(self.rword(i+1,2))
                    break

#       All done
#       --------
        return ch

#........................................................................

    def jcoords ( self ):
        """
        Returns a GaHandle object with attributes defining the
        coordinate variables associated with the current dimension
        environment. This is a Java wrapper; see coords() for 
        additional information.
        """
        qh = self.coords()
        return _toHashMap(qh)

#........................................................................

    def rline ( self, i=None ):
        """Returns the ith line of the most recent GrADS output.
           If *i* is not specified, returns the number of lines available"""
        if i==None:
            return self.nLines
        try:
            return self.Lines[i]
        except IndexError:
            if self.Strict: raise IndexError, "Invalid line number%d"%i
            else:           return '' # no fuss, no muss

    def rword ( self, i, j ):
        """Returns the jth word of the ith line of the most recent
        GrADS output"""
        try:
            return self.Words[i][j-1]
        except IndexError:
            if self.Strict: raise IndexError, \
                                  "Invalid word indices (%d,%d)"%(i,j)
            else:           return '' # no fuss, no muss

#.....................................................................

    def eval ( self, expr ):
        """
        Exports GrADS expression *expr*, returning a flat Python array.

            a = self.eval(expr)

        where

            a  ---  flat Python array (or flat jarray under Jython)
            
        Generalized Expressions
        =======================

        For convenience, *expr* can also be an array.  In such
        cases, the input Field is just returned back. This
        *overloading* is useful for writing high level functions that
        work equaly on a GrADS expression to be exported or on arrays.

        """

#       For convenience, allows calls where expr is not a string, in which
#        case it returns back the input field or raise an exception
#       -------------------------------------------------------------------
        if type(expr) in StringTypes:
            pass # OK, will proceed to export it from GrADS
        elif isinstance(expr,array):
            return expr # just return input
        else:
            raise GrADSError, "input <expr> has invalid type"

#       Retrieve dimension environment
#       ------------------------------
        dh = self.query("dims", Quiet=True) 
        nx, ny, nz, nt, ne = (dh.nx, dh.ny, dh.nz, dh.nt, dh.ne)

#       Tell GrADS to write expression to pipe
#       --------------------------------------
        self.cmd('query gxout', Quiet=True) 
        gxout = self.rword(4,6) # save gxout state
        self.cmd('set gxout fwrite', Quiet=True) 
        self.cmd('set fwrite -', Quiet=True)
 
#       For now, can only handle up to 3 varying dimensions
#       ---------------------------------------------------
        if dh.rank<=2: 
            self.cmd('display %s'%expr, Block=False) # non-blocking
        elif dh.rank==3: # xyz, xyt, xzt, yzt
            if   ne>1: self.cmd('set loopdim e', Quiet=True)
            elif nt>1: self.cmd('set loopdim t', Quiet=True)
            elif nz>1: self.cmd('set loopdim z', Quiet=True)
            self.cmd('display %s'%expr, Block=False) # non-blocking
        else:
            raise GrADSError, 'can only handle 3 varying dimensions'

#       Position stream pointer after <FWRITE> marker
#       ---------------------------------------------
        got = ''
        while got[:8] != '<FWRITE>' :
            got = self.Reader.readline()
            if got[:13] == 'Syntax Error:':
                self.flush()
                raise GrADSError, "Syntax Error - cannot evaluate <%s>"%expr
            elif got == '':
                raise GrADSError, "GrADS terminated while waiting for response"

#       Attempt to read from pipe
#       -------------------------
        a = array('f')
        try:
            n = ne*nt*nz*ny*nx
            a.fromfile(self.Reader,n)
            rc = 0
        except:
            rc = 1

#       Restore gxout settings
#       ----------------------
        self.cmd('disable fwrite')
        self.cmd('set gxout %s'%gxout, Quiet=True) 

#       Something went wrong
#       --------------------
        if rc or len(a) < n:
            raise GrADSError, 'problems evaluating <'+expr+'>'

#       All is well
#       -----------
        else:
            if os.name == 'java':
                if self.byteorder == 'little-endian':
                    a.byteswap() # JVM is always big-endian
                try:
                    return jarray.array(a,'f') # jython 2.3 & older needed this
                except:
                    return a # same as jarray in jython 2.5
            else:
                return a

#.....................................................................

    def setdim (self, dh):
        """
        Sets the dimension environment. On input, *dh* is usually obtained
        as the output of:

            dh = ga.query('dims')
        """
        try:
            self.cmd("set x %d %d"%dh.x,Quiet=True)
            self.cmd("set y %d %d"%dh.y,Quiet=True)
            self.cmd("set z %d %d"%dh.z,Quiet=True)
            self.cmd("set t %d %d"%dh.t,Quiet=True)
            if self.Version[1] is '2':
                self.cmd("set e %d %d"%dh.e,Quiet=True)
        except GrADSError:
            raise GrADSError, 'Cannot restore dimension environment'

#........................................................................

#   This should be private
#   ----------------------
    def _parseReader ( self, Quiet=False ):
        """
        Internal method to parse the GrADS output. Not user callable.
        """

        if Quiet:  Echo = False
        else:      Echo = self.Echo
        Lines = []
        Words = []
        
#       Discard debris
#       --------------
        got = ''
        while got[:5] != '<IPC>' :
            got = self.Reader.readline()
            if got == '':
                raise GrADSError, "GrADS terminated while waiting for response"
            tokens = split(got)

#       Record GrADS command
#       --------------------
        Lines.append(got)
        Words.append(tokens)
###     if Echo: print got

#       Next, record the GrADS output
#       -----------------------------
        rc = -99
        got = self.Reader.readline()
        while got[:6] != '</IPC>':
            tokens = split(got)
            if got[:4] == '<RC>':
                rc = tokens[1]
            else:
                Lines.append(got[:-1])
                Words.append(tokens)
                if Echo: print got[:-1]
            got = self.Reader.readline()
            if got == '':
                raise GrADSError, "GrADS terminated while waiting for response"

#       Save it for later
#       -----------------
        self.Lines = Lines
        self.nLines = len(Lines)-2
        self.Words = Words
        self.rc    = int(rc)

        return self.rc

#........................................................................

# For backward compatibility

GrADS = GaCore

#.....................................................................

def _nint(str):
    """
    Nearest integer, internal use.
    """
#   return int(float(str)+0.5)
    x = float(str)
    if x >= 0: return int(x+0.5) 
    else:      return int(x-0.5)

#.....................................................................

def _toHashMap(qh):
    """
    Convert attributes to a Java HashMap.
    """
    if os.name != 'java':
        raise GrADSError, 'conversion to HashMap requires Java'
    jh = jutil.HashMap();
    for key in qh.__dict__.keys():
        jkey = key
        value = qh.__dict__[key];
        if len(key) == 1:
            jkey = key+key  # Matlab cannot handle single char
        if ( type(value) == TupleType or type(value) == ListType ):
            vtype = type(value[0])
            if vtype == FloatType:
                value = jarray.array(value,'d')
            elif vtype == IntType:
                value = jarray.array(value,'i')
            elif vtype == LongType:
                value = jarray.array(value,'l')
        elif isinstance(value,GaHandle):
            value = _toHashMap(value)
        jh[jkey] = value
    return jh


# --------------------------------------------------------------------------------

__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

def gat2dt(gat):
    """
    Convert grads time to datetime.
    """
    time, date = gat.upper().split('Z')
    if time.count(':') > 0:
        h, m = time.split(":")
    else:
        h = time
        m = '0'
    mmm = date[-7:-4]
    dd, yy = date.split(mmm)
    mm = __Months__.index(mmm) + 1
    dt = datetime(int(yy),int(mm),int(dd),int(h),int(m))
    return dt

def dt2gat(t):
    """
    Convert datetime to grads time.
    """
    gat = "%d:%dZ%d%s%d"%(t.hour,t.minute,t.day,__Months__[t.month-1],t.year)
    return gat

