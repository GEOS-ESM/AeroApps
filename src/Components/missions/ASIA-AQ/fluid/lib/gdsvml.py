import os
import re
import copy

udft = []

class GDSVML(object):

    def __init__(self):

        self.value      = None
        self.macro      = None
        self.macro_name = None

        self.default = { 'parea'   : '0 11 0 8.5',
                         'vpage'   : '0 11 0 8.5',
                         'clip'    : '0 11 0 8.5',
                         'frame'   : 'on',
                         'map'     : '15 1 1',
                         'mpt'     : '0 off',
                         'poli'    : 'on',
                         'mproj'   : 'latlon',
                         'mpdset'  : 'lowres',
                         'mpdraw'  : 'on',
                         'mpvals'  : '',
                         'line'    : '1 1 3',
                         'font'    : '0',
                         'grads'   : 'on',
                         'string'  : '1 bl 6 0',
                         'strsiz'  : '0.2 0.2',
                         'gxout'   : 'contour',
                         'ccolor'  : '0',
                         'cthick'  : '6',
                         'cstyle'  : '1',
                         'clopts'  : '-1 -1 0.09',
                         'ccols'   : '',
                         'clevs'   : '',
                         'csmooth' : 'off',
                         'xlopts'  : '1 5 0.09',
                         'ylopts'  : '1 5 0.09',
                         'rgb'     : {   '1' : '  0   0   0',
                                         '0' : '255 255 255',
                                         '2' : '250  60  60',
                                         '3' : '  0 220   0',
                                         '4' : ' 30  60 255',
                                         '5' : '  0 200 200',
                                         '6' : '240   0 130',
                                         '7' : '230 220  50',
                                         '8' : '240 130  40',
                                         '9' : '160   0 200',
                                        '10' : '160 230  50',
                                        '11' : '  0 160 255',
                                        '12' : '230 175  45',
                                        '13' : '  0 210 140',
                                        '14' : '130   0 220',
                                        '15' : '170 170 170'
                                     },

                         'mpt'    : { '0' : 'off 15 1 1',
                                      '1' : 'off 15 1 1',
                                      '2' : 'off 15 1 1'
                                    }
                       }


        self.fonts  = {}

        self.icolor = 31

        self.colors = {re.sub("\s+"," ",v.strip()):int(k)
                            for k,v in self.default['rgb'].iteritems()}

        self.colors = {}

        self.const  = {}

        self.state  = copy.deepcopy(self.default)

        self.event  = { k:0 for k in self.default }

        self.sticky = [ 'x', 'y', 'z', 't', 'lon', 'lat', 'lev', 'time',
                        'background', 'clip', 'display', 'frame', 'gxout',
                        'line', 'parea', 'vpage', 'xsize', 'font', 'grads',
                        'string', 'strsiz', 'mpdraw', 'mpdset', 'mproj', 'mpt',
                        'map', 'mpvals', 'poli', 'annot', 'black', 'clab',
                        'clopts', 'csmooth', 'cterp', 'cthick', 'rgb', 'rgba',
                        'xlopts', 'ylopts'
                      ]

        self.intrinsic = ['aave','ave','vint','abs','acos','asin','atan2',
                          'cos','exp','log','log10','pow','sin','sqrt','tan',
                          'const','maskout','skip','smth9','hcurl','hdivg',
                          'mag','gr2stn','oacres','re','stnave','stnmin',
                          'stnmax','cdiff','tloop','tvrh2q','tvrh2t',
                          'lat', 'lon', 'x', 'y', 'sum'
                         ]

        self.udf = self.readUDFT()

        self.latlon    = ['set x', 'set y', 'set lon', 'set lat']

        self.dimension = self.latlon + \
                         ['set z', 'set lev', 'set t', 'set time', 'set e']

        self.display = ['d','display']

        self.define = ['define','undefine']

        self.file = ['set dfile','open','sdfopen','xdfopen','close',
                     'reinit']

        self.annotate = ['draw string','draw title','draw xlab','draw ylab',
                         'draw line', 'draw mark', 'draw polyf', 'draw rec',
                         'draw recf', 'draw wxsym', 'draw map', 'draw shp',
                         'draw subtitle', 'draw tmstring', 'draw logo',
                         'draw label', 'draw grid', 'draw cbar', 'draw basemap',
                         'draw arrow', 'draw hilo']

        self.special = ['draw basemap', 'draw line', 'draw polyf',
                        'draw string', 'draw mark', 'draw hilo']

        self.addon   = ['set ARRFILL', 'set SLICE', 'set SKIP', 'set CLIP']

        self.reset = ['reinit', 'c', 'clear', 'reset', 'swap']

        self.rgb = ['set rgb', 'set rgba']

        self.font = ['set font']

        self.run  = ['run']

        self.set = self.dimension + self.addon + \
                   ['set background', 'set clip', 'set cmark', 'set display',
                    'set frame', 'set gridln', 'set gxout', 'set line',
                    'set missconn', 'set parea', 'set vpage', 'set xsize',
                    'set dignum', 'set digsiz', 'set font', 'set grads',
                    'set string', 'set strsiz', 'set timelab', 'set tlsupp',
                    'set mpdraw', 'set mpdset', 'set mproj', 'set mpt',
                    'set mpvals', 'set poli', 'set annot', 'set black',
                    'set ccolor', 'set ccols', 'set cint', 'set clab',
                    'set clevs', 'set clopts', 'set clskip', 'set cmax',
                    'set cmin', 'set csmooth', 'set cstyle', 'set cterp',
                    'set cthick', 'set map', 'set grid', 'set ylpos',
                    'set arrscl', 'set arrowhead']

        self.special_labels = ['label']
        self.header_labels  = ['main', 'subheader', 'header']
        self.trailer_labels = ['xlabel']
        self.label_types    = self.special_labels \
                            + self.header_labels \
                            + self.trailer_labels

    def get_state(self, state=None, default=None):

        if state is None:
            state = self.state

        if default is None:
            default = self.default

        if isinstance(state, dict):

            result = {}

            for key in state:
                d           = default.get(key,state[key])
                result[key] = self.get_state(state[key], d)

            return result
        else:

            list2 = re.sub("\s+"," ",default).split()
            list1 = re.sub("\s+"," ",state).split() + [None]*len(list2)

            if "auto" in list1:
                return ' '.join(list2)
      
            pick  = lambda pair: pair[0] if pair[0] is not None else pair[1]
            return ' '.join(map(pick, zip(list1, list2)))

    def get_color(self, cmd):

        cmd = re.sub("\s+"," ",cmd.strip()).split()
        rgb = ' '.join(cmd[3:])

        if len(cmd) == 4: return int(cmd[3])

        try:
            index = int(cmd[2])
        except:
            index = self.colors.get(rgb, None)
            if not index:
                index = self.icolor
                self.icolor += 1

        self.colors[rgb] = index
        return index

    def readUDFT(self):

        global udft
        if udft: return udft

        file = os.getenv('GA2UDXT', None)
        if not file: return []

        with open(file, 'r') as f:
            lines = [l.strip() for l in f]
            lines = [l.split() for l in lines if re.match(r'^udf', l)]
            udft  = [l[2] for l in lines if re.match(r'^\w+', l[2])]

        return udft

    def is_intrinsic(self, name):

        name.lower()

        if name in self.intrinsic:
            return True

        if name in self.udf:
            return True

        return False

    def is_auto(self, cmd):

        if '--auto' in cmd:
            return True

        return False

    def is_dimension(self, cmd):

        return self.search(cmd, self.dimension)

    def is_file(self, cmd):

        return self.search(cmd, self.file)

    def is_latlon(self, cmd):

        return self.search(cmd, self.latlon)

    def is_display(self, cmd):

        return self.search(cmd, self.display)

    def is_pendown(self, cmd):

        return self.search(cmd, self.display) \
            or self.search(cmd, self.annotate)

    def is_annotate(self, cmd):

        return self.search(cmd, self.annotate)

    def is_special(self, cmd):

        return self.search(cmd, self.special)

    def is_addon(self, cmd):

        return self.search(cmd, self.addon)

    def is_run(self, cmd):

        return self.search(cmd, self.run)

    def is_define(self, cmd):

        return self.search(cmd, self.define)

    def is_data_service(self, cmd):

        return self.is_dimension(cmd) or \
               self.is_file(cmd)      or \
               self.is_display(cmd)   or \
               self.is_define(cmd)

    def is_data(self, cmd):

        return self.is_display(cmd) or \
               self.is_define(cmd)

    def is_reset(self, cmd):

        return self.search(cmd, self.reset)

    def is_set(self, cmd):

        return self.search(cmd, self.set)

    def is_rgb(self, cmd):

        return self.search(cmd, self.rgb)

    def is_open_font(self, cmd):

        if self.search(cmd, self.font):

            cmd = re.sub("\s+"," ",cmd.strip()).split()
            if len(cmd) != 5: return False
            if cmd[3] != 'file': return False

            index = cmd[2]
            name  = cmd[4]

            if self.fonts.get(index, None) == name:
                return True

            self.fonts[index] = name

        return False
            
    def is_action(self, cmd):

        return self.macro is not None

    def is_macro(self, cmd):

        if cmd[0] == '&' and cmd != '&END':
            if self.const.get(cmd[1:], True): return True

        return False

    def is_end(self, cmd):

        if cmd == '&END':
            return True

        return False

    def is_attribute(self, cmd):

        if cmd[0] == '@':
            return True

        return False

    def is_quoted(self, name):

        qRE = r'([\"\'])+.*([\"\'])+$'

        match = re.match('^' + qRE, name)

        if not match:
            return False

        if match.group(1) != match.group(2):
            return False

        return True

    def in_macro(self):

        if self.macro_name:
            return True

        return False

    def register(self, constants):

        self.const.update(constants)

    def restore(self):

        state = {}

        for key in self.sticky:
            if key in self.state:
                state[key] = self.state[key]

        self.state = state

    def eval(self, cmd):

        command   = re.sub("\s+"," ",cmd.strip()).split()

        end_of_plot_event = (self.macro is not None)

        if end_of_plot_event:
            self.restore()

        self.macro = None
        self.value = None

        if self.is_display(cmd) and not self.in_macro():
           
            self.macro = '_'.join(['DISPLAY', self.state['gxout'].upper()])

        elif self.is_annotate(cmd) and not self.in_macro():

            self.value = ' '.join(command[2:])
            self.macro = '_'.join(['DRAW', command[1].upper()])

        elif self.is_run(cmd) and not self.in_macro():

            self.value = ' '.join(command[1:])
            self.macro = 'RUN'

        elif self.is_end(cmd):

            self.macro      = self.macro_name
            self.macro_name = None

        elif self.is_macro(cmd):

            self.macro_name = cmd[1:]

        elif self.is_attribute(cmd):

            self.event[command[0]] = 1
            self.state[command[0]] = ' '.join(command[1:])

        elif self.is_rgb(cmd):

            self.event['rgb'] = 1
            self.state['rgb'][command[2]] = ' '.join(command[3:])

        elif self.is_set(cmd):

            self.event[command[1]] = 1
            self.state[command[1]] = ' '.join(command[2:])

        elif self.is_reset(cmd):
            self.event = { k:0 for k in self.default }
            self.state = copy.deepcopy(self.default)

    def search(self, cmd, list):

        cmd = re.sub("\s+"," ",cmd.strip())

        if any(re.match(r'^'+s+' ', cmd+' ') for s in list):
            return True

        return False

    def __call__(self, key):

        return self.default.get(key,None)
