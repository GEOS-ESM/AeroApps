import re

class Field(object):

    def __init__(self, name, info=None):

        self.name       = name
        self.dexpr      = ''
        self.fileID     = None
        self.expression = None

        if info is None:
            self.expression = name
        else:
            self.__dict__.update(info)
            self.set_dims()

    def update(self, field):

        self.merge(field)

        if self.expression is not None:
            return self

        if self.fileID:
            self.name += '.' + str(self.fileID)

        if self.dexpr:
            self.set_dims()
            self.name += self.dexpr

        return self

    def merge(self, field):

        dRE   = r'\(([A-Za-z0-9=\.,\-]+)\)'
        match = re.match(dRE, field.dexpr)

        if not match:
            return

        f_dexpr = match.group(1)

        if re.match(r'[\(,]' + f_dexpr + r'[\),]', self.dexpr):
            return

        if self.dexpr:
            match      = re.match(dRE, self.dexpr)
            s_dexpr    = match.group(1) 
            self.dexpr = '(' + s_dexpr + ',' + f_dexpr + ')'
        else:
            self.dexpr = field.dexpr

    def set_dims(self):

        if not self.fileID: return

        t_delta = 0
        dout    = []
        dims    = self.dexpr[1:-1].split(',') if self.dexpr else []
        dims_t  = [ d for d in dims if d[0:2] == 't=' ]

        for d in dims:

            if d[0:3] == 'td=':
                t_delta  = int(d[3:])
            else:
                dout.append(d)

        t = self.tm_index + t_delta
        if t != self.tm_index and not dims_t:
            dout.append('t=' + str(t))

        if dout: self.dexpr = '(' + ','.join(dout) + ')'
