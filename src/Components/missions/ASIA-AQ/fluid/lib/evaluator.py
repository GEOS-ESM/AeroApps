import re
import sys
import field
import gdsvml

class Evaluator(object):

    def __init__(self, request, config, dataservice=None):

        self.ds      = dataservice
        self.config  = config
        self.request = request
        self.lang    = gdsvml.GDSVML()
        self.defined = []

        self.quote   = False
        self.qchar   = "'"
        self.qstring = None

        self.file    = None

    def evaluate(self, f):

        if f.expression is None:
            return f.name

        string = ''

        for obj in self.parse(f.expression):

            if type(obj) is field.Field:
                string += self.evaluate(obj.update(f))
            else:
                string += obj

        return string

    def parse(self, expr):

        """Parses Augmented GrADS Expressions

        Python generator to deconstruct augmented GrADS expressions
        and resolve special tokens into runtime values.

        **Args:**
            expr : string : input
                GrADS expression (may have augmented syntax - see notes)

        **Yields:**
            obj : string or Field object.

                string: string fragment from expression.

                Field: Field object representing a new data field name
                with possible dimension expression.

        **Notes:**
            The augmented syntax may include the following constructs:

            1. Fully qualified variable names with stream and collection
               components. For example:

                   slp.inst3_3d_asm_Np
                   slp.GEOSFC.inst3_3d_asm_Np

            2. Time delta dimension expression:

                   slp(td=-1)

            3. Time delta dimension argument(s):

                   sum(slp, td=-3, td=0)
            

        **See Also:**
            field.py

        **Raises:**
            none
        """

        # Declare the regular expression patterns

        sep   = r'\.'
        dARG  = r'([A-Za-z]+=)' # Dimension argument
        dRE   = r'(\([A-Za-z0-9=\.,\-]+\))' # Dimension expression
        fRE   = r'([\"\']*[A-Za-z_][A-Za-z0-9_]*[\"\']*)' # Field name

        while expr:

            match1 = re.match('^' + fRE, expr)
            match2 = re.match('^' + fRE + sep + fRE, expr)
            match3 = re.match('^' + fRE + sep + fRE + sep + fRE, expr)

            match1a = re.match('^' + fRE + dRE, expr)
            match2a = re.match('^' + fRE + sep + fRE + dRE, expr)
            match3a = re.match('^' + fRE + sep + fRE + sep + fRE + dRE, expr)

            match4  = re.match('^' + fRE + sep + r'[1-9]+', expr)

            number = self.number_match(expr)
            td_arg = re.match(r'^(td=\-*[0-9]+)', expr) # Time delta expression
            dimension_arg = re.match('^' + dARG, expr)

            in_quote = self.in_quoted_string(expr[0])

            if in_quote:
                pass
            elif self.qstring:
                obj    = self.qstring
                length = 1
            elif number:
                name   = number.group(1)
                obj    = name
                length = len(name)
            elif td_arg:
                if self.file:
                    dt     = int(td_arg.group(1).split('=')[-1])
                    name   = 't=' + str(dt + self.file['tm_index'])
                else:
                    name   = td_arg.group(1)
                obj    = name
                length = len(td_arg.group(1))
            elif dimension_arg:
                name   = dimension_arg.group(1)
                obj    = name
                length = len(name)
            elif match1 and self.lang.is_intrinsic(match1.group(1)):
                name   = match1.group(1)
                obj    = name
                length = len(name)
            elif match4:
                name   = match4.group(1)
                obj    = name
                length = len(name)
            elif match3a:
                name   = match3a.group(1)
                s      = match3a.group(2)
                c      = match3a.group(3)
                d      = match3a.group(4)
                obj    = self.get_field(name, stream=s, collection=c, dexpr=d)
                length = len('.'.join([name , s, c])+d)
            elif match3:
                name   = match3.group(1)
                s      = match3.group(2)
                c      = match3.group(3)
                obj    = self.get_field(name, stream=s, collection=c)
                length = len('.'.join([name, s, c]))
            elif match2a:
                name   = match2a.group(1)
                c      = match2a.group(2)
                d      = match2a.group(3)
                obj    = self.get_field(name, collection=c, dexpr=d)
                length = len('.'.join([name, c])+d)
            elif match2:
                name   = match2.group(1)
                c      = match2.group(2)
                obj    = self.get_field(name, collection=c)
                length = len('.'.join([name, c]))
            elif match1a:
                name   = match1a.group(1)
                d      = match1a.group(2)
                obj    = self.get_field(name, dexpr=d)
                length = len(name+d)
            elif match1:
                name   = match1.group(1)
                obj    = self.get_field(name)
                length = len(name)
            else:
                obj    = expr[0]
                length = 1

            if not in_quote: yield obj  
    
            if length == len(expr):
                expr = None
            else:
                expr = expr[length:]

    def number_match(self, expr):

        eRE1   = r'([0-9]+[\.]{1}[0-9]*[Ee]{1}[+-]{1}[0-9]+)'
        eRE2   = r'([\.]{1}[0-9]+[Ee]{1}[+-]{1}[0-9]+)'
        eRE3   = r'([0-9]+[Ee]{1}[+-]{1}[0-9]+)'
        rRE1   = r'([0-9]+[\.]{1}[0-9]*)'
        rRE2   = r'([\.]{1}[0-9]+)'
        iRE    = r'([0-9]+)'

        ematch1 = re.match('^' + eRE1, expr)
        ematch2 = re.match('^' + eRE2, expr)
        ematch3 = re.match('^' + eRE3, expr)
        rmatch1 = re.match('^' + rRE1, expr)
        rmatch2 = re.match('^' + rRE2, expr)
        imatch  = re.match('^' + iRE , expr)

        if   ematch1:
            return ematch1
        elif ematch2:
            return ematch2
        elif ematch3:
            return ematch3
        elif rmatch1:
            return rmatch1
        elif rmatch2:
            return rmatch2
        elif imatch:
            return imatch
        else:
            return None

    def in_quoted_string(self, c):

        if not self.quote and (c == "'" or c == '"'):
            self.quote = True
            self.qchar = c
            self.qstring = ''
        elif self.quote and c == self.qchar:
            self.quote = False
        elif self.quote:
            self.qstring += c
        else:
            self.qstring = None

        return self.quote

    def get_field(self, name, stream=None, collection=None, dexpr=None):

        if self.lang.is_quoted(name) or self.is_defined(name):
            return name

        in_stream     = stream
        in_collection = collection

        if stream is None:
            stream = self.request.get('stream','')

        if collection is None:
            collection = self.request.get('collection','')

        if dexpr is None:
            dexpr = ''

        lname = name.lower()

#       Check to see if the referenced field
#       is a derived quantity listed in the 
#       configuration.
#       ===================================

        derived = self.config.fcopy(['field',stream,collection])
        expr    = [lname, collection]

        while len(expr) > 0:

            key = '.'.join(expr)
            if key in derived:
                info = dict(derived[key])
                info['dexpr'] = dexpr
                return field.Field(lname, info)

            del expr[-1]

#       Check to see if the field name references
#       a native field on file.
#       =========================================

        if not self.ds:
            names = [name, in_stream, in_collection]
            name  = '.'.join([arg for arg in names if arg])
            return field.Field(name, {'dexpr': dexpr})

        fh = self.ds.open(self.request,stream=stream,collection=collection)

        native = fh.field.get(lname,None)

        if native:
            fid  = fh.fileinfo.fid
            info = dict(native)
            info['dexpr'] = dexpr
            info.update(self.file_info(fid))
            return field.Field(lname, info)
        else:
            return name

    def define(self, name):

        if name in self.defined:
            return

        self.defined.append(name)

    def is_defined(self, name):

        if name in self.defined:
            return True

        return False

    def file_info(self, fid):

        time       = self.request['time_dt']
        tm_request = time.strftime("%H:%Mz%d%b%Y").lower()

        self.ds.cmd("set dfile %d"%fid)
        self.ds.cmd("set time " + tm_request)
        self.ds.cmd("query time")

        tm_file = self.ds.rword(1,3).lower()
        if ':' not in tm_file: tm_file = tm_file[0:2] + ':00' + tm_file[2:]

        self.ds.cmd("query dims")
        try:
            tm_index = int(self.ds.rword(5,9))
        except:
            tm_index = 1

        self.file = dict(zip(['fileID','tm_request','tm_file','tm_index'],
                             [fid, tm_request, tm_file, tm_index]))

        return self.file
