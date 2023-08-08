import re
import sys
import copy
import yaml
import json

novalue    = object()
config_cache = {}

class Config(dict):

    def __init__(self, *args, **kw):

        self.registry = {}
        super(Config,self).__init__(*args, **kw)

    def find(self, path, name, depth=-1, cfg=None):

        result = []

        if cfg is None:
            cfg = self.follow(path)

        for key in cfg:

            apath = path + [key]

            if key == name:
                result.append(apath)

            if depth == 0:
                continue

            if self.ispartition(cfg[key]):
                result += self.find(apath,name,depth,cfg=cfg[key])
            elif isinstance(cfg[key], dict):
                result += self.find(apath,name,depth-1,cfg=cfg[key])

        return result

    def follow(self, paths):

        cfg   = self
        apath = ''

        for path in paths:

            if path == '/':
                continue

            apath += '/' + str(path)

            if path not in cfg:
                raise KeyError('Config.follow: "' + apath +
                               '" No such file or directory')

            if isinstance(cfg[path],dict):
                cfg = cfg[path]
            else:
                raise KeyError('Config.follow: "' + apath +
                                          '" is not a directory')

        return cfg

    def get_values(self, pathname, default=None):

        return self.get_items(pathname, default).values()

    def get_keys(self, pathname, default=None):

        return self.get_items(pathname, default).keys()

    def get_items(self, pathname, default=None, flat=True, hide=True):

        items     = {}
        pathnames = pathname

        if not isinstance(pathname,list):
            return items

        if not isinstance(pathname[0],list):
            pathnames = [pathname]
        
        for pn in pathnames:

            path = pn[0:-1]
            name = pn[-1]
            key  = pn[-2]
            cfg  = self.follow(path)

            if hide and cfg.get('hide','no') == 'yes': continue

            hash = items
            if not flat: hash = self.mkpath(hash,path[1:-1])

            if name not in cfg:
                hash[key] = default
            else:
                hash[key] = cfg[name]

        return items

    def get_config(self, pathname, default=novalue):

        pathname = self.expand(pathname)

        try:
            cfg = self.follow(pathname[0:-1])
        except KeyError:
            if default is novalue:
                raise
            else:
                return default

        if default is novalue:
            return cfg[pathname[-1]]
        else:
            return cfg.get(pathname[-1],default)

    __call__ = get_config

    def expand(self, paths):

        pathname = []
        if not isinstance(paths,list):
            paths = [paths]

        for path in paths:

            if isinstance(path,basestring):
                pathname += path.split('/')
            else:
                pathname.append(path)

        return pathname

    def mkpath(self, root, path):

        for dir in path:

            if dir not in root:
                root[dir] = {}

            root = root[dir]

        return root

    def read(self, file):

        if file in config_cache: return config_cache[file]

        with open(file, 'r') as ymlfile:
            config_cache[file] = self.copy_yaml(yaml.load(ymlfile))

        return config_cache[file]

    def readJSON(self, file):

        if file in config_cache: return config_cache[file]

        with open(file, 'r') as jsonfile:
            config             = json.load(jsonfile)
            config_cache[file] = self.deserialize(config)

        return config_cache[file]

    def mount(self, cfg, root=None):

        hash = self

        if root is not None:

            for dir in root.split('/'):

                if dir == '/': continue
                if not dir: continue

                if dir not in hash:
                    hash[dir] = {}
                elif not isinstance(hash[dir], dict):
                    hash[dir] = {}

                hash = hash[dir]

        self.overlay(hash,cfg)

    def ispartition(self, dir):

        if not isinstance(dir,dict):
            return False

        result = [key for key in dir.keys() if not isinstance(dir[key],dict)]

        if result:
            return False

        return True

    def fcopy(self, path):

        flat_list = {}
        hash      = self

        for dir in path:
            hash = hash.get(dir,{})
            flat_list.update(hash)
            if dir in flat_list:
                del flat_list[dir]

        return flat_list

    def fdcopy(self, path):
        return copy.deepcopy(self.fcopy(path))

    def overlay(self, hash1, hash2):

        for key2 in hash2:

            if key2 not in hash1:
                if isinstance(hash2[key2], dict):
                    hash1[key2] = copy.deepcopy(hash2[key2])
                else:
                    hash1[key2] = hash2[key2]
            elif isinstance(hash2[key2],dict) and isinstance(hash1[key2],dict):
                self.overlay(hash1[key2], hash2[key2])
            else:
                hash1[key2] = hash2[key2]

    def serialize(self, hash):

        for key in hash:

            if isinstance(hash[key],dict):
                self.serialize(hash[key])

            try:
                json.dumps(hash[key])
            except TypeError:
                hash[key] = hash[key].__module__ + '.' + \
                            hash[key].__class__.__name__ + \
                            ' []'

    def deserialize(self, hash):

        for key in hash:

            if isinstance(hash[key],dict):
                self.deserialize(hash[key])

            if self.is_object_string(hash[key]):

                object      = hash[key]
                module_name = object.split('.')[0]
                class_name  = object.split('.')[1].split()[0]

                module      = __import__(module_name)
                class_      = getattr(module, class_name)
                hash[key]   = class_()

        return hash

    def is_object_string(self, value):

        if not isinstance(value, basestring): return False
        match = re.match(r'\w+\.\w+ \[\]', value)
        if match: return True

        return False

    def copy_yaml(self, hash):

        for key in hash:

            new_key = str(key)

            if not isinstance(key, basestring):
                hash[new_key] = hash[key]
                del hash[key]

            if isinstance(hash[new_key], dict):

                if id(hash[new_key]) in self.registry:
                    hash[new_key] = dict(hash[new_key])
                    self.registry[id(hash[new_key])] = 1
                else:
                    self.registry[id(hash[new_key])] = 1

                self.copy_yaml(hash[new_key])

        return hash

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class UsageError(Error):
    """Exception raise for errors in the input."""
    def __init__(self, msg):
        self.msg = msg
