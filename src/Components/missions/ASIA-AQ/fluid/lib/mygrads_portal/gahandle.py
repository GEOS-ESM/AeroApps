#try:
#    from org.opengrads.interfaces import GaHandleObject
#except:
#    GaHandleObject = object

__version__ = '1.1.0'

class GaHandle(object):
    """
    A simple container class to collect output for query() operations.
    """
    def __init__ (self, name):
        self.name = name
    
    def get (self, key):
        return self.__dict__[key]

    def set (self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        repr  = ''
        repr += object.__repr__(self)
        repr += '\n'
        repr += '==================================\n'
        repr += '  Handle name = %s\n' % self.name
        repr += '==================================\n'
        
        longestKeyLen = max( map(len, self.__dict__.keys()) )
        
        for k, v in self.__dict__.iteritems():
            repr += '%s = %s' % (k.rjust(longestKeyLen), v)
            repr += '\n'
            
        return repr
    
