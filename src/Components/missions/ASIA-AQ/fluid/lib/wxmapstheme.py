
import field
from plot import *
from plotservice import *

class Service(PlotService):

    def __init__(self, *args, **kwargs):

        super(Service,self).__init__(*args, **kwargs)

    def get_capabilities(self, request=None): pass
