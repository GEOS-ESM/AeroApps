import copy

class Encoder(object):

    def __init__(self, service):

        self.service        = copy.copy(service)
        self.service.config = copy.copy(self.service.config)

        keymap = self.service.config.get('keymap',{})
        self.service.config.update(keymap)

    def encode(self, request):

        plots = self.service.get_plot(request, passive=True)

        str = ''
        for plot in plots: str += repr(plot.cmds)

        return str
