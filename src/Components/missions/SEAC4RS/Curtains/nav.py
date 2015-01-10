"""

Classes for implementing DC8 navigation and sampling.

"""

from pyobs import ICARTT

class DC8(ICARTT):
    """
    Tinker ICARTT class for DCO urtain plots.
    """
    def __init__(self,filename,**kwopts):
        """
        Loads ICARTT file and tweak attributes ti make then look like the
        FlighPlan coordinates.
        """
        ICARTT.__init__(self,filename,**kwopts) # load ICART nav file

        self.aircraft = DC8
        self.Altitude = self.GPS_Altitude
        self.Tyme = self.tyme

class ER2(ICARTT):

    def __init__(self,filename,**kwopts):
        """
        Loads ICARTT file and tweak attributes ti make then look like the
        FlighPlan coordinates.
        """
        ICARTT.__init__(self,filename,**kwops) # load ICART nav file

        self.aircraft = ER2
        self.Tyme = self.tyme
