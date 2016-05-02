"""

Classes for implementing DC8 navigation and sampling.

"""

from numpy import array
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
        try:
            self.Altitude = self.GPS_Altitude
        except:
            self.Altitude = self.Alt_ft
            self.Longitude = self.Lon
            self.Latitude = self.Lat
        self.Tyme = self.tyme
        self.Hour = array([t.hour for t in self.tyme])
        
class ER2(ICARTT):

    def __init__(self,filename,**kwopts):
        """
        Loads ICARTT file and tweak attributes ti make then look like the
        FlighPlan coordinates.
        """
        ICARTT.__init__(self,filename,**kwops) # load ICART nav file

        self.aircraft = ER2
        self.Tyme = self.tyme
