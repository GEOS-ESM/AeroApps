"""
Simple climatological albedo interpoleation class.
"""

from grads import GrADS
from numpy import minmax

class ALBEDO(object):

    def __init__(self,filename='albedo_clim.ctl'):
        """Start grads and open albedo file."""
        
        self.ga = GrADS(Echo=0,Window=False)
        self.ga('open %s'%filename)

    def interp(self,lons,lats,gatime):
        """Given date/time, interpolate albedo to lons, lats."""

        # Restrict domain to save memory and I/O
        # --------------------------------------
        self.ga('set lon %f %f'%(min(lons),max(lons)))
        self.ga('set lat %f %f'%(min(lats),max(lats)))

        # Set the closest time
        # --------------------
        self.ga('set time %s'%gatime)
                
        # Interpolate to lat/lon at closest time
        # --------------------------------------
        albedo, levs = ga.interp('albedo', lons, lats)
        return albedo.data
    
#--

if __name__ == "__main__":

    from anet import LAND
    
    a = ALBEDO()

    modl = LAND('SUPER2_combo.Terra.csv',csvVersion=2)


    
