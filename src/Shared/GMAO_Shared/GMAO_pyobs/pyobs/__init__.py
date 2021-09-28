"""
   Python interface to several Level 2 observing systems, mostly from EOS.
"""

# ODS, GFIO dependencies
from bits      import BITS
from aeronet   import AERONET_L2
from man       import MAN
from mxd04     import MxD04_L2
from omaeruv   import OMAERUV_L2
from omso2     import OMSO2_L2
from lidar_l2  import LIDAR_L2
from npz       import NPZ
from hsrl      import HSRL
from oracles   import ORACLES
from hsrl_pbl  import HSRL_PBL
from kde       import *
from igbp      import *
from calipso   import CALIPSO_L1p5

# ODS requirement
try:
    from modis import MODIS
except:
    pass

# MAPL dependency
try:
    from icartt    import ICARTT
    from g5_icartt import G5_ICARTT, G5_AIRCRAFT
    from dragon    import DRAGON
except:
    pass

