#!/usr/bin/env python
"""
Select dee blue dust for a range of years.
"""
import sys
sys.path.insert(1,'../')

from deep import *

#....................................................................

if __name__ == "__main__":

    path = '/nobackup/MODIS/006/Level2'
    years = (2012, 2014)

    grid_dust(path,years[0],years[1],dirn='../DEE')
    
