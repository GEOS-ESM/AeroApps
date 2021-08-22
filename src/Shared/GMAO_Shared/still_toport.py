#!/usr/bin/env python 
from glob import glob
from os.path import basename

def basen(plist):
    blist = []
    for p in plist:
        blist += [basename(p),]
    return blist

    
if __name__ == "__main__":

    py2 = basen(sorted(glob("GMAO_pyobs/pyobs/*.py")))
    py3 = basen(sorted(glob("GMAO_pyobs3/pyobs3/*.py")))

    #print(" 2 --->", py2)
    #print(" 3 --->", py3)

    ignore = [ 'lidar_l2.py', 'oracles.py', 'dragon.py', 'fpl.py',
               'g5_icartt.py', ]
    renamed = [ 'calipso.py', 'calipso_lev2.py' ]

    print("pyobs module")
    for fname in py2:
        if (fname in ignore) or (fname in renamed):
            continue
        if fname not in py3:
            print("   ",fname)
            
        
            
