#!/usr/bin/env python
from model_utils.eval import budget
from model_utils.post import monthly_mod as samplemodis
from optparse   import OptionParser
import os

def evaldriver(expid):
    make_budget(expid)
    make_sample(expid)
    make_sample(expid,sat="MOD04")
    return

def make_sample(expid,path="/home/pcolarco/geos_aerosols/pcolarco",
                collection="inst2d_hwl_x", sat="MYD04",
                obspath="/home/pcolarco/geos_aerosols/pcolarco",
                year=2019):
    varlist = ["TOTEXTTAU550","DUEXTTAU550","SSEXTTAU550","OCEXTTAU550",
               "BCEXTTAU550","SUEXTTAU550","NIEXTTAU550","TOTSCATAU550",
               "TOTANGSTR","PM","PM25"]
    for mm in range(1,13):
        outfile=samplemodis.monthlymodel(year,mm,sat,expid,
                    expdir=f"{path}/{expid}/holding/{collection}/{year:04d}{mm:02d}/",
                    obsdir=f"{obspath}/{sat}/Level3/Y{year:04d}/M{mm:02d}/",
                    varlist=varlist)
    return

def make_budget(expid,path="/home/pcolarco/geos_aerosols/pcolarco",
                collection="tavg2d_aer_x", year=2019):
    filelist = f"{path}/{expid}/{collection}/{expid}.{collection}.monthly.{year}??.nc4"
    # Should have checks that files exist and get 12 full months
    b = budget.AEROSOL(filelist,verbose=True)
#   Make an output directory
    directory_name = expid
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    b.plotAll(expid,verbose=True)
    b.printBudget(expid,verbose=True)
    

    return

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog expid",
                          version='xxx' )
    (options, args) = parser.parse_args()

#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 1:
        expid      = args[0]
    else:
        parser.error("must have 1 argument: expid")

    evaldriver(expid)
    
