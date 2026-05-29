#!/usr/bin/env python
from datetime import datetime, timedelta
import os

with open("fire_table.txt","r") as f:
    data = f.readlines()
f.close()

# Write the profile initialization script
with open("initialize_cases.csh","w") as f:
    f.write("#!/bin/tcsh -f\n")
    for line in data:
        items = line.split()
        if(items[0][0] == "#"):
            continue
        fire  = items[0]
        lat   = items[1]
        lon   = items[2]
        time  = items[3]
        alt   = items[4]
        if(alt != "Y"):
            cmd = f"./initialize_profile.csh {fire}_{time} {lon} {lat} {alt}\n"
            f.write(cmd)
        else:
            cmd = f"./initialize_profile.csh {fire}_{time}_10km {lon} {lat} 10.\n"
            f.write(cmd)
            cmd = f"./initialize_profile.csh {fire}_{time}_12km {lon} {lat} 12.\n"
            f.write(cmd)
            cmd = f"./initialize_profile.csh {fire}_{time}_14km {lon} {lat} 14.\n"
            f.write(cmd)
f.close()
os.system("chmod 755 initialize_cases.csh")



# Write the run script
with open("run_cases.csh","w") as f:
    f.write("#!/bin/tcsh -f\n")
    for line in data:
        items = line.split()
        if(items[0][0] == "#"):
            continue
        fire  = items[0]
        lat   = items[1]
        lon   = items[2]
        time  = items[3]
        alt   = items[4]
        date  = datetime.strptime(time,"%Y-%m-%dT%H:%M:%S")
        date_ = date+timedelta(days=3)
        time2 = date_.isoformat()
        if(alt != "Y"):
            cmd = f"./run_kinematic.csh {fire}_{time} \"{time}\" \"{time2}\"\n"
            f.write(cmd)
        else:
            cmd = f"./run_kinematic.csh {fire}_{time}_10km \"{time}\" \"{time2}\"\n"
            f.write(cmd)
            cmd = f"./run_kinematic.csh {fire}_{time}_12km \"{time}\" \"{time2}\"\n"
            f.write(cmd)
            cmd = f"./run_kinematic.csh {fire}_{time}_14km \"{time}\" \"{time2}\"\n"
            f.write(cmd)
f.close()
os.system("chmod 755 run_cases.csh")

# Write the plot script
with open("plot_cases.csh","w") as f:
    f.write("#!/bin/tcsh -f\n")
    for line in data:
        items = line.split()
        if(items[0][0] == "#"):
            continue
        fire  = items[0]
        lat   = items[1]
        lon   = items[2]
        time  = items[3]
        alt   = items[4]
        if(alt != "Y"):
            cmd = f"./plot_parcel.py parcel_traj.{fire}_{time}.nc\n"
            f.write(cmd)
        else:
            cmd = f"./plot_parcel_forecast.py parcel_traj.{fire}_{time}\n"
            f.write(cmd)
f.close()
os.system("chmod 755 plot_cases.csh")

