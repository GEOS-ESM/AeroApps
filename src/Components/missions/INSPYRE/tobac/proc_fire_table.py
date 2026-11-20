#!/usr/bin/env python
from datetime import datetime, timedelta
import os

with open("fire_table.txt","r") as f:
    data = f.readlines()
f.close()

# Write the JSON driver scripts
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
    date_ = date+timedelta(days=1)
    time_ = date_.isoformat()
    tyme  = time[0:4]+time[5:7]+time[8:10]+time[11:13]+"0000"
    tyme_ = time_[0:4]+time_[5:7]+time_[8:10]+time_[11:13]+"0000"
    sedcmd = f"s/LAT/{lat}/g;s/LON/{lon}/g;s/TYM1/{tyme}/g;s/TYM2/{tyme_}/g;s/FIRE/{fire}/g"
    sedcm2 = f'sed "{sedcmd}" my_events.tmp > {fire}_{tyme}.json'
    r = os.system(f"{sedcm2}")
    print(f"{sedcm2}")

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
        tyme  = time[0:4]+time[5:7]+time[8:10]+time[11:13]+"0000"
        cmd = f"sbatch /discover/nobackup/projects/INSPYRE_25/fliu5/shared/mytool/scripts/submit_mytool.sh --events {fire}_{tyme}.json\n"
        f.write(cmd)
f.close()
os.system("chmod 755 run_cases.csh")
