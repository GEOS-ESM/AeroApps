ncl 'file_in="/discover/nobackup/projects/gmao/nca/firexaq/native/FireWork/Y2019/M08/D27/H00/2019082700_065"' 'file_out="/discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4.tmp1.nc"' /home/dao_ops/jardizzo/FLUID/firex-aq/bin/firework_regrid_latlon.ncl

firework_lv1_2p.py -v -o /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4 /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4.tmp1.nc

firework_lv2_2p.py -v -o /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4.tmp2.nc4 /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4.tmp1.nc

ncks -h -A /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4.tmp2.nc4 /discover/nobackup/projects/gmao/nca/pub/firexaq/data/FireWork/Y2019/M08/D27/H00/FireWork.20190829_17.nc4
