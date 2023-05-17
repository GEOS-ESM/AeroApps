#!/usr/bin/env python3

import idlsave 
#from netCDF4 import Dataset,date2num 
from pyobs import NPZ
import numpy.ma as ma
from numpy import linspace, array, loadtxt, ones
from datetime import datetime
import netCDF4 as nc
MISSING = 1.e+15 


if __name__ == "__main__":

    # input file with lat lon
    #-----------
    coord_file = 'tempo_coords.npz'
    geo = NPZ(coord_file)
    ns, ew = geo.ns_npix, geo.ew_nscan
    ns_e, ew_e = ns+1, ew+1

    # input file with sensor angles
    #-----------
    sensor_file = 'vza_vaza.txt'
    angles = loadtxt(sensor_file)

    view_zen_ang = MISSING * ones((ns,ew))
    view_azi_ang = MISSING * ones((ns,ew))
    for i_ns, i_ew, vza,vaa in zip(angles[:,1],angles[:,0],angles[:,4],angles[:,5]):
           view_zen_ang[i_ns-1,i_ew-1]= vza
           view_azi_ang[i_ns-1,i_ew-1]= vaa  # -1 because indice starts at 1 not 0

    
#def hold():
    # output file
    #------------
    file = nc.Dataset('tempo.la0.invariant_wsens_ang.nc4','w',format='NETCDF4')
    
    # Set global attributes
    # ---------------------
   
    file.institution = 'NASA/Goddard Space Flight Center'
    file.source = 'Global Model and Assimilation Office'
    file.history = 'Created from tempo_geo_footprint_100W.sav'
    file.references = 'n/a'
    file.comment = 'This file contains TEMPO geolocation informations'
    file.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    file.Conventions = 'CF'
    file.sat_lat = geo.satlat
    file.sat_lon = geo.satlon
    file.sat_alt = geo.satalt
    file.Earth_radius = geo.re
    file.Earth_flattening = geo.fe

    # Create dimensions
    #-------------------
    file.createDimension('ns',ns)
    file.createDimension('ew',ew)
    file.createDimension('time',1)

    file.createDimension('ns_e',ns_e)
    file.createDimension('ew_e',ew_e)
    file.createDimension('pix_edge',4)
    

    # create fake lat, lon for GRADS compatibility
    # ----------------
    
    ns = file.createVariable('ns','f8',('ns'),zlib=False)
    ns.long_name = 'pseudo latitude'
    ns.units = 'degrees_north'
    ns[:] = linspace(18,58,2000)   # assumed lat min = 18 and lat max =58 (ns = 2000)

    ew = file.createVariable('ew','f8',('ew',),zlib=False)
    ew.long_name = 'pseudo longitude'
    ew.units = 'degrees_east'
    ew[:] = linspace(-130,-70,1500)  # assumed lon min = -174 and lon max = -70 (ew=1500)

    time = file.createVariable('time','f8',('time'),zlib=False)
    t = datetime(2006,0o6,0o1,12,00)   # just an example
    time_increment          = int('{hours}{minutes}{seconds}'.format(hours='01', minutes='00',
 seconds='00'))
    begin_date              = int(t.strftime('%Y%m%d'))
    begin_time              = int(t.strftime('%H%M%S'))
    time.long_name          = 'time'
    time.units            = 'hours since {:%Y-%m-%d %H:%M:%S}'.format(t)
    time.time_increment   = array(time_increment)
    time.begin_date       = array(begin_date)
    time.begin_time       = array(begin_time)
    time[:] = nc.date2num(t, units=time.units,calendar='gregorian')
 

    
    # Create variables
    # ---------------- 

    pix_size = file.createVariable('pix_size','f8',('pix_edge','ns','ew',),zlib=False)
    pix_size.long_name = 'pixel size --> sizes of each egde in km in clockwise orientation starting bottom'
    pix_size.units = 'unknown'
    pix_size_masked = ma.masked_equal(geo.pixsize,0.0)   # masked 0.0 values
    ma.set_fill_value(pix_size_masked,MISSING)     # replace default fill value by MISSING value
    pix_size[:] = pix_size_masked

    c_lat = file.createVariable('clat','f8',('ns','ew'),zlib=False)
    c_lat.long_name = 'pixel center latitude'
    c_lat.units = 'degrees_north'
    c_lat.missing_value = MISSING
    c_lat_masked = ma.masked_invalid(geo.clat)    # masked nan values
    ma.set_fill_value(c_lat_masked,MISSING)       # replace default fill value by MISSING value
    c_lat[:] = c_lat_masked

    c_lon = file.createVariable('clon','f8',('ns','ew',),zlib=False)
    c_lon.long_name = 'pixel center longitude'
    c_lon.units = 'degrees_east'
    c_lon.missing_value = MISSING
    c_lon_masked = ma.masked_equal(geo.clon,0.0)   # masked 0.0 values
    ma.set_fill_value(c_lon_masked,MISSING)     # replace default fill value by MISSING value
    c_lon[:] = c_lon_masked

    e_lat = file.createVariable('elat','f8',('ns_e','ew_e'),zlib=False)
    e_lat.long_name = 'latitude at pixel edge of southeast corner of pixel'
    e_lat.units = 'degrees_north'
    e_lat.missing_value = MISSING
    e_lat_masked = ma.masked_invalid(geo.elat)    # masked nan values
    ma.set_fill_value(e_lat_masked,MISSING)       # replace default fill value by MISSING value
    e_lat[:] = e_lat_masked

    e_lon = file.createVariable('elon','f8',('ns_e','ew_e',),zlib=False)
    e_lon.long_name = 'longitude at pixel edge of southeast corner of pixel'
    e_lon.units = 'degrees_east'
    e_lon.missing_value = MISSING
    e_lon_masked = ma.masked_equal(geo.elon,0.0)   # masked 0.0 values
    ma.set_fill_value(e_lon_masked,MISSING)     # replace default fill value by MISSING value
    e_lon[:] = e_lon_masked

    vza = file.createVariable('vza','f8',('ns','ew',),zlib=False)
    vza.long_name = 'sensor zenith angle'
    vza.units = 'degrees'
    vza.missing_value = MISSING
    vza_masked = ma.masked_invalid(view_zen_ang)    # masked nan values
    ma.set_fill_value(vza_masked,MISSING)          # replace default fill value by MISSING value
    vza[:] = vza_masked
    
    vaa = file.createVariable('vaa','f8',('ns','ew',),zlib=False)
    vaa.long_name = 'sensor azimuth angle'
    vaa.units = 'degrees'
    vaa.missing_value = MISSING
    vaa_masked = ma.masked_invalid(view_azi_ang)   # masked nan values
    ma.set_fill_value(vaa_masked,MISSING)          # replace default fill value by MISSING value
    vaa[:] = vaa_masked
  

