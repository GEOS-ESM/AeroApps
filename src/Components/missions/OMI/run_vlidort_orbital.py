from MAPL.ShaveMantissa_ import shave32
from   mieobs  import  getAOPvector, getEdgeVars
import sys 
import numpy as np
import h5py
import glob 
import os
import re
import pandas as pd 
from scipy.interpolate import interp1d 
from pyobs.omno2    import OMNO2_L2, granules
from netCDF4 import Dataset
import datetime
from pyobs.nc4ctl   import NC4ctl
from MAPL           import Config, eta
from MAPL.constants import *

VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
VNAMES_NI = ['NO3AN1','NO3AN2','NO3AN3']

AERNAMES = VNAMES_SU  + VNAMES_OC + VNAMES_BC + VNAMES_DU + VNAMES_SS
AERNAMES_NOBC = VNAMES_SU + VNAMES_OC + VNAMES_DU + VNAMES_SS
SDS_AER = ['O3','RH','DELP','PS','H','T'] + AERNAMES
SDS_GEO    = ['clon','clat']


class SampleVar(object):
    """
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name


class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation

def read_l1b(l1b_file,band):
    geo_grp = '/'+band+'_RADIANCE/STANDARD_MODE/GEODATA/'
    obs_grp = '/'+band+'_RADIANCE/STANDARD_MODE/OBSERVATIONS/'
    
    f = Dataset(l1b_file,'r')
    omi_lat = f[geo_grp+'latitude'][0]
    omi_lon = f[geo_grp+'longitude'][0]
    omi_sec = f[obs_grp+'time_TAI93'][0]
    omi_vza = f[geo_grp+'viewing_zenith_angle'][0]
    omi_sza = f[geo_grp+'solar_zenith_angle'][0]
    omi_vaa = f[geo_grp+'viewing_azimuth_angle'][0]
    omi_saa = f[geo_grp+'solar_azimuth_angle'][0]

    omi_time = np.repeat(np.array([datetime.datetime(1993,1,1) + datetime.timedelta(seconds=i) for  i in omi_sec])[:,np.newaxis],60,axis=1)

    return omi_lat, omi_lon, omi_time, omi_sec, omi_sza, omi_vza, omi_saa, omi_vaa

def getVars(rcFile):
    """
    Parse reource file, create variable dictionary with relevant
    metadata.
    """
    cf = Config(rcFile)
    Vars = dict()
    levUnits = 'none'
    levs = []
    for V in cf.keys():
        path = cf(V)
        f = Open(path)
        varList = []
        for v in V.split(','):
            var = f.Vars[v.strip()]
            if var.km>0:
                levUnits = var.levunits
                levs = var.levs
            varList += [var,]
        Vars[path] = varList

    return (Vars, levs, levUnits)

def Open(filename,doNC4ctl=True):
    """
    Uses GFIO or GFIOctl to open either a NetCDF-4 or a control file.
    Very heuristic.

    TO DO: Remove GFIO dependencies, using NC4ctl instead.

    """
    print(filename)
    from gfio import GFIOurl, GFIOctl
    name, ext = os.path.splitext(filename)
    # Open the GFIO dataset
    # ---------------------
    f = None
    if 'HTTP://' == name[:7].upper():
        f = GFIOurl(filename)
        f.lower = True # force variable names to be lower case when sampling.
    elif ext.upper() in ('.NC4','.NC','.HDF','.H5'):
        f = GFIOurl(filename)
        f.lower = False
    else:
        f = GFIOctl(filename)
        f.lower = False
    # Create variable dictionary
    # --------------------------
    Vars = dict()
    for i in range(len(f.vname)):
        if f.lower:
            v = f.vname[i].upper()
        else:
            v = f.vname[i]
        var = SampleVar(v)
        #var.title = f.vtitle[i]
        var.km = f.kmvar[i]
        if var.km>0:
            var.levunits = f.levunits[:]
            var.levs = f.levs[:]
        try:
            var.units = f.vunits[i]
        except:
            var.units = 'unknown'  # opendap currently lacks units
        Vars[v] = var

    f.Vars = Vars

    if doNC4ctl:
        f.nc4 = NC4ctl_(filename)

    return f

class HOLDER(object):
    def __init__(self,aop):
        for sds in aop.SDS_AER:
            print(sds,np.shape(aop.__dict__[sds]))
            if sds == 'PS':
                self.__dict__[sds] = aop.__dict__[sds][:]
            else:
                self.__dict__[sds] = aop.__dict__[sds][:,:]

class MERRA2(object):

    def __init__(self,Vars,sat_lon,sat_lat,sat_time,line_1,line_2,row_1,row_2,channel,rcFile):
        for path in Vars:
            print " <> opening "+path
            g = Open(path)
            for var in Vars[path]:
                
                if g.lower:
                    name = var.name.lower() # GDS always uses lower case
                else:
                    name = var.name


                Z = g.sample(name,sat_lon[line_1:line_2,row_1:row_2].ravel(),sat_lat[line_1:line_2,row_1:row_2].ravel(),sat_time[line_1:line_2,row_1:row_2].ravel(),Transpose=False,squeeze=True,Verbose=False,algorithm='linear')
                self.__dict__[name] = Z.squeeze()


                
                
        self.SDS_AER  = SDS_AER
        self.SDS_GEO  = SDS_GEO
        self.AERNAMES = AERNAMES
        self.AERNAMES_NOBC = AERNAMES_NOBC
        self.VNAMES_BC = VNAMES_BC

        

        aero = HOLDER(self)

        nmom = 300

        tau,ssa,g,pmom = getAOPvector(aero,channel,vnames=self.AERNAMES,Verbose=False,rcfile='Spectral_Scaling.rc',nMom=nmom)
        tot_abs = tau*(1-ssa)
        tot_scat = tau*ssa

        tau_bc,ssa_bc,z,pmom_bc = getAOPvector(aero,channel,vnames=self.VNAMES_BC,Verbose=False,rcfile='Spectral_Scaling.rc',nMom=nmom)
        bc_abs = tau_bc*(1-ssa_bc)
        ssa = tot_scat/(tot_scat+(tot_abs-bc_abs))
        pmom[:,:,:,:,1] = -1.*pmom[:,:,:,:,1]
        pmom[:,:,:,:,3] = -1.*pmom[:,:,:,:,3]

        

        self.__dict__['ssa'] = ssa
        self.__dict__['tau'] = tau
        self.__dict__['pmom'] = pmom
        self.__dict__['g'] = g
        
        for key in self.__dict__.keys():
            print(key,np.shape(self.__dict__[key]))

            if ('SDS' in key) | ('AERNAMES' in key) | ('VNAMES' in key):
                continue
            if np.array(self.__dict__[key]).ndim == 1:
                self.__dict__[key] = np.array(self.__dict__[key]).reshape(line_2-line_1,row_2-row_1)
            elif np.array(self.__dict__[key]).ndim == 2:
                self.__dict__[key] = np.array(self.__dict__[key]).reshape(72,line_2-line_1,row_2-row_1)
            elif np.array(self.__dict__[key]).ndim == 3:
                self.__dict__[key] = np.array(self.__dict__[key]).reshape(72,len(channel),line_2-line_1,row_2-row_1)
            elif np.array(self.__dict__[key]).ndim == 5:
                self.__dict__[key] = np.array(self.__dict__[key]).reshape(72,len(channel),line_2-line_1,row_2-row_1,300,6)

        #self.__dict__['ssa'] = ssa.reshape(72,len(channel),line_2-line_1,row_2-row_1)
        #self.__dict__['tau'] = tau.reshape(72,len(channel),line_2-line_1,row_2-row_1)
        #self.__dict__['pmom'] = pmom.reshape(72,len(channel),line_2-line_1,row_2-row_1,300,6)
        #self.__dict__['g'] = g.reshape(72,len(channel),line_2-line_1,row_2-row_1)



l1brug_files = glob.glob('OML1BRUG/*/*/*/*')


rcFile = 'pace_sampler.rc'
Vars, levs, levUnits = getVars(rcFile)
line_1 = int(sys.argv[1]) -1
line_2 = int(sys.argv[2]) 
row_1 = int(sys.argv[3]) -1
row_2 = int(sys.argv[4])

for filename in l1brug_files:


    date_str = re.findall('[0-9]{4}m[0-9]{4}t[0-9]{4}-o[0-9]{6}',filename)[0]
    yr = date_str[0:4]
    mn = date_str[5:7]
    day = date_str[7:9]

    orb = re.findall('-o[0-9]{6}',filename)[0]

    m2_waves = np.arange(270.,570.,20.)


    omi_lat, omi_lon, omi_time, omi_sec, omi_sza, omi_vza, omi_saa, omi_vaa = read_l1b(filename,'BAND2')

    omi_raa = np.abs(omi_saa-omi_vaa) 
    omi_raa[omi_raa > 180.] = 360.-omi_raa[omi_raa > 180.]
    omi_raa = 180. - omi_raa 

    #omi_lat[1126,50] = 20.991205
    #omi_lon[1126,50] = -157.11382
    #omi_vza[1126,50] = 44.245296
    #omi_sza[1126,50] = 46.23659
    #omi_raa[1126,50] = 127.38314

    m2_data =MERRA2(Vars,omi_lon,omi_lat,omi_time,line_1,line_2,row_1,row_2,m2_waves,rcFile)



    to3_file = glob.glob('OMTO3/'+str(yr)+'/'+str(mn).zfill(2)+'/'+str(day).zfill(2)+'/OMI-Aura*'+orb.replace('-o0','-o')+'*.he5')[0]
    f = h5py.File(to3_file,'r')
    omi_ozone = f['/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/ColumnAmountO3'][:]
    omi_ozone[omi_ozone == -2**100] = np.nan
    f.close()

    m2_ozone = np.nansum((1.0e6)*m2_data.__dict__['O3']*(28.95949/47.998)*0.0078914*m2_data.__dict__['DELP'],axis=0)
    m2_layer_ozone = (1.0e6)*m2_data.__dict__['O3']*(28.95949/47.998)*0.0078914*m2_data.__dict__['DELP']
    ozone_scale = omi_ozone[line_1:line_2,row_1:row_2]/m2_ozone


    scaled_ozone = m2_data.__dict__['O3'] * ozone_scale
    scaled_layer_ozone = m2_layer_ozone* np.repeat(ozone_scale[np.newaxis,:,:],72,axis=0)
    aeruv_ind, = np.where(m2_waves == 390)[0]
    
    aeruv_file = glob.glob('OMAERUV/'+str(yr)+'/'+str(mn).zfill(2)+'/'+str(day).zfill(2)+'/OMI-Aura*'+orb+'.nc')[0]
    f = h5py.File(aeruv_file,'r')
    uv_aod = f['/SCIDATA/FinalAerosolOpticalDepth'][:,:,1]
    f.close()
    aeruv_ind, = np.where(m2_waves == 390)[0]
    uv_aod[uv_aod == -2**100] = np.nan

    print(m2_data.__dict__.keys())
    aod_factor = np.array(uv_aod[line_1:line_2,row_1:row_2]/np.nansum(m2_data.tau[:,aeruv_ind,:,:],axis=0),dtype=np.float)
    aod_factor = np.repeat(aod_factor[np.newaxis,:,:],len(m2_waves),axis=0)
    aod_factor = np.repeat(aod_factor[np.newaxis,:,:,:],72,axis=0)

    
    scaled_aod = m2_data.tau * aod_factor

    gler_file = glob.glob('OMGLER/'+str(yr)+'/'+str(mn).zfill(2)+'/'+str(day).zfill(2)+'/OMI-Aura*'+orb.replace('-o0','-o')+'*.he5')[0]
    f= h5py.File(gler_file,'r')
    omi_wsp = f['/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/WindSpeed'][:]
    omi_chl = f['/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/ChlorophyllConcentration'][:]
    f.close()

    for line in range(line_1,line_2):
        for row in range(row_1,row_2):
            if (omi_sza[line,row] < 0) | (omi_sza[line,row] > 88) | np.isnan(uv_aod[line,row]):
                continue
            
            print('LINE: ',line,' ROW: ',row)
            inp_filename = 'Input_Files/OMI_PACE_UV2_Input_'+str(yr)+'m'+str(mn).zfill(2)+str(day).zfill(2)+orb+'_Line_'+str(line).zfill(4)+'_XTrack_'+str(row).zfill(2)+'.h5'
            f_inp = h5py.File(inp_filename,'w')
            f_inp.create_dataset('Lat',data=omi_lat[line,row],dtype=np.float64)
            f_inp.create_dataset('Lon',data=omi_lon[line,row],dtype=np.float64)
            f_inp.create_dataset('Time',data=omi_sec[line])
            for m2_key in ['DELP','H','PS','T']:
                if np.array(m2_data.__dict__[m2_key]).ndim == 2:
                    f_inp.create_dataset(m2_key,data=m2_data.__dict__[m2_key][line-line_1,row-row_2])
                elif np.array(m2_data.__dict__[m2_key]).ndim == 3:
                    f_inp.create_dataset(m2_key,data=m2_data.__dict__[m2_key][:,line-line_1,row-row_2])
                elif np.array(m2_data.__dict__[m2_key]).ndim == 4:
                    f_inp.create_dataset(m2_key,data=m2_data.__dict__[m2_key][:,:,line-line_1,row-row_2])

            vl_waves = np.arange(310,375,5)
            
            f_inp.create_dataset('LayerOzone',data=scaled_ozone[:,line-line_1,row-row_1])


            f_tau = interp1d(m2_waves,scaled_aod[:,:,line-line_1,row-row_1],axis=1)
            vl_tau_scaled = f_tau(vl_waves)

            f_tau = interp1d(m2_waves,m2_data.tau[:,:,line-line_1,row-row_1],axis=1)
            m2_tau_scaled = f_tau(vl_waves)

            f_ssa = interp1d(m2_waves,m2_data.ssa[:,:,line-line_1,row-row_1],axis=1)
            vl_ssa = f_ssa(vl_waves)

            print(np.shape(m2_waves),np.shape(m2_data.pmom[:,:,line-line_1,row-row_1,:,:]))
            f_pmom = interp1d(m2_waves,m2_data.pmom[:,:,line-line_1,row-row_1,:,:],axis=1)
            vl_pmom = f_pmom(vl_waves)

            f_inp.create_dataset('m2_waves',data=m2_waves)
            f_inp.create_dataset('nwv_merra2',data=len(m2_waves))
            f_inp.create_dataset('tau',data=vl_tau_scaled)
            f_inp.create_dataset('ssa',data=vl_ssa)
            f_inp.create_dataset('pmom',data=vl_pmom)

            f_inp.create_dataset('vl_waves',data=vl_waves)
            f_inp.create_dataset('nwv_vl',data=len(vl_waves))
            f_inp.create_dataset('Year',data=int(yr),dtype=np.int16)
            f_inp.create_dataset('Month',data=int(mn),dtype=np.int16)
            f_inp.create_dataset('Day',data=int(day),dtype=np.int16)
            f_inp.create_dataset('Orbit',data=int(orb.replace('-o','')),dtype=np.int16)
            f_inp.create_dataset('Line',data=line,dtype=np.int16)
            f_inp.create_dataset('Row',data=row,dtype=np.int16)
            f_inp.create_dataset('SZA',data=omi_sza[line,row],dtype=np.float64)
            f_inp.create_dataset('VZA',data=omi_vza[line,row],dtype=np.float64)
            f_inp.create_dataset('RAA',data=omi_raa[line,row],dtype=np.float64)
            f_inp.create_dataset('WSP',data=omi_wsp[line,row],dtype=np.float64)
            f_inp.create_dataset('CHL',data=omi_chl[line,row],dtype=np.float64)


            f_inp.close()

            if 'VIS' in inp_filename:
                #vl_waves = np.arange(350,505,5)
                out_filename = inp_filename.replace('Input_Files/','Output_Files/').replace('OMI_PACE_VIS_Input','OMI_PACE_Vis_V1').replace('.h5','.dat')
                
            elif 'UV2' in inp_filename:
                #continue
                #vl_waves = np.arange(310,375,5)
                out_filename = inp_filename.replace('Input_Files/','Output_Files/').replace('OMI_PACE_UV2_Input','OMI_PACE_UV2_V1').replace('.h5','.dat')


            os.system("./vlidort_test_gas_aerosol_retrieval_single.exe "+inp_filename+" "+out_filename)

            vl_data = pd.read_csv(out_filename,sep='\s*[,]\s*')

            h5_filename = out_filename.replace('.dat','.h5')
            f = h5py.File(h5_filename,'w')
    
            f.create_dataset('BOAFdn_Lw',data=vl_data['ES'],dtype=np.float32)
            f.create_dataset('CHL',data=omi_chl[line,row],dtype=np.float32)
            f.create_dataset('AerosolLayerTau_MERRA2',data=m2_tau_scaled,dtype=np.float32)
            f.create_dataset('AerosolLayerTau_Scaled',data=vl_tau_scaled,dtype=np.float32)
            f.create_dataset('AerosolSSA',data=vl_ssa,dtype=np.float32)
            f.create_dataset('DELP',data=m2_data.__dict__['DELP'][:,line-line_1,row-row_2],dtype=np.float32)
            f.create_dataset('Year',data=int(yr),dtype=np.int)
            f.create_dataset('Month',data=int(mn),dtype=np.int)
            f.create_dataset('Day',data=int(day),dtype=np.int)
            f.create_dataset('Orbit',data=int(orb.replace('-o','')),dtype=np.int)
            f.create_dataset('Line',data=int(line),dtype=np.int)
            f.create_dataset('Row',data=int(row),dtype=np.int)
            f.create_dataset('LayerOzone',data=scaled_layer_ozone[:,line-line_1,row-row_2],dtype=np.float)
            f.create_dataset('Latitude',data=omi_lat[line,row],dtype=np.float32)
            f.create_dataset('Longitude',data=omi_lon[line,row],dtype=np.float32)    
            f.create_dataset('Lw0',data=vl_data['Lw_Bas_0'],dtype=np.float32)
            f.create_dataset('Lw',data=vl_data['Lw_Bas'],dtype=np.float32)
            f.create_dataset('PS',data=m2_data.__dict__['PS'][line-line_1,row-row_2],dtype=np.float32)
            f.create_dataset('SZA',data=omi_sza[line,row],dtype=np.float32)
            f.create_dataset('VZA',data=omi_vza[line,row],dtype=np.float32)
            f.create_dataset('RAA',data=omi_raa[line,row],dtype=np.float32)
            f.create_dataset('TOARadiance_Lw',data=vl_data['Rad'],dtype=np.float32)
            f.create_dataset('Kiso',data=vl_data['Kiso'],dtype=np.float32)
            f.create_dataset('Kdir',data=vl_data['Kdir'],dtype=np.float32)
            f.create_dataset('Trans_Atmos',data=vl_data['TAtm'],dtype=np.float32)
            f.create_dataset('TOARadiance_Lw0',data=vl_data['Rad_0'],dtype=np.float32)
            f.create_dataset('WSP',data=omi_wsp[line,row],dtype=np.float32)
            f.create_dataset('Wavelength',data=vl_data['Wave'],dtype=np.float32)
            f.close()


