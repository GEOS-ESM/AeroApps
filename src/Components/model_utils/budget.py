#!/usr/bin/env python3

"""
    Peter Colarco, July 2025 (Peter.R.Colarco@nasa.gov)
    Utility to compute aerosol budgets from GEOS GOCART output
    Variable mapping is handled in budget.yaml file
    Assumes a full year of model output available in specified
    GrADS control file
"""

import os
import yaml
import xarray as xr
from pyobs import xrctl     as xc
import numpy as np
import math
import sys

import matplotlib.pyplot as plt
from matplotlib import cm #cm is colormap
import matplotlib.ticker as mticker
import matplotlib
matplotlib.use('agg')

import cartopy
from cartopy import util
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Default YAML file mapping GOCART species to budget terms saved in ar_Nx files
# -----------------------------------------------------------------------------
G2G_BudgetMap = """
DU:
  emissions:
    - DUEM001
    - DUEM002
    - DUEM003
    - DUEM004
    - DUEM005
  sedimentation:
    - DUSD001
    - DUSD002
    - DUSD003
    - DUSD004
    - DUSD005
  dry_deposition:
    - DUDP001
    - DUDP002
    - DUDP003
    - DUDP004
    - DUDP005
  wet_deposition:
    - DUWT001
    - DUWT002
    - DUWT003
    - DUWT004
    - DUWT005
  scavenging:
    - DUSV
  load:
    - DUCMASS
  load25:
    - DUCMASS25
  tau:
    - DUEXTTAU550
  ssa:
    - DUSCATAU550
    - DUEXTTAU550
  tauabs:
    - DUSCATAU550
    - DUEXTTAU550
  taufrac:
    - DUEXTTAU550
    - TOTEXTTAU550
SS:
  emissions:
    - SSEM001
    - SSEM002
    - SSEM003
    - SSEM004
    - SSEM005
  sedimentation:
    - SSSD001
    - SSSD002
    - SSSD003
    - SSSD004
    - SSSD005
  dry_deposition:
    - SSDP001
    - SSDP002
    - SSDP003
    - SSDP004
    - SSDP005
  wet_deposition:
    - SSWT001
    - SSWT002
    - SSWT003
    - SSWT004
    - SSWT005
  scavenging:
    - SSSV
  load:
    - SSCMASS
  load25:
    - SSCMASS25
  tau:
    - SSEXTTAU550
  ssa:
    - SSSCATAU550
    - SSEXTTAU550
  tauabs:
    - SSSCATAU550
    - SSEXTTAU550
  taufrac:
    - SSEXTTAU550
    - TOTEXTTAU550
BC:
  emissions:
    - BCEM001
    - BCEM002
  emissions_bioburn:
    - BCEMBB
  emissions_anthro:
    - BCEMAN
  emissions_biofuel:
    - BCEMBF
  prod_hyphil:
    - BCHYPHIL
  sedimentation:
    - BCSD001
    - BCSD002
  dry_deposition:
    - BCDP001
    - BCDP002
  wet_deposition:
    - BCWT002
  scavenging:
    - BCSV
  load:
    - BCCMASS
  tau:
    - BCEXTTAU550
  ssa:
    - BCSCATAU550
    - BCEXTTAU550
  tauabs:
    - BCSCATAU550
    - BCEXTTAU550
  taufrac:
    - BCEXTTAU550
    - TOTEXTTAU550
OC:
  emissions:
    - OCEM001
    - OCEM002
  emissions_bioburn:
    - OCEMBB
  emissions_anthro:
    - OCEMAN
  emissions_biofuel:
    - OCEMBF
  emissions_biogenic:
    - OCEMBG
  prod_hyphil:
    - OCHYPHIL
  prod_SOA:
    - OCPSOA
  sedimentation:
    - OCSD001
    - OCSD002
  dry_deposition:
    - OCDP001
    - OCDP002
  wet_deposition:
    - OCWT002
  scavenging:
    - OCSV
  load:
    - OCCMASS
  tau:
    - OCEXTTAU550
  ssa:
    - OCSCATAU550
    - OCEXTTAU550
  tauabs:
    - OCSCATAU550
    - OCEXTTAU550
  taufrac:
    - OCEXTTAU550
    - TOTEXTTAU550
SU:
  emissions:
    - SUEM003
  emissions_so2:
    - SUEM002
  emissions_dms:
    - SUEM
  prod_so4:
    - SUPSO4G
    - SUPSO4AQ
    - SUPSO4WT
  prod_so2:
    - SUPSO2
  sedimentation:
    - SUSD003
  dry_deposition:
    - SUDP003
  wet_deposition:
    - SUWT003
  scavenging:
    - SUSV
  load:
    - SO4CMASS
  tau:
    - SUEXTTAU550
  ssa:
    - SUSCATAU550
    - SUEXTTAU550
  tauabs:
    - SUSCATAU550
    - SUEXTTAU550
  taufrac:
    - SUEXTTAU550
    - TOTEXTTAU550
BR:
  emissions:
    - BREM001
    - BREM002
  prod_hyphil:
    - BRHYPHIL
  prod_SOA:
    - BRPSOA
  sedimentation:
    - BRSD001
    - BRSD002
  dry_deposition:
    - BRDP001
    - BRDP002
  wet_deposition:
    - BRWT002
  scavenging:
    - BRSV
  load:
    - BRCMASS
  tau:
    - BREXTTAU550
  ssa:
    - BRSCATAU550
    - BREXTTAU550
  tauabs:
    - BRSCATAU550
    - BREXTTAU550
  taufrac:
    - BREXTTAU550
    - TOTEXTTAU550
NI:
  emissions:
    - NH3EM
  prod_aq:
    - NIPNO3AQ
  prod_het:
    - NIHT
    - NIHT002
    - NIHT003
  sedimentation:
    - NISD
    - NISD002
    - NISD003
  dry_deposition:
    - NIDP
    - NIDP002
    - NIDP003
  wet_deposition:
    - NIWT
    - NIWT002
    - NIWT003
  scavenging:
    - NISV
  load:
    - NICMASS
  load25:
    - NICMASS25
  tau:
    - NIEXTTAU550
  ssa:
    - NISCATAU550
    - NIEXTTAU550
  tauabs:
    - NISCATAU550
    - NIEXTTAU550
  taufrac:
    - NIEXTTAU550
    - TOTEXTTAU550
"""


# Compute the area of grid boxes given input lat and lon
# arrays, assuming GEOS regularization
def Area(lon,lat):
    nx = len(lon)
    ny = len(lat)
    dx = 360./nx
    dy = 180./(ny-1)

    rearth = 6370000.
    pi = math.pi

    area = np.zeros((nx,ny))
    #End points
    area[:,ny-1] = 2.*pi*rearth**2.*(1.-np.sin((lat[ny-1]-dy/2.)*pi/180.))*dx/360.
    area[:,0] = area[:,ny-1]
    for j in range(1,ny-1):
        area[:,j] = 2.*pi*rearth**2.*dx/360. *(np.sin((lat[j]+dy/2.)*pi/180.)-np.sin((lat[j]-dy/2.)*pi/180.))

    return area, nx, ny, dx, dy

class AEROSOL(object):

    def __init__ (self,aerFiles,config=None,verbose=False):
        """
        Lazy loads GOCART aerosol budget files (*=aer_Nx) files.

        aerFiles:  str, list, or Dataset with aerosol tracers
        config:    str or YAML file handle with mapping to aerosol budget terms.
                   If None, uses internal default.
        """

#       Load configuration
        if config is None:
            config = G2G_BudgetMap
        elif type(config) is str:
            # get a file handle
            config = open(config)

        self.table = yaml.safe_load(config)

        self.verbose = verbose
        if isinstance(aerFiles,xr.Dataset):
            self.aer = aerFiles
        else:
            self.aer = xc.open_mfdataset(aerFiles)

        try:
            self.lon = self.aer['longitude'].load()
        except:
            self.lon = self.aer['lon'].load()

        try:
            self.lat = self.aer['latitude'].load()
        except:
            self.lat = self.aer['lat'].load()
        lon = self.lon.values
        lat = self.lat.values
        area, nx, ny, dx, dy = Area(lon,lat)
        self.area = np.transpose(area)


    def getField(self,Field=None,Species=None,verbose=False):
        """
        Returns an xarray Dataset with the requested field integrated across
        sub-fields (e.g., return emissions integrated DUEM001+DUEM002+...)
        Also return monthly values
        Fluxes are returned as kg m-2 month-1

        Field:    string
        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of emissions.
        """ 

        months = ['January','February','March','April','May','June',
                  'July','August','September','October','November','December']

        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.table.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.aer   # budget file fields

#       Get a list of variables in the dataset
        vars = list(a.keys())
        
#       Get the dimensions
        lwi = a['LWI'].copy()
        space = lwi.shape
        output = (np.zeros(space))
        output_t = np.zeros(13)  # last entry is total or average

        str_out = []

#       Need mod for leap year...
        nday = [31,28,31,30,31,30,31,31,30,31,30,31]

        simpleaverage = ['tau','load','load25']

        units   = ''
        units_t = ''

        allok = True

#       Levels for plotting in units assigned
        for s in Species: # defaults per species
            if s == 'DU':
                levs = [.1,10,1]
            elif s == 'BC':
                levs = [.001,.1,1]
            else:
                levs = [.1,10,1]

#       Override for specific fields
        if Field == 'tau':
            levs=[0,1,0]
        elif Field == 'taufrac':
            levs = [0,1.,0]
        elif Field == 'tauabs':
            levs = [0,.1,0]


        
#       If SSA assume first is scattering and second is extinction
#       Output_t is area and extinction weighted SSA
        if Field == 'ssa':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]
                if self.verbose:
                    print('   -',Field)

                for fld in Fields:
                    if fld not in vars:
                        allok = False
                        fldnot = fld
                if allok:
                    sca = a[Fields[0]].compute()
                    tau = a[Fields[1]].compute()
                    output_ = sca/tau
                    output  = output_.values
                    for i in range(12):
                        ssa  = np.squeeze(output[i,:,:])
                        tau_ = np.squeeze(tau[i,:,:])
                        output_t[i] = np.sum(tau_*ssa*self.area)/np.sum(tau_*self.area)
                        txt = '{month} ({val:.3f})'
                        str_out.append(txt.format(month=months[i],val=output_t[i]))
                    output_t[12] = np.mean(output_t[0:12])
                    txt = 'Annual Average SSA ({val:.3f})'
                    str_out.append(txt.format(val=output_t[12]))
                    if s == 'BC':
                        levs = [0.2,0.4,0]
                    else:
                        levs = [0.9,1,0]

#       For AOD use simple area weighted averaging in totals
        elif Field == 'tau':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]
                if self.verbose:
                    print('   -',Field)
                for fld in Fields:
                    if fld not in vars:
                        allok = False
                        fldnot = fld
                if allok:
                    output_ = a[Fields[0]].compute()
                    output  = output_.values
                    for i in range(12):
                        val  = np.squeeze(output[i,:,:])
                        output_t[i] = np.sum(val*self.area)/np.sum(self.area)
                        txt = '{month} ({val:.3f})'
                        str_out.append(txt.format(month=months[i],val=output_t[i]))
                    output_t[12] = np.mean(output_t[0:12])
                    txt = 'Annual Average AOD ({val:.3f})'
                    str_out.append(txt.format(val=output_t[12]))

#       For Load use integrated area weighted averaging in totals
        elif Field == 'load' or Field == 'load25':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]
                if self.verbose:
                    print('   -',Field)
                for fld in Fields:
                    if fld not in vars:
                        allok = False
                        fldnot = fld
                if allok:
                    output_ = a[Fields[0]].compute()*1000. # g m-2
                    output  = output_.values
                    for i in range(12):
                        val  = np.squeeze(output[i,:,:])
                        output_t[i] = np.sum(val*self.area)/1.e12
                        txt = '{month} ({val:4.2f} Tg)'
                        str_out.append(txt.format(month=months[i],val=output_t[i]))
                    output_t[12] = np.mean(output_t[0:12])
                    if Field == 'load':
                        txt = 'Annual Average Total Burden ({val:4.2f} Tg)'
                        if (s == 'DU') | (s == 'SS'):
                            levs = [0.01,2,1]
                        else:
                            levs = [0.001,1,1]
                    else:
                        txt = 'Annual Average PM2.5 Burden ({val:4.2f} Tg)'
                        if (s == 'DU') | (s == 'SS'):
                            levs = [0.01,.5,1]
                        else:
                            levs = [0.001,.1,1]
                    str_out.append(txt.format(val=output_t[12]))
                
#       If Tauabs assume first is scattering and second is extinction
#       or if taufrac first is species and second is total
#       and use simple area weighted average in totals
        elif Field == 'tauabs' or Field == 'taufrac':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]
                if self.verbose:
                    print('   -',Field)
                for fld in Fields:
                    if fld not in vars:
                        allok = False
                        fldnot = fld
                if allok:
                    if Field == 'tauabs':
                        output_ = a[Fields[1]].compute()-a[Fields[0]].compute()
                    else:
                        output_ = a[Fields[0]].compute()/a[Fields[1]].compute()
                    output  = output_.values
                    for i in range(12):
                        val  = np.squeeze(output[i,:,:])
                        output_t[i] = np.sum(val*self.area)/np.sum(self.area)
                        txt = '{month} ({val:.4f})'
                        str_out.append(txt.format(month=months[i],val=output_t[i]))
                    output_t[12] = np.mean(output_t[0:12])
                    if Field == 'tauabs':
                        txt = 'Annual Average AAOD ({val:.4f})'
                    else:
                        txt = 'Annual Average Fraction of Total AOD ({val:.4f})'
                    str_out.append(txt.format(val=output_t[12]))

#       Else treat like a flux [kg m-2 mon-1] (total is area integrated sum)
        else:
            for s in Species:   # species

                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]
                for fld in Fields:
                    if fld not in vars:
                        allok = False
                        fldnot = fld
                if allok:
                    for q in Fields:
                        if self.verbose:
                            print('   -',q)
                            output_ = a[q].compute()
                            output += output_.values
#                   Fix sign if scavenging
                    if Field == 'scavenging':
                        output   = -1.*output
#                   Rescale to g m-2 mon-1 and make totals
                    for i in range(12):
                        output[i,:,:] = output[i,:,:]*1000.*nday[i]*86400.
                        val = np.squeeze(output[i,:,:])
                        output_t[i] = np.sum(val*self.area)/1.e12  # Tg
                        txt = '{month} ({val:7.3f} Tg)'
                        str_out.append(txt.format(month=months[i],val=output_t[i]))
                    output_t[12] = np.sum(output_t[0:12])
                    fldname = Field.split("_")
                    if len(fldname) == 2:
                        field_name = fldname[0].capitalize()+' '+fldname[1].capitalize()
                    else:
                        field_name = Field.capitalize()
                    txt = 'Annual Total {field} ({val:.2f} Tg)'
                    str_out.append(txt.format(field=field_name,val=output_t[12]))
                    units = 'g m-2 mon-1'
                    
        if not allok:
            print('%s not found; setting to zero'%(fldnot))
            str_out = ["","","","","","","","","","","","",""]
            
#       Append the last entry of output_t to the list for the Field
        self.table[s][Field].append(output_t[12])
        
        return output, output_t, levs, str_out, units

    def plotField(self,expid,Field=None,Species=None,verbose=False):
        """
        Given a field and species. read the field and make a simple plot
        Assumes 12 months
        """

        months = ['January','February','March','April','May','June',
                  'July','August','September','October','November','December']

        
        val, val_t, levs, str_out, units  = self.getField(Field,Species,verbose=True)

        lon = self.lon.values
        lat = self.lat.values

        map_projection = ccrs.PlateCarree()
        fig, axs = plt.subplots(3,4,figsize=(24,14),layout='constrained',
                                subplot_kw={'projection': ccrs.PlateCarree()})
        if levs[2] == 1:
            a = np.log10(levs[0])
            b = np.log10(levs[1])
            d = (b-a)/50.
            clevs = 10.**(np.arange(a,b,d))
        else:
            a = levs[0]
            b = levs[1]
            d = (b-a)/50.
            clevs = np.arange(a,b,d)
#        for i in range(3):
#            for j in range(4):
#                k = j*3+i
#                print(k)
        k = 0
        for ax in axs.flat:
            pout = np.squeeze(val[k,:,:])
            im = ax.contourf(lon,lat,pout,levels=clevs,cmap='YlOrBr',extend='max')
            ax.set_title(str_out[k],size=24)
            ax.add_feature(cfeature.BORDERS, linestyle='-')
            ax.add_feature(cfeature.COASTLINE, linestyle='-')
            k += 1

        cb = fig.colorbar(im, ax=axs, location='bottom',shrink=.8)
        cb.ax.tick_params(labelsize=24)
        fig.suptitle('%s %s '%(expid,Species)+str_out[12], size=30)
        plt.savefig('%s/%s.%s.%s.png'%(expid,expid,Species,Field))
        plt.close()

    def plotAll(self,expid,verbose=False):
        Species = list(self.table.keys())
        for s in Species:
            Fields = list(self.table[s].keys())
            for fld in Fields:
                self.plotField(expid,Field=fld,Species=s,verbose=True)

    def printBudget(self,expid,verbose=False):
        f = open('%s/%s.budget.txt'%(expid,expid),'w')
        Species = list(self.table.keys())
        for s in Species:
            print(expid, file=f)
            print(s, file=f)
            print('Emissions [Tg]:         ', self.table[s]['emissions'][-1], file=f)
            if s == 'SU':
                print('Emissions (SO2) [Tg]:    ', self.table[s]['emissions_so2'][-1], file=f)
                print('Emissions (DMS) [Tg]:    ', self.table[s]['emissions_dms'][-1], file=f)
                print('Production (SO4) [Tg]:   ', self.table[s]['prod_so4'][-1], file=f)
                print('Production (SO2) [Tg]:   ', self.table[s]['prod_so2'][-1], file=f)
            if s == 'NI':
                print('Production (AQ) [Tg]:    ', self.table[s]['prod_aq'][-1], file=f)
                print('Production (HET) [Tg]:   ', self.table[s]['prod_het'][-1], file=f)
            wet  = self.table[s]['wet_deposition'][-1]
            scav = self.table[s]['scavenging'][-1]
            dep  = self.table[s]['dry_deposition'][-1]
            sed  = self.table[s]['sedimentation'][-1]
            loss = dep + wet + sed + scav
            print('Losses [Tg]:            ', loss, file=f)
            print('-Dry [Tg]:              ',dep+sed, file=f)
            print('-Wet [Tg]:              ',wet+scav, file=f)
            print('Burden [Tg]:            ',self.table[s]['load'][-1], file=f)
            life = self.table[s]['load'][-1] / loss * 365.
            print('Lifetime [days]:        ',life, file=f)
            print('-Wet Removal [1/days]:  ', (wet+scav)/365. /self.table[s]['load'][-1], file=f)
            print(' -Large Scale [1/days]: ', wet/365. /self.table[s]['load'][-1], file=f)
            print(' -Scavenging [1/days]:  ', scav/365. /self.table[s]['load'][-1], file=f)
            print('-Dry Removal [1/days]:  ', (dep+sed)/365. /self.table[s]['load'][-1], file=f)
            print(' -Settling [1/days]:    ', sed/365. /self.table[s]['load'][-1], file=f)
            print(' -Deposition [1/days]:  ', dep/365. /self.table[s]['load'][-1], file=f)
            print('AOT:                    ',self.table[s]['tau'][-1], file=f)
        f.close()    

if __name__ == "__main__":
    expid = 'c180R_v11.8.0_develop'
    aer_Nx = '%s.tavg2d_aer_x.ctl'%(expid)
    b = AEROSOL(aer_Nx,verbose=True)
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
