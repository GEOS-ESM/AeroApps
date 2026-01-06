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
   vars:
    - DUEM001
    - DUEM002
    - DUEM003
    - DUEM004
    - DUEM005
   plotvars:
    - [0,20,0]
    - Total Annual Emissions
  sedimentation:
   vars:
    - DUSD001
    - DUSD002
    - DUSD003
    - DUSD004
    - DUSD005
   plotvars:
    - [0,5,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - DUDP001
    - DUDP002
    - DUDP003
    - DUDP004
    - DUDP005
   plotvars:
    - [0,2,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - DUWT001
    - DUWT002
    - DUWT003
    - DUWT004
    - DUWT005
   plotvars:
    - [0,2,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - DUSV
   plotvars:
    - [0,2,0]
    - Total Annual Scavenging
  load:
   vars:
    - DUCMASS
   plotvars:
    - [0,2,0]
    - Total Annual Average Burden
  load25:
   vars:
    - DUCMASS25
   plotvars:
    - [0,0.5,0]
    - Total Annual Average PM2.5 Burden
  tau:
   vars:
    - DUEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average AOD
  ssa:
   vars:
    - DUSCATAU550
    - DUEXTTAU550
   plotvars:
    - [.91,.93,0]
    - Annual Average SSA
  tauabs:
   vars:
    - DUSCATAU550
    - DUEXTTAU550
   plotvars:
    - [0,.1,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - DUEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average Fraction of AOD
SS:
  emissions:
   vars:
    - SSEM001
    - SSEM002
    - SSEM003
    - SSEM004
    - SSEM005
   plotvars:
    - [0,5,0]
    - Total Annual Emissions
  sedimentation:
   vars:
    - SSSD001
    - SSSD002
    - SSSD003
    - SSSD004
    - SSSD005
   plotvars:
    - [0,2,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - SSDP001
    - SSDP002
    - SSDP003
    - SSDP004
    - SSDP005
   plotvars:
    - [0,.5,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - SSWT001
    - SSWT002
    - SSWT003
    - SSWT004
    - SSWT005
   plotvars:
    - [0,1,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - SSSV
   plotvars:
    - [0,1,0]
    - Total Annual Scavenging
  load:
   vars:
    - SSCMASS
   plotvars:
    - [0,.1,0]
    - Total Annual Average Burden
  load25:
   vars:
    - SSCMASS25
   plotvars:
    - [0,0.01,0]
    - Total Annual Average PM2.5 Burden
  tau:
   vars:
    - SSEXTTAU550
   plotvars:
    - [0,.2,0]
    - Annual Average AOD
  ssa:
   vars:
    - SSSCATAU550
    - SSEXTTAU550
   plotvars:
    - [.99,1,0]
    - Annual Average SSA
  tauabs:
   vars:
    - SSSCATAU550
    - SSEXTTAU550
   plotvars:
    - [0,.01,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - SSEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average Fraction of AOD
BC:
  emissions:
   vars:
    - BCEM001
    - BCEM002
   plotvars:
    - [0,0.04,0]
    - Total Annual Emissions
  emissions_bioburn:
   vars:
    - BCEMBB
   plotvars:
    - [0,0.04,0]
    - Total Annual Biomass Burning Emissions
  emissions_anthro:
   vars:
    - BCEMAN
   plotvars:
    - [0,0.04,0]
    - Total Annual Anthropogenic Emissions
  emissions_biofuel:
   vars:
    - BCEMBF
   plotvars:
    - [0,0.04,0]
    - Total Annual Biofuel Emissions
  prod_hyphil:
   vars:
    - BCHYPHIL
   plotvars:
    - [0,0.05,0]
    - Total Annual Hydrophilic Production
  sedimentation:
   vars:
    - BCSD001
    - BCSD002
   plotvars:
    - [0,0.01,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - BCDP001
    - BCDP002
   plotvars:
    - [0,0.01,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - BCWT002
   plotvars:
    - [0,0.04,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - BCSV
   plotvars:
    - [0,0.04,0]
    - Total Annual Scavenging
  load:
   vars:
    - BCCMASS
   plotvars:
    - [0,0.01,0]
    - Total Annual Average Burden
  tau:
   vars:
    - BCEXTTAU550
   plotvars:
    - [0,0.1,0]
    - Annual Average AOD
  ssa:
   vars:
    - BCSCATAU550
    - BCEXTTAU550
   plotvars:
    - [0.2,0.4,0]
    - Annual Average SSA
  tauabs:
   vars:
    - BCSCATAU550
    - BCEXTTAU550
   plotvars:
    - [0,.1,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - BCEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,0.2,0]
    - Annual Average Fraction of AOD
OC:
  emissions:
   vars:
    - OCEM001
    - OCEM002
   plotvars:
    - [0,0.4,0]
    - Total Annual Emissions
  emissions_bioburn:
   vars:
    - OCEMBB
   plotvars:
    - [0,0.4,0]
    - Total Annual Biomass Burning Emissions
  emissions_anthro:
   vars:
    - OCEMAN
   plotvars:
    - [0,0.2,0]
    - Total Annual Anthropogenic Emissions
  emissions_biofuel:
   vars:
    - OCEMBF
   plotvars:
    - [0,0.4,0]
    - Total Annual Biofuel Emissions
  emissions_biogenic:
   vars:
    - OCEMBG
   plotvars:
    - [0,0.2,0]
    - Total Annual Biogenic Emissions
  prod_hyphil:
   vars:
    - OCHYPHIL
   plotvars:
    - [0,0.1,0]
    - Total Annual Hydrophilic Production
  prod_SOA:
   vars:
    - OCPSOA
   plotvars:
    - [0,0.5,0]
    - Total Annual Secondary Organic Production
  sedimentation:
   vars:
    - OCSD001
    - OCSD002
   plotvars:
    - [0,0.01,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - OCDP001
    - OCDP002
   plotvars:
    - [0,0.1,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - OCWT002
   plotvars:
    - [0,0.2,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - OCSV
   plotvars:
    - [0,0.2,0]
    - Total Annual Scavenging
  load:
   vars:
    - OCCMASS
   plotvars:
    - [0,0.05,0]
    - Total Annual Average Burden
  tau:
   vars:
    - OCEXTTAU550
   plotvars:
    - [0,0.5,0]
    - Annual Average AOD
  ssa:
   vars:
    - OCSCATAU550
    - OCEXTTAU550
   plotvars:
    - [0.95,1,0]
    - Annual Average SSA
  tauabs:
   vars:
    - OCSCATAU550
    - OCEXTTAU550
   plotvars:
    - [0,.05,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - OCEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average Fraction of AOD
SU:
  emissions:
   vars:
    - SUEM003
   plotvars:
    - [0,0.1,0]
    - Total Annual SO4 Emissions
  emissions_so2:
   vars:
    - SUEM002
   plotvars:
    - [0,1,0]
    - Total Annual SO2 Emissions
  emissions_dms:
   vars:
    - SUEM001
   plotvars:
    - [0,0.05,0]
    - Total Annual DMS Emissions
  prod_so4:
   vars:
    - SUPSO4G
    - SUPSO4AQ
    - SUPSO4WT
   plotvars:
    - [0,0.5,0]
    - Total Annual Chemical Production of SO4
  prod_so2:
   vars:
    - SUPSO2
   plotvars:
    - [0,0.05,0]
    - Total Annual Chemical Production of SO2
  sedimentation:
   vars:
    - SUSD003
   plotvars:
    - [0,0.01,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - SUDP003
   plotvars:
    - [0,0.1,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - SUWT003
   plotvars:
    - [0,0.5,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - SUSV
   plotvars:
    - [0,0.2,0]
    - Total Annual Scavenging
  load:
   vars:
    - SO4CMASS
   plotvars:
    - [0,0.02,0]
    - Total Annual Average Burden
  tau:
   vars:
    - SUEXTTAU550
   plotvars:
    - [0,0.5,0]
    - Annual Average AOD
  ssa:
   vars:
    - SUSCATAU550
    - SUEXTTAU550
   plotvars:
    - [0.95,1,0]
    - Annual Average SSA
  tauabs:
   vars:
    - SUSCATAU550
    - SUEXTTAU550
   plotvars:
    - [0,.05,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - SUEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average Fraction of AOD
BR:
  emissions:
   vars:
    - BREM001
    - BREM002
   plotvars:
    - [0,0.5,0]
    - Total Annual Emissions
  prod_hyphil:
   vars:
    - BRHYPHIL
   plotvars:
    - [0,0.5,0]
    - Total Annual Hydrophilic Production
  prod_SOA:
   vars:
    - BRPSOA
   plotvars:
    - [0,0.1,0]
    - Total Annual Secondary Organic Production
  sedimentation:
   vars:
    - BRSD001
    - BRSD002
   plotvars:
    - [0,0.01,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - BRDP001
    - BRDP002
   plotvars:
    - [0,0.2,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - BRWT002
   plotvars:
    - [0,0.2,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - BRSV
   plotvars:
    - [0,0.2,0]
    - Total Annual Scavenging
  load:
   vars:
    - BRCMASS
   plotvars:
    - [0,0.2,0]
    - Total Annual Average Burden
  tau:
   vars:
    - BREXTTAU550
   plotvars:
    - [0,0.75,0]
    - Annual Average AOD
  ssa:
   vars:
    - BRSCATAU550
    - BREXTTAU550
   plotvars:
    - [0.95,1,0]
    - Annual Average SSA
  tauabs:
   vars:
    - BRSCATAU550
    - BREXTTAU550
   plotvars:
    - [0,.02,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - BREXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,1,0]
    - Annual Average Fraction of AOD
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
NI:
  emissions:
   vars:
    - NH3EM
   plotvars:
    - [0,1,0]
    - Total Annual Ammonia Emissions
  prod_aq:
   vars:
    - NIPNO3AQ
   plotvars:
    - [-0.05,0.05,0]
    - Total Annual Aqueous Production of Nitrate
  prod_het:
   vars:
    - NIHT001
    - NIHT002
    - NIHT003
   plotvars:
    - [0,0.5,0]
    - Total Annual Heterogeneous Production of Nitrate
  sedimentation:
   vars:
    - NISD001
    - NISD002
    - NISD003
   plotvars:
    - [0,0.1,0]
    - Total Annual Sedimentation
  dry_deposition:
   vars:
    - NIDP001
    - NIDP002
    - NIDP003
   plotvars:
    - [0,0.1,0]
    - Total Annual Dry Deposition
  wet_deposition:
   vars:
    - NIWT001
    - NIWT002
    - NIWT003
   plotvars:
    - [0,0.2,0]
    - Total Annual Wet Deposition
  scavenging:
   vars:
    - NISV
   plotvars:
    - [0,0.2,0]
    - Total Annual Scavenging
  load:
   vars:
    - NICMASS
   plotvars:
    - [0,0.05,0]
    - Total Annual Average Burden
  load25:
   vars:
    - NICMASS25
   plotvars:
    - [0,0.05,0]
    - Total Annual Average Burden
  tau:
   vars:
    - NIEXTTAU550
   plotvars:
    - [0,0.5,0]
    - Annual Average AOD
  ssa:
   vars:
    - NISCATAU550
    - NIEXTTAU550
   plotvars:
    - [0.95,1,0]
    - Annual Average SSA
  tauabs:
   vars:
    - NISCATAU550
    - NIEXTTAU550
   plotvars:
    - [0,.05,0]
    - Annual Average AAOD
  taufrac:
   vars:
    - NIEXTTAU550
    - TOTEXTTAU550
   plotvars:
    - [0,.5,0]
    - Annual Average Fraction of AOD
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

def printString(f,field,val):
    txt = '{field}{val:.4f}'
    print(txt.format(field=field,val=val), file=f)
    return


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

#       If SSA assume first is scattering and second is extinction
#       Output_t is area and extinction weighted SSA
        if Field == 'ssa':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]['vars']
                levs   = self.table[s][Field]['plotvars'][0]
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

#       For AOD use simple area weighted averaging in totals
        elif Field == 'tau':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]['vars']
                levs   = self.table[s][Field]['plotvars'][0]
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
                    field_name = self.table[s][Field]['plotvars'][1]
                    txt = '{field} ({val:.3f} Tg)'
                    str_out.append(txt.format(field=field_name,val=output_t[12]))

#       For Load use integrated area weighted averaging in totals
        elif Field == 'load' or Field == 'load25':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]['vars']
                levs   = self.table[s][Field]['plotvars'][0]
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
                    field_name = self.table[s][Field]['plotvars'][1]
                    txt = '{field} ({val:.4f} Tg)'
                    str_out.append(txt.format(field=field_name,val=output_t[12]))
                
#       If Tauabs assume first is scattering and second is extinction
#       or if taufrac first is species and second is total
#       and use simple area weighted average in totals
        elif Field == 'tauabs' or Field == 'taufrac':
            for s in Species:   # species
                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]['vars']
                levs   = self.table[s][Field]['plotvars'][0]
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
                    field_name = self.table[s][Field]['plotvars'][1]
                    txt = '{field} ({val:.4f})'
                    str_out.append(txt.format(field=field_name,val=output_t[12]))

#       Else treat like a flux [kg m-2 mon-1] (total is area integrated sum)
        else:
            for s in Species:   # species

                if self.verbose:
                    print('[] working on',s)

                Fields = self.table[s][Field]['vars']
                levs   = self.table[s][Field]['plotvars'][0]
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
                    field_name = self.table[s][Field]['plotvars'][1]
                    txt = '{field} ({val:.2f} Tg)'
                    str_out.append(txt.format(field=field_name,val=output_t[12]))
                    units = 'g m-2 mon-1'
                    
        if not allok:
            print('%s not found; setting to zero'%(fldnot))
            str_out = ["","","","","","","","","","","","",""]

#       Append the last entry of output_t to the list for the Field
        self.table[s][Field]['vars'].append(output_t[12])
                    
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
        return

    def plotAll(self,expid,verbose=False):
        Species = list(self.table.keys())
        for s in Species:
            Fields = list(self.table[s].keys())
            for fld in Fields:
                self.plotField(expid,Field=fld,Species=s,verbose=True)
        return

    def printBudget(self,expid,verbose=False):
        f = open('%s/%s.budget.txt'%(expid,expid),'w')
        Species = list(self.table.keys())
        for s in Species:
            print(expid, file=f)
            print(s, file=f)
            printString(f,'Emissions [Tg]:         ', self.table[s]['emissions']['vars'][-1])
            if s == 'SU':
                printString(f,'Emissions (SO2) [Tg]:   ', self.table[s]['emissions_so2']['vars'][-1])
                printString(f,'Emissions (DMS) [Tg]:   ', self.table[s]['emissions_dms']['vars'][-1])
                printString(f,'Production (SO4) [Tg]:  ', self.table[s]['prod_so4']['vars'][-1])
                printString(f,'Production (SO2) [Tg]:  ', self.table[s]['prod_so2']['vars'][-1])
            if s == 'NI':
                printString(f,'Production (AQ) [Tg]:   ', self.table[s]['prod_aq']['vars'][-1])
                printString(f,'Production (HET) [Tg]:  ', self.table[s]['prod_het']['vars'][-1])
            wet  = self.table[s]['wet_deposition']['vars'][-1]
            scav = self.table[s]['scavenging']['vars'][-1]
            dep  = self.table[s]['dry_deposition']['vars'][-1]
            sed  = self.table[s]['sedimentation']['vars'][-1]
            loss = dep + wet + sed + scav
            printString(f,'Losses [Tg]:            ', loss)
            printString(f,'-Dry [Tg]:              ',dep+sed)
            printString(f,'-Wet [Tg]:              ',wet+scav)
            printString(f,'Burden [Tg]:            ',self.table[s]['load']['vars'][-1])
            life = self.table[s]['load']['vars'][-1] / loss * 365.
            printString(f,'Lifetime [days]:        ',life)
            printString(f,'-Wet Removal [1/days]:  ', (wet+scav)/365. /self.table[s]['load']['vars'][-1])
            printString(f,' -Large Scale [1/days]: ', wet/365. /self.table[s]['load']['vars'][-1])
            printString(f,' -Scavenging [1/days]:  ', scav/365. /self.table[s]['load']['vars'][-1])
            printString(f,'-Dry Removal [1/days]:  ', (dep+sed)/365. /self.table[s]['load']['vars'][-1])
            printString(f,' -Settling [1/days]:    ', sed/365. /self.table[s]['load']['vars'][-1])
            printString(f,' -Deposition [1/days]:  ', dep/365. /self.table[s]['load']['vars'][-1])
            printString(f,'AOT:                    ',self.table[s]['tau']['vars'][-1])
        f.close()
        return

if __name__ == "__main__":
    expid = 'c180R_v11.8.0_newbrcoptics'
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
