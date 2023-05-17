
#!/usr/bin/env python3
"""
	Puts together the GOES-R tiles into 1 file
"""

import os
import glob
from netCDF4 import Dataset
from datetime import datetime, timedelta 
from dateutil.parser import parse
import numpy as np
#---
if __name__ == '__main__':

	# USER INPUTS
	# -------------------------------
	instname   = 'GOES-R' 
	runname    = 'vlidort.vector.iMAIACRTLS.'
	g5name     = 'goes-r-g5nr.lc2.'
	indir      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/GOES-R/DATA/'
	channels   = '2250.00'
	startdate  = '2005-12-31T13:00:00Z'
	enddate    = '2005-12-31T13:00:00Z'
	episode    = 4
	layout     = '41'


	# End User INPUTS
	# -------------------------------
	if episode is not None:
		if episode == 1:
			startdate  = '2005-12-31T00:00:00Z'
			enddate    = '2006-01-01T23:00:00Z'
		elif episode == 2:
			startdate  = '2006-07-27T00:00:00Z'
			enddate    = '2006-08-09T23:00:00Z'
		elif episode == 3:
			startdate  = '2007-03-28T00:00:00Z'
			enddate    = '2007-03-29T23:00:00Z'
		elif episode == 4:
			startdate  = '2007-04-10T00:00:00Z'
			enddate    = '2007-04-11T23:00:00Z'    

	dt          = timedelta(hours=1)
	date        = parse(startdate)
	enddate     = parse(enddate)    
	instnamel   = instname.lower()
	instnameu   = instname.upper()  
	domainvars  = ['time','lev','ew','ns','scanTime','clon','clat']	
	fileatts    = ['angle_inputs','land_inputs','met_inputs','aerosol_inputs']
	attlist     = ['title','institution','source','history', 'grid_inputs'] + \
								fileatts + \
								['surface_inputs','inputs','references','comment','surface_comment','vlidort_comment',
								 'scat_comment','mie_comment','channels','contact','Conventions','Version','Created']
	varatts     = ['long_name','units','positive','axis','missing_value','standard_name','_FillValues']

	if type(channels) is str:
		channels = channels.split()

	# Grid Information
	# -------------------------------
	tyme   = indir + 'LevelG/invariant/'+instnamel+'.lg1.scantime.nc4'
	grid   = indir + 'LevelG/invariant/'+instnamel+'.lg1.invariant.nc4'
	ncGrid = Dataset(grid)
	ncTyme = Dataset(tyme)
	clon = ncGrid.variables['clon'][:,:]
	nX = int(layout[0])
	nY = int(layout[1])
	nNS, nEW = clon.shape		


	while (date <= enddate):
		# Datetime variables
		# -------------------------------
		nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
		hour  = str(date.hour).zfill(2) + 'z'
		nhms   = date.hour * 10000

		for ch in channels:

			varlist = domainvars + ['ref_'+ch,'surf_ref_'+ch,'aod_'+ch]

			ch = ch.replace('.','d')

			indirlong = indir + 'LevelC2/Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2)
			data   = indirlong  + '/'+ g5name + runname + nymd + '_' + hour + '_' + ch + 'nm.'

			outdir = indirlong + '/stitched/'
			
			if (not os.path.exists(outdir)):
				os.makedirs(outdir)

			outfile = outdir + g5name + runname + nymd + '_' + hour + '_' + ch + 'nm.nc4'

			print('data: ', data)
			print('outfile: ', outfile)

			try:
				template = glob.glob(data+'*')[0]
			except:
				template = ''

			if (template):
				# Create outfile
				# -------------------------
				ncOut = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
				
				ncTemplate  = Dataset(template)

				filedic             = {}.fromkeys(varlist)
				filedic['time']     = ncTemplate
				filedic['lev']      = ncTemplate
				filedic['ew']       = ncGrid
				filedic['ns']       = ncGrid
				filedic['scanTime'] = ncTyme
				filedic['clon']     = ncGrid
				filedic['clat']     = ncGrid

				# -- Global attributes
				for att in attlist:
					if att == 'channels':
						setattr(ncOut,att,ch.replace('d','.'))
					elif any(s == att for s in fileatts):
						tempname = getattr(ncTemplate,att)
						tempname = tempname[:-7] + '*.nc4'
						setattr(ncOut,att,tempname)
					else:
						setattr(ncOut,att,getattr(ncTemplate,att))

				# -- Define dimensions
				time  = ncOut.createDimension('time', len(ncTemplate.dimensions['time']))
				lev   = ncOut.createDimension('lev', len(ncTemplate.dimensions['lev']))		
				ew    = ncOut.createDimension('ew', len(ncGrid.dimensions['ew']))
				ns    = ncOut.createDimension('ns', len(ncGrid.dimensions['ns']))		

				# # Read in Data and put fields back together
				# # -------------------------------------------
				def copy_vatts(var,varatts,ncFile,varobj):
					for vatt in varatts:
						if any(vatt == s for s in dir(ncFile.variables[var])):
							setattr(varobj,vatt,getattr(ncFile.variables[var],vatt))


				for var in varlist:

					if any(var == s for s in domainvars):						

						if var == 'time':
							varobj = ncOut.createVariable(var,'i4',('time',))
							copy_vatts(var,varatts,filedic[var],varobj)			
							varobj[:]  = filedic[var].variables[var][:]				

						elif any(var == ss for ss in ['ew','ns','lev']):
							varobj = ncOut.createVariable(var,'f4',(var,))
							copy_vatts(var,varatts,filedic[var],varobj)
							varobj[:] = filedic[var].variables[var][:]

						elif var == 'scanTime':
							varobj = ncOut.createVariable(var,'f4',('ew',))
							copy_vatts(var,varatts,ncTemplate,varobj)
							varobj[:] = filedic[var].variables[var][0,:] 

						else:
							varobj = ncOut.createVariable(var,'f4',('ns','ew',))
							copy_vatts(var,varatts,filedic[var],varobj)
							varobj[:] = filedic[var].variables[var][:]
					else:

						Data = np.ma.masked_all([nNS,nEW])
						for i in np.arange(nX*nY): 
							i_ = i%nX
							j_ = int(i/nX) 
							
							Xstart = i_*nEW/nX
							Xend   = Xstart + nEW/nX
							Ystart = j_*nNS/nY
							Yend   = Ystart + nNS/nY

							datalayout = data + layout+str(i) + '.nc4'
							if os.path.exists(datalayout):
								ncData   = Dataset(datalayout)
								temp     = np.squeeze(ncData.variables[var][:])
								Data[Ystart:Yend,Xstart:Xend] = temp
								ncData.close()

						varobj = ncOut.createVariable(var,'f4',('time','lev','ns','ew',))
						copy_vatts(var,varatts,ncTemplate,varobj)
						varobj[0,0,:,:] = Data

				ncOut.close()
				ncTemplate.close()

		date = date + dt

	ncGrid.close()
