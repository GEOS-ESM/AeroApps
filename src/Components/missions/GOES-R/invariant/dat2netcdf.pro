pro dat2netcdf
	; A program to read in the binary .dat files that describe the GOES-R orbit
	; Save the data to a lg1 netcdf file - trim the domain to correspond to TEMPO
	; 2015 Aug P. Castellanos

	; Read in binary GOES-R lat, lons, and area
	; -----------------------------------------
	dir='/nobackup/GOES-R/LevelG/invariant/Attic/'
	openr,1,dir+'nxny_goes-r'
	readf,1,npixels,nlines
	close,1
	npixels=fix(npixels)
	nlines=fix(nlines)
	print, 'original dimensions', npixels,nlines
	goes_clat=fltarr(npixels,nlines)
	goes_clon=fltarr(npixels,nlines)
	goes_elat=fltarr(npixels+1,nlines+1)
	goes_elon=fltarr(npixels+1,nlines+1)
	area=fltarr(npixels,nlines)

	goes_clat=read_binary(dir+'center_lat_west.dat',data_type=4,data_dims=[npixels,nlines])
	goes_clon=read_binary(dir+'center_lon_west.dat',data_type=4,data_dims=[npixels,nlines])
	goes_elat=read_binary(dir+'edge_lat_west.dat',data_type=4,data_dims=[npixels+1,nlines+1])
	goes_elon=read_binary(dir+'edge_lon_west.dat',data_type=4,data_dims=[npixels+1,nlines+1])
	; Area values don't seem right, there are a lot of zero and missing values....
	area=read_binary(dir+'area_west.dat',data_type=4,data_dims=[npixels,nlines])

	MISSING = -999.9
	f = where(goes_clat eq MISSING)
	goes_clat(f) = !Values.F_NAN
	goes_clon(f) = !Values.F_NAN
	area(f)      = !Values.F_NAN
 	f = where(goes_elat eq MISSING)
 	goes_elat(f) = !Values.F_NAN
 	goes_elon(f) = !Values.F_NAN

	; Read in TEMPO lg1 data
	; ----------------------
	tempofile = '/nobackup/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4'
	ncid = ncdf_open(tempofile, /NOWRITE)

	ncdf_attget, ncid, 'Earth_radius', earth_rad, /GLOBAL
	ncdf_attget, ncid, 'Earth_flattening', flat, /GLOBAL

	varid = ncdf_varid(ncid, 'elat')
	ncdf_varget, ncid, varid, tempo_elat
	varid = ncdf_varid(ncid, 'elon')
	ncdf_varget, ncid, varid, tempo_elon	
	dims    = size(tempo_elon,/dimensions)
	tempo_nlons = dims(0)
	tempo_nlats = dims(1)	

	MISSING = max(tempo_elat)
	f = where(tempo_elat eq MISSING, count)
	if (count gt 0) then tempo_elat(f) = !Values.F_NAN
	f = where(tempo_elon eq MISSING ,count)
	if (count gt 0) then tempo_elon(f) = !Values.F_NAN

	; Trim GOES-R to TEMPO
	; --------------------
	for x=0,npixels-1 do begin
		column = reform(goes_clon(x,*))
		f = where(~finite(column,/NAN), count)
		if (count gt 0) then begin
			column = column(f)
			f = where(column ge max(tempo_elon(0,*),/NAN), countin)
			if (countin eq count) then break
		endif
	endfor
	xmin = x

	for x=npixels-1, 0, -1 do begin
		column = reform(goes_clon(x,*))
		f = where(~finite(column,/NAN), count)
		if (count gt 0) then begin
			column = column(f)
			f = where(column le min(tempo_elon(tempo_nlons-1,*),/NAN), countin)
			if (countin eq count) then break
		endif
	endfor		
	xmax = x

	for y=0,nlines-1 do begin
		row = reform(goes_clat(*,y))
		f = where(~finite(row,/NAN), count)
		if (count gt 0) then begin
			row = row(f)
			f = where(row ge max(tempo_elat(*,0),/NAN), countin)
			if (countin eq count) then break
		endif
	endfor	
	ymin = y

	for y=nlines-1, 0, -1 do begin
		row = reform(goes_clat(*,y))
		f = where(~finite(row,/NAN), count)
		if (count gt 0) then begin
			row = row(f)
			f = where(row le min(tempo_elat(*,tempo_nlats-1),/NAN), countin)
			if (countin eq count) then break
		endif
	endfor		
	ymax = y	

	; eye balling it here....using the hard limits cuts off too much of the domain
	xmin = xmin -1000
	xmax = xmax -1
	ymin = 500 
	goes_clat = goes_clat(xmin:xmax,ymin:ymax)
	goes_clon = goes_clon(xmin:xmax,ymin:ymax)
	area      = area(xmin:xmax,ymin:ymax)
	goes_elat = goes_elat(xmin:xmax+1,ymin:ymax+1)
	goes_elon = goes_elon(xmin:xmax+1,ymin:ymax+1)

	dims = size(goes_clon,/dimensions)
	npixels = dims(0)
	nlines  = dims(1)
	print, 'trimmed dimensions', npixels,nlines

	ew = reform(goes_clon(*, nlines*0.5))
	ns = reform(goes_clat(npixels*0.5, *))

	; Calculate pixel edge lenghts [km] using Vincenty Solution
	; ---------------------------------------------------------
	pix_size = fltarr(npixels,nlines,4)
	pix_size(*,*,*) = !Values.F_NAN
	a = earth_rad
	f = flat	
	radians = !DPI/180.0
	small   = 1e-12
	for x=0,npixels-1 do begin
		for y=0,nlines-1 do begin
			;make sure all 4 corners are defined
			nans = where(finite([goes_elat(x:x+1,y:y+1),goes_elon(x:x+1,y:y+1)],/nan), count)
			if (count gt 0) then begin
				goes_clon(x,y) = !Values.F_NAN
				goes_clat(x,y) = !Values.F_NAN
				area(x,y)      = !Values.F_NAN
			endif else if (count eq 0) then begin
				for e=0,3 do begin
					;bottom
					if (e eq 0) then begin
						phi1 = goes_elat(x,y)   
						phi2 = goes_elat(x+1,y)
						L    = goes_elon(x+1,y) - goes_elon(x,y)
					endif
					;right
					if (e eq 1) then begin
						phi1 = goes_elat(x+1,y)
						phi2 = goes_elat(x+1,y+1)
						L    = goes_elon(x+1,y) - goes_elon(x+1,y+1)
					endif
					;top
					if (e eq 2) then begin
						phi1 = goes_elat(x,y+1)
						phi2 = goes_elat(x+1, y+1)
						L    = goes_elon(x,y+1) - goes_elon(x+1, y+1)
					endif
					;left 
					if (e eq 3) then begin
						phi1 = goes_elat(x,y)
						phi2 = goes_elat(x,y+1)
						L    = abs(goes_elon(x,y) - goes_elon(x,y+1))
					endif
					
					phi1 = phi1*radians
					phi2 = phi2*radians
					L    = L*radians

					tanu1 = (1-f)*tan(phi1)
					tanu2 = (1-f)*tan(phi2)
					cosu1 = 1/sqrt(1+(tanu1^2))
					sinu1 = tanu1*cosu1
					cosu2 = 1/sqrt(1+(tanu2^2))
					sinu2 = tanu2*cosu2				

					lambda = L   ;first approximation
				  dlambda = 0.1
				  i = 0
					while ((dlambda gt small) or (i < 100)) do begin
						sintheta = sqrt((cosu2*sin(lambda))^2 + (cosu1*sinu2 - sinu1*cosu2*cos(lambda))^2)
						costheta = sinu1*sinu2 + cosu1*cosu2*cos(lambda)

						theta = atan(sintheta/costheta)
						sinalpha = cosu1*cosu2*sin(lambda)/sintheta
						cos2alpha = 1-sinalpha^2
						cos2thetam = costheta - 2*sinu1*sinu2/cos2alpha
						C = (f/16) * cos2alpha*(4 + f*(4-3*cos2alpha))
						lambda2 = L + (1-C)*f*sinalpha*(theta + C*sintheta*(cos2thetam + C*costheta*(-1+ 2*cos2thetam*cos2thetam)))
						dlambda = abs(lambda - lambda2)
						lambda  = lambda2
						i = i + 1
					endwhile

					if (e eq 100) then begin
						print, 'I=100'
						stop
					endif

					b = a*(1-f)
					usq = cos2alpha*(a*a - b*b)/(b*b)
					aprime = 1 + (usq/16384)*(4096 + usq*(-768+usq*(320 - 175*usq)))
					bprime = (usq/1024)*(256 + usq*(-128 + usq*(74 - 47*usq)))
					deltatheta = bprime*sintheta*(cos2thetam + 0.25*bprime*(costheta*(-1 + 2*cos2thetam*cos2thetam) - $
											 (bprime/6)*cos2thetam*(-3 + 4*sintheta*sintheta)*(-3 + 4*cos2thetam*cos2thetam)))
					pix_size(x,y,e) = b*aprime*(theta - deltatheta)
					
				endfor
			endif
		endfor
		if (~(x mod 100)) then begin
			print,x,npixels-1
		endif
	endfor

	MISSING = 1e15
	f = where(finite(goes_clon,/nan),count)
	if (count gt 0) then begin
		goes_clon(f) = MISSING
		goes_clat(f) = MISSING
	endif
	f = where(finite(goes_elon,/nan),count)
	if (count gt 0) then begin
		goes_elon(f) = MISSING
		goes_elat(f) = MISSING
	endif
	f = where(finite(pix_size,/nan),count)
	if (count gt 0) then begin
		pix_size(f) = MISSING
	endif
	f = where(finite(ew,/nan),count)
	if (count gt 0) then begin
		ew(f) = MISSING	
	endif
	f = where(finite(ns,/nan),count)
	if (count gt 0) then begin
		ns(f) = MISSING	
	endif

	; Create netcdf file
	; ------------------
	outfile = 'goes-r.lg1.invariant.nc4'
  ncid = ncdf_create(outfile,/clobber,/netcdf4_format)

  ; define dimensions        
  nsdimid = ncdf_dimdef(ncid,'ns',nlines)  
  ewdimid = ncdf_dimdef(ncid,'ew',npixels)
  nsedimid = ncdf_dimdef(ncid,'ns_e',nlines+1)
  ewedimid = ncdf_dimdef(ncid,'ew_e',npixels+1)
  corndimid = ncdf_dimdef(ncid,'pix_corner',4)
  timedimid = ncdf_dimdef(ncid,'time',1)

  ; write global attributes
  ncdf_attput, ncid, /GLOBAL, "institution", "NASA/Goddard Space Flight Center"
  ncdf_attput, ncid, /GLOBAL, "source", "Global Model and Assimilation Office"
  ncdf_attput, ncid, /GLOBAL, "history", "Created from center/edge_lon_west.dat and center/edge_lat_west.dat"
  ncdf_attput, ncid, /GLOBAL, "references", "none"
  ncdf_attput, ncid, /GLOBAL, "comment","This file contains GEOS-R geolocation information for a CONUS image trimmed to the TEMPO domain"
  ncdf_attput, ncid, /GLOBAL, "conventions","CF"
  ncdf_attput, ncid, /GLOBAL, "sat_lat",0.
  ncdf_attput, ncid, /GLOBAL, "sat_lon", -137.
  ncdf_attput, ncid, /GLOBAL, "sat_alt", 35800.
  ncdf_attput, ncid, /GLOBAL, "Earth_radius", earth_rad
  ncdf_attput, ncid, /GLOBAL, "Earth_flattening", flat
  ncdf_attput, ncid, /GLOBAL, "contact", "P. Castellanos<patricia.castellanos@nasa.gov>"

  ; define variables
  clatid = ncdf_vardef(ncid, "clat", [ewdimid, nsdimid], /double)
  ncdf_attput, ncid, clatid, "long_name"    , "pixel center latitude"
  ncdf_attput, ncid, clatid, "units"         , "degrees_north"
  ncdf_attput, ncid, clatid, "missing_value", 1e15

  clonid = ncdf_vardef(ncid, "clon", [ewdimid, nsdimid], /double)
  ncdf_attput, ncid, clonid, "long_name"    , "pixel center longitude"
  ncdf_attput, ncid, clonid, "units"         , "degrees_east"
  ncdf_attput, ncid, clonid, "missing_value", 1e15

  elatid = ncdf_vardef(ncid, "elat", [ewedimid, nsedimid], /double)
  ncdf_attput, ncid, elatid, "long_name"    , "latitude at pixel edge"
  ncdf_attput, ncid, elatid, "units"         , "degrees_north"
  ncdf_attput, ncid, elatid, "missing_value", 1e15

  elonid = ncdf_vardef(ncid, "elon", [ewedimid, nsedimid], /double)
  ncdf_attput, ncid, elonid, "long_name"    , "longitude at pixel edge"
  ncdf_attput, ncid, elonid, "units"         , "degrees_east"
  ncdf_attput, ncid, elonid, "missing_value", 1e15 

  ewid = ncdf_vardef(ncid, "ew", [ewdimid], /double)
  ncdf_attput, ncid, ewid, "long_name"    , "pseudo longitude"
  ncdf_attput, ncid, ewid, "units"         , "degrees_east"

  nsid = ncdf_vardef(ncid, "ns", [nsdimid], /double)
  ncdf_attput, ncid, nsid, "long_name"    , "pseudo latitude"
  ncdf_attput, ncid, nsid, "units"         , "degrees_north"

  sizeid = ncdf_vardef(ncid, "pix_size", [ewedimid, nsedimid, corndimid], /double)
  ncdf_attput, ncid, sizeid, "long_name"    , "pixel size"
  ncdf_attput, ncid, sizeid, "units"        , "km"
  ncdf_attput, ncid, sizeid, "comment"      , 'starts at south edge continues counterclockwise'

  timeid = ncdf_vardef(ncid, "time", [timedimid], /double)
  ncdf_attput, ncid, timeid, "long_name"    , "time"
  ncdf_attput, ncid, timeid, "units"         , "seconds since 2006-06-01 12:00:00"  
  ncdf_attput, ncid, timeid, "time_increment", 10000
  ncdf_attput, ncid, timeid, "begin_date", 20060601
  ncdf_attput, ncid, timeid ,"begin_time", 120000

  ; take file out of define mode
  ncdf_control, ncid, /ENDEF

  ; write data to file
  ;---------------------
  ncdf_varput, ncid, clatid, goes_clat
  ncdf_varput, ncid, clonid, goes_clon
  ncdf_varput, ncid, elatid, goes_elat
  ncdf_varput, ncid, elonid, goes_elon
  ncdf_varput, ncid, ewid, ew
  ncdf_varput, ncid, nsid, ns
  ncdf_varput, ncid, sizeid, pix_size
  ncdf_varput, ncid, timeid, [0]

	; close file
  ncdf_close,ncid

stop
end