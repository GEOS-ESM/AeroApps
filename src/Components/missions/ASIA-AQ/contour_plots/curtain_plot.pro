  pro curtain_plot, caseid, yyyy, mmdd, pickvalue
;;;;;;;;i Huisheng Bian (06/16/2025) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; choose fields (caseid, yyyy, mmdd, pickvalue) to plot ;;;;;;;;;;;;;;;;;;;;
; 1. open idl
; 2. .run curtain_plot.pro
  3. go to line 147, change to your dirgeos directory similar like dirgeos = '/discover/nobackup/projects/gocart/hbian/proposal/2022-AsiaAQ/'+caseid+'/sampled/'
; 4. curtain_plot, 'aaqoxh24','2024','0206','DMS'

  ;caseid = 'aaqoxh24' or 'aaqm24crt'
  caseid = strlowcase(caseid)

  ; study which day
  ;yyyy = '2024'
  ;mmdd = ['0206','0207','0211','0213','0215','0217','0226','0308','0310','0311','0313','0316','0318','0321','0325','0327','0331']

  ; study which species
  ; pickvalue = ['DMS', 'SO2', 'SO4',  'MSA',  'H2O2','OH',  'NH3', 'HNO3','NO3',  'NH4',  'SS',   'BC',   'OA',   'RH']
  pickvalue =strupcase(pickvalue)

  if pickvalue eq 'DMS'  then begin 
     fvarname = 'DMS_TOGA_APEL'                  & unit ='pptv'
  endif   
  if pickvalue eq 'SO2'  then begin 
     fvarname = 'SO2_CIT_CROUNSE'                & unit ='pptv'
  endif   
  if pickvalue eq 'SO4'  then begin
     fvarname = 'Sulfate_PM1.5_AMS_60s_JIMENEZ'  & unit ='ug/m3'
  endif   
  if pickvalue eq 'MSA'  then begin 
     fvarname = 'MSA_PM1.5_AMS_60s_JIMENEZ'      & unit ='ug/m3'
  endif   
  if pickvalue eq 'H2O2' then begin 
     fvarname = 'H2O2_CIT_CROUNSE'               & unit ='pptv'
  endif   
  if pickvalue eq 'NH3'  then begin 
     fvarname = 'NH3_ppbv_LEE'                   & unit ='ppbv'
  endif   
  if pickvalue eq 'HNO3' then begin 
     fvarname = 'HNO3_CIT_CROUNSE'               & unit ='pptv'
  endif   
  if pickvalue eq 'NO3'  then begin 
     fvarname = 'Nitrate_PM1.5_AMS_60s_JIMENEZ'  & unit ='ug/m3'
  endif   
  if pickvalue eq 'NH4'  then begin 
     fvarname = 'Ammonium_PM1.5_AMS_60s_JIMENEZ' & unit ='ug/m3'
  endif   
  if pickvalue eq 'SS'   then begin 
     fvarname = 'Seasalt_PM1.5_AMS_60s_JIMENEZ'  & unit ='ug/m3'
  endif   
  if pickvalue eq 'BC'   then begin 
     fvarname = 'BC_MassConc_LEE'                & unit ='ng/m3'
  endif   
  if pickvalue eq 'OA'   then begin 
     fvarname = 'OA_PM1.5_AMS_60s_JIMENEZ'       & unit ='ug/m3'
  endif   
  if pickvalue eq 'RH'   then begin 
     fvarname = 'RHw_DLH_DISKIN'                 & unit =' '
  endif   

; Get the Oracle observation data 
   close,1
   dirobs = '/discover/nobackup/projects/gmao/iesa/aerosol/campaigns/ASIA-AQ/archive/merges/v20250423/'
   fdatain = file_search(dirobs+'asiaaq-mrg60_dc8_'+yyyy+mmdd+'_RA_20250423T*.ict')
   print, dirobs+fdatain

   headv =['Time_Start','Day_Of_Year_BENNETT','Latitude_BENNETT','Longitude_BENNETT','Altitude_AGL_m_DIGANGI','P','T', $
           'StdtoVol_AMS_60s_JIMENEZ', fvarname]
         nheadv = n_elements(headv)
         print,'nheadv =',nheadv
         locvarv = intarr(nheadv)
         ll = 3000L
         fvarv = 0.d0*fltarr(ll,nheadv)
   close,1
   openr,1,fdatain
   s1=''
       readf,1,s1
       tmpstr = strsplit(s1,',',/extract)
       ntmpstr = n_elements(tmpstr)
       ls = float(tmpstr(0))
       print,'ls = ',ls
   for i=0,ls-2 do begin    
       readf,1,s1
   endfor
               totv = strsplit(s1,' ,',/extract)
               print,'totv =',totv
               if totv(0) eq 'Time_Start' and totv(1) eq 'Time_Stop' then begin
                  ntotv = n_elements(totv)
                  print,'ntotv =',ntotv
                  for i = 0,ntotv-1 do begin
                  for k = 0,nheadv-1 do begin
                      if totv(i) eq headv(k) then begin
                         locvarv(k) = i
                      endif
                  endfor
                  endfor
                  print,'locvarv =',locvarv
               endif

         j = 0
         while(not(eof(1))) do begin 
               readf, 1, s1
               trcv = strsplit(s1,',',/extract)
               for k=0,nheadv-1 do begin
                   fvarv(j,k) = float(trcv(locvarv(k)))
               endfor
               j = j+1L
         endwhile

   print,'size = ',size(fvarv,/dimensions)
   print, 'j =',j

   sizeobs = j
      ;oyyyymmdd = dblarr(sizeobs) & oyyyymmdd(*) = yyyy+mmdd
      print,'sizeobs =',sizeobs
      yyyymmdd = dblarr(sizeobs) & yyyymmdd(*) = 0.0
      utc = fltarr(sizeobs) & utc(*) = 0.0
      fltlat = fltarr(sizeobs) & fltlat(*) = 0.0
      fltlon = fltarr(sizeobs) & fltlon(*) = 0.0
      fltalt = fltarr(sizeobs) & fltalt(*) = 0.0
      fltpress = fltarr(sizeobs) & fltpress(*) = 0.0
      flttemp = fltarr(sizeobs) & flttemp(*) = 0.0
      fltstdtovol = fltarr(sizeobs) & fltstdtovol(*) = 0.0
      fltval = fltarr(sizeobs) & fltval(*) = 0.0

          utc = fvarv(0:j-1,0)+30
          fltlat = fvarv(0:j-1,2)
          fltlon = fvarv(0:j-1,3)
          fltalt = fvarv(0:j-1,4)/1000.  ; from ft to km
          fltpress = fvarv(0:j-1,5)  ;  ; hPa 
          flttemp = fvarv(0:j-1,6)  ;  K 
          fltstdtovol = fvarv(0:j-1,7)
          fltval = fvarv(0:j-1,8)

      fltpress[where(fltpress lt -999.00)] = !values.f_nan
      flttemp[where(flttemp lt -999.00)] = !values.f_nan
      if pickvalue eq 'SO4' or pickvalue eq 'MSA' or pickvalue eq 'NO3' or pickvalue eq 'NH4' or pickvalue eq 'SS' or pickvalue eq 'OA' then begin
         fltval = fltval * fltstdtovol   ; from ug sm-3 to ug m-3
      endif
      fltval[where(fltval le 0 or finite(fltval) ne 1)] = !values.f_nan
      print,'fltval =',max(fltval,/nan),min(fltval,/nan)

  ;pstd = 101325 ; hPa
  ;tstd = 273.15  ;K
  ;vstd = fltpress/pstd*tstd/flttemp

  dirgeos = '/discover/nobackup/projects/gocart/hbian/proposal/2022-AsiaAQ/'+caseid+'/sampled/'
  filemodel = file_search(dirgeos+'asiaaq-mrg60_dc8_'+yyyy+mmdd+'_RA_20250423T*.inst3_aer_Nv.nc4')
  print,'filemodel =',filemodel

  cdfid = ncdf_open(filemodel)
  id = ncdf_varid(cdfid,'time')
  ncdf_varget, cdfid, id, time
  id = ncdf_varid(cdfid,'lon')
  ncdf_varget, cdfid, id, lon
  id = ncdf_varid(cdfid,'lat')
  ncdf_varget, cdfid, id, lat
  id = ncdf_varid(cdfid,'AIRDENS')
  ncdf_varget, cdfid, id, rhoa   ; kg m-3
  id = ncdf_varid(cdfid,'BCPHOBIC')
  ncdf_varget, cdfid, id,bc1 
  id = ncdf_varid(cdfid,'BCPHILIC')
  ncdf_varget, cdfid, id,bc2 
  id = ncdf_varid(cdfid,'BRPHOBIC')
  ncdf_varget, cdfid, id,br1 
  id = ncdf_varid(cdfid,'BRPHILIC')
  ncdf_varget, cdfid, id,br2 
  id = ncdf_varid(cdfid,'DMS')   ; kg kg-1
  ncdf_varget, cdfid, id, dms
  id = ncdf_varid(cdfid,'DU001')  ; radius: 0.1 - 1.0 um
  ncdf_varget, cdfid, id,du001
  id = ncdf_varid(cdfid,'DU002')  ; radious: 1.0 - 1.8 um
  ncdf_varget, cdfid, id,du002
  id = ncdf_varid(cdfid,'H')  ; m
  ncdf_varget, cdfid, id,Hmls
  id = ncdf_varid(cdfid,'HNO3CONC')   ; kg m-3
  ncdf_varget, cdfid, id, hno3
  id = ncdf_varid(cdfid,'MSA')
  ncdf_varget, cdfid, id, msa
  id = ncdf_varid(cdfid,'NH3')   ; kg kg-1
  ncdf_varget, cdfid, id, nh3
  id = ncdf_varid(cdfid,'NH4a')
  ncdf_varget, cdfid, id,nh4 
  id = ncdf_varid(cdfid,'NO3an1')
  ncdf_varget, cdfid, id,no3an1 
  id = ncdf_varid(cdfid,'OCPHOBIC')
  ncdf_varget, cdfid, id,oc1 
  id = ncdf_varid(cdfid,'OCPHILIC')
  ncdf_varget, cdfid, id,oc2 
  id = ncdf_varid(cdfid,'RH')
  ncdf_varget, cdfid, id, rhm  ; K from top to bottom
  id = ncdf_varid(cdfid,'SO2')   ; kg kg-1
  ncdf_varget, cdfid, id, so2
  id = ncdf_varid(cdfid,'SO4')   ; kg kg-1
  ncdf_varget, cdfid, id, so4
  id = ncdf_varid(cdfid,'SS001')  ; radius: 0.03 - 0.1 um
  ncdf_varget, cdfid, id,ss001
  id = ncdf_varid(cdfid,'SS002')  ; radius: 0.1 - 0.5 um
  ncdf_varget, cdfid, id,ss002
  id = ncdf_varid(cdfid,'SS003')  ; radius: 0.5 - 1.5 um
  ncdf_varget, cdfid, id,ss003
  id = ncdf_varid(cdfid,'SUXH2O2')   ; mol mol-1
  ncdf_varget, cdfid, id, h2o2
  id = ncdf_varid(cdfid,'SUXNO3')   ; mol mol-1
  ncdf_varget, cdfid, id, no3
  id = ncdf_varid(cdfid,'SUXOH')   ; mol mol-1
  ncdf_varget, cdfid, id, oh
  id = ncdf_varid(cdfid,'T')
  ncdf_varget, cdfid, id, temp  ; K from top to bottom
  id = ncdf_varid(cdfid,'delp')
  ncdf_varget, cdfid, id, delp
  ncdf_close, cdfid

  h = delp/9.8/rhoa/1000.d0   ; km
     dmsconc  = fltarr(72,sizeobs) & dmsconc  = dms*28.96/78.13*1.e12     ; from kg/kg to pptv
     so2conc  = fltarr(72,sizeobs) & so2conc  = so2*28.96/64.07*1.e12     ; from kg/kg to pptv
     so4conc  = fltarr(72,sizeobs) & so4conc  = so4*rhoa*1.e9             ; kg k-1 to ug m-3
     msaconc  = fltarr(72,sizeobs) & msaconc  = msa*rhoa*1.e9             ; kg k-1 to ug m-3
     h2o2conc = fltarr(72,sizeobs) & h2o2conc = h2o2 *1.e12               ; from mol/mol to pptv
     nh3conc  = fltarr(72,sizeobs) & nh3conc  = nh3*28.96/17.031*1.e9     ; from kg/kg to ppbv
     no3conc  = fltarr(72,sizeobs) & no3conc  = no3an1*rhoa*1.e9          ; kg k-1 to ug m-3
     hno3conc = fltarr(72,sizeobs) & hno3conc = hno3 *1.e12               ; from mol/mol to pptv
     nh4conc  = fltarr(72,sizeobs) & nh4conc  = nh4*rhoa*1.e9             ; kg k-1 to ug m-3
     ssconc   = fltarr(72,sizeobs) & ssconc   = (ss001+ss002)*rhoa*1.e9   ; PM1, kg k-1 to ug m-3
     bcconc   = fltarr(72,sizeobs) & bcconc   = (bc1+bc2)*rhoa*1.e12      ; kg kg-1 to ng m-3
     orgconc  = fltarr(72,sizeobs) & orgconc  = (oc1+oc2+br1+br2)*rhoa*1.e9  ; kg k-1 to ug m-3
     rhconc   = fltarr(72,sizeobs) & rhconc   = rhm * 100.                ; from 0-1 to 0-100 %

     sizemod = size(rhoa,/dimensions)

; make sure obs and mod have the same size
  if sizeobs ne sizemod(1) then begin
     print,'sizeobs, sizemod =', sizeobs, sizemod(1)
     stop
  endif

  ;for io = 0,1 do begin  ; 0 for Pacific and 1 for Atlantic
  xlat = indgen(sizeobs)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  tmpph = fltarr(73,sizeobs)
  tmpph(0,*) = 0.d0
  for k = 0,71 do begin
      tmpph(k+1,*) = tmpph(k,*) + h(71-k,*)
  endfor
  peh = fltarr(73,sizeobs)
  peh = reverse(tmpph) ; now peha from top to bottom in km
  ph = fltarr(72,sizeobs)
  for k = 0,71 do begin
     ph(k,*) = 0.5*(peh(k,*)+peh(k+1,*))
  endfor

  lonv = fltarr(sizeobs) & latv = fltarr(sizeobs) & fltaltv = fltarr(sizeobs)
  lonv = lon(*) & latv = lat(*)+90. & fltaltv = fltalt(*)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; all data from top to bottom
      pval = fltarr(72,sizeobs)
      if pickvalue eq 'DMS'  then pval = dmsconc(*,*)   ; pptv
      if pickvalue eq 'SO2'  then pval = so2conc(*,*)   ; pptv
      if pickvalue eq 'SO4'  then pval = so4conc(*,*)
      if pickvalue eq 'MSA'  then pval = msaconc(*,*)
      if pickvalue eq 'H2O2' then pval = h2o2conc(*,*) ; pptv
      if pickvalue eq 'NH3'  then pval = nh3conc(*,*)   ; ppbv
      if pickvalue eq 'NO3'  then pval = no3conc(*,*)
      if pickvalue eq 'HNO3' then pval = hno3conc(*,*) ; pptv
      if pickvalue eq 'NH4'  then pval = nh4conc(*,*)
      if pickvalue eq 'SS'   then pval = ssconc(*,*)
      if pickvalue eq 'BC'   then pval = bcconc(*,*)
      if pickvalue eq 'OA'   then pval = orgconc(*,*)
      if pickvalue eq 'RH'   then pval = rhconc(*,*)     ; % 

; Make a map of the traj
  set_plot, 'ps'
  fign = caseid+'_'+pickvalue+'_'+mmdd+'.ps'
  !psym=0
  device, file=fign, /color, /helvetica, font_size=12, $
        /landscape
  !p.font=0
          z = [findgen(16),0.]*(!pi*2/16.)
          usersym,cos(z),sin(z),/fill
 !p.color=1

  xx = indgen(sizeobs)
  xtk = [0,50,100,150,200,250,300,350,400,450,500] 
  if mmdd eq '0215' or mmdd eq '0217' or mmdd eq '0226'  then xtk = [0,50,100,150,200,250,300,350,400,450] 
  if mmdd eq '0310' then xtk = [0,50,100,150,200,250,300,350,400] 
  if mmdd eq '0311' then xtk = [0,50,100,150,200,250,300,350,400,450,500,550] 

  nxtkname = n_elements(xtk)
  if xtk(nxtkname-1) ge sizeobs then xtk(nxtkname-1) = sizeobs
  xr1 = 1 & xr2 = xtk(nxtkname-1)
  xnone = replicate(' ',nxtkname)

  cfltlon = strarr(nxtkname)
  cfltlat = strarr(nxtkname)
  for np = 0,nxtkname-1 do begin
      if np le nxtkname-2 then begin
      if fltlon(xtk(np)) lt 0.d0 then cfltlon(np) = string(abs(fltlon(xtk(np))),format='(f5.1)')+'W'
      if fltlon(xtk(np)) ge 0.d0 then cfltlon(np) = string(abs(fltlon(xtk(np))),format='(f5.1)')+'E'
      if fltlat(xtk(np)) lt 0.d0 then cfltlat(np) = string(abs(fltlat(xtk(np))),format='(f4.1)')+'S'
      if fltlat(xtk(np)) ge 0.d0 then cfltlat(np) = string(abs(fltlat(xtk(np))),format='(f4.1)')+'N'
      endif
      if np eq nxtkname-1 then begin
      if fltlon(sizeobs-1) lt 0.d0 then cfltlon(np) = string(abs(fltlon(sizeobs-1)),format='(f5.1)')+'W'
      if fltlon(sizeobs-1) ge 0.d0 then cfltlon(np) = string(abs(fltlon(sizeobs-1)),format='(f5.1)')+'E'
      if fltlat(sizeobs-1) lt 0.d0 then cfltlat(np) = string(abs(fltlat(sizeobs-1)),format='(f4.1)')+'S'
      if fltlat(sizeobs-1) ge 0.d0 then cfltlat(np) = string(abs(fltlat(sizeobs-1)),format='(f4.1)')+'N'
      endif
  endfor
  xtkname = strarr(12)
  xtkname(0) = cfltlon(0)+', '+cfltlat(0)
  xtkname(1) = cfltlon(1)+', '+cfltlat(1)
  xtkname(2) = cfltlon(2)+', '+cfltlat(2)
  xtkname(3) = cfltlon(3)+', '+cfltlat(3)
  xtkname(4) = cfltlon(4)+', '+cfltlat(4)
  xtkname(5) = cfltlon(5)+', '+cfltlat(5)
  xtkname(6) = cfltlon(6)+', '+cfltlat(6)
  xtkname(7) = cfltlon(7)+', '+cfltlat(7)
  xtkname(8) = cfltlon(8)+', '+cfltlat(8)
  if mmdd ne '0310' then xtkname(9) = cfltlon(9)+', '+cfltlat(9)
  if mmdd ne '0215' and mmdd ne '0217' and mmdd ne '0226' and mmdd ne '0310'  then xtkname(10) = cfltlon(10)+', '+cfltlat(10)
  if mmdd eq '0311' then xtkname(11) = cfltlon(11)+', '+cfltlat(11)

; Curtain
  loadct, 39
  yr2 = 4
  if mmdd eq '0215' or mmdd eq '0313' or mmdd eq '0327' then yr2 = 10
  ; draw plot frame
   htit='   ASIAAQ & GEOS '+pickvalue+' '+unit
  print,'max min pval =',max(pval),min(pval)
  contour, transpose(pval), transpose(xlat), transpose(ph), $
   yrange=[0,yr2],xrange=[xr1,xr2], /nodata, /noerase, $
   xtickv=xtk,xtickname=xnone,xticks=nxtkname-1,xminor=0, $
   xticklen=-0.02,yticklen=-0.01, $
   xtitle='', ytitle='Altitude [km]', $
   title= htit, $
   position=[.10,.25,.95,.6]
   
   ypos = Replicate(!Y.Window[0] -0.08, nxtkname)
   xpos = !X.Window[0] - 0.04 + (!X.Window[1] - !X.Window[0]) * $
         Findgen(nxtkname) / (nxtkname-1)
   print,'ypos =',ypos
   print,'xpos =',xpos
   for np = 0, nxtkname-1 do begin
       xyouts,xpos(np),ypos(np),xtkname(np), Alignment=0.0, Orientation=25, /Normal
   endfor

      maxv = max(pval,/NAN) & minv = min(pval,/NAN)
      print,'maxv, minv =',maxv, minv
         levelsb=  [0.001,0.002,0.005,0.01,0.03,0.05,0.1,0.3,0.5,1]

         if pickvalue eq 'DMS' then levelsb=  [0.0001,0.001,0.01,0.1,1,2,5,10,20,50]
         if pickvalue eq 'SO2' then levelsb=  [1,2,5,10,50,100,200,1000,5000,10000]
         if pickvalue eq 'SO4' then begin
            levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]
            if mmdd eq '0313' or mmdd eq '0316' or mmdd eq '0318' or mmdd eq '0321' or mmdd eq '0325' or mmdd eq '0327' then levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*2
         endif   
         if pickvalue eq 'MSA' then begin
            levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*0.025
            if mmdd eq '0211' then levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*0.01 
         endif   
         if pickvalue eq 'H2O2' then begin
            levelsb=  [10,30,50,70,100,200,500,1000,5000,10000]
            if mmdd eq '0206' or mmdd eq '0207' or mmdd eq '0211' or mmdd eq '0213' then begin
               levelsb=  [100,500,700,1000,1300,1500,1700,2000,3000,5000]
            endif
            if mmdd eq '0316' or mmdd eq '0318' or mmdd eq '0321' or mmdd eq '0325' or mmdd eq '0327' then begin
               levelsb=  [100,500,1000,1500,1700,2000,2500,3000,5000,10000]
            endif
         endif



         if pickvalue eq 'NH3' then levelsb=  [0.0001,0.001,0.01,0.1,0.2,0.5,1,5,10,100]
         if pickvalue eq 'HNO3' then begin
            levelsb=  [1,5,10,30,50,100,200,500,1000,2000]
            if mmdd eq '0318' or mmdd eq '0321' or mmdd eq '0325' or mmdd eq '0327' then begin
               levelsb=  [1,50,100,200,500,1000,2000,3000,5000,10000]
            endif
         endif   
         if pickvalue eq 'SS' then levelsb=  [0.001,0.002,0.005,0.01,0.03,0.05,0.1,0.3,0.5,1]
         if pickvalue eq 'SS' then begin
            if mmdd eq '0321' then levelsb=  [0.001,0.002,0.005,0.007,0.01,0.03,0.05,0.1,0.5,1]
            if mmdd eq '0206' or mmdd eq '0207' then levelsb=  [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,10]
            if mmdd eq '0211' or mmdd eq '0213' then levelsb=  [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,10]
            if mmdd eq '0327' then levelsb=  [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,10]
         endif
         if pickvalue eq 'OA' then begin
            levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]
            if mmdd eq '0211' or mmdd eq '0213' or mmdd eq '0215' then levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*2
            if mmdd eq '0217' or mmdd eq '0226' or mmdd eq '0308' or mmdd eq '0310' or mmdd eq '0311' or mmdd eq '0313' then levelsb=  [0.01,0.1,0.25,0.5,1,2,3,4,5,6]*3
            if mmdd eq '0316' or mmdd eq '0318' or mmdd eq '0321' or mmdd eq '0325' or mmdd eq '0327' then levelsb=  [0.01,0.1,0.25,0.5,1,2,3,4,5,6]*5
         endif
         if pickvalue eq 'NO3' then begin
            levelsb=  [0.001,0.01,0.1,0.3,0.5,1,1.5,2,3,4]*0.1
            if  mmdd eq '0217' or mmdd eq '0226' or mmdd eq '0310' or mmdd eq '0311'  then levelsb=  [0.001,0.01,0.1,0.3,0.5,1,1.5,2,3,4]*5
            if  mmdd eq '0308' then levelsb=  [0.001,0.01,0.1,0.3,0.5,1,1.5,2,3,4]*1
            if  mmdd eq '0313' or mmdd eq '0316' or mmdd eq '0318' or mmdd eq '0321' or mmdd eq' 0325' or mmdd eq '0327' then levelsb=  [0.001,0.01,0.1,0.3,0.5,1,1.5,2,3,4]*2
         endif
         if pickvalue eq 'NH4' then begin
            levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*0.5
            if mmdd eq '0217' or mmdd eq '0226' or mmdd eq '0313' or mmdd eq '0316' or mmdd eq '0318' or mmdd eq '0321' or mmdd eq '0325' or mmdd eq '0327' then levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*2.0
            if mmdd eq '0310' or mmdd eq '0311' then levelsb=  [0.01,0.1,0.3,0.5,1,1.5,2,2.5,3,4]*1.0
         endif
         if pickvalue eq 'BC' then begin
            levelsb=  [1,5,10,20,30,50,100,500,1000,2000]
            if mmdd eq '0215' or mmdd eq '0318' then begin
            levelsb=  [1,5,10,20,50,100,1000,5000,10000,50000]
            endif
            if mmdd eq '0321' or mmdd eq '0325' then begin
            levelsb=  [1,10,20,50,100,200,500,1000,2000,5000]
            endif
         endif
         if pickvalue eq 'RH' then begin
            levelsb=  [0.001,0.01,0.1,0.3,0.5,0.6,0.7,0.8,0.9,1]*100.
         endif
         nlevelsb = n_elements(levelsb)
         ;colors = [30,70,90,120,140,160,180,210,254]
         colors = [30,40,70,90,110,170,190,210,254]

  ;; draw modelcontour
  print,'GEOS var =',max(pval),min(pval)
  contour, /over, transpose(pval), transpose(xlat), transpose(ph), $
   levels=levelsb, c_col=colors, /cell, /fill

; Pull the model values at altitude
  xlatflt = fltarr(sizeobs) & xlatflt = xlat(*)
      print,'size fltval =',size(fltval,/dimensions)
      for ip =0L,sizeobs-1 do begin   ; 0927
      for ic = 0, nlevelsb-2 do begin
          if( fltval(ip) ge levelsb(ic) and fltval(ip) le levelsb(ic+1)) then begin
            plots, xlatflt(ip), fltaltv(ip), psym=8, symsize=1.05, color=0
          endif
      endfor
      endfor
      for ip =0L,sizeobs-1 do begin    ; 0927 
      for ic = 0, nlevelsb-2 do begin
          if( fltval(ip) ge levelsb(ic) and fltval(ip) le levelsb(ic+1)) then begin
            plots, xlatflt(ip), fltaltv(ip), psym=8, symsize=0.50, col=colors(ic)
          endif
      endfor
      endfor

  lab = strarr(nlevelsb)
  for i = 0, nlevelsb-1 do begin
      if levelsb(i) ge 0.0001 then lab(i) = string(levelsb(i),format='(f6.4)')
      if levelsb(i) ge 0.001 then lab(i) = string(levelsb(i),format='(f5.3)')
      if levelsb(i) ge 0.01 then lab(i) = string(levelsb(i),format='(f4.2)')
      if levelsb(i) ge 0.1 then lab(i) = string(levelsb(i),format='(f3.1)')
      if levelsb(i) ge 1 then lab(i) = string(levelsb(i),format='(f3.1)')
      if levelsb(i) ge 10 then lab(i) = string(levelsb(i),format='(i3)')
      if levelsb(i) ge 100 then lab(i) = string(levelsb(i),format='(i4)')
      if levelsb(i) ge 1000 then lab(i) = string(levelsb(i),format='(i5)')
  endfor

  makekey, .14, .11, .72, .03, 0, -.035, align=0, $
   colors=make_array(nlevelsb-1,val=255), labels=lab(0:nlevelsb-1)

   ; fill in color
  makekey, .14, .11, .72, .03, 0, -.035, $
   colors=colors, labels=make_array(nlevelsb,val=' ')

  print,'fign =',fign
  device, /close 

end
