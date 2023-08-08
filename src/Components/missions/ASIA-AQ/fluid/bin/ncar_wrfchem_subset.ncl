; extract variables from original wrfout files 
; Usage: ncl 'file_in="name_in.nc"' 'file_out="name_out.nc"' naaps_subset.ncl 
;----------------------------------------------------------------------------

; load in the libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

; set default values for file_in, dir_in, and file_out, if not specified
 if (.not.isvar("file_in")) then
   file_in = "../NCAR-WRFchem/wrfout_hourly_d01_2019-08-13_00:00:00"
 end if
 if (.not.isvar("file_out")) then
   file_out = "../NCAR-WRFchem/wrfout_hourly_d01_2019-08-13_00:00:00.subsetted.nc"
 end if

  f=addfile(file_in, "r")
  global_attnames = getvaratts(f)
  system("rm -f "+file_out) ;-- delete file
  setfileoption("nc","Format","NetCDF4") 
  fo = addfile(file_out,"c") ;-- create new

  ; global attributes 
  att=True
  do i=0,dimsizes(global_attnames)-1
     ;  print(f@$global_attnames(i)$)
      att@$global_attnames(i)$ = f@$global_attnames(i)$ 
  end do

  ;setfileoption(fo,"DefineMode",True)
  fileattdef(fo, att)
  filedimdef(fo,"time",-1,True) 
  

 vars = (/"Times","XLAT", "XLONG", "PSFC", "T", "U", "V", "W", "p", \
          "QVAPOR", "QCLOUD", "CLDFRA", "PBLH", "EXTCOF55",\
          "o3", "no", "no2", "no3", "PM2_5_DRY", \
          "so2", "sulf", "co", "ch4", "P25", "co_fire", "tr17_2", \
          "co_anth", "BC1", "BC2", "OC1", "OC2", "DUST_1", "DUST_2", "DUST_3", \
          "DUST_4", "DUST_5", "SEAS_1", "SEAS_2", "SEAS_3", "SEAS_4", "AOD_550"/) 


  ; write variables
  do i=0,dimsizes(vars)-1
     print("varname:   "+ vars(i))
     
     if (vars(i).eq."p") then 
         va := wrf_user_getvar(f,"p",-1)  ; Pa
        ; printVarSummary(va)
     else 
         v := f->$vars(i)$ 
         if (.not.ismissing(v@stagger) .and. strlen(v@stagger).gt.0)  then
              va := wrf_user_unstagger(v, v@stagger)
         else 
              va:=v 
         end if 
     end if
   
    ;output
    fo->$vars(i)$=va
  end do

  
  ext=flt2dble(f->EXTCOF55)  ; unit: km-1
  ph=f->PH   ;"perturbation geopotential"  (Time, bottom_top_stag, south_north, west_east)
  phb=f->PHB   ; "base state geopotential"  unit: m2 s-2
  d=dimsizes(ph)
  g=9.80665
  gph = flt2dble(ph+phb)/g
  dz= gph(:,1:d(1)-1,:,:)-gph(:,0:d(1)-2,:,:)


  copy_VarMeta(ext,dz)
  dz@description="model layer thickness"
  dz@units = "m"

  Rd=287.05 
  p = wrf_user_getvar(f,"p",-1)
  tk = wrf_user_getvar(f,"tk",-1)
  rho = p/Rd/tk 
  copy_VarMeta(ext,rho)
  rho@description="dry air density"
  rho@units = "kg m-3"
  
  

  ; append rho and dz 
  fo->dz=dz
  fo->rho=rho


