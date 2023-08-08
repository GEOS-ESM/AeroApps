; extract variables from original wrfout files 
; Usage: ncl 'file_in="inname.nc"' 'file_out="outname.nc"' arqi_regrid.ncl 
;----------------------------------------------------------------------------

; load in the libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


 ; default values for file_in, dir_in, and file_out, if not specified
  if (.not.isvar("file_in")) then
    file_in = "../ARQI/2019062912_001"
  end if
  if (.not.isvar("file_out")) then
    file_out = "../ARQI/2019062912_001_grid1.nc"
  end if

  f=addfile(file_in, "r")
  setfileoption("nc","Format","NetCDF4Classic") 
  setfileoption("nc","CompressionLevel",5) ; 0 ~ 9 
  system("rm -f "+file_out) ;-- delete file
  fo = addfile(file_out,"c")


  ; global attributes 
  global_attnames = getvaratts(f)
  att=True
  do i=0,dimsizes(global_attnames)-1
     ;  print(f@$global_attnames(i)$)
      att@$global_attnames(i)$ = f@$global_attnames(i)$ 
  end do
  fileattdef(fo, att)
  filedimdef(fo,"time",-1,True) 

  vars =  getfilevarnames(f)
  vars_skip = (/"lat_2", "lon_2", "rlat2", "rlon2"/)
  lat_1 = f->lat_1
  lon_1 = f->lon_1
  lat_2 = f->lat_2
  lon_2 = f->lon_2
  rlat1 = f->rlat1 
  rlon1 = f->rlon1 


  do i=0,dimsizes(vars)-1

     if (ismissing(str_match( vars_skip, vars(i) )))   then 
         print("processing:   "+ vars(i))
         v := f->$vars(i)$ 
         dims := getvardims(v)
         diml := dimsizes(v)

         if   (  .not.ismissing(str_match(dims,"rlat2")) \
           .and. .not.ismissing(str_match(dims,"rlon2")) )  then 
             ; re-grid from lat_2 lon_2 to lat_1 lon_1 
            
            ;---Set up options for regridding
            wgt_file="/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts/arqi_regrid_wgt.nc"
            Opt                = True
            Opt@WgtFileName    = wgt_file
            Opt@SrcGridLat     = lat_2     ; source grid
            Opt@SrcGridLon     = lon_2
            Opt@DstGridLat     = lat_1     ; destination grid
            Opt@DstGridLon     = lon_1

            Opt@SrcRegional    = True          ; Necessary if grids
            Opt@DstRegional    = True          ; are regional

            Opt@InterpMethod   = "bilinear"    ; "patch", "conserve"
            Opt@ForceOverwrite = True          ; Optional, but recommended.
            Opt@PrintTimings   = True          ; Optional
            Opt@Debug          = False          ; Optional

            ;---Do the regridding
              v_regrid := ESMF_regrid_with_weights(v, wgt_file, Opt)      ; Do the regridding

            delete(v_regrid@lat2d)
            delete(v_regrid@lon2d)
            nd = dimsizes(diml)
            v_regrid!(nd-1)="rlon1"
            v_regrid!(nd-2)="rlat1" 
            v_regrid&rlat1 = rlat1 
            v_regrid&rlon1 = rlon1 

            ; printVarSummary(v)
            ; printVarSummary(v_regrid)
             
            ; write output             
            fo->$vars(i)$ = v_regrid 

            else
        
            fo->$vars(i)$ = v
 
 
         end if 
     end if 
  end do

   
     

  

