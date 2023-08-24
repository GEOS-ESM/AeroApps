; extract variables from original wrfout files 
; Usage: ncl 'file_in="inname.nc"' 'file_out="outname.nc"' firework_regrid.ncl 
;----------------------------------------------------------------------------

; load in the libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


 ; default values for file_in, dir_in, and file_out, if not specified
  if (.not.isvar("file_in")) then
    file_in = "../FireWork/2019070900_000"
  end if
  if (.not.isvar("file_out")) then
    file_out = "../FireWork/2019070900_000_grid1.nc"
  end if

  f=addfile(file_in, "r")
  setfileoption("nc","Format","NetCDF4Classic") 
  setfileoption("nc","CompressionLevel",5) ; 0 ~ 9 
  system("rm -f "+file_out) ;-- delete file
  fo = addfile(file_out,"c")

  latlon_file="/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts/firework.lat_lon_1d.nc"
  f_ll=addfile(latlon_file, "r")
  lat_1d = f_ll->lat
  lon_1d = f_ll->lon 

;-------------------read and write ------------------------------------------
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
  ; vars_skip = (/"lat_2", "lon_2", "rlat2", "rlon2"/)
  vars_subset = (/"time", "pref",  \
                  "a_1","b_1","level1", \
                  "a_2","b_2","level2", \
                  "a_3","b_3","level3", \
                  "a_5","b_5","level5", \
                  "AF","AC","CO","N2","NO","O3","S2","C3H8","GZ_vgrid5", \
                  "HCHO","ISOP", "PAN","NH3","OH", "P0", \
                  "H","RHO","TT_vgrid5","UU_vgrid5","VV_vgrid5" /) 

  lat_1 = f->lat_1
  lon_1 = f->lon_1
  lat_2 = f->lat_2
  lon_2 = f->lon_2
  rlat1 = f->rlat1 
  rlon1 = f->rlon1 


  do i=0,dimsizes(vars)-1

     if (.not.ismissing(ind( vars_subset.eq. vars(i) )))   then 
         print("processing:   "+ vars(i))
         v := f->$vars(i)$ 
         dims := getvardims(v)
         diml := dimsizes(v)

         if   (      (.not.ismissing(str_match(dims,"rlat1")) .and. .not.ismissing(str_match(dims,"rlon1"))) \
               .or.  (.not.ismissing(str_match(dims,"rlat2")) .and. .not.ismissing(str_match(dims,"rlon2")))   )  then 
             ; re-grid to lat lon 
            
            ;---Set up options for regridding
            if  (.not.ismissing(str_match(dims,"rlat1"))) then 
               wgt_file="/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts/firework_regrid_wgt1.nc"
               Opt                = True
               Opt@WgtFileName    = wgt_file
               Opt@SrcGridLat     := lat_1     ; source grid
               Opt@SrcGridLon     := lon_1

            else 
               wgt_file="/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts/firework_regrid_wgt2.nc"
               Opt                = True
               Opt@WgtFileName    = wgt_file
               Opt@SrcGridLat     := lat_2     ; source grid
               Opt@SrcGridLon     := lon_2
            end if 

            Opt@DstGridLat     = lat_1d     ; destination grid
            Opt@DstGridLon     = lon_1d

            Opt@SrcRegional    = True          ; Necessary if grids
            Opt@DstRegional    = True          ; are regional

            Opt@InterpMethod   = "bilinear"    ; "patch", "conserve"
            Opt@ForceOverwrite = True          ; Optional, but recommended.
            Opt@PrintTimings   = True          ; Optional
            Opt@Debug          = False          ; Optional

            ;---Do the regridding
            if (.not.fileexists(wgt_file))  then 
               v_regrid := ESMF_regrid(v, Opt)  
            else 
               v_regrid := ESMF_regrid_with_weights(v, wgt_file, Opt)  
            end if 

            ;delete(v_regrid@lat2d)
            ;delete(v_regrid@lon2d)
            nd = dimsizes(diml)
            v_regrid!(nd-1)="lon"
            v_regrid!(nd-2)="lat" 
            v_regrid&lat = lat_1d 
            v_regrid&lon = lon_1d 

            ; printVarSummary(v)
            ; printVarSummary(v_regrid)
             
            ; write output             
            fo->$vars(i)$ = v_regrid 

            else
        
            fo->$vars(i)$ = v
 
 
         end if 
     end if 
  end do

   
     

  

