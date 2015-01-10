! 
!     Performs test of the I/O interface to HDF
!     Observational Data Stream files.
!
! !REVISION HISTORY: 
!     10May2000 C. Redder  Fixed subscript bug in call to ODS_Info
!                          by predefining n_kt2, n_kx2, and n_qcx2.
!
      implicit none
c      integer     ktmax
c      parameter ( ktmax      = 7 )
c      integer     kxmax
c      parameter ( kxmax      = 87 )

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

c      integer     first_jday
c      parameter ( first_jday = 18000 )
      integer     str_len
      parameter ( str_len     = 50 )
      integer     mobs
      parameter ( mobs        = 65 )
      integer n_kt1,  n_kt2
      integer n_kx1,  n_kx2
      integer n_qcx1, n_qcx2

      character * (str_len) kt_names1  ( Dnkt  + 1 )
      character * (str_len) kt_units1  ( Dnkt  + 1 )
      character * (str_len) kx_names1  ( Dnkx  + 1 )
      character * (str_len) kx_meta1   ( Dnkx  + 1 )
      character * (str_len) qcx_names1 ( Dnqcx + 1 )
      character * (str_len) kt_names2  ( Dnkx  + 1 )
      character * (str_len) kt_units2  ( Dnkt  + 1 )
      character * (str_len) kx_names2  ( Dnkx  + 1 )
      character * (str_len) kx_meta2   ( Dnkx  + 1 )
      character * (str_len) qcx_names2 ( Dnqcx + 1 )
      character * (str_len) filename
      character * (str_len) version
      character * (str_len) mode
      character * (str_len) histag
      character * (str_len) type

      integer nobs1,   nobs2
      integer jday1,   jday2
      integer stime1,  stime2
      integer id1,     id2
      integer nc_id,   ktmaxid

      real    lat1     ( mobs ), lat2     ( mobs )
      real    lon1     ( mobs ), lon2     ( mobs )
      real    level1   ( mobs ), level2   ( mobs )
      integer kx1      ( mobs ), kx2      ( mobs )
      integer kt1      ( mobs ), kt2      ( mobs )
      real    xm1      ( mobs ), xm2      ( mobs )
      integer caldate1 ( mobs ), caldate2 ( mobs )
      integer caltime1 ( mobs ), caltime2 ( mobs )
      integer time1    ( mobs ), time2    ( mobs )
      real    obs1     ( mobs ), obs2     ( mobs )
      real    sigo1    ( mobs ), sigo2    ( mobs )
      integer qc_excl1 ( mobs ), qc_excl2 ( mobs )
      integer qc_hist1 ( mobs ), qc_hist2 ( mobs )

      integer nobs3,   nobs4
      integer jday3,   jday4
      integer stime3,  stime4
      integer id3,     id4

      real    lat3     ( mobs ), lat4     ( mobs )
      real    lon3     ( mobs ), lon4     ( mobs )
      real    level3   ( mobs ), level4   ( mobs )
      integer kx3      ( mobs ), kx4      ( mobs )
      integer kt3      ( mobs ), kt4      ( mobs )
      real    xm3      ( mobs ), xm4      ( mobs ), km4     ( mobs )
      integer caldate3 ( mobs ), caldate4 ( mobs )
      integer caltime3 ( mobs ), caltime4 ( mobs ), julian4 ( mobs )
      integer time3    ( mobs ), time4    ( mobs )
      real    obs3     ( mobs ), obs4     ( mobs )
      real    sigo3    ( mobs ), sigo4    ( mobs )
      integer qc_excl3 ( mobs ), qc_excl4 ( mobs )
      integer qc_hist3 ( mobs ), qc_hist4 ( mobs )

      integer first_jday, first_jday2
      integer first_date, first_date2
      integer ODS_Julian
      integer ODS_CalDat
      integer ODS_StrSize
      integer ODS_NGet
      integer len_str
      integer ierr
      integer ndim
      parameter ( ndim = 3 )

      integer nblocks, iblock 
      integer ibeg, iend
      integer start ( ndim ), count ( ndim )
      integer start_wsp ( ndim ), count_wsp ( ndim )
      integer latest_jday, latest__hour

      histag        = 'test'

      n_kt1  = Dnkt
      n_kx1  = Dnkx
      n_qcx1 = Dnqcx

      call SET_K1   ( n_kt1,  kt_names1, kt_units1,
     .                n_kx1,  kx_names1, kx_meta1,
     .                n_qcx1, qcx_names1 )

      first_date  = 19680523
      first_jday  = ODS_Julian ( first_date )
      first_date2 = ODS_CalDat ( first_jday )

      call ODS_Create ( id1, 'test1.nc','pre_anal', 
     .                  first_jday,
     .                  n_kt1,  kt_names1, kt_units1,
     .                  n_kx1,  kx_names1, kx_meta1,
     .                  n_qcx1, qcx_names1,
     .                  ierr )

      nobs1   = 20
      jday1   = first_jday
      stime1  = 0

      call Set_OBS1 ( nobs1,
     .                lat1, lon1, level1,
     .                jday1,   stime1,
     .                caldate1, caltime1,
     .                kx1, kt1, xm1,
     .                obs1, sigo1, qc_excl1, qc_hist1 )

      call ODS_Cal2Time ( id1, nobs1, caldate1, caltime1,
     .                                time1, ierr )

      call ODS_PutR ( id1,  'lat',     jday1, stime1,
     .                nobs1, lat1,     ierr )
      call ODS_PutR ( id1,  'lon',     jday1, stime1,
     .                nobs1, lon1,     ierr )
      call ODS_PutR ( id1,  'lev',     jday1, stime1,
     .                nobs1, level1,   ierr )
      call ODS_PutI ( id1,  'time',    jday1, stime1,
     .                nobs1, time1,    ierr )
      call ODS_PutI ( id1,  'kt',      jday1, stime1,
     .                nobs1, kt1,      ierr )
      call ODS_PutI ( id1,  'kx',      jday1, stime1,
     .                nobs1, kx1,      ierr )
      call ODS_PutR ( id1,  'xm',      jday1, stime1,
     .                nobs1, xm1,      ierr )
      call ODS_PutR ( id1,  'obs',     jday1, stime1,
     .                nobs1, obs1,     ierr )
c      call ODS_PutR ( id1,  'sigO',    jday1, stime1,
c     .                nobs1, sigo1,    ierr )
      call ODS_PutI ( id1,  'qcexcl',  jday1, stime1,
     .                nobs1, qc_excl1, ierr )
      call ODS_PutI ( id1,  'qchist',  jday1, stime1,
     .                nobs1, qc_hist1, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 11 )
 11   format ( /' Obs set 1 written to file, test1.nc ' )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt1,  kt_names1, kt_units1,
     .                   n_kx1,  kx_names1, kx_meta1,
     .                   n_qcx1, qcx_names1,
     .                   nobs1,
     .                   lat1, lon1, level1,
     .                   caldate1, caltime1,
     .                   kx1, kt1, xm1,
     .                   obs1, sigo1, qc_excl1, qc_hist1 )

      nobs2   = 25
      call Set_OBS2 ( nobs2,
     .                lat2, lon2, level2,
     .                jday1,   stime1,
     .                caldate2, caltime2,
     .                kx2, kt2, xm2,
     .                obs2, sigo2, qc_excl2, qc_hist2 )

      call ODS_Cal2Time ( id1, nobs2, caldate2, caltime2,
     .                                time2, ierr )

      call ODS_Append ( id1, nobs2, ierr )

      call ODS_PutR ( id1,  'lat',     jday1, stime1,
     .                nobs2, lat2,     ierr )
      call ODS_PutR ( id1,  'lon',     jday1, stime1,
     .                nobs2, lon2,     ierr )
      call ODS_PutR ( id1,  'lev',     jday1, stime1,
     .                nobs2, level2,   ierr )
      call ODS_PutI ( id1,  'time',    jday1, stime1,
     .                nobs2, time2,    ierr )
      call ODS_PutI ( id1,  'kt',      jday1, stime1,
     .                nobs2, kt2,      ierr )
      call ODS_PutI ( id1,  'kx',      jday1, stime1,
     .                nobs2, kx2,      ierr )
      call ODS_PutR ( id1,  'xm',      jday1, stime1,
     .                nobs2, xm2,      ierr )
      call ODS_PutR ( id1,  'obs',     jday1, stime1,
     .                nobs2, obs2,     ierr )
c      call ODS_PutR ( id1,  'sigO',    jday1, stime2,
c     .                nobs2, sigo2,    ierr )
      call ODS_PutI ( id1,  'qcexcl',  jday1, stime1,
     .                nobs2, qc_excl2, ierr )
      call ODS_PutI ( id1,  'qchist',  jday1, stime1,
     .                nobs2, qc_hist2, ierr )

      id2     = id1
      jday2   = first_jday
      stime2  = 18

      call Set_OBS2 ( nobs2,
     .                lat2, lon2, level2,
     .                jday2,   stime2,
     .                caldate2, caltime2,
     .                kx2, kt2, xm2,
     .                obs2, sigo2, qc_excl2, qc_hist2 )

      call ODS_Cal2Time ( id2, nobs2, caldate2, caltime2,
     .                                time2, ierr )

      call ODS_PutR ( id2,  'lat',     jday2, stime2,
     .                nobs2, lat2,     ierr )
      call ODS_PutR ( id2,  'lon',     jday2, stime2,
     .                nobs2, lon2,     ierr )
      call ODS_PutR ( id2,  'lev',     jday2, stime2,
     .                nobs2, level2,   ierr )
      call ODS_PutI ( id2,  'time',    jday2, stime2,
     .                nobs2, time2,    ierr )
      call ODS_PutI ( id2,  'kt',      jday2, stime2,
     .                nobs2, kt2,      ierr )
      call ODS_PutI ( id2,  'kx',      jday2, stime2,
     .                nobs2, kx2,      ierr )
      call ODS_PutR ( id2,  'xm',      jday2, stime2,
     .                nobs2, xm2,      ierr )
      call ODS_PutR ( id2,  'obs',     jday2, stime2,
     .                nobs2, obs2,     ierr )
      call ODS_PutI ( id2,  'qcexcl',  jday2, stime2,
     .                nobs2, qc_excl2, ierr )
      call ODS_PutI ( id2,  'qchist',  jday2, stime2,
     .                nobs2, qc_hist2, ierr )

      call ODS_Append ( id2, nobs1, ierr )

      call ODS_PutR ( id2,  'lat',     jday2, stime2,
     .                nobs1, lat1,     ierr )
      call ODS_PutR ( id2,  'lon',     jday2, stime2,
     .                nobs1, lon1,     ierr )
      call ODS_PutR ( id2,  'lev',     jday2, stime2,
     .                nobs1, level1,   ierr )
      call ODS_PutI ( id2,  'time',    jday2, stime2,
     .                nobs1, time1,    ierr )
      call ODS_PutI ( id2,  'kt',      jday2, stime2,
     .                nobs1, kt1,      ierr )
      call ODS_PutI ( id2,  'kx',      jday2, stime2,
     .                nobs1, kx1,      ierr )
      call ODS_PutR ( id2,  'xm',      jday2, stime2,
     .                nobs1, xm1,      ierr )
      call ODS_PutR ( id2,  'obs',     jday2, stime2,
     .                nobs1, obs1,     ierr )
      call ODS_PutI ( id2,  'qcexcl',  jday2, stime2,
     .                nobs1, qc_excl1, ierr )
      call ODS_PutI ( id2,  'qchist',  jday2, stime2,
     .                nobs1, qc_hist1, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 21 )
 21   format ( /' Obs set 2 written to file, test1.nc ' )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt1,  kt_names1, kt_units1,
     .                   n_kx1,  kx_names1, kx_meta1,
     .                   n_qcx1, qcx_names1,
     .                   nobs2,
     .                   lat2, lon2, level2,
     .                   caldate2, caltime2,
     .                   kx2, kt2, xm2,
     .                   obs2, sigo2, qc_excl2, qc_hist2 )

      id3     = id1
      nobs3   = 10
      jday3   = first_jday + 5
      stime3  = 18

      call Set_OBS3 ( nobs3,
     .                lat3, lon3, level3,
     .                jday3,   stime3,
     .                caldate3, caltime3,
     .                kx3, kt3, xm3,
     .                obs3, sigo3, qc_excl3, qc_hist3 )

      call ODS_Cal2Time ( id3, nobs3, caldate3, caltime3,
     .                                time3, ierr )

      call ODS_PutR ( id3,  'lat',     jday3, stime3,
     .                nobs3, lat3,     ierr )
      call ODS_PutR ( id3,  'lon',     jday3, stime3,
     .                nobs3, lon3,     ierr )
      call ODS_PutR ( id3,  'lev',     jday3, stime3,
     .                nobs3, level3,   ierr )
      call ODS_PutI ( id3,  'time',    jday3, stime3,
     .                nobs3, time3,    ierr )
      call ODS_PutI ( id3,  'kt',      jday3, stime3,
     .                nobs3, kt3,      ierr )
      call ODS_PutI ( id3,  'kx',      jday3, stime3,
     .                nobs3, kx3,      ierr )
      call ODS_PutR ( id3,  'xm',      jday3, stime3,
     .                nobs3, xm3,      ierr )
      call ODS_PutR ( id3,  'obs',     jday3, stime3,
     .                nobs3, obs3,     ierr )
      call ODS_PutI ( id3,  'qcexcl',  jday3, stime3,
     .                nobs3, qc_excl3, ierr )
      call ODS_PutI ( id3,  'qchist',  jday3, stime3,
     .                nobs3, qc_hist3, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 31 )
 31   format ( /' Obs set 3 written to file, test1.nc ' )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt1,  kt_names1, kt_units1,
     .                   n_kx1,  kx_names1, kx_meta1,
     .                   n_qcx1, qcx_names1,
     .                   nobs3,
     .                   lat3, lon3, level3,
     .                   caldate3, caltime3,
     .                   kx3, kt3, xm3,
     .                   obs3, sigo3, qc_excl3, qc_hist3 )


      call ODS_Close ( id3, 'test1', ierr )

      call ODS_Open  ( id2, 'test1.nc', 'r', ierr )

      n_kt2  = Dnkt
      n_kx2  = Dnkx
      n_qcx2 = Dnqcx
      call ODS_Info  ( id2, filename, mode, version,
     .                 first_jday, latest_jday, latest__hour,
     .                 n_kt2,  kt_names2, kt_units2,
     .                 n_kx2,  kx_names2, kx_meta2,
     .                 n_qcx2, qcx_names2,
     .                 ierr )

      nobs4   = mobs
      jday3   = first_jday
      stime4  = 18

      call ODS_GetR ( id2,  'lat',     jday2, stime2,
     .                nobs4, lat4,     ierr )
      call ODS_GetR ( id2,  'lon',     jday2, stime2,
     .                nobs4, lon4,     ierr )
      call ODS_GetR ( id2,  'lev',     jday2, stime2,
     .                nobs4, level4,   ierr )
      call ODS_GetI ( id2,  'time',    jday2, stime2,
     .                nobs4, time4,    ierr )
      call ODS_GetI ( id2,  'kt',      jday2, stime2,
     .                nobs4, kt4,      ierr )
      call ODS_GetI ( id2,  'kx',      jday2, stime2,
     .                nobs4, kx4,      ierr )
      call ODS_GetR ( id2,  'xm',      jday2, stime2,
     .                nobs4, xm4,      ierr )
      call ODS_GetR ( id2,  'obs',     jday2, stime2,
     .                nobs4, obs4,     ierr )
c      call ODS_GetR ( id2,  'sigO',    jday2, stime2,
c     .                nobs4, sigo4,    ierr )
      call ODS_GetI ( id2,  'qcexcl',  jday2, stime2,
     .                nobs4, qc_excl4, ierr )
      call ODS_GetI ( id2,  'qchist',  jday2, stime2,
     .                nobs4, qc_hist4, ierr )

      call ODS_Time2Cal ( id2, nobs4, time4,
     .                                caldate4, caltime4, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 41 )
 41   format ( /' Obs set 2 read form file, test1.nc ' )
      len_str = ODS_StrSize ( version )      
      write ( stdout, 42 ) version ( : len_str )
 42   format ( /'   ODS version: ', a )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt2,  kt_names2, kt_units2,
     .                   n_kx2,  kx_names2, kx_meta2,
     .                   n_qcx2, qcx_names2,
     .                   nobs4,
     .                   lat4, lon4, level4,
     .                   caldate4, caltime4,
     .                   kx4, kt4, xm4,
     .                   obs4, sigo4, qc_excl4, qc_hist4 )
      call ODS_CGet ( id2, ':type', type, ierr )
      len_str = ODS_StrSize ( type )
      write ( stdout, 43 ) type ( : len_str)
 43   format ( /'ODS file type is ', a )

      call ODS_Close ( id2, 'test2', ierr )

      call ODS_Open  ( id2, 'sample1.nc', 'r', ierr )

      n_kt2  = Dnkt
      n_kx2  = Dnkx
      n_qcx2 = Dnqcx
      call ODS_Info  ( id2, filename, mode, version,
     .                 first_jday, latest_jday, latest__hour,
     .                 n_kt2,  kt_names2, kt_units2,
     .                 n_kx2,  kx_names2, kx_meta2,
     .                 n_qcx2, qcx_names2,
     .                 ierr )

      nobs4   = mobs

      call ODS_GetR ( id2,  'lat',     jday2, stime2,
     .                nobs4, lat4,     ierr )
      call ODS_GetR ( id2,  'lon',     jday2, stime2,
     .                nobs4, lon4,     ierr )
      call ODS_GetR ( id2,  'lev',     jday2, stime2,
     .                nobs4, level4,   ierr )
      call ODS_GetI ( id2,  'time',    jday2, stime2,
     .                nobs4, time4,    ierr )
      call ODS_GetI ( id2,  'kt',      jday2, stime2,
     .                nobs4, kt4,      ierr )
      call ODS_GetI ( id2,  'kx',      jday2, stime2,
     .                nobs4, kx4,      ierr )
      call ODS_GetR ( id2,  'xm',      jday2, stime2,
     .                nobs4, xm4,      ierr )
      call ODS_GetR ( id2,  'obs',     jday2, stime2,
     .                nobs4, obs4,     ierr )
c      call ODS_GetR ( id2,  'sigO',    jday2, stime2,
c     .                nobs4, sigo4,    ierr )
      call ODS_GetI ( id2,  'qcexcl',  jday2, stime2,
     .                nobs4, qc_excl4, ierr )
      call ODS_GetI ( id2,  'qchist',  jday2, stime2,
     .                nobs4, qc_hist4, ierr )

      call ODS_Time2Cal ( id2, nobs4, time4,
     .                                caldate4, caltime4, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 51 )
 51   format ( /' Obs set 2 read form file, test1_import.nc ' )
      len_str = ODS_StrSize ( version )      
      write ( stdout, 52 ) version ( : len_str )
 52   format ( /'   ODS version: ', a )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt2,  kt_names2, kt_units2,
     .                   n_kx2,  kx_names2, kx_meta2,
     .                   n_qcx2, qcx_names2,
     .                   nobs4,
     .                   lat4, lon4, level4,
     .                   caldate4, caltime4,
     .                   kx4, kt4, xm4,
     .                   obs4, sigo4, qc_excl4, qc_hist4 )

      call ODS_Close ( id2, 'test2', ierr )

      call ODS_Open  ( id2, 'sample2.nc', 'r', ierr )

      first_jday2  = first_jday
      first_jday   = 0
      n_kt2  = Dnkt
      n_kx2  = Dnkx
      n_qcx2 = Dnqcx
      call ODS_Info  ( id2, filename, mode, version,
     .                 first_jday, latest_jday, latest__hour,
     .                 n_kt2,  kt_names2, kt_units2,
     .                 n_kx2,  kx_names2, kx_meta2,
     .                 n_qcx2, qcx_names2,
     .                 ierr )

      nobs4   = 0
      jday2   = first_jday
      nobs4   = ODS_NGet ( id2, jday2, stime2, ierr )

c      call ODS_GetList ( id2, 'kx_meta', n_kx2, kx_meta2, ierr )
      call ClearList ( n_kx2,  kx_meta2 )
      call ClearList ( n_qcx2, qcx_names2 )
      call ClearRVal ( nobs4, xm4 )
      call ClearIVal ( nobs4, qc_hist4 )

      call ODS_GetR ( id2,  'lat',     jday2, stime2,
     .                nobs4, lat4,     ierr )
      call ODS_GetR ( id2,  'lon',     jday2, stime2,
     .                nobs4, lon4,     ierr )
      call ODS_GetR ( id2,  'level',   jday2, stime2,
     .                nobs4, level4,   ierr )
      call ODS_GetI ( id2,  'julian',  jday2, stime2,
     .                nobs4, julian4,  ierr )
      call ODS_GetI ( id2,  'time',    jday2, stime2,
     .                nobs4, time4,    ierr )
      call ODS_GetI ( id2,  'kt',      jday2, stime2,
     .                nobs4, kt4,      ierr )
      call ODS_GetI ( id2,  'kx',      jday2, stime2,
     .                nobs4, kx4,      ierr )
      call ODS_GetR ( id2,  'obs',     jday2, stime2,
     .                nobs4, obs4,     ierr )
      call ODS_GetI ( id2,  'qc_flag', jday2, stime2,
     .                nobs4, qc_excl4, ierr )

      call Time1toTime2 ( nobs4,   first_jday,
     .                    julian4, time4, time4 )
      julian_offset ( id2 ) = first_jday2 - 1
      call ODS_Time2Cal ( id2, nobs4, time4,
     .                    caldate4, caltime4, ierr )

      write ( stdout, '(2x)')
      write ( stdout, 61 )
 61   format ( /' Obs set 2 read form file, test1_v1_02.nc ' )
       len_str = ODS_StrSize ( version )      
       write ( stdout, 62 ) version ( : len_str )
  62   format ( /'   ODS version: ', a )
      write ( stdout, '(2x)')
      call  Print_OBS  ( n_kt2,  kt_names2, kt_units2,
     .                   n_kx2,  kx_names2, kx_meta2,
     .                   n_qcx2, qcx_names2,
     .                   nobs4,
     .                   lat4, lon4, level4,
     .                   caldate4, caltime4,
     .                   kx4, kt4, xm4,
     .                   obs4, sigo4, qc_excl4, qc_hist4 )

      call ODS_Close ( id2, 'test2', ierr )

      stop
      end

*...................................................................

      subroutine Time1toTime2 ( NVal, FirstJDay,
     .                          JDay, Time1, Time2 )

      implicit NONE

      integer  NVal
      integer  FirstJDay
      integer  JDay  ( NVal )
      integer  Time1 ( NVal )
      integer  Time2 ( NVal )

      integer  iVal

      do iVal = 1, NVal
         Time2 ( iVal ) =   Time1 ( iVal )
     .                  + ( JDay  ( iVal ) - FirstJDay ) * 1440

      end do

      return
      end

*...................................................................
 
      subroutine ClearList ( NList, List )

      implicit NONE
      integer                    NList
      character * ( * ) List   ( NList )

      integer  iList

      do iList = 1, NList
         List ( iList ) = ' '

      end do

      return
      end
 
*...................................................................

      subroutine ClearIVal ( NVal, Val )

      implicit   NONE
      integer    NVal
      integer    Val ( NVal )

      integer    iVal

      do iVal = 1, NVal
         Val ( iVal ) = 0

      end do

      return
      end
 
*...................................................................

      subroutine ClearRVal ( NVal, Val )

      implicit   NONE
      integer    NVal
      real       Val ( NVal )

      integer    iVal

      do iVal = 1, NVal
         Val ( iVal ) = 0.0

      end do

      return
      end
 
*...................................................................


      subroutine Print_OBS  ( nkt, ktnames, ktunits,
     .                        nkx, kxnames, kxmeta,
     .                        nqcx, qcxnames,
     .                        nobs,
     .                        lat, lon, level,
     .                        caldate, caltime,
     .                        kx, kt, xm,
     .                        obs, sigo, qc_excl, qc_hist )


      implicit NONE

      integer                      nkt
      character * ( * ) ktnames  ( nkt )
      character * ( * ) ktunits  ( nkt )
      integer                      nkx
      character * ( * ) kxnames  ( nkx )
      character * ( * ) kxmeta   ( nkx )
      integer                      nqcx
      character * ( * ) qcxnames ( nqcx )
      integer                      nobs
      real              lat      ( nobs )
      real              lon      ( nobs )
      real              level    ( nobs )
      integer           caldate  ( nobs )
      integer           caltime  ( nobs )
      integer           kx       ( nobs )
      integer           kt       ( nobs )
      real              xm       ( nobs )
      real              obs      ( nobs )
      real              sigo     ( nobs )
      integer           qc_excl  ( nobs )
      integer           qc_hist  ( nobs )

      integer                      ikt
      integer                      ikx
      integer                      iqcx
      integer                      iob

      include 'ods_stdio.h'

      write ( stdout, 1 )
 1    format ( 2x )

      write ( stdout,  1 )
      write ( stdout, 11 )
      do 10, ikt = 1, nkt
         write ( stdout , 12 ) ikt, ktnames ( ikt ) ( 1 : 40 ),
     .                              ktunits ( ikt ) ( 1 : 10 )
 10   continue
 11   format ( 1x, 'ikt',
     .         2x, 'kt names', 32x,
     .         1x, 'kt units' )
 12   format ( 1x, i3,
     .         2x, a,
     .         1x, a )

      write ( stdout,  1 )
      write ( stdout, 21 )
      do 20, ikx = 1, nkx
         if ( kxnames ( ikx ) ( 1:1 ) .ne. ' ' )
     .   write ( stdout , 22 ) ikx, kxnames ( ikx ) ( 1 : 40 ),
     .                              kxmeta  ( ikx ) ( 1 : 40 )

 20   continue
 21   format ( 1x, 'ikx',
     .         2x, 'kx names', 31x,
     .         2x, 'kx meta' )
 22   format ( 1x, i3,
     .         2x, a,
     .         1x, a )

      write ( stdout,  1 )
      write ( stdout, 31 )
      do 30, iqcx = 1, nqcx
         if ( qcxnames ( iqcx ) ( 1:1 ) .ne. ' ' )
     .   write ( stdout , 32 ) iqcx, qcxnames ( iqcx ) ( 1 : 50 )

 30   continue
 31   format ( 1x, 'ikx',
     .         2x, 'qcx names' )
 32   format ( 1x, i3,
     .         2x, a )

      write ( stdout,  1 )
      write ( stdout, 41 )
      do 40, iob = 1, nobs
         write ( stdout, 42 )
     .                     iob,
     .           lat     ( iob ),
     .           lon     ( iob ),
     .           level   ( iob ),
     .           caldate ( iob ),
     .           caltime ( iob ),
     .           kx      ( iob ),
     .           kt      ( iob ),
     .           xm      ( iob ),
     .           obs     ( iob ),
     .           qc_excl ( iob ),
     .           qc_hist ( iob )
  
 40   continue
 41   format ( 1x,      'iob',
     .         1x, '     lat',
     .         1x, '     lon',
     .         1x, '   level',
     .         1x, '    date',
     .         1x,   '  time',
     .         1x,       'kx',
     .         1x,      ' kt',
     .         1x, '      xm',
     .         1x, '     obs',
     .         1x,  'QC excl',
     .         1x,  'QC hist')

 42   format ( 1x, i3,
     .         1x, f8.2,
     .         1x, f8.2,
     .         1x, f8.2,
     .         1x, i8,
     .         1x, i6,
     .         1x, i2,
     .         1x, i3,
     .         1x, f8.2,
     .         1x, f8.2,
     .         1x, i7,
     .         1x, i7 )
      return
      end

*...................................................................

      subroutine SET_K1 ( nkt,  ktnames, ktunits,
     .                    nkx,  kxnames, kxmeta,
     .                    nqcx, qcxnames )

      implicit NONE

      integer                      nkt
      character * ( * ) ktnames  ( nkt )
      character * ( * ) ktunits  ( nkt )
      integer                      nkx
      character * ( * ) kxnames  ( nkx )
      character * ( * ) kxmeta   ( nkx )
      integer                      nqcx
      character * ( * ) qcxnames ( nqcx )

      integer                      ikt
      integer                      ikx
      integer                      iqcx

      do 10, ikt = 1, nkt
         ktnames ( ikt ) = ' '
         ktunits ( ikt ) = ' '
 10   continue

      do 20, ikx = 1, nkx
         kxnames ( ikx ) = ' '
         kxmeta  ( ikx ) = ' '
 20   continue

      do 30, iqcx = 1, nqcx
         qcxnames ( iqcx ) = ' '
 30   continue

      ktnames   (1) = 'surface pressure $'
      ktunits   (1) = 'HPa $'
      ktnames   (2) = 'surface u-wind component $'
      ktunits   (2) = 'm/s $'
      ktnames   (3) = 'surface v-wind component $'
      ktunits   (3) = 'm/s $'
      ktnames   (4) = 'geopotential height $'
      ktunits   (4) = 'm $'
      ktnames   (5) = 'humidity (dew point) $'
      ktunits   (5) = 'K $'
      ktnames   (6) = 'u-wind component $'
      ktunits   (6) = 'm/s $'
      ktnames   (7) = 'v-wind component $'
      ktunits   (7) = 'm/s $'

      kxnames   (7) = 'rawindsonde $'
      kxmeta    (7) = 'oms_file:rawind.oms $'
      kxnames  (87) = 'tovC $'
      kxmeta   (87) = 'oms_file:tovC.oms $'

      qcxnames ( 1) = 'valid $'
      qcxnames ( 2) = 'warning $'
      qcxnames ( 3) = 'rejected by gross check $'
      qcxnames ( 4) = 'rejected by buddy check $'

      return
      end

*...................................................................


      subroutine SET_OBS  ( nobs,
     .                      lat, lon, level,
     .                      jday,   stime,
     .                      date, time,
     .                      kx, kt, xm,
     .                      obs, sigo, qc_excl, qc_hist )

      implicit NONE

      integer                     nobs
      real              lat     ( nobs )
      real              lon     ( nobs )
      real              level   ( nobs )
      integer           jday
      integer           stime
      integer           date    ( nobs )
      integer           time    ( nobs )
      integer           kx      ( nobs )
      integer           kt      ( nobs )
      real              xm      ( nobs )
      real              obs     ( nobs )
      real              sigo    ( nobs )
      integer           qc_excl ( nobs )
      integer           qc_hist ( nobs )

      integer                     iob
      integer           ODS_CalDat
      integer           sdate

      sdate = ODS_CalDat ( jday )
      do 10, iob = 1, nobs
         lat     ( iob )  = 0
         lon     ( iob )  = 0
         level   ( iob )  = 0
         date    ( iob )  = sdate
         time    ( iob )  = stime * 10000
         kx      ( iob )  = 1
         kt      ( iob )  = 1
         xm      ( iob )  = 0.0
         obs     ( iob )  = 0.0
         sigo    ( iob )  = 0.0
         qc_excl ( iob )  = 0
         qc_hist ( iob )  = 0
 10   continue

      return
      end

*...................................................................


      subroutine SET_OBS3  ( nobs,
     .                       lat, lon, level,
     .                       jday,   stime,
     .                       date, time,
     .                       kx, kt, xm,
     .                       obs, sigo, qc_excl, qc_hist )

      implicit NONE

      integer                     nobs
      real              lat     ( nobs )
      real              lon     ( nobs )
      real              level   ( nobs )
      integer           jday
      integer           stime
      integer           date    ( nobs )
      integer           time    ( nobs )
      integer           kx      ( nobs )
      integer           kt      ( nobs )
      real              xm      ( nobs )
      real              obs     ( nobs )
      real              sigo    ( nobs )
      integer           qc_excl ( nobs )
      integer           qc_hist ( nobs )

      integer                     iob
      integer           ODS_CalDat
      integer           sdate

      sdate = ODS_CalDat ( jday )
      do 10, iob = 1, nobs
         lat     ( iob )  = 0
         lon     ( iob )  = 0
         level   ( iob )  = 0
         date    ( iob )  = sdate
         time    ( iob )  = stime * 10000
         kx      ( iob )  = 1
         kt      ( iob )  = 1
         xm      ( iob )  = 0.0
         obs     ( iob )  = 0.0
         sigo    ( iob )  = 0.0
         qc_excl ( iob )  = 0
         qc_hist ( iob )  = 0
 10   continue

      lat     ( 1  )  = 30.12
      lon     ( 1  )  = -80.90
      level   ( 1  )  = 500.0
      kx      ( 1  )  = 7
      kt      ( 1  )  = 4
      xm      ( 1  )  = 1.5
      obs     ( 1  )  = 5263.0
      sigo    ( 1  )  = 15.0
      qc_excl ( 1  )  = 1
      qc_hist ( 1  )  = 11

      lat     ( 2  )  = 30.12
      lon     ( 2  )  = -80.90
      level   ( 2  )  = 500.0
      kx      ( 2  )  = 7
      kt      ( 2  )  = 5
      xm      ( 2  )  = 2.5
      obs     ( 2  )  = 15.0
      sigo    ( 2  )  = 2.0
      qc_excl ( 2  )  = 2
      qc_hist ( 2  )  = 22

      lat     ( 3  )  = 30.12
      lon     ( 3  )  = -80.90
      level   ( 3  )  = 500.0
      kx      ( 3  )  = 7
      kt      ( 3  )  = 6
      xm      ( 3  )  = 3.5
      obs     ( 3  )  = 11.0
      sigo    ( 3  )  = 2.0
      qc_excl ( 3  )  = 2
      qc_hist ( 3  )  = 33

      lat     ( 4  )  = 30.12
      lon     ( 4  )  = -80.90
      level   ( 4  )  = 300.0
      kx      ( 4  )  = 7
      kt      ( 4  )  = 4
      xm      ( 4  )  = 4.5
      obs     ( 4  )  = 9376.0
      sigo    ( 4  )  = 25.0
      qc_excl ( 4  )  = 5
      qc_hist ( 4  )  = 44

      lat     ( 6  )  = 30.12
      lon     ( 6  )  = -80.90
      level   ( 6  )  = 300.0
      kx      ( 6  )  = 7
      kt      ( 6  )  = 5
      xm      ( 6  )  = 6.5
      obs     ( 6  )  = 17.0
      sigo    ( 6  )  = 2.0
      qc_excl ( 6  )  = 7
      qc_hist ( 6  )  = 66

      lat     ( 7  )  = 30.12
      lon     ( 7  )  = -80.90
      level   ( 7  )  = 300.0
      kx      ( 7  )  = 7
      kt      ( 7  )  = 6
      xm      ( 7  )  = 7.5
      obs     ( 7  )  = 22.0
      sigo    ( 7  )  = 2.0
      qc_excl ( 7  )  = 2
      qc_hist ( 7  )  = 77

      lat     ( 10 )  = 45.82
      lon     ( 10 )  = 100.90
      level   ( 10 )  = 850.0
      kx      ( 10 )  = 7
      kt      ( 10 )  = 4
      xm      ( 10 )  = 10.5
      obs     ( 10 )  = 1367.0
      sigo    ( 10 )  = 2.0
      qc_excl ( 10 )  = 0
      qc_hist ( 10 )  = 111

      return
      end

*...................................................................


      subroutine SET_OBS2 ( nobs,
     .                      lat, lon, level,
     .                      jday,   stime,
     .                      date  , time,
     .                      kx, kt, xm,
     .                      obs, sigo, qc_excl, qc_hist )

      implicit NONE

      integer                     nobs
      real              lat     ( nobs )
      real              lon     ( nobs )
      real              level   ( nobs )
      integer           jday
      integer           stime
      integer           date    ( nobs )
      integer           time    ( nobs )
      integer           kx      ( nobs )
      integer           kt      ( nobs )
      real              xm      ( nobs )
      real              obs     ( nobs )
      real              sigo    ( nobs )
      integer           qc_excl ( nobs )
      integer           qc_hist ( nobs )

      integer                     iob
      integer           ODS_CalDat
      integer           sdate

      sdate = ODS_CalDat ( jday )
      do 10, iob = 1, nobs
         lat     ( iob )  = 0
         lon     ( iob )  = 0
         level   ( iob )  = 0
         date    ( iob )  = sdate
         time    ( iob )  = stime * 10000
         kx      ( iob )  = 1
         kt      ( iob )  = 1
         xm      ( iob )  = 0.0
         obs     ( iob )  = 0.0
         sigo    ( iob )  = 0.0
         qc_excl ( iob )  = 0
         qc_hist ( iob )  = 0
 10   continue

      lat     ( 1  )  = 30.12
      lon     ( 1  )  = -80.90
      level   ( 1  )  = 500.0
      kx      ( 1  )  = 7
      kt      ( 1  )  = 4 
      xm      ( 1  )  = 2.01
      obs     ( 1  )  = 5371.0
      sigo    ( 1  )  = 15.0
      qc_excl ( 1  )  = 0
      qc_hist ( 1  )  = 2

      lat     ( 2  )  = 30.12
      lon     ( 2  )  = -80.90
      level   ( 2  )  = 500.0
      kx      ( 2  )  = 7
      kt      ( 2  )  = 5
      xm      ( 2  )  = 1.2
      obs     ( 2  )  = 12.0
      sigo    ( 2  )  = 2.0
      qc_excl ( 2  )  = 1
      qc_hist ( 2  )  = 7

      lat     ( 3  )  = 30.12
      lon     ( 3  )  = -80.90
      level   ( 3  )  = 500.0
      kx      ( 3  )  = 7
      kt      ( 3  )  = 6
      xm      ( 3  )  = 4.0
      obs     ( 3  )  = 18.0
      sigo    ( 3  )  = 2.0
      qc_excl ( 3  )  = 0
      qc_hist ( 3  )  = 2

      lat     ( 4  )  = 30.12
      lon     ( 4  )  = -80.90
      level   ( 4  )  = 300.0
      kx      ( 4  )  = 7
      kt      ( 4  )  = 4
      xm      ( 4  )  = 7.1
      obs     ( 4  )  = 9485.0
      sigo    ( 4  )  = 25.0
      qc_excl ( 4  )  = 1
      qc_hist ( 4  )  = 23

      lat     ( 6  )  = 30.12
      lon     ( 6  )  = -80.90
      level   ( 6  )  = 300.0
      kx      ( 6  )  = 7
      kt      ( 6  )  = 5
      xm      ( 6  )  = 15.1
      obs     ( 6  )  = 28.0
      sigo    ( 6  )  = 2.0
      qc_excl ( 6  )  = 2
      qc_hist ( 6  )  = 2

      lat     ( 7  )  = 30.12
      lon     ( 7  )  = -80.90
      level   ( 7  )  = 300.0
      kx      ( 7  )  = 7
      kt      ( 7  )  = 6
      xm      ( 7  )  = 21.01
      obs     ( 7  )  = 6.0
      sigo    ( 7  )  = 2.0
      qc_excl ( 7  )  = 0
      qc_hist ( 7  )  = 7

      lat     ( 10 )  = 45.82
      lon     ( 10 )  = 100.90
      level   ( 10 )  = 850.0
      kx      ( 10 )  = 7
      kt      ( 10 )  = 4
      xm      ( 10 )  = 14.0
      obs     ( 10 )  = 1492.0
      sigo    ( 10 )  = 2.0
      qc_excl ( 10 )  = 1
      qc_hist ( 10 )  = 2

      lat     ( 13 )  = 45.82
      lon     ( 13 )  = 100.90
      level   ( 13 )  = 850.0
      kx      ( 13 )  = 7
      kt      ( 13 )  = 5
      xm      ( 13 )  = 66.02
      obs     ( 13 )  = 281.0
      sigo    ( 13 )  = 1.0
      qc_excl ( 13 )  = 1
      qc_hist ( 13 )  = 15

      lat     ( 15 )  = 45.82
      lon     ( 15 )  = 100.90
      level   ( 15 )  = 500.0
      kx      ( 15 )  = 7
      kt      ( 15 )  = 6
      xm      ( 15 )  = 6.2
      obs     ( 15 )  = 21.0
      sigo    ( 15 )  = 2.0
      qc_excl ( 15 )  = 1
      qc_hist ( 15 )  = 22

      lat     ( 16 )  = 45.82
      lon     ( 16 )  = 100.90
      level   ( 16 )  = 500.0
      kx      ( 16 )  = 7
      kt      ( 16 )  = 7
      xm      ( 16 )  = 1.1
      obs     ( 16 )  = 18.0
      sigo    ( 16 )  = 2.0
      qc_excl ( 16 )  = 3
      qc_hist ( 16 )  = 18

      lat     ( 17 )  = 45.82
      lon     ( 17 )  = 100.90
      level   ( 17 )  = 500.0
      kx      ( 17 )  = 7
      kt      ( 17 )  = 5
      xm      ( 17 )  = 13.7
      obs     ( 17 )  = 274.0
      sigo    ( 17 )  = 2.0
      qc_excl ( 17 )  = 1
      qc_hist ( 17 )  = 75

      lat     ( 20 )  = 45.82
      lon     ( 20 )  = 100.90
      level   ( 20 )  = 300.0
      kx      ( 20 )  = 7
      kt      ( 20 )  = 4
      xm      ( 20 )  = 1.0
      obs     ( 20 )  = 9355.0
      sigo    ( 20 )  = 25.0
      qc_excl ( 20 )  = 2
      qc_hist ( 20 )  = 15

      lat     ( 21 )  = 10.51
      lon     ( 21 )  = -71.44
      level   ( 21 )  = 850.0
      kx      ( 21 )  = 7
      kt      ( 21 )  = 4
      xm      ( 21 )  = 4.0
      obs     ( 21 )  = 1655.0
      sigo    ( 21 )  = 7.0
      qc_excl ( 21 )  = 1
      qc_hist ( 21 )  = 19

      lat     ( 23 )  = 10.51
      lon     ( 23 )  = -70.44
      level   ( 23 )  = 700.0
      kx      ( 23 )  = 7
      kt      ( 23 )  = 4
      xm      ( 23 )  = -6.21
      obs     ( 23 )  = 2985.0
      sigo    ( 23 )  = 8.0
      qc_excl ( 23 )  = 4
      qc_hist ( 23 )  = 22

      lat     ( 25 )  = 10.51
      lon     ( 25 )  = -70.44
      level   ( 25 )  = 500.0
      kx      ( 25 )  = 7
      kt      ( 25 )  = 4
      xm      ( 25 )  = 5.11
      obs     ( 25 )  = 5850.0
      sigo    ( 25 )  = 10.0
      qc_excl ( 25 )  = 1
      qc_hist ( 25 )  = 44

 
      return
      end

*...................................................................


      subroutine SET_OBS1 ( nobs,
     .                      lat, lon, level,
     .                      jday,   stime,
     .                      date  , time,
     .                      kx, kt, xm,
     .                      obs, sigo, qc_excl, qc_hist )

      implicit NONE

      integer                     nobs
      real              lat     ( nobs )
      real              lon     ( nobs )
      real              level   ( nobs )
      integer           jday
      integer           stime
      integer           date    ( nobs )
      integer           time    ( nobs )
      integer           kx      ( nobs )
      integer           kt      ( nobs )
      real              xm      ( nobs )
      real              obs     ( nobs )
      real              sigo    ( nobs )
      integer           qc_excl ( nobs )
      integer           qc_hist ( nobs )

      integer                     iob
      integer           ODS_CalDat
      integer           sdate

      sdate = ODS_CalDat ( jday )
      do 10, iob = 1, nobs
         lat     ( iob )  = 0
         lon     ( iob )  = 0
         level   ( iob )  = 0
         date    ( iob )  = sdate
         time    ( iob )  = stime * 10000
         kx      ( iob )  = 1
         kt      ( iob )  = 1
         xm      ( iob )  = 0.0
         obs     ( iob )  = 0.0
         sigo    ( iob )  = 0.0
         qc_excl ( iob )  = 0
         qc_hist ( iob )  = 0
 10   continue

      lat     ( 1  )  = 30.12
      lon     ( 1  )  = -80.90
      level   ( 1  )  = 500.0
      kx      ( 1  )  = 7
      kt      ( 1  )  = 4
      xm      ( 1  )  = 1.1
      obs     ( 1  )  = 5432.0
      sigo    ( 1  )  = 15.0
      qc_excl ( 1  )  = 2
      qc_hist ( 1  )  = 5

      lat     ( 4  )  = 30.12
      lon     ( 4  )  = -80.90
      level   ( 4  )  = 300.0
      kx      ( 4  )  = 7
      kt      ( 4  )  = 4
      xm      ( 4  )  = 2.2
      obs     ( 4  )  = 9522.0
      sigo    ( 4  )  = 25.0
      qc_excl ( 4  )  = 3
      qc_hist ( 4  )  = 2

      lat     ( 6  )  = 30.12
      lon     ( 6  )  = -80.90
      level   ( 6  )  = 300.0
      kx      ( 6  )  = 7
      kt      ( 6  )  = 5
      xm      ( 6  )  = 5.0
      obs     ( 6  )  = 22.0
      sigo    ( 6  )  = 2.0
      qc_excl ( 6  )  = 1
      qc_hist ( 6  )  = 17

      lat     ( 7  )  = 30.12
      lon     ( 7  )  = -80.90
      level   ( 7  )  = 300.0
      kx      ( 7  )  = 7
      kt      ( 7  )  = 5
      xm      ( 7  )  = 22.01
      obs     ( 7  )  = 15.0
      sigo    ( 7  )  = 2.0
      qc_excl ( 7  )  = 4
      qc_hist ( 7  )  = 1

      lat     ( 10 )  = 45.82
      lon     ( 10 )  = 100.90
      level   ( 10 )  = 850.0
      kx      ( 10 )  = 7
      kt      ( 10 )  = 4
      xm      ( 10 )  = -2.1
      obs     ( 10 )  = 1524.0
      sigo    ( 10 )  = 2.0
      qc_excl ( 10 )  = 1
      qc_hist ( 10 )  = 200

      lat     ( 15 )  = 45.82
      lon     ( 15 )  = 100.90
      level   ( 15 )  = 500.0
      kx      ( 15 )  = 7
      kt      ( 15 )  = 6
      xm      ( 15 )  = 33.0
      obs     ( 15 )  = 17.0
      sigo    ( 15 )  = 2.0
      qc_excl ( 15 )  = 2
      qc_hist ( 15 )  = 16

      lat     ( 16 )  = 45.82
      lon     ( 16 )  = 100.90
      level   ( 16 )  = 500.0
      kx      ( 16 )  = 7
      kt      ( 16 )  = 7
      xm      ( 16 )  = 21.0
      obs     ( 16 )  = 13.0
      sigo    ( 16 )  = 2.0
      qc_excl ( 16 )  = 1
      qc_hist ( 16 )  = 4

      lat     ( 17 )  = 45.82
      lon     ( 17 )  = 100.90
      level   ( 17 )  = 500.0
      kx      ( 17 )  = 7
      kt      ( 17 )  = 5
      xm      ( 17 )  = 1.05
      obs     ( 17 )  = 268.0
      sigo    ( 17 )  = 2.0
      qc_excl ( 17 )  = 4
      qc_hist ( 17 )  = 2

      lat     ( 20 )  = 45.82
      lon     ( 20 )  = 100.90
      level   ( 20 )  = 300.0
      kx      ( 20 )  = 7
      kt      ( 20 )  = 4
      xm      ( 20 )  = 0.1
      obs     ( 20 )  = 9522.0
      sigo    ( 20 )  = 25.0
      qc_excl ( 20 )  = 1
      qc_hist ( 20 )  = 8

      return
      end

*....................................................................


      subroutine ODS_Info ( id,        filename,   mode, version,
     .                      FirstJDay, LatestJDay, LatestHour,
     .                      n_kt,      kt_names,   kt_units,
     .                      n_kx,      kx_names,   kx_meta,
     .                      n_qcx,     qcx_names,  ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_Info() --- Returns information about an opened
!                           ODS file
! 
! !DESCRIPTION:
!     Returns the NetCDF and file status parameters
!
! !INTERFACE: 
!     call ODS_Info  ( id,        filename,   mode, version,
!                      FirstJDay, LatestJDay, LatestHour
!                      n_kt,      kt_names,   kt_units,
!                      n_kx,      kx_names,   kx_meta,
!                      n_qcx,     qcx_names,
!                      ierr )
!
!
! !INPUT PARAMETERS:
      integer          id             ! ODS file handle (file id)
                                      !   use this to refer to this
                                      !   file later on.

! !INPUT/OUTPUT PARAMETERS:
      integer          n_kt           ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of GEOS/DAS data
     .                                !   types
      integer          n_kx           ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of GEOS/DAS data
     .                                !   sources
      integer          n_qcx          ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of possible 
                                      !   values for the quality
                                      !   control exclusion mark
                                      !   (qcexcl)
                                      ! note: If the input values
                                      !   for either nkt, kx or
                                      !   nqcx is smaller than
                                      !   the number of values to
                                      !   be read, then the routine
                                      !   exits with a system error
                                      !   message and a code of
                                      !   ODS_DimErr ( = -3 ).
                                      !   Also, nkt, nkx and/or
                                      !   nqcx are set to the
                                      !   minimum value required
                                      !   in order for the routine
                                      !   to execute successfully.
!
! !OUTPUT PARAMETERS:
      character * (*)  filename       ! name of ODS file
      character * (*)  mode           ! mode = 'w' open for writing
                                      ! mode = 'r' open for reading 
      character * (*)  version        ! version tag
      integer          FirstJDay      ! first julian day on file
      integer          LatestJDay     ! latest julian day for which
                                      !   data exists
      integer          LatestHour     ! latest julian hour for which
                                      !   data exists
      character * (*)  kt_names  ( n_kt )
                                      ! name for each data type
      character * (*)  kt_units  ( n_kt )
                                      ! units for each data type
      character * (*)  kx_names  ( n_kx )
                                      ! name for each data source
      character * (*)  kx_meta   ( n_kx )
                                      ! kx specific metadata information
                                      !   Use this to specify the
                                      !   meaning of the parameter "xm"
                                      !   or to specify the name of the
                                      !   external OMS file name.
                                      !   e.g. "oms_file:myinstr.oms"
      character * (*)  qcx_names ( n_qcx ) 
                                      ! information about the meaning
                                      !   of the variable, "qcexcl".
      integer          ierr           ! Error code. If non-zero, an 
                                      !   error has occurred. For a list
                                      !   of possible values see Table 8
                                      !   of da Silva and Redder (1995).
                                      !   If an error has occurred, then
                                      !   file was closed.
!
!     NOTE: No more than id_max ( as defined in header file ods_hdf.h )
!           files can be opened at any one time.  Use the function
!           ODS_Julian to obtain the input parameters, FirstJDay, and
!           LatestJDay.
!
! !SEE ALSO:
!     ODS_Create()   creates the ODS file, sets the dimensions
!                    and saves the text data for kt_names, kt_units
!                    and kx_names. 
!
!
! !REVISION HISTORY:
!     09Apr1998   Redder    Original code.
!
!-------------------------------------------------------------------------

      include  'ods_hdf.h'
      include  'netcdf.inc'
      include  'ods_stdio.h'

*     Function referenced
*     -------------------
      integer   ODS_VerIndex  ! determines the ODS version
                              !  index number

*     Other variables
*     ---------------
      integer   str_len       ! string length for each name
     .                        !  and unit
      integer   ierr_temp     ! temporary storage for the
                              !  returned error code
      integer   version_index ! version index number
      integer   version_2     ! index number for version 2.00

*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Check to determine if the file handle id is valid
*     -------------------------------------------------
      if ( id            .lt. 1      .or.
     .     id            .gt. id_max .or.
     .     IOMode ( id ) .eq. CLOSED ) then
         write ( stderr, 901 )
         ierr = NCEBadID
         return
      end if

*     Get ODS file name
*     -----------------
      call ODS_CGet    ( id, ':filename',   filename,    ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get mode (read or write)
*     ------------------------
      call ODS_CGet    ( id, ':mode',       mode,        ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get ODS version number
*     ----------------------
      call ODS_CGet    ( id, ':version',    version,     ierr )
      if ( ierr .ne. NCNoErr ) return
      version_index = ODS_VerIndex ( version, ierr )

*     Get some NetCDF file dimensions and parameters from ODS
*     -------------------------------------------------------
      version_2     = ODS_VerIndex ( '2.00', ierr )
      if ( ierr .ne. NCNoErr ) return

      if ( version_index .ge. version_2 ) then
         call ODS_IGet    ( id, 'nkt',      n_kt,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'nkx',      n_kx,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'nqcx',     n_qcx,       ierr )
         if ( ierr .ne. NCNoErr ) return

      else
         call ODS_IGet    ( id, 'ktmax',    n_kt,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'kxmax',    n_kx,        ierr )
         if ( ierr .ne. NCNoErr ) return

      end if

      call ODS_IGet    ( id, 'strlen',      str_len,     ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:first_julian_day',  
     .                                      FirstJDay,   ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:latest_julian_day',
     .                                      LatestJDay,  ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:latest_synoptic_hour',
     .                                      LatestHour,  ierr )
      if ( ierr .ne. NCNoErr ) return

*     Read the text labels for data types and sources currently
*     available and quality control codes currently implemented
*     ---------------------------------------------------------
      call ODS_GetList ( id, 'kt_names', n_kt, kt_names, ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_GetList ( id, 'kt_units', n_kt, kt_units, ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_GetList ( id, 'kx_names', n_kx, kx_names, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Read these tables only if the file corre-
*     sponds to one of the later versions
*     -----------------------------------------
      if ( version_index .ge. version_2 ) then
         call ODS_GetList
     .               ( id, 'kx_meta',   n_kx,  kx_meta,   ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_GetList
     .               ( id, 'qcx_names', n_qcx, qcx_names, ierr )
         if ( ierr .ne. NCNoErr ) return

      end if

      return
*     ------

 901  format ( /, ' ODS_Info: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )

      end
*....................................................................
