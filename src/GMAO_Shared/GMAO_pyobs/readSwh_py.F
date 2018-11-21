              subroutine swath_header(file_name,n_obs,n_lvl,n_vrsn)
              character(len=*),intent(in) :: file_name
              integer,intent(out) :: n_obs, n_lvl,n_vrsn


              call readswh_header(file_name,n_obs,n_lvl,n_vrsn)

              return
              end
               
!
!
       subroutine swath_data(file_name, num_obs, ndate, hhmmss,
     *                       ob_lon, ob_lat, ob_typ, ob_swh, ob_wnd,
     *                       ob_qc)
!      subroutine rd_swh (UNIT, n_obs, vrsn)
!
!.............................START PROLOGUE............................
!
! MODULE NAME:  swath_data
!
! DESCRIPTION:  reads the altimeter SWH ocean obs files and produces a
!               report whether you want one or not
!
! NOTES:        the qc probablity of error includes flags indicating
!               swh retrieval in ice covered seas and/or shallow
!               water.  a value of 510 is added to the underlying
!               probability to indicate if the ice concentration
!               exceeds 33% and a value of 512 is added if the
!               bottom depth is less than 5 meters.  note that a
!               composite flag of 510 + 512 can also occur.
!
! PARAMETERS:
!      Name          Type        Usage            Description
!   ----------     ---------     ------    ---------------------------
!   n_obs           integer      input     number swh obs
!   unit            integer      input     FORTRAN unit number
!   vrsn            integer      input     version number data file
!
! SWH VARIABLES:
!     ame       Type                     Description
!   --------  --------    ----------------------------------------------
!   ob_age     real       age of the observation in hours since
!                         January 1, 1992.  provides a continuous
!                         time variable.  reported to the nearest
!                         minute.
!   ob_clm     real       SWH climate (not available)
!   ob_dtg     character  SWH retrieval date time group in the form
!                         year, month, day, hour, minute, second
!                         (YYYYMMDDHHMMSS)
!   ob_glb     real       SWH global FNMOC analysis estimate at the
!                         obs location and sampling time
!   ob_lat     real       SWH retrieval latitude (south negative)
!   ob_lon     real       SWH retrieval longitude (west negative)
!   ob_qc      real       SWH retrieval probability of a gross error
!                         (assumes normal pdf of SWH retrieval errors)
!   ob_rcpt    character  SWH retrieval FNMOC receipt date time group
!                         in the form year, month, day, hour, minute,
!                         second (YYYYMMDDHHMMSS)
!   ob_rgn     real       SWH regional FNMOC analysis estimate at the
!                         obs location and sampling time
!   ob_swh     real       SWH retrieval (m)
!   ob_typ     integer    satellite ID (ERS2, Topex, Jason, GFO,
!                         ENVISAT, Topex Interleaved); see ocn_types.h
!                         for codes
!   ob_xvl     real       cross validation SWH value from QC
!   ob_wnd     real       altimeter colocated wind retrieval (m/s)
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      include 'ocn_types.h'
      character(len=*), intent(in) :: file_name
!
!     ..local array dimension
!
      integer              :: n_obs
      integer, intent(in)  :: num_obs
      integer, intent(out) :: ndate(num_obs)
      integer, intent(out) :: hhmmss(num_obs)
      character * 256   ::  err_msg
!
      integer           :: i, k,n
      real              :: ob_age (num_obs)
      real              :: ob_clm (num_obs)
      character(len=14) :: ob_dtg(num_obs),date
      real              :: ob_glb (num_obs)
      real, intent(out) :: ob_lat (num_obs)
      real, intent(out) :: ob_lon (num_obs)
      real, intent(out) :: ob_qc (num_obs)
      character(len=14) :: ob_rcp(num_obs)
      real              :: ob_rgn (num_obs)
      real, intent(out)              ::  ob_swh (num_obs)
      integer, intent(out)           ::  ob_typ (num_obs)
      real              ::  ob_xvl (num_obs)
      real, intent(out)              ::  ob_wnd (num_obs)
      integer           ::  vrsn,len,n_lvl,nymd,n_vrsn
      logical           ::  exist

!
!          Read header
!
          call readswh_header(file_name,n_obs,n_lvl,n_vrsn)
!
!          Read data
!
          call readswh_data(file_name,num_obs,ndate,hhmmss,ob_lon,
     *                    ob_lat,ob_typ,ob_rgn,ob_swh,ob_wnd,
     *                    ob_glb,ob_xvl,ob_qc)
!

      return
      end
