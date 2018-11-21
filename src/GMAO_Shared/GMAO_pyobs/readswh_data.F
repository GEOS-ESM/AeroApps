!
!
       subroutine readswh_data(file_name,num_obs,ndate,hhmmss,ob_lon,
     *                       ob_lat,ob_typ,ob_rgn,ob_swh,ob_wnd,
     *                       ob_glb,ob_xvl,ob_qc)
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
      integer,intent(in)       ::  num_obs
      integer                  ::  n_obs
      integer                  ::  ndate(num_obs),hhmmss(num_obs)
!
      integer                  ::  i, k,n
      real                     ::  ob_age (num_obs)
      real                     ::  ob_clm (num_obs)
      character(len=14)        ::  ob_rcp(num_obs)
      character(len=14)        ::  ob_dtg(num_obs)
      character*14             ::  date
      character*256            ::  err_msg

      real,intent(out)         ::  ob_glb (num_obs)
      real,intent(out)         ::  ob_lat (num_obs)
      real,intent(out)         ::  ob_lon (num_obs)
      real,intent(out)         ::  ob_qc (num_obs)
      real,intent(out)         ::  ob_rgn (num_obs)
      real,intent(out)         ::  ob_swh (num_obs)
      integer,intent(out)      ::  ob_typ (num_obs)
      real,intent(out)         ::  ob_xvl (num_obs)
      real,intent(out)         ::  ob_wnd (num_obs)
      integer                  ::  UNIT
      integer                  ::  vrsn,len,n_lvl,nymd,n_vrsn
      logical                      exist
      data UNIT /60/ 
!
!          Read header
!
          len = len_trim (file_name)
          inquire (file=file_name(1:len), exist=exist)
          
          if (exist) then
            open (UNIT, file=file_name(1:len), status='old',  
     *                  form='unformatted')
            read (UNIT) n_obs, n_lvl, n_vrsn
!
!...............................executable..............................
!
!     ..read swh variables
!
            read (UNIT) (ob_glb(i), i = 1, n_obs)
            read (UNIT) (ob_lat(i), i = 1, n_obs)
            read (UNIT) (ob_lon(i), i = 1, n_obs)
            read (UNIT) (ob_age(i), i = 1, n_obs)
            read (UNIT) (ob_clm(i), i = 1, n_obs)
            read (UNIT) (ob_qc(i),  i = 1, n_obs)
            read (UNIT) (ob_typ(i), i = 1, n_obs)
            read (UNIT) (ob_rgn(i), i = 1, n_obs)
            read (UNIT) (ob_swh(i), i = 1, n_obs)
            read (UNIT) (ob_wnd(i), i = 1, n_obs)
            read (UNIT) (ob_xvl(i), i = 1, n_obs)
            read (UNIT) (ob_dtg(i), i = 1, n_obs)
            read (UNIT) (ob_rcp(i), i = 1, n_obs)

            nymd   = 0
            hhmmss = 0

            do n = 1,n_obs
             date = ob_dtg(n)
             read(date,'(i8,i6)') ndate(n),hhmmss(n)
      
!
!          ..produce swh report
!
!             k = 1
!             write (6, '(''  reporting skip factor: '', i10)') k
!             write (6, '(11x,''dtg'', 11x,''rcpt'', 5x,''lat'', 5x,''lon'',   
!    *               4x,''type'', 5x,''swh'', 4x,''wind'', 4x,''clim'',  
!    *               4x,''glbl'', 4x,''regn'', 4x,''xval'', 8x,''qc'')')

!              do i = 1, n_obs, k
!                write (6, '(a,1x,a,2f8.2,i8,6f8.1,2x,f8.3,2x,a)')
!    *                  ob_dtg(i), ob_rcp(i), ob_lat(i), ob_lon(i),  
!    *                   ob_typ(i), ob_swh(i), ob_wnd(i), ob_clm(i),  
!    *                    ob_glb(i), ob_rgn(i), ob_xvl(i), ob_qc(i),   
!    *                    data_lbl(ob_typ(i))
!              enddo
!
            enddo

          else
             write (err_msg, '(''file "'', a, ''" does not exist'')')
     *              file_name(1:len)
             n_obs = -1
          endif
      close (UNIT)
      return
      end
