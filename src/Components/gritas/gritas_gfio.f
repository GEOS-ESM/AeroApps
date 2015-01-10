!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_Output --- Write del grids to IEEE GrADS file
! 
! !INTERFACE:
!
       subroutine GFIO_Output ( grid_obasen, 
     &                          var_list, var_descr, listsz, timeinc,
     &                          l2d, list_2d, bias_2d, stdv_2d, nobs_2d,
     &                          l3d, list_3d, bias_3d, stdv_3d, nobs_3d,
     &                          nstat, yyyymmdd, hhmmss, fid )
!
! !USES:
!
      use m_gritas_grids, only : im, jm, km
c      use m_gritas_grids, only : glon
c      use m_gritas_grids, only : glat
c      use m_gritas_grids, only : AMISS
c      use m_gritas_grids, only : nlevs
c      use m_gritas_grids, only : plevs

      implicit NONE
c      include 'gritas.h'       ! lmax
!
! !INPUT PARAMETERS: 
!

      character(len=*)  grid_obasen      ! output file basen name

      integer       listsz               ! size of the list (< lmax)
      character*11  var_list(listsz)     ! variable name for output 
      character*256 var_descr(listsz)    ! ... and its description
      character*11  var_num     ! variable name of number of obs for output 

      integer         l2d          ! actual number of 2d grids
      integer         l3d          ! actual number of 3d grids

      integer         list_2d(l2d) ! list item for 2d grids
      integer         list_3d(l3d) ! list item for 3d grids

                                   ! Surface (2D) grids:
      real  bias_2d(im,jm,l2d)     !   time mean
      real  stdv_2d(im,jm,l2d)     !   stdv
      real  nobs_2d(im,jm,l2d)     !   number of obs per grid box

                                   ! Upper-air (3D) grids:
      real  bias_3d(im,jm,km,l3d)  !   time mean
      real  stdv_3d(im,jm,km,l3d)  !   stdv
      real  nobs_3d(im,jm,km,l3d)  !   number of obs per grid box
      integer nstat                ! No of stats to write to file:
                                   ! nstat = 1   ...  bias
                                   ! nstat = 2   ...  bias, stdv
                                   ! nstat = 3   ...  bias, stdv, nobs
                                   ! nstat < 0   ...  write each stat in
                                   !   separate files each with the base
                                   !   name grid_obasen // '.' // stat
                                   !   where stat = 'bias', 'stdv', or 'nobs' 
      integer        yyyymmdd      ! Year-month-day, e.g., 19971003
      integer        hhmmss        ! Hour-minute-second, e.g., 120000
  
      integer        timeinc        ! output frequency (HHMMSS)
!
! !OUTPUT PARAMETERS: 
!
      integer        fid           ! GFIO file (for closing the file later)
!
! !DESCRIPTION: Write del-grids to HDF file using the GFIO interface.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  17OCt97  G. P. Lou  Modifications. (RT: What modifications?)
!  26Mar98  da Silva   Added nstat parameter for GrISAS's sake.
!  10Jun04  Todling    Bug fix: fid missing from argument list.
!  05Jun09  Redder     Added code to enable routine to write each stat 
!                      in separate files (by setting n < 0).  Fixed bug
!                      by generalizing the algorithm for increment time.
!  20Jul09  Redder     Added the input parameter, var_list
!
!EOP
!-------------------------------------------------------------------------
!BOC
c      integer  i, j, k, l, nymd, nhms


!     GFIO work space
!     ---------------
c      character*80    fname, title, source, contact
c      character*12    vunits(lmax)
c      integer         kmvar(lmax)
c      real            valid_range(2,lmax), packing_range(2,lmax)
c      integer         rc, prec

c      save fname

c      character*255 PrevBasen
c      save          PrevBasen

c      data          PrevBasen / '!@#$%^&**()_+|' /   ! just garbage


      character*255 BaseN
      integer nymd, nhms, nstat2

      nstat2 = abs ( nstat )
      if ( nstat2 .lt. 1 .or. nstat2 .gt. 3 ) then
         call gritas_perror ( 'GFIO_Output: invalid nstat' )
      endif

      BaseN = grid_obasen
      nymd = yyyymmdd
      nhms = hhmmss

      if ( nstat .lt. 0 ) then
         BaseN = trim ( grid_obasen ) // '.bias'
      end if
      call    GFIO_TGroup ( BaseN,
     &                      var_list, var_descr, listsz, timeinc,
     &                      l2d,  list_2d, bias_2d,
     &                      l3d,  list_3d, bias_3d,
     &                      nymd, nhms, fid )

      if ( nstat2 .ge. 2 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.stdv'
         else
            call GFIO_TIncr ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup ( BaseN, 
     &                      var_list, var_descr, listsz, timeinc,
     &                      l2d,  list_2d, stdv_2d,
     &                      l3d,  list_3d, stdv_3d,
     &                      nymd, nhms, fid )
      end if

      if ( nstat2 .ge. 3 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.nobs'
         else
            call GFIO_TIncr ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup ( BaseN, 
     &                      var_list, var_descr, listsz, timeinc,
     &                      l2d,  list_2d, nobs_2d,
     &                      l3d,  list_3d, nobs_3d,
     &                      nymd, nhms, fid )
      end if

!     All done
!     --------
      return
      end

!...............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_TGroup --- Write del grids to IEEE GrADS file for one time group
! 
! !INTERFACE:
!
       subroutine GFIO_TGroup ( grid_obasen, 
     &                          var_list, var_descr, listsz, timeinc,
     &                          l2d, list_2d, data_2d,
     &                          l3d, list_3d, data_3d,
     &                          yyyymmdd, hhmmss, fid )
!
! !USES:
!
      use m_gritas_grids, only : im, jm, km
      use m_gritas_grids, only : glon
      use m_gritas_grids, only : glat
      use m_gritas_grids, only : AMISS
      use m_gritas_grids, only : nlevs
      use m_gritas_grids, only : plevs

      implicit NONE
      include 'gritas.h'       ! lmax
!
! !INPUT PARAMETERS: 
!

      character(len=*)  grid_obasen      ! output file basen name

      integer       listsz               ! size of the list (< lmax)
      character*11  var_list(listsz)     ! variable name for output 
      character*256 var_descr(listsz)    ! ... and its description
      character*11  var_num     ! variable name of number of obs for output 

      integer         l2d          ! actual number of 2d grids
      integer         l3d          ! actual number of 3d grids

      integer         list_2d(l2d) ! list item for 2d grids
      integer         list_3d(l3d) ! list item for 3d grids

                                   ! Surface (2D) grids:
      real  data_2d(im,jm,l2d)     !   time mean
                                   ! Upper-air (3D) grids:
      real  data_3d(im,jm,km,l3d)  !   time mean
      integer        yyyymmdd      ! Year-month-day, e.g., 19971003
      integer        hhmmss        ! Hour-minute-second, e.g., 120000
  
      integer        timeinc       ! output frequency (HHMMSS)
!
! !OUTPUT PARAMETERS: 
!
      integer        fid           ! GFIO file (for closing the file later)
!
! !DESCRIPTION: Write del-grids to HDF file using the GFIO interface.
!
! !REVISION HISTORY: 
!
!  04Jun2009  Redder  Initial code.
!  20Jul2009  Redder  Added the input parameter, var_descr.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      integer  i, j, k, l, nymd, nhms

!     GFIO work space
!     ---------------
      character*80    fname, title, source, contact
      character*12    vunits(lmax)
      integer         kmvar(lmax)
      real            valid_range(2,lmax), packing_range(2,lmax)
      integer         rc, prec

      save fname

      character*255 PrevBasen
      save          PrevBasen

      data          PrevBasen / '!@#$%^&**()_+|' /   ! just garbage

!     New outputfile, created it
!     --------------------------
      if ( trim(PrevBasen) .ne. trim(grid_obasen) ) then            

         PrevBasen = grid_obasen


!        Create the GFIO file
!        --------------------
         fname = trim(grid_obasen)//'.hdf'
         
         print *
         print *, ' [] Writing to new HDF/GFIO file ', trim(fname),
     &            ' on ', yyyymmdd, hhmmss/10000, 'Z'

         title = 'Gridded O-F, O-A or Obs Values'
         source = 'From GEOS/DAS DEL or ODS files'
         contact = 'data@gmao.gsfc.nasa.gov'
         do i = 1, listsz
            vunits(i) = 'none'  ! for now
            valid_range(1,i)   = amiss      
            valid_range(2,i)   = amiss      
            packing_range(1,i) = amiss      
            packing_range(2,i) = amiss      
         end do
         kmvar = -1
         do i = 1, l2d
            l = list_2d(i)
            kmvar(l) = 0
         end do
         do i = 1, l3d
            l = list_3d(i)
            kmvar(l) = nlevs
         end do
         prec = 0               ! 32 bits

         call GFIO_Create ( trim(fname), title, source, contact, amiss,
     &                      im, jm, nlevs, glon, glat, plevs, 'hPa', 
     &                      yyyymmdd, hhmmss, timeinc,
     &                      listsz, var_list, var_descr, vunits, kmvar,
     &                      valid_range, packing_range, prec,
     &                      fid, rc )

         if ( rc .ne. 0 ) call
     & gritas_perror ( 'GFIO_TGroup: could not create HDF/GFIO file.' )

      else

         print *
         print *, ' [] Writing to existing HDF/GFIO file ', 
     &       trim(fname), ' on ', yyyymmdd, hhmmss/10000, 'Z'

      end if


!     Write the data to file
!     ----------------------
      print *

!     Means (uses GrADS t=1 slot)
!     ---------------------------
      nymd = yyyymmdd
      nhms = hhmmss
      do i = 1, l2d
         l = list_2d(i)
         call GFIO_PutVar ( fid, var_list(l), nymd, nhms,
     &                      im, jm, 0, 1, data_2d(1,1,i),  
     &                      rc )
         if ( rc .ne. 0 ) then
            print *, 'rc = ', rc, 'var = ', var_list(l)
            call gritas_perror ( 
     &           'GFIO_TGroup: cannot write 2d bias' )
         end if
      end do
      do i = 1, l3d
         l = list_3d(i)
         do k = 1, nlevs
         call GFIO_PutVar ( fid, var_list(l), nymd, nhms,
     &                      im, jm, k, 1, data_3d(1,1,k,i),  
     &                      rc )
         if ( rc .ne. 0 ) then 
            print *, 'rc = ', rc, 'var = ', var_list(l)
            call gritas_perror ( 
     &          'GFIO_TGroup: cannot write 3d bias' )
         end if
      end do
      end do
      
!     All done
!     --------
      return
      end
!...............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_TIncr --- Increment time
! 
! !INTERFACE:
!
      subroutine GFIO_TIncr ( Date, Time, TInc, Date2, Time2 ) 
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
      integer  Date          ! Date (in YYYYMMDD format )
      integer  Time          ! Time (in   HHMMSS format )
      integer  TInc          ! TInc (in time increment )
!
! !OUTPUT PARAMETERS: 
      integer  Date2         ! Output date
      integer  Time2         ! Output time
!
!...............................................................

      integer JInc, JTime, ODS_Julian, ODS_CalDat
      integer JTime_New, JDay_New, Date_New, Time_New

      JInc  =       ( abs ( TInc ) / 10000 ) *  60   * 60
     .      + ( mod ( abs ( TInc ),  10000 ) / 100 ) * 60
     .      +   mod ( abs ( TInc ),    100 )
      if ( TInc .lt. 0 ) JInc = - JInc  
      JTime =             ( Time   / 10000 ) * 60    * 60 
     .      + ( mod (       Time,    10000 ) / 100 ) * 60
     .      +   mod (       Time,      100 )
      JTime_New = modulo ( JTime + JInc,  24 * 60 * 60 )
      JDay_New  = ODS_Julian ( Date ) + ( JTime + JInc - JTime_New )
     .                                / ( 24 * 60 * 60 )
      
      Date_New  = ODS_CalDat    ( JDay_New )
      Time_New  = 10000 *       ( JTime_New / ( 60 * 60 ))
     .          +   100 * ( mod ( JTime_New,    60 * 60 ) / 60 )
     .          +           mod ( JTime_New,    60 )
      Date2    = Date_New
      Time2    = Time_New 
      return
      end
!...............................................................

!EOC


