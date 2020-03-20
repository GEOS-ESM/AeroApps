!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_gritas_grids --- GRITAS dimension and grid parameters
!
! !INTERFACE:

      MODULE m_gritas_grids

! !USES:

      Implicit none

! !DESCRIPTION: 
!
!  Defines dimensions and other relevant parameters.
!
! !DEFINED PARAMETERS: 

!     Grid dimension
!     --------------
      integer, save ::     im                        ! no. zonal      gridpoints
      integer, save ::     jm                        ! no. meridional gridpoints
      integer, save ::     km                        ! max.   no. vertical levels
      integer, save ::     nlevs                     ! actual no. vertical levels

!     Longitudinal grid
!     -----------------
      real, parameter         :: LonMin = -180.      ! = lambda(1)
      real, save              :: dLon                ! mesh size (deg)
      real, allocatable, save :: glon(:)             ! zonal gridpoints

!     real, parameter :: dLon = 360. / im    ! mesh size (deg)
   

!     Latitudinal grid
!     ----------------
      real, parameter         :: LatMin = -90.       ! = phi(1)
      real, save              :: dLat                ! mesh size (deg)
      real, allocatable, save :: glat(:)             ! neridional gridpoints

!     Vertical grid
!     -------------
      real, allocatable, save :: plevs(:)            ! vertical levels in hPa

!     Missing value
!     -------------
      real, parameter :: AMISS = 1.0E15


! !REVISION HISTORY: 
!
!  05sep97   da Silva   Original code.
!  09Jun06   Todling    Converted to module.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_Gritas_Grids'

      interface grids_init  ; module procedure init_        ; end interface
      interface grids_clean ; module procedure clean_       ; end interface
      interface GrITAS_Pint ; module procedure GrITAS_Pint_ ; end interface
      interface GrITAS_Norm ; module procedure GrITAS_Norm_ ; end interface
      interface GrITAS_Accum; module procedure GrITAS_Accum_; end interface

      CONTAINS

      subroutine init_ (im_def,jm_def,km_def,stat)

      implicit none
      integer, intent(in)  :: im_def
      integer, intent(in)  :: jm_def
      integer, intent(in)  :: km_def
      integer, intent(out) :: stat

      integer  i, j, ier

      stat = 0

      im = im_def
      jm = jm_def
      km = km_def

      allocate ( glon(im), glat(jm), plevs(km), stat=ier )
         if(ier/=0) then
            stat = 1
            return
         endif

      dLon = 360. / im
      dLat = 180. / ( jm - 1 )

!     Create horizontal grids for GFIO's sake
!     ---------------------------------------
      do i = 1, im
         glon(i) = LonMin + (i-1) * dLon
      end do
      do j = 1, jm
         glat(j) = LatMin + (j-1) * dLat
      end do

      end subroutine init_

      subroutine clean_(stat)
      implicit none
      integer, intent(out) :: stat
      integer ier
      stat = 0
      deallocate ( glon, glat, plevs, stat=ier )
         if(ier/=0) then
            stat = 1
            return
         endif
      end subroutine clean_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Pint_ --- Set pressure Intervals for gridding
! 
! !INTERFACE:
!
      subroutine GrITAS_Pint_ ( p1, p2 )

! !INPUT PARAMETERS: 
 
! !OUTPUT PARAMETERS:
!
      real       p1(nlevs), p2(nlevs)   ! pressure intervals for gridding
!
! !DESCRIPTION: Determine the pressure level ranges for binning
!               off-mandatory level data. Observations with pressure
!               p1(k) <= p < p2(k) will be gridded at level plevs(k). 
!               The binnings are log-pressure linear. 
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  17Oct97  G. P. Lou  Modification.
!  25Mar98  da Silva   Removed tighning of vertical level interval
!                      introduced by G. P. Lou.
!  24Jun04  Todling    Changed handling of 1-level case.
!
!EOP
!-------------------------------------------------------------------------


      character(len=*), parameter :: myname_ = myname//'GrITAS_Pint_'

      integer i, j, k

!     Check number of levels
!     ----------------------
      if ( nlevs .le. 0 ) call gritas_perror
     &   ( 'Pint: must have at least 1 vertical levels' )

!     Special single level case
!     -------------------------
      if ( nlevs .eq. 1 ) then
           p1(1) = plevs(1)
           p2(1) = plevs(1)
           print *, '       p1         p        p2'
           print 100, p1(1), plevs(1), p2(1)
           print *
           return
      end if


!     Check monotonicity
!     ------------------
      do k = 1, nlevs-1
         if ( plevs(k) .le. plevs(k+1) ) call gritas_perror
     &      ('Pint: pressure levels do not decrease monotonically')
      end do

      if ( plevs(nlevs) .eq. 0. ) call gritas_perror
     &   ( 'Pint: cannot handle 0 hPa' )

      k = 1
      p1(1) = 2000. 
      p2(1) = exp ( 0.5 * ( log(plevs(k+1)) + log(plevs(k)) ) )

      do k = 2, nlevs-1
         p1(k) = exp ( 0.5 * ( log(plevs(k-1)) + log(plevs(k)) ) )
         p2(k) = exp ( 0.5 * ( log(plevs(k+1)) + log(plevs(k)) ) )
      end do

      k = nlevs
      p1(k) = exp ( 0.5 * ( log(plevs(k-1)) + log(plevs(k)) ) )
      p2(k) = 0.
       
!     Echo parameters
!     ---------------
      print *, '       p1         p        p2'
      do k = 1, nlevs
         print 100, p1(k), plevs(k), p2(k)
      end do
 100  format(1x,3F10.1)

      print *

      end subroutine GrITAS_Pint_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GrITAS_Accum_ --- Accumulate O-F on del grids
! 
! !INTERFACE:
!
      subroutine GrITAS_Accum_ ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                           kxkt_table, p1, p2,
     &                           ob_flag,
     &                           l2d, bias_2d, rtms_2d, xcov_2d, nobs_2d,
     &                           l3d, bias_3d, rtms_3d, xcov_3d, nobs_3d )

! !USES:
!

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lat(nobs)     ! latitude (degrees)
      real,    intent(in) ::       lon(nobs)     ! longitude (degrees)
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index

      
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      real,    intent(in) ::           p1(nlevs)
      real,    intent(in) ::           p2(nlevs)     ! pressure intervals for gridding

      integer,    intent(in) ::         l2d          ! actual number of 2d grids
      integer,    intent(in) ::         l3d          ! actual number of 3d grids

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer,    intent(inout) :: ob_flag(nobs)     ! On input, ob is rejected if
                                                     !    ob_flag is other than 0
                                                     ! On output, ob_flag = 0 if ob
                                                     !    was used to accumulate the
                                                     !    bias/RMS
                                                     ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(im,jm,l2d,nr)    !   time means
      real, intent(inout) ::  rtms_2d(im,jm,l2d,nr)    !   root time mean square
      real, intent(inout) ::  xcov_2d(im,jm,l2d)       !   x (co)variance
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                     ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(im,jm,km,l3d,nr) !  time means
      real, intent(inout) ::  rtms_3d(im,jm,km,l3d,nr) !   root time mean square
      real, intent(inout) ::  xcov_3d(im,jm,km,l3d)    !   root time mean square
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box

!
! !DESCRIPTION: Given the O-F for a given synoptic time, this routine
!               accumulates the bias and the RMS at the nearest gridpoint.
!
! !REMARKS:
!
!   1.  xrtms is redundant when nr=1, i.e., xrtms=rtms
!   2.  when nr=2, assumes residuals are from same experiment, that is,
!       no need to distinguish between nobs
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  25Mar98  da Silva   Removed "goto 10" statement; now all must be
!                      assigned a level or else a fatal error will occur
!                      (this was the original design).
!  26Mar98  da Silva   Included periodic boundary conditions in longitude.
!  14Jul99  da SIlva   Added land-surface kt's
!  28oct06  Todling    Generalized to handle cross-terms; note xrtms is 
!                      redundant when nr=1, i.e., xrtms=rtms
!  28May09  Redder     Changed code to define negative value in
!                      kxkt_table as indices for surface variables.
!  07Jul09  Redder     Added parameters ob_flag
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_Accum_'

      integer n, i, j, k, kk, nx, l, ir
      logical surface_var, reject_ob

      if ( nr > 2 ) call gritas_perror ('Accum: could not handle '
     &                               // 'more then 2 residuals')
      nx = max(1,nr)

!     Loop over observations
!     ----------------------
      do 10 n = 1, nobs
         reject_ob = ob_flag ( n ) .ne. 0
         if ( reject_ob  ) go to 10                    ! reject unwanted obs

         l = kxkt_table(kx(n),kt(n))
         surface_var = kxkt_table (kx(n),kt(n)) .lt. 0
         l = abs (l)

         if ( l .le. 0 ) then
            ob_flag ( n ) = 1
            go to 10                                   ! forget it
         end if
         
!        Horizontal gridbox: gridpoint is assumed at center of box
!        ---------------------------------------------------------
         i = 1 + nint ( (lon(n) - LonMin) / dLon )
         if ( i .eq. (im+1) ) i = 1                      ! periodic bc
         if ( i .eq.   0    ) i = im                     ! periodic bc
         j = 1 + nint ( (lat(n) - LatMin) / dLat )

         i = min(im,max(1,i))
         j = min(jm,max(1,j))
        

 
!        Surface (2D) grid
!        -----------------
         if ( surface_var )  then  ! surface variables 

!             Accumulate
!             ----------              
              do ir = 1, nr
                 bias_2d(i,j,l,ir) = bias_2d(i,j,l,ir) + del(n,ir)
                 rtms_2d(i,j,l,ir) = rtms_2d(i,j,l,ir) + del(n,ir)*del(n,ir)
              end do
              xcov_2d(i,j,l)       = xcov_2d(i,j,l)    + del(n,1) *del(n,nx)
              nobs_2d(i,j,l)       = nobs_2d(i,j,l)    + 1.


!        Upper-air (3D) grid
!        -------------------
         else

!             Find vertical level
!             -------------------
              do kk = 1, nlevs
                 if ( lev(n) .le. p1(kk) .AND. 
     &                lev(n) .gt. p2(kk) ) then
                      k = kk
                      go to 21
                 end if
              end do

              print *, 'Accum: obs level is ', lev(n), ' hPa'
              call gritas_perror ('Accum: could not find level')

 21           continue

!             Accumulate
!             ----------              
              do ir = 1, nr
                 bias_3d(i,j,k,l,ir) = bias_3d(i,j,k,l,ir) + del(n,ir)
                 rtms_3d(i,j,k,l,ir) = rtms_3d(i,j,k,l,ir) + del(n,ir)*del(n,ir)
              end do
              xcov_3d(i,j,k,l)       = xcov_3d(i,j,k,l)    + del(n,1) *del(n,nx)
              nobs_3d(i,j,k,l)       = nobs_3d(i,j,k,l)    + 1.

           end if

 10     continue



        end subroutine GrITAS_Accum_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Norm_ --- Normalize time mean statistics
! 
! !INTERFACE:
!
       subroutine GrITAS_Norm_ ( nr, nrmlz,
     &                           l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                           l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,
     &                           ntresh )
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::      nr           ! number of distinct residual series
      logical, intent(in) ::      nrmlz        ! normalize or not
      integer, intent(in) ::      l2d          ! actual number of 2d grids
      integer, intent(in) ::      l3d          ! actual number of 3d grids
      integer, OPTIONAL, intent(in) :: ntresh  ! tresh hold to compute means

! !INPUT/OUTPUT PARAMETERS:
!
                                                      ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(im,jm,l2d,nr)    !   time means
      real, intent(inout) ::  rtms_2d(im,jm,l2d,nr)    !   root time mean square
      real, intent(inout) ::  xcov_2d(im,jm,l2d,nr)    !   x-(co)variance
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                      ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(im,jm,km,l3d,nr) !   time means
      real, intent(inout) ::  rtms_3d(im,jm,km,l3d,nr) !   root time mean square
      real, intent(inout) ::  xcov_3d(im,jm,km,l3d,nr) !   x-(co)variance
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box
!
! !OUTPUT PARAMETERS:
!
      real, intent(inout) ::  stdv_2d(im,jm,l2d,nr)    ! 2D root time mean square
      real, intent(inout) ::  stdv_3d(im,jm,km,l3d,nr) ! 3D root time mean square

!
! !DESCRIPTION:  Normalizes the time mean statistics (such as
!                dividing the accumulated O-F by the number of
!  observations), set to 'missing value' the gridboxes with no
!  data.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  26Mar98  da Silva   Change normalization: n=1 is now OK for means.
!                      I needed this for being able to use the same
!                      code for synoptic time averages (GrISAS).
!  24Jun04  Todling    Added protection for nlev<2.
!  12oct05  da Silva   Added ntresh to impose a min number of obs
!  28Oct06  Todling    Generalized to handle cross-terms.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_Norm_'

      integer i, j, k, l, m, nx, ir
      real    tmp(nr), xtmp, n, ntresh_
   

      if ( nr > 2 ) call gritas_perror ('Accum: could not handle more then 2 residuals')
      nx = max(1,nr)

      ntresh_ = 0.                  ! by default no effect
      if ( present(ntresh) ) then
           ntresh_ = ntresh
      end if

!     Surface (2D) grids
!     ------------------
      m=1
      do l = 1, l2d
         do j = 1, jm
            do i = 1, im

               n = nobs_2d(i,j,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = bias_2d(i,j,l,ir) / m                                ! mean(x)
                       tmp(ir)           = rtms_2d(i,j,l,ir) / m                                ! mean(x**2)
                       rtms_2d(i,j,l,ir) = sqrt(tmp(ir))
                       tmp(ir)           = abs(tmp(ir)  - bias_2d(i,j,l,ir) *bias_2d(i,j,l,ir)) ! mean(x'**2)
                    end do
                    xtmp             = xcov_2d(i,j,l,1) / m                                    ! <x,y>
                    xcov_2d(i,j,l,1) = xtmp                                                    ! save xRMS in slot 1
                    xtmp             = xtmp  - bias_2d(i,j,l,1) *bias_2d(i,j,l,nx)             ! <x',y'>
               else
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = AMISS
                       rtms_2d(i,j,l,ir) = AMISS
                    end do
                       xcov_2d(i,j,l,1)  = AMISS
               end if

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr 
                       stdv_2d(i,j,l,ir)  = sqrt( m * tmp(ir) / (m-1) )
                    end do
                       xcov_2d(i,j,l,nx)  =       m * xtmp    / (m-1) 
               else
                    do ir = 1, nr 
                       stdv_2d(i,j,l,ir)  = AMISS
                    end do
                       xcov_2d(i,j,l,nx)  = AMISS
               end if

            end do
         end do
      end do

      if( nlevs<2 ) return

!     Upper-air (3D) grids
!     --------------------

      m=1
      do l = 1, l3d
         do k = 1, nlevs
            do j = 1, jm
               do i = 1, im
                  
               n = nobs_3d(i,j,k,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = bias_3d(i,j,k,l,ir) / m                        ! mean(x) and mean(y)
                       tmp(ir)             = rtms_3d(i,j,k,l,ir) / m                        ! mean(x**2) and mean(y**2)
                       rtms_3d(i,j,k,l,ir) = sqrt(tmp(ir))
                       tmp(ir)  = abs(tmp(ir)  - bias_3d(i,j,k,l,ir) *bias_3d(i,j,k,l,ir))  ! mean(x'**2) and mean(y'**2)
                    end do
                    xtmp               = xcov_3d(i,j,k,l,1) / m                             ! <x,y>
                    xcov_3d(i,j,k,l,1) = xtmp                                               ! save xRMS in slot 1 
                    xtmp               = xtmp  - bias_3d(i,j,k,l,1) *bias_3d(i,j,k,l,nx)    ! <x',y'>
               else
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = AMISS
                       rtms_3d(i,j,k,l,ir) = AMISS
                    end do
                       xcov_3d(i,j,k,l,1)  = AMISS
               endif

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    stdv_3d(i,j,k,l,1:nr)  = sqrt( m * tmp(1:nr) / (m-1) )
                    xcov_3d(i,j,k,l,nx)    =       m * xtmp      / (m-1)
               else
                    stdv_3d(i,j,k,l,1:nr) = AMISS
                    xcov_3d(i,j,k,l,nx)   = AMISS
               end if

               end do
            end do
         end do
      end do

      end subroutine GrITAS_Norm_

      end module m_gritas_grids
