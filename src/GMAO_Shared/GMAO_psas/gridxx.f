**************************************************************************
*                                                                        *
*                  NASA/Goddard Space Flight Center                      *
*                     Laboratory for Atmospheres                         *
*                      Data Assimilation Office                          *
*                            Code 910.3                                  *
*                                                                        *
*            PHYSICAL-SPACE STATISTICAL ANALYSIS SYSTEM (PSAS)           *
*                                                                        *
**************************************************************************


*........................... ROUTINE PROLOGUE .............................
* !BOP
*
* !FILE: gridxx.f
*
* !ROUTINE: gridxx0
* 
* !DESCRIPTION: Initialize analysis grid parameters.
*
* !CALLING SEQUENCE: 
*
*         call GRIXX0
*
* !INPUT PARAMETERS: none
*
* !OUTPUT PARAMETERS: none
*
* !SYSTEM ROUTINES: none.
*
* !FILES USED: none.
*
* !WRITTEN BY: Jim Pf, 02jun93 
* 
* !REVISION HISTORY:
!
!	21Oct99	- Jing Guo
!		. Replaced qea.h dependency with m_qea(.F90).
!		. Removed stdio.h dependency by removing write()
!		  statements.
!	01Oct99	- Jing Guo
!		. Removed gridxx.h, since it is not used in here.
!
*  18Jul95 - Jing G.	- Splited gridxx.h and qea.h for variable grid
*			  settings
*  16nov94   A. da Silva  Implemented quasi-equal area grid option.
*                         Initialization of mandatory levels moved
*                         to gridxx.h for portability.
*
* !EOP
*.........................................................................

      subroutine GRIDXX0

	use m_inpak90, only : lablin,fltget
	use m_qea,     only : eaytresh
	implicit none

*     Parameters and table for analysis grid.
*     --------------------------------------
      call LABLIN ( 'latitude_treshold_for_equal_area_grid:' )
      eaytresh = abs(FLTGET ( 90. ))

      return
      end

*........................... ROUTINE PROLOGUE .............................
* !BOP
*
* !FILE: gridxx.f
*
* !ROUTINE: eagrid_set
* 
* !DESCRIPTION: 
*
*     Returns quasi-equal area grid points on the sphere poleward
*  of ytresh, regular lat/lon grid elsewhere. If ytresh = 90, then
*  a regular lat/lon grid is returned.
*
* !CALLING SEQUENCE:
*
*      call EAGRID_set ( ier, glon, glat, ndecl, ngrid, 
*    .               idim, jdim )
*
* !INPUT PARAMETERS: 
*
*      integer  ndecl     ! declared dimension in calling program
*      integer  idim      ! no. of longitudes at equator
*      integer  jdim      ! no. of equally space latitudes
*      integer  npole     ! no. zonal grid points at pole
*      real     ytresh    ! only change grid poleward of ytresh
*
* !OUTPUT PARAMETERS:
*
*      integer  ier           ! 0 if OK, .ne. 0 otherwise
*      integer  ngrid         ! total number of grid points in
*                             ! quasi equal area grid.
*      real     glon(ndecl)   ! longitude in degree
*      real     glat(ndecl)   ! latitude in degree
*
* !SEE ALSO: n/a
*
* !SYSTEM ROUTINES: none. 
*
* !FILES USED: none.
*
* !WRITTEN BY: A. da Silva, 16nov94.
* 
* !REVISION HISTORY:
!
!	08Nov99	- Jing Guo
!		. renamed EAGRID() to EAGRID_set, since the original
!		  procedure name conflicts with the module datatype
!		  name of m_EAGrid(.F90).
!	01Oct99	- Jing Guo
!		. Removed gridxx.h, since it is not used here.
*
* !EOP
*.........................................................................

      subroutine EAGRID_set ( ier, glon, glat, ndecl, ngrid, 
     .                    idim, jdim )

	use m_qea,only : qea_init
	use m_qea,only : lea_beg,lea_len,ea_lon
	use m_qea,only : npole,eaytresh
	use m_qea,only : j_south,j_north

	implicit none
      
*     On output
*     ---------      
      integer  ier           ! 0 if OK, .ne. 0 otherwise
      integer  ngrid         ! total number of grid points in
                             ! quasi equal area grid.

      integer  ndecl     ! declared dimension in calling program

      real     glon(ndecl)   ! longitude in degree
      real     glat(ndecl)   ! latitude in degree

*     On input:
*     --------
      integer  idim      ! no. of longitudes at equator
      integer  jdim      ! no. of equally space latitudes

! Local variables

	integer :: i,j,l
	integer :: nlon,neq
	real :: dlon,dlat,ylat,alat
	real :: pi,ytresh,d2r
*
*     Returns quasi-equal area grid points on the sphere.
*     This version defines a peculiar quasi-equal area grid at
*     the polar caps (northward of ytresh below). To simplify the
*     interpolation, the poles have a full set of longitudes
*     (usually the longitudinal mesh size decreases as the poles
*     are approached.)

	call qea_init(idim,jdim)

      ytresh = eaytresh
      pi = 4. * atan ( 1. )
      d2r = pi / 180.

*     Generate quasi-equal area grid
*     ------------------------------
      neq  = idim - npole
      dlat = 180. / ( jdim - 1 )
      l = 0
      do 10 j = 1, jdim

         ylat = -90. + (j-1) * dlat
         alat = d2r * ylat

*        Decide how many long. gridpoints for this latitude
*        --------------------------------------------------
         if ( abs(ylat) .gt. ytresh ) then
              nlon = nint ( npole + neq * cos(alat) )
         else
              nlon = idim
         end if
         if (j.eq.1 .or. j.eq.jdim) nlon = idim ! no interp at poles
         dlon = 360. / nlon

         lea_beg(j) = l + 1    ! pointer for interpolation
         lea_len(j) = nlon

*        Create gridpoints for this latitude
*        -----------------------------------
         do 20 i = 1, nlon
            l = l + 1
            if ( l .gt. ndecl ) then
               ier = 1
               return
            end if
            glon(l) = ( i - 1 ) * dlon - 180.  ! this will be sorted
            glat(l) = ylat                     ! this will be sorted
            ea_lon(l) = glon(l)                ! this will be kept as is

 20      continue

 10   continue


*     No. of q.e.a. grid points
*     -------------------------
      ngrid = l


*     Compute beginning and end latitudinal index of
*      regular lat/lon grid (based on ytresh) needed
*      for interpolation
*     ----------------------------------------------
      do 50 j = 1, jdim
         ylat = -90. + (j-1) * dlat
         j_south = j
         if ( abs(ylat) .le. ytresh ) go to 51
 50   continue
 51   continue
      do 60 j = jdim, 1, -1
         ylat = -90. + (j-1) * dlat
         j_north = j
         if ( abs(ylat) .le. ytresh ) go to 61
 60   continue
 61   continue

      ier = 0

*     All done.
*     ---------
      ier = 0
      return
      end


*........................... ROUTINE PROLOGUE .............................
* !BOP
*
* !FILE: gridxx.f
*
* !ROUTINE: ea2ll
* 
* !DESCRIPTION: Converts from quasi-equal area to a lat/lon grid.
*
* !CALLING SEQUENCE:
*
* !INPUT PARAMETERS:
*
* !OUTPUT PARAMETERS:
*
* !SEE ALSO:
*
* !SYSTEM ROUTINES: 
*
* !FILES USED: 
*
* !WRITTEN BY: A. da Silva, 16nov94.
* 
* !REVISION HISTORY:
!
!	01Oct99	- Jing Guo
!		. Removed gridxx.h, since it is not used here.
*
* !EOP
*.........................................................................

      subroutine EA2LL ( All, idim, jdim, Aea, ngrid )
	use m_qea,only : lea_beg,lea_len,ea_lon
	use m_qea,only : eaytresh
	use m_qea,only : j_south,j_north
	use m_stdio,only : stderr
	use m_die,only : die
	implicit none

	integer,intent(in) :: idim
	integer,intent(in) :: jdim
      
*     Output:
*     ------
      real,intent(out) :: All(idim,jdim)      ! array on lat/lon grid

*     Input:
*     -----
      integer,intent(in) :: ngrid ! No. of quasi-equal area grid points
      real,   intent(in) :: Aea(ngrid) ! array on quasi-equal area grid

      character(len=*), parameter :: myname='ea2ll'

*     Local work space
*     ----------------
	integer nw,ier
	real, allocatable :: b(:),c(:),d(:),x(:),y(:)	! spline coef.
	real :: dlon,dlat,ytresh,u
	integer :: i,j,l,n

	real,external :: SEVAL

	nw=idim+1

      ytresh = eaytresh

*     For latitudes equatorward of ytresh no interpolation is needed
*     --------------------------------------------------------------
      l = lea_beg(j_south)
      do j = j_south, j_north
         do i = 1, idim
            All(i,j) = Aea(l)
            l = l + 1
         end do
      end do


*     No interpolation needed if entire grid is lat/lon
*     -------------------------------------------------

      if ( ytresh .ge. 90.-90./(jdim-1) ) return

	allocate( b(nw), c(nw), d(nw), x(nw), y(nw), stat=ier)
	if(ier.ne.0) then
	  write(stderr,'(2a,i5)') myname,
     &	    ': allocate() error, stat =',ier
	  call die(myname)
	endif

*     Southern polar cap
*     ------------------
      dlon = 360. / idim
      do j = 2, j_south - 1                 ! no south pole

c-       ylat = -90. + (j-1)*dlat

         l = lea_beg(j)
         n = lea_len(j)

*        Copy q.e.a. data to local vector repeating first
*         longitudinal point to ensure periodicity
*        ------------------------------------------------
         do i = 1, n
            x(i) = ea_lon(l+i-1)
            y(i) = Aea(l+i-1)
         end do
         n = n + 1
         x(n) = +180.      ! -180 = + 180 - date line.
         y(n) = Aea(l)

*        Compute spline coefficients
*        ---------------------------
         call SPLINE ( b, c, d, n, x, y )

*        Compute interpolated values
*        ---------------------------
         do i = 1, idim
            u = -180. + (i-1) * dlon
            All(i,j) =  SEVAL ( u, x, y, n, b, c, d )
         end do

      end do


*     Northern polar cap
*     ------------------
      do j = j_north+1, jdim-1                  ! no north pole

c-       ylat = -90. + (j-1)*dlat

         l = lea_beg(j)
         n = lea_len(j)

*        Copy q.e.a. data to local vector repeating first
*         longitudinal point to ensure periodicity
*        ------------------------------------------------
         do i = 1, n
            x(i) = ea_lon(l+i-1)
            y(i) = Aea(l+i-1)
         end do
         n = n + 1
         x(n) = +180.      ! -180 = + 180 - date line.
         y(n) = Aea(l)

*        Compute spline coefficients
*        ---------------------------
         call SPLINE ( b, c, d, n, x, y )

*        Compute interpolated values
*        ---------------------------
         do i = 1, idim
            u = -180. + (i-1) * dlon
            All(i,j) =  SEVAL ( u, x, y, n, b, c, d )
         end do

      end do

*     No interpolation needed at poles
*     --------------------------------
      do i = 1, idim
         All(i,1)    = Aea(i)                 ! south pole
         All(i,jdim) = Aea(ngrid-idim+i)      ! north pole
      end do

      deallocate(b,c,d,x,y)

      end






