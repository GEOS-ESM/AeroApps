C  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
C			  to follow the guide lines.
c  modification for dynamic storage on CRAY                     05/28/93
c  03oct94 - Implemented CRAY specifics with IFDEFs (A. da S.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Compute (x,y,z) coordinates from (rlons,rlats).

      subroutine ll2xyz( rlons, rlats, length, x, y, z, ierr )

	use m_stdio,only : stderr
	implicit none
c.......................................................................
c.... Argument declarations.

      integer      length
      real         rlons(length)
      real         rlats(length)
      real         x(length)
      real         y(length)
      real         z(length)
      integer      ierr
c.......................................................................
c.... Local storage.

	real,allocatable :: coslat(:)
	real :: deg
	real :: rlon,clon,slon
	real :: rlat,clat,slat
	integer :: n

	character*6 myname
	parameter  (myname='ll2xyz')

c.......................................................................
c.... Statement function.

      logical      llbad
      llbad(n) = (abs(rlons(n)).gt.180.0).or.(abs(rlats(n)).gt.90.0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   ..Consistency check

	allocate(coslat(length),stat=ierr)
	if(ierr.ne.0) then
	  write(stderr,'(2a,i5)') myname,
     &		': allocate() error, stat =',ierr
	  return
	endif

      ierr = 0

c.......................................................................
c.... Loop over the points.

      do 100 n = 1, length

c.......................................................................
c....... Check that longitude and latitude are in range.

         if( llbad(n) ) then
            print *, '  ll2xyz:  latitude or longitude out of range'
            print *, '  ll2xyz:  n = ', n, '  rlats(n) = ', rlats(n),
     $                                     '  rlons(n) = ', rlons(n)
            ierr = 1
	    deallocate(coslat)
            return
         endif

  100 continue

c.......................................................................
c.... Compute sines and cosines of longitude and latitude.

      !call qtrig( length, rlons, x, y )
      !call qtrig( length, rlats, coslat, z )

c.......................................................................
c.... Compute x and y coordinates, ( z = sin(lat) ).

      deg=4.*atan(1.)/180.
      do 200 n = 1, length
	 rlon=rlons(n)*deg
	 rlat=rlats(n)*deg
	 clon=cos(rlon)
	 slon=sin(rlon)
	 clat=cos(rlat)
	 slat=sin(rlat)
	 x(n)=clat*clon
	 y(n)=clat*slon
	 z(n)=slat
         !x(n) = coslat(n)*x(n)
         !y(n) = coslat(n)*y(n)
  200 continue

      deallocate(coslat)
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
