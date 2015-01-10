!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_kt_CorU - correlation matrix generators for uncorrel.Obs.Err.
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_kt_CorU
      implicit none
      private	! except

      public :: uHHcor1
      public :: sHHmcxpy
      public :: sHH1cxpy

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Added new interface that avoids storing the block
!                 but rather applies elements to right/left vectors as
!                 individual elements are constructed.  Significant
!                 performance improvement.
!	29Aug00	- Jing Guo
!		- initial prototype/prolog/code
!		- combined uHHcor1.F into this module
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_kt_CorU'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: uHHcor1 - forms a packed symmetric <H,H> correlation matrix
!
! !DESCRIPTION:
!
!	uHHcor1 forms the upper triangle of a symmetric <height,height>
!	correlation matrix.  The returned matrix is in a packed form in
!	the order of ((i=1,j),j=1,len).  The matrix contains only the
!	uncorrelated components.
!
! !INTERFACE:

    subroutine uHHcor1(len,kx,ks,kl, Mtyp, corrU)
      implicit none

      integer,intent(in) :: len	! dimension of the matrix
      integer,dimension(len),intent(in) :: kx
      integer,dimension(len),intent(in) :: ks
      integer,dimension(len),intent(in) :: kl
      character(len=1),intent(inout) :: Mtyp		! "I" or "U"
      real,dimension(len*(len+1)/2),intent(out) :: corrU

! !REVISION HISTORY:
! 	29Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- Unified argument style and modular style with other
!		  corrlation matrix generator subroutines.
!
! 	09Nov95 - J. Guo	- rewrite the code and the prolog.  The
!		recoding 1) separated forecast error correlation matrix
!		from observation error correlation matrix;  2) used a
!		new structure in hope that it will be optimized better;
!		and 3) getting ready for the new changes (see BUGS:)
!		that may effect the results.
!
!  29Sep95  - Jing G.	- added #ifdef for dynamic dimension check.
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  27Mar95  - Jing G.	- Modified for text file based data tables.
!  19Jan95  - Jing G.	- Added wobs table to replace pindx2().  One
!			  could use rlevs for the same purpose to reduce
!			  the overhead, since rlevs has no real purpose
!			  in this subroutine.
!
!   7oct94  Meta S.  rearranged observation error correlation test
!                     to avoid index problem
!  10may94  Meta S.  change observation error correlation tabulation
!                     to use cos(dist/rade)
!  28apr94  Meta S.  change forecast error correlation tabulation
!                     to use cos(dist/rade), added prologue
!  02feb94  Meta S.  change radiosonde obs covariance criterion to specify
!                    -matching- lat/lon and alter 'coslim' definition
!                    (now included in "hHHotab.h", "hHHotab.h")
!  15dec93  Meta S.  add observation covariances, covariance tables
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::uHHcor1'

!-----------------------------------------------------------------------
	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj
	integer kxj,ksj
	integer jvtype
	real vercor
	logical correlated

!=======================================================================

	if(len.le.0) return

	corrU(1)=1.

	do j=2,len

	  ij=j*(j-1)/2

	  ksj=ks(j)
	  kxj=kx(j)

		! is there an uncorrelated part for the class?

	  jvtype=i_voecHu(kxj)
	  correlated=jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (ks).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      correlated=ks(i).eq.ksj.and.kx(i).eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then
	    Mtyp='U'		! a packed matrix, not "I" any more

	    kj = kl(j)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      corrU(l)=0.

			! ksi.eq.ksj may or may not implies kxi.eq.kxj.
			! Errors in this correlation model are
			! correlated only if they are from the same
			! observation platform.  e.g. the same rawind
			! sond.

	      if(ks(i).eq.ksj.and.kx(i).eq.kxj) then

		ki=kl(i)

		vercor=voecHH(ki,kj,jvtype)

		corrU(l)=vercor
	      endif	! ks(i).eq.ksj.and.kx(i).eq.kxj
	    end do		! i=1,j-1

	  else
	    do i=1,j-1
	      corrU(ij+i)=0.
	    end do
	  endif

	  corrU(ij+j)=1.	! where i=j
	end do		! j=2,len

end subroutine uHHcor1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sHHcxpy - applies a symmetric <H,H> correlation matrix
!
! !DESCRIPTION:
!
!	sHHcxpy applies a symmetric <height,height>
!	correlation matrix to right and left vectors.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sHHmcxpy(len,kx,ks,kl, nvecs, x, y)
      implicit none

      integer,intent(in) :: len	! dimension of the matrix
      integer,dimension(len),intent(in) :: kx
      integer,dimension(len),intent(in) :: ks
      integer,dimension(len),intent(in) :: kl
      integer, intent(in) :: nvecs
      real,dimension(nvecs,len), intent(in)    :: x
      real,dimension(nvecs,len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sHHmcxpy'

!-----------------------------------------------------------------------
	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj
	integer kxj,ksj
	integer jvtype
	real vercor
	logical correlated
	real corrU

!=======================================================================

	if(len.le.0) return

	y(:,1) = y(:,1) + x(:,1)

	do j=2,len

	  ij=j*(j-1)/2

	  ksj=ks(j)
	  kxj=kx(j)

		! is there an uncorrelated part for the class?

	  jvtype=i_voecHu(kxj)
	  correlated=jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (ks).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      correlated=ks(i).eq.ksj.and.kx(i).eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    kj = kl(j)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

			! ksi.eq.ksj may or may not implies kxi.eq.kxj.
			! Errors in this correlation model are
			! correlated only if they are from the same
			! observation platform.  e.g. the same rawind
			! sond.

	      if(ks(i).eq.ksj.and.kx(i).eq.kxj) then

		ki=kl(i)

		vercor=voecHH(ki,kj,jvtype)

		corrU=vercor

		y(:,i) = y(:,i) + corrU * x(:,j)
		y(:,j) = y(:,j) + corrU * x(:,i)

	      endif	! ks(i).eq.ksj.and.kx(i).eq.kxj
	    end do		! i=1,j-1

	  else
	     ! no contributions
	  endif

	  y(:,j) = y(:,j) + x(:,j) ! where i=j
	end do		! j=2,len

      end subroutine sHHmcxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sHHcxpy - applies a symmetric <H,H> correlation matrix
!
! !DESCRIPTION:
!
!	sHHcxpy applies a symmetric <height,height>
!	correlation matrix to right and left vectors.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sHH1cxpy(len,kx,ks,kl, nvecs, x, y)
      implicit none

      integer,intent(in) :: len	! dimension of the matrix
      integer,dimension(len),intent(in) :: kx
      integer,dimension(len),intent(in) :: ks
      integer,dimension(len),intent(in) :: kl
      integer, intent(in) :: nvecs
      real,dimension(len), intent(in)    :: x
      real,dimension(len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sHH1cxpy'

!-----------------------------------------------------------------------
	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj
	integer kxj,ksj
	integer jvtype
	real vercor
	logical correlated
	real corrU
!=======================================================================

	if(len.le.0) return

	y(1) = y(1) + x(1)

	do j=2,len

	  ij=j*(j-1)/2

	  ksj=ks(j)
	  kxj=kx(j)

		! is there an uncorrelated part for the class?

	  jvtype=i_voecHu(kxj)
	  correlated=jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (ks).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      correlated=ks(i).eq.ksj.and.kx(i).eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    kj = kl(j)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

			! ksi.eq.ksj may or may not implies kxi.eq.kxj.
			! Errors in this correlation model are
			! correlated only if they are from the same
			! observation platform.  e.g. the same rawind
			! sond.

	      if(ks(i).eq.ksj.and.kx(i).eq.kxj) then

		ki=kl(i)

		vercor=voecHH(ki,kj,jvtype)

		corrU=vercor

		y(i) = y(i) + corrU * x(j)
		y(j) = y(j) + corrU * x(i)

	      endif	! ks(i).eq.ksj.and.kx(i).eq.kxj
	    end do		! i=1,j-1

	  else
	     ! no contributions
	  endif

	  y(j) = y(j) + x(j) ! where i=j
	end do		! j=2,len

      end subroutine sHH1cxpy

end module m_kt_CorU
