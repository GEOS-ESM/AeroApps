!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corUxpy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_block_corUxpy
      implicit none
      private	! except

      public :: sCxpy	! for a block symmetric block matrix
      interface sCxpy; module procedure sCxpy_; end interface

! !REVISION HISTORY:
!       23Feb01 - Tom Clune <clune@sgi.com>
!               . significant redesign to directly call
!                 low level routines that comput elements
!                 The BLAS call is "fused" into those routines
!                 thereby eliminating the need to store the entire
!                 block.  (Improves cache usage!)
!	21Aug00	- Jing Guo
!		. created from m_block_symCxpy.F90 and
!		  m_block_recCxpy.F90
! 	25Mar99 - Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_block_corUxpy'

  include "kind_covs.h"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sCxpy_ - y=C*x+y, with covariance sub-matrix block C.
!
! !DESCRIPTION:
!
!	sCxpy_() operates as a real symmetric covariance matrix
!	multiply one or more vectors.
!
! !INTERFACE:

    subroutine sCxpy_(kind_cov, kr,kt,ln,kx,ks,kl,nvecs,x,y,ierr)

      use m_kt_corU, only : sHHmcxpy, sHH1cxpy
      use config, only : ktHH
      use m_die,only : perr
      implicit none

      integer,intent(in) :: kind_cov	! kind of cov. matrix
      integer,intent(in) :: kr		! region index
      integer,intent(in) :: kt		! variable type
      integer,intent(in) :: ln		! size

      integer,intent(in) :: ks(ln)	! instrument indices
      integer,intent(in) :: kx(ln)	! instrument indices
      integer,intent(in) :: kl(ln)	! level references

      integer,intent(in) :: nvecs	! number of input vectors
      real   ,intent(in) :: x(nvecs,ln)	! multiplier vectors
      real   ,intent(out):: y(nvecs,ln)	! product C times x

      integer,intent(out):: ierr	! status

! !REVISION HISTORY:
!       23Feb01 - Tom Clune <clune@sgi.com>
!               . significant redesign to directly call
!                 low level routines that comput elements
!                 The BLAS call is "fused" into those routines
!                 thereby eliminating the need to store the entire
!                 block.  (Improves cache usage!)
!	21Aug00	- Jing Guo
!		. customized sCxpy() as a special interface for CorF
!
! 	25Mar99 - Jing Guo <guo@dao.gsfc.nasa.gov>
!		- renamed from Cprod1() to sCxpy()
!
! 	22Mar99 - Jing Guo <guo@dao.gsfc.nasa.gov>
!		- changed from Cx to Cxpy
!		- Switched the 2 dimensions of both x and y
!
! 	04Dec95 - J. Guo
!		- implemented new correlation function interface
!		- converted to FORTRAN 90
!		- revised prolog
!
!  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
!			  to follow the guide lines.
!  19Jan95  - Jing G.	- Added wobs tables to pass pindx2() values to
!			  ??cor1() and ??corx() routines.  One could use
!			  rlevs for the same purpose to reduce the over-
!			  head, since rlevs has no real purpose in this
!			  subroutine and subsequent routines.
!  03oct94  - A. da Silva - Implemented CRAY specifics with IFDEFS.
!              Use of exit() instead of STOP.
!  22jun94  - Jim Pf.  - Added prologue, BLAS matrix call
!  07jan94  - Meta S.  - Added pass of trig lat/lon
!  28may93  - J. Searl - Modification for dynamic storage on CRAY
!  20feb93  - Jim Pf.  - Original program
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sCxpy_'

!-----------------------------------------------------------------------
!..local vars.
	integer ivec
	character(len=1) :: Mtyp
	character(len=len('diagcor___'))	:: diagcorX

!_______________________________________________________________________
	! allow zero-sized output

	ierr=0
	if(ln<=0 .or. nvecs<=0) return

!_______________________________________________________________________
	! Argument checking, kind_cov

	select case(kind_cov)
	case (kind_covU)
	  diagcorX='diagcorR_u'
	case default
	  call perr(myname_,'unexpected kind_cov',kind_cov)
	  ierr=-1
	  return
	end select
!_______________________________________________________________________

	If (kt == ktHH) Then
	  select case(nvecs)
	  case(1)
	    Call sHH1cxpy(ln, kx, ks, kl, nvecs, x, y )

	  case default
	    Call sHHmcxpy(ln, kx, ks, kl, nvecs, x, y)
	  end select

	Else
	   y(:,:) = y(:,:) + x(:,:)
	End If

!-----------------------------------------------------------------------

end subroutine sCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diagperr_ - error messenger
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine diagperr_(where,corx,kt,mtyp,ierr)
      use m_stdio,only : stderr
      implicit none
      character(len=*),intent(in) :: where
      character(len=*),intent(in) :: corx
      integer,intent(in) :: kt
      character(len=*),intent(in) :: Mtyp
      integer,intent(in) :: ierr


! !REVISION HISTORY:
! 	05Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diagperr_'
  character(len=16) :: ckt
  character(len=16) :: cerr

  write(ckt ,'(i16)') kt
  write(cerr,'(i16)') ierr

  write(stderr,*) where,			&
	': diagcor("',corx,'") error',		&
	', kt = ',trim(adjustl(ckt)),		&
	', Mtyp = "',Mtyp,			&
	'", stat =',trim(adjustl(cerr))

end subroutine diagperr_

end module m_block_corUxpy
!.
