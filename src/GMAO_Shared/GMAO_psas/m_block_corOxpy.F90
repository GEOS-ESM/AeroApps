!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corOxpy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_block_corOxpy
      implicit none
      private	! except

      public :: sCxpy	! for a block symmetric block matrix
      public :: xCxpy	! for a symmetric pair of block matrices
      public :: rCxpy	! for a generic block matrix

!!$      interface sCxpy; module procedure sCxpy_; end interface
!!$      interface xCxpy; module procedure xCxpy_; end interface
!!$      interface rCxpy; module procedure rCxpy_; end interface

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

  character(len=*),parameter :: myname='m_block_corOxpy'

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

    subroutine sCxpy(kind_cov, kr,kt,ln,kx,qr,kl,nvecs,x,y,ierr)

      use m_kt_corO, only : sHHmcxpy, sHH1cxpy
      use config, only : ktHH
      use m_die,only : perr
      implicit none

      integer,intent(in) :: kind_cov	! kind of cov. matrix
      integer,intent(in) :: kr		! region index
      integer,intent(in) :: kt		! variable type
      integer,intent(in) :: ln		! size

      integer,intent(in) :: kx(  ln)	! instrument indices
      real   ,intent(in) :: qr(3,ln)	! unit vectors of locations
      integer,intent(in) :: kl(  ln)	! level references

      integer,intent(in)    :: nvecs	! number of input vectors
      real   ,intent(in)    :: x(nvecs,ln)	! multiplier vectors
      real   ,intent(inout) :: y(nvecs,ln)	! product C times x

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
	character(len=len('diagcor___'))	:: diagcorX

!_______________________________________________________________________
	! allow zero-sized output

	ierr=0
	if(ln<=0 .or. nvecs<=0) return

!_______________________________________________________________________
	! Argument checking, kind_cov

	select case(kind_cov)
	case (kind_covC)
	  diagcorX='diagcorR_c'
	case default
	  call perr(myname_,'unexpected kind_cov',kind_cov)
	  ierr=-1
	  return
	end select
!_______________________________________________________________________

	select case (kt)
		!--------------------------------
	case (ktHH)
	  select case(nvecs)
	  case(1)
	   call sHH1cxpy(ln, kx, qr, kl, nvecs, x, y )

	  case default
	   call sHHmcxpy(ln, kx, qr, kl, nvecs, x, y )
	  end select

	Case Default
	   ! block is identity matrix
	   y(:,:) = y(:,:) + x(:,:)

		!--------------------------------
	end select

!-----------------------------------------------------------------------

end subroutine sCxpy


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xCxpy_ - performs both C(ai,aj)*xj+yi and C(aj,ai)*xi+yj
!
! !DESCRIPTION:
!
! 	xCxpy_() performs both C'*xi+yj and C*xj+yi in the same
!	call.  C is a block matrix of a symmetric covariance matrix.
!
! !INTERFACE:

    subroutine xCxpy(kind_cov,					&
		kri,kti,lni,kxi,qri,kli,			&
		krj,ktj,lnj,kxj,qrj,klj,			&
		nvecs,xi,yj,					&
		      xj,yi,					&
		ierr)

      use m_die  , only : perr
      use m_kt_corO, only : xHHmcxpy, xHH1cxpy
      use config, only : ktHH

      implicit none

      integer,intent(in) :: kind_cov	! kind of the cov. matrix

      integer,intent(in) :: kri		! region index of i-vector(row)
      integer,intent(in) :: kti		! variable type of i-vector
      integer,intent(in) :: lni		! size of i-vector

      integer,intent(in) :: kxi(  lni)	! instrument indices
      real   ,intent(in) :: qri(3,lni)	! unit vectors of locations
      integer,intent(in) :: kli(  lni)	! level references

      integer,intent(in) :: krj		! region index of j-vector(col.)
      integer,intent(in) :: ktj		! variable type of j-vector
      integer,intent(in) :: lnj		! size of j-vector

      integer,intent(in) :: kxj(  lnj)	! instrument indices
      real   ,intent(in) :: qrj(3,lnj)	! unit vectors of locations
      integer,intent(in) :: klj(  lnj)	! level references

      integer,intent(in) :: nvecs	! number of input vectors

      real   ,intent(in)    :: xi(nvecs,lni)	! i_vectors
      real   ,intent(inout) :: yj(nvecs,lnj)	! C'*i_vectors + y_j

      real   ,intent(in)    :: xj(nvecs,lnj)	! j_vectors
      real   ,intent(inout) :: yi(nvecs,lni)	! C *j_vectors + y_i

      integer,intent(out) :: ierr		! status

! !REVISION HISTORY:
!       23Feb01 - Tom Clune <clune@sgi.com>
!               . significant redesign to directly call
!                 low level routines that comput elements
!                 The BLAS call is "fused" into those routines
!                 thereby eliminating the need to store the entire
!                 block.  (Improves cache usage!)
!	21Aug00	- Jing Guo
!		. customized xCxpy() as a special interface for CorF
!
! 	25Mar99 - Jing Guo <guo@thunder>
!		- renamed from Cprodx() to xCxpy()
!
! 	22Mar99 - Jing Guo <guo@thunder>
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
!  22jun94  - Jim Pf.  - Added prologue, BLAS matrix routine
!  07jan94  - Meta S.  - Added pass of trig lat/lon
!  28may93  - J. Searl - Modification for dynamic storage on CRAY
!  20feb93  - Jim Pf.  - original program
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::xCxpy'

!   Local vars.
!  =============
	character(len=len('offdcor___'))	:: offdcorX


!-----------------------------------------------------------------------
	! If any vector is zero-sized, both yi and yj are untouched.

	ierr=0
	if(nvecs <= 0) return
	if(lni <=0 .or. lnj <=0) return

!_______________________________________________________________________
	! Argument checking, kind_cov

	select case(kind_cov)
	case (kind_covC)
	  offdcorX='offdcorR_c'
	case default
	  call perr(myname_,'unexpected kind_cov',kind_cov)
	  ierr=-1
	  return
	end select

!-----------------------------------------------------------------------

	If (kti == ktHH .and. ktj == ktHH) Then

	  select case(nvecs)
	  case(1)
	    Call xHH1cxpy(lni,kxi,qri,kli,	&
		lnj,kxj,qrj,klj,nvecs,xi,yj,xj,yi)

	  case default
	    Call xHHmcxpy(lni,kxi,qri,kli,	&
		lnj,kxj,qrj,klj,nvecs,xi,yj,xj,yi)
	  end select

	Else

	   ! sparsity should have taken care of normal cases
	   ierr=-1

	End If

!_______________________________________________________________________

end subroutine xCxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rCxpy_ - performs y=Cx+y
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rCxpy( kind_cov,		&
		kri,kti,lni,kxi,qri,kli,	&
		krj,ktj,lnj,kxj,qrj,klj,	&
		nvecs,x,y, ierr)

	use m_kt_corO, only : rHHmcxpy, rHH1cxpy
	use config, only : ktHH
	use m_die,only : perr

	implicit none

	integer, intent(in) :: kind_cov	! kind of the cov. matrix

		! Attrbutes of rows (i) or y

	integer, intent(in) :: kri	! variable type of i-vector
	integer, intent(in) :: kti	! variable type of i-vector
	integer, intent(in) :: lni	! size of i-vector
	integer, intent(in) :: kxi(  lni)	! instrument indices
	real,    intent(in) :: qri(3,lni)	! unit vectors
	integer, intent(in) :: kli(  lni) ! indices to hcor/vcor tables

		! Attributes of columns (j) or x

	integer, intent(in) :: krj	! variable type of j-vector
	integer, intent(in) :: ktj	! variable type of j-vector
	integer, intent(in) :: lnj	! size of j-vector
	integer, intent(in) :: kxj(  lnj)	! instrument indices
	real,    intent(in) :: qrj(3,lnj)	! unit vectors
	integer, intent(in) :: klj(  lnj) ! indices to hcor/vcor tables

	integer, intent(in) :: nvecs	! number of input vectors

        real,intent(in)    :: x(nvecs,lnj) ! column vectors
	real,intent(inout) :: y(nvecs,lni) ! row vectors

	integer, intent(out):: ierr	! status

! !REVISION HISTORY:
!       23Feb01 - Tom Clune <clune@sgi.com>
!               . significant redesign to directly call
!                 low level routines that comput elements
!                 The BLAS call is "fused" into those routines
!                 thereby eliminating the need to store the entire
!                 block.  (Improves cache usage!)
!	21Aug00	- Jing Guo
!		. customized rCxpy() as a special interface for CorF
!
! 	21Dec98 - Jing Guo <guo@thunder>
!		- renamed from gCprodx() to rCxpy()
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
!  22jun94  - Jim Pf.  - Added prologue, BLAS matrix routine
!  07jan94  - Meta S.  - Added pass of trig lat/lon
!  28may93  - J. Searl - Modification for dynamic storage on CRAY
!  20feb93  - Jim Pf.  - original program
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rCxpy_'

!   Local vars.
!  =============
	character(len=len('offdcor___'))	:: offdcorX

!_______________________________________________________________________
	! Argument checking, kind_cov

	select case(kind_cov)
	case (kind_covC)
	  offdcorX='offdcorR_c'
	case default
	  call perr(myname_,'unexpected kind_cov',kind_cov)
	  ierr=-1
	  return
	end select
!-----------------------------------------------------------------------
	! Special cases

	ierr=0
	if(nvecs <= 0) return
	if(lni <= 0) return	! a case of zero-row matrix
	if(lnj <= 0) return	! a case of zero-column matrix

!-----------------------------------------------------------------------

	If (kti == ktHH .and. ktj == ktHH) Then
	  select case(nvecs)
	  case(1)
	    Call rHH1cxpy(lni,kxi,qri,kli,lnj,kxj,qrj,klj,nvecs,x,y)

	  case default
	    Call rHHmcxpy(lni,kxi,qri,kli,lnj,kxj,qrj,klj,nvecs,x,y)
	  end select

	Else
	   ! sparsity should have taken care of normal cases
	   ierr=-1
	End If
!_______________________________________________________________________
end subroutine rCxpy


end module m_block_corOxpy
!.
