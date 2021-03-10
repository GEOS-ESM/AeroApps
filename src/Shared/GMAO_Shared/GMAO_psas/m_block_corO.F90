!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corO - block correlation operators of Fcst.Err.
!
! !DESCRIPTION:
!
!	m_block_corO() is a module of PSAS' correlation function driver
!	subroutines.  It defines the interfaces of the following
!	correlation function driver subroutines,
!
!	  diagcor : forming a diagonal block of the forecast error
!		    correlation matrix;
!	  offdcor : forming a off-diagonal block of the forecast error
!		    correlation matrix;
!
! !INTERFACE:

    module m_block_corO

      use config, only	: ktpm => ktus	! dp/dm gradient type
      use config, only	: ktpl => ktvs	! dp/dl gradient type
      use config, only	: ktHm => ktuu	! dH/dm gradient type
      use config, only	: ktHl => ktvv	! dH/dl gradient type
      use config, only	: ktHH,ktslp,ktqq

      implicit none
      private	! except

      public :: diagCor
      public :: offdCor

      interface diagCor; module procedure diagcor_; end interface
      interface offdCor; module procedure offdcor_; end interface

! !REVISION HISTORY:
! 	16Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- redefined to a true module and customized for corF
!
!	25Jan95 - J. Guo	- added diagcorU() interface.  Some
!		corrections also made according to the changes in other
!		include subroutines.
!
! 	15Dec95 - J. Guo	- programed and added the prolog
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_block_corO'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diagcor_ - form a packed diagonal block correlation matrix
!
! !DESCRIPTION:
!
!	diagcor_() forms a diagonal block of a correlation matrix for
!	a given data type `kt'.  The output matrix `corM' is an error
!	correlation matrix stored in a packed form, representing the
!	upper triangle of the matrix with row being the first running
!	index.  `Mtyp' is either a "U" for the upper triangle or "I"
!	for an identity matrix, both with `istat'=0.
!
! !INTERFACE:

    subroutine diagcor_(kind_cov,kt,len, kx,qr,kl, Mtyp,corM,istat)

      use m_kt_corO,only : oHHcor1

      implicit none

      integer,intent(in) :: kind_cov	! which correlation matrix
      integer,intent(in) :: kt		! the variable type
      integer,intent(in) :: len		! dimension of the matrix

      integer,intent(in) :: kx(  len)	! level references
      real   ,intent(in) :: qr(3,len)	! unit vectors at locations
      integer,intent(in) :: kl(  len)	! level references

      character(len=*),intent(out) :: Mtyp	! type of the matrix
      real,intent(out) :: corM(len*(len+1)/2)	! a packed matrix
      integer,intent(out) :: istat		! status

! !REVISION HISTORY:
!
! 	18Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- redesigned with new data structure for unit vectors
!		- redesigned as a module procedure
!		- reformatted the prolog
!
!	25Jan96  - Jing G.	- adopted partially the interface from
!		diagcor() with a new structure, functionalities, and a
!		new prolog.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diagcor_'
  integer :: i

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   For different variable types, different functions are used.

	istat=0

	Mtyp='I'	! assuming an identity matrix.  The actual type
			! is to be determined by oHHcor1() under certain
			! conditions.

	if(len<=0) return	! allow zero-sized vector/matrix

	select case(kt)
	case(ktHH)	! <H,H> correlation

	  call oHHcor1(len,kx,qr,kl, Mtyp, corM)

	case default

	  corM=0.		! all zeroes, but ...
	  do i=1,len
	    corM(i*(i+1)/2)=1.
	  end do

	end select
!_______________________________________________________________________
end subroutine diagcor_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: offdcor_ - form an off-diagonal block correlation matrix
!
! !DESCRIPTION:
!
!	offdcor_ forms a off-diagonal block of a correlation matrix for
!	the given variable types (kti vs. ktj).  The returned matrix
!	in corM may be zeros (Mtyp='Z'), or normal form (Mtyp='N'), or
!	in error (Mtyp='E', istat=-1)
!
! !INTERFACE:

    subroutine offdcor_(kind_cov,		&
			kti,leni,kxi,qri,kli,	&
			ktj,lenj,kxj,qrj,klj,	&
			Mtyp, corM,	istat	)

      use m_kt_corO,only : oHHcorx

      implicit none

      integer,intent(in) :: kind_cov	! which correlation matrix

      integer,intent(in) :: kti		! variable type of the rows
      integer,intent(in) :: leni	! row dimension of the matrix
      integer,intent(in) :: kxi(  leni)	! instrument types
      real   ,intent(in) :: qri(3,leni)	! unit vectors
      integer,intent(in) :: kli(  leni)	! indices to hcor/vcor tables

      integer,intent(in) :: ktj		! variable type of the rows
      integer,intent(in) :: lenj	! row dimension of the matrix
      integer,intent(in) :: kxj(  lenj)	! instrument types
      real   ,intent(in) :: qrj(3,lenj)	! unit vectors
      integer,intent(in) :: klj(  lenj)	! indices to hcor/vcor tables

      character(len=1),intent(out) :: Mtyp	! if is a "Transpose"
      real   ,intent(out) :: corM(leni,lenj)	! the correlation matrix
      integer,intent(out) :: istat		! status

! !REVISION HISTORY:
! 	16Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- redesigned with new data structure for unit vectors
!		- redesigned as a module procedure
!		- reformatted the prolog
!
! 	28Nov95 - J. Guo	- added the prolog
!
!  19Jan95  - Jing G.	- Added wobs tables to pass pindx2() values to
!			  ??cor1() and ??corx() routines.  One could use
!			  rlevs for the same purpose to reduce the over-
!			  head, since rlevs has no real purpose in this
!			  subroutine and subsequent routines.
!  25jul94  - Meta S.  - Added QQ correlation
!  22jun94  - Jim Pf.  - Added header
!  18jan94  - Meta S.  - Added sea level pressure correlation
!  04jan94  - Meta S.  - Passing trig values as arguments to UU- VVcor
!  19dec93  - Meta S.  - Added UU and VV wind correlations
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::offdcor_'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   For different variable types, different functions are used.

!-----------------------------------------------------------------------
	istat=0
	Mtyp='N'

	if(leni<=0) return	! allow zero-sized vector/matrix
	if(lenj<=0) return	! allow zero-sized vector/matrix

	if( (kti.eq.ktHH .and.ktj.eq.ktHH ) ) then

!		 case of <H,H> correlation
!		---------------------------
	  Mtyp='N'		! the actual type is to be determined by
				! oHHcorx().

	  call oHHcorx(	leni,kxi,qri,kli,	&
			lenj,kxj,qrj,klj,	&
			Mtyp,corM		)

	else
		! sparsity should have taken care of normal cases

	  Mtyp='E'
	  istat=-1
	endif
!_______________________________________________________________________
end subroutine offdcor_
end module m_block_corO
