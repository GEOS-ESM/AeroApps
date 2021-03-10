!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corU - uncorrelated obs. err. matrix generators
!
! !DESCRIPTION:
!
!	m_block_corU() is a module of PSAS' correlation function driver
!	subroutines.  It defines the interfaces of the following
!	correlation function driver subroutines,
!
!	  diagcor : forming a diagonal block of the observation error
!		    correlation matrix with uncorrelated horizontal
!		    correlation functions;
!	  offdcor : (TBD) forming a off-diagonal block of the
!		    observation error correlation matrix without
!		    horizontal correlation functions.
!
! !INTERFACE:

    module m_block_corU
      use config, only	: ktHH
      implicit none
      private	! except

      public :: diagCor

      interface diagCor; module procedure diagcor_; end interface

#ifdef TODO
      public :: offdCor
      interface offdCor; module procedure offdcor_; end interface
#endif

! !REVISION HISTORY:
! 	16Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- redefined to a true module and customized for corU
!
!	25Jan95 - J. Guo	- added diagcorU() interface.  Some
!		corrections also made according to the changes in other
!		include subroutines.
!
! 	15Dec95 - J. Guo	- programed and added the prolog
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_block_corU'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diagcor_ - form a packed diagonal block correlation matrix
!
! !DESCRIPTION:
!
!	diagcor_() forms a diagonal block of a correlation matrix for a
!	given variable type `kt'.  The output matrix `corU' is of
!	observation error correlation with only horizontal correlation
!	between observations in the same profile (ks).  The matrix is
!	stored in a packed form, representing the upper triangle of the
!	matrix with rows being the first running index.  `Mtyp' is
!	either "U" for the upper triangle or "I" for any unimplemented
!	`kt'. `istat' is always zero.
!
! !INTERFACE:

    subroutine diagcor_(kind_cov,kt,len, kx, ks, kl, Mtyp, corU, istat)
      use m_kt_corU,only : uHHcor1
      implicit none

      integer,intent(in) :: kind_cov	! which correlation matrix

      integer,intent(in) :: kt	! the variable type
      integer,intent(in) :: len	! dimension of the matrix
      integer,intent(in) :: kx(  len)	! instrument indices
      integer,intent(in) :: ks(  len)	! sounding indices
      integer,intent(in) :: kl(  len)	! hcor/vcor table reference

      character(len=*),intent(out) :: Mtyp	! matrix type
      real,intent(out) :: corU(len*(len+1)/2)	! a packed matrix
      integer,intent(out) :: istat		! status

! !REVISION HISTORY:
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
!   For different variable types, different functions are used.  It may
!   need to be changed in the near future.

	istat=0
	Mtyp='I'	! assuming an identity matrix.  The actual type
			! is to be determined by uHHcor1() under certain
			! conditions.

	if(len<=0) return	! allow zero-sized vector/matrix

	select case(kt)
	case(ktHH)	! <H,H> correlation

	  call uHHcor1(len,kx,ks,kl,Mtyp,corU)

	case default

	  corU=0.		! all zeroes, but ...
	  do i=1,len
	    corU(i*(i+1)/2)=1.
	  end do

	end select
!_______________________________________________________________________
end subroutine diagcor_

#ifdef TODO
.. offdcor_(..)
#endif

end module m_block_corU
