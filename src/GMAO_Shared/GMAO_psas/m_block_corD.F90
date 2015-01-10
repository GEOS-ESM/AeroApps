!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corD - block diagonal correlation of vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_block_corD
      implicit none
      private	! except

      public :: diagcorD

	interface diagcorD; module procedure	&
	  diagcorD_Phi,	&
	  diagcorD_Oth
	end interface

! !REVISION HISTORY:
! 	06Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_block_corD'
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: diagcorD_Phi - form a packed diagonal block cor. matrix
!
! !INTERFACE: @interface

    subroutine diagcorD_Phi(kind_cov,kt,n_x,	&
	ktab,jtab,qr,qm,ql,kl, Mtyp, corD, istat)

      use config, only	: ktus, ktvs, ktuu, ktvv
      use config, only	: kind_covF
      use m_kt_corF,only : fDDcor1,fDDcorx
      use FEalpha_imat,only : Aum_imat,Aul_imat
      use FEalpha_imat,only : Avm_imat,Avl_imat
      use m_die,only : perr

      implicit none

      integer,intent(in) :: kind_cov	! which correlation matrix
      integer,intent(in) :: kt		! the variable type
      integer,intent(in) :: n_x		! dimension of the matrix

      integer,intent(in) :: ktab(n_x)	! indices to Aum tables
      integer,intent(in) :: jtab(n_x)	! indices to Aum tables

      real,   intent(in) :: qr(3,n_x)	! x,y,z of qr
      real,   intent(in) :: qm(3,n_x)	! x,y,z of d(qr)/dm
      real,   intent(in) :: ql(3,n_x)	! x,y,z of d(qr)/dl
      integer,intent(in) :: kl(  n_x)	! indices to hfecH/vfecH tables

      character(len=1),intent(out) :: Mtyp	! type of the matrix
      real   ,intent(out) :: corD(n_x*(n_x+1)/2) ! a packed cor. matrix
      integer,intent(out) :: istat		! status

! !DESCRIPTION:
!	diagcorD() forms a diagonal block of a correlation matrix for
!	a given data type `kt'.  The output matrix `corD' is of forecast
!	error correlation, stored in a packed form, representing the
!	upper triangle of the matrix with row being the first running
!	index.  `Mtyp' is either a "U" for the upper triangle with
!	`istat'=0 or "E" for an unimplemented `kt' with `istat'=-1.
!
!	Sigma matrix and the transformation matrix are considered
!	commutative, based on the assumption that the variances for the
!	two components of the gradiant d/dm and d/dl are the same at
!	every location.  Notice that the sigma matrix is the variances
!	of gradiant components, not of the wind components u and v.
!
!	Further more
!
! !SEE ALSO: fHHcor1(), fDDcor1(), and fQQcor1().
!	
! !BUGS:
!	- derived data types are expected to be used in future version.
!
! !REVISION HISTORY:
!	22Feb96  - Jing G.	- adopted from diagcorF() for mass field
!		balanced <u,u> and <v,v> correlation matrix
!	25Jan96  - Jing G.	- adopted partially the interface from
!		diagcor() with a new structure, functionalities, and a
!		new prolog.
!_______________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diagcorD_Phi'

!-----------------------------------------------------------------------
!   Parameters
!  ============

	real,allocatable :: corMM(:),corLL(:),corML(:,:)

	integer i,j,ij

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   For different variable types, different functions are used.

	Mtyp='E'
	istat=-1
	if(kind_cov/=kind_covF) return

	if( kt.ne.ktuu.and.kt.ne.ktvv .and.	&
	    kt.ne.ktus.and.kt.ne.ktvs		) return

	istat=0
	Mtyp='U'

	  allocate(corMM(n_x*(n_x+1)/2),		&
		corLL(n_x*(n_x+1)/2),corML(n_x,n_x), stat=istat	)
	  if(istat.ne.0) then
	    call perr(myname_,'allocate()',istat)
	    return
	  endif

	  call fDDcor1(kind_covF,n_x,qr,qm,kl,corMM)
	  call fDDcorx(kind_covF,n_x,qr,qm,kl,n_x,qr,ql,kl,corML)
	  call fDDcor1(kind_covF,n_x,qr,ql,kl,corLL)

	  if(kt.eq.ktuu.or.kt.eq.ktus) then

	    do j=1,n_x
	      ij=j*(j-1)/2	! where the last column ends

			! ~ (corD(i,j) = ..., i=1,j)
	      do i=1,j

		corD(ij+i) =	&
		  Aum_imat(ktab(i),jtab(i)) *	&
		    ( corMM(ij+i)*Aum_imat(ktab(j),jtab(j)) +	&
		      corML(i, j)*Aul_imat(ktab(j),jtab(j)) ) + &
		  Aul_imat(ktab(i),jtab(i)) *	&
		    ( corML(j, i)*Aum_imat(ktab(j),jtab(j)) +	&
		      corLL(ij+i)*Aul_imat(ktab(j),jtab(j)) )

	      end do	! i=1,j
	    end do	! j=1,n_x

	  elseif(kt.eq.ktvv.or.kt.eq.ktvs) then

	    do j=1,n_x
	      ij=j*(j-1)/2	! where the last column ends

			! ~ (corD(i,j) = ..., i=1,j)
	      do i=1,j

		corD(ij+i) =	&
		  Avm_imat(ktab(i),jtab(i)) *	&
		    ( corMM(ij+i)*Avm_imat(ktab(j),jtab(j)) +	&
		      corML(i, j)*Avl_imat(ktab(j),jtab(j)) ) +	&
		  Avl_imat(ktab(i),jtab(i)) *	&
		    ( corML(j, i)*Avm_imat(ktab(j),jtab(j)) +	&
		      corLL(ij+i)*Avl_imat(ktab(j),jtab(j)) )

	      end do	! i=1,j
	    end do	! j=1,n_x

	  endif

	  deallocate(corMM,corLL,corML,stat=istat)
		if(istat/=0) then
		  call perr(myname_,'deallocate()',istat)
		  return
		endif

!_______________________________________________________________________
	end subroutine diagcorD_Phi

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: diagcorD_Oth - form a packed diagonal block cor. matrix
!
! !INTERFACE:

    subroutine diagcorD_Oth(kind_cov,kt,n_x,	&
	qr,qm,ql,kl, Mtyp, corD, istat	)

      use config, only	: ktus, ktvs, ktuu, ktvv
      use config, only	: kind_covS, kind_covV
      use m_kt_corF,only : fDDcor1
      implicit none

      integer, intent(in) :: kind_cov	! which correlation matrix
      integer, intent(in) :: kt		! the variable type
      integer, intent(in) :: n_x	! dimension of the matrix

      real,    intent(in) :: qr(3,n_x)	! x,y,z of qr
      real,    intent(in) :: qm(3,n_x)	! x,y,z of d(qr)/dm
      real,    intent(in) :: ql(3,n_x)	! x,y,z of d(qr)/dl
      integer, intent(in) :: kl(  n_x)	! indices to hfecH/vfecH tables

      character(len=1),intent(out) :: Mtyp	! type of the matrix
      real,    intent(out) :: corD(n_x*(n_x+1)/2) ! a packed cor. matrix
      integer, intent(out) :: istat		! status

! !DESCRIPTION:
!	diagcorD() forms a diagonal block of a correlation matrix for
!	a given data type `kt'.  The output matrix `corD' is of forecast
!	error correlation, stored in a packed form, representing the
!	upper triangle of the matrix with row being the first running
!	index.  `Mtyp' is either a "U" for the upper triangle with
!	`istat'=0 or "E" for an unimplemented `kt' with `istat'=-1.
!
!	Sigma matrix and the transformation matrix are considered
!	commutative, based on the assumption that the variances for the
!	two components of the gradiant d/dm and d/dl are the same at
!	every location.  Notice that the sigma matrix is the variances
!	of gradiant components, not of the wind components u and v.
!
!	Further more
!
! !SEE ALSO: fHHcor1(), fDDcor1(), and fQQcor1().
!	
! !BUGS:
!	- derived data types are expected to be used in future version.
!
! !REVISION HISTORY:
!	22Feb96  - Jing G.	- adopted from diagcorF() for mass field
!		balanced <u,u> and <v,v> correlation matrix
!	25Jan96  - Jing G.	- adopted partially the interface from
!		diagcor() with a new structure, functionalities, and a
!		new prolog.
!_______________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diagcorD_Oth'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   For different variable types, different functions are used.

	Mtyp='E'
	istat=-1
	if(.not. (kind_cov==kind_covS	.or.	&
		  kind_cov==kind_covV)	) return

	if( kt.ne.ktuu.and.kt.ne.ktvv .and.	&
	    kt.ne.ktus.and.kt.ne.ktvs		) return

	istat=0
	Mtyp='U'

	select case(kind_cov)
	case(kind_covS)

	  select case(kt)
	  case(ktuu,ktus)

	    call fDDcor1(kind_covS, n_x, qr,qm,kl, corD)

	  case(ktvv,ktvs)

	    call fDDcor1(kind_covS, n_x, qr,ql,kl, corD)

	  end select

	case(kind_covV)

	  select case(kt)
	  case(ktuu,ktus)

	    call fDDcor1(kind_covV, n_x, qr,ql,kl, corD)

	  case(ktvv,ktvs)

	    call fDDcor1(kind_covV, n_x, qr,qm,kl, corD)

	  end select

	end select
!_______________________________________________________________________
	end subroutine diagcorD_Oth

end module m_block_corD
