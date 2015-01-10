subroutine diagcorD(kind_cov,kt,n_x, qr_x,qr_y,qr_z,		&
			    qm_x,qm_y,qm_z, ql_x,ql_y,		&
			    ktab,jtab,				&
		    Mtyp, corD, istat				)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: diagcorD - form a packed diagonal block correlation matrix
!
! !INTERFACE: @interface
!
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

use cormats
use config, only	: ktus, ktvs, ktuu, ktvv
use config, only	: kind_covF, kind_covS, kind_covV
use m_stdio,only	: stderr
use FEalpha_imat
implicit none

!-----------------------------------------------------------------------
!   Arguments
!  ===========
!@interface
integer, intent(in)	:: kind_cov	! which correlation matrix
integer, intent(in)	:: kt		! the variable type
integer, intent(in)	:: n_x		! dimension of the matrix

  real,    intent(in)	:: qr_x(n_x)	! x of qr
  real,    intent(in)	:: qr_y(n_x)	! y
  real,    intent(in)	:: qr_z(n_x)	! z
  real,    intent(in)	:: qm_x(n_x)	! x of d(qr)/dm
  real,    intent(in)	:: qm_y(n_x)	! y
  real,    intent(in)	:: qm_z(n_x)	! z
  real,    intent(in)	:: ql_x(n_x)	! x of d(qr)/dl
  real,    intent(in)	:: ql_y(n_x)	! y

  integer, intent(in)	:: ktab(n_x)	! indices to vfecH tables
  integer, intent(in)	:: jtab(n_x)	! weights of vfecH tables (tmp)

character*1, intent(out) :: Mtyp	! type of the matrix
real,    intent(out)	:: corD(n_x*(n_x+1)/2)	! a packed cor. matrix
integer, intent(out)	:: istat		! status
!@end/interface
!-----------------------------------------------------------------------
!   Parameters
!  ============
	character(len=*), parameter	:: myname='diagcorD'

	real,allocatable :: ql_z(:)
	real,allocatable :: corMM(:),corLL(:),corML(:,:)

	integer i,j,ij

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   For different variable types, different functions are used.

	Mtyp='E'
	istat=-1
	if( kt.ne.ktuu.and.kt.ne.ktvv .and.	&
	    kt.ne.ktus.and.kt.ne.ktvs		) return

	istat=0
	Mtyp='U'

	select case(kind_cov)
	case(kind_covF)

	  allocate(ql_z(n_x),corMM(n_x*(n_x+1)/2),		&
		corLL(n_x*(n_x+1)/2),corML(n_x,n_x), stat=istat	)
	  if(istat.ne.0) then
	    write(stderr,'(2a,i4)') myname,		&
		': allocate() error, stat =',istat
	    return
	  endif

	  ql_z(:)=0.

	  call fDDcor1(kind_covF,			&
	    n_x, qr_x,qr_y,qr_z,qm_x,qm_y,qm_z,ktab,	&
	    corMM					)

	  call fDDcorx(kind_covF,			&
	    n_x, qr_x,qr_y,qr_z,qm_x,qm_y,qm_z,ktab,	&
	    n_x, qr_x,qr_y,qr_z,ql_x,ql_y,ql_z,ktab,	&
	    corML					)

	  call fDDcor1(kind_covF,			&
	    n_x, qr_x,qr_y,qr_z,ql_x,ql_y,ql_z,		&
	    ktab, corLL					)

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

	  deallocate(ql_z,corMM,corLL,corML)

	case(kind_covS)

	  if(kt.eq.ktuu.or.kt.eq.ktus) then

	    call fDDcor1(kind_covS,			&
		n_x, qr_x,qr_y,qr_z,qm_x,qm_y,qm_z,	&
		ktab, corD				)

	  elseif(kt.eq.ktvv.or.kt.eq.ktvs) then
	    allocate(ql_z(n_x), stat=istat)
	    if(istat.ne.0) then
	      write(stderr,'(2a,i4)') myname,		&
		': allocate() error, stat =',istat
	      return
	    endif
	    ql_z(1:n_x)=0.

	    call fDDcor1(kind_covS,			&
		n_x, qr_x,qr_y,qr_z,ql_x,ql_y,ql_z,	&
		ktab, corD				)

	    deallocate(ql_z)
	  endif

	case(kind_covV)

	  if(kt.eq.ktuu.or.kt.eq.ktus) then
	    allocate(ql_z(n_x), stat=istat)
	    if(istat.ne.0) then
	      write(stderr,'(2a,i4)') myname,		&
		': allocate() error, stat =',istat
	      return
	    endif
	    ql_z(1:n_x)=0.

	    call fDDcor1(kind_covV,			&
		n_x, qr_x,qr_y,qr_z,ql_x,ql_y,ql_z,	&
		ktab, corD				)

	    deallocate(ql_z)
	  elseif(kt.eq.ktvv.or.kt.eq.ktvs) then

	    call fDDcor1(kind_covV,			&
		n_x, qr_x,qr_y,qr_z,qm_x,qm_y,qm_z,	&
		ktab, corD				)

	  endif

	end select
!_______________________________________________________________________
	end
!.
