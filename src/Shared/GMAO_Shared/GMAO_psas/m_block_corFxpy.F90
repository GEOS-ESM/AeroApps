!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_corFxpy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_block_corFxpy
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

  character(len=*),parameter :: myname='m_block_corFxpy'


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

  subroutine sCxpy(kind_cov, kr,kt,ln,qr,qd,kl,nvecs,x,y,ierr)

    use config, only	: ktpm => ktus	! dp/dm gradient type
    use config, only	: ktpl => ktvs	! dp/dl gradient type
    use config, only	: ktHm => ktuu	! dH/dm gradient type
    use config, only	: ktHl => ktvv	! dH/dl gradient type
    use config, only	: ktHH,ktslp,ktqq

    use m_mvcorF_matx, only : sHHmcxpy, sHH1cxpy
    use m_mvcorF_matx, only : sDDmcxpy, sDD1cxpy
    use m_kt_uvcorF, only : sQQmcxpy, sQQ1cxpy
    use m_die,only : perr

    implicit none
    
    integer,intent(in) :: kind_cov	! kind of cov. matrix
    integer,intent(in) :: kr		! region index
    integer,intent(in) :: kt		! variable type
    integer,intent(in) :: ln		! size
    
    real   ,intent(in) :: qr(3,ln)	! unit vectors of locations
    real   ,intent(in) :: qd(3,ln,*)	! unit vectors of directions
    integer,intent(in) :: kl(  ln)	! level references
    
    integer,intent(in) :: nvecs	! number of input vectors
    real   ,intent(in) :: x(nvecs,ln,*)	! multiplier vectors
    real   ,intent(inout):: y(nvecs,ln,*)	! product C times x
    
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
    case (kind_covF)
       diagcorX='diagcorPhi'
    case (kind_covS)
       diagcorX='diagcorPsi'
    case (kind_covV)
       diagcorX='diagcorChi'
    case default
       call perr(myname_,'unexpected kind_cov',kind_cov)
       ierr=-1
       return
    end select

       ! Determine variable for consideration "p", "H", or "Q"

    Select Case (kt)
    Case (ktHm, ktpm, ktHl, ktpl)
      select case(nvecs)
      case(1)
	call sDD1cxpy(kind_cov,	ln, qr, qd, kl, nvecs, x, y )

      case default
	call sDDmcxpy(kind_cov,	ln, qr, qd, kl, nvecs, x, y )
      end select

    Case (ktHH, ktslp)
      select case(nvecs)
      case(1)
	call sHH1cxpy(kind_cov, ln, qr, kl, nvecs, x, y )

      case default
	call sHHmcxpy(kind_cov, ln, qr, kl, nvecs, x, y )
      end select

    Case (ktqq)
      select case(nvecs)
      case(1)
	call sQQ1cxpy(ln, qr, kl, nvecs, x, y )

      case default
	call sQQmcxpy(ln, qr, kl, nvecs, x, y )
      end select

    Case Default
       ierr = -1
    End Select

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

  subroutine xCxpy(kind_cov,      &
		kri,kti,lni,qri,qdi,kli, &
		krj,ktj,lnj,qrj,qdj,klj, &
		nvecs,xi,yj,xj,yi,	 &
		ierr)

    use m_mvcorF_matx, only : xHHmcxpy, xHH1cxpy
    use m_mvcorF_matx, only : xHDmcxpy, xHD1cxpy
    use m_mvcorF_matx, only : xDDmcxpy, xDD1cxpy
    use m_kt_uvcorF, only : xQQmcxpy, xQQ1cxpy
    use m_die  , only : perr
    implicit none
    
    integer,intent(in) :: kind_cov	! kind of the cov. matrix
    
    integer,intent(in) :: kri		! region index of i-vector(row)
    integer,intent(in) :: kti		! variable type of i-vector
    integer,intent(in) :: lni		! size of i-vector
    
    real   ,intent(in) :: qri(3,lni)	! unit vectors of locations
    real   ,intent(in) :: qdi(3,lni,*)	! unit vectors of directions
    integer,intent(in) :: kli(  lni)	! level references
    
    integer,intent(in) :: krj		! region index of j-vector(col.)
    integer,intent(in) :: ktj		! variable type of j-vector
    integer,intent(in) :: lnj		! size of j-vector
    
    real   ,intent(in) :: qrj(3,lnj)	! unit vectors of locations
    real   ,intent(in) :: qdj(3,lnj,*)	! unit vectors of directions
    integer,intent(in) :: klj(  lnj)	! level references
    
    integer,intent(in) :: nvecs	! number of input vectors
    
    real   ,intent(in)    :: xi(nvecs,lni,*)	! i_vectors
    real   ,intent(inout) :: yj(nvecs,lnj,*)	! C'*i_vectors + y_j
    
    real   ,intent(in)    :: xj(nvecs,lnj,*)	! j_vectors
    real   ,intent(inout) :: yi(nvecs,lni,*)	! C *j_vectors + y_i
    
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

    character(len=*),parameter :: myname_=myname//'::xCxpy_'

!   Local vars.
!  =============
    integer ivec
    
    character(len=len('offdcor___'))	:: offdcorX
    
    Character(Len=1) :: var_i, var_j
    Character(Len=2) :: HD_Select
    
!_______________________________________________________________________
	! Argument checking, kind_cov

    select case(kind_cov)
    case (kind_covF)
       offdcorX='offdcorPhi'
    case (kind_covS)
       offdcorX='offdcorPsi'
    case (kind_covV)
       offdcorX='offdcorChi'
    case default
       call perr(myname_,'unexpected kind_cov',kind_cov)
       ierr=-1
       return
    end select
    
!-----------------------------------------------------------------------
	! If any vector is zero-sized, both yi and yj are untouched.
    ierr=0
    if(nvecs <= 0) return
    if(lni <=0 .or. lnj <=0) return

!-----------------------------------------------------------------------
        ! Determine variable for consideration "p", "H", or "Q"
    var_i = Var(kti, ierr)
       If (ierr /= 0) Return
    var_j = Var(ktj, ierr)
       If (ierr /= 0) Return
	   
	   ! Verify consistency
       If (var_i /= var_j) Then
	  ierr = -1
	  Return
       End If

    Select Case (var_i)

    Case ('Q')
      select case(nvecs)
      case(1)
       call xQQ1cxpy(lni,qri,kli, lnj,qrj,klj, nvecs, xi, yj, xj, yi )

      case default
       call xQQmcxpy(lni,qri,kli, lnj,qrj,klj, nvecs, xi, yj, xj, yi )
      end select

    Case ('p','H')

           ! Determine whether to use gradient ("D") or var ("H")
       HD_Select(1:1) = H_or_D(kti)
       HD_Select(2:2) = H_or_D(ktj)
       
       Select Case (HD_Select)
	  
       Case ('HH') ! case of <H,H> correlation
	 select case(nvecs)
	 case(1)
	   call xHH1cxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,klj, nvecs, xi, yj, xj, yi )

	 case default
	   call xHHmcxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,klj, nvecs, xi, yj, xj, yi )
	 end select

       Case ('HD') ! case of <H,M> or <H,L> correlation
	 select case(nvecs)
	 case(1)
	   call xHD1cxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,qdj,klj, nvecs, xi, yj, xj, yi )

	 case default
	   call xHDmcxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,qdj,klj, nvecs, xi, yj, xj, yi )
	 end select

       Case ('DH') ! case of <M,H> or <L,H> correlation
	 select case(nvecs)
	 case(1)
	   call xHD1cxpy(kind_cov,	&
		lnj,qrj,klj, lni,qri,qdi,kli, nvecs, xj, yi, xi, yj )

	 case default
	   call xHDmcxpy(kind_cov,	&
		lnj,qrj,klj, lni,qri,qdi,kli, nvecs, xj, yi, xi, yj )
	 end select

       Case ('DD') ! case of <M,M>, <M,L>, <L,M>, or <L,L> correlation
	 select case(nvecs)
	 case(1)
	   call xDD1cxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,qdj,klj, nvecs, xi, yj, xj, yi )

	 case default
	   call xDDmcxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,qdj,klj, nvecs, xi, yj, xj, yi )
	 end select
       End Select
    End Select

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

  subroutine rCxpy( kind_cov,                   &
                     kri,kti,lni,qri,qdi,kli,	&
                     krj,ktj,lnj,qrj,qdj,klj,	&
		     nvecs, x, y, ierr)

    use m_mvcorF_matx, only : rHHmcxpy, rHH1cxpy
    use m_mvcorF_matx, only : rHDmcxpy, rHD1cxpy
    use m_mvcorF_matx, only : rDHmcxpy, rDH1cxpy
    use m_mvcorF_matx, only : rDDmcxpy, rDD1cxpy
    use m_kt_uvcorF, only : rQQmcxpy, rQQ1cxpy
    use m_die  , only : perr
    implicit none

    integer, intent(in) :: kind_cov	! kind of the cov. matrix

      ! Attrbutes of rows (i) or y
      
    integer, intent(in) :: kri	! variable type of i-vector
    integer, intent(in) :: kti	! variable type of i-vector
    integer, intent(in) :: lni	! size of i-vector
    real,    intent(in) :: qri(3,lni)	! unit vectors
    real,    intent(in) :: qdi(3,lni,*)	! unit vectors
    integer, intent(in) :: kli(  lni)   ! indices to hcor/vcor tables
    
      ! Attributes of columns (j) or x
      
    integer, intent(in) :: krj	! variable type of j-vector
    integer, intent(in) :: ktj	! variable type of j-vector
    integer, intent(in) :: lnj	! size of j-vector
    real,    intent(in) :: qrj(3,lnj)	! unit vectors
    real,    intent(in) :: qdj(3,lnj,*)	! unit vectors
    integer, intent(in) :: klj(  lnj)   ! indices to hcor/vcor tables
    
    integer, intent(in) :: nvecs	! number of input vectors
    
    real,intent(in)    :: x(nvecs,lnj,*) ! column vectors
    real,intent(inout) :: y(nvecs,lni,*) ! row vectors
    
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
    integer ivec
    character(len=len('offdcor___'))	:: offdcorX
    character(len=1) :: var_i, var_j
    character(len=2) :: HD_Select

!_______________________________________________________________________
	! Argument checking, kind_cov

    select case(kind_cov)
    case (kind_covF)
       offdcorX='offdcorPhi'
    case (kind_covS)
       offdcorX='offdcorPsi'
    case (kind_covV)
       offdcorX='offdcorChi'
    case default
       ierr=-1
       call perr(myname_,'unexpected kind_cov',kind_cov)
       return
    end select
!-----------------------------------------------------------------------
	! Special cases

    ierr=0
    if(nvecs <= 0) return
    if(lni <= 0) return	! a case of zero-row matrix
    if(lnj <= 0) return	! a case of zero-column matrix
    
!-----------------------------------------------------------------------

	   ! Determine variable for consideration "p", "H", or "Q"
    var_i = Var(kti, ierr)
       If (ierr /= 0) Return
    var_j = Var(ktj, ierr)
       If (ierr /= 0) Return
	   
	   ! Verify consistency
       If (var_i /= var_j) Then
	  ierr = -1
	  Return
       End If


    Select Case (var_i)

    Case ('Q')
      select case(nvecs)
      case(1)
	call rQQ1cxpy( lni, qri, kli, lnj, qrj, klj, nvecs, x, y )
      
      case default
	call rQQmcxpy( lni, qri, kli, lnj, qrj, klj, nvecs, x, y )
      end select

    Case ('p','H')

           ! Determine whether to use gradient ("D") or var ("H")
       HD_Select(1:1) = H_or_D(kti)
       HD_Select(2:2) = H_or_D(ktj)
       
       Select Case (HD_Select)
	  
       Case ('HH') ! case of <H,H> correlation
	 select case(nvecs)
	 case(1)
	   call rHH1cxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,klj, nvecs, x, y )

	 case default
	   call rHHmcxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,klj, nvecs, x, y )
	 end select

       Case ('HD') ! case of <H,M> or <H,L> correlation
	 select case(nvecs)
	 case(1)
	   call rHD1cxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,qdj,klj, nvecs, x, y )

	 case default
	   call rHDmcxpy(kind_cov,	&
		lni,qri,kli, lnj,qrj,qdj,klj, nvecs, x, y )
	 end select

       Case ('DH') ! case of <M,H> or <L,H> correlation
	 select case(nvecs)
	 case(1)
	   call rDH1cxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,klj, nvecs, x, y )

	 case default
	   call rDHmcxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,klj, nvecs, x, y )
	 end select

       Case ('DD') ! case of <M,M>, <M,L>, <L,M>, or <L,L> correlation
	 select case(nvecs)
	 case(1)
	   call rDD1cxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,qdj,klj, nvecs, x, y )

	 case default
	   call rDDmcxpy(kind_cov,	&
		lni,qri,qdi,kli, lnj,qrj,qdj,klj, nvecs, x, y )
	 end select
       End Select
    End Select
!_______________________________________________________________________
  end subroutine rCxpy

Character(Len=1) Function Var(kt, ierr)
  use config, only	: ktpm => ktus	! dp/dm gradient type
  use config, only	: ktpl => ktvs	! dp/dl gradient type
  use config, only	: ktHm => ktuu	! dH/dm gradient type
  use config, only	: ktHl => ktvv	! dH/dl gradient type
  use config, only	: ktHH,ktslp,ktqq

  Implicit None
  Integer, Intent(In) :: kt
  Integer, Intent(Out) :: ierr

  
  Select Case (kt)
  Case (ktHH, ktHm, ktHl)
     Var = 'H'
     ierr = 0
  Case (ktslp, ktpm, ktpl)
     Var = 'p'
     ierr = 0
  Case (ktqq)
     Var = 'Q'
  Case Default
     ierr = -1   
  End Select
  Return
End Function Var

Character(Len=1) Function H_or_D(kt)
  use config, only	: ktpm => ktus	! dp/dm gradient type
  use config, only	: ktpl => ktvs	! dp/dl gradient type
  use config, only	: ktHm => ktuu	! dH/dm gradient type
  use config, only	: ktHl => ktvv	! dH/dl gradient type
  use config, only	: ktHH,ktslp,ktqq

  Implicit None
  Integer, Intent(In) :: kt
  
  Select Case (kt)
  Case (ktHH, ktslp)
     H_or_D = 'H'
  Case (ktHm,ktHl,ktpm,ktpl)
     H_or_D = 'D'
  End Select
  Return

End Function H_or_D

end module m_block_corFxpy
!.
