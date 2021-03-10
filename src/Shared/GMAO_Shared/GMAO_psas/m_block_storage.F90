!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_block_storage
!
! !DESCRIPTION:
!
! !INTERFACE:

#include "assert.H"

module m_block_storage
  Use m_die,only : assert_
  implicit none
  private	! except
  
  Public :: reserve_storage
  Public :: get_handle
  Public :: release_handle
  Public :: clean

  Interface get_handle
     Module Procedure get_handle_diag_
!!$     Module Procedure get_handle_offdiag_
  End Interface
  Interface release_handle
     Module Procedure release_handle_diag_
!!$     Module Procedure release_handle_offdiag_
  End Interface

! !REVISION HISTORY:
!       29Oct01 - Tom Clune
!               - Added mall_mci/mall_mco calls
!	17Feb01	- Tom Clune
!		- initial prototype/prolog/code

!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_block_storage'

#include "kind_covs.h"
#include "kind_mats.h"

!!$  Integer, Parameter :: TOT_STORAGE  =  20 * 10**6 ! elements
!!$  Integer, Parameter :: TOT_STORAGE  =  5 * 10**6 ! elements
  Integer, Parameter :: TOT_STORAGE  =  3 * 10**6 ! elements
 
  Integer :: MAX_DIAG_STORAGE = 0
  Integer :: MAX_REGN_STORAGE = 0
  Integer :: MAX_GLOB_STORAGE = 0

  Integer :: Excess_Storage   = TOT_STORAGE

  Type ptr_storage
     Real, Pointer :: ptr(:,:)		! => Null()
  End Type ptr_storage

  Type diag_ptr_storage
     Real, Pointer :: ptr(:)		! => Null()
  End Type diag_ptr_storage

#include "ktmax.h"
  Integer, Parameter :: n_kind_cov = 3
  Integer, Parameter :: MAX_DIAG_SEGS = 10000
  Type (diag_ptr_storage), Save :: diag_ptr_corr(MAX_DIAG_SEGS)
  logical,save :: associated_diag_ptr_corr(MAX_DIAG_SEGS) = .false.

  Integer :: n_elem_GLOB = 0
  Integer :: n_elem_REGN = 0
  Integer :: n_elem_DIAG = 0

Contains

  Function reserve_storage(kmat, n_elem) Result (succeed)
    use m_die, only : warn
    Integer, Intent(In) :: kmat
    Integer, Intent(In) :: n_elem
    Logical :: succeed

    Character(len=*),parameter :: myname_=myname//'::reserve_storage'
    Integer, Dimension(4), Parameter :: allowed_kmat =	&
	(/ kind_Umat, kind_Rmat, kind_4mat, kind_5mat /)


    ASSERT(Any(kmat==allowed_kmat))

    Select Case (kmat)
    Case (kind_Umat) ! Diagonal - highest priority
       If (n_elem <= excess_storage) Then
	  MAX_DIAG_STORAGE = MAX_DIAG_STORAGE + n_elem
	  excess_storage = excess_storage - n_elem
	  succeed = .true.
       ElseIf (n_elem <= MAX_GLOB_STORAGE+excess_storage) Then
	  MAX_DIAG_STORAGE = MAX_DIAG_STORAGE + n_elem
	  MAX_GLOB_STORAGE = MAX_GLOB_STORAGE - (n_elem - excess_storage)
	  excess_storage = 0
	  succeed = .true.
	  Call warn(myname_,'DIAG borrowed from GLOB')
       ElseIf (n_elem <= MAX_REGN_STORAGE+MAX_GLOB_STORAGE+excess_storage) Then
	  MAX_DIAG_STORAGE = MAX_DIAG_STORAGE + n_elem
	  MAX_REGN_STORAGE = MAX_REGN_STORAGE - (n_elem - (MAX_GLOB_STORAGE + excess_storage))
	  excess_storage = 0
	  MAX_GLOB_STORAGE = 0
	  succeed = .true.
	  Call warn(myname_,'DIAG borrowed from REGN')
       Else
	  MAX_DIAG_STORAGE = MAX_DIAG_STORAGE + MAX_REGN_STORAGE + MAX_GLOB_STORAGE + excess_storage
	  excess_storage   = 0
	  MAX_GLOB_STORAGE = 0
	  MAX_REGN_STORAGE = 0
	  succeed = .false.
       End If

    Case (kind_Rmat) ! Regional-Diagonal - medium priority
       If (n_elem <= excess_storage) Then
	  MAX_REGN_STORAGE = MAX_REGN_STORAGE + n_elem
	  excess_storage = excess_storage - n_elem
	  succeed = .true.
       ElseIf (n_elem <= MAX_GLOB_STORAGE+excess_storage) Then
	  MAX_REGN_STORAGE = MAX_REGN_STORAGE + n_elem
	  MAX_GLOB_STORAGE = MAX_GLOB_STORAGE - (n_elem - excess_storage)
	  excess_storage = 0
	  succeed = .true.
	  Call warn(myname_,'REGN borrowed from GLOB')
       Else
	  MAX_REGN_STORAGE = MAX_REGN_STORAGE + MAX_GLOB_STORAGE + excess_storage
	  excess_storage   = 0
	  MAX_GLOB_STORAGE = 0
	  succeed = .false.
       End If
    Case (kind_4mat,kind_5mat) ! Global - low priority
       If (n_elem <= excess_storage) Then
	  MAX_GLOB_STORAGE = MAX_GLOB_STORAGE + n_elem
	  excess_storage = excess_storage - n_elem
	  succeed = .true.
       Else
	  MAX_GLOB_STORAGE = MAX_GLOB_STORAGE + excess_storage
	  excess_storage   = 0
	  MAX_GLOB_STORAGE = 0
	  succeed = .false.
       End If

    End Select

    Return
    
  End Function reserve_storage

  Subroutine release_handle_diag_(handle)
    use m_mall, only : mall_mco, mall_ison
    use m_die,  only : perr
    use m_mpout, only : mpout_log
    Implicit None

    Integer, Intent(In) :: handle

    Integer :: n_elem, kcov, ierr
    Character(Len=*), Parameter :: myname_ = myname // '::release_handle_diag_'

       If (mall_ison()) Call mall_mco(diag_ptr_corr(handle)%ptr, myname)
    Deallocate(diag_ptr_corr(handle)%ptr)
       nullify(diag_ptr_corr(handle)%ptr)
       associated_diag_ptr_corr(handle)=.false.

  End Subroutine release_handle_diag_
!_______________________________________________________________________

  Function get_handle_diag_(handle,nn,compute,release,ierr) Result (ptr)
    use m_die,only : perr
    use m_mall, only : mall_ison, mall_mci
    use m_mpout, only : mpout_log
    use m_stdio,only : stderr
    Implicit None

    Integer, Intent(In) :: handle
    Integer, Intent(In) :: nn

    Logical, Intent(Out) :: compute
    Logical, Intent(Out) :: release
    Integer, Intent(InOut) :: ierr

    Real, Pointer        :: ptr(:)

    Integer :: n_elem, kcov
    Character(Len=*), Parameter :: myname_ = myname // '::get_handle_diag_'

    If (handle > MAX_DIAG_SEGS) Then
       write(stderr,*) myname_, &
            & 'too many segments'
       ierr = -1
       return
    End If

    If (associated_diag_ptr_corr(handle)) Then

       ptr  => diag_ptr_corr(handle)%ptr
!!$       mtyp => mtyp_storage(1,1,1,0)
       compute = .false.
       release = .false.

    Else
       ! Check to see if space remains
       compute = .true.
       n_elem = nn

       If (n_elem <= MAX_DIAG_STORAGE - n_elem_DIAG) Then
	  release = .false.
	  compute = .true.
	  n_elem_DIAG = n_elem_DIAG + n_elem
       Else
	  release = .true. ! no permanent space available
	  compute = .true.
       End If

       Allocate(diag_ptr_corr(handle)%ptr(n_elem), STAT = ierr)
          ALWAYS_ASSERT(ierr == 0)
          If (mall_ison()) Call mall_mci(diag_ptr_corr(handle)%ptr, myname)
       associated_diag_ptr_corr(handle)=.true.
       ptr  => diag_ptr_corr(handle)%ptr
!!$       mtyp => mtyp_storage(1,1,1,0)

    End If

  End Function get_handle_diag_

  Function get_handle_offdiag_(kti,ktj,kri,krj,kind_cov,ni,nj,mtyp,compute,release) Result (ptr)
    use m_die,only : perr
    use m_mall, only : mall_mci, mall_ison
    use m_Spherical_Partition, Only : NumberOfRegions
    Implicit None

    Integer, Intent(In) :: kti
    Integer, Intent(In) :: ktj
    Integer, Intent(In) :: kri
    Integer, Intent(In) :: krj
    Integer, Intent(In) :: kind_cov
    Integer, Intent(In) :: ni,nj

    Character*1, Pointer :: mtyp
    Logical, Intent(Out) :: compute
    Logical, Intent(Out) :: release
    Real, Pointer :: ptr(:,:)

    Integer :: n_elem, kcov, kr, dkr, ierr
    Character(Len=*), Parameter :: myname_ = myname // '::get_handle_offdiag_'

    select case(kind_cov)
    case (kind_covF)
       kcov = 1
    case (kind_covS)
       kcov = 2
    case (kind_covV)
       kcov = 3
    case default
       call perr(myname_,'unexpected kind_cov',kind_cov)
       ierr=-1
       return
    end select

    kr  = kri
    dkr = mod(krj-kri+NumberOfRegions(),NumberOfRegions())


    release = .true.
    compute = .true.

    ! ptr  => Null()
    nullify(ptr)
!!$    mtyp => mtyp_storage(kti,ktj,kr,dkr)
    Allocate(ptr(ni,nj), STAT = ierr)
    ALWAYS_ASSERT(ierr == 0)
    If (mall_ison()) Call mall_mci(ptr, myname)


!!$    If (Associated(ptr_corr(kti,ktj,kr,dkr,kcov)%ptr)) Then
!!$       ptr  => ptr_corr(kti,ktj,kr,dkr,kcov)%ptr
!!$       mtyp => mtyp_storage(kti,ktj,kr,dkr,kcov)
!!$       compute = .false.
!!$       release = .false.
!!$    Else
!!$       compute = .true.
!!$       ! Check to see if space remains
!!$       n_elem = ni * nj
!!$       If (((dkr == 0) .and. (n_elem <= MAX_REGN_STORAGE - n_elem_REGN)) .or. &
!!$	    & ((dkr > 0) .and. (n_elem <= MAX_GLOB_STORAGE - n_elem_GLOB))) Then
!!$	  release = .false.
!!$	  Allocate(ptr_corr(kti,ktj,kr,dkr,kcov)%ptr(ni,nj), STAT = ierr)
!!$	  ALWAYS_ASSERT(ierr == 0)
!!$	  ptr  => ptr_corr(kti,ktj,kr,dkr,kcov)%ptr
!!$	  mtyp => mtyp_storage(kti,ktj,kr,dkr,kcov)
!!$
!!$	  If (dkr == 0) Then
!!$	     n_elem_REGN = n_elem_REGN - n_elem
!!$	  Else
!!$	     n_elem_GLOB = n_elem_GLOB - n_elem
!!$	  End If
!!$
!!$       Else
!!$
!!$	  release = .true.
!!$	  Allocate(ptr(ni,nj), STAT = ierr)
!!$	  ALWAYS_ASSERT(ierr == 0)
!!$	  mtyp => mtyp_storage(kti,ktj,kr,dkr,kcov)
!!$
!!$       End If
!!$    End If

  End Function get_handle_offdiag_
  
  Subroutine Clean()
    use m_mall, only : mall_ison, mall_mco
    Implicit None
    Integer :: i,j,k

!!$    Do k = 1, n_kind_cov
       Do j = 1, MAX_DIAG_SEGS
	     If (associated_diag_ptr_corr(j)) then
                   If (mall_ison()) Call mall_mco(diag_ptr_corr(j)%ptr, myname)
                Deallocate(diag_ptr_corr(j)%ptr)
                associated_diag_ptr_corr(j)=.false.
	     endif
	  End Do
!!$    End Do
!!$    mtyp_storage  = 'E'


    MAX_DIAG_STORAGE = 0
    MAX_REGN_STORAGE = 0
    MAX_GLOB_STORAGE = 0

    n_elem_DIAG = 0
    n_elem_REGN = 0
    n_elem_GLOB = 0

    excess_storage = TOT_STORAGE

  End Subroutine Clean

End module m_block_storage
