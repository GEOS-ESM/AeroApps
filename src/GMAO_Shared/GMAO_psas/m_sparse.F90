Module m_sparse
  use config, only : ktmax,ktHH,ktslp,ktqq
  use config, only : kind_covU,kind_covC,kind_covF,kind_covS,kind_covV
  use config, only : kind_Umat,kind_Rmat,kind_3mat,kind_4mat,kind_5mat
  implicit none
  private

  public :: sparse
  public :: sparse_init
  public :: sparse_clean

  Interface sparse_init
     module procedure init_sp_
     module procedure init_gp_
  end Interface
  Interface sparse_clean
     module procedure clean_
  end Interface

! Module data

  logical,save :: initialized_ = .false.
  logical,save :: save_ = .false.

! Precomputed, to save time (see init_) sparse mask.

  integer,parameter :: NBPW=4	! assigned number of bytes/word
  Logical(Kind=1),Save,Allocatable :: sparse_mask(:,:)
  integer,Save,Allocatable,dimension(:) :: kr_base80
  integer,save :: msize_=4		! mask size in words

  character(len=*),parameter :: myname='m_sparse::'

  include "kttabl.h"

Contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sparse - determine sparsity of a block matrix
!
! !DESCRIPTION:
!	sparse returns sparsity between kr/kt regions using a pre-
!   calculated sparse_mask.  Note this should be interchangeable with
!   the ealier function sparse() which didn't precalculate sepang
!
! !INTERFACE:

logical function sparse(kind_mat,kind_cov, kr_i,kt_i, kr_j,kt_j)

  use m_die , only : die,perr
  implicit none

  integer, intent(in)	:: kind_mat	! which matrix
  integer, intent(in)	:: kind_cov	! which covariance
  integer, intent(in)	:: kr_i,kt_i	! row block indices
  integer, intent(in)	:: kr_j,kt_j	! column block indices

! !REVISION HISTORY:
!	30Jan02	- Jing Guo
!		. Changed sepang_mask(:,:,:) to sparse_mask(:,:) and
!		  many other implementations.
!        1Dec00 - T. Clune and P. Lyster modified to use Clune's sepang_mask rather than cossepang
!       23Oct00 - P. Lyster modified to allow for precomputed cossepang_matrix
! 	20Feb96 - J. Guo	- (to do)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'sparse'

  logical	:: sparse_kr, sparse_kt

	! bitwise-or to create a bit mask

  integer, parameter :: kind_covPf=	&
	ior(ior(kind_covF,kind_covS),kind_covV)

  integer  :: i,j,n
!-----------------------------------------------------------------------
  sparse=.true.		! default value if not a known matrix
!-----------------------------------------------------------------------
	! A sparse matrix block is:
	!
	! 	( a sparse "uncorrelated" observation
	!	  error covariance matrix block		)	.AND.
	! 	( a sparse "correlated" observation
	!	  error covariance matrix block		)	.AND.
	! 	( a sparse "decoupled" wind forecast
	!	  error covariance matrix block		)	.AND.
	! 	( a sparse "mass-wind balanced" forecast
	!	  error covariance matrix block		)

	! In details:
	! 	( a sparse "uncorrelated" observation
	!	  error covariance matrix block		)	.AND.

  if(sparse .and. (iand(kind_cov,kind_covU) /= 0) ) then

    sparse = kt_i.ne.kt_j			! sparse_kt() .or.
    						! sparse_kr()
    if(.not.sparse) sparse = kr_Base80(kr_i).ne.kr_Base80(kr_j)

  endif	! kind_covU

!-----------------------------------------------------------------------
	! 	( a sparse "correlated" observation
	!	  error covariance matrix block		)	.AND.

  if(sparse .and. (iand(kind_cov,kind_covC) /= 0) ) then

    sparse = kt_i .ne. kt_j			! sparse_kt() .or.
    if(.not.sparse) sparse=kt_i/=ktHH			! .or.
    if(.not.sparse) sparse=sparse_kr_(kr_i,kr_j)	! sparse_kr()

  endif	! kind_covC

!-----------------------------------------------------------------------
	! 	( a sparse "decoupled" wind forecast
	!	  error covariance matrix block		)	.AND.

  if(sparse .and. (iand(kind_cov,ior(kind_covS,kind_covV)) /= 0) ) then

    sparse = sparse_kr_(kr_i,kr_j)		! sparse_kr() .or.
							! sparse_kt()
    if(.not. sparse) then
      sparse = .not.ktmvar(kt_i,kt_j)	! sparse_kt() -- by user .or.
						! .not. wind-wind
      if(.not. sparse) then
	sparse = kt_i.eq.ktHH.or.kt_i.eq.ktslp.or.kt_i.eq.ktqq .or. &
		 kt_j.eq.ktHH.or.kt_j.eq.ktslp.or.kt_j.eq.ktqq
      endif
    endif

  endif

!-----------------------------------------------------------------------
	! 	( a sparse "mass-wind balanced" forecast
	!	  error covariance matrix block		)

  if(sparse .and. (iand(kind_cov,kind_covF) /= 0) ) then

    sparse = sparse_kr_(kr_i,kr_j)		! sparse_kr() .or.
							! sparse_kt()
    if(.not. sparse) then
      sparse= .not.ktmvar(kt_i,kt_j)	! sparse_kt() -- by user .or.
						! .not. "univariate"

		! Note: "univariate" is considered only for kind_Umat
		! and winds are considered as the univariate.

      if(.not.sparse .and. kind_mat.eq.kind_Umat) then
	sparse = kt_i.ne.kt_j .and.				&
	  ( kt_i.eq.ktHH.or.kt_i.eq.ktslp.or.kt_i.eq.ktqq.or.	&
	    kt_j.eq.ktHH.or.kt_j.eq.ktslp.or.kt_j.eq.ktqq)
      endif

    endif

  endif

!=======================================================================

contains

!-----------------------------------------------------------------------
  logical function sparse_kr_(kr_i,kr_j)
    implicit none
    integer, intent(in) :: kr_i,kr_j

      sparse_kr_ = .true.
      if(kind_mat.eq.kind_Umat.or.kind_mat.eq.kind_Rmat) then

        sparse_kr_ = kr_Base80(kr_i) .ne. kr_Base80(kr_j)

      else

! Sparse_mask tells if the regions are .not.correlated.  sparse_mask is
! .false. if sepang(kr_i,kr_j) .lt. seplim(kind_mat), where sepang is
! the minimum angle distance between regions and seplim is the cutoff
! that's read in from the resource file.

        sparse_kr_ = kr_i .ne. kr_j .and. sparse_mask(kr_i,kr_j)

      endif

  end function sparse_kr_

end function sparse
!=======================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_sp_ - initialize sparse_mask matrix.
!
! !DESCRIPTION:
!
!	Sparse_mask is .true. if two regions are .not.correlated, or if
! sepang(kr_i,kr_j) .lt. seplim(kind_mat), where sepang is the minimum
! angle distance between regions and seplim is the cutoff that's read in
! from the resource file.
!
! !INTERFACE:

 subroutine init_sp_(SphPart,comm,root)
  use m_Spherical_Partition, only : Spherical_Partition
  use m_Spherical_Partition, only : NumberOfRegions, BuildMask
  use m_Spherical_Partition, only : Base80Region
  use m_redwin, only : redwin_initialized
  use m_redwin, only : redwin_seplim
  use m_mall , only : mall_ison, mall_ci, mall_mci
  use m_mpout, only : mpout,mpout_log
  use m_die  , only : die
  Implicit None
  type(Spherical_Partition),intent(in) :: SphPart
  integer,intent(in) :: comm	! to be used in future.
  integer,intent(in) :: root	! to be used in future.

! !REVISION HISTORY:
!       19Oct00 - Peter Lyster <lys@dao.gsfc.nasa.gov> init sparse matrix
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::init_sp_'
  
   integer :: ier=0
   integer :: nReg
   integer :: i,j,n

! initialize the hardwired masking matrix.  Mostly
! this is done here for optimization purposes (i.e., so you only have to
! calculate sparse_mask once per run).
!
! Calculate sparse_mask once per run may be fine with specific
! applications, but problematic in a general case, since a Spherical_
! Partition might be changed in the same run.  Therefore, a mechenism
! (save_) is added to allow initialization only once; and force_ is
! added to allow forced clean() if needed.  (J.G.)

   If (initialized_) then
     If (.not.save_) call die(myname_,'multiple object definition')
     Return
   endif

   If (.not.redwin_initialized())	&
   	call die(myname_,'redwin not initialized')

	! This module will be initialized once if the same
	! GlobalPartition is set to save_.

   initialized_ = .true.

   nReg=NumberOfRegions(partition=SphPart)

   allocate(sparse_mask(nReg,nReg),kr_base80(nReg),stat=ier)
      If(ier /= 0) call die(myname_,'allocate(sparse_mask)',ier)

      msize_=(size(sparse_mask)+NBPW-1)/NBPW

      If (mall_ison()) then
	Call mall_ci(msize_,myname) ! 4 logicals/word
	Call mall_mci(kr_base80,myname)
      endif
 
! sparse_mask is true if sepang(kr_i,kr_j) .lt. seplim(kind_mat), where
! sepang is the minimum angle distance between regions and seplim is the
! cutoff that's read in from the resource file.

   Call BuildMask(redwin_seplim,sparse_mask(:,:),partition=SphPart)

   sparse_mask(:,:)=.not.sparse_mask(:,:)

   do i=1,nReg
     kr_base80(i)=Base80Region(i,partition=SphPart)
   end do
   
   return
  end subroutine init_sp_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_gp_ - initialize sparse_mask with a GlobalPartition.
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine init_gp_(gpart,comm,root)
    use m_GlobalPartition,only : GlobalPartition
    use m_GlobalPartition,only : ptr_Partition
    Implicit None

    type(GlobalPartition),intent(in) :: gpart
    integer,intent(in) :: comm	! to be used in future.
    integer,intent(in) :: root	! to be used in future.

! !REVISION HISTORY:
!	20Aug02	- Jing Guo
!		. Modified.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_gp_'
  
  call init_sp_(ptr_Partition(gpart),comm,root)
end subroutine init_gp_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean up sparse_mask (don't do if you still
!                         intend to use sparse)
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine clean_(force)
   use m_mall, only : mall_ison, mall_co, mall_mco
   use m_die,  only : die
   implicit none
   logical,optional,intent(in) :: force

! !REVISION HISTORY:
!       19Oct00 - Peter Lyster <lys@dao.gsfc.nasa.gov> init sparse matrix
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::clean_'
  logical :: force_
  integer :: ier

  if(.not.initialized_) call die(myname_,'object undefined')

	! The object will not be cleaned, if it is initialized to
	! "save", unless it is asked to be "force"d to clean().

  if(save_) then
    force_=.false.
    if(present(force)) force_=force
    if(.not.force_) return
  endif

	! If the object is not to be saved, clean all of its
	! components.

      If (mall_ison()) then
	Call mall_co(msize_,myname) ! 4 logicals/word
	Call mall_mco(kr_base80,myname)
      endif
   deallocate(sparse_mask,kr_base80,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  initialized_=.false.
  save_=.false.
  msize_=0
   return

  end subroutine clean_

end Module m_sparse
