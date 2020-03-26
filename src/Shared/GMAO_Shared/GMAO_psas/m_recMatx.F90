!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_recMatx - blocked matrix attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_recMatx
      use m_AttrVect, only : recMatx => AttrVect, clean

      implicit none
      private	! except

      public :: recMatx			! The class data structure
      public :: recMatx_schedule	! define a (recMatx)
      public :: recMatx_clean		! clean a (recMatx)

      public :: recMatx_irow		! index of "irow"
      public :: recMatx_icol		! index of "icol"
      public :: recMatx_ierr		! index of "ierr"
#ifdef _OPENMP
    public :: recMatx_ithd		! index of "ithd"
#endif
      public :: recMatx_cost		! index of "cost"
      public :: recMatx_nels		! index of "nels"
      public :: recMatx_mcst		! index of "mcst"

      interface recMatx_schedule;module procedure	&
	schedule_
      endinterface

      interface recMatx_clean; module procedure clean_; end interface

! !REVISION HISTORY:
! 	15Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	21Jul99 - J.W. Larson <jlarson@dao> - extension of recMatx type
!                 to include element counts and measured cost (i.e. actual
!                 load) tracking.  Included interface to simple, tunable 
!                 cost function defined in module m_cost.
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.  This change
!                 affects only the routine Cxpy_().
!       29Oct01 - T. Clune <clune@sgi.com> - modified to use
!                 Sparse communication patterns.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_recMatx'

! The recMatx data type is an AttrVect data type (see m_AttrVect for 
! details).  This data type has two storage arrays, one for integer 
! attributes (recMatx%iAttr) and one for real attributes (recMatx%rAttr).
! The set of attributes in each storage array are defined by a List; 
! recMatx%iList for integer attributes, and recMatx%rList for real 
! attributes.  Below are defined the attribute lists, and the indices 
! that reference them:
!
! recMatx_iList components:
!    irow : row index
!    icol : column index
!    ierr : error flag
#ifdef _OPENMP
  character(len=*),parameter :: recMatx_iList='irow:icol:ierr:ithd'
#else
  character(len=*),parameter :: recMatx_iList='irow:icol:ierr'
#endif
  integer,parameter :: recMatx_irow=1
  integer,parameter :: recMatx_icol=2
  integer,parameter :: recMatx_ierr=3
#ifdef _OPENMP
  integer,parameter :: recMatx_ithd=4
#endif

! recMatx_rList components:
!    cost : cost used to partition matrix
!    nels : number of elements in this block (used to compute
!           global number of elements and sparsity)
!    mcst : measured cost to compute/apply block (for 
!           use in cost monitoring and load re-balancing)
  character(len=*),parameter :: recMatx_rList='cost:nels:mcst'
  integer,parameter :: recMatx_cost=1
  integer,parameter :: recMatx_nels=2
  integer,parameter :: recMatx_mcst=3

#include "assert.H"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: schedule_ - Schedule the works for the local PE.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine schedule_(mNav, spGath, spScat,		&
        row_collNav, col_collNav,       &
	row_Nav,row_krNav,row_ktNav,	&
	col_Nav,col_krNav,col_ktNav,	&
	comm, kind_C,kind_M,sparse,TotalCost,myCost)

      use m_AttrVect ,only : AttrVect_init => init
      use m_AttrVect ,only : AttrVect_clean=>clean
      use m_AttrVect ,only : AttrVect_lsize=>lsize
      use m_AttrVect ,only : ptr_iAttr
      use m_AttrVect ,only : ptr_rAttr

      use m_Collector,only : Collector
      use m_Collector, only : getcoll => get

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
!!$      use m_PartMatx,only : PartMatx, PartMatx_RECTANG, ChooseMethod
      use m_PartMatx,only : PartMatx, PartMatx_SCATTER, ChooseMethod
      use m_SortingTools,only : IndexSet
      use m_SortingTools,only : IndexSort
      use m_recurPart, only : recurPart

      use m_die, only : perr_die
      use m_die,only  : die
      use m_die, only : assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      use m_SparseComm, only : SparseComm
      use m_SparseComm, only : Analyze
      use m_SparseComm, only : DescribeComm
      use m_SparseComm, only : GATHER, REDUCE_SCATTER ! operations

      use m_mpout, only : mpout_log, mpout

      implicit none
      type(recMatx), intent(out):: mNav
      type(SparseComm), intent(inout) :: spGath
      type(SparseComm), intent(inout) :: spScat
      type(Collector), intent(in) :: row_CollNav
      type(Collector), intent(in) :: col_CollNav
	
      type(Navigator),intent(in) :: row_Nav
	integer,dimension(:),intent(in) :: row_krNav
	integer,dimension(:),intent(in) :: row_ktNav
      type(Navigator),intent(in) :: col_Nav
	integer,dimension(:),intent(in) :: col_krNav
	integer,dimension(:),intent(in) :: col_ktNav
      integer,intent(in) :: comm
      integer,intent(in) :: kind_C	! kind of correlation
      integer,intent(in) :: kind_M	! kind of the matrix
      Interface 
         logical function sparse(kind_mat,kind_cov, kr_i,kt_i, kr_j,kt_j)
           use m_Spherical_Partition, only : Base80Region
           integer, intent(in)	:: kind_mat	! which matrix
           integer, intent(in)	:: kind_cov	! which covariance
           integer, intent(in)	:: kr_i,kt_i	! row block indices
           integer, intent(in)	:: kr_j,kt_j	! column block indices
         End function sparse
      End Interface
      real,optional,intent(out) :: TotalCost
      real,optional,intent(out) :: myCost

! !REVISION HISTORY:
! 	15Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::schedule_'

  type(recMatx) :: rMatx
  integer,dimension(:),allocatable :: iMatx
  integer :: mBlox
  integer :: i,l,ier

  integer,pointer,dimension(:,:) ::  mNav_iAttr
  real   ,pointer,dimension(:,:) ::  mNav_rAttr
  integer,pointer,dimension(:,:) :: rMatx_iAttr
  real   ,pointer,dimension(:,:) :: rMatx_rAttr

! variables used by Sparse communication analysis
  logical, allocatable :: colseg_in(:), colseg_out(:)
  logical, allocatable :: rowseg_in(:), rowseg_out(:)
  Integer :: lc_row, le_row, lc_col, le_col, nseg_row, nseg_col, lcseg, leseg
  Integer :: irow, icol, iseg
  Integer :: method

#ifdef _OPENMP
  Integer :: i_thread, n_threads, mthrdBlox
  Integer, Allocatable :: ithrdMatx(:,:)
  Real    :: thrdCost
  Integer, External :: OMP_GET_MAX_THREADS
#endif
      include 'ktmax.h'
      include "kind_covs.h"

	! Create a (recMatx) for the full matrix

  call schedAll_(rMatx,kind_C,kind_M,sparse,comm,	&
	row_Nav,row_krNav,row_ktNav,			&
	col_Nav,col_krNav,col_ktNav)

	! Allocate a workspace

  mBlox=AttrVect_lsize(rMatx)
	call mpout_log(myname_,'mBlox = ',mBlox)

  allocate(iMatx(mBlox),stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(iMatx,myname)

	! The blocks of rMatx is handled through iMatx.  Therefore,
	! one may arrange mNav in any order through indexed sorting
	! of rMatx.

  call IndexSet(iMatx)

	! Partition rMatx.  iMatx contains the indices to only the
	! local matrix blocks.

        method = ChooseMethod('REC_SCHEDULE',PartMatx_SCATTER)
  call PartMatx(rMatx, iMatx,mBlox, comm,	&
	TotalCost=TotalCost,myCost=myCost,method=method)

	! Create a (recMatx) for the local partition of the matrix,
	! by keeping only the blocks handled _locally_.

  call AttrVect_init(mNav,rMatx,mBlox)

	 mNav_iAttr => ptr_iAttr(mNav)
	 mNav_rAttr => ptr_rAttr(mNav)
	rMatx_iAttr => ptr_iAttr(rMatx)
	rMatx_rAttr => ptr_rAttr(rMatx)

  do i=1,mBlox
    l=iMatx(i)
    mNav_iAttr(:,i)=rMatx_iAttr(:,l)
    mNav_rAttr(:,i)=rMatx_rAttr(:,l)
  end do
	
#ifdef _OPENMP
!   Determine thread assignment in schedule
!   ---------------------------------------

  n_threads = OMP_GET_MAX_THREADS()
  Allocate(ithrdMatx(mBlox,0:n_threads-1), STAT = ier)
     if(ier /= 0) call die(myname_,'allocate()',ier)
     Call mall_mci(ithrdMatx,myname)

  ithrdMatx = Spread(iMatx(1:mBlox),2,n_threads) ! duplicate
  mNav_iAttr(recMatx_ithd,1:mBlox) = -1
!$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC), PRIVATE(i_thread,thrdCost,mthrdBlox), &
!$OMP&  SHARED(mblox, myCost, n_threads, mNav_rAttr, ithrdMatx, mNav_iAttr)
  Do i_thread = 0, n_threads - 1
     thrdCost = Sum(mNav_rAttr(recMatx_cost,1:mBlox))
     mthrdBlox = mBlox
     call IndexSet(ithrdMatx(:,i_thread))
     call IndexSort(mthrdBlox,ithrdMatx(:,i_thread),mNav_rAttr(recMatx_cost,1:mBlox),descend=.true.)
     mthrdBlox = mBlox
     call RecurPart(mNav_rAttr(recMatx_cost,1:mblox),ithrdMatx(:,i_thread),mthrdBlox,thrdCost, &
          & i_thread, 0, n_threads-1)
     mNav_iAttr(recMatx_ithd,ithrdMatx(1:mthrdBlox,i_thread)) = i_thread
  End Do
!$OMP END PARALLEL DO

     Call mall_mco(ithrdMatx,myname)
  Deallocate(ithrdMatx, STAT = ier)
     if(ier /= 0) call die(myname_,'deallocate()',ier)
  If (any(mNav_iAttr(recMatx_ithd,:) < 0)) Call die(myname_,': Error in thread assignment',ier)
#ifdef VERIFY_
  Write(mpout,*)'rec schedule: ',mNav_iAttr
#endif
#endif

	nullify(rMatx_iAttr)
	nullify(rMatx_rAttr)
	nullify( mNav_iAttr)
	nullify( mNav_rAttr)

	! Clean up working variables

	if(mall_ison()) call mall_mco(iMatx,myname)
  deallocate(iMatx,stat=ier)
	if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

  call AttrVect_clean(rMatx)

!-------------------------------------------------------------
  ! analyze sparse communication for efficient global operations

  nseg_col = lsize(col_Nav)
  nseg_row = lsize(row_Nav)
  Allocate(colseg_in(nseg_col), colseg_out(nseg_col), &
           rowseg_in(nseg_row), rowseg_out(nseg_row), STAT=ier)
        ALWAYS_ASSERT(ier==0)
	if(mall_ison()) Then
	   call mall_mci(colseg_in,myname)
	   call mall_mci(colseg_out,myname)
	   call mall_mci(rowseg_in,myname)
	   call mall_mci(rowseg_out,myname)
	endif

  colseg_in  = .false.
  colseg_out = .false.
  rowseg_in  = .false.
  rowseg_out = .false.

  Call getcoll(col_CollNav, lbound=lc_col,ubound=le_col)
  colseg_in(lc_col:le_col) = .true.

  Call getcoll(row_CollNav, lbound=lc_row,ubound=le_row)
  rowseg_out(lc_row:le_row) = .true.

  ! Actually these should be left around from before.
  mNav_iAttr => ptr_iAttr(mNav)
  mblox = size(mNav_iAttr,2)

  Do i = 1, mblox
     irow = mNav_iattr(recMatx_irow,i)
     icol = mNav_iattr(recMatx_icol,i)
     ASSERT(irow <= nseg_row)
     ASSERT(icol <= nseg_col)
     colseg_out(icol) = .true.
     rowseg_in(irow)  = .true.
         If ((row_ktnav(irow) == ktuu .or. row_ktnav(irow) == ktus) .and. &
	       (kind_c /= kind_covU .and. kind_c /= kind_covC)) rowseg_in(irow+1) = .true. ! u-v symmetry kludge
        If ((col_ktnav(icol) == ktuu .or. col_ktnav(icol) == ktus) .and. &
	       (kind_c /= kind_covU .and. kind_c /= kind_covC)) colseg_out(icol+1) = .true. ! u-v symmetry kludge

  End Do

! Special handling for zero-sized local vectors (cases where
! lc_row > le_row and lc_col > le_col) is made inside subroutine
! Analyze().  See m_SparseComm_xxx for details.

  Call Analyze(spGath, colseg_in, colseg_out, col_Nav, lc_col, GATHER, comm)
  Call Analyze(spScat, rowseg_in, rowseg_out, row_Nav, lc_row, REDUCE_SCATTER, comm)

	if(mall_ison()) Then
	   call mall_mco(colseg_in,myname)
	   call mall_mco(colseg_out,myname)
	   call mall_mco(rowseg_in,myname)
	   call mall_mco(rowseg_out,myname)
	endif
  Deallocate(colseg_in, colseg_out, rowseg_in, rowseg_out)

#ifndef NDEBUG
  Call DescribeComm(spGath,'Gather')
  Call DescribeComm(spScat,'Scatter')
#endif


end subroutine Schedule_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a recMatx
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(rMat)
      use m_AttrVect,only : AttrVect_clean => clean
      implicit none
      type(recMatx),intent(inout) :: rMat

! !REVISION HISTORY:
! 	01Apr99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call AttrVect_clean(rMat)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: schedAll_ - initialize a full (recMatx) data structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine schedAll_(mNav, kind_cor,kind_mat,sparse,comm,	&
	row_Nav,row_krNav,row_ktNav,	&
	col_Nav,col_krNav,col_ktNav)

      use m_AttrVect, only : AttrVect_init
      use m_AttrVect, only : clean
      use m_AttrVect, only : resize
      use m_AttrVect, only : ptr_iAttr
      use m_AttrVect, only : ptr_rAttr

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get

      use m_Collector,only : Collector
      use m_Collector,only : Collector_init
      use m_Collector,only : clean
      use m_CollectorComm,only : allgatherv

      use m_costs,    only : mb_cost
      use m_die,      only : die
      use m_die,      only : MP_die

      use m_mpif90,   only : MP_comm_size
      use m_mpif90,   only : MP_comm_rank

      implicit none

      type(recMatx),intent(out) :: mNav

      integer,      intent(in)  :: kind_cor
      integer,      intent(in)  :: kind_mat
      Interface 
         logical function sparse(kind_mat,kind_cov, kr_i,kt_i, kr_j,kt_j)
           use m_Spherical_Partition, only : Base80Region
           integer, intent(in)	:: kind_mat	! which matrix
           integer, intent(in)	:: kind_cov	! which covariance
           integer, intent(in)	:: kr_i,kt_i	! row block indices
           integer, intent(in)	:: kr_j,kt_j	! column block indices
         End function sparse
      End Interface
      integer,      intent(in)  :: comm

      type(Navigator),intent(in) :: row_Nav
      integer,dimension(:),intent(in) :: row_krNav
      integer,dimension(:),intent(in) :: row_ktNav

      type(Navigator),intent(in) :: col_Nav
      integer,dimension(:),intent(in) :: col_krNav
      integer,dimension(:),intent(in) :: col_ktNav

! !REVISION HISTORY:
!	23Aug00	- Jing Guo
!		. Implemented message-passing calculations for the
!		  sparsity tests.  See m_symMatx for Tom Clune's
!		  initial algorithm.
!		. Simplified out the first loop.  It was implemented
!		  to count the exact number of matrix blocks.  An
!		  estimated local size is used now with a new resize()
!		  functionality of m_AttrVect.
!
! 	15Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::schedAll_'

  include 'ktmax.h'

  integer,parameter :: irow=recMatx_irow
  integer,parameter :: icol=recMatx_icol
  integer,parameter :: ierr=recMatx_ierr
  integer,parameter :: icst=recMatx_cost
  integer,parameter :: inel=recMatx_nels

  integer :: mcol,mrow
  integer :: iBlox,nBlox
  integer :: kri,kti,lni
  integer :: krj,ktj,lnj
  integer :: i,j
  real    :: cost
  integer :: nPEs,myID
  integer :: ier
  integer :: resid,count,displ

  type(recMatx)   :: tNav
  type(Collector) :: coll

  integer,pointer,dimension(:,:) :: iAttr
  real   ,pointer,dimension(:,:) :: rAttr

    mrow=lsize(row_Nav)	! number of row_blocks
    mcol=lsize(col_Nav)	! number of col_blocks

	! Estimate a size.  Make it roomy

    call MP_comm_size(comm,nPEs,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
    call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

    resid=mod(mrow,nPEs)
    count=    mrow/nPEs
    if(myID< resid) count=count+1
    displ=count*myID
    if(myID>=resid) displ=displ+resid

    nBlox=count*mcol

	! Create a directory data structure for sub-matrix blocks.

    call AttrVect_init(tNav,iList=recMatx_iList,	&
			    rList=recMatx_rList, lsize=nBlox	)

	! Need to accessing private data of tNav

    iAttr=>ptr_iAttr(tNav)
    rAttr=>ptr_rAttr(tNav)

	! Assign the attributes of each blocks to the (recMatx).
	! Note that for this implementation, the order of the blocks
	! are sort of column-major.

    iBlox=0
    do i=displ+1,displ+count	! only if it is my row.

      call get(row_Nav,i,ln=lni)
      if(lni <= 0) cycle

	kri=row_krNav(i)
	kti=row_ktNav(i)
	If (kti == ktvs .or.	&
	    kti == ktvv) Cycle ! Exploit u-v symmetry


      do j=1,mcol

        call get(col_Nav,j,ln=lnj)
	if(lnj <= 0) cycle

		krj=col_krNav(j)
		ktj=col_ktNav(j)
		If (ktj == ktvs .or.	&
		    ktj == ktvv) Cycle ! Exploit u-v symmetry
	
	if(sparse(kind_mat,kind_cor, kri,kti, krj,ktj)) cycle

			! it seems count
	iBlox=iBlox+1

	iAttr(irow,iBlox)=i
	iAttr(icol,iBlox)=j
	iAttr(ierr,iBlox)=0

!jwl Simple cost function for rectangular matrix block found in m_costs.
!jwl This cost function returns the cost for the _entire matrix block_:

        cost = mb_cost( kind_cor, kri,kti,lni, krj,ktj,lnj )

		! Remove any negative cost values

	rAttr(icst,iBlox) = 0.
	if(cost>0.) rAttr(icst,iBlox) = cost

		! Compute the total number of elements in a matrix
		! block.  Floating point numbers are used in case that
		! the magnitude of the product becomes too large.

	rAttr(inel,iBlox) = real(lni)*real(lnj)

      end do
    end do

			! Nullify temporary pointers to protect the
			! innocents.
	nullify(iAttr)
	nullify(rAttr)

		! Set it to it's real size

    call resize(tNav,lsize=iBlox)

		! Populate the information

    call Collector_init(coll,iBlox,comm)

    call allgatherv(tNav,mNav, coll,comm,stat=ier)
	if(ier/=0) call die(myname_,'allgatherv()',ier)

    call clean(coll)
    call clean(tNav)

end subroutine schedAll_

end module m_recMatx
!.

