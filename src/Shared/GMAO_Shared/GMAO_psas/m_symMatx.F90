!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_symMatx - blocked matrix attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_symMatx
      use m_AttrVect, only : symMatx => AttrVect, clean
      implicit none
      private	! except

      public :: symMatx			! The class data structure
      public :: symMatx_schedule	! Define a (symMatx)
      public :: symMatx_clean		! Clean a (symMatx)

      public :: symMatx_irow		! index of "irow"
      public :: symMatx_icol		! index of "icol"
      public :: symMatx_ierr		! index of "ierr"
#ifdef _OPENMP
      public :: symMatx_ithd		! index of "OpenMP thread"
#endif
      public :: symMatx_cost		! index of "cost"
      public :: symMatx_nels		! index of "nels"
      public :: symMatx_mcst		! index of "mcst"


      interface symMatx_schedule;module procedure	&
	schedule_
      endinterface

      interface symMatx_clean; module procedure clean_; end interface

! !REVISION HISTORY:
!	24Aug00	- Jing Guo
!		. Restructured m_symMatx, with cxpy_() moved to a newly
!		  from module m_CorMatxF (possiblly similar procedures
!		  in other similar modules too).
!		. Reconciled with Tom Clune's optimization changes.
!       06Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.  This 
!                 change affects only the routine Cxpy_().
! 	15Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_symMatx'

! The symMatx data type is an AttrVect data type (see m_AttrVect for 
! details).  This data type has two storage arrays, one for integer 
! attributes (symMatx%iAttr) and one for real attributes (symMatx%rAttr).
! The set of attributes in each storage array are defined by a List; 
! symMatx%iList for integer attributes, and symMatx%rList for real 
! attributes.  Below are defined the attribute lists, and the indices 
! that reference them:
!
! symMatx_iList components:
!    irow : row index
!    icol : column index
!    ierr : error flag

#ifdef _OPENMP
  character(len=*),parameter :: symMatx_iList='irow:icol:ierr:ithd'
#else
  character(len=*),parameter :: symMatx_iList='irow:icol:ierr'
#endif

  integer,parameter :: symMatx_irow=1
  integer,parameter :: symMatx_icol=2
  integer,parameter :: symMatx_ierr=3
#ifdef _OPENMP
  integer,parameter :: symMatx_ithd=4
#endif

! symMatx_rList components:
!    cost : cost used to partition matrix
!    nels : number of elements in this block (used to compute
!           global number of elements and sparsity)
!    mcst : measured cost to compute/apply block (for 
!           use in cost monitoring and load re-balancing)

  character(len=*),parameter :: symMatx_rList='cost:nels:mcst'
  integer,parameter :: symMatx_cost=1
  integer,parameter :: symMatx_nels=2
  integer,parameter :: symMatx_mcst=3

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

    subroutine schedule_(mNav,spGath, spScat, collnav,  &
	 aNav,krNav,ktNav,              		&
	 comm,kind_C,kind_M,sparse,TotalCost,myCost)

      use m_SparseComm, only : SparseComm
      use m_SparseComm, only : Analyze
      use m_SparseComm, only : DescribeComm
      use m_SparseComm, only : GATHER, REDUCE_SCATTER ! operations

      use m_Collector, only : Collector
      use m_Distribution, only : ptr_coll
      use m_Collector, only : getcoll => get

      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init  => init
      use m_AttrVect,only : AttrVect_lsize => lsize
      use m_AttrVect,only : AttrVect_clean => clean
      use m_AttrVect,only : ptr_iAttr
      use m_AttrVect,only : ptr_rAttr

      use m_Navigator ,only : Navigator
      use m_Navigator, only : lsize, get

      use m_PartMatx,only : PartMatx, PartMatx_SCATTER, ChooseMethod
      use m_SortingTools,only : IndexSet
      use m_SortingTools,only : IndexSort
      use m_die       ,only : die
      use m_mall      ,only : mall_ison,mall_mci,mall_mco

      use m_recurPart, only : recurPart

      use m_mpout, only : mpout, mpout_ison, mpout_flush, mpout_log
      use m_die, only : assert_

      implicit none
      type(symMatx),  intent(out) :: mNav
      type(SparseComm), intent(inout) :: spGath
      type(SparseComm), intent(inout) :: spScat
      type(Collector), intent(in) :: collnav

      type(Navigator),intent(in) :: aNav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav
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

  type(symMatx) :: sMatx
  integer,dimension(:),allocatable :: iMatx
  integer :: mBlox
  integer :: i,l,ier

	! Aliases for mNav and sMatx components

  integer,pointer,dimension(:,:) ::  mNav_iAttr
  real   ,pointer,dimension(:,:) ::  mNav_rAttr
  integer,pointer,dimension(:,:) :: sMatx_iAttr
  real   ,pointer,dimension(:,:) :: sMatx_rAttr

  logical, allocatable :: segin(:), segmid(:)
  Integer :: lc, le, nseg, lcseg, leseg
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

	! Create a (symMatx) of the full matrix

  call schedAll_(sMatx,kind_C,kind_M,sparse,comm,aNav,krNav,ktNav)

!________________________________________

	! Allocate a workspace

  mBlox=AttrVect_lsize(sMatx)

  allocate(iMatx(mBlox),stat=ier)
	if(ier /= 0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(iMatx,myname)
!________________________________________

	! The blocks of sMatx is handled through iMatx.  Therefore,
	! one may arrange mNav in any order through indexed sorting
	! of sMatx.

  call IndexSet(iMatx)

	! Partition sMatx.  iMatx contains the indices to only the
	! local matrix blocks.  Note that the logical size of iMatx
	! (or the value of mBlox) will be reduced at the return.

        method = ChooseMethod('SYM_SCHEDULE',PartMatx_SCATTER)

  call PartMatx(sMatx, iMatx, mBlox, comm,	&
	TotalCost=TotalCost,myCost=myCost,method=method)

		! Private components of sMatx are accessed through
		! aliases.  These aliases will be nullified later.

	sMatx_iAttr => ptr_iAttr(sMatx)
	sMatx_rAttr => ptr_rAttr(sMatx)

	! Incremental post-partitioning sorting for the local set.
	! Note that while iMatx(1:mBlox) contains only a subset of all
	! matrix blocks, their locations in sMatx could be anywhere.
	! Therefore, the whole sMatx storage is passed for the keys.
	! This is from Tom's version.  The purpose was not stated.
	! However, it could be for the cache optimization or to
	! reorder the blocks after SchedAll_() becomes message-passing
	! parallel in a inter-leaving approach.

	! Create a (symMatx) of the local partition of the matrix,
	! by keeping only _local_ blocks.

  call AttrVect_init(mNav,sMatx,mBlox)

		! The private components of mNav are accessed through
		! aliases.  These aliases will be nullified later.

	 mNav_iAttr => ptr_iAttr(mNav)
	 mNav_rAttr => ptr_rAttr(mNav)

  do i=1,mBlox
    l=iMatx(i)
    mNav_iAttr(:,i)=sMatx_iAttr(:,l)
    mNav_rAttr(:,i)=sMatx_rAttr(:,l)
  end do

#ifdef _OPENMP
!   Determine thread assignment in schedule
!   ---------------------------------------
  n_threads = OMP_GET_MAX_THREADS()
  Allocate(ithrdMatx(mBlox,0:n_threads-1), STAT = ier)
     if(ier /= 0) call die(myname_,'allocate()',ier)
     Call mall_mci(ithrdMatx,myname)

  ithrdMatx = Spread(iMatx(1:mBlox),2,n_threads) ! duplicate
  mNav_iAttr(symMatx_ithd,1:mBlox) = -1

!$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC), PRIVATE(i_thread,thrdCost,mthrdBlox), &
!$OMP&  SHARED(mBlox, n_threads,mNav_rAttr,ithrdMatx,mNav_iAttr)
  Do i_thread = 0, n_threads - 1
     thrdCost = Sum(mNav_rAttr(symMatx_cost,1:mBlox))
     mthrdBlox = mBlox
     call IndexSet(ithrdMatx(:,i_thread))
     call IndexSort(mthrdBlox,ithrdMatx(:,i_thread),mNav_rAttr(symMatx_cost,1:mBlox),descend=.true.)
     call RecurPart(mNav_rAttr(symMatx_cost,1:mblox),ithrdMatx(:,i_thread),mthrdBlox,thrdCost, &
          & i_thread, 0, n_threads-1)
     mNav_iAttr(symMatx_ithd,ithrdMatx(1:mthrdBlox,i_thread)) = i_thread
  End Do
!$OMP END PARALLEL DO
     Call mall_mco(ithrdMatx,myname)
  Deallocate(ithrdMatx, STAT = ier)
     if(ier /= 0) call die(myname_,'deallocate()',ier)
  If (any(mNav_iAttr(symMatx_ithd,:) < 0)) Call die(myname_,': Error in thread assignment',ier)
#ifndef NDEBUG  
!!!! #ifdef VERIFY_
  Write(mpout,*)'sym schedule: ',mNav_iAttr
#endif
#endif

#ifndef NDEBUG  

  if (mpout_ison()) Then 
     write(mpout,*)' '
     write(mpout,*)'****  load balance summary ***'
     write(mpout,*)'   kind_M = ',kind_M,'; kind_C = ',kind_C, '; mblox = ',mblox,'; nblox = ',size(iMatx)
     write(mpout,'(a16,f8.4)')'   local cost = ',Sum(mNav_rattr(symMatx_cost,:))
     write(mpout,'(a16,f8.4)')'   total cost = ',Sum(sMatx_rAttr(symMatx_cost,:))
     Write(mpout,'(a16,f8.4)')'   est. spdup = ',Sum(sMatx_rAttr(symMatx_cost,:))/(1.e-10 + Sum(mNav_rattr(symMatx_cost,:)))
     write(mpout,*)' '
     Call flush(mpout)
  End if

#endif

!-------------------------------------------------------------
  ! analyze sparse communication for efficient global operations

  nseg = lsize(aNav)
  Allocate(segin(nseg), segmid(nseg), STAT=ier)
	if(mall_ison()) Then
	   call mall_mci(segin,myname)
	   call mall_mci(segmid,myname)
	endif
  ALWAYS_ASSERT(ier==0)

  segin  = .false. ! default
  segmid = .false. ! default

  Call getcoll(collnav, lbound=lc,ubound=le)
  segin(lc:le) = .true.

#ifndef NDEBUG
  Write(mpout,*)'Begin analysis ',kind_C, kind_M
  Write(mpout,*)'Checking on ob_dstr: (lc,le) = ',lc,le,nseg
  Call flush(mpout)
#endif

  ! Actually these should be left around from before.
  mNav_iAttr => ptr_iAttr(mNav)
  mblox = size(mnav_iattr,2)

  Do i = 1, mblox
     irow = mNav_iattr(symMatx_irow,i)
     icol = mNav_iattr(symMatx_icol,i)
     ASSERT(irow <= nseg)
     ASSERT(icol <= nseg)
     segmid(irow) = .true.
         If ((ktnav(irow) == ktuu .or. ktnav(irow) == ktus) .and. &
	       (kind_c /= kind_covU .and. kind_c /= kind_covC)) segmid(irow+1) = .true. ! u-v symmetry kludge
     segmid(icol) = .true.
        If ((ktnav(icol) == ktuu .or. ktnav(icol) == ktus) .and. &
	       (kind_c /= kind_covU .and. kind_c /= kind_covC)) segmid(icol+1) = .true. ! u-v symmetry kludge

  End Do

! Special handling for zero-sized local vectors (cases where lc>le) is
! made inside subroutine Analyze().  See m_SparseComm_xxx for details.

  Call Analyze(spGath,  segin, segmid, aNav, lc, GATHER, comm)
  Call Analyze(spScat, segmid,  segin, aNav, lc, REDUCE_SCATTER, comm)

	if(mall_ison()) Then
	   call mall_mco(segin,myname)
	   call mall_mco(segmid,myname)
	endif
  Deallocate(segin, segmid)

#ifndef NDEBUG
  Call DescribeComm(spGath,'Gather')
  Call DescribeComm(spScat,'Scatter')
#endif

#ifndef NDEBUG
if(mpout_ison()) then
  write(mpout,*)' '
  write(mpout,*)'****  load balance summary ***'
  write(mpout,*)'   kind_M = ',kind_M,'; kind_C = ',kind_C, '; mblox = ',mblox,'; nblox = ',size(iMatx)
  write(mpout,*)'   local cost = ',Sum(mNav_rattr(symMatx_cost,:)),'; total cost = ',Sum(sMatx_rAttr(symMatx_cost,:)),&
       & '; ratio = ',Sum(sMatx_rAttr(symMatx_cost,:))/(1.e-10 + Sum(mNav_rattr(symMatx_cost,:)))
  write(mpout,*)' '
endif
#endif

!-------------------------------------------------------------
! Clean up working variables

	nullify( mNav_iAttr)	! nullify aliases of mNav
	nullify( mNav_rAttr)

	nullify(sMatx_iAttr)	! nullify aliases of sMatx
	nullify(sMatx_rAttr)



	if(mall_ison()) call mall_mco(iMatx,myname)
  deallocate(iMatx,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate()',ier)

  call AttrVect_clean(sMatx)

end subroutine Schedule_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a symMatx
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(sMat)
      use m_AttrVect,only : AttrVect_clean => clean
      implicit none
      type(symMatx),intent(inout) :: sMat

! !REVISION HISTORY:
! 	01Apr99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call AttrVect_clean(sMat)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: schedAll_ - initialize a full (symMatx) data structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine schedAll_(mNav, kind_cor,kind_mat,sparse,comm,	&
	vNav,krNav,ktNav)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init   => init
      use m_AttrVect, only : resize
      use m_AttrVect, only : clean
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
      use m_die,      only : die,MP_die
      Use m_mpif90,   only : MP_comm_size
      Use m_mpif90,   only : MP_comm_rank
      use m_mpout, only : mpout_log, mpout, mpout_ison

      implicit none

      type(symMatx),intent(out) :: mNav

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
      integer,       intent(in) :: comm

      type(Navigator),intent(in)  :: vNav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

! !REVISION HISTORY:
! 	15Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	19Jul99 - J.W. Larson <jlarson@dao> - extension of symMatx type
!                 to include element counts and measured cost (i.e. actual
!                 load) tracking.  Included interface to simple, tunable 
!                 cost function defined in module m_cost.
!
!EOP ___________________________________________________________________

      include 'ktmax.h'

  integer,parameter :: irow=symMatx_irow
  integer,parameter :: icol=symMatx_icol
  integer,parameter :: ierr=symMatx_ierr
  integer,parameter :: icst=symMatx_cost
  integer,parameter :: inel=symMatx_nels

  character(len=*),parameter :: myname_=myname//'::schedAll_'

!  integer,parameter :: ROOT=0
!
  include "kind_covs.h"

  integer :: mcol,mrow,mband
  integer :: iBlox,nBlox,jb
  integer :: kri,kti,lni
  integer :: krj,ktj,lnj
  integer :: i,j,jm
  real    :: cost, cost_max, cost_tot, cost_n

  integer :: myid, nprocs
  integer :: resid,count,displ
  integer :: ier

  type (symMatx)   :: mNav_tmp
  Type (Collector) :: coll

  integer,pointer,dimension(:,:) :: iAttr
  real   ,pointer,dimension(:,:) :: rAttr

#ifndef NDEBUG
    call mpout_log(myname_,'kind_mat = ',kind_mat)
    call mpout_log(myname_,'kind_cov = ',kind_cor)
#endif

    mrow=lsize(vNav)	! number of row_blocks
    mcol=mrow		! number of colume_blocks (symmetric)
    mband=mcol/2	! number of diagonal bands (0:mband)


    cost_tot = 0
    cost_n   = 0
    cost_max = -1

	! Assign the attributes of each blocks to the (symMatx).
	! Note that for this implementation, the order of the blocks
	! are in a sort of column-major.

    ! This procedure now has several stages to allow parallelization
    ! and the requisite syncronization ...

    ! 1) Create a temporary workspace for the local set with a maximum
    !   estimated size.
    ! 2) Loop through the reduced (local) subset of blocks and compute
    !   sparsity and costs.
    ! 3) all-gather all "local" subsets to get replicated "global"
    !   copies

	! Determine the processor group size, and the local rank:

  call MP_COMM_SIZE( comm, nprocs, ier ) 
	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

  call MP_COMM_RANK( comm, myid, ier )
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  If (kind_mat == 1) call MP_die(myname_,'cannot do univariate here',ier)

	! Determine the location of local workloads

  resid=mod(mrow,nprocs)

  count=mrow/nprocs
  if(myid< resid) count=count+1

  displ=count*myid
  if(myid>=resid) displ=displ+resid

	! Determine a size for the local workspace to avoid two
	! two passes.

  nBlox=count*(mband+1)
!_______________________________________________________________________

	! Create a directory data structure for sub-matrix blocks.
	! This is a temporary attribute vector for local PE only.

  call AttrVect_init(mNav_tmp,iList=symMatx_iList,	&
       rList=symMatx_rList, lsize=nBlox)

       iAttr=>ptr_iAttr(mNav_tmp)
       rAttr=>ptr_rAttr(mNav_tmp)

	! The following loop is now parallelized in a blocking form.
	! This approach can preserve the original block sequence from
	! a serial algorithm.  The earlier parallel algorithm was in
	! an inter-leaving form.

    iBlox = 0
    do i=displ+1,displ+count	! Count only if it is my row.

      call get(vNav,i,ln=lni)
      if(lni <= 0) cycle

	kri=krNav(i)
	kti=ktNav(i)
	If ((kti == ktvs .or. kti == ktvv) .and.       &
	     (kind_cor /= kind_covU .and. kind_cor /= kind_covC)) Cycle ! Exploit u-v symmetry

		! If mrow is even, the last _band_ is half
      jm=mband
      if( mod(mrow,2)==0 .and. i>mband) jm=jm-1

      do jb=0,jm	! check all columns for sparsity

	j=i+jb
	if(j > mcol) j=j-mcol

	call get(vNav,j,ln=lnj)
	if(lnj <= 0) cycle

	krj=krNav(j)
	ktj=ktNav(j)
	If ((ktj == ktvs .or. ktj == ktvv) .and.      &
	     (kind_cor /= kind_covU .and. kind_cor /= kind_covC)) Cycle ! Exploit u-v symmetry

	if(sparse(kind_mat,kind_cor, kri,kti, krj,ktj)) cycle

			! This block seems count

	iBlox = iBlox + 1

	iAttr(irow,iBlox)=i
	iAttr(icol,iBlox)=j
	iAttr(ierr,iBlox)=0

#ifdef SIZECOSTF

	cost=lni*lnj
	if(i==j) cost=lni*(lni+1)/2

#else

	! A fancy cost function may be built and used here
	! as this:
	!
	!   = cost(kri,kti,lni, krj, ktj,lnj, diag=(i==j) )

!jwl this cost function returns the cost for the _entire matrix block_:

        cost = mb_cost( kind_mat, kind_cor, kri,kti,lni, krj,ktj,lnj )
	if (cost > cost_max) cost_max = cost
	cost_tot = cost_tot + cost
	cost_n = cost_n + 1
#endif

	rAttr(icst,iBlox)=0.
	if(cost>0.) rAttr(icst,iBlox)=cost

 		! Floating point numbers are used to compute the
 		! product, in case the product is too big for integer.

 	rAttr(inel,iBlox)=real(lni)*real(lnj)

      end do
    end do

#ifndef SIZECOSTF
#ifndef NDEBUG
  if(mpout_ison()) then
    write(mpout,*)' load balance info? ', cost_n
    write(mpout,*)'    max: ',cost_max
    write(mpout,*)'    avg: ',cost_tot/max(1.,cost_n)
  endif
#endif
#endif


	nullify(iAttr)
	nullify(rAttr)

	
	! Fix the size

    call resize(mNav_tmp,lsize=iBlox)
!_______________________________________________________________________

	! Create a directory data structure for sub-matrix blocks.
	! This is the real one with a real size

    call Collector_init(coll,iBlox,comm)

    call allgatherv(mNav_tmp,mNav,coll,comm,stat=ier)
	if(ier/=0) call die(myname_,'allgatherv()',ier)

    call clean(coll)
    call clean(mNav_tmp)

end subroutine schedAll_

end module m_symMatx
!.
