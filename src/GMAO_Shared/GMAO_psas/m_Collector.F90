!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Collector - a size description of a distributed array
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Collector
      implicit none
      private	! except

      public :: Collector		! The class data structure

      public :: Collector_init ,init	! initialize a Collector
      public :: Collector_clean,clean	! clean (finish) a Collector

      public :: new		! create an object (for a pointer)
      public :: delete		! delete an object (for a pointer)

      public :: localSize	! find the local size
      public :: globalSize	! find the global size
      public :: getRange	! get the range of a segment
      public :: get		! similar to getRange()
      public :: getRank		! locate a given element

		! Pointer functions below are dangeous!  They should
		! be used for information output only.

      public :: ptr_counts	! refering to %counts
      public :: ptr_displs	! refering to %displs

    type Collector
      private
      integer :: myPE
      integer :: gsize				! the Global size
      integer :: lsize				! my local size
      integer :: lleft				! my left index
      integer,pointer,dimension(:) :: counts	! all local sizes
      integer,pointer,dimension(:) :: displs	! PE ordered locations
    end type Collector

    interface Collector_init ; module procedure	&
	initd_,	&	! initialize from all PEs (allgather(counts))
	initr_, &	! initialize from the root (bcast(counts))
	inita_		! every PE initialize their own (no comm.)
    end interface
    interface init ; module procedure	&
	initd_,	&	! initialize from all PEs (allgather(counts))
	initr_, &	! initialize from the root (bcast(counts))
	inita_		! every PE initialize their own (no comm.)
    end interface

    interface Collector_clean; module procedure clean_; end interface
    interface clean; module procedure clean_; end interface

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

    interface globalSize; module procedure gsize_; end interface

    interface localsize; module procedure lsize_; end interface

    interface getRange; module procedure	&
	range_		! get the segment range of the local or a given
    end interface	! PE

    interface get; module procedure get_; end interface
    interface getRank; module procedure rank_ ; end interface

    interface ptr_counts; module procedure ptr_counts_; end interface
    interface ptr_displs; module procedure ptr_displs_; end interface

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Collector'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - define the map from distributed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initd_(GMap,ln,comm)
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Collector),intent(out) :: GMap
      integer,intent(in) :: ln	! the local size
      integer,intent(in) :: comm

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initd_'
  integer :: nPEs,myID,ier,l,i
  integer :: mesgType

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(GMap%counts,myname)
	  call mall_mci(GMap%displs,myname)
	endif

  mesgType=MP_type(ln)

  call MPI_allgather(ln,1,mesgType,GMap%counts,1,mesgType,comm,ier)
  if(ier/=0) call MP_die(myname_,'MPI_allgather()',ier)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%myPE =myID
  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ initialize the map from the root
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initr_(GMap,lns,root,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_comm_size,MP_comm_rank
      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Collector),intent(out) :: GMap
      integer,dimension(:),intent(in) :: lns	! the distributed sizes
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	29May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: nPEs,myID,ier,l,i
  integer :: mesgType

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(GMap%counts,myname)
	  call mall_mci(GMap%displs,myname)
	endif

  if(myID == root) then
    if(size(lns(:)) /= nPEs) call die(myname_,	&
	'size(lns)',size(lns),'nPEs',nPEs)

    GMap%counts(:)=lns(:)
  endif

  mesgType=MP_type(GMap%counts)

  call MPI_bcast(GMap%counts,nPEs,mesgType,root,comm,ier)
  if(ier/=0) call MP_die(myname_,'MPI_bcast()',ier)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%myPE =myID
  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: inita_ - define the map from known counts
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine inita_(GMap,lns,comm)
      use m_mpif90,only : MP_comm_size,MP_comm_rank
      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Collector),    intent(out) :: GMap
      integer,dimension(0:),intent(in)  :: lns	! all sizes
      integer,              intent(in)  :: comm	! which PE I am

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from initd_()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inita_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(nPEs/=size(lns))   &
        call die(myname_,'nPEs',nPEs,'size(lns)',size(lns))

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
        if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(GMap%counts,myname)
	  call mall_mci(GMap%displs,myname)
	endif

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    GMap%counts(i)=lns(i)
    l=l+lns(i)
  end do

  GMap%myPE =myID
  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

end subroutine inita_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(GMap)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(Collector),intent(inout) :: GMap

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

	if(mall_ison()) then
	  call mall_mco(GMap%counts,myname)
	  call mall_mco(GMap%displs,myname)
	endif

  deallocate(GMap%counts,GMap%displs,stat=ier)
  if(ier /= 0) call die(myname_,'deallocate()',ier)

  GMap%myPE =-1
  GMap%lsize=0
  GMap%lleft=0
  GMap%gsize=0

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - find the local size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(GMap)
      implicit none
      type(Collector),intent(in) :: GMap
      integer :: lsize_

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=GMap%lsize

end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - find the global size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    function gsize_(GMap)
      implicit none
      type(Collector),intent(in) :: GMap
      integer :: gsize_

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=GMap%gsize

end function gsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - rank (which PE) of a given global index
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rank_(GMap,i_g,rank)
      implicit none
      type(Collector),intent(in) :: GMap
      integer, intent(in)  :: i_g	! a global index
      integer, intent(out) :: rank

! !REVISION HISTORY:
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank_'
  integer :: i,ilc,ile

  rank=-1	! if nowhere fits
  do i=0,size(GMap%displs)-1
    ilc=GMap%displs(i)
    ile=ilc+GMap%counts(i)

		! If i_g in (ilc,ile].  Note that i_g := [1:..]

    if(ilc < i_g .and. i_g <= ile) then
      rank=i
      return
    endif
  end do

end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: range_ - the range of global indices on the local PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine range_(GMap,lbound,ubound,iPE)
      implicit none
      type(Collector),intent(in) :: GMap
      integer,intent(out) :: lbound,ubound
      integer,optional,intent(in) :: iPE

! !REVISION HISTORY:
!	10Apr00	- Jing Guo
!		. Added optional iPE
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::range_'

  integer :: iPE_

  iPE_=GMap%myPE
  if(present(iPE)) iPE_=iPE

  lbound=GMap%displs(iPE_)+1
  ubound=GMap%displs(iPE_)+GMap%counts(iPE_)

end subroutine range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get reference information on local or a given PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(GMap,displ,count,lbound,ubound,iPE)
      implicit none
      type(Collector),intent(in) :: GMap
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lbound
      integer,optional,intent(out) :: ubound
      integer,optional,intent(in) :: iPE

! !REVISION HISTORY:
!	10Apr00	- Jing Guo
!		. Added optional iPE
!	24Mar00	- Jing Guo
!		. initial prototype/prolog/code from m_GlobalMap
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'

  integer :: iPE_

  iPE_=GMap%myPE
  if(present(iPE)) iPE_=iPE

  if(present(lbound)) lbound=GMap%displs(iPE_)+1
  if(present(ubound)) ubound=GMap%displs(iPE_)+GMap%counts(iPE_)

  if(present(displ )) displ =GMap%displs(iPE_)
  if(present(count )) count =GMap%counts(iPE_)

end subroutine get_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_counts_ - get the reference to %counts
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_counts_(gMap)
      implicit none
      type(Collector),intent(in) :: gMap
      integer,pointer,dimension(:) :: ptr_counts_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_counts_'

  ptr_counts_ => gMap%counts

end function ptr_counts_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_displs_ - get the reference to %displs
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_displs_(gMap)
      implicit none
      type(Collector),intent(in) :: gMap
      integer,pointer,dimension(:) :: ptr_displs_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_displs_'

  ptr_displs_ => gMap%displs

end function ptr_displs_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Collector),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Collector),pointer :: new_

! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Collector),pointer :: coll
  integer :: ier

  if(present(stat)) stat=0
  allocate(coll,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'allocate()',ier)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => coll
  nullify(coll)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(coll,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Collector),pointer :: coll
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(coll,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine delete_

end module m_Collector
