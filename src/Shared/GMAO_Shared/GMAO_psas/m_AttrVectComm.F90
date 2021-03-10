!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVectComm - Communications of AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AttrVectComm
      use m_AttrVect, only : AttrVect
      implicit none
      private	! except

      public :: gather		! gather all local vectors to the root
      public :: scatter		! scatter from the root to all PEs
      public :: bcast		! bcast from root to all PEs

    interface gather ; module procedure gather_ ; end interface
    interface scatter; module procedure scatter_; end interface
    interface bcast  ; module procedure bcast_  ; end interface

! !REVISION HISTORY:
! 	10Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.  This 
!                 change affects gather_(), scatter_(), and bcast_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVectComm'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gather_ - gather a vector according to a given _map_
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gather_(iV,oV,Mp,root,comm,stat)
      use m_stdio
      use m_die
      use m_GlobalMap
      use m_mpif90
      implicit none
      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	15Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gather_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _scatterd_ vector

  niV=lsize(Mp)
  noV=lsize_(iV)
  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, lsize(Mp) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV=gsize(Mp)		! the gathered local size, as for the output
  if(myID /= root) nov=0
  call initv_(oV,iV,noV)

  niV=lsize(Mp)		! the scattered local size, as for the input

  nIA=nIAttr_(oV)	! number of INTEGER attributes
  nRA=nRAttr_(oV)	! number of REAL attributes

  call MPI_gatherv(iV%iAttr,niV*nIA,MP_INTEGER,			&
	oV%iAttr,Mp%counts*nIA,Mp%displs*nIA,MP_INTEGER,	&
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_gatherv(iV%rAttr,niV*nRA,MP_type(iV%rAttr(1,1)),	     &
	oV%rAttr,Mp%counts*nRA,Mp%displs*nRA,MP_type(oV%rAttr(1,1)), &
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine gather_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatter_ - scatter a vecter according to a given _map_
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatter_(iV,oV,Mp,root,comm,stat)
      use m_stdio
      use m_die
      use m_GlobalMap
      use m_mpif90
      implicit none
      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatter_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _gathered_ vector

  niV=gsize(Mp)		! the _gathered_ local size
  if(myID /= root) niV=0

  noV=lsize(iV)
  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, rsize(Mp) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV=lsize(Mp)		! the _scatterd_ local size
  call initv_(oV,iV,noV)

  nIA=nIAttr_(iV)	! number of INTEGER attributes
  nRA=nRAttr_(iV)	! number of REAL attributes

  call MPI_scatterv(iV%iAttr(1,1),Mp%counts*nIA,	&
	Mp%displs*nIA,MP_INTEGER,			&
	oV%iAttr(1,1),noV*nRA,MP_INTEGER,root,comm,ier )
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_scatterv(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_scatterv(iV%rAttr(1,1),Mp%counts*nRA,	&
	Mp%displs*nRA,MP_type(iV%rAttr(1,1)),		&
	oV%rAttr(1,1),noV*nRA,MP_type(oV%rAttr(1,1)),   &
        root, comm,ier )
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_scatterv(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast from the root to all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine bcast_(aV,root,comm,stat)
      use m_die, only : die, perr
      use m_mpif90
      use m_String, only : String,bcast,char
      use m_List, only : get
      implicit none
      type(AttrVect) :: aV	! (IN) on the root, (OUT) elsewhere
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'
  type(String) :: iLStr,rLStr
  integer :: nIA,nRA,lsize
  integer :: myID
  integer :: ier

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Convert the two Lists to two Strings

  if(myID == root)	&
    call get(iLStr,aV%iList)

  call bcast(iLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(iLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root)	&
    call get(rLStr,aV%rList)

  call bcast(rLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(rLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root) lsize=lsize_(aV)

	! Set lsize for all PEs

  call MPI_bcast(lsize,1,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(lsize)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID /= root)	&
	call init_(aV,iList=char(iLStr),rList=char(rLStr),lsize=lsize)

  nIA=nIAttr_(aV)
  nRA=nRAttr_(aV)

  call MPI_bcast(aV%iAttr,nIA*lsize,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_bcast(aV%rAttr,nRA*lsize,MP_type(aV%rAttr(1,1)), &
                 root, comm, ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine bcast_
end module m_AttrVectComm
!.
