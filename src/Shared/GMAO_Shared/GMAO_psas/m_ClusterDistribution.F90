!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ClusterDistribution - Distribution by "clusters"
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ClusterDistribution
      implicit none
      private	! except

      public :: ClusterDistribution_init

      interface ClusterDistribution_init; module procedure	&
	init_; end interface

! !REVISION HISTORY:
!       02Jul01 - Tom Clune <clune@sgi.com>
!               . Minor improvement to distribution heuristic
! 	08Jun01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ClusterDistribution'

  Real, Parameter :: DEFAULT_GRANULARITY = 1. ! default granularity

#ifdef _TUNING_
  Character(Len=*), Parameter :: infile = 'granularity'
#endif

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Define (Distribution) by "clusters"
!
! !DESCRIPTION:
!
!	A "cluster" is a data unit between KR and KS.
!
! !INTERFACE:

    subroutine init_(dstr,attr,man,attr0,root,comm,symmetric,	&
    	wpower,gpart)
      use m_die  ,only : die,MP_die, assert_
      use m_mpout,only : mpout_ison,mpout_log,mpout_flush
      use m_mall ,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_type,MP_MAX,MP_SUM
      use m_mpif90,only : MP_comm_size,MP_comm_rank

      use m_Distribution,only : Distribution, Distribution_Init
      use m_Distribution,only : clean

      use m_Attributes,only : Attributes
      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : lsize
      use m_Attributes,only : permute
      use m_Attributes,only : unpermute
      use m_Attributes,only : clean
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_AttributesComm,only : trans
      use m_AttributesMAN ,only : build_fromAttr

      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_MultiAccessNavigator,only : nType
      use m_MultiAccessNavigator,only : nProf
      use m_MultiAccessNavigator,only : getType
      use m_MultiAccessNavigator,only : getProf
      use m_MultiAccessNavigator,only : ptr_krType
      use m_MultiAccessNavigator,only : ptr_krProf

      use m_Dictionary, only : Dictionary
      use m_Dictionary, only : Dictionary_clean
      use m_DictionaryTable, only : DictionaryTable_Init

      use m_GlobalPartition,only : GlobalPartition
      use m_GlobalPartition,only : ptr_partition
      use m_GlobalPartition,only : inquire

      use m_Spherical_Partition, only : NumberOfRegions
      use m_Spherical_Partition, only : LatLon_to_Region, Get

      use m_ktList,only : ktus,ktvs,ktUU,ktVV
      use m_ioutil,  only : luavail, opntext, clstext

#include "assert.H"
      
      implicit none
			! Distributions as "clusters"

      type(Distribution)        ,intent(out) :: dstr	! target
      type(Attributes)          ,intent(out) :: attr	! received
      type(MultiAccessNavigator),intent(out) :: man	! navigators

			! attr0 is in-out because a pointer association

      type(Attributes)          ,intent(inout) :: attr0
      integer                   ,intent(in) :: root
      integer                   ,intent(in) :: comm

      logical, optional         ,intent(in) :: symmetric
      real   , optional         ,intent(in) :: wpower
      type(GlobalPartition),optional,intent(in) :: gpart

! !REVISION HISTORY:
! 	08Jun01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

		! The distribution-triple (dstr,attr,man) is first
		! defined using region numbers (kr).  If the maximum
		! block size is too large, it will be re-defined using
		! "cluster IDs" as pseudo-region numbers.

  integer :: mtyp
  integer :: mxkr,mxkt
  integer :: ikr,ikt
  integer,dimension(1:2) :: mxk		! a message buffer

  integer :: npes, n_threads, nCPUs
  Integer :: unit
  Integer :: myid
  integer :: ier

  integer :: kcid
  integer :: iType,iProf
  integer :: ln,lc,le

  integer,allocatable,dimension(:,:) :: loc_pop
  integer,allocatable,dimension(:,:) :: all_pop

  integer,pointer,dimension(:) :: krp, ktp

  Type (Dictionary) :: dictTbl
  Type (Attributes) :: attrTbl

  real    :: total_work, thread_work, final_granularity
  integer :: max_seg_size
  integer :: max_cluster, n_cluster, n_bad
  logical :: symmetric_
  integer :: i
  integer :: nregs, level, max_level
  Real, Pointer :: slat(:), slon(:)
  Real :: wpower_, granularity

#ifdef _OPENMP
  Integer, External :: OMP_GET_MAX_THREADS
#endif
!________________________________________
! Basic configuration parameters

	    ! assume that work scales as the square of a segment size
	    ! then the total work is sum(all_pop**2).  we wish
	    ! to make the largest block at most 1/npes * total_work.
	    ! use _real_ quantities to prevent overflow of 32 bit
	    ! integers.

  symmetric_ = .false.
  If (Present(symmetric)) symmetric_ = symmetric

  	! If symmetric is not set, default to linear computational
	! complicity.

  wpower_=1.

	! Load balance is most critical for level 0 preconditioners.
	! These cost roughly the square of their dimension to evaluate
	! on a block by block basis.  In this case, a squired
	! computational complicity is assumed.

  if(symmetric_) wpower_=2.
  if(present(wpower)) wpower_=wpower

  n_threads = 1
#ifdef _OPENMP
  n_threads = OMP_GET_MAX_THREADS()
#endif

	call MP_comm_rank(comm,myid,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
	call MP_comm_size(comm,npes,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

  nCPUs = npes*n_threads
!________________________________________

  granularity=DEFAULT_GRANULARITY		! the default value

#ifdef _TUNING_

	! This input option is turned off to ensure that all sources of
	! input are explicit to a system integrator or operator.  This
	! configuration may be left as default if all configuration
	! data files are integrated into a single database.

	   ! Determine the target granularity (cluster size)

	If (myid == 0) Then ! read granularity file if it exists
	   unit = luavail()
	   call opntext(unit,infile,"old",ier)
	   If (ier == 0) Then
	      Read(unit,*) granularity
	      Call clstext(unit, ier)
	   Else
	      granularity = 1.0
	   End If
           If (mpout_ison()) Call mpout_log(myname_,	&
	   	'Target granularity = ',granularity)
	End If

	! send granularity everywhere
	Call mpi_bcast(granularity, 1, MP_TYPE(granularity),	&
	  root, comm, ier)
        	if(ier/=0) call die(myname_,'mpi_bcast()',ier)
#endif
!________________________________________

  	! NumberOfRegions() is the number of regions at the deepest
	! level (i.e. refinement-level).  At this level, regions are
	! likely to be more than currently in use.

  if(present(gpart)) then
    nregs=NumberOfRegions(ptr_partition(gpart))	! max region number
  else
    nregs=NumberOfRegions()
  endif

	mxk(1)=nregs			! max region number
	mxk(2)=maxval(ptr_kt(attr0))	! max type number
	mtyp=MP_type(mxk)

  call MPI_allreduce((mxk),mxk,size(mxk),mtyp,MP_MAX,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_ALLgather()',ier)

	mxkr=mxk(1)
	mxkt=mxk(2)

		! Verify if all PEs have the same gpart

	if(mxkr/=nregs)	&
	  call die(myname_,'mxkr',mxkr,'nregs',nregs)

  	! If the GlobalPartition (gpart) seems to be the same on all
	! PEs, get its parameters.

  level = 0					! baselevel
  Call Get(refinement_level = max_level)	! max refinement

  if(present(gpart))	&
    call inquire(gpart,baselevel=level,refinementlevel=max_level)

!________________________________________

  Allocate(loc_pop(mxkt,mxkr), all_pop(mxkt,mxkr), STAT = ier)
       	if(ier/=0) call die(myname_,'allocate(pop)',ier)

	if(mall_ison()) then
	  call mall_mci(loc_pop, myname_)
	  call mall_mci(all_pop, myname_)
	endif

		! Local counts for all regions and types.  The
		! additions (+ln) may be totally unnecessary.

  loc_pop(:,:)=0

  krp => ptr_kr(attr0)
  ktp => ptr_kt(attr0)

  Do i = 1, Size(krp)
     ikr = krp(i)
     ikt = ktp(i)
     Select Case(ikt)
     Case (ktvs,ktvv)
        ikt = ikt - 1
     End Select
     loc_pop(ikt,ikr) = loc_pop(ikt,ikr) + 1
  End Do
!________________________________________

		! Global counts for all regions and types.  These
		! numbers are used to estimate total work and an
		! optimal maximum "cluster" size.

  Call MPI_ALLreduce(loc_pop, all_pop, size(loc_pop),	&
  	MP_TYPE(loc_pop), MP_SUM, comm, ier)
	if(ier/=0) call MP_die(myname_,'MPI_ALLreduce()',ier)

  total_work  = Sum( Real(all_pop) ** wpower_)
  thread_work = (total_work/(nCPUs*granularity))**(1./wpower_)
  max_seg_size = max(1,Floor(thread_work))

  n_bad = Count(all_pop > max_seg_size)	! segments marked to refine
!________________________________________

	    ! Output analysis to mpout

  If (mpout_ison()) then
    Call mpout_log(myname_,'initial nregs   = ', nregs)
    Call mpout_log(myname_,'initial nsegs   = ', Count(all_pop > 0))	
    Call mpout_log(myname_,'largest segment = ', Maxval(all_pop))
    Call mpout_log(myname_,'target seg size = ', max_seg_size)
    Call mpout_log(myname_,'total work      = ', total_work)
    Call mpout_log(myname_,'thread work     = ', thread_work)
    Call mpout_log(myname_,'base level      = ', level)
    Call mpout_log(myname_,'max refinement  = ', max_level)
    Call mpout_log(myname_,'seg to refine   = ', n_bad)
  endif
!________________________________________

                ! For regions that are "overpopulated" refine
                ! until there are no such regions, or refinement is exausted
  slat => ptr_lat(attr0)
  slon => ptr_lon(attr0)

  Do
     If (level >= max_level) Exit
     
     loc_pop = 0
     Do i = 1, Size(krp)
        ikr = krp(i)
        ikt = ktp(i)
        Select Case(ikt)
        Case (ktvs,ktvv) ! Operate alongside ktus and ktuu
           ikt = ikt - 1
        End Select
        If (all_pop(ikt,ikr) > max_seg_size) Then

	  if(present(gpart)) then
            Call LatLon_to_Region(1,slat(i:i),slon(i:i),	&
	    	krp(i:i),ier,atlevel = level+1,			&
		partition=ptr_partition(gpart))
	  else
             Call LatLon_to_Region(1,slat(i:i),slon(i:i),	&
	     	krp(i:i),ier,atlevel = level+1)
	  endif
              if(ier/=0) call die(myname_,'LatLon_to_Region()',ier)

        End If
        loc_pop(ikt,krp(i)) = loc_pop(ikt,krp(i)) + 1
     End Do

	!________________________________________
	! Global counts for all regions and types.  These
	! numbers are used to estimate total work and an
	! optimal maximum "cluster" size.

     Call MPI_ALLreduce(loc_pop, all_pop, size(loc_pop),	&
     	MP_TYPE(loc_pop), MP_SUM, comm, ier)
	if(ier/=0) call MP_die(myname_,'MPI_ALLreduce()',ier)

	!________________________________________
	! assume that work scales as the square of a segment size
	! then the total work is sum(all_pop**2).  we wish
	! to make the largest block at most 1/npes * total_work.
	! use _real_ quantities to prevent overflow of 32 bit
	! integers.

     total_work = Sum( Real(all_pop) ** wpower_)
     thread_work= (total_work/(nCPUs*granularity))**(1./wpower_)
     max_seg_size = max(1,Floor(thread_work))

	!________________________________________
	! max_cluster is the estimated maximum number of "cluster"
	! segments that any region segment is divided into.

        max_cluster = 1 + (MaxVal(all_pop)-1)/max_seg_size

     If (mpout_ison()) then
       Call mpout_log(myname_,'max re-seg size = ', max_seg_size)
       Call mpout_log(myname_,'max_cluster     = ', max_cluster )
     endif
	!________________________________________
        ! increase the level and repeat if necessary
     level = level + 1
  End Do

  nullify(slat)
  nullify(slon)
!________________________________________

  If (mpout_ison()) then
    final_granularity=total_work/(nCPUs*Real(MaxVal(all_pop))**wpower_)
    Call mpout_log(myname_,'at level        = ',level)
    Call mpout_log(myname_,'n_segments      = ',Count(all_pop >0))
    Call mpout_log(myname_,'final gran.     = ',final_granularity)
  endif
!________________________________________

    ! Create pseudo region numbers for binning
    Do i = 1, Size(krp)
       
       ikt = ktp(i)
       Select Case(ikt)
       Case (ktvs,ktvv)
          ikt = ikt - 1
       End Select
       krp(i) = mxkt * (krp(i)-1) + ikt
    End Do
		!________________________________________
		! Use pseudo-region numbers to re-distribute the
		! source data, by first creating a dictionary based
		! upon all_pop(:,:), global attributes counts.  The
		! Corresponding attrTbl%kr(:) are defined by expression
		! k=1,size(all_pop) (or k=1,mxkr*mxkt).

    call DictionaryTable_init(dictTbl,attrTbl,	&
	pseudo_population = Reshape(all_pop, (/ size(all_pop) /)))
		!________________________________________
    		! Distribute based on a user defined Dictionary.  At
		! this point, attr0%kr(:) are also defined in the range
		! of k=1,mxkr*mxkt.

    call Distribution_init(dstr, attr, man, attr0,	&
             dictTbl,attrTbl, comm,wpower=wpower_)

	     	! At this point, attr%kr(:) are defined in the range
		! of k=1,mxkr*mxkt.

		!________________________________________
		! Clean workspaces

    call Dictionary_clean(dictTbl)
    call clean(attrTbl)

		if(mall_ison()) then
		  call mall_mco(loc_pop,myname_)
		  call mall_mco(all_pop,myname_)
		endif
	deallocate(loc_pop,all_pop,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

!________________________________________
		! Reset the region numbers of the source data, attr0.
		! After this, attr0%kr(:) are defined in the range of
		! kr=1,mxkr.

    krp(:) = ((krp(:)-1)/mxkt) + 1
	ALWAYS_ASSERT(All(krp <= nregs))
	ALWAYS_ASSERT(All(krp > 0))
	nullify(krp)
	nullify(ktp)

		! Reset the region numbers of the distributed data.
		! After this, attr%kr(:) are defined in the range of
		! kr=1,mxkr.

	krp => ptr_kr(attr)
    krp(:) = ((krp(:)-1)/mxkt) + 1
	ALWAYS_ASSERT(All(krp <= nregs))
	ALWAYS_ASSERT(All(krp > 0))

		! After this, man%krType(:) are defined in the range of
		! kr=1,mxkr.

	krp => ptr_krType(man)
    krp(:) = ((krp(:)-1)/mxkt) + 1
	ALWAYS_ASSERT(All(krp <= nregs))
	ALWAYS_ASSERT(All(krp > 0))

		! After this, man%krProf(:) are defined in the range of
		! kr=1,mxkr.

	krp => ptr_krProf(man)
    krp(:) = ((krp(:)-1)/mxkt) + 1
	ALWAYS_ASSERT(All(krp <= nregs))
	ALWAYS_ASSERT(All(krp > 0))

	nullify(krp)

end subroutine init_

end module m_ClusterDistribution
