!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalPartition - A Global Spherical Partition Data Holder
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_GlobalPartition
      use m_Spherical_Partition,only : Spherical_Partition
      implicit none
      private	! except

		! Interfaces

      public :: GlobalPartition		! class name
      public :: GlobalPartition_init	! initialize
      public :: clean	! clean
      public :: setkr
      public :: inquire
      public :: ptr_Partition

      type GlobalPartition

		! "thePartition" is the main component.  Rests are
		! parameters controlling the partition.

	type(Spherical_Partition),pointer :: thePartition

		! BaseLevel is used to define the basic blocking that
		! may be used to define the preconditioners and matrix
		! related quantities.

	integer :: BaseLevel

		! RefinementLevel is used to define the _capacity_ of
		! the Spherical_Partition component of this
		! GlobalPartition object.

	integer :: RefinementLevel

	! Alternative axes of a partition are not implemented
	! logical :: alteraxes
	! real    :: zaxis_lon,zaxis_lat
	! real    :: xaxis_lon,xaxis_lat
      end type GlobalPartition

      interface GlobalPartition_init ; module procedure	&
	init_ ; end interface
      interface clean; module procedure	&
	clean_; end interface
      interface inquire; module procedure	&
	inquire_; end interface
      interface ptr_Partition ; module procedure	&
	ptr_; end interface

      interface setkr ; module procedure	&
	setbyks_,	&
	setbyll_
      end interface
!________________________

    public :: new	! create an object as a pointer
    public :: delete	! delete an object as a pointer

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

! !REVISION HISTORY:
! 	29Jan02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- changed from m_Regioner
!       31Oct01 - Tom Clune <clune@sgi.com>
!               . Modified to use m_Spherical_Partition
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalPartition'

!_______________________________________________________________________
!
! Additional configuration: A GlobalPartition can be configured
! to based on a cartision coordinate system different from the
! default North-South pole oriented spherical partition, by
! changing the position of the z-axis only or both the z-axis
! and the x-axis.
!
! 1) The default configuration assumes:
!
!	z-axis = (-90, 90)
!	x-axis = (  0,  0)
!	y-axis = ( 90,  0)
!
!   The values are in (degree-longutide,degree-latitude), for
!   the corresponding location of a given axis projected on the
!   sphere.
!
! 2) The user can set
!
!	z-axis = (lambda    ,   phi)
!
!   The derived axes are:
!
!	x-axis = (lambda+ 90, 0    )
!	y-axis = (lambda+180,90-phi)
!
! 3) The user can also set both
!
!	z-axis = (lambda_z    ,   phi_z)
!	x-axis = (lambda_x    ,   phi_x)
!
!   First, the orthoganality between x- and z-axis will be
!   verified through
!
!	x-axis .dot. z-axix ~= 0.
!
!   where the notation "-axis" are used as the corresponding
!   vector defined by the location on the sphere with respect to
!   the origin (the center of the sphere).  Note that the
!   meaning of "~= 0." can be implemented rather loosely,
!   because the actual x-axis will be redefined.
!
!   If the dot-product passes the test, new axes will be derived
!   from the two axes:
!
!	y-axis = unit(z-axis .cross. x-axis)
!	x-axis = unit(y-axis .cross. z-axis)
!
!   Note that x-axis is redefined to ensure the orthoganality.
!   For a given non-zero vector v, its unit vector is defined by
!
!	unit(v) = v/|v|
!
! *) This is now just an idea.  No implementation has been down.

#include "assert.H"
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Initialize the definition of the "regions"
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(obj,comm,root,rsrc)
      use m_Spherical_Partition, only : Initialize,new

      use m_inpak90,only : i90_loadf,i90_release
      use m_inpak90,only : i90_label,i90_gint,i90_gtoken

      use m_die   ,only : die,MP_die
      use m_chars ,only : uppercase
      use m_mpif90,only : MP_comm_rank,MP_type

      use m_mpout, only : mpout_log

      implicit none
      type(GlobalPartition),intent(out) :: obj
      integer,intent(in) :: comm
      integer,intent(in) :: root
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  character(len=*),parameter ::	&
	rsrcBaseLevel      ="GlobalPartition_BaseLevel:"
  character(len=*),parameter ::	&
	rsrcRefinementLevel="GlobalPartition_RefinementLevel:"

  integer,dimension(2) :: ibuf
  integer :: ier
  integer :: myID
!________________________________

	call MP_comm_rank(comm,myID,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
	! Set the default configuration before trying to read the user
	! configuration from a resource file.

    obj%RefinementLevel=1
    obj%BaseLevel=0

    if(present(rsrc)) then
	! Redefine the configuration from a resource file.

      call i90_loadf(rsrc,ier)
	if(ier/=0) call die(myname_,'i90_loadf("'//trim(rsrc)//'")',ier)

	call mpout_log(myname_, 'using "'//		&
		trim(rsrc)//'" for the runtime resource input')
	!________________________________

      call i90_label(rsrcRefinementLevel,ier)
	if(ier==0) then
			! Parsing error is a fatal error, iff its label
			! exists.

	  obj%RefinementLevel=i90_gint(ier)
	  if(ier/=0) call die(myname_,	&
		'i90_gint("'//rsrcRefinementLevel//'")',ier)
    
	  call mpout_log(myname_,	&
		'Using "'//rsrcRefinementLevel//'" = ', &
		obj%RefinementLevel)
	else
          call mpout_log(myname_,	&
		'Using default "'//rsrcRefinementLevel//'" = ', &
		obj%RefinementLevel)
	endif
	!________________________________

      call i90_label(rsrcBaseLevel,ier)
	if(ier==0) then
			! Parsing error is a fatal error, iff its label
			! exists.

	  obj%BaseLevel=i90_gint(ier)
	  if(ier/=0) call die(myname_,	&
		'i90_gint("'//rsrcBaseLevel//'"',ier)


	  call mpout_log(myname_,	&
		'Using "'//rsrcBaseLevel//'" = ',obj%BaseLevel)
	else
    	  call mpout_log(myname_,	&
	    'Using default "'//rsrcBaseLevel//'" = ',obj%BaseLevel)
	endif
	!________________________________

      call i90_release(ier)
	if(ier/=0) call die(myname_,	&
		'i90_release("'//trim(rsrc)//'")',ier)
    endif

  endif
!________________________________

	! Broadcast the input data

	ibuf(1)=obj%RefinementLevel
	ibuf(2)=obj%BaseLevel

  call MPI_Bcast(ibuf,size(ibuf),MP_type(ibuf),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_Bcast()',ier)

	obj%RefinementLevel=ibuf(1)
	obj%BaseLevel      =ibuf(2)
!________________________________

	! Because of the use of "atlevel=" argument in calling routine
	! LatLon_to_Region(), for atlevel values different from the
	! initial level of the Spherical_Partition (n_levels=
	! Region_Refinement_Level), argument compress should always be
	! set to .false., to ensure m_sparse() to be useful for
	! different atlevel values for the same Spherical_Partition.

	obj%thePartition => new(obj%thePartition)
  Call Initialize(n_levels =obj%RefinementLevel,	&
		  partition=obj%thePartition,compress = .false.)

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - In case the region definition becoming complicate
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(obj)
      use m_Spherical_Partition,only : clean,delete
      use m_die, only : die
      implicit none
      type(GlobalPartition),intent(inout) :: obj

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  obj%RefinementLevel=1
  obj%BaseLevel=0

  call clean(partition=obj%thePartition)
  call delete(obj%thePartition)
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ - module data reference
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_(obj)
      use m_Spherical_Partition,only : Spherical_Partition
      use m_die,only : die
      implicit none
      type(GlobalPartition),target,intent(in) :: obj
      type(Spherical_Partition),pointer :: ptr_

! !REVISION HISTORY:
! 	29Jan02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_'

  ptr_ => obj%thePartition
end function ptr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setbyks_ - set region IDs for soundings
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setbyks_(kr,obj,kx,ks,rlat,rlon,comm)
      use m_die ,only : die,assert_
      implicit none

      integer,dimension(:) ,intent(out) :: kr

      type(GlobalPartition),intent(in)  :: obj
      integer,dimension(:) ,intent(in)  :: kx
      integer,dimension(:) ,intent(in)  :: ks
      real   ,dimension(:) ,intent(in)  :: rlat
      real   ,dimension(:) ,intent(in)  :: rlon
      integer,intent(in) :: comm

! !REVISION HISTORY:
!	29Jan02	- Jing Guo
!		. Changed to merge the changes between MPI PSAS,
!		  MPI-"optimized" PSAS, "reduced" PSAS, and OOp PSAS.
! 	23Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setbyks_'
  integer :: n

  n=size(kr)
  if(n<=0) return

	ASSERT(n==size(kx))
	ASSERT(n==size(ks))
	ASSERT(n==size(rlat))
	ASSERT(n==size(rlon))
!________________________________________
!
!	This algorithm needs major fixes.  This is because it assumes
!   that data from the same sounding are on the same local processor
!   and ordered together.  This is not true in general.  A correct
!   algorithm should 1) define globally a list of sounding locations,
!   2)  assign the partition indices (kr) for the soundings based on
!   their locations, then 3) assign every partition index to all data
!   in the same sounding.

  call setbyll_(kr,obj,rlat,rlon,comm)

end subroutine setbyks_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setbyll_ - set region IDs for lat-long datasets
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setbyll_(kr,obj,rlat,rlon,comm)
      use m_Spherical_Partition, only : LatLon_to_region
      use m_die,only : die,assert_
      implicit none
      integer,dimension(:) ,intent(out) :: kr

      type(GlobalPartition),intent(in)  :: obj
      real   ,dimension(:) ,intent(in)  :: rlat
      real   ,dimension(:) ,intent(in)  :: rlon
      integer,intent(in) :: comm

! !REVISION HISTORY:
!	29Jan02	- Jing Guo
!		. Changed to merge the changes between MPI PSAS,
!		  MPI-"optimized" PSAS, "reduced" PSAS, and OOp PSAS.
! 	23Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setbyll_'
  integer :: n
  integer :: ier

  n=size(kr)
  if(n<=0) return

	ASSERT(n==size(rlat))
	ASSERT(n==size(rlon))

  call LatLon_to_Region(n,rlat(1:n),rlon(1:n),kr(1:n),ier,	&
	partition=obj%thePartition,atlevel=obj%BaseLevel)

	if(ier/=0) call die(myname_,'LatLon_to_Region()',ier)

end subroutine setbyll_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: inquire_ - inquiring configuration information
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine inquire_(obj,baselevel,refinementlevel)
      use m_die,only : die
      implicit none
      type(GlobalPartition),intent(in) :: obj
      integer,optional,intent(out) :: baselevel
      integer,optional,intent(out) :: refinementlevel

! !REVISION HISTORY:
! 	30Jan02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inquire_'

  if(present(baselevel)) baselevel=obj%Baselevel
  if(present(refinementlevel)) refinementlevel=obj%Refinementlevel
end subroutine inquire_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(GlobalPartition),pointer :: mold
      integer,optional,intent(out) :: stat
      type(GlobalPartition),pointer :: new_

! !REVISION HISTORY:
! 	07Feb02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(GlobalPartition),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => obj
  nullify(obj)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(obj,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(GlobalPartition),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	07Feb02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine delete_
end module m_GlobalPartition
