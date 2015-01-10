!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_redwin 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_redwin
      implicit none
      private	! except

      public :: redwin		! data structure
      public :: qxWtb		! unit length in tau
      public :: mxWtb		! size of the table.
      public :: redwin_seplim	! combined seplim of redwin and bands

      public :: redwin_init
      public :: redwin_clean
      public :: redwin_initialized

      interface redwin_init ; module procedure init_ ; end interface
      interface redwin_clean; module procedure clean_; end interface
      interface redwin_initialized; module procedure	&
      	inited_; end interface

! !REVISION HISTORY:
! 	08Nov02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_redwin'

  	! Currently, it is one window for all corralation components.

  real,allocatable,dimension(:),save :: redwin
  real   ,save :: qxWtb
  integer,save :: mxWtb
  real   ,save :: redwin_seplim	! support of redwin in deg-angles
	!________________________________________

	integer,parameter :: NOWINDOW = 0
	integer,parameter :: GASPCOHN = 1

					! etc.

  integer,parameter :: DEFAULT_FUNCTION = NOWINDOW
  real   ,parameter :: DEFAULT_CUTDISTN = 3000.
  integer,parameter :: DEFAULT_MXTBSIZE = 15000
  real   ,parameter :: BAND_SEPLIM = 58.25	! == "band level 5"
	!________________________________________

! Usecase:
!
!	For every analysis, a redwin module needs to be defined to
!   support both the sparsity function (m_sparse) and correlation
!   matrix operator modules (m_mvcorF_nsep_bmop, m_kt_uvcorF, and
!   m_kt_corO).
!
!	call redwin_init(comm,root,rsrc)
!
!		! in m_sparse
!
!		  use m_redwin,only : redwin_initialized
!		  use m_redwin,only : redwin_seplim
!		  ASSERT(redwin_initialized())
!		  ..
!		  call BuildMask(redwin_seplim,sparse_mask,..)
!		  ..
!
!		! in m_mvcorF_nsep_bmop, m_kt_uvcorF, and m_kt_corO
!
!		  use m_redwin,only : redwin_initialized
!		  use m_redwin,only : redwin,qxWtb,mxWtb
!		  ASSERT(redwin_initialized())
!		  ..
!		  ..qxWtb
!		  ..mxWtb
!		  ..redwin(..)
!		  ..
!	  
!	call redwin_clean()

	!________________________________________

  logical,save :: initialized_=.false.
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(comm,root,rsrc)

    use m_psasrc ,only : psasrc_open,psasrc_close
    use const,only : Radius_of_EARTH	! in meters

    use m_simplePart,only : simplePart_comm

    use m_Collector,only : Collector
    use m_Collector,only : Collector_init
    use m_Collector,only : clean
    use m_CollectorComm,only : allgatherv

    use m_inpak90,only : i90_label,i90_gfloat,i90_gint
    use m_mpout  ,only : mpout_log,mpout
    use m_die    ,only : die,perr,MP_die,assert_
    use m_mall   ,only : mall_ison,mall_mci
    use m_mpif90 ,only : MP_comm_rank,MP_type

    implicit none

    integer,intent(in) :: comm
    integer,intent(in) :: root
    character(len=*),optional,intent(in) :: rsrc

    include "kind_mats.h"	! use only kind_5mat

! !REVISION HISTORY:
! 	08Nov02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  character(len=*),parameter :: rsrcCUTDISTN="redwin_cutdistance:"
  character(len=*),parameter :: rsrcFUNCTION="redwin_function:"

! This is a phased out parameter.  It is defined here to issue a
! warning only.

  character(len=*),parameter ::	&
	rsrcLevelOfBanded="level_for_banded_approximation:"

  integer :: myID
  integer :: ier		! will this be re-initialized?
  integer :: lc,le,ln,ld
  integer :: i,i_support

  real    :: dtau
  real    :: PI
  real    :: redwin_cutdistance 
  real    :: rad_support
  integer :: redwin_function

  type(collector) :: coll

  if(initialized_) call die(myname_,'module already initialized')
!_______________________________________________________________________

	call MP_comm_rank(comm,myID,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  PI=4.*atan(1.)

if(myID==root) then

	! Loading the resource handle to initialize configuration
	! variables.

  call psasrc_open(rsrc=rsrc,stat=ier)
	if(ier/=0) call die(myname_,'psasrc_open()',ier)
!________________________________________

! Each init is responsible for its own error detection and response

  redwin_cutdistance = DEFAULT_CUTDISTN
  redwin_function    = DEFAULT_FUNCTION

	! Define cut-distance, if a user has defined it.

  call i90_label(rsrcCUTDISTN,ier)
  if(ier==0) then
    redwin_cutdistance=i90_gfloat(ier)
    if(ier/=0) call die(myname_,'i90_gfloat("'//rsrcCUTDISTN//'")',ier)
    if(redwin_cutdistance<=0.)		&
	call die(myname_,'redwin_cutdistance <= 0.')
  endif

  if(redwin_cutdistance > .5 * Radius_of_Earth/1000. * PI)	&
    call die(myname_,'redwin_cutdistnace exceeding 1/4 of the sphere')

	! Redefine the form of the generating function

  call i90_label(rsrcFUNCTION,ier)
  if(ier==0) then
    redwin_function=i90_gint(ier)
    if(ier/=0) call die(myname_,'i90_gint("'//rsrcFUNCTION//'")',ier)
  endif

	! Release PSASRC

  call psasrc_close(stat=ier)
	if(ier /= 0) call die(myname_,'psasrc_close()',ier)
endif

call MPI_bcast(redwin_function,1,MP_type(redwin_function),root,comm,ier)
    if(ier/=0) call MP_die(myname_,'MPI_bcast(redwin_function)',ier)

call MPI_bcast(redwin_cutdistance,1,MP_type(redwin_cutdistance), &
	root,comm,ier)
    if(ier/=0) call MP_die(myname_,'MPI_bcast(redwin_cutdistance)',ier)

!_______________________________________________________________________

initialized_=.true.

  mxWtb=DEFAULT_MXTBSIZE

  	allocate(redwin(mxWtb),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(redwin,myname)

  call simplePart_comm(mxWtb,comm,count=ln,displ=ld)
  lc=ld+1
  le=ld+ln

  rad_support=redwin_cutdistance*2.*1000./Radius_of_Earth
  dtau=(1.-cos(rad_support))/max(1,mxWtb-1)
  qxWtb=1./dtau

  redwin_seplim=rad_support*180./PI	! rad_support in degree

! Define the default configuration

  if(ln>0) then
    select case(redwin_function)
    case(NOWINDOW)
      redwin(lc:le)=1.

    case(GASPCOHN)
      call gaspari_cohn_(redwin_cutdistance,lc,le,dtau,redwin(lc:le))

    case default
      call die(myname_,'Unknown redwin_function',redwin_function)
    end select
  endif

  	! Populate the data.

  call Collector_init(coll,ln,comm)
  call allgatherv((redwin(lc:le)),redwin,coll,comm,stat=ier)
  	if(ier/=0) call die(myname_,'allgatherv()',ier)
  call clean(coll)

  call mpout_log(myname_,'lc',lc)
  call mpout_log(myname_,'le',le)
  call mpout_log(myname_,'mxWtb',mxWtb)
  call mpout_log(myname_,'redwin(1)',redwin(1))
  call mpout_log(myname_,'redwin(lc)',redwin(lc))
  call mpout_log(myname_,'redwin(le)',redwin(le))
  call mpout_log(myname_,'redwin(mxWtb)',redwin(mxWtb))
  call mpout_log(myname_,'dtau',dtau)
  call mpout_log(myname_,'rad_support',rad_support)
  call mpout_log(myname_,'redwin_seplim',redwin_seplim)

  i_support=mxWtb	! where |redwin| < 1.e-6
  do i=mxWtb-1,1,-1
    if(abs(redwin(i)) > 1.E-6) exit
    i_support=i
  end do
  redwin_seplim=180./PI * acos(1.- (i_support-1)*dtau)
  redwin_seplim=min(redwin_seplim,BAND_SEPLIM)	! to be removed later
 
  	! List
  select case(redwin_function)
  case(NOWINDOW)
    call mpout_log(myname_,'redwin_function = NOWINDOW')

  case(GASPCOHN)
    call mpout_log(myname_,'redwin_function = GASPCOHN')

  case default
    call die(myname_,'Unknown redwin_function',redwin_function)
  end select
  call mpout_log(myname_,'redwin_cutdistance = ',redwin_cutdistance)
  call mpout_log(myname_,'adjusted redwin_seplim = ',redwin_seplim)
end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none

! !REVISION HISTORY:
! 	12Nov02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  	if(.not.initialized_) call die(myname_,'module uninitialized')

  initialized_=.false.

  	if(mall_ison()) call mall_mco(redwin,myname)
  deallocate(redwin,stat=ier)
  	if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: inited_ - status of initialization
!
! !DESCRIPTION:
!
! !INTERFACE:

    function inited_()
      implicit none
      logical :: inited_

! !REVISION HISTORY:
! 	12Nov02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inited_'

  inited_=initialized_
end function inited_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gaspari_cohn_ - initialized an equal interval array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gaspari_cohn_(cut,lc,le,dtau,redwin)
      use hcorfuns,only : i_GASPARI_COHN,corfun
      use m_realkinds,only : kind_r8
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_die ,only : die
      implicit none
      real   ,intent(in) :: cut		! in km, half of the support
      integer,intent(in) :: lc,le	! range of ctaus
      real   ,intent(in) :: dtau	! the interval
      real   ,dimension(lc:),intent(out) :: redwin  ! window function

! !REVISION HISTORY:
! 	12Nov02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gaspari_cohn_'
  integer :: i,ier
  real    :: pars(2)
  real(kind_r8),allocatable,dimension(:) :: ctaus,hcors
  real(kind_r8) :: dctau

	allocate(ctaus(lc:le),hcors(lc:le),stat=ier)
		if(ier/=0) call die(myname_,'allocate(ctaus,...)',ier)
		if(mall_ison()) then
		  call mall_mci(ctaus,myname)
		  call mall_mci(hcors,myname)
		endif

  do i=lc,le
    ctaus(i)=(i-1)*dtau
  end do

  pars(1)=cut*2.			! expected support
  pars(2)=cut/sqrt(10.d0/3.d0)		! as a correlation length

  call corfun(le-lc+1,ctaus,0.0_kind_r8,hcors,	&
    i_GASPARI_COHN,size(pars),pars,ier)
  	if(ier/=0) call die(myname_,'gaspari_cohn()',ier)

  redwin(lc:le)=hcors(lc:le)

		if(mall_ison()) then
		  call mall_mco(ctaus,myname)
		  call mall_mco(hcors,myname)
		endif
	deallocate(ctaus,hcors,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ctaus,...)',ier)

end subroutine gaspari_cohn_
end module m_redwin
