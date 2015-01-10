!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sigmaPhi -
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sigmaPhi
      use m_Navigator,only : Navigator
      implicit none
      private	! except

      public :: sigmaPhi	! the class data structure
      public :: sigmaPhi_init
      public :: clean
      public :: sigmaPhi_matvec
      public :: sigmaPhi_adjvec

      public :: col_size
      public :: row_size

      public :: ptr_sigma

      public :: col_Navigator
      public :: col_krNav
      public :: col_ktNav

      public :: col_rlon
      public :: col_rlat
      public :: col_rlev

		! For now, it is not clear how to define
		! diag(sigmaPhi(DCD')sigmaPhi').  It really should be
		! done right sometime.

      type sigmaPhi
	private

	type(Navigator),pointer :: nav
	integer,pointer,dimension(:) :: krNav
	integer,pointer,dimension(:) :: ktNav

	integer :: nrow
	integer :: ncol
	real,pointer,dimension(:) :: rlon
	real,pointer,dimension(:) :: rlat
	real,pointer,dimension(:) :: rlev
	real,pointer,dimension(:) :: sigma
      end type sigmaPhi

      interface sigmaPhi_init ; module procedure	&
	init_ ; end interface
      interface clean; module procedure	&
	clean_; end interface

      interface sigmaPhi_matvec; module procedure	&
	matvec_; end interface
      interface sigmaPhi_adjvec; module procedure	&
	adjvec_; end interface

      interface col_size; module procedure	&
	col_size_; end interface
      interface row_size; module procedure	&
	row_size_; end interface

      interface ptr_sigma; module procedure	&
	ptr_sigma_; end interface

      interface col_Navigator; module procedure	&
	col_Navigator_; end interface
      interface col_krNav; module procedure	&
	col_krNav_; end interface
      interface col_ktNav; module procedure	&
	col_ktNav_; end interface

      interface col_rlon; module procedure	&
	col_rlon_; end interface
      interface col_rlat; module procedure	&
	col_rlat_; end interface
      interface col_rlev; module procedure	&
	col_rlev_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sigmaPhi'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a sigmaPhi operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(sigPhi,nav,krNav,ktNav,rlon,rlat,rlev)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci
      use m_xOp_sigFi,only : intp_sigFi

      implicit none

      type(sigmaPhi),intent(out) :: sigPhi

      type(Navigator),target,intent(in) :: nav
      integer,target,dimension(:),intent(in) :: krNav
      integer,target,dimension(:),intent(in) :: ktNav

      real,target,dimension(:),intent(in) :: rlon
      real,target,dimension(:),intent(in) :: rlat
      real,target,dimension(:),intent(in) :: rlev

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nsize
  integer :: ier

	! For this operator, the column attributes are the same as the
	! row attributes, except the types of the variables.  The types
	! of the column variables are scaled row variables.

  sigPhi%rlon => rlon
  sigPhi%rlat => rlat
  sigPhi%rlev => rlev

  sigPhi%nav   => nav
  sigPhi%krNav => krNav
  sigPhi%ktNav => ktNav

  nsize=size(rlon)

	allocate(sigPhi%sigma(nsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(sigPhi%sigma,myname)

  call intp_sigFi(nav,krNav,ktNav,rlat,rlon,rlev,sigPhi%sigma)

  sigPhi%nrow=nsize
  sigPhi%ncol=nsize

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a sigmaPhi operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(sigPhi)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(sigmaPhi),intent(inout) :: sigPhi

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

		if(mall_ison()) call mall_mco(sigPhi%sigma,myname)
	deallocate(sigPhi%sigma,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  nullify(sigPhi%rlon)
  nullify(sigPhi%rlat)
  nullify(sigPhi%rlev)

  nullify(sigPhi%nav)
  nullify(sigPhi%krNav)
  nullify(sigPhi%ktNav)

  sigPhi%nrow=0
  sigPhi%ncol=0

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: matvec_ -  
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine matvec_(sigPhi, nvecs,Cvec,Rvec		)

      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die

      implicit none

      type(sigmaPhi),intent(in) :: sigPhi

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Cvec
      real   ,dimension(:,:),intent(out) :: Rvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code to replace mv_diag()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::matvec_'

  integer :: ivec
  integer :: lc,le
  integer :: inav
  integer :: ncol,nrow

!-----------------------------------------------------------------------
  ncol=size(Cvec,1)
  nrow=size(Rvec,1)

	! loop over all Rvec partitions

  do inav=1,lsize(sigPhi%nav)
    call get(sigPhi%nav,inav,lc=lc,le=le)

    if(lc<1   ) call die(myname_,'invalid lc',lc)
    if(le>ncol) call die(myname_,'invalid le',le)
    if(le>nrow) call die(myname_,'invalid le',le)

    do ivec=1,nvecs
      Rvec(lc:le,ivec)=sigPhi%sigma(lc:le)*Cvec(lc:le,ivec)
    end do
  end do

end subroutine matvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: adjvec_ -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine adjvec_(sigPhi, nvecs,Rvec,Cvec		)

      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die

      implicit none

      type(sigmaPhi),intent(in) :: sigPhi

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Rvec
      real   ,dimension(:,:),intent(out) :: Cvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::adjvec_'

  integer :: ivec
  integer :: lc,le
  integer :: inav
  integer :: ncol,nrow

!-----------------------------------------------------------------------
  ncol=size(Cvec,1)
  nrow=size(Rvec,1)

	! loop over all Rvec partitions

  do inav=1,lsize(sigPhi%nav)
    call get(sigPhi%nav,inav,lc=lc,le=le)

    if(lc<1   ) call die(myname_,'invalid lc',lc)
    if(le>ncol) call die(myname_,'invalid le',le)
    if(le>nrow) call die(myname_,'invalid le',le)

    do ivec=1,nvecs
      Cvec(lc:le,ivec)=sigPhi%sigma(lc:le)*Rvec(lc:le,ivec)
    end do
  end do

end subroutine adjvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_size_ - the size of the columns
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_size_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      integer :: col_size_

! !REVISION HISTORY:
! 	15Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_size_'

  col_size_ = sigPhi%ncol

end function col_size_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_size_ - the size of the rows
!
! !DESCRIPTION:
!
! !INTERFACE:

    function row_size_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      integer :: row_size_

! !REVISION HISTORY:
! 	15Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::row_size_'

  row_size_ = sigPhi%nrow

end function row_size_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_sigma_ - referencing %sigma
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_sigma_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      real,pointer,dimension(:) :: ptr_sigma_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_sigma_'

  ptr_sigma_ => sigPhi%sigma

end function ptr_sigma_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_Navigator_ - referencing the column navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_Navigator_(sigPhi)
      use m_Navigator,only : Navigator
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      type(Navigator),pointer :: col_Navigator_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_Navigator_'

  col_Navigator_ => sigPhi%nav

end function col_Navigator_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_krNav_ - referencing the column krNav
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_krNav_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      integer,pointer,dimension(:) :: col_krNav_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_krNav_'

  col_krNav_ => sigPhi%krNav

end function col_krNav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_ktNav_ - referencing the column ktNav
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_ktNav_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      integer,pointer,dimension(:) :: col_ktNav_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_ktNav_'

  col_ktNav_ => sigPhi%ktNav

end function col_ktNav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_rlon_ - referencing the column rlon
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_rlon_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      real,pointer,dimension(:) :: col_rlon_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_rlon_'

  col_rlon_ => sigPhi%rlon

end function col_rlon_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_rlat_ - referencing the column rlat
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_rlat_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      real,pointer,dimension(:) :: col_rlat_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_rlat_'

  col_rlat_ => sigPhi%rlat

end function col_rlat_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: col_rlev_ - referencing the column rlev
!
! !DESCRIPTION:
!
! !INTERFACE:

    function col_rlev_(sigPhi)
      implicit none
      type(sigmaPhi),intent(in) :: sigPhi
      real,pointer,dimension(:) :: col_rlev_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::col_rlev_'

  col_rlev_ => sigPhi%rlev

end function col_rlev_

end module m_sigmaPhi
