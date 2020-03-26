!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !MODULE: m_GrADSfiles - Read-only files in the GrADS format
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_GrADSfiles

      implicit none
      private

      public	:: GrADSfiles	! the class data stucture
      public	:: GrADS_open	! open a GrADSfiles instance
      public	:: GrADS_input	! read a GrADS 2d/3d field
      public	:: GrADS_close	! close a GrADSfiles instance
      public	:: GrADS_getdims! get the dimensions of a variable
      public	:: GrADS_zdef	! get the levels of a variable

		! The current version only has limited supports to a
		! GrADS file, with restrictions on the grid, time,
		! and the file format, etc.

        type GrADSfiles

	  integer		:: nvar	! number of variables
	  integer		:: idim		! longitude dimension
	  integer		:: jdim		! latitude dimension
	  integer		:: kdim		! level dimension
	  integer		:: ldim		! time dimension

	  real,dimension(:),pointer :: zdef	! size(zdef)=kdim

	  integer		:: irec		! current location
	  integer		:: iacc		! access control
	  integer		:: ilen		! record length
	  real			:: undef	! missing value flag

	  integer		:: nblock
	  character(len=12),dimension(:),pointer :: vname
	  integer,          dimension(:),pointer :: n_rec
	  integer,          dimension(:),pointer :: i_rec

	  integer	     :: lu	! logical unit if already opened
	  character(len=256) :: file	! the filename for input

	  real*4,dimension(:,:),pointer :: dbuf
        end type GrADSfiles

	! Interface definitions

      interface GrADS_open;  module procedure open_;  end interface
      interface GrADS_close; module procedure close_; end interface
      interface GrADS_input; module procedure	&
	input3d_,	&
	input2d_
      end interface

      interface GrADS_getdims; module procedure getdims_; end interface
      interface GrADS_zdef;    module procedure zdef_;    end interface

! !EXAMPLES: (to do)
! !BUGS: (to do)
!
!   Output interfaces should be added soon.  See aio_grads.f for more
!   information.
!
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	16Jul96 - J. Guo	- (to do)
!_______________________________________________________________________

!=======================================================================
!
! !REVISION HISTORY
! 	grads.h - last change: Wed Jul 20 21:00:31 EDT 1994 (ams)
!				- Original source from A. da Silva
!	01Dec94 - Jing G. -	added zdef_gr for small values
!	16Jul96 - Jing G. -	combined to form a module
! file: grads.f - last change: Wed Jul 20 21:00:31 EDT 1994 (ams)
!
!  Routines to read in GrADS like files.
!......................................................................

	! access methods of a Fortran unformatted file
	!---------------------------------------------

  integer, parameter :: iacc_DIRECT = 1
  integer, parameter :: iacc_SEQUENTIAL = 2

  integer, parameter :: stat_DEFINED = 1
  integer, parameter :: stat_UNDEF   = 0

  character(len=*),parameter :: myname='m_GrADSfiles'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !IROUTINE: open_ - open an input GrADS "control" file for input
!
! !DESCRIPTION: (to do)
! !INTERFACE:

    subroutine open_( gs, ctl_file, stat)

      use m_inpak90, only : i90_loadf
      use m_inpak90, only : i90_release
      use m_inpak90, only : i90_label
      use m_inpak90, only : i90_gtoken
      use m_inpak90, only : i90_gint
      use m_inpak90, only : i90_gfloat
      use m_inpak90, only : i90_gline
      use m_chars,   only : uppercase
      use m_stdio,   only : stderr
      use config,    only : FILL
      use m_die,     only : perr
      use m_die,     only : die
      use m_mall,    only : mall_mci,mall_ison
      use m_ioutil,  only : luavail

      implicit none

      type(GrADSfiles), intent(out) :: gs
      character(len=*), intent(in)  :: ctl_file	! filename
      integer, optional,intent(out) :: stat	! status

! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
!	21Jan00	- Jing Guo
!		. Added "direct" access to open_()
! 	16Jul96 - J. Guo	- modified as a Fortran 90 module.
!	01Dec94 - Jing G.	- added zdef_gr for small values.
!				- Original source from A. da Silva
!_______________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::open_'

! Local variables:

	character(len=64) :: str
	integer  i,k,l,lu

	integer :: idim,jdim,kdim,ldim,nvar
	integer :: ierr
	logical :: formatdefined

	if(present(stat)) stat=0

!----------------------------------------
! Use m_inpak90 read the table file
!----------------------------------------
	call i90_loadf(ctl_file,ierr)
	if(ierr /= 0 ) then
	  write(stderr,'(4a,i3)') myname_,			&
	    ': i90_loadf("',trim(ctl_file),'") error, ierr =', ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif

!----------------------------------------
!  Mandatory GrADS settings:
!
!	dset xdef ydef zdef tdef vars
!----------------------------------------
! DSET
	call i90_label('DSET',ierr)
	if(ierr /= 0) call i90_label('dset',ierr)
	if(ierr == 0) call i90_gtoken(gs%file, ierr)
	if(ierr /= 0) then
	  write(stderr,'(4a,i3)') myname_,			&
	    ': DSET error with "',trim(ctl_file),'", ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	if(gs%file(1:1) == '^') then
	  gs%file=gs%file(2:)
	  i=index(ctl_file,'/',back=.true.)
	  if(i > 0) gs%file=ctl_file(1:i)//gs%file
	endif
!----------------------------------------
! XDEF
	gs%idim = -1
	call i90_label('XDEF',ierr)
	if(ierr /= 0) call i90_label('xdef',ierr)
	if(ierr == 0) gs%idim = i90_gint(ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,': XDEF error, ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	idim=gs%idim
!----------------------------------------
! YDEF
	gs%jdim = -1
	call i90_label('YDEF',ierr)
	if(ierr /= 0) call i90_label('ydef',ierr)
	if(ierr == 0) gs%jdim = i90_gint(ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,': YDEF error, ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	jdim=gs%jdim
!----------------------------------------
! ZDEF
	gs%kdim = -1
	call i90_label('ZDEF',ierr)
	if(ierr /= 0) call i90_label('zdef',ierr)
	if(ierr == 0) gs%kdim = i90_gint(ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,': ZDEF error, ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	kdim=gs%kdim

	if(kdim>0) then
	  call i90_gtoken(str,ierr)
	  if(ierr /=0) then
	    write(stderr,'(2a,i3)') myname_,	&
		': ZDEF type error, ierr = ',ierr
	    if(.not.present(stat)) call die(myname_)
	    stat=ierr
	    return
	  endif

	  str=uppercase(str)
	  select case(str)
	  case ('LEVELS')

	    allocate(gs%zdef(1:kdim),stat=ierr)
	    if(ierr /= 0) then
	      write(stderr,'(2a,i4)') myname_,	&
		': allocate(gs%zdef) error, stat = ',ierr
	      if(.not.present(stat)) call die(myname_)
	      stat=ierr
	      return
	    endif

		if(mall_ison()) call mall_mci(gs%zdef,myname)

	    ierr=0
	    do k=1,kdim
	      if(ierr.eq.0) gs%zdef(k)=i90_gfloat(ierr)
	    end do
	    if(ierr /= 0) then
	      write(stderr,'(2a,i3)') myname_,	&
		': ZDEF level error, ierr = ',ierr
	      if(.not.present(stat)) call die(myname_)
	      stat=ierr
	      return
	    endif

	  case default
	    write(stderr,'(4a)') myname_,	&
		': unknown ZDEF type, "',trim(str),'"'
	    if(.not.present(stat)) call die(myname_)
	    stat=1
	    return
	  end select

	endif
!----------------------------------------
! TDEF
	gs%ldim=-1
	call i90_label('TDEF',ierr)
	if(ierr /= 0) call i90_label('tdef',ierr)
	if(ierr == 0) gs%ldim=i90_gint(ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,': TDEF error, ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	ldim=gs%ldim
!----------------------------------------
! VARS -- ENDVARS
	gs%nvar=-1
	call i90_label('VARS',ierr)
	if(ierr /= 0) call i90_label('vars',ierr)
	if(ierr == 0) gs%nvar=i90_gint(ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,': VARS error, ierr = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif
	nvar=gs%nvar

	allocate( gs%vname(nvar),gs%n_rec(nvar),	&
		  gs%i_rec(nvar), stat=ierr)

	if(ierr /= 0) then
	  write(stderr,'(2a,i4)') myname_,	&
	    ': allocate(VARS) error, stat = ',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif

		if(mall_ison()) then
		  call mall_mci(gs%vname,myname)
		  call mall_mci(gs%n_rec,myname)
		  call mall_mci(gs%i_rec,myname)
		endif

!     Get variable names and labels
!     -----------------------------
	ierr=0
	do k=1,nvar
	  if(ierr == 0) call i90_gline(ierr)
	  if(ierr == 0) call i90_gtoken(str,ierr)
	  if(ierr == 0) l=i90_gint(ierr)
	  if(ierr == 0) then
	    gs%vname(k)=str
	    gs%n_rec(k)=l
	  endif
	end do
	if(ierr /= 0) then
	  write(stderr,'(2a,i3)') myname_,	&
	    ': VARS entry error, line ',k+1
	  if(.not.present(stat)) call die(myname_)
	  stat=2
	  return
	endif

!	if(ierr == 0 .and. (.upper.trim(str) == 'ENDVARS') exit
!
!	if(k /= nvar) then
!	  write(stdout,'(2a)') myname_,	&
!	    ': mismatched VARS <nvar> and ENDVARS?'
!	  gs%nvar=min(nvar,k)
!	  nvar=gs%nvar
!	  return
!	endif

	gs%i_rec(1)=1
	do k=2,nvar
	  gs%i_rec(k)=gs%i_rec(k-1)+max(gs%n_rec(k-1),1)
	end do
	gs%nblock = gs%i_rec(nvar)+max(gs%n_rec(nvar),1) -1

!----------------------------------------
! Optional settings
!
!	options
!----------------------------------------
! OPTIONS 
	formatdefined=.false.

	gs%iacc = iacc_DIRECT
	call i90_label('OPTIONS',ierr)
	if(ierr /= 0) call i90_label('options',ierr)
	if(ierr == 0) then

	  call i90_gtoken(str,ierr)
	  if(ierr == 0) then

	    str=uppercase(str)
	    select case(str)
	    case('SEQUENTIAL')
	      formatdefined=.true.
	      gs%iacc = iacc_SEQUENTIAL

	    case('DIRECT')
	      formatdefined=.true.
	      gs%iacc = iacc_DIRECT

	    case default
	      write(stderr,'(4a)') myname_,	&
	        ': unsupported option, "',trim(str),'"'
	      if(.not.present(stat)) call die(myname_)
	      stat=3
	      return
	    end select

	  endif
	endif

!----------------------------------------
! Optional settings
!
!	format
!----------------------------------------
! FORMAT
	if(.not. formatdefined) then
	  gs%iacc = iacc_DIRECT
	  call i90_label('FORMAT',ierr)
	  if(ierr /= 0) call i90_label('format',ierr)
	  if(ierr == 0 ) then
	    call i90_gtoken(str,ierr)
	    if(ierr == 0) then

	      str=uppercase(str)
	      if(str /= 'DIRECT' .and. str /= 'SEQUENTIAL') then
	        write(stderr,'(4a)') myname_,	&
	          ': invalid FORMAT option, "',trim(str),'"'
	        if(.not.present(stat)) call die(myname_)
	        stat=3
	        return
	      endif

	      if(str == 'SEQUENTIAL') gs%iacc = iacc_SEQUENTIAL
	    endif
	  endif
	endif

!----------------------------------------
! Optional settings
!----------------------------------------
! UNDEF
	gs%undef=FILL
	call i90_label('UNDEF',ierr)
	if(ierr /= 0) call i90_label('undef',ierr)
	if(ierr == 0) then
	  gs%undef=i90_gfloat(ierr)
	  if(ierr /= 0) then
	    write(stderr,'(2a,i3)') myname_,': UNDEF error, ierr =',ierr
	    if(.not.present(stat)) stat=0
	    stat=3
	    return
	  endif
	endif

!----------------------------------------
	call i90_release(ierr)
	  if(ierr /= 0) then
	    call perr(myname_,'i90_release()',ierr)
	    if(.not.present(stat)) call die(myname_)
	    stat=ierr
	    return
	  endif
!----------------------------------------

	lu=gs%lu
	if(lu >= 0) close(lu)
	gs%lu=-1

	!--------------------------------------------------------
		! allocate the input buffer

	allocate(gs%dbuf(gs%idim,gs%jdim), stat=ierr)
	if(ierr /= 0) then
	  write(stderr,'(2a,i5)') myname_,	&
	    ': allocate(gs%dbuf) error, stat =',ierr
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
	endif

		if(mall_ison()) call mall_mci(gs%dbuf,myname)

	gs%lu = luavail()
	gs%irec=1
	call opendset_(gs%lu,gs%file,gs%iacc,gs%dbuf,gs%ilen,ierr)
		if(ierr/=0) then
		  call perr(myname_,'opendset_()',ierr)
		  if(.not.present(stat)) call die(myname_)
		  stat=ierr
		  return
		endif

end subroutine open_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: opendset_ - open a DSET file
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine opendset_(lu,name,iacc,dbuf,ilen,ierr)
      use m_ioutil,only : opnieee
      use m_die,   only : perr
      implicit none

      integer,         intent(in)  :: lu
      character(len=*),intent(in)  :: name
      integer,         intent(in)  :: iacc
      real*4,dimension(:,:),intent(in) :: dbuf

      integer,         intent(out) :: ilen
      integer,         intent(out) :: ierr

! !REVISION HISTORY:
! 	21Jan00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::opendset_'
  character(len=16) :: clen


  ilen=0
  inquire(iolength=ilen) dbuf

  select case(iacc)
  case(iacc_SEQUENTIAL)

    ilen=0		! reset %ilen to avoid confusion
    call opnieee(lu,trim(name),'old',ierr)
	  if(ierr.ne.0) then
	    call perr(myname_,'opnieee("'//	&
		trim(name)//'")',ierr		)
	    return
	  endif

  case(iacc_DIRECT)

    call opnieee(lu,trim(name),'old',ierr,recl=ilen)
	  if(ierr.ne.0) then
	    clen='****************'
	    write(clen,'(i16)',iostat=ierr) ilen
	    clen=adjustl(clen)
	    call perr(myname_,'opnieee("'//	&
		trim(name)//'",recl='//		&
		trim(clen)//')',ierr		)
	    return
	  endif
	
  case default

	  call perr(myname_,'unknown iacc',iacc)
	  return
  end select

	!--------------------------------------------------------
end subroutine opendset_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !IROUTINE: close_ - close a GrADSfiles variable
!
! !DESCRIPTION: (to do)
!
! !INTERFACE:

    subroutine close_(gs,stat)
      use m_die, only : perr,die
      use m_mall,only : mall_mco,mall_ison
      use m_ioutil,only : clsieee
      implicit none

      type(GrADSfiles), intent(inout) :: gs
      integer,optional, intent(out)   :: stat

! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	18Mar97 - Jing Guo <guo@eramus> - initial prototyping and coding
!_______________________________________________________________________
  character(len=*), parameter :: myname_=myname//'::close_'

  integer :: lu
  integer :: ier

  if(present(stat)) stat=0

  lu=gs%lu
  call clsieee(lu,ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) then
	  call mall_mco(gs%zdef ,myname)
	  call mall_mco(gs%vname,myname)
	  call mall_mco(gs%n_rec,myname)
	  call mall_mco(gs%i_rec,myname)
	  call mall_mco(gs%dbuf ,myname)
	endif

  deallocate(gs%zdef,gs%vname,gs%n_rec,gs%i_rec,gs%dbuf,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  gs%nvar = -1

  gs%idim = -1
  gs%jdim = -1
  gs%kdim = -1
  gs%ldim = -1

  gs%nblock = -1

  gs%irec = 0
  gs%iacc = iacc_DIRECT
  gs%ilen = -1

  gs%file = ' '
  gs%lu=-1

end subroutine close_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: input2d_ - input a 2-d field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine input2d_(gs,vnam,llev,klev,vfld,stat)

      use m_stdio, only : stderr
      use m_die,   only : die

      implicit none

      type(GrADSfiles), intent(inout) :: gs	! the input
      character(len=*), intent(in)    :: vnam	! what variable?
      integer,          intent(in)    :: llev	! what time?
      integer,          intent(in)    :: klev	! which level?
      real,dimension(:,:),intent(out) :: vfld	! a 2-d gridded field
      integer, optional,intent(out)   :: stat

! !REVISION HISTORY:
! 	18Mar97 - Jing Guo <guo@eramus> - initial prototyping and coding
! 	08Dec98 - Jing Guo <guo@thunder> - modified from read_ with
!			dbuf for portable input
!EOP ___________________________________________________________________

  character(len=*), parameter :: myname_=myname//'::input2d_'
  integer :: i,ierr,nrec,nskp,ivar,lu
  logical :: no_buffer

  if(present(stat)) stat=0

	! Sanity checks

		! Check the file status

  if( gs%file == ' ' .or. gs%lu   <  0 .or.	&
      gs%nvar <= 0   .or. gs%idim <= 0 .or.	&
      gs%jdim <= 0   .or. gs%ldim <= 0		) then

    write(stderr,'(2a)') myname_,': uninitialized type(GrADSfiles)?'
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

		! Check the buffer dimensions

  if(size(vfld,1) /= gs%idim .or. size(vfld,2) /= gs%jdim) then

    write(stderr,'(2a,$)') myname_,': invalid arguments'
    write(stderr,'(a,2i6,a,$)') ', shape(vfld) = (',shape(vfld),')'
    write(stderr,'(a,2i6,a,$)') ', gs%[ij]dim = (',gs%idim,gs%jdim,')'
    write(stderr,*)
    if(.not.present(stat)) call die(myname_)
    stat=2
    return
  endif

		! Check/index the requested variable

  ivar=lindex_(gs%nvar,gs%vname,vnam)
  if(ivar <= 0) then
    write(stderr,'(4a)') myname_,': unknown variable "',trim(vnam),'"'
    if(.not.present(stat)) call die(myname_)
    stat=3
    return
  endif

		! Check the requested time

  if(llev < 1 .or. llev > gs%ldim) then
    write(stderr,'(2a,$)') myname_,': invalid time request'
    write(stderr,'(2(a,i3))') ', llev =',llev,', gs%ldim =',gs%ldim
    if(.not.present(stat)) call die(myname_)
    stat=4
    return
  endif

		! Check the requested level

  if(klev < 1 .or. klev > gs%n_rec(ivar) .or. klev > gs%kdim) then
    write(stderr,'(2a,$)') myname_,': invalid level request'
    write(stderr,'(a,i3,3a,i3,a,i3)') ', klev =',klev,		&
	', gs%n_rec("',trim(vnam),'") =', gs%n_rec(ivar),	&
	', gs%kdim =',gs%kdim
    if(.not.present(stat)) call die(myname_)
    stat=5
    return
  endif

	!--------------------------------------------------------
	! Compute the record number

  nrec = (llev-1) * gs%nblock + gs%i_rec(ivar) + klev-1

	!--------------------------------------------------------
	! Read the nrec-th record.  The current position may be
	! taken into account if the file is sequentially accessed.
	! See read_() for details.

  call read_(gs%lu,gs%iacc,nrec,gs%irec,vfld,gs%dbuf,ierr)

  if(ierr /= 0) then
    write(stderr,'(2a,i5)') myname_,	&
	': read_() error, ierr =',ierr
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  endif

  gs%irec=nrec+1	! current record position

	!--------------------------------------------------------
end subroutine input2d_
!=======================================================================

  function lindex_(nlst,lsts,entr)
    use m_chars, only : uppercase
    implicit none
    integer,                      intent(in) :: nlst
    character(len=*),dimension(:),intent(in) :: lsts
    character(len=*),             intent(in) :: entr

    integer :: lindex_	! the result

	!--------------------------------------------------------
    integer :: i

	!--------------------------------------------------------
    lindex_=0
    do i=1,nlst
      if(uppercase(entr) == uppercase(lsts(i))) then
	lindex_=i
	return
      endif
    end do
  end function lindex_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read_ - read the n-th record
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read_(lu,iacc,nrec,irec,vfld,dbuf,ierr)
      implicit none
      integer,intent(in)  :: lu		! the input unit
      integer,intent(in)  :: iacc	! access
      integer,intent(in)  :: nrec	! which to read
      integer,intent(in)  :: irec	! where it is now
      real,  dimension(:,:),intent(out) :: vfld
      real*4,dimension(:,:),intent(out) :: dbuf
      integer,intent(out) :: ierr

! !REVISION HISTORY:
! 	22Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read_'
  logical :: no_buffer
  integer :: nskp
  integer :: i
  integer :: nx,ny

  no_buffer =	kind(dbuf)   == kind(vfld)	.and.	&
		size(dbuf,1) == size(vfld,1)	.and.	&
		size(dbuf,2) == size(vfld,2)

  ierr=-1
  select case(iacc)
  case(iacc_SEQUENTIAL)

			! Sequential skip
    ierr=0
    if(nrec < irec) then
      rewind(lu)		! can we trust backspace()?
      nskp=nrec-1
    else
      nskp=nrec-irec
    endif
    do i=1,nskp
      if(ierr == 0) read(lu,iostat=ierr)
    end do

			! Sequential read()
    if(ierr == 0) then
      if(no_buffer) then
        read(lu,iostat=ierr) vfld
      else
        read(lu,iostat=ierr) dbuf
	if(ierr == 0) then
	  nx=min(size(vfld,1),size(dbuf,1))
	  ny=min(size(vfld,2),size(dbuf,2))
	  vfld(1:nx,1:ny)=dbuf(1:nx,1:ny)
	endif
      endif
    endif

  case(iacc_DIRECT)
			! Direct read()
    if(ierr == 0) then
      if(no_buffer) then
	read(lu,rec=nrec,iostat=ierr) vfld
      else
	read(lu,rec=nrec,iostat=ierr) dbuf
	if(ierr == 0) then
	  nx=min(size(vfld,1),size(dbuf,1))
	  ny=min(size(vfld,2),size(dbuf,2))
	  vfld(1:nx,1:ny)=dbuf(1:nx,1:ny)
	endif
      endif
    endif
  end select
end subroutine read_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getdims_ get the dimensions of a given variable
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getdims_(gs,vnam,nlon,nlat,nlev,udef,stat)
      use m_stdio, only : stderr
      use m_die,   only : die
      implicit none
      character(len=*), intent(in) :: vnam
      type(GrADSfiles), intent(in) :: gs
      integer,intent(out) :: nlon
      integer,intent(out) :: nlat
      integer,intent(out) :: nlev
      real,   optional,intent(out) :: udef
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getdims_'
  integer :: ivar

  if(present(stat)) stat=0

  nlon=gs%idim
  nlat=gs%jdim

  ivar=lindex_(gs%nvar,gs%vname,vnam)
  if(ivar==0) then
    write(stderr,'(4a)') myname_,	&
	': unknown variable name, "',vnam,'"'
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  nlev=gs%n_rec(ivar)
  if(present(udef)) udef=gs%undef

	! A post-condition?

  if(nlev < 0 .or. nlev > gs%kdim) then
    write(stderr,'(4a,i4)') myname_,	&
	': improper number of records, gs%n_rec("',vnam,'") =', nlev
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

end subroutine getdims_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: zdef_ - return the leading levels of ZDEF
!
! !DESCRIPTION:
!
! !INTERFACE:

    function zdef_(gs,nlev)
      use m_stdio, only : stderr
      use m_die,   only : die
      implicit none
      type(GrADSfiles),intent(in) :: gs
      integer,intent(in) :: nlev
      real,dimension(nlev) :: zdef_

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::zdef_'

  if(nlev > size(gs%zdef) ) then
    write(stderr,'(2a,i4)') myname_,	&
	': improper number of "ZDEF" levels, nlev =', nlev
    call die(myname_)
  endif

  zdef_(1:nlev)=gs%zdef(1:nlev)

end function zdef_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: input3d_ - input a 3-d field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine input3d_(gs,vnam,llev,vfld,stat)
      use m_stdio, only : stderr
      use m_die,   only : die
      implicit none
      type(GrADSfiles),intent(inout) :: gs
      character(len=*),intent(in) :: vnam
      integer,         intent(in) :: llev
      real,            intent(out):: vfld(:,:,:)
      integer,optional,intent(out):: stat

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::input3d_'

  integer :: nlon,nlat,nlev
  integer :: k,ier

  if(present(stat)) stat=0

  call getdims_(gs,vnam,nlon,nlat,nlev)

  do k=1,nlev
    call input2d_(gs,vnam,llev,k, vfld(:,:,k), stat=ier)

    if(ier /= 0) then
      write(stderr,'(2a,i4)') myname_,	&
	': input2d_() error, stat =',ier
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
  end do

end subroutine input3d_

end module m_GrADSfiles
!.
