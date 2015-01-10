!#define _TRACE	!
subroutine set_OEclas
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_OEclas - set information from observation error tables
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	18Jan96 - J. Guo	- (to do)
!_______________________________________________________________________
use config, only : ktmax,kxmax,lvmax_oe,ktHH
use m_stdio,only : stderr,stdout
use m_die,  only : die
use OEclass_tbl
use m_output, only : output_ObsErrCov
implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  include "kttabl.h"
  include "kxtabl.h"

!  Locals
! ========
  character(len=len('set_OEclas')) myname
  parameter( myname='set_OEclas' )

  integer iClass
  integer iStat
  integer i,l
  integer kt,kc,kx
  integer lc,ln
  character*17 rc_tmp
  character*132 line
  logical swapUC

!  Working arrays
! ================

  integer mxOEt
  parameter(mxOEt=2*ktmax)

  integer iOEt,n_OEt
  character*16 name_OEt(mxOEt)
  character*16 eName
  character*1  eType
  integer lc_OEt(mxOEt)
  integer ln_OEt(mxOEt)
  real    sig_OEt(lvmax_oe,mxOEt)
  integer ihcHHc,ivcHHc,ivcHHu

  integer listvals,lstins
  external listvals,lstins

!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! Tabulate unique and sorted OE classes

#ifdef	_TRACE
	_TRACE write(stdout,'(2a)') myname,': tabSlist()'
#endif
  nOEclas=0
  call tabSlist(kxmax,kxclas,kxmax,nOEclas,OEclas)
#ifdef	_TRACE
	_TRACE write(stdout,'(2a,i3)') myname,': nOEclas=',nOEclas
#endif

	!------------------------------------------
	! Remove ' ' and '-' from the OEclas list

  iClass=0
  do i=1,nOEclas
    if(OEclas(i).ne.' '.and.OEclas(i).ne.'-') then
      iClass=iClass+1
      if(iClass.ne.i) OEclas(iClass)=OEclas(i)
    endif
  end do
  nOEclas=iClass	! the new size
#ifdef	_TRACE
	_TRACE write(stdout,'(2a,i3)') myname,': nOEclas=',nOEclas
#endif

  if(nOEclas.le.0) call die(myname,'no valid ObsErr class')

	!------------------------------------------
	! Index all kxclas to OEclas table.

#ifdef	_TRACE
	_TRACE write(stdout,'(2a)') myname,': calling inxSlist()'
#endif
  call inxSlist(nOEclas,OEclas,kxmax,kxclas,i_kxclas)
#ifdef	_TRACE
	_TRACE write(stdout,'(2a,2i3)') myname,': i_kxclas=',	&
	_TRACE	i_kxclas(1),i_kxclas(kxmax)
#endif
!-----------------------------------------------------------------------
	! Search the resource file for the OE level table.

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling rdlevels()'
#endif
  call rdlevels(RC_OEplev,lvmax_oe,nlev_oe,plev_oe,levtype,iStat)
  if(iStat.ne.0)	&
	call die(myname,'rdlevels("'//trim(RC_OEplev)//'")',iStat)

	!------------------------------------------
	! Compare level type with "pres~sure", not case sensitive

  select case (levtype)
  case ('pres','Pres','pressure','Pressure')
	! all valid names, do nothing

  case default
    write(stderr,'(4a)') myname,': unexpected level type, "',levtype,'"'
    call die(myname)
  end select

	!------------------------------------------
	! Echo the level input
if(output_ObsErrCov) then
  l=max(listvals(nlev_oe,plev_oe,line),1)
  write(stdout,'(/5a)') RC_OEplev,' ',levtype,' ',line(1:l)
endif

!-----------------------------------------------------------------------
	! Initialize all components of oetabl.h

  voecH_u(1:nOEclas)=' '
  loc_sigOu(1:ktmax,1:nOEclas)=0
  len_sigOu(1:ktmax,1:nOEclas)=0
  sigOu(1:lvmax_oe,1:ktmax,1:nOEclas)=-1.

  hoecH_c(1:nOEclas)=' '
  voecH_c(1:nOEclas)=' '
  loc_sigOc(1:ktmax,1:nOEclas)=0
  len_sigOc(1:ktmax,1:nOEclas)=0
  sigOc(1:lvmax_oe,1:ktmax,1:nOEclas)=-1.

!   ..For all tabulated classes, process the resource file to set the
!   obs. err. information for all variable types and both 'c' and 'u'
!   error types, including sigO (sqrt(variances)) and i_hcHH/i_vcHH
!   (function types).

  do iClass=1,nOEclas
!============

	! Creating a resource name for the class

    l=len_trim(OEclas(iClass))
    rc_tmp='ObsErr*'//OEclas(iClass)(1:l)//'::'

	! Search the resource file and return the information

#ifdef	_TRACE
    	_TRACE	l=max(len_trim(rc_tmp),1)
	_TRACE	write(stdout,'(5a)') myname,': rdoetbl("',	&
	_TRACE	  rc_tmp(1:l),'")', OEclas(iClass)
#endif

    call rdoetbl(rc_tmp, mxOEt,n_OEt,name_OEt,		&
      lvmax_oe,lc_OEt,ln_OEt,sig_OEt,			&
      ivcHHu,ivcHHc,ihcHHc, iStat			)

    if(iStat.lt.0) then		! is there a mistake?
      l=max(len_trim(rc_tmp),1)
      write(stderr,'(4a)') myname,': table not found, "',rc_tmp(1:l),'"'
      call die(myname)
    endif

    if(iStat.gt.0) then		! input error
      l=max(len_trim(rc_tmp),1)
      write(stderr,'(4a,i3)') myname,				&
	': from rdoetbl() with ',rc_tmp(1:l),', iStat =',iStat
      call die(myname)
    endif

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a,i2,3i5)') myname,		&
	_TRACE	  ': from rdoetbl(), with n_OEt,ihc,ivc,ivu=',	&
	_TRACE	  n_OEt,ihcHHc,ivcHHc,ivcHHu
#endif

    if((ihcHHc.eq.0.or.ihcHHc.eq.-1).and.ivcHHu.eq.0) then
	! entries with only vCor_HH or hCor_HH=-1, such as
	! for kx=7,8,9 in the earlier resource format

	! to be changed later to direct name them
      hoecH_c(iClass)=' '
      voecH_c(iClass)=' '
      voecH_u(iClass)=' '

      if(ivcHHc.ne.0) write(voecH_u(iClass),'(i1)') ivcHHc
      swapUC=.true.

    else
	! regular entries as designed

	! to be changed later to direct name them
      hoecH_c(iClass)=' '
      voecH_c(iClass)=' '
      voecH_u(iClass)=' '

      if(ihcHHc.gt.0) write(hoecH_c(iClass),'(i1)') ihcHHc
      if(ivcHHc.gt.0) write(voecH_c(iClass),'(i1)') ivcHHc
      if(ivcHHu.gt.0) write(voecH_u(iClass),'(i1)') ivcHHu
      swapUC=.false.
    endif

    ALL_TYPES: do iOEt=1,n_OEt
  !===========
     
      eName=name_OEt(iOEt)

		! bracket the variable type

      l=index(eName,'.')-1
      if(l.le.0) l=max(len_trim(eName),1)

		! index the variable type
      kt=0
      if(l.gt.0) kt=lstins(ktmax,ktname,eName(1:l),.true.)

		! bracket the error type

      eType='u'
      if(kt.eq.ktHH.and..not.swapUC) eType='c'
      if(l.gt.0. .and. l+2.le.len_trim(eName)) eType=eName(l+2:l+2)

		! check for possible errors

      l=max(len_trim(eName),1)
      if(kt.eq.0) then
	write(stdout,'(4a)') myname,				&
		  ': WARNING - unknown variable type for "',	&
		  eType(1:l),'", entry ignored'

	cycle ALL_TYPES
      endif

      if(kt.ne.ktHH.and.eType.ne.'u') then
	write(stdout,'(4a)') myname,				&
		  ': WARNING - unknown error type for "',	&
		  eType(1:l),'", entry ignored'

	cycle ALL_TYPES
      endif

      if( lc_OEt(iOEt) .eq. 0 .or. ln_OEt(iOEt) .eq.0 ) then

	write(stdout,'(4a)') myname,				&
		  ': WARNING - no valid data value for "',	&
		  eType(1:l),'", entry ignored'

	cycle ALL_TYPES
      endif

		! mark this variable type (permitted in analysis)

      select case (eType)
      case ('u')
	loc_sigOu(kt,iClass)=lc_OEt(iOEt)
	len_sigOu(kt,iClass)=ln_OEt(iOEt)
	sigOu(:,kt,iClass)=sig_OEt(:,iOEt)

      case ('c')
	loc_sigOc(kt,iClass)=lc_OEt(iOEt)
	len_sigOc(kt,iClass)=ln_OEt(iOEt)
	sigOc(:,kt,iClass)=sig_OEt(:,iOEt)

      case default
	l=max(len_trim(eName),1)
	write(stdout,'(4a)') myname,				&
		  ': WARNING - unknown error type ignored, "',	&
		  eType(1:l),'"'
      end select

    end do ALL_TYPES

    do kt=1,ktmax

	! Make sure the two numbers are consistent

      if( len_sigOc(kt,iClass).eq.0 ) loc_sigOc(kt,iClass)=0
      if( len_sigOu(kt,iClass).eq.0 ) loc_sigOu(kt,iClass)=0

	! Verify the validity of each variable type based on if any
	! of two error types for a variable type has been specified.

      KTclas(kt,iClass)=loc_sigOc(kt,iClass).gt.0 .or.		&
			loc_sigOu(kt,iClass).gt.0
    end do

  end do

	! pass the datatype mask(logical) for the class to the mask
	! (integer) for the data source.  This is needed because both
	! F90 and F77 code do not share the same module/common

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': kxmax loop'
#endif
  do kx=1,kxmax
    kc=i_kxclas(kx)
    do kt=1,ktmax
      kxtmask(kx,kt)=0
      if(kc.gt.0) then
         if(KTclas(kt,kc)) kxtmask(kx,kt)=1
      end if
    end do
  end do

  if(.not. output_ObsErrCov) return
	! Echo the input, OE tables

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': echo loop'
#endif
  do iClass=1,nOEclas

	! Creating a resource name for the class

    l=len_trim(OEclas(iClass))
    rc_tmp='ObsErr*'//OEclas(iClass)(1:l)//'::'
    l=len_trim(rc_tmp)

    write(stdout,'(/a)') rc_tmp(1:l)

    do kt=1,ktmax

      if(KTclas(kt,iClass)) then
	lc=loc_sigOc(kt,iClass)
	ln=len_sigOc(kt,iClass)

	if(ln.gt.0) then		! a valid var.c type
          l=len_trim(ktname(kt))
	  eName=ktname(kt)(1:l)//'.c'

	  write(stdout,'(2x,a,$)') eName
	  do l=1,lc-1
	    write(stdout,'(x,a5,$)') '   - '
	  end do
	  do l=lc,lc+ln-1
	    write(stdout,'(x,f5.1,$)') sigOc(l,kt,iClass)
	  end do
	  write(stdout,*)
	endif

	lc=loc_sigOu(kt,iClass)
	ln=len_sigOu(kt,iClass)

	if(ln.gt.0) then		! a valid var.u type
          l=len_trim(ktname(kt))
	  eName=ktname(kt)(1:l)//'.u'

	  write(stdout,'(2x,a,$)') eName
	  do l=1,lc-1
	    write(stdout,'(x,a5,$)') '   - '
	  end do
	  do l=lc,lc+ln-1
	    write(stdout,'(x,f5.1,$)') sigOu(l,kt,iClass)
	  end do
	  write(stdout,*)
	endif
      endif
    end do

    if(KTclas(ktHH,iClass)) then
      if(hoecH_c(iClass).ne.' ')				&
	write(stdout,'(2x,a,2x,a)') 'hCor_HH.c',hoecH_c(iClass)
      if(voecH_c(iClass).ne.' ')				&
	write(stdout,'(2x,a,2x,a)') 'vCor_HH.c',voecH_c(iClass)
      if(voecH_u(iClass).ne.' ')				&
	write(stdout,'(2x,a,2x,a)') 'vCor_HH.u',voecH_u(iClass)
    endif
      
    write(stdout,'(a)') '::'
  end do
end
!.
