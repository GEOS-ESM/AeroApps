subroutine rdPars(rc_Pars,					&
	ParType,ParDesc, mlev,nlev,plev,mPar,nPar,Pars, istat	)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: rdPars - read a multi-level parameter table from a file
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	05Mar96 - J. Guo	- (to do)
!_______________________________________________________________________
use m_inpak90, only : lablin, rdnext, str2rn, getwrd, getstr
use m_stdio, only : stderr	! error message output unit
implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*), intent(in)	:: rc_Pars	! the resource name

  character(len=*), intent(out)	:: ParType	! model type
  character(len=*), intent(out)	:: ParDesc	! model description

  integer, intent(in)		:: mlev		! possible level count
  integer, intent(out)		:: nlev		! actual level count
  real,    intent(out)		:: plev(mlev)	! level values
  integer, intent(in)		:: mPar		! max. no. of Pars
  integer, intent(out)		:: nPar		! actual number of Pars
  real,    intent(out)		:: Pars(mPar,mlev)	! parameters

  integer, intent(out)		:: istat	! input status
!_______________________________________________________________________

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Local vars.
  integer ios,ln,ll,ipar
  character(len=len(ParType)) snum,slev
  real val
!-----------------------------------------------------------------------
!   Parameters
  character(len=*), parameter	:: myname='rdPars'
!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  istat=0			! default return if everything OK
  npar=mpar			! most possible number of paras.

!-----------------------------------------------------------------------
!  Look for the parameter table resource
  call lablin(trim(rc_Pars))	! locate the resource block
  call rdnext(ios)	! get the next line
  if(ios.eq.2) then
    ln=max(len_trim(rc_Pars),1)
    write(stderr,'(4a)') myname,': table not found, "',rc_Pars(1:ln),'"'
    istat=-1
    return	! bad news to the parent
  endif

!-----------------------------------------------------------------------
!   Process the parameter table line by line
  ParType=' '
  nlev=0
  do while(ios.eq.0.and.nlev.lt.mlev)

	! token #1 is either a level (pressure) or the class

    call getwrd(ios,slev)
    if(ios /= 0) then
      ln=max(len_trim(slev),1)
      write(stderr,'(4a)') myname,				&
	': unexpected error from getwrd(), "',slev(1:ln),'"'
      istat=1
      return	! bad news to the parent
    endif

	! rest tokens are dependent on the type of the first one

    ll=max(len_trim(slev),1)
    ios=verify(slev(1:ll),'0123456789.')	! a f90 function
    if(ios /= 0) then
		! it can not be a valid level value, therefore,
		! possiblly a correlation class name (ParType)

      if(ParType /= ' ') then	! ParType is already defined
	write(stderr,'(4a)') myname,				&
	  ': unexpected string found, "',			&
	  slev(1:ll),'"'
	istat=2
	return
      endif

		! the rest is treated as a string

      ParType=slev		! it is a class name
      call getstr(ios,ParDesc)	! get the description
      if(ios /= 0) ParDesc='?'	! ignore error

    else
		! it might be a number, therefore, a pressure level
		! value

      snum=slev

      val=str2rn(snum,ios)
      if(ios /= 0) then		! a format problem
	write(stderr,'(6a)') myname,': level ',slev(1:ll),	&
	  ', expected a number but found "',snum(1:ll),'"'
	istat=3
	return
      endif

      nlev=nlev+1		! Count as a new level
      plev(nlev)=val

      ipar=0
      call getwrd(ios,snum)
      do while(ios == 0 .and. ipar < npar)

	val=str2rn(snum,ios)
	if(ios /= 0) then	! there is a format problem
	  ln=max(len_trim(snum),1)
	  write(stderr,'(6a)') myname,': level ',slev(1:ll),	&
	    ', expected a number but found "',snum(1:ln),'"'
	  istat=3
	  return
	endif

	ipar=ipar+1
        Pars(ipar,nlev)=val

	call getwrd(ios,snum)
      end do	! while(ios==0 .and. ipar < npar)

	! Verify the line of input

      if(ios == 0) then
		! There is more token than expected (npar)
	ln=max(len_trim(snum),1)
	write(stderr,'(6a)') myname,': WARNING - level ',	&
	  slev(1:ll),', more tokens than expected, "',snum(1:ln),'"'

      else
		! it might be alright, except ...
	if(nlev > 1 .and. ipar < npar) then
		! that means the number of tokens is not the same as
		! the first line.
	  ln=max(len_trim(snum),1)
	  write(stderr,'(4a)') myname,': level ',slev(1:ll),	&
	    ': expected more, but last token is "',snum(1:ln),'"'
	  istat=4
	  return
	endif
      endif

      npar=ipar		! the actual number of tokens
    endif

    call rdnext(ios)
  end do	! ios == 0 .and. nlev < mlev
!
		! If the table is not finished, tell the parent.
  if(ios.eq.0) then
    write(stderr,'(2a,i2,a)') myname,				&
	': buffer is full (mxlev=',nlev,'), but more in the table'
    istat=5
    return
  endif

		! If the table is empty, tell the parent.
  if(nlev.eq.0) then
    ln=max(len_trim(rc_Pars),1)
    write(stderr,'(4a)') myname,				&
	': empty table, "',rc_Pars(1:ln),'"'
    istat=6
    return
  endif
end
!.
