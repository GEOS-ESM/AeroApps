subroutine setbox(nbmax, nboxes, boxes)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: setbox - read the resource file data selection box settings
!
! !INTERFACE:
!	<@interface
!
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	09May96 - J. Guo	- changed the label data structure to
!				table data structure.
!				- some documentation changes.
!				- now coded in Fortran 90.
!				- interface changes
!
!	06Jan95	- Jing G.	- Extended the boxes from three spatial
!				  dimensions to 6 hyper dimensions,
!				  including kx and kt.
!
!  04oct94 A. da S. - Changed box definition from hardwired DATA
!	statement to input from resource file with inpak77. Notice that
!	calling convention has changed.
!
!  04feb93 ?	- code modifications based on: analy.f 
!_______________________________________________________________________

use m_inpak90, only : lablin, rdnext, str2rn, getwrd
use config, only : kxmax, ktmax
use m_stdio,only : stdout, stderr
use m_die,  only : die
implicit none

integer, intent(in)	:: nbmax	! maximum number of boxes
integer, intent(out)	:: nboxes	! number of output boxes
real,    intent(out)	:: boxes(2,6,nbmax)	! box information

!=======================================================================
!   Parameters
  character(len=*), parameter	:: myname='setbox'
  character(len=*), parameter	:: RC_boxes='DataBoxes::'

!-----------------------------------------------------------------------
!   Local vars.
  integer ios,ln,ibox,iset
  integer i1,i2
  real    x1,x2
  character(len=16) snum
  real    val
  logical valid

!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nboxes=0
!-----------------------------------------------------------------------
!  Look for the table resource

  call lablin(trim(RC_boxes))	! locate the resource block
  call rdnext(ios)	! get the next line
  if(ios.eq.2) then
    ln=max(len_trim(RC_boxes),1)
    write(stderr,'(4a)') myname,': table not found, "',	&
      &	RC_boxes(1:ln),'"'
    call die(myname)
  endif

!-----------------------------------------------------------------------
!   Process the parameter table line by line
  ibox=0
  do while(ios.eq.0.and.ibox.lt.nbmax)

    ibox=ibox+1

    do iset=1,6
	!--------------------------------
		! the first number of the set

      call getwrd(ios,snum)
      if(ios == 0) val=str2rn(snum,ios)
      if(ios /= 0) then
        ln=max(len_trim(snum),1)
        write(stderr,'(4a)') myname,				&
	  ': expected a number but found "',snum(1:ln),'"'
        call die(myname)
      endif

      boxes(1,iset,ibox)=val

	!--------------------------------
		! the second number

      call getwrd(ios,snum)
      if(ios == 0) val=str2rn(snum,ios)
      if(ios /= 0) then
        ln=max(len_trim(snum),1)
        write(stderr,'(4a)') myname,				&
	  ': expected a number but found "',snum(1:ln),'"'
        call die(myname)
      endif

      boxes(2,iset,ibox)=val
    end do

    call rdnext(ios)
  end do	! ios == 0 .and. ibox < nbmax
!-----------------------------------------------------------------------
	! If the table is not finished, report the error

  if(ios.eq.0) then
    write(stderr,'(2a,i2,a)') myname,				&
	': buffer is full (nbmax=',nbmax,'), but more in the table'
    call die(myname)
  endif

  nboxes=ibox	! number of actual boxes
!=======================================================================
	! verify the inputs

  valid=.true.
  write(stdout,'(a)') RC_boxes
  do ibox=1,nboxes

    write(stdout,'(2i4,x,2i3,x,2f4.0,x,2f5.0,x,2f9.3,x,2f5.0)') &
	& nint(boxes(1,1,ibox)),nint(boxes(2,1,ibox)),		&
	& nint(boxes(1,2,ibox)),nint(boxes(2,2,ibox)),		&
	&      boxes(1,3,ibox) ,     boxes(2,3,ibox) ,		&
	&      boxes(1,4,ibox) ,     boxes(2,4,ibox) ,		&
	&      boxes(1,5,ibox) ,     boxes(2,5,ibox) ,		&
	&      boxes(1,6,ibox) ,     boxes(2,6,ibox)

	! parameter set 1 is kx

    i1=nint(boxes(1,1,ibox))
    i2=nint(boxes(2,1,ibox))
    if(i1 < 1 .or. i1 > kxmax .or. i2 < 1 .or. i2 > kxmax) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 1'
      valid=.false.
    endif

	! parameter set 2 is kt

    i1=nint(boxes(1,2,ibox))
    i2=nint(boxes(2,2,ibox))
    if(i1 < 1 .or. i1 > kxmax .or. i2 < 1 .or. i2 > kxmax) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 2'
      valid=.false.
    endif

	! parameter set 3 is latitude

    x1=boxes(1,3,ibox)
    x2=boxes(2,3,ibox)
    if(x1 < -90. .or. x1 > +90. .or. x2 < -90. .or. x2 > +90.) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 3'
      valid=.false.
    endif

	! parameter set 4 is longitude

    x1=boxes(1,4,ibox)
    x2=boxes(2,4,ibox)
    if(x1 < -180. .or. x1 > +180. .or. x2 < -180. .or. x2 > +180.) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 4'
      valid=.false.
    endif

	! parameter set 5 is pressure

    x1=boxes(1,5,ibox)
    x2=boxes(2,5,ibox)
    if(x1 <= 0. .or. x2 <= 0.) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 5'
      valid=.false.
    endif

	! parameter set 6 is longitude

    x1=boxes(1,6,ibox)
    x2=boxes(2,6,ibox)
    if(x1 < -180. .or. x1 > +180. .or. x2 < -180. .or. x2 > +180.) then
      write(stderr,'(2a,i3,a)') myname,': invalid box, ibox =',	&
	& ibox,' and iset = 6'
      valid=.false.
    endif

  end do
  write(stdout,'(a)') '::'

  if(.not.valid) then
    write(stderr,'(4a)') myname,': invalid table, "',RC_boxes,'"'
    call die(myname)
  endif
end
!.
