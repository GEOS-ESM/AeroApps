	subroutine savevect(prefix,kr,kt,ln,ibuf)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: savevect - save a packed upper matrix to a file
!
! !SYNOPSIS:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	28Dec95 - J. Guo	- added the prolog
!_______________________________________________________________________
	implicit none

	character(len=*),intent(in) :: prefix
	integer, intent(in) :: kr,kt,ln
	integer, intent(in) :: ibuf(ln)

! Locals.

	integer  luavail,lnblnk
	external luavail,lnblnk

	integer lumat,i,j,ios

	character*80 fnmat
	
		! Set the filename

	  j=max(lnblnk(prefix),1)
	  write(fnmat,'(a,2(a,i2.2),a,i5.5,a)')
     &	    prefix(1:j),'.',kr,'.v',kt,'.',ln,'.mat'
	  do i=1,j
	    if(fnmat(i:i).eq.' ') fnmat(i:i)='_'
	  end do

		! Check the existense

	  lumat=luavail()
	  open(lumat,file=fnmat,status='new',iostat=ios,
     &	    form='unformatted',access='direct',recl=ln)

		! Save the file

	  if(ios.eq.0) then
	    write(lumat,rec=1,iostat=ios) ibuf
	    close(lumat)
	  endif

	end
!.
