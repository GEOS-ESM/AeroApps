	subroutine savemat1(prefix,kr,kt,ln,cor)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: savemat1 - save a packed upper matrix to a file
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
	real,    intent(in) :: cor(ln*(ln+1)/2)

! Locals.
	integer  luavail,lnblnk
	external luavail,lnblnk

	integer lumat,i,j,ios

	character*80 fnmat
	
	integer lx
	lx(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

		! Set the filename

	  j=max(lnblnk(prefix),1)
	  write(fnmat,'(a,4(a,i2.2),2(a,i5.5),a)') prefix(1:j),
     &	    '.',kr,'-',kr,'.v',kt,'-',kt,'.',ln,'x',ln,'.mat'
	  do i=1,j
	    if(fnmat(i:i).eq.' ') fnmat(i:i)='_'
	  end do

		! Check existense

	  lumat=luavail()
	  open(lumat,file=fnmat,status='new',iostat=ios,
     &	    form='unformatted',access='direct',recl=ln)

		! Save the file

	  if(ios.eq.0) then

	    do j=1,ln
	      write(lumat,rec=j,iostat=ios) (cor(lx(i,j)),i=1,ln)
	    end do

	    close(lumat)
	  endif

	end
!.
