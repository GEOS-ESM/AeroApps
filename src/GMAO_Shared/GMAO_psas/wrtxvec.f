	subroutine wrtxvec(luvec,nvec,nymd,nhms,
     &	  rlats,rlons,rlevs,kx,kt,Xvec,sigO,sigF,rtims,rqflg,ierr)
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: wrtxvec
! 
! !DESCRIPTION: 
!	Creating a del2 format output file for the given vector
!
! !CALLING SEQUENCE: 
!	call opnieee(luvec,...)	! see also opnieee()
!	call wrtxvec(luvec,nvec,nymd,nhms,
!    &	  rlats,rlons,rlevs,rtims,kx,kt,xvec,sigF,sigO)
!
! !INPUT PARAMETERS: 
!	integer		luvec		! logical unit for the output
!	integer		nvec		! size of xvec
!	integer		nymd		! yymmdd date of the data
!	integer		nhms		! hhmmss time of the data
!	real		rlats(nvec)	! latitudes
!	real		rlons(nvec)	! longitudes
!	real		rlevs(nvec)	! levels(pressure)
!	integer		kx(nvec)	! data sources
!	integer		kt(nvec)	! data types
!	real		Xvec(nvec)	! data values
!	real		sigO(nvec)	! observation errors
!	real		sigF(nvec)	! forecast errors
!	real 		rtims(nvec)	! times(in minutes)
!	real		rqflg(nvec)	! quality flags
!
! !OUTPUT PARAMETERS:
!	integer		ierr		! status, !=0 in case of error
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
! 	22Aug95 - J. Guo	- added the prolog
!
!-----------------------------------------------------------------------
	implicit none
c-----------------------------------------------------------------------
c   ..Input arguments
	integer		luvec		! logical unit for the output
	integer		nvec		! size of xvec
	integer		nymd		! yymmdd date of the data
	integer		nhms		! hhmmss time of the data
	real		rlats(nvec)	! latitudes
	real		rlons(nvec)	! longitudes
	real		rlevs(nvec)	! levels(pressure)
	integer		kx(nvec)	! data sources
	integer		kt(nvec)	! data types
	real		Xvec(nvec)	! data values
	real		sigO(nvec)	! observation errors
	real		sigF(nvec)	! forecast errors
	real 		rtims(nvec)	! times(in minutes)
	real		rqflg(nvec)	! quality flags

c-----------------------------------------------------------------------
c   ..Output arguments
	integer		ierr		! status, !=0 in case of error

c-----------------------------------------------------------------------
c   ..Local vars.

	integer n

	integer lenrecs		! size (words) of a record
	parameter(lenrecs=10)
c		=sizeof(rlats(1))+sizeof(rlons(1))+...+sizeof(rqflg(1))

c------------------------------------------------------
c   ..Creating a header record
	write(luvec,iostat=ierr) nvec,lenrecs,nymd,nhms

c   ..Write records.
	n=0
	do while(ierr.eq.0.and.n.lt.nvec)
	  n=n+1
	  write(luvec,iostat=ierr) rlats(n),rlons(n),rlevs(n),
     &	    kx(n),kt(n),Xvec(n),sigO(n),sigF(n),rtims(n),rqflg(n)
	end do

	end
c.
