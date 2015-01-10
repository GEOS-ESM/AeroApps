!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: slintab - set indices and linear weights respect to a table
!
! !INTERFACE:
  subroutine slintab(nearest,ntab,vtab, ndat,vdat,indx,wght)
    implicit none
    logical, intent(in)		:: nearest	! nearest entries?
    integer, intent(in)		:: ntab		! size of vtab
      real,    intent(in)	:: vtab(ntab)	! a lookup table

    integer, intent(in)		:: ndat		! size of vdat
      real,    intent(in)	:: vdat(ndat)	! a list of values

    integer, intent(out)	:: indx(ndat)	! indices to vtab
    real,    intent(out)	:: wght(ndat)	! weights of vtab
! end subroutine slintab
!
! !DESCRIPTION:
!	*slintab()* sets indices and linear weights for a given data
!	list *vdat(ndat)* respect to a given lookup table *vtab(ntab)*.
!	If *nearest* is .false., indices *indx* point to *vtab(ntab)*
!	entries of nearest (in a linear sense) values, with *wght(ndat)*
!	defined within the range [-.5,+.5).  If *nearest* is .true.,
!	indices *indx* are lower entries of table intervals, with *wght*
!	defined withing [0,1).
!
!	*nearest*=.true. may be used to create an index table for
!	direct table reference, while *nearest*=.false. may be used to
!	create an index table with weights for linear interpolating
!	functions.
!
!	Notice that *vtab(ntab)* must be in a sorted order.
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
!	18May96 - J. Guo	- Fixed a problem with out-of-ranged
!				  data
! 	05Dec95 - J. Guo	- initial coding
!_______________________________________________________________________
!_______________________________________________________________________

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! local vars.
	integer i,j
	real v
!_______________________________________________________________________

	! extreme condition.  Treated as normal

	if(ntab.le.0) return

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ..Index the array to the lower index entry of an intervals in vtab.
	call vindex(ntab,vtab,ndat,vdat,indx)

!  ..Compute the linear weights, and adjust the indices if needed.
	do i=1,ndat
	  j=indx(i)

	  v=0.
	  if(j.lt.ntab) v=(vdat(i)-vtab(j))/(vtab(j+1)-vtab(j))
	  v=min(max(0.,v),1.)

	  if(nearest) then
		! indx(i) is the index to the nearest entry in vtab,
		! with wght(i) is within [-.5,+.5).  It is likely to
		! be used for direct table reference

	    indx(i)=j + nint(v)
	    wght(i)=v - nint(v)

	  else
		! indx(i) is the lower index of an interval in vtab,
		! with wght(i) is in [0,1).
	    wght(i)=v
	  endif
	end do
!_______________________________________________________________________
	end
!.
