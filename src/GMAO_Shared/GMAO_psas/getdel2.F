!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: getdel2 - Get the quality controled deltas in "del2" format
!
! !SYNOPSIS:

      subroutine GETDEL2 ( del2file,nymd,nhms,nobsmax,nobs,
     &	rlats,rlons,rlevs,tstamp,kx,kt,del,qcflag,exclude,
     &	verbose,idelprb,idelpre,idelpri,ierr )

! !INPUT PARAMETERS:
!		logical verbose		! output control
!		integer idelprb		! for selected output begin from
!		integer idelpre		! for selected output end at
!		integer idelpri		! for selected output step by
!		integer nobsmax		! size of the output arrays
!		integer nymd		! for this YYMMDD only,
!		integer nhms		! for this HHMMSS only
!
! !OUTPUT PARAMETERS:
!		real rlats(nobsmax)	! latitudes of del2
!		real rlons(nobsmax)	! longitudes of del2
!		real rlevs(nobsmax)	! levels(pressures) of del2
!		integer kx(nobsmax)	! data sources of del2
!		integer kt(nobsmax)	! variable types of del2
!		real   del(nobsmax)	! values of del2
!		real  sigO(nobsmax)	! sqrt of obs. error variances
!		real  sigF(nobsmax)	! sqrt of fcst. error variances
!		real tstamp(nobsmax)	! time of del2
!		integer qcflag(nobsmax)	! QC flags of del2
!		integer nymd		! actual YYMMDD, see above
!		integer nhms		! actual HHMMSS, see above
!		integer nobs		! actual size of del2 arrays
!		character*(*) del2file	! filename of del2
!		integer ierr		! input status
!
! !DESCRIPTION:
!	"getdel2()" reads a "deltas" data file in ".del2." format.  It
!	gets the input filename from a "psas.rc" ASCII file that has
!	been opened before getdel2() is called.
!
!	On the return, getdel2() has ierr=0 if the input is considered
!	to be successful.
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
! 	19Sep95 - J. Guo	- added the prolog
!   02Feb95 - Jing G.	- Changed open() to opnieee() to remove the
!			  explicit system I/O dependency on number
!			  representations.
!   22Dec94 - Jing G.	- Fixed a typo of ludel ("rewind ludelf").  Also
!			  made the unit number determined by luavail().
!			  One more reason why "implicit none" should be
!			  used.
!			- luavail() is used to seek any logical unit
!			  which may become available at the time.
!   04oct94 - A. da S.	- major changes: change of input parameters,
!                    file name now comes from resource file and
!                    eliminated calls to opndel/clsdel.
!
!   05nov93 - ???	- code modified to filter out observations which
!   			  did not pass Quality Control
!   04feb93 - ???	- getdel2.f
!_______________________________________________________________________
	use m_inpak90, onlY : lablin,strget
	use m_ioutil,  only : luavail,opnieee
	use m_stdio,   only : stderr,stdout

	use m_stdio
	implicit none

*     Argument declarations.
*     ---------------------
      logical      verbose
      integer      idelprb
      integer      idelpre
      integer      idelpri
      integer      nobsmax

      real         rlats(nobsmax)
      real         rlons(nobsmax)
      real         rlevs(nobsmax)
      integer      kx(nobsmax)
      integer      kt(nobsmax)
      real         del(nobsmax)
      real         sigO(nobsmax)
      real         sigF(nobsmax)
      real         tstamp(nobsmax)
      integer      qcflag(nobsmax)	! no longer used in the solver
      logical      exclude		! if qcflag(i)==0 only

      character*(*) del2file

      integer      nymd
      integer      nhms
      integer      nobs
      integer      ierr

      integer lf,ludel
      integer mqcf
      integer irec,iobs,iymd,iy2kmd,ihms,lenrecs
      integer ikx,ikt
#ifdef _UNICOS
      real xlat,xlon,xlev,xdel,xsigO,xsigF,xtim,xqcf
#else
      real*4 xlat,xlon,xlev,xdel,xsigO,xsigF,xtim,xqcf
#endif

*.................................................................

      character*80 filen

	character*7 myname
	parameter(myname='getdel2')

      filen=del2file
      lf=len_trim(filen)

*     Open the oiVer2.0 del dataset for input.
*     ----------------------------------------
	ludel = luavail()
	call opnieee(ludel,filen(1:lf),'old',ierr)

	if(ierr.ne.0) then
	  write(stderr,'(4a,i4)') myname,': error return from opnieee(',
     &	    filen(1:lf),'), iostat =',ierr
	  ierr=1
	  return
	endif
	write(stdout,'(4a)') myname,': filename = "',filen(1:lf),'"'

*     Read header record from del dataset.
*     ------------------------------------

	rewind(ludel)
	read(ludel,iostat=ierr) iobs,lenrecs,iymd,ihms
	if(ierr.eq.0) then
	  if(nymd.eq.0.and.nhms.eq.0) then
	    iy2kmd=iymd
	    if(iymd.lt.100) iy2kmd=1900+iymd
	    nymd=iy2kmd
	    nhms=ihms
	  endif
	endif

	do while(ierr.eq.0 .and. nymd.gt.0.and.nhms.ge.0 .and.
     &	  iobs.ge.0.and.iy2kmd.ge.0.and.ihms.ge.0 .and.
     &	  (iy2kmd.lt.nymd .or. (iy2kmd.eq.nymd.and.ihms.lt.nhms)))

	  write(stdout,'(2a,i8.8,a,i6.6,a,i6)') myname,
     &	    ': skipping nymd = ',iy2kmd,',nhms = ',ihms,', nobs = ',iobs

	  do irec=1,iobs		! skip all record
	    if(ierr.eq.0) read(ludel,iostat=ierr)
	  end do

	  if(ierr.ne.0) then
	    write(stderr,'(4a,i4)') myname,
     &        ': corrupted input "',filen(1:lf),'", iostat =',ierr
	    ierr=2
	    return

	  else
	    read(ludel,iostat=ierr) iobs,lenrecs,iymd,ihms
	    if(ierr.eq.0) then
	      iy2kmd=iymd
	      if(iymd.lt.100) iy2kmd=1900+iymd
	    endif
	  endif
	end do

c	..Found the dataset?
c	=====================
	if(ierr.ne.0) then

	  write(stderr,'(3a,i8.8,a,i6.6,a,i4)') myname,
     &	    ': not found expected date-time',
     &	    ', nymd = ',nymd,', nhms = ',nhms,', iostat = ',ierr
	  ierr=-1
	  return

	elseif(iobs.lt.0.or.iy2kmd.lt.0.or.ihms.lt.0) then
	  write(stderr,'(4a,i4)') myname,
     &      ': corrupted input "',filen(1:lf),'", iostat =',ierr
	  ierr=2
	  return

	elseif(iy2kmd.gt.nymd.or.(iy2kmd.eq.nymd.and.ihms.gt.nhms)) then

	  write(stderr,'(2a,2(a,i8.8,a,i6.6))') myname,
     &	    ': not found expected date-time',
     &	    ', nymd = ',nymd,', nhms = ',nhms,
     &	    ', iymd = ',iy2kmd,', ihms = ',ihms
	  ierr=2
	  return

	elseif(nymd.le.0.and.nhms.le.0) then
	  nymd=iy2kmd
	  nhms=ihms

	endif

	write(stdout,'(2a,i8.8,a,i6.6,a,i6)') myname,
     &	  ': processing nymd = ',nymd,',nhms = ',nhms,', nobs = ',iobs

*     Read the data from the del dataset.  One rec. for one obs.
*     ----------------------------------------------------------
	ierr=0
	irec=0
	nobs=0
	do while(ierr.eq.0.and.irec.lt.iobs)

c	  ..Read a record
	  irec=irec+1
	  read(ludel,iostat=ierr)
     &	    xlat,xlon,xlev,ikx,ikt,xdel,xsigO,xsigF,xtim,xqcf

	  mqcf=0
	  if(xqcf<1.) mqcf=floor((xqcf-9.)*10000.)

c	  ..Preliminary selections
	  if(.not.exclude .or. mqcf==0) then
	    nobs=nobs+1
	    if(nobs.le.nobsmax) then
	      rlats(nobs)=xlat
	      rlons(nobs)=xlon
	      rlevs(nobs)=xlev
	      kx(nobs)=ikx
	      kt(nobs)=ikt
	      del(nobs)=xdel
	      tstamp(nobs)=xtim
	      qcflag(nobs)=mqcf
	    endif
	  endif
	end do

c	..Input error?

	if(ierr.ne.0) then
	  write(stderr,'(4a,i4)') myname,
     &	    ': corrupted input "',filen(1:lf),'", iostat =',ierr
	  ierr=2
	  return
	endif

c	..Out of the memory?

	if(nobs.gt.nobsmax) then
	  write(stderr,'(2a,2(i6,a))') myname,
     &	    ': insufficient workspace, (nobs=',
     &	    nobs,') > (nobsmax=',
     &	    nobsmax,')'
	  ierr=3
	  return
	endif

c	..A summary of the input
	if(nobs.lt.iobs) write(stdout,'(2a,2(i6,a))') myname,
     &	    ': selected ',nobs,' observations out of ',iobs,' records'

*     Print selected records from del dataset.
*     ----------------------------------------
      if( (verbose).and.(idelpri.ne.0) ) then
        do irec = idelprb, idelpre, idelpri
	  if(irec.ge.1.and.irec.le.nobs) then
            write( stdout, 9001 ) irec,
     $	      rlats(irec), rlons(irec), rlevs(irec), kx(irec), kt(irec),
     $	      del(irec),tstamp(irec),qcflag(irec)
	  endif
	end do
 9001   format(1x,i6,3f8.2,2i3,1x,f9.4,f8.2,i6)
      endif


*     Close the oiVer2.0 del dataset.
*     ------------------------------
      close(ludel)


*     All done
*     --------
      end
