cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Restrict the choice to data within a union of lat/lon boxes
c.... i.e. we want to set kl(n) = kl(n) and (box1.or.box2.or....or.boxk)

      subroutine llboxes( verbose, luverb, nboxes, boxes,
     &	nobs, rlons, rlats, rlevs, kx, kt, tstamp, kl )

	use config, only : ktslp,ktus,ktvs
c.......................................................................
c.... Argument declarations.

      logical      verbose
      integer      luverb
      integer      nboxes
      real         boxes(2,6,nboxes)
      integer      nobs
      real         rlons(nobs)
      real         rlats(nobs)
      real         rlevs(nobs)
      integer	   kx(nobs)
      integer	   kt(nobs)
      real	   tstamp(nobs)
      logical      kl(nobs)

c.......................................................................
c.... Local storage

      logical      inunion

c.......................................................................
c.... Statement functions

      real rlondn,rlonup
      real rlatdn,rlatup
      real rlevdn,rlevup
      integer kxdn,kxup
      integer ktdn,ktup
      real tmdn,tmup

      kxdn(k)   = nint(boxes(1,1,k))	! Data source
      kxup(k)   = nint(boxes(2,1,k))	! Data source
      ktdn(k)   = nint(boxes(1,2,k))	! Data type
      ktup(k)   = nint(boxes(2,2,k))	! Data type
      rlatdn(k) = boxes(1,3,k)		! South
      rlatup(k) = boxes(2,3,k)		! North
      rlondn(k) = boxes(1,4,k)		! West
      rlonup(k) = boxes(2,4,k)		! East
      rlevdn(k) = boxes(2,5,k)		! Low (high pressure)
      rlevup(k) = boxes(1,5,k)		! High (low pressure)
      tmdn(k)	= boxes(1,6,k)		! ealiest time stamp
      tmup(k)	= boxes(2,6,k)		! latest time stamp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if(nboxes.eq.0) then
	  write(luverb,'(a)') 'llboxes: no box, all obs. are included'
	  return	! if no box, default is all in
	endif

	kin=0
	do n=1,nobs
	  if(kl(n)) kin=kin+1
	end do

	do n = 1, nobs

          if( kl(n) ) then

	    k=0
	    inunion=.false.
	    do while(.not.inunion.and.k.lt.nboxes)
	      k=k+1
	      inunion=
     &		rlondn(k).le. rlons(n).and. rlons(n).le.rlonup(k) .and.
     &		rlatdn(k).le. rlats(n).and. rlats(n).le.rlatup(k) .and.
     &		  kxdn(k).le.    kx(n).and.    kx(n).le.  kxup(k) .and.
     &		  ktdn(k).le.    kt(n).and.    kt(n).le.  ktup(k) .and.
     &		  tmdn(k).le.tstamp(n).and.tstamp(n).le.  tmup(k)

	      if(inunion .and. kt(n).ne.ktslp .and. kt(n).ne.ktus .and.
     &			       kt(n).ne.ktvs ) inunion =
     &		rlevup(k).le. rlevs(n).and. rlevs(n).le.rlevdn(k)
	    end do

	    kl(n)=inunion	! update kl
	  endif

	end do

	kout=0
	do n=1,nobs
	  if(kl(n)) kout=kout+1
	end do

	if(verbose) write(luverb,920) kin, kout
  920 format('  llboxes: from ',I6,' observations, ',I6,' remain'/
     $       '  llboxes:')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
