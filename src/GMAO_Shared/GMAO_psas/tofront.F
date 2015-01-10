!  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
!			  to follow the guide lines.
!  03oct94 - Implemented CRAY specifics with IFDEFs.
!  modification for dynamic storage on CRAY                     05/28/93
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.... Slide the selected data to the front of the data arrays.

      subroutine TOFRONT ( nobs,  kx,    kt,    kl,   
     $                     rlats, rlons, rlevs,
     $                     del,   sigU,   sigO,  sigF,  tstamp, newnr )

	use m_die,only : die
!.......................................................................
!.... Argument declarations.
	implicit none

      integer      nobs
      integer      kx(nobs)
      integer      kt(nobs)
      logical      kl(nobs)
      real         rlats(nobs)
      real         rlons(nobs)
      real         rlevs(nobs)
      real         del(nobs)
      real         sigU(nobs)
      real         sigO(nobs)
      real         sigF(nobs)
      real         tstamp(nobs)
!-    real         qcflag(nobs) 	! no longer used
      integer      newnr

!.......................................................................
!.... Local storage

	character(len=*), parameter :: myname='tofront'
	real flip
	integer istat,n

	integer,allocatable :: iperm(:)
	real,   allocatable :: rsort(:)

	allocate(iperm(nobs),rsort(nobs), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

!.......................................................................
!.... Initialize the vector to be sorted

      newnr = 0

      do 100 n = 1, nobs

         if( kl(n) ) then
            newnr    = newnr + 1
            rsort(n) = - float(newnr)
         else
            rsort(n) = float(n)
         endif

  100 continue

      if( newnr.gt.0 ) then

!.......................................................................
!....... Perform the sort.

         flip = - float(newnr)
         do 110 n = 1, nobs
            if( kl(n) ) rsort(n) = flip - rsort(n)
  110    continue

         call INDEXXR ( nobs, rsort, iperm )

!.......................................................................
!....... Apply the sorting permutation to all arrays.

         call PERMUTI ( kx,    iperm, nobs, kx    )
         call PERMUTI ( kt,    iperm, nobs, kt    )
         call PERMUTL ( kl,    iperm, nobs, kl    )
         call PERMUTR ( rlats, iperm, nobs, rlats )
         call PERMUTR ( rlons, iperm, nobs, rlons )
         call PERMUTR ( rlevs, iperm, nobs, rlevs )
         call PERMUTR ( del,   iperm, nobs, del   )
         call PERMUTR ( sigU,  iperm, nobs, sigU  )
         call PERMUTR ( sigO,  iperm, nobs, sigO  )
         call PERMUTR ( sigF,  iperm, nobs, sigF  )
         call PERMUTR ( tstamp,iperm, nobs, tstamp)
!-       call PERMUTR ( qcflag,iperm, nobs, qcflag)

      endif

      deallocate(iperm,rsort)
      return

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
