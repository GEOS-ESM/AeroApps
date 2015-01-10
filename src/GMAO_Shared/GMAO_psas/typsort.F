C  16Nov01  - Todling   - Added TYPSORT_ALL
C  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
C			  to follow the guide lines.
c  03oct94 - Implemented CRAY specifics with IFDEFs.
c  modification for dynamic storage on CRRAY                    05/28/93
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Sort data in regions by kt-kx-rlats-rlons-rlevs and set ityplen.

      subroutine TYPSORT ( verbose, luverb,  nnobs,   kx,    kt,
     $                    rlats,   rlons,   rlevs,   del,
     &			  sigU,    sigO,    sigF,
     $                    tstamp,  maxreg,  iregbeg, ireglen, 
     $                    ktmax, ityplen                               )

	use m_die,only : die
c.......................................................................
c.... Argument declarations.
	implicit none

      logical      verbose
      integer      luverb
      integer      nnobs
      integer      kx(nnobs)
      integer      kt(nnobs)
      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      real         del(nnobs)
      real         sigU(nnobs)
      real         sigO(nnobs)
      real         sigF(nnobs)
      real         tstamp(nnobs)
c-    real         qcflag(nnobs)	! no longer used
      integer      maxreg
      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ktmax
      integer      ityplen(ktmax,maxreg)

c.......................................................................
c.... Local storage.

	character(len=*), parameter :: myname='typsort'

	integer,allocatable	:: iperm(:)	! (nnobs)
	integer istat,ireg,ktt,ibeg,ilen,n

!-----------------------------------------------------------------------

	allocate(iperm(nnobs), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

c.......................................................................
c.... Loop over the regions.

      do 500 ireg = 1, maxreg

c.......................................................................
c....... Initialize array ityplen.

         do 100 ktt = 1, ktmax
            ityplen(ktt,ireg) = 0
  100    continue

c.......................................................................
c....... Check whether there is data in the current region.

         if( ireglen(ireg).eq.0 ) goto 500

c.......................................................................
c....... Count the data of each type in the current region.

         ibeg = iregbeg(ireg)
         ilen = ireglen(ireg)

         do 200 n = ibeg, ibeg+ilen-1
            ityplen(kt(n),ireg) = ityplen(kt(n),ireg) + 1
  200    continue

c.......................................................................
c....... Sort the data in the current region by type.

         call INDEXXI ( ilen, kt(ibeg), iperm )

c.......................................................................
c....... Apply the sorting permutation to all data items for the region.

         call PERMUTI ( kx(ibeg),    iperm, ilen, kx(ibeg)    )
         call PERMUTI ( kt(ibeg),    iperm, ilen, kt(ibeg)    )
         call PERMUTR ( rlats(ibeg), iperm, ilen, rlats(ibeg) )
         call PERMUTR ( rlons(ibeg), iperm, ilen, rlons(ibeg) )
         call PERMUTR ( rlevs(ibeg), iperm, ilen, rlevs(ibeg) )
         call PERMUTR ( del(ibeg),   iperm, ilen, del(ibeg)   )
         call PERMUTR ( sigU(ibeg),  iperm, ilen, sigU(ibeg)  )
         call PERMUTR ( sigO(ibeg),  iperm, ilen, sigO(ibeg)  )
         call PERMUTR ( sigF(ibeg),  iperm, ilen, sigF(ibeg)  )
         call PERMUTR ( tstamp(ibeg),iperm, ilen, tstamp(ibeg))
c-       call PERMUTR ( qcflag(ibeg),iperm, ilen, qcflag(ibeg))

c.......................................................................
c....... Loop over the data types.

         do 300 ktt = 1, ktmax

c.......................................................................
c.......... Check whether there is any type = kt data in current region.

            ilen = ityplen(ktt,ireg)
            if( ilen.eq.0 ) goto 300

c.......................................................................
c.......... Sort the data by kx-lat-lon-lev.

            call INDEX3R ( ilen, kx(ibeg), rlats(ibeg), rlons(ibeg), 
     $                                    rlevs(ibeg), iperm        )

c.......................................................................
c.......... Apply the sorting permutation to the data items.

            call PERMUTI ( kx(ibeg),    iperm, ilen, kx(ibeg)    )
            call PERMUTR ( rlats(ibeg), iperm, ilen, rlats(ibeg) )
            call PERMUTR ( rlons(ibeg), iperm, ilen, rlons(ibeg) )
            call PERMUTR ( rlevs(ibeg), iperm, ilen, rlevs(ibeg) )
            call PERMUTR ( del(ibeg),   iperm, ilen, del(ibeg)   )
            call PERMUTR ( sigU(ibeg),  iperm, ilen, sigU(ibeg)  )
            call PERMUTR ( sigO(ibeg),  iperm, ilen, sigO(ibeg)  )
            call PERMUTR ( sigF(ibeg),  iperm, ilen, sigF(ibeg)  )
            call PERMUTR ( tstamp(ibeg),iperm, ilen, tstamp(ibeg))
c-          call PERMUTR ( qcflag(ibeg),iperm, ilen, qcflag(ibeg))

            ibeg = ibeg + ilen

  300    continue

  500 continue

      deallocate(iperm)
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
      Subroutine typsort_all(verbose,luverb,nnobs,kx,kt,rlats,rlons,
     .                       rlevs,
     .                       kid,maxreg,iregbeg,ireglen,ktmax,ityplen)
c******************************************************************************
c*                                                                            *
c*                                                                            *
c*                                                                            *
c*                                                                            *
c******************************************************************************
      implicit none

      integer luverb
      integer nnobs
      integer maxreg
      integer ktmax

      real rlats(nnobs)
      real rlons(nnobs)
      real rlevs(nnobs)

      integer kx(nnobs)
      integer kt(nnobs)
      integer kid(nnobs)
      integer iperm(nnobs)
      integer iregbeg(maxreg)
      integer ireglen(maxreg)
      integer ityplen(ktmax,maxreg)

      logical verbose


      integer ktt, ireg, n, ibeg, ilen

c******************************
c*    Loop over the regions   *
c******************************
cc$omp parallel do private(ireg,ktt,ibeg,ilen,n,iperm)
      do 500 ireg=1,maxreg

c-----Initialize array ityplen.
      do 100 ktt=1,ktmax
      ityplen(ktt,ireg)=0
  100 continue

c-----Check whether there is data in the current region.
      if(ireglen(ireg).eq.0 ) goto 500

c-----Count the data of each type in the current region.
      ibeg=iregbeg(ireg)
      ilen=ireglen(ireg)

      do 200 n=ibeg,ibeg+ilen-1
      ityplen(kt(n),ireg)=ityplen(kt(n),ireg)+1
  200 continue

c-----Sort the data in the current region by type.
      call INDEXXI(ilen,kt(ibeg),iperm)

c-----Apply the sorting permutation to all data items for the region.
      call PERMUTI(kid(ibeg),iperm,ilen,kid(ibeg))
      call PERMUTI(kx(ibeg),iperm,ilen,kx(ibeg))
      call PERMUTI(kt(ibeg),iperm,ilen,kt(ibeg))
      call PERMUTR(rlats(ibeg),iperm,ilen,rlats(ibeg))
      call PERMUTR(rlons(ibeg),iperm,ilen,rlons(ibeg))
      call PERMUTR(rlevs(ibeg),iperm,ilen,rlevs(ibeg))

c-----Loop over the data types.
      do 300 ktt=1,ktmax

c-----Check whether there is any type = kt data in current region.
      ilen=ityplen(ktt,ireg)
      if(ilen.eq.0) goto 300

c-----Sort the data by kx-lat-lon-lev.
      call INDEX3R(ilen,kx(ibeg),rlats(ibeg),rlons(ibeg),
     .             rlevs(ibeg),iperm)

c-----Apply the sorting permutation to the data items.
      call PERMUTI(kid(ibeg),iperm,ilen,kid(ibeg))
      call PERMUTI(kx(ibeg),iperm,ilen,kx(ibeg))
      call PERMUTR(rlats(ibeg),iperm,ilen,rlats(ibeg))
      call PERMUTR(rlons(ibeg),iperm,ilen,rlons(ibeg))
      call PERMUTR(rlevs(ibeg),iperm,ilen,rlevs(ibeg))
      ibeg = ibeg+ilen
  300 continue

  500 continue

      return
      end
