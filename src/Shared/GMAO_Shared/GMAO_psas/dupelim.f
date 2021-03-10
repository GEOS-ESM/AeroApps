c                                                   01/13/94 - dupelim.f
c  07oct94 - A. da S. - introduced tolerance (TOL) to compare lat/lon/lev
c                       this fix needed because CRAY was missing duplicates
c                       because of its IEEE conversion to 64 bits.
c  05oct94 - A. da S. - eliminated 'nobs' from parameter list; now
c                       nnobs is reset if duplicates are found.
c
c  1/23 - pass 'kl' flag to main routine -mes-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Eliminate duplicate observations and readjust region pointers

      subroutine DUPELIM (verbose, luverb,  nnobs,   kx,    kt,   kl,
     $                    rlats,   rlons,   rlevs,   del,sigU,sigO,sigF,
     $                    tstamp,  maxreg,  iregbeg, ireglen, 
     $                    ktmax, ityplen )

	use m_stdio, only : stdout

*     Argument declarations.
*     ----------------------
      logical      verbose
      integer      luverb
      integer      nnobs
      integer      kx(nnobs)
      integer      kt(nnobs)
      logical      kl(nnobs)
      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      real         del(nnobs)
      real         sigU(nnobs)
      real         sigO(nnobs)
      real         sigF(nnobs)
      real         tstamp(nnobs)
c-    real         qcflag(*)		! no longer used
      integer      maxreg
      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ktmax
      integer      ityplen(ktmax,maxreg)


*     Tolerance to define a duplicate
*     NOTE: This is important because of CRAY's IEEE conversion
*     ---------------------------------------------------------
      parameter ( tol = 1.E-4 )


*     Local storage.
*     --------------
      logical      kflag

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  input at this pt is sorted by region, type, kx, lat, lon, and level
c  scan through all regions 
c  mark kl as .true. if it does not duplicate previous ob.
c  mark duplicates as .false., and decrement 'ityplen', 'ireglen'
c  (ireglen of current region, ityplen of type kx of that region)
c
c  when completed, call tofront to shift data around
c  then redo 'iregbeg' count (code from regsort)
c.......................................................................

      itot = 0
      nobs = nnobs

*     Loop over the regions.
*     ----------------------
      do 500 ireg = 1, maxreg


*        Check whether there is data in the current region.
*        --------------------------------------------------
         if( ireglen(ireg).eq.0 ) goto 500

         ibeg = iregbeg(ireg)
         ilen = ireglen(ireg)

         kl(ibeg) = .true.

         do 200 n = ibeg+1, ibeg+ilen-1
            n1 = n - 1
            kflag = ( (abs(rlevs(n)-rlevs(n1)) .gt. TOL) .or.
     $                (abs(rlons(n)-rlons(n1)) .gt. TOL) .or.
     $                (abs(rlats(n)-rlats(n1)) .gt. TOL) .or.
     $                kx(n)    .ne. kx(n1)    .or.
     $                kt(n)    .ne. kt(n1) )

            if ( .not. kflag ) then
               ireglen(ireg) = ireglen(ireg) - 1
               ityplen(kt(n),ireg) = ityplen(kt(n),ireg) - 1
               itot = itot + 1
            endif
            kl(n) = kflag

  200    continue

  500 continue


*     Report number of duplicate found
*     --------------------------------
      if ( itot .eq. 0 )  then
          return
      end if


*     Slide the non-duplicate data to the front of the arrays.
*     -------------------------------------------------------
      call TOFRONT ( nnobs, kx,    kt,    kl,   
     $               rlats, rlons, rlevs,
     $               del,   sigU,sigO,  sigF,  tstamp, nobs  )


      if ( verbose ) then
          write(luverb, *) '   dupelim: ', itot,' duplicate obs found'
          write(luverb, *) '   dupelim: nobs reset from ', nnobs,
     .                     ' to ', nobs
          if ( luverb .ne. stdout ) then
            write(stdout, *) '   dupelim: ', itot,' duplicate obs found'
            write(stdout, *) '   dupelim: nobs reset from ', nnobs,
     .                     ' to ', nobs
          end if
      end if


*     Reset number of observations
*     ----------------------------
      nnobs = nobs


*     Update region pointers
*     ----------------------
      ksofar = 0
      do 600 ireg = 1, maxreg
         if( ireglen(ireg).gt.0 ) then
            iregbeg(ireg) = ksofar + 1
            ksofar = ksofar + ireglen(ireg)
         endif
  600 continue

*     All done.
*     ---------
      return
      end
