
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Print summary of data in selected region.

      subroutine prtregn( verbose, luverb,  idreg,   prtdata,
     $                    maxreg,  iregbeg, ireglen, ityplen,
     $                    nnobs,   kx,      kt, 
     $                    rlats,   rlons,   rlevs,
     $                    del,     sigO,    sigF,    tstamp )

c.......................................................................
c.... Argument declarations.

      include     "ktmax.h"

      logical      verbose
      integer      luverb
      integer      idreg
      logical      prtdata
      integer      maxreg
      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ityplen(ktmax,maxreg)
      integer      nnobs
      integer      kx(nnobs)
      integer      kt(nnobs)
      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      real         del(nnobs)
      real         sigO(nnobs)
      real         sigF(nnobs)
      real         tstamp(nnobs)
c-    real         qcflag(*)		! no longer used

c.......................................................................
c.... Data source names.

      include     "kxmax.h"
      include     "kxtabl.h"

c.......................................................................
c.... Data type names.

      include     "kttabl.h"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if( verbose ) then

         write(luverb,910) idreg, iregbeg(idreg), ireglen(idreg),
     $                            (ityplen(k,idreg),k=1,ktmax)
  910    format('  prtregn:  region= ',I2,', iregbeg=',I8,
     $                                              ', ireglen=',I8,/
     $         ('  prtregn: ityplen=',8I8))

         if( prtdata ) then

            ibeg = iregbeg(idreg)
            do 200 k = 1, ktmax

               len = ityplen(k,idreg)
               if( len.eq.0 ) goto 200

               write(luverb,920)  ktname(k)
  920          format('  prtregn: ',A8)

	       do n=ibeg,ibeg+len-1
                 write(luverb,930) kxdesc(kx(n)), rlats(n), rlons(n),
     $                             rlevs(n), del(n), sigO(n), sigF(n)
c-   $                             tstamp(n), qcflag(n)
  930            format('  prtregn: ',9X,A25,F6.2,F7.2,F9.2,3G13.4,F7.1,
     $                F7.4)
	       end do

               ibeg = ibeg + len

  200       continue

         endif

         write(luverb,900)
  900    format('  prtregn:')

      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
