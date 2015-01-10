
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Print out selected data.

      subroutine prtdata( verbose, luverb, nobs,    kx,    kt,    kl,
     $                    rlats,   rlons,  rlevs,   del,   sigO,  sigF,
     $                    tstamp                                       )

c.......................................................................
c.... Argument declarations.

      logical      verbose
      integer      luverb
      integer      nobs
      integer      kx(nobs)
      integer      kt(nobs)
      logical      kl(nobs)
      real         rlats(nobs)
      real         rlons(nobs)
      real         rlevs(nobs)
      real         del(nobs)
      real         sigO(nobs)
      real         sigF(nobs)
      real         tstamp(nobs)
c-    real         qcflag(*)		! not used

c.......................................................................
c.... Data source names.

      include     "ktmax.h"
      include     "kxmax.h"
      include     "kxtabl.h"

c.......................................................................
c.... Data type names.

      include     "kttabl.h"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if( verbose ) then

         kount = 0
         do 100 n = 1, nobs

            if( kl(n) ) then
               kount = kount + 1
               write(luverb,920) kount, kxdesc(kx(n)), ktname(kt(n)),
     $                           rlats(n), rlons(n), rlevs(n), del(n),
     $                           sigO(n),  sigF(n)
					!!!! , tstamp(n), qcflag(n)
  920          format('  prtdata:',I7,A26,A9,F7.2,F8.2,F8.2,3G13.5,F7.1,
     $                F7.4)
            endif

  100    continue

         write(luverb,930)
  930    format('  prtdata:')

      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
