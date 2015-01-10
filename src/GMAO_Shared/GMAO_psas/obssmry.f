
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Print a summary of all the good observations.

      subroutine OBSSMRY ( lu, nobs, kx, ktype )
      use m_ioutil,only : luflush

c.......................................................................
c.... Argument declarations.

      integer      lu
      integer      nobs
      integer      kx(nobs)
      integer      ktype(nobs)

c.......................................................................
c.... Data source names.

      include     "ktmax.h"
      include     "kxmax.h"
      include     "kxtabl.h"

c.......................................................................
c.... Data type names.

      include     "kttabl.h"

c.......................................................................
c.... Local storage for counters.

      integer      ksums(kxmax,ktmax)
      integer      ktots(ktmax)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.......................................................................
c.... Initialize counters.

      do 110 j = 1, ktmax
         ktots(j) = 0
         do 100 i = 1, kxmax
            ksums(i,j) = 0
  100    continue
  110 continue

c.......................................................................
c.... Count up the data by source and type.

      do 200 n = 1, nobs
         ksums(kx(n),ktype(n)) = ksums(kx(n),ktype(n)) + 1
         ktots(ktype(n))       = ktots(ktype(n))       + 1
  200 continue

c.......................................................................
c.... Write the observation summary sub-tables for each data source.

      do 300 i = 1, kxmax

         ksum = 0
         do 250 j = 1, ktmax
            ksum = ksum + ksums(i,j)
  250    continue

         if( ksum.gt.0 ) then

            write(lu,900) i, kxdesc(i)
  900       format('  obssmry:  kx    = ',I3,'  kxdesc = ',A25)

            do 275 j = 1, ktmax

               if( ksums(i,j).ne.0 ) then
                  write(lu,910) j, ktname(j), ksums(i,j)
  910             format('  obssmry:  ktype = ',I3,'  ktname = ',A9,I10)
               endif

  275       continue
            write(lu,920)
  920       format('  obssmry:')

         endif

  300 continue

c.......................................................................
c.... Write the observation summary table for all data sources.

      kgrand = 0
      write(lu,900) 9999, 'ALL KXs '

      do 325 j = 1, ktmax

         kgrand = kgrand + ktots(j)
         write(lu,910) j, ktname(j), ktots(j)

  325 continue

      write(lu,920)
      write(lu,910) 9999, 'ALL KTs ', kgrand
      write(lu,920)
      call luflush(lu)

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
