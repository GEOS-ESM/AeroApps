
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Indexes an array arrin of length n, i.e. outputs the array indx
c.... such that arrin(indx(j)) is in ascending order for j=1,2,3,...,n.
c.... The input quantities arrin and n are not changed.

c.... ref: Numerical recipes p.233

c History:
c	02Feb95 - Jing G.	- Added if(n.eq.0) ..

      subroutine indexxi( n, arrin, indx )

c.......................................................................
c.... Argument declarations.

      integer      n
      integer      arrin(n)
      integer      indx(n)

c.......................................................................
c.... Local storage.

      integer      q

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(n.eq.0) return		! fool proved condition checking

      do 11 j = 1, n
         indx(j) = j
   11 continue

      if( n.eq.1 ) return

      L = n/2 + 1
      ir = n

   10 continue

         if( L.gt.1 ) then
            L = L - 1
            indxt = indx(L)
            q = arrin(indxt)

         else
            indxt = indx(ir)
            q = arrin(indxt)
            indx(ir) = indx(1)
            ir = ir - 1
            if( ir.eq.1 ) then
               indx(1) = indxt
               return
            endif

         endif

         i = L
         j = L + L

   20    if( j.le.ir ) then

            if( j.lt.ir ) then
               if( arrin(indx(j)).lt.arrin(indx(j+1)) ) j = j + 1
            endif

            if( q.lt.arrin(indx(j)) ) then
               indx(i) = indx(j)
               i = j
               j = j + j

            else
               j = ir + 1

            endif

            goto 20

         endif

         indx(i) = indxt

      goto 10

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
