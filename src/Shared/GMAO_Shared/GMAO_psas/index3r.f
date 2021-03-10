
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Indexes an array arrin of length n, i.e. outputs the array indx
c.... such that arrin(indx(j)) is in ascending order for j=1,2,3,...,n.
c.... The input quantities arrin and n are not changed.

c.... ref: Numerical recipes p.233

c History:
c	02Feb95 - Jing G.	- Added if(n.eq.0) ..

      subroutine index3r( n, ir1, ar1, ar2, ar3, indx )

c.......................................................................
c.... Argument declarations.

      integer      n
      integer      ir1(n)
      real         ar1(n)
      real         ar2(n)
      real         ar3(n)
      integer      indx(n)

c.......................................................................
c.... Local storage.

      integer      i1
      real         q1
      real         q2
      real         q3

c.......................................................................
c.... Statement function definition.

      logical altb

      altb(ia1,a1,a2,a3,ib1,b1,b2,b3) =  
     $      (ia1.lt.ib1).or.
     $     ((ia1.eq.ib1).and.(a1.lt.b1)).or.
     $     ((ia1.eq.ib1).and.(a1.eq.b1).and.(a2.lt.b2)).or.
     $     ((ia1.eq.ib1).and.(a1.eq.b1).and.(a2.eq.b2).and.(a3.lt.b3))

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
            i1 = ir1(indxt)
            q1 = ar1(indxt)
            q2 = ar2(indxt)
            q3 = ar3(indxt)

         else

            indxt = indx(ir)
            i1 = ir1(indxt)
            q1 = ar1(indxt)
            q2 = ar2(indxt)
            q3 = ar3(indxt)
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
               ma = indx(j)
               mb = indx(j+1)
               if( altb(ir1(ma),ar1(ma),ar2(ma),ar3(ma),
     $                  ir1(mb),ar1(mb),ar2(mb),ar3(mb)) ) j = j + 1
            endif

            mb = indx(j)

            if( altb(i1,q1,q2,q3,ir1(mb),ar1(mb),ar2(mb),ar3(mb)) ) then
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
