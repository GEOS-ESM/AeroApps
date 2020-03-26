
c  Sets pa(i)=a(iperm(i)) i=1,2,...,n for real a and pa.
c  The vectors a and pa can share the same storage locations

      subroutine permutr( a, iperm, n, pa )

      real       a(n), pa(n), asave
      integer    iperm(n)

      if( n.gt.0 ) then

         do 100 i = 1, n
            iperm(i) = -iperm(i)
  100    continue

         kbeg = 0
  200    kbeg = kbeg + 1

         if( kbeg.le.n ) then

            if( iperm(kbeg).lt.0 ) then

               asave = a(kbeg)
               k = kbeg
  300          l = -iperm(k)

               if( l.eq.kbeg ) then
                  pa(k) = asave
                  iperm(k) = l
                  goto 200
               else
                  pa(k) = a(l)
                  iperm(k) = l
                  k = l
                  goto 300
               endif

            else
               goto 200

            endif

         endif

      endif

      return
      end
