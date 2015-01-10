
c  subroutine smex( a, n, ii, m, b )
c  Extracts the symmetric m by m submatrix b of the symmetric n by n 
c  matrix a located in rows ii through ii+m-1 and columns ii through
c  ii+m-1. Both a and b are stored in compressed form.

      subroutine smex( a, n, ii, m, b )

c  Tested with symslvt.f

      real              a(n*(n+1)/2)
      real              b(m*(m+1)/2)

      iim = ii - 1
      ib = 1
      ia = ii*(ii+1)/2 - 1
      do 110 j = 1, m
         do 100 i = 1, j
            b(ib) = a(ia+i)
            ib = ib + 1
  100    continue
         ia = ia + iim + j
  110 continue

      return
      end
