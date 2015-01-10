
c    Shift the diagonal elements of "a" by a constant "del", and
c    renormalized.  "a" is a symmentric matrix in compressed form.

      subroutine smexsh( a, m, del )

c  Tested with symslvt.f

      real              a(m*(m+1)/2)
      real              del

      denom=1.+del
      ia=0
      do j=1,m
	do i=1,j-1
	  ia=ia+1
	  a(ia)=a(ia)/denom
	end do
	ia=ia+1
      end do

      end
