        subroutine QVMV ( a, b, c, n )

        real a(n), b(n), c(n)

*
*       Sets a(i) = b(i) * c(i), i = 1, n
*

        do 10 i = 1, n
           a(i) = b(i) * c(i)
10      continue


        return
        end


*.......................................................................

        subroutine QVDV ( a, b, c, n )

        real a(n), b(n), c(n)

*
*       Sets a(i) = b(i) / c(i), i = 1, n
*

        do 10 i = 1, n
           a(i) = b(i) / c(i)
10      continue


        return
        end

	subroutine QVES ( a, n, s )
	real a(n)

c	..Set a(i)=s

	do i=1,n
	  a(i)=s
	end do
	end
