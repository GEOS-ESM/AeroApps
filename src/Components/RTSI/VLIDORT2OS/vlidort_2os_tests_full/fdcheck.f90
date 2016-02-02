program fdcheck
implicit none
integer :: n, nd
real*8  :: grad,pert,bas
open(34,file='fort.34',status='old')
open(35,file='fort.35',status='old')
do n  = 1, 60
   read(34,*)nd,grad,bas ; read(35,*)nd,pert
   write(*,*)nd,grad,(pert-bas)/0.001d0
enddo
close(34) ; close(35)
stop
end program fdcheck
