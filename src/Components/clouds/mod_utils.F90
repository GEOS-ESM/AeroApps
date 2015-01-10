! file: mod_utils.f90
! various basic utilities

module mod_utils

  implicit none

  logical :: exception = .FALSE.

contains

  subroutine myError (message)
    character(*), intent(in) :: message
    write(*,'(a)') 'Error: ', message
    exception = .TRUE.
  end subroutine myError

  ! given a lower trangular matrix A(N,N) representing a symmetrix matrix
  ! copy the off-diagonal lower elements to the upper off-diagonal elements
  ! to make A actually symmetric.

  subroutine L2U(N,A)
    integer, intent(in)    :: N
    real*4,  intent(inout) :: A(N,N)
    integer :: I, J
    do J = 1, N-1
      do I = J+1, N
        A(J,I) = A(I,J)
      end do
    end do
  end subroutine L2U

  ! given a upper trangular matrix A(N,N) representing a symmetrix matrix
  ! copy the off-diagonal upper elements to the lower off-diagonal elements
  ! to make A actually symmetric.

  subroutine U2L(N,A)
    integer, intent(in)    :: N
    real*4,  intent(inout) :: A(N,N)
    integer :: I, J
    do J = 1, N-1
      do I = J+1, N
        A(I,J) = A(J,I)
      end do
    end do
  end subroutine U2L

end module mod_utils
