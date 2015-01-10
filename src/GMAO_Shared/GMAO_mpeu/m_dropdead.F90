!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_dropdead - An abort() with a style
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_dropdead
      implicit none
      private	! except

      public	:: die	! terminate a program with a condition

      interface die; module procedure	&
	die_,	&
	diex_
      end interface

! !REVISION HISTORY:
! 	20Feb97 - Jing Guo <guo@eramus> - defined template
!EOP
!_______________________________________________________________________

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: die_ - Clean up and raise an exception to the OS
!
! !DESCRIPTION:
!
!   A call to die() exits the program with minimum information for
!   both the user and the operating system.
!
! !INTERFACE:

    subroutine die_(where)
      use m_stdio, only : stderr
      use m_mpif90,only : MP_comm_world
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_abort
      implicit none
      character(len=*),intent(in) :: where	! where it is called

! !REVISION HISTORY:
! 	20Feb97 - Jing Guo <guo@eramus> - defined template
!
!EOP
!_______________________________________________________________________

  character(len=*),parameter :: myname_='die.'
  integer :: myrank,ier

	!-------------------------------------------------
	! MPI_ should have been initialized for this call
	!-------------------------------------------------

    call MP_comm_rank(MP_comm_world,myrank,ier)

	! a message for the users:

    write(stderr,'(z3.3,5a)') myrank,'.',myname_,	&
      ': from ',trim(where),'()'

	! raise a condition to the OS

    call MP_abort(MP_comm_world,2,ier)

end subroutine die_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: diex_ - Clean up and raise an exception to the OS
!
! !DESCRIPTION:
!
!   A call to die() exits the program with minimum information for
!   both the user and the operating system.  This implementation,
!   however, may be used in conjunction with with a source preprocessor
!   to produce more detailed location information.
!
! !INTERFACE:

    subroutine diex_(where,fnam,line)
      use m_stdio, only : stderr
      use m_mpif90,only : MP_comm_world
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_abort
      implicit none
      character(len=*),intent(in) :: where	! where it is called
      character(len=*),intent(in) :: fnam
      integer,intent(in) :: line

! !REVISION HISTORY:
! 	20Feb97 - Jing Guo <guo@eramus> - defined template
!
!EOP
!_______________________________________________________________________

  character(len=*),parameter :: myname_='die.'
  integer :: myrank,ier
  character(len=16) :: lineno

	!-------------------------------------------------
	! MPI_ should have been initialized for this call
	!-------------------------------------------------

    call MP_comm_rank(MP_comm_world,myrank,ier)

	! a message for the users:

    write(lineno,'(i16)') line

    write(stderr,'(z3.3,9a)') myrank,'.',myname_,	&
      ': from ',trim(where),'()',	&
      ', line ',trim(adjustl(lineno)),	&
      ' of file ',fnam

	! raise a condition to the OS

    call MP_abort(MP_comm_world,2,ier)

end subroutine diex_
!=======================================================================
end module m_dropdead
!.
