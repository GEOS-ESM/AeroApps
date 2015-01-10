!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_gdstat - interface
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_gdstat
      implicit none
      private	! except

      public :: gdstat		! The class data structure
      public :: lvstat		! The class data structure

      interface lvstat
	subroutine LVSTAT (lu,mx,my,a,h,atype,htype,amiss,flag)
	implicit none
	integer,intent(in) :: lu		! Output unit
	integer,intent(in) :: mx,my		! Array sizes
	real,intent(in) :: a(mx,my)		! The array
	real,intent(in) :: h			! The argument
	character*4,intent(in) :: atype	! Type of the variable(array)
	character*4,intent(in) :: htype	! Type of the level(argument)
	real,intent(in) :: amiss	! missing value flag of a
	character*(*),intent(in) :: flag
	end subroutine LVSTAT
      end interface

      interface gdstat
	subroutine GDSTAT (lu,mx,my,mz,a,h,atype,htype,amiss,header,inc)
	implicit none
	integer,intent(in) :: lu		! Output unit
	integer,intent(in) :: mx,my,mz	! Array sizes
	real,intent(in) :: a(mx,my,mz)	! The array
	real,intent(in) :: h(mz)		! The argument(levels)
	character*4,intent(in) :: atype	! Type of the variable
	character*4,intent(in) :: htype	! Typf of the levels
	real,intent(in) :: amiss	! missing value flag of a
	character*(*),intent(in) :: header	! A header message
	integer,intent(in) :: inc		! order of the listing
	end subroutine GDSTAT
      end interface

! !REVISION HISTORY:
! 	05Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_gdstat'

end module m_gdstat
