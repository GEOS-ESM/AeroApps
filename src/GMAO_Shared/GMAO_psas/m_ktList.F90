!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ktList - a list of externally defined KT values
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ktList
      implicit none
      private	! except

      public :: ktUs
      public :: ktVs
      public :: ktslp
      public :: ktUU
      public :: ktVV
      public :: ktHH
      public :: ktQQ

      public :: KTMAX

      public :: ktList_dim	! check the class of variable type

      interface ktList_dim; module procedure dim_; end interface

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- converted from "ktmax.h"
!	24Mar95 - Jing G.
!		- New version for text based data table.
!	02/05/93 - ktmax.h
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ktList'

  integer,parameter :: ktUs  = 1	! sea-level wind
  integer,parameter :: ktVs  = 2	! sea-level wind
  integer,parameter :: ktslp = 3	! sea-level pressure
  integer,parameter :: ktUU  = 4	! upper-air wind
  integer,parameter :: ktVV  = 5	! upper-air wind
  integer,parameter :: ktHH  = 6	! upper-air hight
  integer,parameter :: ktQQ  = 7	! upper-air mixing ratio

  integer,parameter :: KTMAX = 7	! Total number of data types

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dim_ - find if the variable has a "third" demension
!
! !DESCRIPTION:
!
! !INTERFACE:

    function dim_(kt)
      implicit none
      integer,intent(in) :: kt
      integer :: dim_

! !REVISION HISTORY:
! 	08Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dim_'

  select case(kt)
  case(ktUs,ktVs,ktslp)
    dim_=2
  case(ktUU,ktVV,ktHH,ktQQ)
    dim_=3
  case default
    dim_=-1		! an unexpected variable type
  end select

end function dim_

end module m_ktList
