!
! This is intented as a OMS template for a generic instrument XXXX.
! Usually, you would replace XXXX with the name of your instrument,
! and create the OMS type specific for your instrument.
!

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_mxxxx --- Implements XXXX Observation Metadata Stream 
! 
! !INTERFACE:
!

    module  m_mxxxx

! !USES:

   implicit none
   
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  mxxxx_init
   PUBLIC  mxxxx_clean
   PUBLIC  mxxxx_get
   PUBLIC  mxxxx_put
!
! !PUBLIC DATA MEMBERS:
!
!
! !DESCRIPTION:
! \label{MXXX:Mod}
!  This module implements the XXXX Observation Metadata Stream,
!
! !REVISION HISTORY: 
!
!  15Nov2001  da Silva  First crack.
!  24Sep2002  C. Redder  Redefined protex label to prevent conflicts with
!                        labels in other modules
!
!EOP
!-------------------------------------------------------------------------

!  XXXX metadata type
!  ------------------
   type XXXX_OMS
        integer           :: nsta       ! No. of stations
        integer           :: nvct       ! Max size of vectors
   end type XXXX_OMS

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MXXXX_Init --- Allocates necessary memory
! 
! !INTERFACE:

  subroutine MXXXX_Init ( this, nsta, rc )

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS: 
!
    integer, intent(in)            :: nsta  ! number of stations

! !INPUT/OUTPUT PARAMETERS:
!
    type(xxxx_oms), intent(inout)  ::  this   ! XXXX OMS

    integer, intent(out)           ::  rc     ! Error return code:
                                              ! = 0    - all is well
!
! !DESCRIPTION: This routine allocates memory for hold XXXX OMS.
!
! !REVISION HISTORY: 
!
!  07Apr2002   da Silva   First crack.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mxxxx_init'

! All done
! --------
  return


end subroutine MXXXX_Init

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MXXXX_Clean --- Deallocates memory used by XXXX OMS
! 
! !INTERFACE:

  subroutine MXXXX_Clean ( this, rc )

!
! !USES
!
    Implicit NONE

! !INPUT/OUTPUT PARAMETERS:
!
    type(xxxx_oms), intent(inout)  ::  this   ! XXXX OMS

! !OUTPUT PARAMETERS:
!
    integer, intent(out)           ::  rc     ! Error return code:
                                              ! = 0    - all is well
!
! !DESCRIPTION: This routine allocates memory for hold XXXX OMS.
!
! !REVISION HISTORY: 
!
!  07Apr2002   da Silva   First crack.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mxxxx_clean'

! All done
! --------
  return

end subroutine MXXXX_Clean

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MXXXX_Get --- Reads XXXX OMS for a given synoptic time
! 
! !INTERFACE:

  subroutine MXXXX_Get ( fname, nymd, nhms, this, rc )

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! OMS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

!
! !OUTPUT PARAMETERS:
!
    type(xxxx_oms), intent(inout)  :: this    ! XXXX OMS

    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well

! !DESCRIPTION:
! \label{MXXX:Get}
!  This routine reads metadata from a XXXX OMS file, allocates the necessary
!  memory for the OMS vector, and loads the data for the synoptic time
!  (nymd,nhms). Usually the OMS file is opened, the data is read in and
!  the file is closed. 
!
! !REVISION HISTORY: 
!
!  07apr2002  da Silva   Initial code.
!  24Sep2002  C. Redder  Redefined protex label to prevent conflicts with
!                        labels in other modules
!  11Dec2002  Herdies    Made intent of this inout
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mxxxx_get'

  end subroutine mxxxx_get

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MXXXX_Put --- Writes XXXX OMS for a given synoptic time
! 
! !INTERFACE:

  subroutine MXXXX_Put ( fname, nymd, nhms, this, rc )

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! OMS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000
    type(xxxx_oms), intent(in)     :: this    ! XXXX OMS

!
! !OUTPUT PARAMETERS:
!
    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well

! !DESCRIPTION:
!  This routine writes metadata to a XXXX OMS file on a given synoptic
!  time. Usually the OMS file is opened, the metadata is written to
!  the file is closed. 
!
! !REVISION HISTORY: 
!
!  07apr2002  da Silva   Initial code.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mxxxx_put'

  end subroutine mxxxx_put

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MXXXX_Map --- Maps ODS xm into OMS index
! 
! !INTERFACE:

  subroutine MXXXX_Map ( nobs, xm, this, ioms, rc )

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS: 
!
    integer, intent(in)         :: nobs     ! size of XM vector
    real, intent(in)            :: xm(nobs) ! ODS metadata index
    type(xxxx_oms), intent(in)  :: this     ! XXXX OMS

!
! !OUTPUT PARAMETERS:
!
    integer, intent(out)        :: ioms(nobs) ! size of XM vector
    integer, intent(out)        :: rc         ! Error return code:
                                              !  0 - all is well

! !DESCRIPTION:
!  This routine translates the ODS metadata index XM into the OMS 
!  metadata index IOMS. With this index, one can address the OMS
!  vector directory, say, this%id(ioms(i)).
!
! !REVISION HISTORY: 
!
!  07apr2002  da Silva   Initial code.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mxxxx_map'

  end subroutine mxxxx_map

end module m_mxxxx
