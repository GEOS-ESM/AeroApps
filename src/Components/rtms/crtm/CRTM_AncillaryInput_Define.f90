!
! CRTM_AncillaryInput_Define
!
! Module defining a container for CRTM ancillary input.
!
!
! CREATION HISTORY:
!       Modifed by:     Paul van Delst, 08-Dec-2009
!                       paul.vandelst@noaa.gov


MODULE CRTM_AncillaryInput_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE SSU_Input_Define
  USE Zeeman_Input_Define
  ! Disable implicit typing
  IMPLICIT NONE
  
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PRIVATE, PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_AncillaryInput_Define.f90,v 1.2 2010/06/09 10:10:22 rtodling Exp $'
  
  !--------------------
  ! Structure defintion
  !--------------------
  !:tdoc+:
  TYPE :: CRTM_AncillaryInput_type
    TYPE(SSU_Input_type)    :: SSU
    TYPE(Zeeman_Input_type) :: Zeeman
  END TYPE CRTM_AncillaryInput_type
  !:tdoc-:

END MODULE CRTM_AncillaryInput_Define
