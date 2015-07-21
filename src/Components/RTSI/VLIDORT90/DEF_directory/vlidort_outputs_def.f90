
! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

      MODULE VLIDORT_Outputs_def

!  This module contains the following structures:

!                 VLIDORT_Main_Outputs  nested in VLIDORT_Outputs
!           VLIDORT_Exception_Handling  nested in VLIDORT_Outputs
!     VLIDORT_Input_Exception_Handling  Intent(Out) from Input settings
!                      VLIDORT_Outputs  Intent(Out)

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Main_Outputs


!  Fourier values

!      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
!        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_STOKES_F

!  Fourier-summed values

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_STOKES

!  Results for mean-value output:

!  Complete Actinic and Regular Fluxes (including Direct terms)

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_MEAN_STOKES

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_STOKES

!  Direct Fluxes only

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES ) :: TS_MEAN_DIRECT

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES ) :: TS_FLUX_DIRECT

!  Contribution functions (TOA Upwelling only)

!  Fourier component of MS component of Diffuse Field

!      REAL(fpk), DIMENSION ( MAX_USER_VZANGLES, MAX_SZANGLES, &
!        MAXSTOKES, MAXLAYERS  ) :: TS_MS_CONTRIBS_F

!  Fourier-summed values of Whole Diffuse Field

!      REAL(fpk), DIMENSION ( MAX_GEOMETRIES, MAXSTOKES, &
!        MAXLAYERS ) :: TS_CONTRIBS

!  Fourier numbers used (bookkeeping)

      INTEGER, dimension ( MAX_SZANGLES )  :: TS_FOURIER_SAVED

!  Number of geometries (bookkeeping output)

      INTEGER                              :: TS_N_GEOMETRIES

!  Offsets for indexing geometries

      INTEGER :: TS_SZA_OFFSETS ( MAX_SZANGLES )
      INTEGER :: TS_VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      REAL(fpk), dimension ( MAXBEAMS ) :: TS_SOLARBEAM_BOATRANS

!  Contribution function (TOA Upwelling only)

!      REAL(fpk), dimension ( MAX_GEOMETRIES, MAXSTOKES, &
!        MAXLAYERS ) :: TS_SS_CONTRIBS


      END TYPE VLIDORT_Main_Outputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Exception_Handling


!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTCHECK
      INTEGER      :: TS_NCHECKMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_CHECKMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_ACTIONS

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER             :: TS_STATUS_CALCULATION
      CHARACTER (Len=120) :: TS_MESSAGE, TS_TRACE_1, TS_TRACE_2, TS_TRACE_3


      END TYPE VLIDORT_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Input_Exception_Handling


!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTREAD
      INTEGER      :: TS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTACTIONS


      END TYPE VLIDORT_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Outputs


      TYPE(VLIDORT_Main_Outputs)       :: Main
      TYPE(VLIDORT_Exception_Handling) :: Status


      END TYPE VLIDORT_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Main_Outputs, &
                VLIDORT_Exception_Handling, &
                VLIDORT_Input_Exception_Handling, &
                VLIDORT_Outputs

      END MODULE VLIDORT_Outputs_def
