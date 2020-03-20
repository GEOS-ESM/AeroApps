
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

      MODULE VLIDORT_LinSup_SS_def

!  This module contains the following structures,
!  with intents :

!       VLIDORT_LinSup_SS_Col    nested in VLIDORT_LinSup_SS_InOut
!      VLIDORT_LinSup_SS_Prof    nested in VLIDORT_LinSup_SS_InOut
!      VLIDORT_LinSup_SS_Surf    nested in VLIDORT_LinSup_SS_InOut
!     VLIDORT_LinSup_SS_InOut    Intent(InOut)

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Col


!  SS atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_COLUMNWF_SS

!  DB atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES ) :: TS_COLUMNWF_DB


      END TYPE VLIDORT_LinSup_SS_Col

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Prof


!  SS atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_PROFILEWF_SS

!  DB atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES ) :: TS_PROFILEWF_DB


      END TYPE VLIDORT_LinSup_SS_Prof

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Surf


!  SS surface weighting functions

!      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
!        MAXSTOKES, MAX_DIRECTIONS ) :: TS_SURFACEWF_SS

!  DB surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES ) :: TS_SURFACEWF_DB


      END TYPE VLIDORT_LinSup_SS_Surf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS


      TYPE(VLIDORT_LinSup_SS_Col)   :: Col
      TYPE(VLIDORT_LinSup_SS_Prof)  :: Prof
      TYPE(VLIDORT_LinSup_SS_Surf)  :: Surf


      END TYPE VLIDORT_LinSup_SS

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_LinSup_SS_Col,&
                VLIDORT_LinSup_SS_Prof,&
                VLIDORT_LinSup_SS_Surf,&
                VLIDORT_LinSup_SS

      END MODULE VLIDORT_LinSup_SS_def
