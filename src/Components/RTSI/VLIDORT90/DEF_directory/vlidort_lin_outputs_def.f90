
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

      MODULE VLIDORT_LinOutputs_def

!  This module contains the following VLIDORT output structures,
!  with intents :

!          VLIDORT_LinCol    nested in VLIDORT_LinOutputs
!         VLIDORT_LinProf    nested in VLIDORT_LinOutputs
!        VLIDORT_LinAtmos    nested in VLIDORT_LinOutputs
!         VLIDORT_LinSurf    nested in VLIDORT_LinOutputs
!      VLIDORT_LinOutputs    Intent(Out)

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinCol


!  Atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_COLUMNWF

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_MINT_COLUMNWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_COLUMNWF

!  Direct Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_MINT_COLUMNWF_DIRECT

!  Direct Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_FLUX_COLUMNWF_DIRECT


      END TYPE VLIDORT_LinCol

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinProf


!  Atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_PROFILEWF

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_MINT_PROFILEWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_PROFILEWF

!  Direct Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_MINT_PROFILEWF_DIRECT

!  Direct Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_FLUX_PROFILEWF_DIRECT


      END TYPE VLIDORT_LinProf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinAtmos

!  Atmospheric thermal emission: Blackbody weighting functions (BBWFS)
!     26 March 2014, Introduced for Version 2.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
        0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)           :: TS_ABBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, &
        0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)           :: TS_ABBWFS_FLUXES

!  These Jacobians were superseded in Version 2.7 by the above ABBWFS arrays
!      LTE linearization - Temperature weighting functions for BB functions (RT Solutions Use Only)
!      REAL(fpk), dimension ( 0:MAXLAYERS, MAX_USER_LEVELS, &
!        MAX_USER_VZANGLES, MAX_DIRECTIONS ) :: TS_LTE_ATMOSWF

      END TYPE VLIDORT_LinAtmos

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSurf

!  Surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_SURFACEWF

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_MINT_SURFACEWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_SURFACEWF

!  Surface thermal emission: Blackbody weighting functions (BBWFS)
!     26 March 2014, Introduced for Version 2.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                  MAXSTOKES, MAX_DIRECTIONS )  :: TS_SBBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, &
                  MAXSTOKES, MAX_DIRECTIONS )  :: TS_SBBWFS_FLUXES

      END TYPE VLIDORT_LinSurf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinOutputs

      TYPE(VLIDORT_LinCol)   :: Col
      TYPE(VLIDORT_LinProf)  :: Prof
      TYPE(VLIDORT_LinAtmos) :: Atmos
      TYPE(VLIDORT_LinSurf)  :: Surf

      END TYPE VLIDORT_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_LinCol,&
                VLIDORT_LinProf,&
                VLIDORT_LinAtmos,&
                VLIDORT_LinSurf,&
                VLIDORT_LinOutputs

      END MODULE VLIDORT_LinOutputs_def
