
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

      MODULE VLIDORT_LinWork_def

!  This module contains the following VLIDORT linearized work structures:
!  (n/a = not available at present)

!         VLIDORT_LinWork_RadTran    n/a
!            VLIDORT_LinWork_User    n/a
!           VLIDORT_LinWork_Solar    n/a
!         VLIDORT_LinWork_Thermal
!        VLIDORT_LinWork_Geometry    n/a
!    VLIDORT_LinWork_Miscellanous
!      VLIDORT_LinWork_Multiplier
!     VLIDORT_LinWork_Corrections    n/a
!      VLIDORT_LinWork_FirstOrder    n/a
!         VLIDORT_LinWork_Control    n/a
!         VLIDORT_LinWork_Optical    n/a
!            VLIDORT_LinWork_BRDF    n/a
!          VLIDORT_LinWork_SLEAVE    n/a
!          VLIDORT_LinWork_Output    n/a

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinWork_Thermal


      DOUBLE PRECISION :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )


      END TYPE VLIDORT_LinWork_Thermal

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinWork_Miscellanous


      DOUBLE PRECISION :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_GREEKMAT_TOTAL &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: L_DELTAU_SLANT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      LOGICAL ::          DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LP_AVERAGE_SECANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

      DOUBLE PRECISION :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)

      END TYPE VLIDORT_LinWork_Miscellanous

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinWork_Multiplier


      DOUBLE PRECISION :: LP_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

      END TYPE VLIDORT_LinWork_Multiplier

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: & !VLIDORT_LinWork_RadTran, &
                !VLIDORT_LinWork_User, &
                !VLIDORT_LinWork_Solar, &    
                VLIDORT_LinWork_Thermal, &    
                !VLIDORT_LinWork_Geometry, &    
                VLIDORT_LinWork_Miscellanous, &    
                VLIDORT_LinWork_Multiplier!, &    
                !VLIDORT_LinWork_Corr_Nadir, &
                !VLIDORT_LinWork_Corr_OutGo, &
                !VLIDORT_LinWork_Corr_DB, &
                !VLIDORT_LinWork_FirstOrder, &    
                !VLIDORT_LinWork_Control, &    
                !VLIDORT_LinWork_Optical, &    
                !VLIDORT_LinWork_BRDF, &    
                !VLIDORT_LinWork_SLEAVE, &    
                !VLIDORT_LinWork_Output, & 
                !VLIDORT_LinWork_Corrections

   END MODULE VLIDORT_LinWork_def
