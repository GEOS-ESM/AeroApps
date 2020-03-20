
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

      MODULE VLIDORT_Work_def

!  This module contains the following VLIDORT work structures:
!  (n/a = not available at present)

!         VLIDORT_Work_RadTran    n/a
!            VLIDORT_Work_User    n/a
!           VLIDORT_Work_Solar    n/a
!         VLIDORT_Work_Thermal
!        VLIDORT_Work_Geometry    n/a
!    VLIDORT_Work_Miscellanous
!      VLIDORT_Work_Multiplier
!     VLIDORT_Work_Corrections
!      VLIDORT_Work_FirstOrder
!         VLIDORT_Work_Control    n/a
!         VLIDORT_Work_Optical    n/a
!            VLIDORT_Work_BRDF
!          VLIDORT_Work_SLEAVE
!          VLIDORT_Work_Output    n/a

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Thermal


      DOUBLE PRECISION :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )


      END TYPE VLIDORT_Work_Thermal

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Miscellanous


      DOUBLE PRECISION :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      INTEGER ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )


      END TYPE VLIDORT_Work_Miscellanous

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Multiplier


      LOGICAL ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )


      END TYPE VLIDORT_Work_Multiplier

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Corr_Nadir


      DOUBLE PRECISION :: ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: TMS ( MAXLAYERS )
      DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )
      DOUBLE PRECISION :: SS_CUMSOURCE_UP &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: SS_CUMSOURCE_DN &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )


      END TYPE VLIDORT_Work_Corr_Nadir

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Corr_OutGo


      DOUBLE PRECISION :: UP_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: UP_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: BOA_ATTN ( MAX_GEOMETRIES )


      END TYPE VLIDORT_Work_Corr_OutGo

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Corr_DB


      DOUBLE PRECISION :: DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: EXACTDB_SOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: ATTN_DB_SAVE ( MAX_GEOMETRIES )


      END TYPE VLIDORT_Work_Corr_DB

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_FirstOrder


      DOUBLE PRECISION :: FO_STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES_DTA &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DTS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )


      END TYPE VLIDORT_Work_FirstOrder

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_BRDF


      DOUBLE PRECISION ::   EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )


      END TYPE VLIDORT_Work_BRDF

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_SLEAVE


      DOUBLE PRECISION ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   SLTERM_USERANGLES &
          ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )


      END TYPE VLIDORT_Work_SLEAVE 

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Corrections


      TYPE(VLIDORT_Work_Corr_Nadir) :: Nadir
      TYPE(VLIDORT_Work_Corr_OutGo) :: OutGo
      TYPE(VLIDORT_Work_Corr_DB)    :: DB


      END TYPE VLIDORT_Work_Corrections

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: & !VLIDORT_Work_RadTran, &
                !VLIDORT_Work_User, &
                !VLIDORT_Work_Solar, &    
                VLIDORT_Work_Thermal, &    
                !VLIDORT_Work_Geometry, &    
                VLIDORT_Work_Miscellanous, &    
                VLIDORT_Work_Multiplier, &    
                VLIDORT_Work_Corr_Nadir, &
                VLIDORT_Work_Corr_OutGo, &
                VLIDORT_Work_Corr_DB, &
                VLIDORT_Work_FirstOrder, &    
                !VLIDORT_Work_Control, &    
                !VLIDORT_Work_Optical, &    
                VLIDORT_Work_BRDF, &    
                VLIDORT_Work_SLEAVE, &    
                !VLIDORT_Work_Output, & 
                VLIDORT_Work_Corrections

      END MODULE VLIDORT_Work_def
