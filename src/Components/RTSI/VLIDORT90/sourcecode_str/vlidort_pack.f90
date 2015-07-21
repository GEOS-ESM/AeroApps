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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_PACK_MISC                                #
! #            VLIDORT_PACK_THERM                               #
! #            VLIDORT_PACK_MULT                                #
! #                                                             #
! ###############################################################

      MODULE vlidort_pack

      PRIVATE
      PUBLIC :: VLIDORT_PACK_MISC, &
                VLIDORT_PACK_THERM, &
                VLIDORT_PACK_MULT

      CONTAINS

      SUBROUTINE VLIDORT_PACK_MISC ( &
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_VZANGLES, N_SZANGLES,& ! Input 
        NMOMENTS,                                                             & ! Input
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, OMEGA_GREEK,                  & ! Input
        LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,          & ! Input
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                       & ! Input
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, & ! Input
        CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                  & ! Input
        Misc )                                                                  ! Output

      USE VLIDORT_PARS
      USE VLIDORT_Work_def

!  Input
!  -----

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS

!  From VLIDORT_MISCSETUPS

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

      !DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT (IN) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Output
!  ------

      TYPE(VLIDORT_Work_Miscellanous), INTENT (OUT) :: Misc

!  Pack structure

      Misc%DELTAU_VERT ( 1:NLAYERS ) = &
           DELTAU_VERT ( 1:NLAYERS )
      Misc%PARTAU_VERT ( 1:N_PARTLAYERS ) = &
           PARTAU_VERT ( 1:N_PARTLAYERS )
      Misc%DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES ) = &
           DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )

      !Misc%TRUNC_FACTOR ( 1:NLAYERS ) = &
      !     TRUNC_FACTOR ( 1:NLAYERS )
      !Misc%FAC1 ( 1:NLAYERS ) = &
      !     FAC1 ( 1:NLAYERS )

      Misc%OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES ) = &
           OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES )
      Misc%LAYER_PIS_CUTOFF ( 1:N_SZANGLES ) = &
           LAYER_PIS_CUTOFF ( 1:N_SZANGLES )
      Misc%TRANS_SOLAR_BEAM ( 1:N_SZANGLES ) = &
           TRANS_SOLAR_BEAM ( 1:N_SZANGLES )
      Misc%DO_REFLECTED_DIRECTBEAM ( 1:N_SZANGLES ) = &
           DO_REFLECTED_DIRECTBEAM ( 1:N_SZANGLES )

      Misc%T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS ) = &
           T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS )
      Misc%T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS )  = &
           T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS )
      Misc%T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS )  = &
           T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS )

      Misc%T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES ) = &
           T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES )
      Misc%T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
           T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES )
      !Misc%T_UTUP_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
      !     T_UTUP_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES )

      Misc%T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES ) = &
           T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES )
      Misc%T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES ) = &
           T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES )
      Misc%T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES ) = &
           T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES )

      Misc%CUMTRANS ( 1:NLAYERS, 1:N_USER_VZANGLES ) = &
           CUMTRANS ( 1:NLAYERS, 1:N_USER_VZANGLES )
      Misc%INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES ) = &
           INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES )
      !Misc%ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      !     ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )
      Misc%AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES ) = &
           AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES )
      Misc%LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES ) = &
           LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES )

      END SUBROUTINE VLIDORT_PACK_MISC

!

      SUBROUTINE VLIDORT_PACK_THERM ( &
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input
        Therm )                                                     ! Output

      USE VLIDORT_PARS
      USE VLIDORT_Work_def

!  Input
!  -----

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  From THERMAL_SETUP

      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Output
!  ------

      TYPE(VLIDORT_Work_Thermal), INTENT (OUT) :: Therm

!  Pack structure

      Therm%THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
            THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS )
      Therm%DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
            DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS )
      Therm%XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS ) = &
            XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS )
      Therm%TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
            TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS )

      Therm%T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS ) = &
            T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS )
      Therm%T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS ) = &
            T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS )

      Therm%T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS ) = &
            T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS )
      Therm%T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS ) = &
            T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS )

      END SUBROUTINE VLIDORT_PACK_THERM

!

      SUBROUTINE VLIDORT_PACK_MULT ( &
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input
        Mult )                                                ! Output

      USE VLIDORT_PARS
      USE VLIDORT_Work_def

!  Input
!  -----

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  From EMULT_MASTER or EMULT_MASTER_OBSGEO

      !LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Output
!  ------

      TYPE(VLIDORT_Work_Multiplier), INTENT (OUT) :: Mult

!  Pack structure

      !Mult%EMULT_HOPRULE ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      !     EMULT_HOPRULE ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )

      !Mult%SIGMA_M ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      !     SIGMA_M ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )
      !Mult%SIGMA_P ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      !     SIGMA_P ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )

      Mult%EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES ) = &
           EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES )
      Mult%EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES ) = &
           EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES )

      Mult%UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )  = &
           UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )
      Mult%UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )  = &
           UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )

      END SUBROUTINE VLIDORT_PACK_MULT

      END MODULE vlidort_pack
