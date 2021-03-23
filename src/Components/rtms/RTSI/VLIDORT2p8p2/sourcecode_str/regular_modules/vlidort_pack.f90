
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p2                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_PACK_MISC                                #
! #            VLIDORT_PACK_THERM                               #
! #            VLIDORT_PACK_MULT                                #
! #                                                             #
! ###############################################################

!  4/15/20. Version 2.8.2. No changes

      MODULE vlidort_pack_m

      PUBLIC :: VLIDORT_PACK_MISC,  &
                VLIDORT_PACK_THERM, &
                VLIDORT_PACK_MULT

      CONTAINS

      SUBROUTINE VLIDORT_PACK_MISC ( &
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_VZANGLES, N_SZANGLES, NMOMENTS, & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                                          & ! Input Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,                             & ! Input Optical
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,                          & ! Input Solar
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                                  & ! Input Trans. D.O.
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,            & ! Input Trans solar/User
        CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                             & ! Input Av. secant.
        Misc )                                                                             ! Output structure to be packed

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, MAXSTOKES, MAXBEAMS
      USE VLIDORT_Work_def_m

!  Input
!  -----

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS

!  From VLIDORT_MISCSETUPS, optical variables
!mick mod 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to input

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT  ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )

!  Solar

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Solar beam Transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  User-angle transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Average-secant parameterization

      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA   ( 0:MAXLAYERS, MAXBEAMS )

!  Output
!  ------

      TYPE(VLIDORT_Work_Miscellanous), INTENT (OUT) :: Misc

!  Pack structure. Not including FAC1 and TRUNC_FACTOR

      Misc%DELTAU_VERT ( 1:NLAYERS )      = &
           DELTAU_VERT ( 1:NLAYERS )
      Misc%PARTAU_VERT ( 1:N_PARTLAYERS ) = &
           PARTAU_VERT ( 1:N_PARTLAYERS )
      Misc%DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES ) = &
           DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )

      Misc%LEVELS_SOLARTRANS ( 0:NLAYERS, 1:N_SZANGLES )        = &
           LEVELS_SOLARTRANS ( 0:NLAYERS, 1:N_SZANGLES )
      Misc%PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
           PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:N_SZANGLES )

      Misc%OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES ) = &
           OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES )

      !Misc%TRUNC_FACTOR ( 1:NLAYERS ) = TRUNC_FACTOR ( 1:NLAYERS )
      !Misc%FAC1 ( 1:NLAYERS )         = FAC1 ( 1:NLAYERS )

      Misc%BEAM_CUTOFF ( 1:N_SZANGLES ) = &
           BEAM_CUTOFF ( 1:N_SZANGLES )
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

      Misc%AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES ) = &
           AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES )
      Misc%LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES ) = &
           LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES )

      !Misc%ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      !     ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )

      END SUBROUTINE VLIDORT_PACK_MISC

!

      SUBROUTINE VLIDORT_PACK_THERM ( &
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input Numbers
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input Thermal setup
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input Thermal solutions
        Therm )                                                     ! Output

      USE VLIDORT_PARS_m, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS
      USE VLIDORT_Work_def_m

!  Input
!  -----

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  From THERMAL_SETUP

      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Thermal direct solutions, also from THERMAL_SETUP

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

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
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input numbers
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input multipliers
        Mult )                                                ! Output

      USE VLIDORT_PARS_m, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS
      USE VLIDORT_Work_def_m

!  Input
!  -----

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  Multipliers from EMULT_MASTER or EMULT_MASTER_OBSGEO

      !LOGICAL, INTENT (IN) ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

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

      END MODULE vlidort_pack_m

