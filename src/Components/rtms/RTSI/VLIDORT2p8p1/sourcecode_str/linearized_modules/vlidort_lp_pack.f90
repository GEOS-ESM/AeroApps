
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p1                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr                          #
! #                Matt Christi                                 #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p1                             #
! #  Release Date :   31 August 2019                            #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
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
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.1 of the VLIDORT software library.          #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2019.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of VLIDORT_2p8p1 ( Version 2.8.1 )            #
! #                                                                 #
! # VLIDORT_2p8p1 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p1 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p1  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_PACK_LAP_MISC                            #
! #            VLIDORT_PACK_LP_MULT                             #
! #                                                             #
! ###############################################################

      MODULE vlidort_lp_pack_m

      PRIVATE
      PUBLIC :: VLIDORT_PACK_LAP_MISC, &
                VLIDORT_PACK_LP_MULT

      CONTAINS

      SUBROUTINE VLIDORT_PACK_LAP_MISC ( &
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,                       & ! Input
        N_SZANGLES, N_USER_VZANGLES, NMOMENTS, N_TOTALPROFILE_WFS,      & ! Input
        L_DELTAU_VERT, L_OMEGA_GREEK,                                   & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,           & ! Input
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                 & ! Input
        LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                            & ! Input
        LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                               & ! Input
        LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS,                   & ! Input
        LAP_Misc )                                                        ! Output

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAXSTREAMS, &
                                 MAX_ATMOSWFS, MAX_PARTLAYERS, MAX_USER_STREAMS

      USE VLIDORT_LinWork_def_m

!  Input control
!  -------------

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_TOTALPROFILE_WFS

!  From VLIDORT_LAP_MISCSETUPS
!  ---------------------------

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL  ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      !DOUBLE PRECISION, INTENT (IN) :: L_TRUNC_FACTOR  ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized Transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Pseudo-spherical average-secant linearization
!mick mod 9/19/2017 - added LP_LEVELS_SOLARTRANS & LP_PARTIALS_SOLARTRANS to input

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output
!  ------

      TYPE(VLIDORT_LinWork_Miscellanous), INTENT (OUT) :: LAP_Misc

!  Pack structure

      !LAP_Misc%L_OMEGA_TOTAL ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
      !         L_OMEGA_TOTAL ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )
      !LAP_Misc%L_GREEKMAT_TOTAL ( 1:N_TOTALPROFILE_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ ) = &
      !         L_GREEKMAT_TOTAL ( 1:N_TOTALPROFILE_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ )
      !LAP_Misc%L_DELTAU_SLANT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )  = &
      !         L_DELTAU_SLANT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )
      !LAP_Misc%DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALPROFILE_WFS ) = &
      !         DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALPROFILE_WFS )
      !LAP_Misc%L_TRUNC_FACTOR ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
      !         L_TRUNC_FACTOR ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )

      LAP_Misc%L_DELTAU_VERT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
               L_DELTAU_VERT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )
      LAP_Misc%L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALPROFILE_WFS ) = &
               L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALPROFILE_WFS )

      LAP_Misc%L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALPROFILE_WFS ) = &
               L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS ) = &
               L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS ) = &
               L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS )

      LAP_Misc%L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )

      LAP_Misc%LP_AVERAGE_SECANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_AVERAGE_SECANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%LP_INITIAL_TRANS ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_INITIAL_TRANS ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%LP_T_DELT_MUBAR ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_T_DELT_MUBAR ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%LP_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%LP_LEVELS_SOLARTRANS ( 0:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_LEVELS_SOLARTRANS ( 0:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LAP_Misc%LP_PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
               LP_PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

      END SUBROUTINE VLIDORT_PACK_LAP_MISC

!

      SUBROUTINE VLIDORT_PACK_LP_MULT ( &
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, N_TOTALPROFILE_WFS, & ! Input
        LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,               & ! Input
        LP_Mult )                                                                 ! Output

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_ATMOSWFS

      USE VLIDORT_LinWork_def_m

!  Input control
!  -------------

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_TOTALPROFILE_WFS

!  From LP_EMULT_MASTER or LP_EMULT_MASTER_OBSGEO
!  ----------------------------------------------

      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output
!  ------

      TYPE(VLIDORT_LinWork_Multiplier), INTENT (OUT) :: LP_Mult

!  Pack structure

      LP_Mult%LP_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
              LP_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LP_Mult%LP_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
              LP_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

      LP_Mult%LP_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
              LP_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
      LP_Mult%LP_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
              LP_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

      END SUBROUTINE VLIDORT_PACK_LP_MULT

      END MODULE vlidort_lp_pack_m
