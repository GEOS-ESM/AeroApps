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
! #            VLIDORT_PACK_LAC_MISC                            #
! #            VLIDORT_PACK_LC_MULT                             #
! #                                                             #
! ###############################################################

      MODULE vlidort_lc_pack

      PRIVATE
      PUBLIC :: VLIDORT_PACK_LAC_MISC, &
                VLIDORT_PACK_LC_MULT

      CONTAINS

      SUBROUTINE VLIDORT_PACK_LAC_MISC ( &
        NSTOKES, NSTOKES_SQ, NSTREAMS, NLAYERS, N_PARTLAYERS,           & ! Input
        N_SZANGLES, N_USER_VZANGLES, NMOMENTS, N_TOTALCOLUMN_WFS,       & ! Input
        L_DELTAU_VERT, L_OMEGA_GREEK,                                   & ! Input
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                            & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,           & ! Input
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,                               & ! Input
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                 & ! Input
        LAC_Misc )                                                        ! Output

      USE VLIDORT_PARS
      USE VLIDORT_LinWork_def

!  Input
!  -----

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS

!  From VLIDORT_LAC_MISCSETUPS

      !DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL &
      !    ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_SLANT &
      !    ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )

      !LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION &
      !    ( MAXLAYERS, MAX_ATMOSWFS )
      !DOUBLE PRECISION, INTENT (IN) :: L_TRUNC_FACTOR &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)

!  Output
!  ------

      TYPE(VLIDORT_LinWork_Miscellanous), INTENT (OUT) :: LAC_Misc

!  Pack structure

      !LAC_Misc%L_OMEGA_TOTAL ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS ) = &
      !         L_OMEGA_TOTAL ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS )
      LAC_Misc%L_DELTAU_VERT ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS ) = &
               L_DELTAU_VERT ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS )
      !LAC_Misc%L_GREEKMAT_TOTAL ( 1:N_TOTALCOLUMN_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ ) = &
      !         L_GREEKMAT_TOTAL ( 1:N_TOTALCOLUMN_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ )
      !LAC_Misc%L_DELTAU_SLANT ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES ) = &
      !         L_DELTAU_SLANT ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )

      !LAC_Misc%DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALCOLUMN_WFS ) = &
      !         DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALCOLUMN_WFS )
      !LAC_Misc%L_TRUNC_FACTOR ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS ) = &
      !         L_TRUNC_FACTOR ( 1:N_TOTALCOLUMN_WFS, 1:NLAYERS )
      LAC_Misc%L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALCOLUMN_WFS ) = &
               L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALCOLUMN_WFS )

      LAC_Misc%L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALCOLUMN_WFS )

      LAC_Misc%L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALCOLUMN_WFS )

      LAC_Misc%LC_AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               LC_AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%LC_INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               LC_INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%LC_T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               LC_T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )
      LAC_Misc%LC_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
               LC_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )

      END SUBROUTINE VLIDORT_PACK_LAC_MISC

!

      SUBROUTINE VLIDORT_PACK_LC_MULT ( &
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,       & ! Input
        N_TOTALCOLUMN_WFS,                                        & ! Input
        LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN, & ! Input
        LC_Mult )                                                   ! Output

      USE VLIDORT_PARS
      USE VLIDORT_LinWork_def

!  Input
!  -----

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS

!  From LC_EMULT_MASTER or LC_EMULT_MASTER_OBSGEO

      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output
!  ------

      TYPE(VLIDORT_LinWork_Multiplier), INTENT (OUT) :: LC_Mult

!  Pack structure

      LC_Mult%LC_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
              LC_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )
      LC_Mult%LC_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
              LC_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )

      LC_Mult%LC_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
              LC_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )
      LC_Mult%LC_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS ) = &
              LC_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES, 1:N_TOTALCOLUMN_WFS )

      END SUBROUTINE VLIDORT_PACK_LC_MULT

      END MODULE vlidort_lc_pack
