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
! #            VLIDORT_UNPACK_L_THERM                           #
! #                                                             #
! ###############################################################

      MODULE vlidort_l_unpack

      PRIVATE
      PUBLIC :: VLIDORT_UNPACK_L_THERM

      CONTAINS

      SUBROUTINE VLIDORT_UNPACK_L_THERM ( L_Therm,                                  & ! Input
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, N_TOTALATMOS_WFS, & ! Input
        L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                       & ! Output
        L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )            ! Output

      USE VLIDORT_PARS
      USE VLIDORT_LinWork_def

!  Input
!  -----

      TYPE(VLIDORT_LinWork_Thermal), INTENT (IN) :: L_Therm

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_TOTALATMOS_WFS

!  Output
!  ------

!  From THERMAL_SETUP_PLUS

      DOUBLE PRECISION, INTENT (OUT) :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Unpack structure

              L_THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )

              L_T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS )
              L_T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS )

              L_T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS )
              L_T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS )

      END SUBROUTINE VLIDORT_UNPACK_L_THERM

      END MODULE vlidort_l_unpack
