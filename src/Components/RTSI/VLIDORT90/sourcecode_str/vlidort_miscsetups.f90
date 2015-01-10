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
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
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
! #            VLIDORT_MISCSETUPS (master, calling:)            #
! #              VLIDORT_DELTAMSCALE                            #
! #              VLIDORT_SSALBINIT                              #
! #              VLIDORT_QSPREP                                 #
! #              VLIDORT_PREPTRANS                              #
! #                                                             #
! #            VLIDORT_SURFACE_DIRECTBEAM                       #
! #                                                             #
! #            VLIDORT_CHAPMAN                                  #
! #              BEAM_GEOMETRY_PREPARE                          #
! #                                                             #
! ###############################################################


      MODULE vlidort_miscsetups_module

      PRIVATE
      PUBLIC :: VLIDORT_MISCSETUPS, &
                VLIDORT_SURFACE_DIRECTBEAM, &
                VLIDORT_CHAPMAN

      CONTAINS

      SUBROUTINE VLIDORT_MISCSETUPS ( &
        DO_DELTAM_SCALING, NSTOKES, &
        NLAYERS, N_USER_LEVELS, &
        OMEGA_TOTAL_INPUT, &
        DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        NMOMENTS, NBEAMS, &
        PARTLAYERS_OUTFLAG, DO_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        TAUGRID_INPUT, CHAPMAN_FACTORS, &
        MUELLER_INDEX, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        DO_SPECIALIST_OPTION_3, NLAYERS_CUTOFF, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        DO_SOLUTION_SAVING, NSTREAMS, &
        DO_TOA_CONTRIBS, QUAD_STREAMS, &
        N_USER_STREAMS, DO_USER_STREAMS, &
        USER_SECANTS, N_PARTLAYERS, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        OMEGA_TOTAL, DELTAU_VERT, &
        PARTAU_VERT, GREEKMAT_TOTAL, &
        TAUGRID, DELTAU_SLANT, &
        TRUNC_FACTOR, FAC1, &
        OMEGA_GREEK, LAYER_PIS_CUTOFF, &
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
        INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_DISORDS_UTDN, T_DELT_MUBAR, &
        T_UTDN_MUBAR, T_UTUP_MUBAR, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, CUMTRANS, ITRANS_USERM )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NBEAMS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_3
      INTEGER, INTENT (IN) ::          NLAYERS_CUTOFF
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (OUT) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )
      LOGICAL, INTENT (OUT) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_MUBAR &
          ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_MUBAR &
          ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  miscellaneous setup operations, Master routine

!  Performance set-up is in VLIDORT_DERIVE_INPUTS (vlidort_inputs.f)
!      CALL VLIDORT_PERFORMANCE_SETUP

!  Delta-m scaling of input quantities

      CALL VLIDORT_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, &
        NLAYERS, N_USER_LEVELS, &
        OMEGA_TOTAL_INPUT, &
        DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        NMOMENTS, NBEAMS, &
        PARTLAYERS_OUTFLAG, DO_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        TAUGRID_INPUT, CHAPMAN_FACTORS, &
        OMEGA_TOTAL, DELTAU_VERT, &
        PARTAU_VERT, GREEKMAT_TOTAL, &
        TAUGRID, DELTAU_SLANT, &
        TRUNC_FACTOR, FAC1 )

!  initialise single scatter albedo terms
!    GREEKMAT x OMEGA

      CALL VLIDORT_SSALBINIT ( &
        NSTOKES, NLAYERS, &
        NMOMENTS, MUELLER_INDEX, &
        OMEGA_TOTAL, GREEKMAT_TOTAL, &
        OMEGA_GREEK )

!  Prepare quasi-spherical attenuation

      CALL VLIDORT_QSPREP ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        NLAYERS, &
        DO_SPECIALIST_OPTION_3, NLAYERS_CUTOFF, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        NBEAMS, DELTAU_VERT, &
        TAUGRID, DELTAU_SLANT, &
        LAYER_PIS_CUTOFF, &
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        LOCAL_CSZA )

!  Transmittances and Transmittance factors

      CALL VLIDORT_PREPTRANS ( &
        DO_SOLUTION_SAVING, NSTREAMS, &
        NLAYERS, DO_TOA_CONTRIBS, &
        LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, &
        USER_SECANTS, N_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        DELTAU_VERT, PARTAU_VERT, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_DISORDS_UTDN, T_DELT_MUBAR, &
        T_UTDN_MUBAR, T_UTUP_MUBAR, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, CUMTRANS, &
        ITRANS_USERM )

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_MISCSETUPS

!

      SUBROUTINE VLIDORT_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, &
        NLAYERS, N_USER_LEVELS, &
        OMEGA_TOTAL_INPUT, &
        DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        NMOMENTS, NBEAMS, &
        PARTLAYERS_OUTFLAG, DO_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        TAUGRID_INPUT, CHAPMAN_FACTORS, &
        OMEGA_TOTAL, DELTAU_VERT, &
        PARTAU_VERT, GREEKMAT_TOTAL, &
        TAUGRID, DELTAU_SLANT, &
        TRUNC_FACTOR, FAC1 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NBEAMS
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTLAYERS_VALUES   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )

!  local variables

      DOUBLE PRECISION :: FDEL, FAC2, DNL1, FDNL1
      DOUBLE PRECISION :: DNM1, DELS, DT, XTD
      INTEGER          :: K, K1, K2, N, N1, L, UT, UTA, NM1, IB

!  Indexing

      INTEGER, DIMENSION(4) :: KTYPE1, KTYPE2

!mick - singularity buster output
      INTEGER   :: I
      LOGICAL   :: SBUST(6)

!  Setup the two ktype arrays

      KTYPE1 = (/ 1, 6, 11, 16 /)
      KTYPE2 = (/ 2, 5, 12, 15 /)

!  First do the slant input
!  ------------------------

!  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!  debug

!      do n = 1, nlayers
!       dels = height_grid(n-1)-height_grid(n)
!       write(*,*)n,dels,height_grid(n)
!       write(44,'(i4,101f10.5)')n,
!     &  (CHAPMAN_FACTORS(N,K,1)*(height_grid(k-1)-height_grid(k))
!     &           ,k=1,nlayers)
!      enddo
!      pause

!mick fix - initialise
      GREEKMAT_TOTAL = ZERO

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

!  New section by R. Spurr, RT Solutions Inc.
!   Based in part on the code in VDISORT.

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DBLE(2*NM1+1)

!  Scaling for layer input
!  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

!  overall truncation factor

          FDEL = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

!  Scale Greek Matrix entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = &
                  ( GREEKMAT_TOTAL_INPUT(L,N,K1) - FDNL1 ) / FAC2
                GREEKMAT_TOTAL(L,N,K2) = &
                    GREEKMAT_TOTAL_INPUT(L,N,K2) / FAC2
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = &
                  ( GREEKMAT_TOTAL_INPUT(L,N,1) - FDNL1 ) / FAC2
            ENDIF
          ENDDO

!  Maintain phase function normalization

          GREEKMAT_TOTAL(0,N,1) = ONE

!  scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

!  end layer loop

        ENDDO

!  Scaling for user-defined off-grid optical depths
!  ------------------------------------------------

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled)

        IF ( DO_PARTLAYERS ) THEN
          UT = 0
          DO UTA = 1, N_USER_LEVELS
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT  = UT + 1
              N   = PARTLAYERS_LAYERIDX(UT)
              DT  = PARTLAYERS_VALUES(UT)
              XTD = DELTAU_VERT_INPUT(N) * DT
              PARTAU_VERT(UT) = XTD * FAC1(N)
            ENDIF
          ENDDO
        ENDIF

!  Scale layer path thickness values

        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          FAC1(N)        = ONE
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
          TAUGRID(N)     = TAUGRID_INPUT(N)
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N)
          DO L = 0, NMOMENTS
            DO K = 1, 4
              K1 = KTYPE1(K)
              K2 = KTYPE2(K)
              GREEKMAT_TOTAL(L,N,K1) = &
                    GREEKMAT_TOTAL_INPUT(L,N,K1)
              GREEKMAT_TOTAL(L,N,K2) = &
                    GREEKMAT_TOTAL_INPUT(L,N,K2)
            ENDDO
          ENDDO
        ENDDO

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled)

        IF ( DO_PARTLAYERS ) THEN
          UT = 0
          DO UTA = 1, N_USER_LEVELS
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT  = UT + 1
              N   = PARTLAYERS_LAYERIDX(UT)
              DT  = PARTLAYERS_VALUES(UT)
              XTD = DELTAU_VERT_INPUT(N) * DT
              PARTAU_VERT(UT) = XTD
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!mick fix 2/13/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if
!        necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        !(1) Divide by 2L+1 where L = 1 to get the asym par
        GREEKMAT_TOTAL(1,N,1) = GREEKMAT_TOTAL(1,N,1)/3.0D0
        !(2) Modify the asym par if necessary
        IF (GREEKMAT_TOTAL(1,N,1) > 0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (GREEKMAT_TOTAL(1,N,1) < -0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) >= 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) < 1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) < 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) > -1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = -1.0D-9
          SBUST(6) = .true.
        END IF
        !(3) Reconstruct the 1st-order phase func moment
        GREEKMAT_TOTAL(1,N,1) = 3.0D0*GREEKMAT_TOTAL(1,N,1)

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO
!READ(*,*)

!  Finish module

      RETURN
      END SUBROUTINE VLIDORT_DELTAMSCALE

!

      SUBROUTINE VLIDORT_SSALBINIT ( &
        NSTOKES, NLAYERS, &
        NMOMENTS, MUELLER_INDEX, &
        OMEGA_TOTAL, GREEKMAT_TOTAL, &
        OMEGA_GREEK )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  local variables

      INTEGER ::          N, L, O1, O2, OM
      DOUBLE PRECISION :: SUM

!  phase matrix-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          DO O1 = 1, NSTOKES
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              SUM = OMEGA_TOTAL(N)*GREEKMAT_TOTAL(L,N,OM)
              OMEGA_GREEK(L,N,O1,O2)   = SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_SSALBINIT

!

      SUBROUTINE VLIDORT_QSPREP ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        NLAYERS, &
        DO_SPECIALIST_OPTION_3, NLAYERS_CUTOFF, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        NBEAMS, DELTAU_VERT, &
        TAUGRID, DELTAU_SLANT, &
        LAYER_PIS_CUTOFF, &
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        LOCAL_CSZA )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_3
      INTEGER, INTENT (IN) ::           NLAYERS_CUTOFF
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::           NBEAMS
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

      INTEGER, INTENT (OUT) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )
      LOGICAL, INTENT (OUT) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: N, K, IB
      DOUBLE PRECISION :: S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)

      DOUBLE PRECISION :: TAUSLANT ( 0:MAXLAYERS, MAXBEAMS )

!  Specialist code
!  ---------------

!  Cutting off the solar source arbitrarily

      IF ( DO_SPECIALIST_OPTION_3 ) THEN
        DO IB = 1, NBEAMS
          LAYER_PIS_CUTOFF(IB) = NLAYERS_CUTOFF
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          LAYER_PIS_CUTOFF(IB) = NLAYERS
        ENDDO
      ENDIF

!  TOA  cosines

      DO IB = 1, NBEAMS
        LOCAL_CSZA(0,IB) = COS_SZANGLES(IB)
      ENDDO

!  plane-parallel case
!  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       DO IB = 1, NBEAMS
        SEC0 = ONE / COS_SZANGLES(IB)
!        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_CSZA(N,IB)      = COS_SZANGLES(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
            LOCAL_CSZA(N,IB)      = ZERO
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

!  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = ZERO
        DO N = 1, NLAYERS
          TAU = ZERO
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

!  set up the average secant formulation

        S_T_1 = ZERO
        S_T_0 = ONE
!        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) &
                                     / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

!  Set the Local solar zenith cosines
!  Distinguish between the refractive and non-refractive cases.

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              LOCAL_CSZA(N,IB) = SUN_SZA_COSINES(N,IB)
            ELSE
              LOCAL_CSZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              LOCAL_CSZA(N,IB) = COS_SZANGLES(IB)
            ELSE
              LOCAL_CSZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ENDIF

       ENDDO
      ENDIF

!  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        SOLAR_BEAM_OPDEP(IB) = ZERO
        DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        IF ( .NOT.DO_SPECIALIST_OPTION_3 ) THEN
          IF ( TAU_SOLAR(IB) .LT. MAX_TAU_SPATH ) THEN
            SOLAR_BEAM_OPDEP(IB) = DEXP( - TAU_SOLAR(IB) )
            DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
          ENDIF
        ENDIF
      ENDDO

! debug
!      do n = 1, nlayers
!        write(*,*)AVERAGE_SECANT(N,1),  INITIAL_TRANS(N,1)
!      enddo
!      pause

!  finish

      RETURN
      END SUBROUTINE VLIDORT_QSPREP

!

      SUBROUTINE VLIDORT_PREPTRANS ( &
        DO_SOLUTION_SAVING, NSTREAMS, &
        NLAYERS, DO_TOA_CONTRIBS, &
        LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, &
        USER_SECANTS, N_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        DELTAU_VERT, PARTAU_VERT, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_DISORDS_UTDN, T_DELT_MUBAR, &
        T_UTDN_MUBAR, T_UTUP_MUBAR, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, CUMTRANS, &
        ITRANS_USERM )

!  Prepare transmittances and transmittance factors

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_SOLUTION_SAVING
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_MUBAR &
          ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_MUBAR &
          ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER          :: N, UT, UM, IB, I
      DOUBLE PRECISION :: XT, SPHER, HELP

!  Transmittance factors for discrete ordinate streams
!  ===================================================

!  New code by R. Spurr, RT Solutions, 12 April 2005
!    Required for the solution saving option
!  Partial layer code added 30 December 2005 by R. Spurr

      IF ( DO_SOLUTION_SAVING ) THEN

!  whole layers

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            SPHER = DELTAU_VERT(N) / QUAD_STREAMS(I)
            IF ( SPHER .GT. MAX_TAU_QPATH ) THEN
              T_DELT_DISORDS(I,N) = ZERO
            ELSE
              T_DELT_DISORDS(I,N) = DEXP ( - SPHER )
            ENDIF
          ENDDO
        ENDDO

!  Atmosphere Partial layers

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = PARTAU_VERT(UT)
          DO I = 1, NSTREAMS
            HELP =  XT / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTDN(I,UT) = ZERO
            ELSE
              T_DISORDS_UTDN(I,UT) = DEXP(-HELP)
            ENDIF
            HELP = ( DELTAU_VERT(N) - XT ) / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTUP(I,UT) = ZERO
            ELSE
              T_DISORDS_UTUP(I,UT) = DEXP(-HELP)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!  Transmittance factors for average secant stream
!  ===============================================

!  start solar loop

      DO IB = 1, NBEAMS

!  Whole layer Transmittance factors
!  ---------------------------------

!  layer transmittance

       DO N = 1, NLAYERS
        IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
          T_DELT_MUBAR(N,IB) = ZERO
        ELSE
          SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
          IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
            T_DELT_MUBAR(N,IB) = ZERO
          ELSE
            T_DELT_MUBAR(N,IB) = DEXP ( - SPHER )
          ENDIF
        ENDIF
       ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  -----------------------------------------------------------------

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
           T_UTDN_MUBAR(UT,IB) = ZERO
           T_UTUP_MUBAR(UT,IB) = ZERO
         ELSE
           XT = PARTAU_VERT(UT)
           SPHER = XT * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTDN_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTDN_MUBAR(UT,IB) = DEXP ( - SPHER )
           ENDIF
           SPHER = ( DELTAU_VERT(N) - XT ) * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTUP_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTUP_MUBAR(UT,IB) = DEXP ( - SPHER )
           ENDIF
         ENDIF
        ENDDO

!  end solar beam loop

      ENDDO

!  Transmittances for User Streams
!  ===============================

!  return if not flagged

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  Initial transmittances divided by user streams
!  ----------------------------------------------

      DO IB = 1, NBEAMS
       DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
         ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
        ENDDO
       ENDDO
      ENDDO

!  Whole Layer transmittances

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        DO UM = 1, N_USER_STREAMS
         SPHER = DELTAU_VERT(N) * USER_SECANTS(UM)
         IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
          T_DELT_USERM(N,UM) = ZERO
         ELSE
          T_DELT_USERM(N,UM) = DEXP ( - SPHER )
         ENDIF
        ENDDO
       ENDIF
      ENDDO

!  Cumulative tranmsittances (TOA contribution functions)
!     New section, 27 Janaury 2010

      if ( DO_TOA_CONTRIBS ) THEN
        DO UM = 1, N_USER_STREAMS
          CUMTRANS(1,UM) = ONE
          DO N = 1, NLAYERS - 1
            CUMTRANS(N+1,UM) = CUMTRANS(N,UM) * T_DELT_USERM(N,UM)
          ENDDO
        ENDDO
      endif

!  Partial Layer transmittances for off-grid optical depths

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)
        XT = PARTAU_VERT(UT)
        DO UM = 1, N_USER_STREAMS
          SPHER = XT * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTDN_USERM(UT,UM) = ZERO
          ELSE
            T_UTDN_USERM(UT,UM) = DEXP ( - SPHER )
          ENDIF
          SPHER = ( DELTAU_VERT(N) - XT ) * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTUP_USERM(UT,UM) = ZERO
          ELSE
            T_UTUP_USERM(UT,UM) = DEXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PREPTRANS

!

      SUBROUTINE VLIDORT_SURFACE_DIRECTBEAM ( &
        DELTA_FACTOR, FOURIER_COMPONENT, &
        DO_REFRACTIVE_GEOMETRY, NSTOKES, &
        NSTREAMS, NLAYERS, &
        FLUX_FACTOR, &
        SZA_LOCAL_INPUT, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        BRDF_F_0, USER_BRDF_F, &
        USER_BRDF_F_0, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  &
        SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0, &
        FLUXVEC, COS_SZANGLES, &
        NBEAMS, N_USER_STREAMS, &
        MUELLER_INDEX, DO_USER_STREAMS, &
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
        ATMOS_ATTN, DIRECT_BEAM, &
        USER_DIRECT_BEAM )

!  Generalized to include the BRDF case as well
!   RT Solutions Inc., 26 July 2010

!    New Surface-Leaving arguments added, 17 May 2012

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) ::  DELTA_FACTOR
      INTEGER, INTENT (IN)          ::  FOURIER_COMPONENT
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )

!  New Surface-Leaving stuff 17 May 2012
      LOGICAL, INTENT (IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::            DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::  SOLAR_BEAM_OPDEP ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: X0_FLUX, X0_BOA, ATTN, REFLEC, SUM, SL
      INTEGER          :: I, UI, O1, O2, IB, M, OM

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            DIRECT_BEAM(I,IB,O1) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              USER_DIRECT_BEAM(UI,IB,O1) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Fourier component

      M = FOURIER_COMPONENT

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
        ELSE
         X0_BOA = COS_SZANGLES(IB)
        ENDIF

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!       X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR / DELTA_FACTOR

        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4             ! New
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

!  Lambertian case
!  ---------------

!  Set value for Fourier = 0
!  Natural light only, so only the first Stokes component is nonzero

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          REFLEC = ATTN * LAMBERTIAN_ALBEDO
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            DO I = 1, NSTREAMS
              DIRECT_BEAM(I,IB,1) = REFLEC
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                USER_DIRECT_BEAM(UI,IB,1) = REFLEC
              ENDDO
            ENDIF
          ENDIF
        ENDIF

!  Non-Lambertian case
!  -------------------

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN

!  Solar beam reflected into quad directions

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              SUM = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                SUM = SUM + FLUXVEC(O2) * BRDF_F_0(M,OM,I,IB)
              ENDDO
              DIRECT_BEAM(I,IB,O1) = ATTN * SUM

            ENDDO
          ENDDO

!  Solar beam reflected into User directions

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                SUM = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  SUM = SUM + FLUXVEC(O2)*USER_BRDF_F_0(M,OM,UI,IB)
                ENDDO
                USER_DIRECT_BEAM(UI,IB,O1) = ATTN * SUM
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ START
!  New Surface-Leaving stuff 17 May 2012
        IF ( DO_SURFACE_LEAVING ) THEN
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            DO O1 = 1, NSTOKES
              SL = SLTERM_ISOTROPIC(O1)
              DO I = 1, NSTREAMS
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_DIRECT_BEAM(UI,IB,O1) = USER_DIRECT_BEAM(UI,IB,O1) + SL
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                SL = SLTERM_F_0(M,O1,I,IB)
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  SL = USER_SLTERM_F_0(M,O1,UI,IB)
                  USER_DIRECT_BEAM(UI,IB,O1) = USER_DIRECT_BEAM(UI,IB,O1) + SL
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END 

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
      END SUBROUTINE VLIDORT_SURFACE_DIRECTBEAM

!

      SUBROUTINE VLIDORT_CHAPMAN ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        NLAYERS, N_SZANGLES, &
        SZANGLES, &
        EARTH_RADIUS, RFINDEX_PARAMETER, &
        HEIGHT_GRID, PRESSURE_GRID, &
        TEMPERATURE_GRID, FINEGRID, &
        SZA_LOCAL_INPUT, CHAPMAN_FACTORS, &
        FAIL, MESSAGE, TRACE )

!  This is the internal Chapman function calculation of the slant path
!  Chapman Factors required for slant-path optical thickness values.

!  This module calculates CHAPMAN_FACTORS internally inside VLIDORT,
!  saving you the job of doing it yourself, though you can input these
!  quantities, if the DO_CHAPMAN_FACTORS flag is not set.

!  The following options apply:

!   1. If the plane-parallel flag is on, no further inputs are required

!   2. If the plane-parallel flag is off, then Pseudo-spherical:
!       (a) Straight line geometry, must specify
!               Earth_radius, height grid
!       (b) Refractive geometry, must specify
!               Earth_radius, height grid
!               pressure grid, temperature grid

!  The logic will be checked before the module is called.

!  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

!    Based round a call to a pure geometry module which returns slant
!    path distances which was adapted for use in the Radiant model by
!    R. Spurr during an OCO L2 intensive April 24-29, 2005.

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) ::  SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) ::  RFINDEX_PARAMETER
      DOUBLE PRECISION, INTENT (IN) ::  HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER, INTENT (IN) ::           FINEGRID ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (OUT) :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  number of iterations (refractive case only)
!      This is debug output

      INTEGER          ::  ITERSAVE(MAXLAYERS)

!  other local variables

      INTEGER          :: IB
      DOUBLE PRECISION :: SUN0

!  get spherical optical depths
!  ----------------------------

!  start beam loop

!mick - added initialization
      CHAPMAN_FACTORS = ZERO
      SZA_LOCAL_INPUT = ZERO

      DO IB = 1, N_SZANGLES

        SUN0 = SZANGLES(IB)

        CALL BEAM_GEOMETRY_PREPARE &
           ( MAX_SZANGLES, MAXLAYERS, IB, NLAYERS, FINEGRID, &
             SUN0, EARTH_RADIUS, RFINDEX_PARAMETER, &
             DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
             HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, &
             CHAPMAN_FACTORS, SZA_LOCAL_INPUT, &
             ITERSAVE, FAIL, MESSAGE )

!  return if failed

        IF ( FAIL ) THEN
          TRACE = 'Geometry failure in VLIDORT_CHAPMAN'
          RETURN
        ENDIF

!  end beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHAPMAN

!

      SUBROUTINE BEAM_GEOMETRY_PREPARE ( &
        MAXBEAMS, MAXLAYERS, IBEAM, NLAYERS, FINEGRID, &
        SZA_GEOM_TRUE, REARTH, RFINDEX_PARAMETER, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        HEIGHTS, PRESSURES, TEMPERATURES, &
        CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT, &
        ITERSAVE, FAIL, MESSAGE )

      IMPLICIT NONE

!  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Coarse layering is input to the module. Values of Z, P, T  are
!  given at the layer boundaries, with the first value (index 0) at TOA.

!  The refractive geometry is assumed to start at the TOA level.

!  We also require the earth radius and the refractive index parameter
!   (For the Born-Wolf approximation)

!  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
!  In this case we do not require pressure and temperature information.
!  The calculation will then be for geometric rays.

!  The plane parallel Flag and the refractive geometry flag should not
!  both be true - this should be checked outside.

!  In the refracting case, fine-gridding of pressure and temperature is
!  done internally, temperature is interpolated linearly with height,
!  and pressure log-linearly. The refraction uses Snell's law rule.
!  Finelayer gridding assumes equidistant heights within coarse layers
!  but the number of fine layers can be varied

!  Output is specified at coarse layer boundaries

!  Module is stand-alone.

!  Reprogrammed for the OCO L2 algorithm
!   R. Spurr, RT Solutions, Inc.   April 27, 2005

!  Intended use in LIDORT and Radiant RT models.

!  Input arguments
!  ===============

!  input dimensioning

      INTEGER, INTENT (IN) ::          MAXLAYERS, MAXBEAMS

!  Beam index

      INTEGER, INTENT (IN) ::          IBEAM

!  number of coarse layers

      INTEGER, INTENT (IN) ::          NLAYERS

!  number of fine layers within coarse layers

      INTEGER, INTENT (IN) ::          FINEGRID(MAXLAYERS)

!  True solar zenith angle (degrees)

      DOUBLE PRECISION, INTENT (IN) :: SZA_GEOM_TRUE

!  Earth radius (km)

      DOUBLE PRECISION, INTENT (IN) :: REARTH

!  Refractive index parametaer (Born-Wolf approximation)

      DOUBLE PRECISION, INTENT (IN) :: RFINDEX_PARAMETER

!  flag for plane parallel case

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  flag for refractive geometry

      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY

!  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION, INTENT (IN) :: HEIGHTS (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: PRESSURES (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: TEMPERATURES (0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

!mick
!      DOUBLE PRECISION, INTENT (OUT) :: &
!        CHAPMAN_FACTORS (MAXLAYERS, MAXLAYERS, MAXBEAMS)
      DOUBLE PRECISION, INTENT (INOUT) :: &
        CHAPMAN_FACTORS (MAXLAYERS, MAXLAYERS, MAXBEAMS)

!  solar zenith angles at nadir

!mick
!      DOUBLE PRECISION, INTENT (OUT) :: &
!        SZA_LEVEL_OUTPUT (0:MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (INOUT) :: &
        SZA_LEVEL_OUTPUT (0:MAXLAYERS,MAXBEAMS)

!  number of iterations (refractive case only)
!   This is debug output

      INTEGER, INTENT (OUT)         :: ITERSAVE(MAXLAYERS)

!  output status

      LOGICAL, INTENT (OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE

!  Local variables
!  ===============

!  local dimensioning

      INTEGER, PARAMETER :: &
        LOCAL_MAXLAYERS = 141, &
        LOCAL_MAXFINELAYERS = 20

!  fine layer gridding for refraction

      DOUBLE PRECISION :: ZRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION :: PRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION :: TRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)

!  local height arrays

      DOUBLE PRECISION :: H(0:LOCAL_MAXLAYERS)
      DOUBLE PRECISION :: DELZ(LOCAL_MAXLAYERS)

!  help variables

      INTEGER          :: N, J, NRFINE, K, ITER, MAXF,IB
      LOGICAL          :: LOOP
      DOUBLE PRECISION :: GM_TOA, TH_TOA, MU_TOA, MU_NEXT
      DOUBLE PRECISION :: Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q
      DOUBLE PRECISION :: FU, FL, DEG_TO_RAD

      DOUBLE PRECISION :: LAYER_DIST, &
            MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK, &
            PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO, &
            RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL

!  Standard temperature (K) and pressure (mbar).

      DOUBLE PRECISION, PARAMETER :: &
        T_STANDARD = 273.16D0, &
        P_STANDARD = 1013.25D0, &
        STP_RATIO  = T_STANDARD / P_STANDARD

!  Loschmidt's number (particles/cm2/km).

      DOUBLE PRECISION, PARAMETER :: &
        RHO_STANDARD = 2.68675D+24

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0,IB) = 0.0D0
      DO N = 1, NLAYERS
        SZA_LEVEL_OUTPUT(N,IB) = 0.0D0
        DO K = 1, NLAYERS
          CHAPMAN_FACTORS(N,K,IB) = 0.0D0
        ENDDO
      ENDDO
      FAIL    = .FALSE.
      MESSAGE = ' '

!  check local dimensioning

      IF ( LOCAL_MAXLAYERS .LT. NLAYERS ) THEN
        MESSAGE = 'local coarse layer dimensioning insufficient'
        FAIL = .TRUE.
        RETURN
      ENDIF

!  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

!  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )
      STH2D  = 0.0D0

!  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!mick fix - moved inside if block from above
        !  Check fine layers do not exceed local dimensions assigned
        MAXF = 0
        DO N = 1, NLAYERS
          MAXF = MAX(MAXF,FINEGRID(N))
        ENDDO
        IF ( LOCAL_MAXFINELAYERS .LT. MAXF ) THEN
          MESSAGE = 'local fine layer dimensioning insufficient'
          FAIL = .TRUE.
          RETURN
        ENDIF

        Z0 = HEIGHTS(0)
        P0 = PRESSURES(0)
        T0 = TEMPERATURES(0)
        Q0 = DLOG(P0)
        DO N = 1, NLAYERS
          NRFINE = FINEGRID(N)
          LOCAL_SUBTHICK = DELZ(N) / DBLE(NRFINE)
          P1 = PRESSURES(N)
          Z1 = HEIGHTS(N)
          T1 = TEMPERATURES(N)
          Q1 = DLOG(P1)
          ZRFINE(N,0) = Z0
          PRFINE(N,0) = P0
          TRFINE(N,0) = T0
          DO J = 1, NRFINE - 1
            Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
            FL = ( Z0 - Z ) / DELZ(N)
            FU = 1.0d0 - FL
            Q  = FL * Q1 + FU * Q0
            T  = FL * T0 + FU * T1
            PRFINE(N,J) = DEXP (  Q )
            TRFINE(N,J) = T
            ZRFINE(N,J) = Z
          ENDDO
          PRFINE(N,NRFINE) = P1
          TRFINE(N,NRFINE) = T1
          ZRFINE(N,NRFINE) = Z1
!              write(*,'(i3,11F10.4)')N,(PRFINE(N,J),J=0,NRFINE)
          Z0 = Z1
          P0 = P1
          T0 = T1
          Q0 = Q1
        ENDDO
      ENDIF

!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Refractive Geometry case
!  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  Start value of SZA cosine

        MU_PREV = MU_TOA

!  start layer loop

        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1 = DASIN(SINTH1)
          PHI_0 = TH_TOA - STH1
          NRFINE = FINEGRID(N)

!  iteration loop

          ITER = 0
          LOOP = .TRUE.
          DO WHILE (LOOP.AND.ITER.LT.100)
            ITER = ITER + 1
            PHI_CUM = 0.0D0
            RE_UPPER = ZRFINE(1,0) + REARTH
            RATIO  = PRFINE(1,0) * STP_RATIO / TRFINE(1,0)
            REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
            SNELL = REFRAC * RE_UPPER * SINTH1
            DO K = 1, N
              LAYER_DIST = 0.0D0
              LOCAL_SUBTHICK = DELZ(K) / DBLE(NRFINE)
              DO J = 0, NRFINE - 1
                RATIO  = PRFINE(K,J) * STP_RATIO / TRFINE(K,J)
                REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
                RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                SINTH2 = SNELL/ (REFRAC * RE_UPPER )
                IF ( SINTH2.GT.1.0D0 ) SINTH2 = 1.0D0
                STH2 = DASIN(SINTH2)
                SINTH2D = RE_UPPER * SINTH2 / RE_LOWER
                IF ( SINTH2D .GT. 1.0D0 ) THEN
                  MESSAGE = 'refraction yields angles > 90 some levels'
                  FAIL = .TRUE.
                  RETURN
                ENDIF
                STH2D = DASIN(SINTH2D)
                PHI = STH2D - STH2
                SINPHI = DSIN(PHI)
                PHI_CUM = PHI_CUM + PHI
                DIST = RE_UPPER * SINPHI / SINTH2D
                LAYER_DIST = LAYER_DIST +  DIST
                RE_UPPER = RE_LOWER
              ENDDO
              CHAPMAN_FACTORS(N,K,IB) = LAYER_DIST / DELZ(K)
            ENDDO

!  examine convergence

            DELPHI = PHI_0 - PHI_CUM
            LOOP = (DABS(DELPHI/PHI_CUM).GT.1.0D-4)

!  Fudge factors to speed up the iteration

            IF ( SZA_GEOM_TRUE .GT. 88.7D0 ) THEN
              STH1 = STH1 + 0.1 * DELPHI
              PHI_0 = TH_TOA - STH1
            ELSE IF ( SZA_GEOM_TRUE .LT. 80.0D0 ) THEN
              PHI_0 = PHI_CUM
              STH1 = TH_TOA - PHI_0
            ELSE
              STH1 = STH1 + 0.3 * DELPHI
              PHI_0 = TH_TOA - STH1
            ENDIF
            SINTH1 = DSIN(STH1)

          ENDDO

!  failure

          IF ( LOOP ) THEN
            MESSAGE = 'refractive iteration not converged'
            FAIL = .TRUE.
            RETURN
          ENDIF

!  Update and save angle output

          MU_NEXT = DCOS(STH2D)
          MU_PREV = MU_NEXT
          SZA_LEVEL_OUTPUT(N,IB) = DACOS(MU_NEXT) / DEG_TO_RAD
          ITERSAVE(N) = ITER

        ENDDO

!  Straight line geometry
!  ======================

      ELSE

        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1   = DASIN(SINTH1)
          RE_UPPER = H(0)

!  solar zenith angles are all the same = input value

          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N

          DO K = 1, N

!  sine-rule; PHI = earth-centered angle

            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = DASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = DSIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

!  re-set

            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2

          ENDDO

!  finish main layer loop

        ENDDO

!  Finish

      ENDIF

!  end of routine

      RETURN
      END SUBROUTINE BEAM_GEOMETRY_PREPARE

      END MODULE vlidort_miscsetups_module

