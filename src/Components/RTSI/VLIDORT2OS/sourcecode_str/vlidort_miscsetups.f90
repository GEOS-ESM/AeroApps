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
! #            VLIDORT_MISCSETUPS (master, calling:)            #
! #              VLIDORT_DELTAMSCALE                            #
! #              VLIDORT_SSALBINIT                              #
! #              VLIDORT_QSPREP                                 #
! #              VLIDORT_PREPTRANS                              #
! #                                                             #
! #            VLIDORT_DIRECTBEAM                               #
! #            VLIDORT_PIMATRIX_SETUP                           #
! #            VLIDORT_PIMATRIX_SETUP_OMP  (Version 2.7)        #
! #                                                             #
! ###############################################################


      MODULE vlidort_miscsetups_module

      PRIVATE
      PUBLIC :: VLIDORT_MISCSETUPS, &
                VLIDORT_DIRECTBEAM, &
                VLIDORT_PIMATRIX_SETUP, &
                VLIDORT_PIMATRIX_SETUP_OMP

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
        DO_OBSERVATION_GEOMETRY, &
        USER_SECANTS, N_PARTLAYERS, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        OMEGA_TOTAL, DELTAU_VERT, &
        PARTAU_VERT, GREEKMAT_TOTAL, &
        TAUGRID, DELTAU_SLANT, SOLARBEAM_BOATRANS, &   ! Rob fix 11/17/14 added argument
        TRUNC_FACTOR, FAC1, &
        OMEGA_GREEK, LAYER_PIS_CUTOFF, &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
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
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
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
      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
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

!  Rob fix 11/17/14. Added Argument

      DOUBLE PRECISION, INTENT (OUT) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

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
        TAUGRID, DELTAU_SLANT, SOLARBEAM_BOATRANS, &   ! Rob fix 11/17/14 added argument
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
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        LOCAL_CSZA )

!  Transmittances and Transmittance factors

      CALL VLIDORT_PREPTRANS ( &
        DO_SOLUTION_SAVING, &
        NSTREAMS, NLAYERS, DO_TOA_CONTRIBS, &
        LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
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
        TAUGRID, DELTAU_SLANT, SOLARBEAM_BOATRANS, &   ! Rob fix 11/17/14 added argument
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

!  Rob fix 11/17/14. Added Argument for Diagnostic output

      DOUBLE PRECISION, INTENT (OUT) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

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

!mick fix 7/23/2014 - initialized for packing
      DELTAU_SLANT = ZERO

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!rob fix 11/17/2014 - initialize Solar Beam Transmittance for output, then calculate
!    Must be UNSCALED, so you need it now.

      SOLARBEAM_BOATRANS = ZERO
      DO IB = 1, NBEAMS
        DELS = SUM(DELTAU_SLANT(NLAYERS,1:NLAYERS,IB))
        IF ( DELS.LT. MAX_TAU_SPATH ) SOLARBEAM_BOATRANS(IB) = EXP ( - DELS )
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

!Rob Fix for 2OS correction - initialize
      TRUNC_FACTOR = zero
      OMEGA_TOTAL  = zero
      DELTAU_VERT  = zero

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
!mick fix 1/21/2013 -  added if structure and else section
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = &
                      GREEKMAT_TOTAL_INPUT(L,N,K1)
                GREEKMAT_TOTAL(L,N,K2) = &
                      GREEKMAT_TOTAL_INPUT(L,N,K2)
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = &
                    GREEKMAT_TOTAL_INPUT(L,N,1)
            ENDIF
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
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
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
      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
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
        TRANS_SOLAR_BEAM(IB) = ZERO
        DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        IF ( .NOT.DO_SPECIALIST_OPTION_3 ) THEN
          IF ( TAU_SOLAR(IB) .LT. MAX_TAU_SPATH ) THEN
            TRANS_SOLAR_BEAM(IB) = DEXP( - TAU_SOLAR(IB) )
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
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
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
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
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

      INTEGER          :: N, UT, UM, IB, I, LUM
      DOUBLE PRECISION :: XT, SPHER, HELP

!  Local user index

      LUM = 1

!  Transmittance factors for discrete ordinate streams
!  ===================================================

!  New code by R. Spurr, RT Solutions, 12 April 2005
!    Required for the solution saving option
!  Partial layer code added 30 December 2005 by R. Spurr

!mick fix 7/23/2014 - initialized for packing
        T_DELT_DISORDS = ZERO
        T_DISORDS_UTDN = ZERO
        T_DISORDS_UTUP = ZERO

        ITRANS_USERM = ZERO
        T_DELT_USERM = ZERO
        T_UTDN_USERM = ZERO
        T_UTUP_USERM = ZERO

        CUMTRANS = ZERO

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

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
          DO UM = 1, N_USER_STREAMS
           ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
          ENDDO
         ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
           ITRANS_USERM(N,LUM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(IB)
         ENDDO
        ENDDO
      ENDIF

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

      SUBROUTINE VLIDORT_DIRECTBEAM ( &
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
        DO_OBSERVATION_GEOMETRY, &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
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
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC &
          ( MAXSTOKES, MAXBEAMS )
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
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      DOUBLE PRECISION, INTENT (IN) ::  TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: X0_FLUX, X0_BOA, ATTN, REFLEC, SUM, SL, Help
      INTEGER          :: I, UI, O1, O2, IB, M, OM, LUI

!  Initialize
!  ----------

!  Local user index

      LUI = 1

!  Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            DIRECT_BEAM(I,IB,O1) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UI = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                USER_DIRECT_BEAM(UI,IB,O1) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              USER_DIRECT_BEAM(LUI,IB,O1) = ZERO
            ENDDO
          ENDIF
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
        ATTN           = X0_FLUX * TRANS_SOLAR_BEAM(IB)
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
              IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_DIRECT_BEAM(UI,IB,1) = REFLEC
                ENDDO
              ELSE
                USER_DIRECT_BEAM(LUI,IB,1) = REFLEC
              ENDIF
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
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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
            ELSE
              DO O1 = 1, NSTOKES
                SUM = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  SUM = SUM + FLUXVEC(O2)*USER_BRDF_F_0(M,OM,LUI,IB)
                ENDDO
                USER_DIRECT_BEAM(LUI,IB,O1) = ATTN * SUM
              ENDDO
            ENDIF
          ENDIF

        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ START
!  New Surface-Leaving stuff 17 May 2012

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

        IF ( DO_SURFACE_LEAVING ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            DO O1 = 1, NSTOKES
              SL = SLTERM_ISOTROPIC(O1,IB) * HELP
              DO I = 1, NSTREAMS
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                  DO UI = 1, N_USER_STREAMS
                    USER_DIRECT_BEAM(UI,IB,O1) = USER_DIRECT_BEAM(UI,IB,O1) + SL
                  ENDDO
                ELSE
                  USER_DIRECT_BEAM(LUI,IB,O1) = USER_DIRECT_BEAM(LUI,IB,O1) + SL
                ENDIF
              ENDIF
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                SL = SLTERM_F_0(M,O1,I,IB) * HELP
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                  DO UI = 1, N_USER_STREAMS
                    SL = USER_SLTERM_F_0(M,O1,UI,IB)* HELP
                    USER_DIRECT_BEAM(UI,IB,O1) = USER_DIRECT_BEAM(UI,IB,O1) + SL
                  ENDDO
                ELSE
                  SL = USER_SLTERM_F_0(M,O1,LUI,IB)* HELP
                  USER_DIRECT_BEAM(LUI,IB,O1) = USER_DIRECT_BEAM(LUI,IB,O1) + SL
                ENDIF
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
      END SUBROUTINE VLIDORT_DIRECTBEAM

!

      SUBROUTINE VLIDORT_PIMATRIX_SETUP ( &
        FOURIER, &
        DO_REFRACTIVE_GEOMETRY, NSTOKES, &
        NSTREAMS, NLAYERS, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        QUAD_STREAMS, NMOMENTS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, USER_STREAMS, &
        MUELLER_INDEX, DMAT, &
        PI_XQP, PI_XQM, &
        PI_XUP, PI_XUM, &
        PI_X0P, &
        PI_XQM_POST, PI_XQM_PRE, &
        PI_XQP_PRE, PI_XUM_POST, &
        PI_XUP_PRE )

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )

!  Indices and Dmatrix from VLIDORT_DERIVE_INPUTS

      INTEGER         , INTENT (IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)  :: DMAT ( MAXSTOKES, MAXSTOKES )

!  Output (generalized spherical functions)

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI &
          ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PHI ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )

      DOUBLE PRECISION, SAVE :: DF_L ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_2LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LSQM4 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: ZHELP ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PLEG20 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: M2X_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: XUMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X_RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_KM ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X2LP1 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  local variables

!     INTEGER          :: LOCAL_NSTOKES
      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, M1, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

!  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

!  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = DSQRT(PHI(L,M))
      ENDDO

!  Additional saved quantities (to commons in VLIDORT_PISETUP.VARS)
!  ----------------------------------------------------------------

!  Only for the first fundamental harmonic

      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DBLE(L)
          DF_LP1(L)  = DBLE(L+1)
          DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = DSQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = DSQRT(DBLE((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

      DO I = 1, NSTREAMS
        XM(I) = QUAD_STREAMS(I)
      ENDDO
      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          XM(I+NSTREAMS) = USER_STREAMS(I)
        ENDDO
      ENDIF

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO N = 1, NLAYERS
!          DO IB = 1, NBEAMS
!            IP = NS + NBEAMS * (N-1) + IB
!            XM(IP) = SUN_SZA_COSINES(N,IB)
!          ENDDO
!        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          IP = IB + NS
          XM(IP) = COS_SZANGLES(IB)
        ENDDO
      ENDIF

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

      IF ( M .EQ. 0 )  THEN
        DO I = 1, NSA
          XSQ = XM(I) * XM(I)
          PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          UMXSQ(I)         = ONE - XSQ
          XUMXSQ(I)        = UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           RT_UMXSQ(I)      = ZERO
           X_RT_UMXSQ(I)    = ZERO
           UPXSQ_D_UMXSQ(I) = ZERO
           M2X_D_UMXSQ(I)   = ZERO
          ELSE
           RT_UMXSQ(I)      = DSQRT(UMXSQ(I))
           X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
           UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
           M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * DF_2LP1(L)
          ENDDO
        ENDDO

!  D-matrix and Mueller indices. NOW DONE in VLIDORT_DERIVE_INPUTS
!        DO SK1 = 1, MAXSTOKES
!          DO SK2 = 1, MAXSTOKES
!            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
!            DMAT(SK1,SK2) = ZERO
!          ENDDO
!          IF ( SK1.GT.2) THEN
!            DMAT(SK1,SK1) = -ONE
!          ELSE
!            DMAT(SK1,SK1) = ONE
!          ENDIF
!          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
!        ENDDO

      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          XLM_DIAG(L,1) = DF_LP1(L)
          XLM_DIAG(L,2) = DF_RT_LP3XLM1(L)
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            XLM_DIAG(L,SK1) = ONE/XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          YLM_DIAG(L,1) = DF_L(L)
          YLM_DIAG(L,2) = DF_LSQM4(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            XLM_DIAG(L,2) = ZERO
          ELSE
            XLM_DIAG(L,2) = DF_RT_LP3XLM1(L) * DF_LP1MM / DF_LP1(L)
            XLM_DIAG(L,2) = ONE / XLM_DIAG(L,2)
          ENDIF
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------

!  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!mick - alternate PI initialization
      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PIMM_11(I) =   THREE     * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) =   RCONST_21
            ELSE
              PIMM_KM(I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = DF_2LP1(M-1)
          H2 = DF_LP1(M-1)*H1/DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
            PIMM_11(I) = H1 * RT_UMXSQ(I) * PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) = ZERO
            ELSE
              PIMM_KM(I) = H2 * RT_UMXSQ(I) * PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF

        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PIMM_KM(I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * PI_XQP(L,I,SK3,SK2)
                H2 = H2 + PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              PI_XQP_PRE(L,I,SK1,SK2)  = H1
              PI_XQM_PRE(L,I,SK1,SK2)  = H1 * PISIGN(L,M)
              PI_XQM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              PI_XQM(L,I,SK1,SK2)  = H1 * PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  .. at positive user_defined angles

        IF ( DO_USER_STREAMS ) THEN

          DO I = 1, N_USER_STREAMS
            I1 = I + NSTREAMS

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_XUP(L,I,SK1,SK2) = PI(L,I1,SK1,SK2) * PINORM(L,M)
              ENDDO
            ENDDO

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES

                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PI_XUP(L,I,SK1,SK3)
                ENDDO
                PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PI_XUM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PI_XUM(L,I,SK1,SK2) = H1 * PISIGN(L,M)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

!  .. at solar zenith angles

!  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              DO SK1 = 1, NSTOKES
!                DO SK2 = 1, NSTOKES
!                  PI_X0P(L,IB,N,SK1,SK2) =
!     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = &
                     PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PIMATRIX_SETUP

!

      SUBROUTINE VLIDORT_PIMATRIX_SETUP_OMP ( &
        FOURIER, &
        DO_REFRACTIVE_GEOMETRY, NSTOKES, &
        NSTREAMS, NLAYERS, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        QUAD_STREAMS, NMOMENTS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, USER_STREAMS, &
        MUELLER_INDEX, DMAT, &
        PIMM_11, PIMM_KM, &
        PI_XQP, PI_XQM, &
        PI_XUP, PI_XUM, &
        PI_X0P, &
        PI_XQM_POST, PI_XQM_PRE, &
        PI_XQP_PRE, PI_XUM_POST, &
        PI_XUP_PRE )

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )

!  Indices and Dmatrix from VLIDORT_DERIVE_INPUTS

      INTEGER         , INTENT (IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)  :: DMAT ( MAXSTOKES, MAXSTOKES )

!  InOut variables

      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Output (generalized spherical functions)

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI &
          ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PHI ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )

      DOUBLE PRECISION, SAVE :: DF_L ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_2LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LSQM4 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: ZHELP ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PLEG20 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: M2X_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: XUMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X_RT_UMXSQ ( MAX_ALLSTRMS_P1 )
!mick fix 7/28/2014 - these two defined InOut for now to make VLIDORT threadsafe
      !DOUBLE PRECISION :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      !DOUBLE PRECISION :: PIMM_KM ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X2LP1 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  local variables

!     INTEGER          :: LOCAL_NSTOKES
      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, M1, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

!  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

!  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = DSQRT(PHI(L,M))
      ENDDO

!  Additional saved quantities
!  ---------------------------

!  Only for the first fundamental harmonic

!mick fix 7/28/2014 - do for each fourier right now to make VLIDORT threadsafe
      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DBLE(L)
          DF_LP1(L)  = DBLE(L+1)
          DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = DSQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = DSQRT(DBLE((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

      DO I = 1, NSTREAMS
        XM(I) = QUAD_STREAMS(I)
      ENDDO
      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          XM(I+NSTREAMS) = USER_STREAMS(I)
        ENDDO
      ENDIF

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO N = 1, NLAYERS
!          DO IB = 1, NBEAMS
!            IP = NS + NBEAMS * (N-1) + IB
!            XM(IP) = SUN_SZA_COSINES(N,IB)
!          ENDDO
!        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          IP = IB + NS
          XM(IP) = COS_SZANGLES(IB)
        ENDDO
      ENDIF

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

!mick fix 7/28/2014 - do for each fourier right now to make VLIDORT threadsafe
      IF ( M .EQ. 0 )  THEN
        DO I = 1, NSA
          XSQ = XM(I) * XM(I)
          PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          UMXSQ(I)         = ONE - XSQ
          XUMXSQ(I)        = UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           RT_UMXSQ(I)      = ZERO
           X_RT_UMXSQ(I)    = ZERO
           UPXSQ_D_UMXSQ(I) = ZERO
           M2X_D_UMXSQ(I)   = ZERO
          ELSE
           RT_UMXSQ(I)      = DSQRT(UMXSQ(I))
           X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
           UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
           M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * DF_2LP1(L)
          ENDDO
        ENDDO

!  D-matrix and Mueller indices. NOW DONE in VLIDORT_DERIVE_INPUTS
!        DO SK1 = 1, MAXSTOKES
!          DO SK2 = 1, MAXSTOKES
!            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
!            DMAT(SK1,SK2) = ZERO
!          ENDDO
!          IF ( SK1.GT.2) THEN
!            DMAT(SK1,SK1) = -ONE
!          ELSE
!            DMAT(SK1,SK1) = ONE
!          ENDIF
!          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
!        ENDDO

      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          XLM_DIAG(L,1) = DF_LP1(L)
          XLM_DIAG(L,2) = DF_RT_LP3XLM1(L)
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            XLM_DIAG(L,SK1) = ONE/XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          YLM_DIAG(L,1) = DF_L(L)
          YLM_DIAG(L,2) = DF_LSQM4(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            XLM_DIAG(L,2) = ZERO
          ELSE
            XLM_DIAG(L,2) = DF_RT_LP3XLM1(L) * DF_LP1MM / DF_LP1(L)
            XLM_DIAG(L,2) = ONE / XLM_DIAG(L,2)
          ENDIF
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------

!  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!mick - alternate PI initialization
      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PIMM_11(I) =   THREE * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) =   RCONST_21
            ELSE
              PIMM_KM(I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = DF_2LP1(M-1)
          H2 = DF_LP1(M-1)*H1/DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
            PIMM_11(I) = H1 * RT_UMXSQ(I) * PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) = ZERO
            ELSE
              PIMM_KM(I) = H2 * RT_UMXSQ(I) * PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF

        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PIMM_KM(I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * PI_XQP(L,I,SK3,SK2)
                H2 = H2 + PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              PI_XQP_PRE(L,I,SK1,SK2)  = H1
              PI_XQM_PRE(L,I,SK1,SK2)  = H1 * PISIGN(L,M)
              PI_XQM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              PI_XQM(L,I,SK1,SK2)  = H1 * PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  .. at positive user_defined angles

        IF ( DO_USER_STREAMS ) THEN

          DO I = 1, N_USER_STREAMS
            I1 = I + NSTREAMS

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_XUP(L,I,SK1,SK2) = PI(L,I1,SK1,SK2) * PINORM(L,M)
              ENDDO
            ENDDO

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES

                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PI_XUP(L,I,SK1,SK3)
                ENDDO
                PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PI_XUM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PI_XUM(L,I,SK1,SK2) = H1 * PISIGN(L,M)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

!  .. at solar zenith angles

!  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              DO SK1 = 1, NSTOKES
!                DO SK2 = 1, NSTOKES
!                  PI_X0P(L,IB,N,SK1,SK2) =
!     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = &
                     PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PIMATRIX_SETUP_OMP

      END MODULE vlidort_miscsetups_module

