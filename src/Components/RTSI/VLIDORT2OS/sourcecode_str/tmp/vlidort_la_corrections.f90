! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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
! # --------------------------                                  #
! #                                                             #
! #              VLIDORTSS_L_ZMATRICES                          #
! #              SSCORR_OUTGOING_L_ZMATRIX                      #
! #                                                             #
! #      Version 2.4R ------- Notes -------                     #
! #                                                             #
! #            Additional correction to ZMATRIX and FMATRIX     #
! #            routines for Azimuth > 180. HvdM, Eqs(94/95)     #
! #            Implemented by V. Natraj and R. Spurr, 5/1/09    #
! #                                                             #
! ###############################################################


      MODULE vlidort_la_corrections

      PRIVATE
      PUBLIC :: VLIDORTSS_L_ZMATRICES, &
                SSCORR_OUTGOING_L_ZMATRIX

      CONTAINS

      SUBROUTINE VLIDORTSS_L_ZMATRICES ( &
        DO_SUNLIGHT, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NGREEKMOMS, NV, NV_PARAMETERS, &
        N_GEOMETRIES, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
        DO_SCATMAT_VARIATION, LAYER_MAXMOMENTS, GREEKMAT, L_GREEKMAT, &
        DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL, VSIGN, VZA_OFFSETS, &
        CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, &
        L_ZMAT )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  -----

!  CONTROL

      LOGICAL, INTENT (IN) ::          DO_SUNLIGHT
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES, NGREEKMOMS, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

!  +/- SIGN, OFFSETS

      DOUBLE PRECISION, INTENT (IN) :: VSIGN
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAXBEAMS, MAX_USER_STREAMS )

!  SCATTERING INPUT INFORMATION

      LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )
      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  TRUNCATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION, INTENT (IN) :: SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  ZENITH ANGLE COSINES/SINES, AZIMUTH ANGLE COSINES

      DOUBLE PRECISION, INTENT (IN) :: CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN) :: STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN) :: CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT (IN) :: SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT (IN) :: CPHI   (MAX_USER_RELAZMS)

!  AZIMUTH ANGLES. ADDED V2.4R. FOR PHI > 180 CASE

      DOUBLE PRECISION, INTENT (IN) :: PHI    (MAX_USER_RELAZMS)

!  OUTPUT
!  ------

!  Z-MATRIX LINEARIZED

      DOUBLE PRECISION, INTENT (INOUT) :: L_ZMAT &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 4, 4 )

!  LOCAL
!  -----

!  ROTATION ANGLE COSINES/SINES

      DOUBLE PRECISION :: C1 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: S1 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: C2 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: S2 (MAX_GEOMETRIES)

!  LOCAL F-MATRIX

      DOUBLE PRECISION :: L_FMAT ( MAX_ATMOSWFS, MAX_GEOMETRIES, 6 )

!  HELP VARIABLES

      INTEGER          :: IB, UM, IA, V, Q, GREEKMAT_INDEX(6)
      DOUBLE PRECISION :: COSSCAT, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION :: UUU(MAX_GEOMETRIES)

!  GENERALIZED SPHERICAL FUNCTIONS
!    P2P2 AND P2M2 ADDED 20 MARCH 2006

      DOUBLE PRECISION :: P00(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P02(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P2P2(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P2M2(MAX_GEOMETRIES,2)

      INTEGER          :: K, L, LNEW, LOLD, ITMP
      INTEGER          :: INDEX_11, INDEX_12, INDEX_34
      INTEGER          :: INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION :: DL, QROOT6, FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION :: TMP1, TMP2, SUM23, DIF23
      DOUBLE PRECISION :: FL2, FLL1, PERLL4, WFACT, QFAC, DL1
      DOUBLE PRECISION :: HELP1, HELP2, FDNL1, DNL1, FACT
      DOUBLE PRECISION :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

      DOUBLE PRECISION :: GK11, L_GK11(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK12, L_GK12(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK44, L_GK44(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK34, L_GK34(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK22, L_GK22(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK33, L_GK33(MAX_ATMOSWFS)

!  INDEXING KEY

!      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
!      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
!      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
!      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
!      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
!      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      SQL41 = ZERO

!  GEOMETRICAL QUANTITIES
!  ----------------------

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN

        DO IB = 1, NBEAMS
          DO UM = 1, N_USER_STREAMS
            DO IA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + IA

!  COSINE SCATTER ANGLE (THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING

              COSSCAT = VSIGN * CTHETA(NV,IB) * CALPHA(UM) + &
                                STHETA(NV,IB) * SALPHA(UM) * CPHI(IA)
              UUU(V)  = COSSCAT

!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

              HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
              IF ( HELP_SINSCAT.LE.ZERO ) THEN
                SINSCAT = 1.0D-12
              ELSE
                SINSCAT = DSQRT ( HELP_SINSCAT )
              ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

              IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
                CSIG1 = ZERO
                CSIG2 = ZERO
              ELSE
                IF ( STHETA(1,IB) .EQ. ZERO ) THEN
                  CSIG1 = -  CPHI(IA)
                ELSE
                  CSIG1   = ( - VSIGN * CALPHA(UM) + CTHETA(NV,IB)*COSSCAT ) &
                            / SINSCAT / STHETA(NV,IB)
                ENDIF
                IF ( SALPHA(UM) .EQ. ZERO ) THEN
                  CSIG2 = -  CPHI(IA)
                ELSE
                  CSIG2   = ( - CTHETA(NV,IB) + VSIGN*CALPHA(UM)*COSSCAT ) &
                             / SINSCAT / SALPHA(UM)
                ENDIF
              ENDIF
!      WRITE(7,'(4I5,1P4E15.7)')V,IB,UM,IA,CSIG1,CSIG2,SINSCAT

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

              IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
              IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

              IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
              IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)
!    ROTATION SINES AND COSINES

              CSIG1_2 = TWO * CSIG1
              CSIG2_2 = TWO * CSIG2
              IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
                SSIG1 = ZERO
              ELSE
                SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
              ENDIF
              IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
                SSIG2 = ZERO
              ELSE
                SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
              ENDIF

!  FOR RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1 AND S2
!  SEE H/VDM, EQS. 94-95. V. NATRAJ AND R. SPURR, 01 MAY 2009.

              C1(V) = CSIG1_2 * CSIG1 - ONE
              C2(V) = CSIG2_2 * CSIG2 - ONE

              IF (PHI(IA) .LE. 180.D0) THEN
                S1(V) = CSIG1_2 * SSIG1
                S2(V) = CSIG2_2 * SSIG2
              ELSE
                S1(V) = -CSIG1_2 * SSIG1
                S2(V) = -CSIG2_2 * SSIG2
              ENDIF

!  END GEOMETRY LOOPS

            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO V = 1, N_GEOMETRIES

!  COSINE SCATTER ANGLE (THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING

          COSSCAT = VSIGN * CTHETA(NV,V) * CALPHA(V) + &
                            STHETA(NV,V) * SALPHA(V) * CPHI(V)
          UUU(V)  = COSSCAT

!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

          HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
          IF ( HELP_SINSCAT.LE.ZERO ) THEN
            SINSCAT = 1.0D-12
          ELSE
            SINSCAT = DSQRT ( HELP_SINSCAT )
          ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

          IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
            CSIG1 = ZERO
            CSIG2 = ZERO
          ELSE
            IF ( STHETA(1,V) .EQ. ZERO ) THEN
              CSIG1 = -  CPHI(V)
            ELSE
              CSIG1   = ( - VSIGN * CALPHA(V) + CTHETA(NV,V)*COSSCAT ) &
                        / SINSCAT / STHETA(NV,V)
            ENDIF
            IF ( SALPHA(V) .EQ. ZERO ) THEN
              CSIG2 = -  CPHI(V)
            ELSE
              CSIG2   = ( - CTHETA(NV,V) + VSIGN*CALPHA(V)*COSSCAT ) &
                         / SINSCAT / SALPHA(V)
            ENDIF
          ENDIF
!    WRITE(7,'(4I5,1P4E15.7)')V,IB,UM,IA,CSIG1,CSIG2,SINSCAT

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

          IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
          IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

          IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
          IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)
!    ROTATION SINES AND COSINES

          CSIG1_2 = TWO * CSIG1
          CSIG2_2 = TWO * CSIG2
          IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
            SSIG1 = ZERO
          ELSE
            SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
          ENDIF
          IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
            SSIG2 = ZERO
          ELSE
            SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
          ENDIF

!  FOR RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1 AND S2
!  SEE H/VDM, EQS. 94-95. V. NATRAJ AND R. SPURR, 01 MAY 2009.

          C1(V) = CSIG1_2 * CSIG1 - ONE
          C2(V) = CSIG2_2 * CSIG2 - ONE

          IF (PHI(V) .LE. 180.D0) THEN
            S1(V) = CSIG1_2 * SSIG1
            S2(V) = CSIG2_2 * SSIG2
          ELSE
            S1(V) = -CSIG1_2 * SSIG1
            S2(V) = -CSIG2_2 * SSIG2
          ENDIF

!  END GEOMETRY LOOP

        ENDDO

      ENDIF

!  LINEARIZED F-MATRICES
!  ---------------------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! INITIALISE F-MATRIX

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO Q = 1, NV_PARAMETERS
            L_FMAT(Q,V,K) = ZERO
          END DO
        END DO
      END DO

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  SET THE LOCAL GREEK MATRIX ELEMENTS THAT YOU NEED
!   44 AND 34 ARE NOT REQUIRED WITH NATURAL SUNLIGHT (DEFAULT HERE)
!   22 AND 33 REQUIRED FOR NON-MIE SPHEROIDAL PARTICLES

        DO Q = 1, NV_PARAMETERS

         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( DO_SSCORR_TRUNCATION ) THEN

           DNL1  = DBLE(2*L + 1 )
           FDNL1 = SSFDEL(NV) * DNL1
           FACT  = ONE - SSFDEL(NV)

!  (1,1) AND (1,2) COMPONENTS

           GK11 = ( GREEKMAT(L,NV,1) - FDNL1 ) / FACT
           HELP1 = L_GREEKMAT(Q,L,NV,1) * GREEKMAT(L,NV,1)
           HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(NV,Q)
           L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
           IF ( NSTOKES.GT.1 ) THEN
             GK12 =   GREEKMAT(L,NV,2) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
             HELP2 = GK12 * L_SSFDEL(NV,Q)
             L_GK12(Q) = ( HELP1 + HELP2 ) / FACT
           ENDIF

!  NON-SUNLIGHT CASES. (4,4) AND (3,4)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
!             GK44 = ( GREEKMAT(L,NV,16) - FDNL1 ) / FACT
!             GK34 =   GREEKMAT(L,NV,12) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
!             HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
!             HELP2 = GK34 * L_SSFDEL(NV,Q)
!             L_GK34(Q) = ( HELP1 + HELP2 ) / FACT
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
             GK44 = ( GREEKMAT(L,NV,16) - FDNL1 ) / FACT
             GK34 =   GREEKMAT(L,NV,12) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
             HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
             HELP2 = GK34 * L_SSFDEL(NV,Q)
             L_GK34(Q) = ( HELP1 + HELP2 ) / FACT
           ELSE
             L_GK44(Q) = ZERO
             L_GK34(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (2,2) AND (3,3)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT  ) THEN
!             GK22 = ( GREEKMAT(L,NV,6)  - FDNL1 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,6) * GREEKMAT(L,NV,6)
!             HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
!             GK33 = ( GREEKMAT(L,NV,11) - FDNL1 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
!             HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK33(Q) = ( HELP1 + HELP2 ) / FACT
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT  ) THEN
             GK22 = ( GREEKMAT(L,NV,6)  - FDNL1 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,6) * GREEKMAT(L,NV,6)
             HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
             GK33 = ( GREEKMAT(L,NV,11) - FDNL1 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
             HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK33(Q) = ( HELP1 + HELP2 ) / FACT
           ELSE
             L_GK22(Q) = ZERO
             L_GK33(Q) = ZERO
           ENDIF

!  USUAL CASE WITH NO TRUNCATION

          ELSE

!  (1,1) AND (1,2) COMPONENTS

           L_GK11(Q) = L_GREEKMAT(Q,L,NV,1) * GREEKMAT(L,NV,1)

!mick fix
!           IF ( NSTOKES.GT.1 ) THEN
!             L_GK12(Q) = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
!           ENDIF
           IF ( NSTOKES.GT.1 ) THEN
             L_GK12(Q) = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
           ELSE
             L_GK12(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (4,4) AND (3,4)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
!             L_GK44(Q) = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
!             L_GK34(Q) = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
             L_GK44(Q) = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
             L_GK34(Q) = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
           ELSE
             L_GK44(Q) = ZERO
             L_GK34(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (2,2) AND (3,3)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT  ) THEN
!             L_GK22(Q) = L_GREEKMAT(Q,L,NV,6)  * GREEKMAT(L,NV,6)
!             L_GK33(Q) = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT  ) THEN
             L_GK22(Q) = L_GREEKMAT(Q,L,NV,6)  * GREEKMAT(L,NV,6)
             L_GK33(Q) = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
           ELSE
             L_GK22(Q) = ZERO
             L_GK33(Q) = ZERO
           ENDIF

          ENDIF
         ENDIF
        ENDDO

        IF ( L .EQ. 0 ) THEN

!  ADDING PAPER EQS. (76) AND (77) WITH M=0
!   ADDITIONAL FUNCTIONS P2M2 AND P2P2 ZERO FOR M = 0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = ONE
            P00(V,LNEW) = ZERO
            P02(V,LOLD) = ZERO
            P02(V,LNEW) = ZERO
            P2P2(V,LOLD) = ZERO
            P2P2(V,LNEW) = ZERO
            P2M2(V,LOLD) = ZERO
            P2M2(V,LNEW) = ZERO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = DL1/DL

! ADDING PAPER EQ. (81) WITH M=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = FAC1*UUU(V)*P00(V,LNEW) - FAC2*P00(V,LOLD)
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = QROOT6*(ONE-UUU(V)*UUU(V))
            P02(V,LNEW) = ZERO
          END DO
          SQL41 = ZERO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          DO V = 1, N_GEOMETRIES
            P2P2(V,LOLD)= 0.25D0*(ONE+UUU(V))*(ONE+UUU(V))
            P2M2(V,LOLD)= 0.25D0*(ONE-UUU(V))*(ONE-UUU(V))
          ENDDO

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = TMP1*UUU(V)*P02(V,LNEW) - TMP2*P02(V,LOLD)
          END DO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          FL2    = TWO * DL - ONE
          FLL1   = DL * DL1
          PERLL4 = ONE/(DL1*SQL41**2)
          QFAC   = DL  * ( DL1*DL1 - FOUR)
          DO V = 1, N_GEOMETRIES
           WFACT = FL2 * ( FLL1 * UUU(V) - FOUR )
           P2P2(V,LOLD) = (WFACT*P2P2(V,LNEW) - &
                QFAC*P2P2(V,LOLD)) * PERLL4
           WFACT = FL2 * ( FLL1 * UUU(V) + FOUR )
           P2M2(V,LOLD) = (WFACT*P2M2(V,LNEW) - &
                QFAC*P2M2(V,LOLD)) * PERLL4
          ENDDO

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

! SECTION FOR RANDOMLY-ORIENTED SPHEROIDS, ADDED 20 MARCH 2006
!  R. SPURR AND V. NATRAJ

        DO Q = 1, NV_PARAMETERS

         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( L.LE.LAYER_MAXMOMENTS(NV) ) THEN
           DO V = 1, N_GEOMETRIES
            L_FMAT(Q,V,INDEX_11) = &
                  L_FMAT(Q,V,INDEX_11) + L_GK11(Q) * P00(V,LNEW)
            L_FMAT(Q,V,INDEX_12) = &
                  L_FMAT(Q,V,INDEX_12) + L_GK12(Q) * P02(V,LNEW)

!  THIS SECTION ADDS F-MATRIX TERMS NEEDED FOR RANDOMLY-ORIENTED SPHEROIDS
            SUM23 = L_GK22(Q) + L_GK33(Q)
            DIF23 = L_GK22(Q) - L_GK33(Q)
            L_FMAT(Q,V,INDEX_22) = &
                  L_FMAT(Q,V,INDEX_22) + SUM23 * P2P2(V,LNEW)
            L_FMAT(Q,V,INDEX_33) = &
                  L_FMAT(Q,V,INDEX_33) + DIF23 * P2M2(V,LNEW)
!  END OF NEW SECTION---------------------------------------------------

            L_FMAT(Q,V,INDEX_44) = &
                  L_FMAT(Q,V,INDEX_44) + L_GK44(Q) * P00(V,LNEW)
            L_FMAT(Q,V,INDEX_34) = &
                  L_FMAT(Q,V,INDEX_34) + L_GK34(Q) * P02(V,LNEW)

           ENDDO
          ENDIF
         ENDIF
        END DO

!  END MOMENT LOOP

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP.

      DO V = 1, N_GEOMETRIES
        DO Q = 1, NV_PARAMETERS
          L_FMAT(Q,V,INDEX_22) = HALF * &
             ( L_FMAT(Q,V,INDEX_22) + L_FMAT(Q,V,INDEX_33) )
          L_FMAT(Q,V,INDEX_33) = &
             ( L_FMAT(Q,V,INDEX_22) - L_FMAT(Q,V,INDEX_33) )
        END DO
      END DO

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE.
!   R. SPURR AND V. NATRAJ, 20 MARCH 2006

!      DO V = 1, N_GEOMETRIES
!        DO Q = 1, NV_PARAMETERS
!          L_FMAT(Q,V,INDEX_22) = L_FMAT(Q,V,INDEX_11)
!          L_FMAT(Q,V,INDEX_33) = L_FMAT(Q,V,INDEX_44)
!        END DO
!      END DO

!  CREATE LINEARIZED Z MATRICES
!  ----------------------------

!mick fix 1/16/2013 - initialize matrix entries for the current NV
        DO Q = 1, NV_PARAMETERS
          DO V = 1, N_GEOMETRIES
            L_ZMAT(Q,V,NV,:,:) = ZERO
          ENDDO
        ENDDO

!  COPY FOR STOKES = 1

      IF ( NSTOKES .EQ. 1 ) THEN
        DO Q = 1, NV_PARAMETERS
          DO V = 1, N_GEOMETRIES
            L_ZMAT(Q,V,NV,1,1) = L_FMAT(Q,V,INDEX_11)
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  FOR POLARIZED CASE, SUNLIGHT ONLY
!    ** ONLY NEED THE FIRST COLUMN OF THE Z-MATRIX

      IF ( DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            L_ZMAT(Q,V,NV,1,1) =   L_FMAT(Q,V,INDEX_11)
            L_ZMAT(Q,V,NV,2,1) =  -L_FMAT(Q,V,INDEX_12) * C2(V)
            L_ZMAT(Q,V,NV,3,1) =   L_FMAT(Q,V,INDEX_12) * S2(V)
            L_ZMAT(Q,V,NV,4,1) =   ZERO
          ENDDO
        ENDDO
      ENDIF

!  FOR POLARIZED CASE, NON-SUNLIGHT GENERAL CASE

      IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            L_ZMAT(Q,V,NV,1,1) =   L_FMAT(Q,V,INDEX_11)
            L_ZMAT(Q,V,NV,2,1) =  -L_FMAT(Q,V,INDEX_12) * C2(V)
            L_ZMAT(Q,V,NV,3,1) =   L_FMAT(Q,V,INDEX_12) * S2(V)
            L_ZMAT(Q,V,NV,4,1) =   ZERO
!  THIS CODE UNTESTED
            HELP2C1 = L_FMAT(Q,V,INDEX_22) * C1(V)
            HELP2S1 = L_FMAT(Q,V,INDEX_22) * S1(V)
            HELP3C1 = L_FMAT(Q,V,INDEX_33) * C1(V)
            HELP3S1 = L_FMAT(Q,V,INDEX_33) * S1(V)
            L_ZMAT(Q,V,NV,1,2) =   L_FMAT(Q,V,INDEX_12) * C1(V)
            L_ZMAT(Q,V,NV,1,3) = - L_FMAT(Q,V,INDEX_12) * S1(V)
            L_ZMAT(Q,V,NV,1,4) =   ZERO
            L_ZMAT(Q,V,NV,2,2) =   C2(V) * HELP2C1 - S2(V) * HELP3S1
            L_ZMAT(Q,V,NV,2,3) = - C2(V) * HELP2S1 - S2(V) * HELP3C1
            L_ZMAT(Q,V,NV,2,4) = - L_FMAT(Q,V,INDEX_34) * S2(V)
            L_ZMAT(Q,V,NV,3,2) =   S2(V) * HELP2C1 + C2(V) * HELP3S1
            L_ZMAT(Q,V,NV,3,3) = - S2(V) * HELP2S1 + C2(V) * HELP3C1
            L_ZMAT(Q,V,NV,3,4) =   L_FMAT(Q,V,INDEX_34) * C2(V)
            L_ZMAT(Q,V,NV,4,2) = - L_FMAT(Q,V,INDEX_34) * S1(V)
            L_ZMAT(Q,V,NV,4,3) = - L_FMAT(Q,V,INDEX_34) * C1(V)
            L_ZMAT(Q,V,NV,4,4) =   L_FMAT(Q,V,INDEX_44)
          ENDDO
        ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE VLIDORTSS_L_ZMATRICES

!

      SUBROUTINE SSCORR_OUTGOING_L_ZMATRIX ( &
        DO_SSCORR_TRUNCATION, DO_SUNLIGHT, DO_LINEAR, &
        LAYER, NSTOKES, NGREEKMOMS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        LAYER_MAXMOMENTS, GREEKMATRIX, SSFDEL, &
        DO_SCATMAT, L_GREEKMATRIX, L_SSFDEL, &
        CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
        PHI, COSSCAT, VSIGN, &
        ZMAT, L_ZMAT, FMAT )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  -----

!  LINEARIZATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_LINEAR
      LOGICAL, INTENT (IN) ::          DO_SCATMAT ( MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER

!  CONTROL INTEGERS AND SUNLIGHT FLAG

      LOGICAL, INTENT (IN) ::          DO_SUNLIGHT
      INTEGER, INTENT (IN) ::          LAYER, NSTOKES, NGREEKMOMS

!  SCATTERING INPUT INFORMATION

      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS
      DOUBLE PRECISION, INTENT (IN) :: GREEKMATRIX &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMATRIX &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )

!  TRUNCATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION, INTENT (IN) :: SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SSFDEL ( MAXLAYERS,MAX_ATMOSWFS )

!  ZENITH ANGLE, AZIMUTH ANGLE, SCATTER ANGLE (SINES/COSINES)

      DOUBLE PRECISION, INTENT (IN) :: CTHETA
      DOUBLE PRECISION, INTENT (IN) :: STHETA
      DOUBLE PRECISION, INTENT (IN) :: CALPHA
      DOUBLE PRECISION, INTENT (IN) :: SALPHA
      DOUBLE PRECISION, INTENT (IN) :: CPHI
      DOUBLE PRECISION, INTENT (IN) :: COSSCAT
      DOUBLE PRECISION, INTENT (IN) :: VSIGN

!  AZIMUTH ANGLE. ADDED FOR VERSION 2.4R, REQUIRED FOR PHI > 180.
!    R. SPURR AND V. NATRAJ, 01 MAY 2009

      DOUBLE PRECISION, INTENT (IN) :: PHI

!  OUTPUT
!  ------

!  Z-MATRIX, FIRST COLUMN ONLY. F-MATRIX ( SUNLIGHT CASE )
!  GENERALIZED, 05 OCTOBER 2010, NON-TESTED

      DOUBLE PRECISION, INTENT (OUT) :: ZMAT ( 4, 4 )
      DOUBLE PRECISION, INTENT (OUT) :: L_ZMAT ( 4, 4, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: FMAT ( 6 )

!  LOCAL
!  -----

      DOUBLE PRECISION :: C1, S1, C2, S2, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION :: DL, QROOT6, UUU, PHIDEG, DNL1, FDNL1, FACT
      DOUBLE PRECISION :: FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION :: TMP1, TMP2, SUM23, DIF23, HELP1, HELP2
      DOUBLE PRECISION :: GK11, GK12, GK34, GK44, GK22, GK33
      DOUBLE PRECISION :: L_GK11(MAX_ATMOSWFS),L_GK12(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_GK44(MAX_ATMOSWFS),L_GK34(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_GK22(MAX_ATMOSWFS),L_GK33(MAX_ATMOSWFS)
      DOUBLE PRECISION :: FL2, FLL1, PERLL4, QFAC, WFACT, DL1
      DOUBLE PRECISION :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

      DOUBLE PRECISION :: L_FMAT ( MAX_ATMOSWFS, 6 )

      INTEGER ::          GREEKMAT_INDEX(6), K, L, LNEW, LOLD, ITMP, N, Q
      INTEGER ::          INDEX_11, INDEX_12, INDEX_34
      INTEGER ::          INDEX_22, INDEX_33, INDEX_44

!  GENERALIZED SPHERICAL FUNCTIONS
!    P2P2 AND P2M2 ADDED, 05 OCTOBER 2010

      DOUBLE PRECISION :: P00(2), P02(2), P2P2(2), P2M2(2)

!  INDEXING KEY

!      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
!      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
!      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
!      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
!      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
!      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      N = LAYER

!      IF ( V.EQ.42.AND.N.EQ.1.AND.VSIGN.GT.0.0D0 ) THEN
!         WRITE(*,*)'JJJ'
!         WRITE(*,*)CTHETA, STHETA, CALPHA, SALPHA, CPHI,
!     I    PHI, COSSCAT, VSIGN
!      ENDIF

!  GEOMETRICAL QUANTITIES
!  ----------------------

!  COSINE SCATTER ANGLE (THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  SHOULD BE THE SAME AS THE INPUT VALUE CSA

      UUU  = COSSCAT

!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

      HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
      IF ( HELP_SINSCAT.LE.ZERO ) THEN
        SINSCAT = 1.0D-12
      ELSE
        SINSCAT = DSQRT ( HELP_SINSCAT )
      ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

      IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
        CSIG1 = ZERO
        CSIG2 = ZERO
      ELSE
        IF ( STHETA .EQ. ZERO ) THEN
          CSIG1 = -  CPHI
        ELSE
          CSIG1 = (-VSIGN*CALPHA+CTHETA*COSSCAT)/SINSCAT/STHETA
        ENDIF
        IF ( SALPHA .EQ. ZERO ) THEN
          CSIG2 = -  CPHI
        ELSE
          CSIG2 = (-CTHETA+VSIGN*CALPHA*COSSCAT)/SINSCAT/SALPHA
        ENDIF
      ENDIF

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

      IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
      IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

      IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
      IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

      CSIG1_2 = TWO * CSIG1
      CSIG2_2 = TWO * CSIG2
      IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
        SSIG1 = ZERO
      ELSE
        SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
      ENDIF
      IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
        SSIG2 = ZERO
      ELSE
        SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
      ENDIF

!  FOR RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1 AND S2
!  SEE H/VDM, EQS. 94-95

      C1 = CSIG1_2 * CSIG1 - ONE
      C2 = CSIG2_2 * CSIG2 - ONE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  PHI IS INPUT IN RADIANS HERE. MUST CONVERT TO DEGREES
!     WE WERE USING DEGREES FOR THE S1/S2 SIGN REVERSLA

      PHIDEG = PHI / DEG_TO_RAD
!      IF ( PHI    .LE. 180.0D0 ) THEN            ! OLD
      IF ( PHIDEG .LE. 180.0D0 ) THEN             ! NEW
        S1 = CSIG1_2 * SSIG1
        S2 = CSIG2_2 * SSIG2
      ELSE
        S1 = -CSIG1_2 * SSIG1
        S2 = -CSIG2_2 * SSIG2
      ENDIF

 !  B G CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  F-MATRICES AND LINEARIZED F-MATRICES
!  ------------------------------------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! INITIALISE F-MATRIX AND ITS LINEARIZATION

      DO K = 1, 6
        FMAT(K) = 0.0D0
      END DO
      IF ( DO_LINEAR ) THEN
       IF ( LAYER_VARY_FLAG ) THEN
        DO Q = 1, LAYER_VARY_NUMBER
         DO K = 1, 6
          L_FMAT(Q,K) = 0.0D0
         ENDDO
        ENDDO
       ENDIF
      ENDIF

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  SET THE LOCAL GREEK MATRIX ELEMENTS THAT YOU NEED
!   44 AND 34 ARE NOT REQUIRED WITH NATURAL SUNLIGHT (DEFAULT HERE)
!   22 AND 33 REQUIRED FOR NON-MIE SPHEROIDAL PARTICLES

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DNL1  = DBLE(2*L + 1 )
          FDNL1 = SSFDEL(N) * DNL1
          FACT  = ONE - SSFDEL(N)

!  REGULAR COEFFICIENTS

          GK11 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(1)) - FDNL1 ) / FACT
          GK12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
           GK34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
           GK22 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(2)) - FDNL1 ) / FACT
           GK33 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(4)) - FDNL1 ) / FACT
          ENDIF

!  LINEARIZED COEFFICIENTS

          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN

!  (1,1) AND (1,2) CASES

               HELP1 = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(N,Q)
               L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
               HELP1 = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
               HELP2 = GK12 * L_SSFDEL(N,Q)
               L_GK12(Q) = ( HELP1 + HELP2 ) / FACT

!  OTHER NON-SUNLIGHT CASES

               IF ( .NOT. DO_SUNLIGHT ) THEN

                 HELP1 = L_GREEKMATRIX(Q,L,N,16) * GREEKMATRIX(L,N,16)
                 HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
                 HELP1 = L_GREEKMATRIX(Q,L,N,12) * GREEKMATRIX(L,N,12)
                 HELP2 = GK34 * L_SSFDEL(N,Q)
                 L_GK34(Q) = ( HELP1 + HELP2 ) / FACT

                 HELP1 = L_GREEKMATRIX(Q,L,N,6)  * GREEKMATRIX(L,N,6)
                 HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
                 HELP1 = L_GREEKMATRIX(Q,L,N,11) * GREEKMATRIX(L,N,11)
                 HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK33(Q) = ( HELP1 + HELP2 ) / FACT

               ENDIF

!  NO SCAT-MAT VARIATION

             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = ZERO
                L_GK34(Q) = ZERO
                L_GK22(Q) = ZERO
                L_GK33(Q) = ZERO
               ENDIF
             ENDIF

!  END LINEARIZATION CONTROL

            ENDDO
           ENDIF
          ENDIF

!  NO SSCORR TRUNCATION, JUST COPY

        ELSE

!  REGULAR COEFFICIENTS

          GK11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
          GK12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
           GK34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
           GK22 = GREEKMATRIX(L,N,GREEKMAT_INDEX(2))
           GK33 = GREEKMATRIX(L,N,GREEKMAT_INDEX(4))
          ENDIF

! LINAERIZED COEFFICIENTS

          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN
               L_GK11(Q) = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               L_GK12(Q) = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = L_GREEKMATRIX(Q,L,N,16)*GREEKMATRIX(L,N,16)
                L_GK34(Q) = L_GREEKMATRIX(Q,L,N,12)*GREEKMATRIX(L,N,12)
                L_GK22(Q) = L_GREEKMATRIX(Q,L,N,6) *GREEKMATRIX(L,N,6)
                L_GK33(Q) = L_GREEKMATRIX(Q,L,N,11)*GREEKMATRIX(L,N,11)
               ENDIF
             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = ZERO
                L_GK34(Q) = ZERO
                L_GK22(Q) = ZERO
                L_GK33(Q) = ZERO
               ENDIF
             ENDIF
            ENDDO
           ENDIF
          ENDIF
        ENDIF

!  FIRST MOMENT

        IF ( L .EQ. 0 ) THEN

!  ADDING PAPER EQS. (76) AND (77) WITH M=0
!   ADDITIONAL FUNCTIONS P2M2 AND P2P2 ZERO FOR M = 0

          P00(LOLD) = ONE
          P00(LNEW) = ZERO
          P02(LOLD) = ZERO
          P02(LNEW) = ZERO
          P2P2(LOLD) = ZERO
          P2P2(LNEW) = ZERO
          P2M2(LOLD) = ZERO
          P2M2(LNEW) = ZERO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = DL1/DL

! ADDING PAPER EQ. (81) WITH M=0

          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)

        END IF

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          P02(LOLD) = QROOT6*(ONE-UUU*UUU)
          P02(LNEW) = ZERO
          SQL41 = ZERO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          P2P2(LOLD)= 0.25D0*(ONE+UUU)*(ONE+UUU)
          P2M2(LOLD)= 0.25D0*(ONE-UUU)*(ONE-UUU)

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-4.0D0)
          TMP1  = (2.0D0*DL-1.0D0)/SQL41
          TMP2  = SQL4/SQL41
          P02(LOLD) = TMP1*UUU*P02(LNEW) - TMP2*P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          FL2 = TWO * DL - ONE
          FLL1 = DL * DL1
          PERLL4=ONE/(DL1*SQL41**2)
          QFAC  = DL  * ( DL1*DL1 - FOUR)
          WFACT = FL2 * ( FLL1 * UUU - FOUR )
          P2P2(LOLD) = (WFACT*P2P2(LNEW) - QFAC*P2P2(LOLD)) * PERLL4
          WFACT = FL2 * ( FLL1 * UUU + FOUR )
          P2M2(LOLD) = (WFACT*P2M2(LNEW) - QFAC*P2M2(LOLD)) * PERLL4

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  DO NOT NEED THESE WITH NATURAL SUNLIGHT:
!            FMAT(INDEX_34)AND FMAT(INDEX_44)

! SECTION FOR RANDOMLY-ORIENTED SPHEROIDS, ADDED 05 OCTOBER 2010
!  R. SPURR AND V. NATRAJ

        IF ( L.LE.LAYER_MAXMOMENTS ) THEN

          FMAT(INDEX_11) = FMAT(INDEX_11) + GK11 * P00(LNEW)
          FMAT(INDEX_12) = FMAT(INDEX_12) + GK12 * P02(LNEW)

          IF ( .NOT. DO_SUNLIGHT ) THEN
            SUM23 = GK22 + GK33
            DIF23 = GK22 - GK33
            FMAT(INDEX_22) = FMAT(INDEX_22) + SUM23 * P2P2(LNEW)
            FMAT(INDEX_33) = FMAT(INDEX_33) + DIF23 * P2M2(LNEW)
            FMAT(INDEX_44) = FMAT(INDEX_44) + GK44 * P00(LNEW)
            FMAT(INDEX_34) = FMAT(INDEX_34) + GK34 * P02(LNEW)
          ENDIF

        ENDIF

!  LINEARIZATION

        IF ( DO_LINEAR ) THEN
         IF ( LAYER_VARY_FLAG ) THEN
          DO Q = 1, LAYER_VARY_NUMBER
           IF ( DO_SCATMAT(N,Q) ) THEN
            IF ( L.LE.LAYER_MAXMOMENTS ) THEN

             L_FMAT(Q,INDEX_11) = L_FMAT(Q,INDEX_11)+L_GK11(Q)*P00(LNEW)
             L_FMAT(Q,INDEX_12) = L_FMAT(Q,INDEX_12)+L_GK12(Q)*P02(LNEW)
             IF ( .NOT. DO_SUNLIGHT ) THEN
              SUM23 = L_GK22(Q) + L_GK33(Q)
              DIF23 = L_GK22(Q) - L_GK33(Q)
              L_FMAT(Q,INDEX_22) = L_FMAT(Q,INDEX_22)+SUM23*P2P2(LNEW)
              L_FMAT(Q,INDEX_33) = L_FMAT(Q,INDEX_33)+DIF23*P2M2(LNEW)
              L_FMAT(Q,INDEX_44) =L_FMAT(Q,INDEX_44)+L_GK44(Q)*P00(LNEW)
              L_FMAT(Q,INDEX_34) =L_FMAT(Q,INDEX_34)+L_GK34(Q)*P02(LNEW)
             ENDIF

            ENDIF
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  END MOMENTS

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP. (NON-SUNLIGHT ONLY)

      IF ( .NOT. DO_SUNLIGHT ) THEN

        FMAT(INDEX_22) = HALF * ( FMAT(INDEX_22) + FMAT(INDEX_33) )
        FMAT(INDEX_33) = FMAT(INDEX_22) - FMAT(INDEX_33)

        IF ( DO_LINEAR ) THEN
          IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
              IF ( DO_SCATMAT(N,Q) ) THEN
                L_FMAT(Q,INDEX_22) = HALF * &
                   ( L_FMAT(Q,INDEX_22) + L_FMAT(Q,INDEX_33) )
                L_FMAT(Q,INDEX_33) = &
                   ( L_FMAT(Q,INDEX_22) - L_FMAT(Q,INDEX_33) )
              ENDIF
            ENDDO
          ENDIF
        ENDIF

      ENDIF

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE. R. SPURR AND V. NATRAJ,, 20 MARCH 2006

!      IF ( .NOT. DO_SUNLIGHT ) THEN
!        FMAT(INDEX_22) = FMAT(INDEX_11)
!        FMAT(INDEX_33) = FMAT(INDEX_44)
!        IF ( DO_LINEAR ) THEN
!          IF ( LAYER_VARY_FLAG ) THEN
!            DO Q = 1, LAYER_VARY_NUMBER
!              IF ( DO_SCATMAT(N,Q) ) THEN
!                L_FMAT(Q,INDEX_22) = L_FMAT(Q,INDEX_11)
!                L_FMAT(Q,INDEX_33) = L_FMAT(Q,INDEX_44)
!              ENDIF
!            ENDDO
!          ENDIF
!        ENDIF
!      ENDIF

!  Z-MATRIX CALCULATION
!  ====================

!  INITIALIZE Z-MATRIX + LINEARIZATION

      ZMAT   = ZERO
      L_ZMAT = ZERO

!  SCALAR CASE

      IF ( NSTOKES .EQ. 1 ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
      ENDIF

!  NATURAL LIGHT, FIRST COLUMN ONLY

      IF ( DO_SUNLIGHT .AND. NSTOKES.GT.1 ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
        ZMAT(2,1) =  -FMAT(INDEX_12) * C2
        ZMAT(3,1) =   FMAT(INDEX_12) * S2
        ZMAT(4,1) =   ZERO
      ENDIF

!  OTHER SOURCES, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.GT.1 ) THEN
         HELP2C1 = FMAT(INDEX_22) * C1
         HELP2S1 = FMAT(INDEX_22) * S1
         HELP3C1 = FMAT(INDEX_33) * C1
         HELP3S1 = FMAT(INDEX_33) * S1
         ZMAT(1,2) =   FMAT(INDEX_12) * C1
         ZMAT(1,3) = - FMAT(INDEX_12) * S1
         ZMAT(1,4) =   ZERO
         ZMAT(2,2) =   C2 * HELP2C1 - S2 * HELP3S1
         ZMAT(2,3) = - C2 * HELP2S1 - S2 * HELP3C1
         ZMAT(2,4) = - FMAT(INDEX_34) * S2
         ZMAT(3,2) =   S2 * HELP2C1 + C2 * HELP3S1
         ZMAT(3,3) = - S2 * HELP2S1 + C2 * HELP3C1
         ZMAT(3,4) =   FMAT(INDEX_34) * C2
         ZMAT(4,2) = - FMAT(INDEX_34) * S1
         ZMAT(4,3) = - FMAT(INDEX_34) * C1
         ZMAT(4,4) =   FMAT(INDEX_44)
      ENDIF

!  RETURN IF NO LINEARIZATION

      IF ( .NOT. DO_LINEAR ) RETURN

!  LINEARIZED Z-MATRIX
!  -------------------

!  COPY FOR STOKES = 1
!  FOR POLARIZED CASE, SUNLIGHT ONLY !!!!!!!
!   FOR SUNLIGHT, ONLY NEED THE FIRST COLUMN OF THE Z-MATRIX

      DO Q = 1, LAYER_VARY_NUMBER

        IF ( NSTOKES .EQ. 1 ) THEN
          L_ZMAT(1,1,Q) = L_FMAT(Q,INDEX_11)
        ENDIF

        IF ( DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
          L_ZMAT(1,1,Q) =   L_FMAT(Q,INDEX_11)
          L_ZMAT(2,1,Q) =  -L_FMAT(Q,INDEX_12) * C2
          L_ZMAT(3,1,Q) =   L_FMAT(Q,INDEX_12) * S2
          L_ZMAT(4,1,Q) =   ZERO
        ENDIF

!  OTHER SOURCES, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

        IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
           HELP2C1 = L_FMAT(Q,INDEX_22) * C1
           HELP2S1 = L_FMAT(Q,INDEX_22) * S1
           HELP3C1 = L_FMAT(Q,INDEX_33) * C1
           HELP3S1 = L_FMAT(Q,INDEX_33) * S1
           L_ZMAT(1,2,Q) =   L_FMAT(Q,INDEX_12) * C1
           L_ZMAT(1,3,Q) = - L_FMAT(Q,INDEX_12) * S1
           L_ZMAT(1,4,Q) =   ZERO
           L_ZMAT(2,2,Q) =   C2 * HELP2C1 - S2 * HELP3S1
           L_ZMAT(2,3,Q) = - C2 * HELP2S1 - S2 * HELP3C1
           L_ZMAT(2,4,Q) = - L_FMAT(Q,INDEX_34) * S2
           L_ZMAT(3,2,Q) =   S2 * HELP2C1 + C2 * HELP3S1
           L_ZMAT(3,3,Q) = - S2 * HELP2S1 + C2 * HELP3C1
           L_ZMAT(3,4,Q) =   L_FMAT(Q,INDEX_34) * C2
           L_ZMAT(4,2,Q) = - L_FMAT(Q,INDEX_34) * S1
           L_ZMAT(4,3,Q) = - L_FMAT(Q,INDEX_34) * C1
           L_ZMAT(4,4,Q) =   L_FMAT(Q,INDEX_44)
        ENDIF

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE SSCORR_OUTGOING_L_ZMATRIX

      END MODULE vlidort_la_corrections

