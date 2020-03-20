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
! #                                                             #
! #              VLIDORT_LP_QBEAM_SOLUTION                      #
! #              VLIDORT_LP_UBEAM_SOLUTION                      #
! #                                                             #
! ###############################################################


      MODULE vlidort_lp_solutions

      PRIVATE
      PUBLIC :: VLIDORT_LP_QBEAM_SOLUTION, &
                VLIDORT_LP_UBEAM_SOLUTION

      CONTAINS

      SUBROUTINE VLIDORT_LP_QBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, &
        DOVARY, N_PARAMETERS, &
        DO_PLANE_PARALLEL, NSTOKES, &
        NSTREAMS, FLUX_FACTOR, &
        LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        NMOMENTS, NSTREAMS_2, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        DMAT, DFLUX, &
        AVERAGE_SECANT, &
        PI_XQP, PI_X0P, PI_XQM_POST, &
        SAB, DAB, &
        QSUMVEC_SAVE, QDIFVEC_SAVE, &
        QVEC_SAVE, QDIF_SAVE, &
        QMAT_SAVE, QPIVOT, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        L_SAB, L_DAB, L_EIGENMAT, &
        L_OMEGA_GREEK, LP_AVERAGE_SECANT, &
        LP_BVEC, &
        STATUS, MESSAGE, TRACE )

!  linearized values of the classical particular solution.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

!  regular values (strictly in)

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM
      LOGICAL, INTENT (IN) ::          DOVARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Regular values (strictly in)

      DOUBLE PRECISION, INTENT (IN) :: SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (IN) :: QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (IN) :: QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER, INTENT (IN) ::          QPIVOT ( MAXSTRMSTKS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Linearized  values (strictly in)

      DOUBLE PRECISION, INTENT (IN) :: L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EIGENMAT &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly out)

      DOUBLE PRECISION, INTENT (INOUT) ::  LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Exception

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Saved vector

      DOUBLE PRECISION ::  QSAVE(MAXSTRMSTKS)

!  linearization arrays

      DOUBLE PRECISION ::  L_QVEC(MAXSTRMSTKS,MAX_ATMOSWFS)

      DOUBLE PRECISION ::  L_QSUMVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION ::  L_QDIFVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION ::  HELP_TEMPO(0:MAXMOMENTS,MAXSTOKES)

!  help variables

      DOUBLE PRECISION ::  HELP, WSUM, WDIF, SECBAR
      DOUBLE PRECISION ::  TM1, TM2, TM3, L_WVEC_D, L_WVEC_U
      DOUBLE PRECISION ::  SUM, GSUM, S_TPA, S_TMA, TPA, TMA, F1

      INTEGER          ::  I, J, I1, L, N, M, INFO, IB, LVARY
      INTEGER          ::  O1, O2, O3, IR, JC, IROW, JCOL
      INTEGER          ::  Q, K, K_PARAMETERS
      LOGICAL          ::  DO_FIRST
      CHARACTER (LEN=3) :: CN, CI

!  Start of code
!  -------------

!  initialise Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Set layer, fourier, beam index

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  safety

      SECBAR = ZERO

!  set up driving vector (using saved results)
!  This must be done regardless of whether layer N is varying or not.

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          QSAVE(IROW) = QDIFVEC_SAVE(I,O1) + &
                       TWO * AVERAGE_SECANT(N,IB) * QVEC_SAVE(IROW)
        ENDDO
      ENDDO

!  Linearization for layer N is in two parts:
!    1A. Linearization due to variations in Layer N itself (Profiles)
!    1B. Linearization due to columns
!    2. Linearization due to variations in layers K < N (Profiles)

!  Part 1.
!  =======

!  Set layer to vary

      LVARY = N

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. &
                      DO_LAYER_SCATTERING(M,N)

!  No particular solution beyond the cutoff layer.
!  Or no scattering in this layer
!    [ Zero the boundary layer values and start Part 2 )

      IF (.NOT. DOVARY .OR. .NOT. DO_FIRST  ) THEN
        DO I = 1, NSTREAMS_2
         DO O1 = 1, NSTOKES
          DO Q = 1, N_PARAMETERS
            LP_BVEC(I,O1,N,LVARY,Q) = ZERO
          ENDDO
         ENDDO
        ENDDO
        GO TO 2222
      ENDIF

!  solar zenith cosine for this layer

      SECBAR = AVERAGE_SECANT(N,IB)

!  For each varying parameter

      DO Q = 1, N_PARAMETERS

!  Auxiliary matrix for Q functions

        DO L = M, NMOMENTS
          DO O1 = 1, NSTOKES
            SUM = ZERO
            DO O2 = 1, NSTOKES
              GSUM = ZERO
              DO O3 = 1, NSTOKES
                GSUM = GSUM + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
              ENDDO
              SUM = SUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * GSUM
            ENDDO
            HELP_TEMPO(L,O1) = SUM
          ENDDO
        ENDDO

!  Set up linearized sum and difference vectors for Beam source terms

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            S_TPA = ZERO
            S_TMA = ZERO
            DO L = M, NMOMENTS
              TPA = ZERO
              TMA = ZERO
              DO O2 = 1, NSTOKES
                TPA = TPA + PI_XQP(L,I,O1,O2)      * HELP_TEMPO(L,O2)
                TMA = TMA + PI_XQM_POST(L,I,O1,O2) * HELP_TEMPO(L,O2)
              ENDDO
              S_TPA = S_TPA + TPA
              S_TMA = S_TMA + TMA
            ENDDO
            L_QSUMVEC(I,O1,Q) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
            L_QDIFVEC(I,O1,Q) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
          ENDDO
        ENDDO

!   setup linearized RHS vector
!  ( use results from the original solution )
!    WARNING. BE CAREFUL OF SIGNS on TM1, TM2, TM3
!     slgihtly different from the scalar model case.

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
           JC = NSTOKES*(J-1)
           DO O2 = 1, NSTOKES
            JCOL = JC + O2
            TM1 = L_EIGENMAT(IROW,JCOL,N,Q) * QVEC_SAVE(JCOL)
            TM2 =   DAB(I,J,O1,O2,N)   * L_QSUMVEC(J,O2,Q)
            TM3 = L_DAB(I,J,O1,O2,N,Q) *   QSUMVEC_SAVE(J,O2)
            HELP = HELP - TM1 + TM2 + TM3
!            HELP = HELP - TM1 - TM2 - TM3
           ENDDO
          ENDDO
          L_QVEC(IROW,Q)  = HELP + L_QDIFVEC(I,O1,Q) * SECBAR
         ENDDO
        ENDDO

!  end parameter loop

      ENDDO

!  additional terms for the quasi-spherical case
!  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                HELP = QSAVE(IROW) * LP_AVERAGE_SECANT(N,LVARY,IB,Q)
                L_QVEC(IROW,Q) = L_QVEC(IROW,Q) + HELP
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS &
              ('N',NSTKS_NSTRMS,N_PARAMETERS,QMAT_SAVE, &
                MAXSTRMSTKS,QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

!  Exception handling 1

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CN, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call # 1 in VLIDORT_LP_QBEAM_SOLUTION,'// &
                  ' Beam P.I. linearization (N,N), layer # '//CN
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  assign solutions for the quasi-spherical case, N > 1
!    Note change of sign with HELP, this is because SAB as
!    defined in VLIDORT is the negative of SAB in the scalar model.

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          TM3 = - QDIF_SAVE(IROW) / SECBAR
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
             TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
             HELP = HELP + TM1 + TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           TM2 = ( HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           WDIF = LP_AVERAGE_SECANT(N,LVARY,IB,Q) * TM3 + TM2
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           LP_BVEC(I,O1,N,LVARY,Q)  = L_WVEC_D
           LP_BVEC(I1,O1,N,LVARY,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) LP_BVEC(I1,O1,N,LVARY,Q) = - L_WVEC_U
!  Use this debug code
!         if(q.eq.1.and.i.eq.4.and.n.eq.24) write(39,'(4i3,1p2e20.10)')
!     &           m,ibeam,i,O1,BVEC(I,O1,N),LP_BVEC(I,O1,N,LVARY,Q)
          ENDDO
         ENDDO
        ENDDO

!  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
             TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
             HELP = HELP - TM1 - TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           WDIF = ( - HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           LP_BVEC(I,O1,N,LVARY,Q)  = L_WVEC_D
           LP_BVEC(I1,O1,N,LVARY,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) LP_BVEC(I1,O1,N,LVARY,Q) = - L_WVEC_U
          ENDDO
         ENDDO
        ENDDO

      ENDIF


!  Old code: Toggle to LP_BVEC

!      DO Q = 1, N_PARAMETERS
!       DO I = 1, NSTREAMS
!        I1 = I + NSTREAMS
!        DO O1 = 1, NSTOKES
!         LP_BVEC(I,O1,N,LVARY,Q) = L_WVEC(I,O1,N,LVARY,Q)
!         HELP = ZERO
!         DO O2 = 1, NSTOKES
!          HELP = HELP + DMAT(O1,O2) * L_WVEC(I1,O2,N,LVARY,Q)
!         ENDDO
!         LP_BVEC(I1,O1,N,LVARY,Q) = HELP
!        ENDDO
!       ENDDO
!      ENDDO


!  Part 2.
!  =======

!  Continuation point

 2222 CONTINUE

!  Only for the pseudo-spherical case, profile weighting functions
!  Also not required for the column weighting functions

      IF ( DO_PLANE_PARALLEL )       RETURN

!  No particular solution beyond the cutoff layer.
!    [ Zero the boundary layer values and exit )

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
              DO Q = 1, K_PARAMETERS
                LP_BVEC(I,O1,N,K,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Loop over layers K above N

      DO K = 1, N - 1

!  If there is a varying layer

        IF ( LAYER_VARY_FLAG(K) ) THEN

!  number of varying parameters for this layer

          K_PARAMETERS = LAYER_VARY_NUMBER(K)

!  Set up vector

          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, K_PARAMETERS
             L_QVEC(IROW,Q) = QSAVE(IROW) * LP_AVERAGE_SECANT(N,K,IB,Q)
            ENDDO
           ENDDO
          ENDDO

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

          CALL DGETRS &
              ('N',NSTKS_NSTRMS,K_PARAMETERS,QMAT_SAVE, &
                MAXSTRMSTKS,QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

!  Exception handling 2

          IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGETRS call # 2 in VLIDORT_LP_QBEAM_SOLUTION,'// &
                      ' Beam P.I. linearization (N,K), layer # '//CN
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  assign linearized solutions for layer N due to variations in layer K
!    Note change of sign with HELP, this is because SAB as
!    defined in VLIDORT is the negative of SAB in the scalar model.

          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           I1 = I + NSTREAMS
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            TM3 = - QDIF_SAVE(IROW) / SECBAR
            DO Q = 1, K_PARAMETERS
             HELP = ZERO
             DO J = 1, NSTREAMS
              JC = NSTOKES*(J-1)
              DO O2 = 1, NSTOKES
               JCOL = JC + O2
               HELP = HELP + SAB(I,J,O1,O2,N) * L_QVEC(JCOL,Q)
              ENDDO
             ENDDO
             WSUM = L_QVEC(IROW,Q)
             WDIF = LP_AVERAGE_SECANT(N,K,IB,Q) * TM3 + (HELP/SECBAR)
             L_WVEC_D = HALF * ( WSUM + WDIF )
             L_WVEC_U = HALF * ( WSUM - WDIF )
             LP_BVEC(I,O1,N,K,Q)  = L_WVEC_D
             LP_BVEC(I1,O1,N,K,Q) = L_WVEC_U
             IF ( O1 .GT. 2 ) LP_BVEC(I1,O1,N,K,Q) = - L_WVEC_U
            ENDDO
           ENDDO
          ENDDO

!  Old code.Toggle to LP_BVEC

!          DO Q = 1, K_PARAMETERS
!           DO I = 1, NSTREAMS
!            I1 = I + NSTREAMS
!            DO O1 = 1, NSTOKES
!             LP_BVEC(I,O1,N,K,Q) = L_WVEC(I,O1,N,K,Q)
!             HELP = ZERO
!             DO O2 = 1, NSTOKES
!              HELP = HELP + DMAT(O1,O2) * L_WVEC(I1,O2,N,K,Q)
!             ENDDO
!             LP_BVEC(I1,O1,N,K,Q) = HELP
!            ENDDO
!           ENDDO
!          ENDDO

!  end K-layer loop

        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LP_QBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_LP_UBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, &
        DO_VARY, N_PARAMETERS, DO_PLANE_PARALLEL, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NSTREAMS, FLUX_FACTOR, &
        LAYER_PIS_CUTOFF, QUAD_HALFWTS, &
        NMOMENTS, N_USER_STREAMS, &
        DO_LAYER_SCATTERING, LOCAL_UM_START, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        DFLUX, OMEGA_GREEK, &
        PI_XQP, PI_XUP, PI_XUM, PI_X0P, &
        PI_XQM_PRE, HELPSTOKES_BEAM, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        LP_BVEC, L_OMEGA_GREEK, &
        L_UPAR_DN_1, L_UPAR_UP_1, &
        LP_UPAR_DN_2, LP_UPAR_UP_2 )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Regular input (strictly in)

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM
      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: HELPSTOKES_BEAM &
          ( 0:MAXMOMENTS, MAXSTOKES )

!  Linearized (strictly in)

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
         ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized (Ostensibly in)

      DOUBLE PRECISION, INTENT (INOUT) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized (Ostensibly OUTPUT)

      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LP_UPAR_DN_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LP_UPAR_UP_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER, DO_FIRST
      INTEGER ::          IB, UM, L, N, M, O1, O2, O3, J, J1
      INTEGER ::          Q, K, K_PARAMETERS, LVARY

      DOUBLE PRECISION :: SUM1, SUM2, T1, T2, H, HH, F1
      DOUBLE PRECISION :: SP, SM, SPS, SPT, SGW(MAXSTOKES)

      DOUBLE PRECISION :: L_HELP_Q1(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_HELP_Q2(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_BEAM(MAXSTOKES)

      SGW = (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Layer and Fourier

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      DO_LAYER = .FALSE.
      LVARY = 0
      F1 = FLUX_FACTOR / PI4

!  Part 1. Linearizations for quantities varying in Layer N
!  ========================================================

!  Profiles and Columns

      LVARY = N

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. &
                      DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO Q = 1, N_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                L_UPAR_UP_1(UM,O1,N,Q)        = ZERO
                LP_UPAR_UP_2(UM,O1,N,LVARY,Q) = ZERO
              ENDDO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                L_UPAR_DN_1(UM,O1,N,Q)        = ZERO
                LP_UPAR_DN_2(UM,O1,N,LVARY,Q) = ZERO
              ENDDO
            ENDDO
           ENDIF
          ENDDO
        ELSE
          DO Q = 1, N_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO O1 = 1, NSTOKES
              L_UPAR_UP_1(IB,O1,N,Q)        = ZERO
              LP_UPAR_UP_2(IB,O1,N,LVARY,Q) = ZERO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO O1 = 1, NSTOKES
              L_UPAR_DN_1(IB,O1,N,Q)        = ZERO
              LP_UPAR_DN_2(IB,O1,N,LVARY,Q) = ZERO
            ENDDO
           ENDIF
          ENDDO
        ENDIF
        GO TO 2222
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. &
                   STERM_LAYERMASK_DN(N) )

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  first function
!    (linearization of primary scattering of beam)

       IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SP = ZERO
            DO O3 = 1, NSTOKES
              SP = SP + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SPS = SPS + L_OMEGA_GREEK(L,N,O1,O2,Q) * SP
          ENDDO
          L_HELP_Q1(L,O1) = SPS
         ENDDO
        ENDDO
       ENDIF

!  For each moment, do inner sum over computational angles
!   (Linearization of Diffuse scattering of beam)

       IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = ZERO
           SM = ZERO
           H = QUAD_HALFWTS(J)
           DO O2 = 1, NSTOKES
            HH = H * SGW(O2)
            SP = SP + H  * PI_XQP    (L,J,O1,O2)*LP_BVEC(J,O2,N,LVARY,Q)
            SM = SM + HH * PI_XQM_PRE(L,J,O1,O2)*LP_BVEC(J1,O2,N,LVARY,Q)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          L_HELPSTOKES_BEAM(O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SPT = ZERO
          DO O2 = 1, NSTOKES
           SPT = SPT + OMEGA_GREEK(L,N,O1,O2)   * L_HELPSTOKES_BEAM(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) * HELPSTOKES_BEAM(L,O2)
          ENDDO
          L_HELP_Q2(L,O1) = SPT
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integrated contributions

       IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            T1 = ZERO
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM1 = ZERO
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUM(L,UM,O1,O2)
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,UM,O1,O2)
             ENDDO
             T1 = T1 + SUM1
             T2 = T2 + SUM2
            ENDDO
            L_UPAR_UP_1(UM,O1,N,Q)  = F1 * T1
            LP_UPAR_UP_2(UM,O1,N,LVARY,Q) = T2
           ENDDO
          ENDDO
        ELSE
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUM(L,IB,O1,O2)
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,IB,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          L_UPAR_UP_1(IB,O1,N,Q)  = F1 * T1
          LP_UPAR_UP_2(IB,O1,N,LVARY,Q) = T2
         ENDDO
        ENDIF
       ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integrated contributions

       IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            T1 = ZERO
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM1 = ZERO
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUP(L,UM,O1,O2)
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,UM,O1,O2)
             ENDDO
             T1 = T1 + SUM1
             T2 = T2 + SUM2
            ENDDO
            L_UPAR_DN_1(UM,O1,N,Q)  = F1 * T1
            LP_UPAR_DN_2(UM,O1,N,LVARY,Q) = T2
           ENDDO
          ENDDO
        ELSE
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUP(L,IB,O1,O2)
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,IB,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          L_UPAR_DN_1(IB,O1,N,Q)  = F1 * T1
          LP_UPAR_DN_2(IB,O1,N,LVARY,Q) = T2
         ENDDO
        ENDIF
       ENDIF

!  end parameter loop

      ENDDO

!  Part 2.
!  =======

!  Continuation point

 2222 CONTINUE

!  Only these extra variations when the following conditions
!  are not satisfied (pseudo-spherical, layer > 1)
!  Profiles only

      IF ( N .EQ. 1 )                RETURN
      IF ( DO_PLANE_PARALLEL )       RETURN

!  No particular solution beyond the cutoff layer.
!   Or no scattering in this layer
!    [ Zero the boundary layer values and exit )

      IF ( .NOT. DO_FIRST ) THEN
       DO K = 1, N - 1
        K_PARAMETERS = LAYER_VARY_NUMBER(K)

        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO Q = 1, K_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LP_UPAR_UP_2(UM,O1,N,K,Q) = ZERO
              ENDDO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LP_UPAR_DN_2(UM,O1,N,K,Q) = ZERO
              ENDDO
            ENDDO
           ENDIF
          ENDDO
        ELSE
          DO Q = 1, K_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO O1 = 1, NSTOKES
              LP_UPAR_UP_2(IB,O1,N,K,Q) = ZERO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO O1 = 1, NSTOKES
              LP_UPAR_DN_2(IB,O1,N,K,Q) = ZERO
            ENDDO
           ENDIF
          ENDDO
        ENDIF

       ENDDO
       RETURN
      ENDIF

!  start loop over all layers above N

      DO K = 1, N - 1

!  Start parameter loop -- only do if layer K has some variation

       IF ( LAYER_VARY_FLAG(K) ) THEN
        K_PARAMETERS = LAYER_VARY_NUMBER(K)
        DO Q = 1, K_PARAMETERS

!  For each moment, do inner sum over computational angles
!   (Linearization of Diffuse scattering of beam)

         IF ( DO_LAYER ) THEN
          DO L = M, NMOMENTS
           DO O1 = 1, NSTOKES
            SPS = ZERO
            DO J = 1, NSTREAMS
             J1 = J + NSTREAMS
             SP = ZERO
             SM = ZERO
             H = QUAD_HALFWTS(J)
             DO O2 = 1, NSTOKES
              HH = H * SGW(O2)
              SP = SP + H  * PI_XQP    (L,J,O1,O2)*LP_BVEC(J,O2,N,K,Q)
              SM = SM + HH * PI_XQM_PRE(L,J,O1,O2)*LP_BVEC(J1,O2,N,K,Q)
             ENDDO
             SPS = SPS + SP + SM
            ENDDO
            L_HELPSTOKES_BEAM(O1) = SPS
           ENDDO
           DO O1 = 1, NSTOKES
            SPT = ZERO
            DO O2 = 1, NSTOKES
             SPT = SPT + OMEGA_GREEK(L,N,O1,O2) * L_HELPSTOKES_BEAM(O2)
            ENDDO
            L_HELP_Q2(L,O1) = SPT
           ENDDO
          ENDDO
         ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Diffuse (integrated) contributions only

         IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO O1 = 1, NSTOKES
              T2 = ZERO
              DO L = M, NMOMENTS
               SUM2 = ZERO
               DO O2 = 1, NSTOKES
                SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,UM,O1,O2)
               ENDDO
               T2 = T2 + SUM2
              ENDDO
              LP_UPAR_UP_2(UM,O1,N,K,Q) = T2
             ENDDO
            ENDDO
          ELSE
           DO O1 = 1, NSTOKES
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,IB,O1,O2)
             ENDDO
             T2 = T2 + SUM2
            ENDDO
            LP_UPAR_UP_2(IB,O1,N,K,Q) = T2
           ENDDO
          ENDIF
         ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Diffuse (integrated) contributions only

         IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO O1 = 1, NSTOKES
              T2 = ZERO
              DO L = M, NMOMENTS
               SUM2 = ZERO
               DO O2 = 1, NSTOKES
                SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,UM,O1,O2)
               ENDDO
               T2 = T2 + SUM2
              ENDDO
              LP_UPAR_DN_2(UM,O1,N,K,Q) = T2
             ENDDO
            ENDDO
          ELSE
           DO O1 = 1, NSTOKES
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,IB,O1,O2)
             ENDDO
             T2 = T2 + SUM2
            ENDDO
            LP_UPAR_DN_2(IB,O1,N,K,Q) = T2
           ENDDO
          ENDIF
         ENDIF

!  end parameter and K-layer loop

        ENDDO
       ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LP_UBEAM_SOLUTION

      END MODULE vlidort_lp_solutions

