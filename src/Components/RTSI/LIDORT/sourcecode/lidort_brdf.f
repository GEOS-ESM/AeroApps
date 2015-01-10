C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LIDORT_BRDF_MASTER (master), calling             #
C #              LIDORT_BRDF_MAKER                              #
C #              BRDF_QUADRATURE                                #
C #              BRDF_QUADRATURE_TPZ  (not used)                #
C #                                                             #
C #            LIDORT_BRDF_FOURIER (called by FOURIER)          #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_BRDF_MASTER 

C  Prepares the bidirectional reflectance functions
C  necessary for LIDORT.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Input arguments controlling type of surface
C     output arguments also go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  BRDF functions
C  --------------

      EXTERNAL       ROSSTHIN_FUNCTION
      EXTERNAL       ROSSTHICK_FUNCTION
      EXTERNAL       LISPARSE_FUNCTION
      EXTERNAL       LIDENSE_FUNCTION
      EXTERNAL       HAPKE_FUNCTION
      EXTERNAL       ROUJEAN_FUNCTION
      EXTERNAL       RAHMAN_FUNCTION
      EXTERNAL       COXMUNK_FUNCTION
      EXTERNAL       COXMUNK_FUNCTION_DB

C  Hapke old uses exact DISORT code
C      EXTERNAL       HAPKE_FUNCTION_OLD

C  local arguments
C  ---------------

      INTEGER          K, N
      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )

C  BRDF quadrature
C  ---------------

C  Save these quantities for efficient coding

      CALL BRDF_QUADRATURE

C  Fill BRDF arrays
C  ----------------

      DO K = 1, N_BRDF_KERNELS

C  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO N = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(N) = BRDF_PARAMETERS(K,N)
        ENDDO

C  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, ROSSTHIN_FUNCTION, ROSSTHIN_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, ROSSTHICK_FUNCTION, ROSSTHICK_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, LISPARSE_FUNCTION, LISPARSE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, LIDENSE_FUNCTION, LIDENSE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, HAPKE_FUNCTION, HAPKE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Rahman kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, RAHMAN_FUNCTION, RAHMAN_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, ROUJEAN_FUNCTION, ROUJEAN_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Cox-Munk kernel: (2 free parameters).
C    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_GLITTER_DBMS ) THEN
            CALL LIDORT_BRDF_MAKER
     I          ( K, COXMUNK_FUNCTION, COXMUNK_FUNCTION_DB, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, COXMUNK_FUNCTION, COXMUNK_FUNCTION, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

      ENDDO

C  Finish

      RETURN
      END


      SUBROUTINE BRDF_QUADRATURE

C  Exactly the same as the scalar version

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of inputs.
C     Require arguments controlling type of surface

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  output arguments go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  local variables
C  ---------------

      INTEGER             I, I1, K

C  BRDF quadrature (Gauss-Legendre)
C  ---------------

C  Save these quantities for efficient coding

      NBRDF_HALF = NSTREAMS_BRDF / 2
      CALL GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
        DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
          X_BRDF(I1) = - X_BRDF(I)
          A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = DSQRT(ONE-X_BRDF(I)*X_BRDF(I))
        ENDDO
        DO I = 1, NSTREAMS_BRDF
       X_BRDF(I) = PIE * X_BRDF(I)
       CX_BRDF(I) = DCOS ( X_BRDF(I) )
       SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

C  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_SURFACE_EMISSION .AND.
     &       .NOT. DO_LAMBERTIAN_SURFACE) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BRDF_QUADRATURE_TPZ

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of inputs.
C     Require arguments controlling type of surface

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  output arguments go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  local variables
C  ---------------

      INTEGER             I, I1
      DOUBLE PRECISION DF1, DEL

C  BRDF quadrature (Trapezium)
C  ---------------

C  Save these quantities for efficient coding

      DF1 = DBLE(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
        DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DBLE(I1) * DEL - PIE
        X_BRDF(I) = DBLE(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
        CXE_BRDF(I) = CX_BRDF(I)
        SXE_BRDF(I) = SX_BRDF(I)
      ENDDO
        DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

C  Half space cosine-weight arrays (emission only, non-Lambertian)

C      IF ( DO_SURFACE_EMISSION .AND.
C     &       .NOT. DO_LAMBERTIAN_SURFACE) THEN
C        DO K = 1, NBRDF_HALF
C          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
C        ENDDO
C      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_BRDF_MAKER
     I       ( BRDF_INDEX, BRDF_FUNCTION, BRDF_FUNCTION_DB,
     I         LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )

C  Prepares the bidirectional reflectance scatter matrices

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (model/control/surface) 

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  output arguments go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Module arguments
C  ----------------

C  Index

      INTEGER          BRDF_INDEX

C  BRDF functions (external calls)

      EXTERNAL         BRDF_FUNCTION
      EXTERNAL         BRDF_FUNCTION_DB

C  Local number of parameters and local parameter array

      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )
      
C  local variables
C  ---------------

      INTEGER          I, UI, J, K, KE, IB
      DOUBLE PRECISION MUX, SZASURCOS(MAXBEAMS),SZASURSIN(MAXBEAMS)
      DOUBLE PRECISION PHIANG, COSPHI, SINPHI

C  Usable solar beams

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          MUX =  DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          SZASURCOS(IB) = X0(IB)
          SZASURSIN(IB) = SX0(IB)
        ENDDO
      ENDIF

C  Exact DB calculation
C  --------------------

      IF ( DO_DBCORRECTION ) THEN
        DO K = 1, N_USER_RELAZMS
          PHIANG = USER_RELAZMS(K)
          COSPHI = DCOS(PHIANG*DEG_TO_RAD)
          SINPHI = DSIN(PHIANG*DEG_TO_RAD)
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              CALL BRDF_FUNCTION_DB
     C          ( MAX_BRDF_PARAMETERS, 
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G            SZASURCOS(IB), SZASURSIN(IB),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            PHIANG, COSPHI, SINPHI,
     O            EXACTDB_BRDFUNC(BRDF_INDEX,UI,K,IB) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Exit for SSFULL calculation

      IF ( DO_SSFULL ) RETURN

c      pause

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam

      DO IB = 1, NBEAMS 
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
         CALL BRDF_FUNCTION
     C        ( MAX_BRDF_PARAMETERS, 
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G          SZASURCOS(IB), SZASURSIN(IB), X(I), SX(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          BRDFUNC_0(BRDF_INDEX,I,IB,K) )
        ENDDO
       ENDDO
      ENDDO

C  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION
     C          ( MAX_BRDF_PARAMETERS, 
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G            X(J), SX(J), X(I), SX(I),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            BRDFUNC(BRDF_INDEX,I,J,K) )
          ENDDO
        ENDDO
      ENDDO

C  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        NBRDF_HALF = NSTREAMS_BRDF / 2
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION
     C            ( MAX_BRDF_PARAMETERS,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G              CXE_BRDF(KE), SXE_BRDF(KE), X(I), SX(I),
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              EBRDFUNC(BRDF_INDEX,I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  User-streams outgoing directions
C  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

C  Incident Solar beam

        DO IB = 1, NBEAMS
         DO UI = 1, N_USER_STREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION
     C          ( MAX_BRDF_PARAMETERS, 
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G            SZASURCOS(IB), SZASURSIN(IB),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            USER_BRDFUNC_0(BRDF_INDEX,UI,IB,K) )
          ENDDO
         ENDDO
        ENDDO

C  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION
     C            ( MAX_BRDF_PARAMETERS, 
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I              X(J), SX(J), USER_STREAMS(UI), USER_SINES(UI),
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              USER_BRDFUNC(BRDF_INDEX,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

C  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION
     C            ( MAX_BRDF_PARAMETERS,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I              CXE_BRDF(KE), SXE_BRDF(KE),
     I              USER_STREAMS(UI), USER_SINES(UI), 
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              USER_EBRDFUNC(BRDF_INDEX,UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_BRDF_FOURIER
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER, DELFAC )

C  Prepares Fourier components of the bidirectional reflectance functions

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (model/control/surface) 

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Input and output arguments go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Module arguments
C  ----------------

      LOGICAL             DO_INCLUDE_SURFACE
      LOGICAL             DO_INCLUDE_SURFEMISS
      INTEGER             FOURIER
      DOUBLE PRECISION    DELFAC

C  local variables

      INTEGER             I, UI, J, K, KPHI, KL, IB
      DOUBLE PRECISION    SUM, REFL, HELP, HELP_A

C  Weighted azimuth factor
C  ( Results are stored in commons )

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
       IF ( FOURIER .NE. 0 ) THEN
        DO K = 1, NSTREAMS_BRDF
          BRDF_AZMFAC(K) = A_BRDF(K) * DCOS ( FOURIER * X_BRDF(K) )
        ENDDO
       ELSE
        DO K = 1, NSTREAMS_BRDF
          BRDF_AZMFAC(K) = A_BRDF(K)
        ENDDO
       ENDIF
      ENDIF

C  surface factor

      HELP = HALF * DELFAC

C  Skip Fourier components if Lambertian

      IF ( DO_LAMBERTIAN_SURFACE ) GO TO 778

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam (direct beam reflections)

      DO KL = 1, N_BRDF_KERNELS
       IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
        DO IB = 1, NBEAMS
         IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          DO I = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
             SUM  = SUM + BRDFUNC_0(KL,I,IB,K)*BRDF_AZMFAC(K)
            ENDDO
            BIREFLEC_0(KL,I,IB) = SUM * HELP
          ENDDO
         ENDIF
        ENDDO
       ENDIF
      ENDDO

C  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO KL = 1, N_BRDF_KERNELS
          IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + BRDFUNC(KL,I,J,K) * BRDF_AZMFAC(K)
                ENDDO
                BIREFLEC(KL,I,J) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

C  debug information

c      IF ( DO_DEBUG_WRITE ) THEN
c        WRITE(555,'(A)')'BRDF_1 Fourier 0 quad values'
c        IF ( FOURIER .EQ. 0 ) THEN
c          DO I = 1, NSTREAMS
c          WRITE(555,'(1PE12.5,3x,1P10E12.5)')
c     &     BIREFLEC_0(1,I,1),(BIREFLEC(1,I,J),J=1,NSTREAMS)
c         ENDDO
c        ENDIF
c      ENDIF

C  albedo check, always calculate the spherical albedo.
C   (Plane albedo calculations are commented out)

      IF ( FOURIER .EQ. 0 ) THEN
        DO KL = 1, N_BRDF_KERNELS
         IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

C  ............................... Plane albedo calculations
C          SUM = ZERO
C          DO I = 1, NSTREAMS
C            SUM = SUM + BIREFLEC_0(KL,I) * AX(I)
C          ENDDO
C          SUM = SUM * TWO
C          write(*,*)0,BRDF_FACTORS(KL),sum
C          DO J = 1, NSTREAMS
C            SUM = ZERO
C            DO I = 1, NSTREAMS
C              SUM = SUM + BIREFLEC(KL,I,J) * AX(I)
C            ENDDO
C            SUM = SUM * TWO
C            write(*,*)j,BRDF_FACTORS(KL),sum
C          ENDDO

C  ............................... Spherical albedo calculation

          HELP_A = ZERO
          DO I = 1, NSTREAMS
            SUM = ZERO
            DO J = 1, NSTREAMS
               SUM = SUM + BIREFLEC(KL,I,J) * AX(J)
            ENDDO
            HELP_A = HELP_A + SUM * AX(I)
          ENDDO
          SPHERICAL_ALBEDO(KL) = HELP_A*FOUR
         ENDIF
        ENDDO
      ENDIF

C  User-streams outgoing directions
C  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

C  Incident Solar beam (direct beam reflections)

       DO KL = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
         DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
           DO UI = 1, N_USER_STREAMS
             SUM = ZERO
             DO K = 1, NSTREAMS_BRDF
              SUM = SUM+USER_BRDFUNC_0(KL,UI,IB,K)*BRDF_AZMFAC(K)
             ENDDO
             USER_BIREFLEC_0(KL,UI,IB) = SUM * HELP
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO

C  incident quadrature directions (surface multiple reflections)

       IF ( DO_INCLUDE_SURFACE ) THEN
        DO KL = 1, N_BRDF_KERNELS
         IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
          DO UI = 1, N_USER_STREAMS
           DO J = 1, NSTREAMS
             SUM = ZERO
             DO K = 1, NSTREAMS_BRDF
              SUM = SUM + USER_BRDFUNC(KL,UI,J,K)*BRDF_AZMFAC(K)
             ENDDO
             USER_BIREFLEC(KL,UI,J) = SUM * HELP
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDIF

      ENDIF

C  Continuation point

 778  CONTINUE

C  Emissivity
C  ----------

C  Assumed to exist only for the total intensity
C        (first element of Stokes Vector) - is this right ??????

      IF ( DO_INCLUDE_SURFEMISS ) THEN

C  Initialise

        DO I = 1, NSTREAMS
          EMISSIVITY(I) = ONE
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            USER_EMISSIVITY(UI) = ONE
          ENDDO
        ENDIF

C  Start loop over kernels

        DO KL = 1, N_BRDF_KERNELS

C  Lambertian case

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

            DO I = 1, NSTREAMS
              ALBEDO_EMISSIVITY(I,KL) = BRDF_FACTORS(KL)
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                ALBEDO_USER_EMISSIVITY(UI,KL) = BRDF_FACTORS(KL)
              ENDDO
            ENDIF

C  bidirectional reflectance

          ELSE

C  Quadrature polar directions

            DO I = 1, NSTREAMS
              REFL = ZERO
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = ZERO
                DO K = 1, NBRDF_HALF
                  SUM = SUM + EBRDFUNC(KL,I,K,KPHI) * BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              ALBEDO_EMISSIVITY(I,KL) = REFL * BRDF_FACTORS(KL)
            ENDDO

C   user-defined polar directions

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                REFL = ZERO
                DO KPHI= 1, NSTREAMS_BRDF
                  SUM = ZERO
                  DO K = 1, NBRDF_HALF
                    SUM = SUM+USER_EBRDFUNC(KL,UI,K,KPHI)*BAX_BRDF(K)
                  ENDDO
                  REFL = REFL + A_BRDF(KPHI) * SUM
                ENDDO
                ALBEDO_USER_EMISSIVITY(UI,KL) = REFL*BRDF_FACTORS(KL)
              ENDDO
            ENDIF

          ENDIF

C  end loop over kernels

        ENDDO

C  Total emissivities

        DO KL = 1, N_BRDF_KERNELS
          DO I = 1, NSTREAMS
            EMISSIVITY(I) = EMISSIVITY(I) -
     &               ALBEDO_EMISSIVITY(I,KL)
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_EMISSIVITY(UI) = USER_EMISSIVITY(UI) -
     &               ALBEDO_USER_EMISSIVITY(UI,KL)
            ENDDO
          ENDIF
        ENDDO

C  end emissivity clause

      ENDIF

C  Finish

      RETURN
      END
