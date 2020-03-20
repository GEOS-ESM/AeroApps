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
C #            LIDORT_L_BRDF_MASTER (master)                    #
C #              LIDORT_BRDF_MAKER_PLUS                         #
C #                                                             #
C #            LIDORT_L_BRDF_FOURIER                            #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_BRDF_MASTER 

C  Prepares (linearizations of) the bidirectional reflectance functions
C  necessary for LIDORT.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Include files of LIDORT linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Input arguments controlling type of surface
C     output arguments also go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  BRDF functions
C  --------------

C  ordinary BRDF without derivatives

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

C  BRDFs with derivatives

      EXTERNAL       LISPARSE_FUNCTION_PLUS
      EXTERNAL       LIDENSE_FUNCTION_PLUS
      EXTERNAL       HAPKE_FUNCTION_PLUS
      EXTERNAL       RAHMAN_FUNCTION_PLUS
      EXTERNAL       COXMUNK_FUNCTION_PLUS
      EXTERNAL       COXMUNK_FUNCTION_PLUS_DB

C  local arguments
C  ---------------

      INTEGER          K, Q
      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL          LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

C  BRDF quadrature
C  ---------------

C  Save these quantities for efficient coding

      CALL BRDF_QUADRATURE

C  Fill BRDF arrays
C  ----------------

      DO K = 1, N_BRDF_KERNELS

C  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO Q = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(Q) = BRDF_PARAMETERS(K,Q)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO Q = 1, MAX_BRDF_PARAMETERS
            LOCAL_BRDF_DERIVS(Q) = DO_KERNEL_PARAMS_WFS(K,Q)
          ENDDO
        ENDIF

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
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, LISPARSE_FUNCTION_PLUS, LISPARSE_FUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, LISPARSE_FUNCTION,  LISPARSE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, LIDENSE_FUNCTION_PLUS, LIDENSE_FUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, LIDENSE_FUNCTION, LIDENSE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, HAPKE_FUNCTION_PLUS, HAPKE_FUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, HAPKE_FUNCTION, HAPKE_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER
     I          ( K, ROUJEAN_FUNCTION, ROUJEAN_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Rahman kernel: (3 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, RAHMAN_FUNCTION_PLUS, RAHMAN_FUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, RAHMAN_FUNCTION, RAHMAN_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
         IF (  DO_GLITTER_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, COXMUNK_FUNCTION_PLUS, COXMUNK_FUNCTION_PLUS_DB,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, COXMUNK_FUNCTION, COXMUNK_FUNCTION_DB,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS
     I          ( K, COXMUNK_FUNCTION_PLUS, COXMUNK_FUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL LIDORT_BRDF_MAKER
     I          ( K, COXMUNK_FUNCTION, COXMUNK_FUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ENDIF
        ENDIF

      ENDDO

C  Finish

      RETURN
      END

c

      SUBROUTINE LIDORT_BRDF_MAKER_PLUS
     I       ( BRDF_INDEX, BRDF_FUNCTION_PLUS, BRDF_FUNCTION_PLUS_DB,
     I         LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I         LOCAL_BRDF_DERIVS )

C  Prepares the bidirectional reflectance function derivatives

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of LIDORT weighting function input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include file of basic reflectance values

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  output arguments go in here.

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  Module arguments
C  ----------------

C  Index

      INTEGER          BRDF_INDEX

C  BRDF functions (external calls)

      EXTERNAL         BRDF_FUNCTION_PLUS
      EXTERNAL         BRDF_FUNCTION_PLUS_DB

C  Local number of parameters and local parameter array

      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL          LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )
            
C  local variables

      INTEGER          I, UI, J, K, KE, Q, B, IB
      DOUBLE PRECISION FUNC, DFUNC ( MAX_BRDF_PARAMETERS )
      DOUBLE PRECISION MUX, SZASURCOS(MAXBEAMS),SZASURSIN(MAXBEAMS)
      DOUBLE PRECISION PHIANG, COSPHI, SINPHI

C  Usable solar beams. Repeat code in the  standard module.

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

C  shorthand

      B = BRDF_INDEX

C  Exact DB calculation
C  --------------------

      IF ( DO_DBCORRECTION ) THEN
        DO K = 1, N_USER_RELAZMS
          PHIANG = USER_RELAZMS(K)
          COSPHI = DCOS(PHIANG*DEG_TO_RAD)
          SINPHI = DSIN(PHIANG*DEG_TO_RAD)
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              CALL BRDF_FUNCTION_PLUS_DB
     C        ( MAX_BRDF_PARAMETERS,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          SZASURCOS(IB), SZASURSIN(IB),
     G          USER_STREAMS(UI), USER_SINES(UI),
     G          PHIANG, COSPHI, SINPHI,
     O          FUNC, DFUNC )
              EXACTDB_BRDFUNC(B,UI,K,IB) = FUNC
              DO Q = 1, LOCAL_BRDF_NPARS
                D_EXACTDB_BRDFUNC(B,Q,UI,K,IB)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Exit for SSFULL calculation

      IF ( DO_SSFULL ) RETURN

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam

      DO IB = 1, NBEAMS
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
          CALL BRDF_FUNCTION_PLUS
     C        ( MAX_BRDF_PARAMETERS,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          SZASURCOS(IB), SZASURSIN(IB), X(I), SX(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          FUNC, DFUNC )
          BRDFUNC_0(B,I,IB,K) = FUNC
          DO Q = 1, LOCAL_BRDF_NPARS
            D_BRDFUNC_0(B,Q,I,IB,K)  = DFUNC(Q)
          ENDDO
        ENDDO
       ENDDO
      ENDDO

C  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
          CALL BRDF_FUNCTION_PLUS
     C        ( MAX_BRDF_PARAMETERS,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          X(J), SX(J), X(I), SX(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          FUNC, DFUNC )
            BRDFUNC(B,I,J,K) = FUNC
            DO Q = 1, LOCAL_BRDF_NPARS
              D_BRDFUNC(B,Q,I,J,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        NBRDF_HALF = NSTREAMS_BRDF / 2
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS
     C          ( MAX_BRDF_PARAMETERS,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C            LOCAL_BRDF_DERIVS,
     G            CXE_BRDF(KE), SXE_BRDF(KE), X(I), SX(I),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            FUNC, DFUNC )
              EBRDFUNC(B,I,KE,K) = FUNC
              DO Q = 1, LOCAL_BRDF_NPARS
                D_EBRDFUNC(B,Q,I,KE,K)  = DFUNC(Q)
              ENDDO
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
            CALL BRDF_FUNCTION_PLUS
     C          ( MAX_BRDF_PARAMETERS,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C            LOCAL_BRDF_DERIVS,
     G            SZASURCOS(IB), SZASURSIN(IB),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            FUNC, DFUNC )
            USER_BRDFUNC_0(B,UI,IB,K) = FUNC
            DO Q = 1, LOCAL_BRDF_NPARS
              D_USER_BRDFUNC_0(B,Q,UI,IB,K)  = DFUNC(Q)
            ENDDO
          ENDDO
         ENDDO
        ENDDO

C  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS
     C          ( MAX_BRDF_PARAMETERS,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C            LOCAL_BRDF_DERIVS,
     I            X(J), SX(J), USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            FUNC, DFUNC )
              USER_BRDFUNC(B,UI,J,K) = FUNC
              DO Q = 1, LOCAL_BRDF_NPARS
                D_USER_BRDFUNC(B,Q,UI,J,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

C  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_FUNCTION_PLUS
     C            ( MAX_BRDF_PARAMETERS,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C              LOCAL_BRDF_DERIVS,
     I              CXE_BRDF(KE), SXE_BRDF(KE),
     I              USER_STREAMS(UI), USER_SINES(UI), 
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O             FUNC, DFUNC )
                USER_EBRDFUNC(B,UI,KE,K) = FUNC
                DO Q = 1, LOCAL_BRDF_NPARS
                  D_USER_EBRDFUNC(B,Q,UI,KE,K)  = DFUNC(Q)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_BRDF_FOURIER
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER_COMPONENT,
     I      SURFACE_FACTOR )

C  Prepares Fourier components of the bidirectional reflectance derivatives

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include file of LIDORT weighting function input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include files of basic reflectance values

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  output arguments go in here.

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  Module arguments
C  ----------------

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables

      INTEGER          I, UI, J, K, KPHI, Q, B, IB
      DOUBLE PRECISION SUM, REFL, HELP

C  factor

      HELP = HALF * SURFACE_FACTOR

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam (direct beam reflections)

      DO B = 1, N_BRDF_KERNELS
       IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
        DO IB = 1, NBEAMS
         IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          DO Q = 1, N_BRDF_PARAMETERS(B)
           IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
            DO I = 1, NSTREAMS
             SUM = ZERO
             DO K = 1, NSTREAMS_BRDF
              SUM = SUM + D_BRDFUNC_0(B,Q,I,IB,K) * BRDF_AZMFAC(K)
             ENDDO
             D_BIREFLEC_0(B,Q,I,IB) = SUM * HELP
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDIF
      ENDDO
      
C  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
       DO B = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
         DO Q = 1, N_BRDF_PARAMETERS(B)
          IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
           DO I = 1, NSTREAMS
            DO J = 1, NSTREAMS
             SUM = ZERO
             DO K = 1, NSTREAMS_BRDF
              SUM = SUM + D_BRDFUNC(B,Q,I,J,K) * BRDF_AZMFAC(K)
             ENDDO
             D_BIREFLEC(B,Q,I,J) = SUM * HELP
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDIF

C  User-streams outgoing directions
C  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

C  Incident Solar beam (direct beam reflections)

       DO B = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
         DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
           DO Q = 1, N_BRDF_PARAMETERS(B)
            IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
             DO UI = 1, N_USER_STREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM +
     &            D_USER_BRDFUNC_0(B,Q,UI,IB,K)*BRDF_AZMFAC(K)
              ENDDO
              D_USER_BIREFLEC_0(B,Q,UI,IB) = SUM * HELP
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO

C  incident quadrature directions (surface multiple reflections)

       IF ( DO_INCLUDE_SURFACE ) THEN
        DO B = 1, N_BRDF_KERNELS
         IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
          DO Q = 1, N_BRDF_PARAMETERS(B)
           IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
            DO UI = 1, N_USER_STREAMS
             DO J = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM + D_USER_BRDFUNC(B,Q,UI,J,K)*BRDF_AZMFAC(K)
              ENDDO
              D_USER_BIREFLEC(B,Q,UI,J) = SUM * HELP
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDIF

      ENDIF

C  Emissivity
C  ----------

      IF ( DO_INCLUDE_SURFEMISS ) THEN

C  Start loop over kernels

       DO B = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
         DO Q = 1, N_BRDF_PARAMETERS(B)
          IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN

C  Quadrature polar directions

           DO I = 1, NSTREAMS
            REFL = ZERO
            DO KPHI= 1, NSTREAMS_BRDF
             SUM = ZERO
             DO K = 1, NBRDF_HALF
              SUM = SUM + D_EBRDFUNC(B,Q,I,K,KPHI) * BAX_BRDF(K)
             ENDDO
             REFL = REFL + A_BRDF(KPHI) * SUM
            ENDDO
            D_EMISSIVITY(B,Q,I) = - REFL * BRDF_FACTORS(B)
           ENDDO

C   user-defined polar directions

           IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
             REFL = ZERO
             DO KPHI= 1, NSTREAMS_BRDF
              SUM = ZERO
              DO K = 1, NBRDF_HALF
               SUM = SUM + D_USER_EBRDFUNC(B,Q,UI,K,KPHI) * BAX_BRDF(K)
              ENDDO
              REFL = REFL + A_BRDF(KPHI) * SUM
             ENDDO
             D_USER_EMISSIVITY(B,Q,UI) = - REFL*BRDF_FACTORS(B)
            ENDDO
           ENDIF

C  parameter and kernel loop

          ENDIF
         ENDDO
        ENDIF
       ENDDO

C  end emissivity clause

      ENDIF

C  Finish

      RETURN
      END
