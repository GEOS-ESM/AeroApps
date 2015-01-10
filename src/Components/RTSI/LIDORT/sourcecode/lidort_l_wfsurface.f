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

C  Subroutines in this Module

C            KFACTOR_BRDFWF (master)
C            KPARAMS_BRDFWF (master)

C              BOA_KFACTOR_SURFACEWF
C              BOA_KPARAMS_SURFACEWF
C              UPUSER_SURFACEWF
C              DNUSER_SURFACEWF
C              MIFLUX_SURFACEWF

C ###############################################################

      SUBROUTINE KFACTOR_BRDFWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          DO_INCLUDE_DIRECTBEAM,
     I          FOURIER_COMPONENT, IBEAM,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER,
     I          BRDF_INDEX, WF_INDEX,
     O          STATUS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution and setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of Linearized solution variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier surface factor, Fourier number, Beam index, Flux multiplier

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Weighting function index number

      INTEGER          WF_INDEX

C  BRDF index (indicates which BRDF kernel is active)

      INTEGER          BRDF_INDEX

C  Output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Other local variables

      INTEGER          N, I, I1, K
      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Initialise status

       STATUS = LIDORT_SUCCESS

C  Avoid boundary value problem for thermal-transmittance

       IF ( DO_THERMAL_TRANSONLY ) GO TO 699

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'
C    Include flag here controlled by DO_REFLECTED_DIRECTBEAM

      CALL KFACTOR_WF_COLSETUP
     I       (  DO_REFLECTED_DIRECTBEAM(IBEAM),
     I          DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR,
     I          IBEAM, BRDF_INDEX )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (multilayer) in LIDORT_KFACTOR_BRDFWF'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            NCON_ALB(I,N) = COL2_WFALB(LCONMASK(I,N),1)
            PCON_ALB(I,N) = COL2_WFALB(MCONMASK(I,N),1)
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2_WFALB, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (one layer) in LIDORT_KFACTOR_BRDFWF'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          NCON_ALB(I,N) = SCOL2_WFALB(I,1)
          PCON_ALB(I,N) = SCOL2_WFALB(I1,1)
        ENDDO

C  end clause
 
      ENDIF

C  debug------------------------------------------
        if ( do_debug_write.and.fourier_component.eq.0 ) then
         DO N = 1, NLAYERS
          DO K = 1, NSTREAMS
           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
     &                LCON(K,N), MCON(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N)
          ENDDO
         ENDDO
        ENDIF

C  Continuation point (thermal transmittance only)

 699    continue

C  Get the weighting functions
C  ---------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE)

        CALL BOA_KFACTOR_SURFACEWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT,
     I          IBEAM, BRDF_INDEX, 
     O          L_BOA_SOURCE,
     O          L_BOA_THTONLY_SOURCE )

C  Upwelling Albedo weighting functions

        CALL UPUSER_SURFACEWF
     I      ( DO_INCLUDE_MVOUTPUT,
     I        FLUX_MULTIPLIER, IBEAM,
     I        WF_INDEX, BRDF_INDEX,
     I        L_BOA_SOURCE,
     I        L_BOA_THTONLY_SOURCE )

      ENDIF

C  Downwelling Albedo weighting functions

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF
     I       ( DO_INCLUDE_MVOUTPUT,
     I         FLUX_MULTIPLIER, IBEAM,
     I         WF_INDEX, BRDF_INDEX )
      ENDIF

C  mean value output

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL MIFLUX_SURFACEWF ( IBEAM, WF_INDEX )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE KPARAMS_BRDFWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          DO_INCLUDE_DIRECTBEAM,
     I          FOURIER_COMPONENT, IBEAM,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER,
     I          BRDF_INDEX, PAR_INDEX, WF_INDEX,
     O          STATUS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution and setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of Linearized solution variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier surface factor, Fourier number, Beam index, Flux multiplier

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Weighting function index number

      INTEGER          WF_INDEX

C  BRDF index (indicates which BRDF kernel)

      INTEGER          BRDF_INDEX

C  Parameter index

      INTEGER          PAR_INDEX

C  Output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Other local variables

      INTEGER          N, I, I1, K
      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Initialise status

       STATUS = LIDORT_SUCCESS

C  Avoid boundary value problem for thermal-transmittance

       IF ( DO_THERMAL_TRANSONLY ) GO TO 699

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'
C    Include flag here controlled by DO_REFLECTED_DIRECTBEAM

      CALL KPARAMS_WF_COLSETUP
     I       (  DO_REFLECTED_DIRECTBEAM(IBEAM),
     I          DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR,
     I          IBEAM, BRDF_INDEX, PAR_INDEX )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (multilayer) in LIDORT_KPARAMS_BRDFWF'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            NCON_ALB(I,N) = COL2_WFALB(LCONMASK(I,N),1)
            PCON_ALB(I,N) = COL2_WFALB(MCONMASK(I,N),1)
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2_WFALB, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (one layer) in LIDORT_KPARAMS_BRDFWF'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          NCON_ALB(I,N) = SCOL2_WFALB(I,1)
          PCON_ALB(I,N) = SCOL2_WFALB(I1,1)
        ENDDO

C  end clause
 
      ENDIF

C  debug------------------------------------------
        if ( do_debug_write.and.fourier_component.eq.0 ) then
         DO N = 1, NLAYERS
          DO K = 1, NSTREAMS
           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,1,K,
     &                LCON(K,N), MCON(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N)
          ENDDO
         ENDDO
        ENDIF

C  Continuation point (thermal transmittance only)

 699    continue

C  Get the weighting functions
C  ---------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE) for these weighting functions

        CALL BOA_KPARAMS_SURFACEWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, PAR_INDEX,
     O          L_BOA_SOURCE,
     O          L_BOA_THTONLY_SOURCE )

C  Upwelling Surface weighting functions (w.r.t Kernel amplitudes)

        CALL UPUSER_SURFACEWF
     I      ( DO_INCLUDE_MVOUTPUT,
     I        FLUX_MULTIPLIER,
     I        IBEAM, WF_INDEX, BRDF_INDEX,
     I        L_BOA_SOURCE,
     I        L_BOA_THTONLY_SOURCE )

      ENDIF

C  Downwelling Surface weighting functions (w.r.t Kernel amplitudes)

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF
     I       ( DO_INCLUDE_MVOUTPUT,
     I         FLUX_MULTIPLIER, IBEAM,
     I         WF_INDEX, BRDF_INDEX )
      ENDIF

C  mean value output

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL MIFLUX_SURFACEWF ( IBEAM, WF_INDEX )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BOA_KFACTOR_SURFACEWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, 
     O          L_BOA_SOURCE,
     O          L_BOA_THTONLY_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of Linearized solution variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  Surface inclusion flag

      LOGICAL          DO_INCLUDE_SURFACE

C  directbeam inclusion flag

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  MV output inclusion flag

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM
C  kernel index

      INTEGER          BRDF_INDEX

C  output

      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          UM, I, J, AA, KL, N, IB
      DOUBLE PRECISION INTEGRAND(MAXSTREAMS), REFLEC, SUM
      DOUBLE PRECISION FACTOR_BRDF, H1, H2

C  Initialise
C  ----------

      N  = NLAYERS
      IB = IBEAM

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  initialise Derivative of BOA source function

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_BOA_SOURCE(UM) = ZERO
        ENDDO
      ENDIF

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          L_BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

C  Return if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Skip diffuse-field variation for thermal transmittance-only

      IF ( DO_THERMAL_TRANSONLY. OR. .NOT. DO_USER_STREAMS ) GO TO 599

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

      DO I = 1, NSTREAMS
        SUM = ZERO
        DO AA = 1, NSTREAMS
          H1  = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
          H2  = PCON_ALB(AA,N)*XNEG(I,AA,N) 
          SUM = SUM + H1 + H2
        ENDDO
        INTEGRAND(I) = SUM * AX(I)
      ENDDO

C  integrated BRDF reflectance term (loop over all kernels)
C  --------------------------------------------------------

      DO KL = 1, N_BRDF_KERNELS

C  .. Either, integrate Lambertian case, same for all user-streams

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) .AND.
     &                 FOURIER_COMPONENT.EQ.0 ) THEN

          FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)
          SUM = ZERO
          DO J = 1, NSTREAMS
             SUM = SUM + INTEGRAND(J)
          ENDDO
          REFLEC = FACTOR_BRDF * SUM
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + REFLEC
          ENDDO

C  .. Or, integrate with BRDF kernel function at user angles

        ELSE
  
          FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            SUM = ZERO
            DO J = 1, NSTREAMS
              SUM = SUM + INTEGRAND(J) * USER_BIREFLEC(KL,UM,J)
            ENDDO
            REFLEC = FACTOR_BRDF * SUM
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + REFLEC
          ENDDO

        ENDIF

C  End loop over all kernels

      ENDDO

C  Continuation point for avoiding diffuse field computation

 599  continue

C  Contributions due to direct variation of kernel factor
C  ------------------------------------------------------

C  Add linearization due to kernel factor variation, diffuse term

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM)
     &              + ALBEDO_BOA_SOURCE  ( UM, BRDF_INDEX )
        ENDDO
      ENDIF

C  Thermal tranmsittance only (additional quadrature terms if flagged)

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
             L_BOA_THTONLY_SOURCE(I) = L_BOA_THTONLY_SOURCE(I) +
     &           ALBEDO_BOA_THTONLY_SOURCE(I,BRDF_INDEX)
        ENDDO
      ENDIF

C  Due to factor variation of direct beam (only if flagged)

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM)
     &              + A_USER_DIRECT_BEAM ( UM, IB, BRDF_INDEX )
        ENDDO
      ENDIF

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)
C   Needs testing. Checks out. 2 March 2007

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_USER_STREAMS ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) - 
     &           SURFBB * ALBEDO_USER_EMISSIVITY(UM,BRDF_INDEX)
          ENDDO
        ENDIF
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
             L_BOA_THTONLY_SOURCE(I) = L_BOA_THTONLY_SOURCE(I) -
     &           SURFBB * ALBEDO_EMISSIVITY(I,BRDF_INDEX)
          ENDDO
        ENDIF
      ENDIF

C  debug

c      if ( fourier_component.le.4.and.ibeam.eq.1) then
c         write(99,'(1p4e17.9)')fourier_component,L_BOA_SOURCE(1)
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BOA_KPARAMS_SURFACEWF
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, PAR_INDEX,
     O          L_BOA_SOURCE,
     O          L_BOA_THTONLY_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of Linearized solution and reflectance variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  Module arguments
C  ----------------

C  Surface inclusion flag

      LOGICAL          DO_INCLUDE_SURFACE

C  directbeam inclusion flag

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  MV output inclusion flag

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM
C  kernel index

      INTEGER          BRDF_INDEX

C  parameter index

      INTEGER          PAR_INDEX

C  output

      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          UM, I, J, AA, KL, N, Q, B, IB
      DOUBLE PRECISION INTEGRAND(MAXSTREAMS), REFLEC, SUM
      DOUBLE PRECISION FACTOR_BRDF, H1, H2, PAR, REFL_ATTN

C  Initialise
C  ----------

C  Local indices

      N  = NLAYERS
      IB = IBEAM
      B = BRDF_INDEX
      Q = PAR_INDEX

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  initialise Derivative of BOA source function

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_BOA_SOURCE(UM) = ZERO
        ENDDO
      ENDIF

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          L_BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

C  Return if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Skip diffuse-field variation for thermal transmittance-only

      IF ( DO_THERMAL_TRANSONLY. OR. .NOT. DO_USER_STREAMS ) GO TO 599

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

      DO I = 1, NSTREAMS
        SUM = ZERO
        DO AA = 1, NSTREAMS
          H1  = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
          H2  = PCON_ALB(AA,N)*XNEG(I,AA,N) 
          SUM = SUM + H1 + H2
        ENDDO
        INTEGRAND(I) = SUM * AX(I)
      ENDDO

C  integrated BRDF reflectance term (loop over all kernels)
C  --------------------------------------------------------

      DO KL = 1, N_BRDF_KERNELS

C  .. Either, integrate Lambertian case, same for all user-streams

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) .AND.
     &                 FOURIER_COMPONENT.EQ.0 ) THEN

          FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)
          SUM = ZERO
          DO J = 1, NSTREAMS
             SUM = SUM + INTEGRAND(J)
          ENDDO
          REFLEC = FACTOR_BRDF * SUM
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + REFLEC
          ENDDO

C  .. Or, integrate with BRDF kernel function at user angles

        ELSE
  
          FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            SUM = ZERO
            DO J = 1, NSTREAMS
              SUM = SUM + INTEGRAND(J) * USER_BIREFLEC(KL,UM,J)
            ENDDO
            REFLEC = FACTOR_BRDF * SUM
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + REFLEC
          ENDDO

        ENDIF

C  End loop over all kernels

      ENDDO

C  Continuation point for avoiding diffuse field computation

 599  continue

C  Contributions due to direct variation of kernel parameter
C  ---------------------------------------------------------

C  variation of the diffuse reflectance function

      PAR = BRDF_PARAMETERS(B,Q)
      FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(B) * PAR
      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SUM = ZERO
          DO J = 1, NSTREAMS
            SUM = SUM + IDOWNSURF(J) * D_USER_BIREFLEC(B,Q,UM,J)
          ENDDO
          REFLEC = FACTOR_BRDF * SUM
          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + REFLEC
        ENDDO
      ENDIF
      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          SUM = ZERO
          DO J = 1, NSTREAMS
            SUM = SUM + IDOWNSURF(J) * D_BIREFLEC(B,Q,I,J)
          ENDDO
          REFLEC = FACTOR_BRDF * SUM
          L_BOA_THTONLY_SOURCE(I) = L_BOA_THTONLY_SOURCE(I) + REFLEC
        ENDDO
      ENDIF

C  variation of direct beam reflectance (only if flagged)

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        REFL_ATTN = BRDF_FACTORS(B) * ATMOS_ATTN(IB)
        FACTOR_BRDF = PAR * REFL_ATTN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) +
     &                     FACTOR_BRDF * D_USER_BIREFLEC_0(B,Q,UM,IB)
        ENDDO
      ENDIF

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)
C  Emissivity derivative is straight-up (not normalized) --> multiply by PAR
C  Emissivity derivative has correct sign already worked in
C   Linearization code tested with Thermal supplement, 02 March 2007.

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_USER_STREAMS ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) +
     &                  SURFBB * PAR * D_USER_EMISSIVITY(B,Q,UM)
          ENDDO
        ENDIF
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
            L_BOA_THTONLY_SOURCE(I) = L_BOA_THTONLY_SOURCE(I) +
     &                  SURFBB * PAR * D_EMISSIVITY(B,Q,I)
          ENDDO
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE UPUSER_SURFACEWF
     I      ( DO_INCLUDE_MVOUTPUT,
     I        FLUX_MULTIPLIER,
     I        IBEAM, WF_INDEX, BRDF_INDEX,
     I        L_BOA_SOURCE,
     I        L_BOA_THTONLY_SOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution, multipliers variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of Linearized solution variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux multiplier

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local control flag

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Beam index

      INTEGER          IBEAM

C  weighting function index

      INTEGER          WF_INDEX

C  BRDF kernel index 

      INTEGER          BRDF_INDEX

C  derivatives of reflected surface upwelling intensity

      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, NL
      INTEGER          UTA, UM, IUM, IQD, I, I1, UT, AA, IB, K
      DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM, THELP

C  Initial section
C  ---------------

C  Local index

      IB = IBEAM

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_OUT_USERTAUS
          DO UM = 1, n_USER_STREAMS
            IUM = USEROUTPUT_INDEX(UM)
            SURFACEWF_F(WF_INDEX,UTA,IUM,IB,UPIDX) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_CUMUL_SOURCE(UM) = L_BOA_SOURCE(UM)
        ENDDO
      ENDIF

C  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1

C  Thermal transmittance only

            IF ( DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
               L_LAYER_SOURCE(UM) = ZERO
             ENDDO
            ENDIF

C  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO AA = 1, NSTREAMS
                NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XPOS(UM,AA,N)
                PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XNEG(UM,AA,N)
              ENDDO
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                H2 = PCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              L_LAYER_SOURCE(UM) = SHOM
             ENDDO
            ENDIF

C  Recursion

            DO UM = LOCAL_UM_START, N_USER_STREAMS
              L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM)
            ENDDO

C  end layer loop

          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output, mean value + flux output

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
           IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_ALB(AA,N)*XPOS(I1,AA,N)*T_UTDN_EIGEN(AA,UT)
                H2 = PCON_ALB(AA,N)*XNEG(I1,AA,N)*T_UTUP_EIGEN(AA,UT)
                SHOM = SHOM + H1 + H2
              ENDDO
              QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) = 
     &                       FLUX_MULTIPLIER * SHOM
            ENDDO
           ELSE
            DO I = 1, NSTREAMS
              THELP = L_BOA_THTONLY_SOURCE(I)
              DO K = NLAYERS, N+1, -1
                THELP = THELP*T_DELT_DISORDS(I,K)
              ENDDO
              THELP = THELP*T_DISORDS_UTUP(I,UT)
              QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) = 
     &                       FLUX_MULTIPLIER * THELP
            ENDDO
           ENDIF
          ENDIF

C  Copy results if quadrature output required (unlikely)

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
               SURFACEWF_F(WF_INDEX,UTA,IQD,IB,UPIDX) = 
     &                     QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX)
            ENDDO
          ENDIF

C  User-defined stream output

          IF ( DO_USER_STREAMS ) THEN

C  Thermal transmittance only

            IF ( DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
               L_LAYER_SOURCE(UM) = ZERO
             ENDDO
            ENDIF

C  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO AA = 1, NSTREAMS
                NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XPOS(UM,AA,N)
                PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XNEG(UM,AA,N)
              ENDDO
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_HELP(UM,AA) * UT_HMULT_UD(AA,UM,UT)
                H2 = PCON_HELP(UM,AA) * UT_HMULT_UU(AA,UM,UT)
                SHOM = SHOM + H1 + H2
              ENDDO
              L_LAYER_SOURCE(UM) = SHOM
             ENDDO
            ENDIF

C  assign final albedo weighting function

            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
               SURFACEWF_F(WF_INDEX,UTA,IUM,IB,UPIDX) = 
     &                        FLUX_MULTIPLIER * L_FINAL_SOURCE
            ENDDO

          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output, mean value + flux output

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

            NL = NLEVEL
            N = NL + 1

C  For the lowest level

            IF ( NLEVEL .EQ. NLAYERS ) THEN

              IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 =NCON_ALB(AA,NL)*XPOS(I1,AA,NL)*T_DELT_EIGEN(AA,NL)
                  H2 =PCON_ALB(AA,NL)*XNEG(I1,AA,NL)
                  SHOM = SHOM + H1 + H2
                ENDDO
                 QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) =
     &                          FLUX_MULTIPLIER * SHOM
               ENDDO
              ELSE
               DO I = 1, NSTREAMS
                 THELP = L_BOA_THTONLY_SOURCE(I)
                 QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) =
     &                          FLUX_MULTIPLIER * THELP
               ENDDO
              ENDIF

C  For other levels in the atmosphere

            ELSE

              IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_ALB(AA,N)*XPOS(I1,AA,N)
                  H2 = PCON_ALB(AA,N)*XNEG(I1,AA,N)*T_DELT_EIGEN(AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                 QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) = 
     &                          FLUX_MULTIPLIER * SHOM
               ENDDO
              ELSE
               DO I = 1, NSTREAMS
                 THELP = L_BOA_THTONLY_SOURCE(I)
                 DO K = NLAYERS, N, -1
                   THELP = THELP*T_DELT_DISORDS(I,K)
                 ENDDO
                 QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX) =
     &                          FLUX_MULTIPLIER * THELP
               ENDDO
              ENDIF

            ENDIF

          ENDIF

C  Copy results if quadrature output required (unlikely)

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
               SURFACEWF_F(WF_INDEX,UTA,IQD,IB,UPIDX) = 
     &                   QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX)
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
               SURFACEWF_F(WF_INDEX,UTA,IUM,IB,UPIDX) = 
     &                      FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM)
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  debug

c      do i = 1, nstreams
c        write(89,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     & (QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,UPIDX),UTA=1,N_OUT_USERTAUS)
c      enddo

c      if ( wf_index.eq.3.and.ibeam.eq.1) then
c        write(99,*)ibeam,
c     &        SURFACEWF_F(WF_INDEX,2,1,IB,UPIDX) 
c      endif

C  Finish

      RETURN
      END

C

      SUBROUTINE DNUSER_SURFACEWF
     I         ( DO_INCLUDE_MVOUTPUT,
     I           FLUX_MULTIPLIER,
     I           IBEAM, WF_INDEX, BRDF_INDEX )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution, multipliers variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of Linearized solution variables

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  BRDF and weighting function index

      INTEGER          BRDF_INDEX, WF_INDEX

C  Beam index

      INTEGER          IBEAM

C  Flux multiplier

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local inclusion flag for mean-value ouptut

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, NL
      INTEGER          UTA, UM, IUM, IQD, I, UT, AA, IB
      DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION L_TOA_SOURCE  (MAX_USER_STREAMS)
      DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM

C  Initial section
C  ---------------

C  Local indices

      IB = IBEAM

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
       DO UTA = 1, N_OUT_USERTAUS
         DO UM = 1, N_USER_STREAMS
           IUM = USEROUTPUT_INDEX(UM)
           SURFACEWF_F(WF_INDEX,UTA,IUM,IB,DNIDX) = ZERO
         ENDDO
       ENDDO
      ENDIF

C  get the TOA source term
C  -----------------------

C  initialise toa source function and recursion

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          L_TOA_SOURCE(UM) = ZERO
          L_CUMUL_SOURCE(UM) = L_TOA_SOURCE(UM)
        ENDDO
      ENDIF

C  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT

C  Thermal transmittance only

            IF ( DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
               L_LAYER_SOURCE(UM) = ZERO
             ENDDO
            ENDIF

C  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO AA = 1, NSTREAMS
                NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XNEG(UM,AA,N)
                PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XPOS(UM,AA,N)
              ENDDO
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                H2 = PCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              L_LAYER_SOURCE(UM) = SHOM
             ENDDO
            ENDIF

C  recursion

            DO UM = LOCAL_UM_START, N_USER_STREAMS
              L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM)
            ENDDO

C  End cumulative layer loop

          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output if flagged

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
           IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
            DO I = 1, NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_UTDN_EIGEN(AA,UT)
                H2 = PCON_ALB(AA,N)*XNEG(I,AA,N)*T_UTUP_EIGEN(AA,UT)
                SHOM = SHOM + H1 + H2
              ENDDO
              QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX) =
     &                       FLUX_MULTIPLIER * SHOM
            ENDDO
           ELSE
            DO I = 1, NSTREAMS
              QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX) = ZERO
            ENDDO
           ENDIF
          ENDIF

C  Copy results if quadrature output required (unlikely)

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
               SURFACEWF_F(WF_INDEX,UTA,IQD,IB,DNIDX) = 
     &                  QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX)
            ENDDO
          ENDIF

C  User-defined stream output

          IF ( DO_USER_STREAMS ) THEN

C  Thermal transmittance only

            IF ( DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
               L_LAYER_SOURCE(UM) = ZERO
             ENDDO
            ENDIF

C  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO AA = 1, NSTREAMS
                NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XNEG(UM,AA,N)
                PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XPOS(UM,AA,N)
              ENDDO
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON_HELP(UM,AA) * UT_HMULT_DD(AA,UM,UT)
                H2 = PCON_HELP(UM,AA) * UT_HMULT_DU(AA,UM,UT)
                SHOM = SHOM + H1 + H2
              ENDDO
              L_LAYER_SOURCE(UM) = SHOM
             ENDDO
            ENDIF

C  assign final albedo weighting function

            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
               SURFACEWF_F(WF_INDEX,UTA,IUM,IB,DNIDX) = 
     &                 FLUX_MULTIPLIER * L_FINAL_SOURCE
            ENDDO

          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output, mean value + flux output

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

            NL = NLEVEL
            N = NL

C  For the highest level

            IF ( NLEVEL .EQ. 0 ) THEN

              DO I = 1, NSTREAMS
                QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX) = ZERO
              ENDDO        

C  For other levels in the atmosphere

            ELSE

              IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO I = 1, NSTREAMS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
                  H2 = PCON_ALB(AA,N)*XNEG(I,AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX) =
     &                          FLUX_MULTIPLIER * SHOM
               ENDDO
              ELSE
                DO I = 1, NSTREAMS
                  QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX) = ZERO
                ENDDO        
              ENDIF

            ENDIF

          ENDIF

C  Copy results if quadrature output required (unlikely)

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
               SURFACEWF_F(WF_INDEX,UTA,IQD,IB,DNIDX) = 
     &                     QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,DNIDX)
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
               SURFACEWF_F(WF_INDEX,UTA,IUM,IB,DNIDX) = 
     &                 FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM)
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  debug

c      do i = 1, nstreams
c        write(56,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     & (QUAD_SURFACEWF(WF_INDEX,UTA,I,IB,dnIDX),UTA=1,N_OUT_USERTAUS)
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE MIFLUX_SURFACEWF (IBEAM,WF_INDEX)

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  arguments

      INTEGER          IBEAM, WF_INDEX

C  local variables

      INTEGER          I, IDIR, WDIR, UTA
      DOUBLE PRECISION SUM_MI, SUM_FX

C  mean intensity and flux
C  -----------------------

C  direction loop

      DO IDIR = 1, N_DIRECTIONS

        WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

        DO UTA = 1, N_OUT_USERTAUS

C  integrations

          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            SUM_MI = SUM_MI +
     &                A(I)  * QUAD_SURFACEWF(WF_INDEX,UTA,I,IBEAM,WDIR)
            SUM_FX = SUM_FX + 
     &                AX(I) * QUAD_SURFACEWF(WF_INDEX,UTA,I,IBEAM,WDIR)
          ENDDO
          MINT_SURFACEWF(WF_INDEX,UTA,IBEAM,WDIR) = SUM_MI * HALF
          FLUX_SURFACEWF(WF_INDEX,UTA,IBEAM,WDIR) = SUM_FX * PI2

C  end loops

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

