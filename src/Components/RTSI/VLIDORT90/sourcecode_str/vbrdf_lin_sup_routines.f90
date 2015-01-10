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
! #              VBRDF_LIN_MAKER                                #
! #              VBRDF_LIN_FOURIER                              #
! #                                                             #
! ###############################################################


      MODULE vbrdf_LinSup_routines_m

      PRIVATE
      PUBLIC :: VBRDF_LIN_MAKER, &
                VBRDF_LIN_FOURIER

      CONTAINS

! 

      SUBROUTINE VBRDF_LIN_MAKER &
         ( BRDF_VFUNCTION_PLUS, BRDF_VFUNCTION_DB_PLUS, &
           DO_EXACTONLY, DO_MSRCORR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
           DO_USER_STREAMS, DO_SURFACE_EMISSION,  n_muquad, n_phiquad, &
           NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
           NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
           QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
           SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
           CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, BRDF_PARS, BRDF_DERIVS, &
           X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
           DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
           BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
           D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
           D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices
!        And their linearizations w.r.t surface parameters

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_VFUNCTION_PLUS
      EXTERNAL         BRDF_VFUNCTION_DB_PLUS

!  Exact only flag (no Fourier term calculations)

      LOGICAL ::          DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      INTEGER ::          MSRCORR_ORDER
      LOGICAL ::          DO_MSRCORR_EXACTONLY

!  Local flags

      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFACE_EMISSION

!  Number of Azimuth waudrature streams

      INTEGER ::          NSTREAMS_BRDF
      INTEGER ::          NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER ::          NSTOKESSQ

!  Local number of Kernel parameters

      INTEGER ::          BRDF_NPARS

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local MSR quadrature control

      INTEGER          :: n_muquad, n_phiquad

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  Local parameter array

      DOUBLE PRECISION :: BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL          :: BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CXE_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Local arrays for MSR quadrature

      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      DOUBLE PRECISION :: EBRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                                     MAXSTREAMS, MAXSTREAMS, &
                                     MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                                     MAXSTREAMS, MAXBEAMS, &
                                     MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: D_USER_BRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_USER_BRDFUNC_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: D_DBKERNEL_BRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      DOUBLE PRECISION :: D_EBRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, &
                       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_USER_EBRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTHALF_BRDF, &
                       MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER          :: I, UI, J, K, KE, IB, ORDER

!  Exact DB calculation
!  --------------------

      DO K = 1, N_USER_RELAZMS
         DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
               if(MSRCORR_ORDER.eq.2)write(*,*)'Doing SZA/VZA = ',IB,UI
               IF ( DO_MSRCORR .or. DO_MSRCORR_EXACTONLY ) THEN
                 CALL BRDF_VFUNCTION_DB_PLUS &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, MSRCORR_ORDER, &
                    NSTOKESSQ, n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),           &
                    USER_STREAMS(UI), USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),      &
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,        &
                    DBKERNEL_BRDFUNC(:,UI,K,IB), D_DBKERNEL_BRDFUNC(:,:,UI,K,IB) )
               ELSE
                 CALL BRDF_VFUNCTION_PLUS &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,   &
                    NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                    USER_SINES(UI), PHIANG(K),  COSPHI(K), SINPHI(K),          &
                    DBKERNEL_BRDFUNC(:,UI,K,IB), D_DBKERNEL_BRDFUNC(:,:,UI,K,IB) )
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!  Return if this is all you require

      IF ( DO_EXACTONLY ) RETURN

!  Local

      ORDER = MSRCORR_ORDER

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
         IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
          CALL BRDF_VFUNCTION_DB_PLUS &
            ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,    &
              NSTOKESSQ, n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),      &
              QUAD_STREAMS(I), QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
              X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,   &
              BRDFUNC_0(:,I,IB,K), D_BRDFUNC_0(:,:,I,IB,K) )
         ELSE
          CALL BRDF_VFUNCTION_PLUS &
            ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,  &
              NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I), &
              QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         &
              BRDFUNC_0(:,I,IB,K), D_BRDFUNC_0(:,:,I,IB,K) )
         ENDIF
        ENDDO
       ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
           IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
            CALL BRDF_VFUNCTION_DB_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,     &
                 NSTOKESSQ, n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J),     &
                 QUAD_STREAMS(I), QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),  &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,    &
                 BRDFUNC(:,I,J,K), D_BRDFUNC(:,:,I,J,K)  )
           ELSE
            CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,    &
                 NSTOKESSQ, QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I), &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),           &
                 BRDFUNC(1,I,J,K), D_BRDFUNC(:,:,I,J,K)  )
           ENDIF
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
             IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
              CALL BRDF_VFUNCTION_DB_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,    &
                   NSTOKESSQ, n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE),        &
                   QUAD_STREAMS(I), QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,   &
                   EBRDFUNC(:,I,KE,K), D_EBRDFUNC(:,:,I,KE,K)  )
             ELSE
              CALL BRDF_VFUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, &
                   NSTOKESSQ, CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),   &
                   QUAD_SINES(I), X_BRDF(K),  CX_BRDF(K), SX_BRDF(K),        &
                   EBRDFUNC(:,I,KE,K), D_EBRDFUNC(:,:,I,KE,K) )
             ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
         DO UI = 1, N_USER_STREAMS
          DO K = 1, NSTREAMS_BRDF
           IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
            CALL BRDF_VFUNCTION_DB_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,      &
                 NSTOKESSQ, n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),        &
                 USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                 USER_BRDFUNC_0(:,UI,IB,K), D_USER_BRDFUNC_0(:,:,UI,IB,K) )
           ELSE
            CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,   &
                 NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         &
                 USER_BRDFUNC_0(:,UI,IB,K), D_USER_BRDFUNC_0(:,:,UI,IB,K) )
           ENDIF
          ENDDO
         ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
             IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
              CALL BRDF_VFUNCTION_DB_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,     &
                   NSTOKESSQ, n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J),      &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_BRDFUNC(:,UI,J,K), D_USER_BRDFUNC(:,:,UI,J,K) )
             ELSE
              CALL BRDF_VFUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,     &
                   NSTOKESSQ, QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI), &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),           &
                   USER_BRDFUNC(:,UI,J,K), D_USER_BRDFUNC(:,:,UI,J,K) )
             ENDIF
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_EXACTONLY ) THEN
                CALL BRDF_VFUNCTION_DB_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, ORDER,      &
                   NSTOKESSQ, n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE),          &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_EBRDFUNC(:,UI,KE,K), D_USER_EBRDFUNC(:,:,UI,KE,K) )
               ELSE
                CALL BRDF_VFUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, &
                   NSTOKESSQ, CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI), &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                   USER_EBRDFUNC(1,UI,KE,K), D_USER_EBRDFUNC(:,:,UI,KE,K) )
               ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_LIN_MAKER

!

      SUBROUTINE VBRDF_LIN_FOURIER &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION, LAMBERTIAN_FLAG, M,   &
           NSTOKES, NSTOKESSQ, NBEAMS, NSTREAMS, N_USER_STREAMS,       &
           NSTREAMS_BRDF, NBRDF_HALF, BRDF_NPARS, BRDF_DERIVS, DELFAC, &
           FACTOR, BRDF_COSAZMFAC, BRDF_SINAZMFAC, A_BRDF, BAX_BRDF,   &
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,     &
           D_EBRDFUNC, D_USER_EBRDFUNC, D_LOCAL_BRDF_F,                  &
           D_LOCAL_BRDF_F_0, D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0, &
           D_LOCAL_EMISSIVITY, D_LOCAL_USER_EMISSIVITY )

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Prepares linearizations of Fourier component of the bidirectional reflectance functions

!  Input arguments
!  ===============

!  Control

      LOGICAL ::          LAMBERTIAN_FLAG
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFACE_EMISSION

!  Local numbers

      INTEGER ::          M, NSTOKES, NSTOKESSQ
      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

!  linearization Control

      INTEGER ::          BRDF_NPARS
      LOGICAL ::          BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Surface factors

      DOUBLE PRECISION :: DELFAC, FACTOR

!  Azimuth cosines/sines and weights

      DOUBLE PRECISION :: BRDF_COSAZMFAC ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDF_SINAZMFAC ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: A_BRDF         ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BAX_BRDF       ( MAXSTHALF_BRDF  )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                                     MAXSTREAMS, MAXSTREAMS, &
                                     MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                                     MAXSTREAMS, MAXBEAMS, &
                                     MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: D_USER_BRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_USER_BRDFUNC_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Values for Emissivity

      DOUBLE PRECISION :: D_EBRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, &
                       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: D_USER_EBRDFUNC &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTHALF_BRDF, &
                       MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: D_LOCAL_BRDF_F &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_BRDF_F_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      DOUBLE PRECISION :: D_LOCAL_USER_BRDF_F &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_USER_BRDF_F_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      DOUBLE PRECISION :: D_LOCAL_EMISSIVITY &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_USER_EMISSIVITY &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES, &
                       MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER ::          I, UI, J, K, KPHI, IB, Q, O1, O2, W
      DOUBLE PRECISION :: SUM, REFL, HELP, EMISS(16)
      INTEGER ::          COSSIN_MASK(16)

      COSSIN_MASK = (/ 1,1,2,0,1,1,2,0,2,2,1,0,0,0,0,1 /)

!  surface factor

      HELP = HALF * DELFAC

!  Zeroing

      D_LOCAL_BRDF_F        = ZERO
      D_LOCAL_BRDF_F_0      = ZERO
      D_LOCAL_USER_BRDF_F   = ZERO
      D_LOCAL_USER_BRDF_F_0 = ZERO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO W = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(W) ) THEN
            DO IB = 1, NBEAMS
              DO I = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  SUM = ZERO
                  IF ( COSSIN_MASK(Q) .EQ. 1 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC_0(W,Q,I,IB,K) * &
                             BRDF_COSAZMFAC(K)
                    ENDDO
                  ELSE IF ( COSSIN_MASK(Q) .EQ. 2 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC_0(W,Q,I,IB,K) * &
                             BRDF_SINAZMFAC(K)
                    ENDDO
                  ENDIF
                  D_LOCAL_BRDF_F_0(W,Q,I,IB) = SUM * HELP
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO W = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(W) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  SUM = ZERO
                  IF ( COSSIN_MASK(Q) .EQ. 1 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC(W,Q,I,J,K) * &
                             BRDF_COSAZMFAC(K)
                    ENDDO
                  ELSE IF ( COSSIN_MASK(Q) .EQ. 2 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC(W,Q,I,J,K) * &
                             BRDF_SINAZMFAC(K)
                    ENDDO
                  ENDIF
                  D_LOCAL_BRDF_F(W,Q,I,J) = SUM * HELP
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  DO Q = 1, NSTOKESSQ
                    SUM = ZERO
                    IF ( COSSIN_MASK(Q) .EQ. 1 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC_0(W,Q,UI,IB,K) * &
                               BRDF_COSAZMFAC(K)
                      ENDDO
                    ELSE IF ( COSSIN_MASK(Q) .EQ. 2 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC_0(W,Q,UI,IB,K) * &
                               BRDF_SINAZMFAC(K)
                      ENDDO
                    ENDIF
                    D_LOCAL_USER_BRDF_F_0(W,Q,UI,IB) = SUM * HELP
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO UI = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  DO Q = 1, NSTOKESSQ
                    SUM = ZERO
                    IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC(W,Q,UI,J,K) * &
                               BRDF_COSAZMFAC(K)
                      ENDDO
                    ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC(W,Q,UI,J,K) * &
                               BRDF_SINAZMFAC(K)
                      ENDDO
                    ENDIF
                    D_LOCAL_USER_BRDF_F(W,Q,UI,J) = SUM * HELP
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Azimuth independent contribution, from Kirchhoff's law
!  Assumed to exist only for the total intensity
!        (first element of Stokes Vector) - is this right ??????

      IF ( DO_SURFACE_EMISSION .AND. M .EQ. 0 ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG .AND. M .EQ. 0 ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO I = 1, NSTREAMS
                D_LOCAL_EMISSIVITY(W,1,I) = ZERO
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  D_LOCAL_USER_EMISSIVITY(W,1,UI) = ZERO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!  bidirectional reflectance

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN

!  Inserted Polarization sum here.   Still to be checked.....!!!!!!!

!  Quadrature polar directions

!mick fix
EMISS = ZERO

          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  REFL = ZERO
                  DO KPHI= 1, NSTREAMS_BRDF
                    SUM = ZERO
                    DO K = 1, NBRDF_HALF
                      SUM = SUM + D_EBRDFUNC(W,Q,I,K,KPHI) * BAX_BRDF(K)
                    ENDDO
                    REFL = REFL + A_BRDF(KPHI) * SUM
                  ENDDO
                  EMISS(Q) = REFL
                ENDDO
                DO O1 = 1, NSTOKES
                  REFL = ZERO
                  DO O2 = 1, NSTOKES
                    Q = 4 * ( O1 - 1 ) + O2
                    REFL = REFL + EMISS(Q)
                  ENDDO
                  D_LOCAL_EMISSIVITY(W,O1,I) = REFL * FACTOR
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO W = 1, BRDF_NPARS
              IF ( BRDF_DERIVS(W) ) THEN
                DO UI = 1, N_USER_STREAMS
                  DO Q = 1, NSTOKESSQ
                    REFL = ZERO
                    DO KPHI= 1, NSTREAMS_BRDF
                      SUM = ZERO
                      DO K = 1, NBRDF_HALF
                        SUM = SUM + D_USER_EBRDFUNC(W,Q,UI,K,KPHI) * &
                              BAX_BRDF(K)
                      ENDDO
                      REFL = REFL + A_BRDF(KPHI) * SUM
                    ENDDO
                    EMISS(Q) = REFL
                  ENDDO
                  DO O1 = 1, NSTOKES
                    REFL = ZERO
                    DO O2 = 1, NSTOKES
                      Q = 4 * ( O1 - 1 ) + O2
                      REFL = REFL + EMISS(Q)
                    ENDDO
                    D_LOCAL_USER_EMISSIVITY(W,O1,UI) = REFL * FACTOR
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Not lambertian

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_LIN_FOURIER

!  End module

      END MODULE vbrdf_LinSup_routines_m

