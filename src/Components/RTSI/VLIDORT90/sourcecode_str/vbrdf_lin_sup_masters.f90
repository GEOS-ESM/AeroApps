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
! #            VBRDF_LIN_INPUTMASTER                            #
! #                                                             #
! #            VBRDF_LIN_MAINMASTER (master), calling           #
! #              VBRDF_LIN_MAKER                                #
! #              BRDF_QUADRATURE_Gaussian                       #
! #              BRDF_QUADRATURE_Trapezoid (not used)           #
! #              VBRDF_LIN_FOURIER                              #
! #                                                             #
! ###############################################################


MODULE vbrdf_LinSup_masters_m

      PRIVATE
      PUBLIC :: VBRDF_LIN_INPUTMASTER, &
                VBRDF_LIN_MAINMASTER

      CONTAINS

      SUBROUTINE VBRDF_LIN_MAINMASTER (  &
        DO_DEBUG_RESTORATION, & ! Inputs
        NMOMENTS_INPUT,       & ! Inputs
        VBRDF_Sup_In,         & ! Inputs
        VBRDF_LinSup_In,      & ! Inputs
        VBRDF_Sup_Out,        & ! Outputs
        VBRDF_LinSup_Out )      ! Outputs

!  Prepares the bidirectional reflectance functions
!  necessary for VLIDORT.

      USE VLIDORT_PARS

      USE vbrdf_sup_inputs_def
      USE vbrdf_linsup_inputs_def

      USE vbrdf_sup_outputs_def
      USE vbrdf_linsup_outputs_def

      USE vbrdf_sup_aux_m, only : BRDF_GAULEG,              &
                                  BRDF_QUADRATURE_Gaussian, &
                                  BRDF_QUADRATURE_Trapezoid

      USE vbrdf_sup_kernels_m
      USE vbrdf_linsup_kernels_m

      USE vbrdf_sup_routines_m
      USE vbrdf_linsup_routines_m

      IMPLICIT NONE

!  Inputs
!  ------

!  Debug flag for restoration

      LOGICAL, INTENT(IN) ::                 DO_DEBUG_RESTORATION

!  Input number of moments (only used for restoration debug)

      INTEGER, INTENT(IN) ::                 NMOMENTS_INPUT

!  Input structure
!  ---------------

      TYPE(VBRDF_Sup_Inputs)    , INTENT(IN)  :: VBRDF_Sup_In
      TYPE(VBRDF_LinSup_Inputs) , INTENT(IN)  :: VBRDF_LinSup_In

!  Output structure
!  ----------------

      TYPE(VBRDF_Sup_Outputs)   , INTENT(OUT) :: VBRDF_Sup_Out
      TYPE(VBRDF_LinSup_Outputs), INTENT(OUT) :: VBRDF_LinSup_Out

!  VLIDORT local variables
!  ++++++++++++++++++++++

!  Input arguments
!  ===============

!  User stream Control

      LOGICAL ::          DO_USER_STREAMS

!  Surface emission

      LOGICAL ::          DO_SURFACE_EMISSION

!  number of Stokes components

      INTEGER ::          NSTOKES

!   Number and index-list of bidirectional functions

      INTEGER ::          N_BRDF_KERNELS
      INTEGER ::          WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER ::          N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      DOUBLE PRECISION :: BRDF_PARAMETERS &
          ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF names

      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Lambertian Surface control

      LOGICAL ::          LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      DOUBLE PRECISION :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER ::          NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL ::          DO_SHADOW_EFFECT

!  Exact only flag (no Fourier term calculations)

      LOGICAL ::          DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      INTEGER ::          MSRCORR_ORDER
      LOGICAL ::          DO_MSRCORR_EXACTONLY

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL ::          DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL ::          DO_KERNEL_PARAMS_WFS &
          ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL ::          DO_KPARAMS_DERIVS ( MAX_BRDF_KERNELS )

!  number of surfaceweighting functions

      INTEGER ::          N_SURFACE_WFS
      INTEGER ::          N_KERNEL_FACTOR_WFS
      INTEGER ::          N_KERNEL_PARAMS_WFS

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Angles

      DOUBLE PRECISION :: BEAM_SZAS (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  BRDF External functions
!  =======================

!  lambertian

!      EXTERNAL         LAMBERTIAN_VFUNCTION

!  Modis-type kernels

!      EXTERNAL         ROSSTHIN_VFUNCTION
!      EXTERNAL         ROSSTHICK_VFUNCTION
!      EXTERNAL         LISPARSE_VFUNCTION
!      EXTERNAL         LIDENSE_VFUNCTION
!      EXTERNAL         HAPKE_VFUNCTION
!      EXTERNAL         ROUJEAN_VFUNCTION
!      EXTERNAL         RAHMAN_VFUNCTION

!  Cox-munk types

!      EXTERNAL         COXMUNK_VFUNCTION
!      EXTERNAL         COXMUNK_VFUNCTION_DB
!      EXTERNAL         GISSCOXMUNK_VFUNCTION
!      EXTERNAL         GISSCOXMUNK_VFUNCTION_DB

!  GCM CRI is not an external call
!      EXTERNAL         GCMCRI_VFUNCTION
!      EXTERNAL         GCMCRI_VFUNCTION_DB

!  new for Version 2.4R, introduced 30 April 2009, 6 May 2009
!    2009 function is final Kernel supplied by Breon, May 5 2009.

!      EXTERNAL         BPDF2009_VFUNCTION

!  Local BRDF functions
!  ====================

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

!  DB Kernel values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      DOUBLE PRECISION :: EBRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF)
      DOUBLE PRECISION :: USER_EBRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF)

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

!  Linearized Exact DB values

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

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_WEIGHTS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_STRMWTS(MAXSTREAMS)

!  Viewing zenith streams

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER ::          NBRDF_HALF
      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      DOUBLE PRECISION :: BAX_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: CXE_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Azimuth factors

      DOUBLE PRECISION :: BRDF_COSAZMFAC(MAXSTREAMS_BRDF)
      DOUBLE PRECISION :: BRDF_SINAZMFAC(MAXSTREAMS_BRDF)

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: LOCAL_BRDF_F &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_BRDF_F_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      DOUBLE PRECISION :: LOCAL_USER_BRDF_F &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_USER_BRDF_F_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      DOUBLE PRECISION :: LOCAL_EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Local Derivative-kernel Fourier components
!  ==========================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: D_LOCAL_BRDF_F &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, &
                       MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_BRDF_F_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, &
                       MAXBEAMS )

!  at user-defined stream directions

      DOUBLE PRECISION :: D_LOCAL_USER_BRDF_F &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_USER_BRDF_F_0 &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, &
                       MAX_USER_STREAMS, MAXBEAMS )

!  emissivities

      DOUBLE PRECISION :: D_LOCAL_EMISSIVITY &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: D_LOCAL_USER_EMISSIVITY &
                     ( MAX_BRDF_PARAMETERS, MAXSTOKES, &
                       MAX_USER_STREAMS )

!  Other local variables
!  =====================

!  Spherical albedo

      DOUBLE PRECISION :: SPHERICAL_ALBEDO (MAX_BRDF_KERNELS)

!  help

      INTEGER          :: WOFFSET ( MAX_BRDF_KERNELS)
      INTEGER          :: K, B, I, I1, J, IB, UI, UM, IA, M, O1, Q, P, W
      INTEGER          :: BRDF_NPARS, NMOMENTS, NSTOKESSQ, N_phiquad_HALF
      DOUBLE PRECISION :: PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL          :: DERIVS ( MAX_BRDF_PARAMETERS )
      DOUBLE PRECISION :: MUX, DELFAC, HELP_A, SUM, ARGUMENT, XM, FF
      LOGICAL          :: ADD_FOURIER, LOCAL_MSR

!  default, use Gaussian quadrature

      LOGICAL, PARAMETER :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Copy from input structure
!  -------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = VBRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = VBRDF_Sup_In%BS_DO_BRDF_SURFACE
      DO_SURFACE_EMISSION = VBRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  Copy Geometry results

      NSTOKES           = VBRDF_Sup_In%BS_NSTOKES
      NSTREAMS          = VBRDF_Sup_In%BS_NSTREAMS
      NBEAMS            = VBRDF_Sup_In%BS_NBEAMS
      BEAM_SZAS         = VBRDF_Sup_In%BS_BEAM_SZAS
      N_USER_RELAZMS    = VBRDF_Sup_In%BS_N_USER_RELAZMS
      USER_RELAZMS      = VBRDF_Sup_In%BS_USER_RELAZMS
      N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
      USER_ANGLES_INPUT = VBRDF_Sup_In%BS_USER_ANGLES_INPUT

!  Copy BRDF inputs

      N_BRDF_KERNELS         = VBRDF_Sup_In%BS_N_BRDF_KERNELS
      BRDF_NAMES             = VBRDF_Sup_In%BS_BRDF_NAMES
      WHICH_BRDF             = VBRDF_Sup_In%BS_WHICH_BRDF
      N_BRDF_PARAMETERS      = VBRDF_Sup_In%BS_N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = VBRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG
      BRDF_FACTORS           = VBRDF_Sup_In%BS_BRDF_FACTORS
      NSTREAMS_BRDF          = VBRDF_Sup_In%BS_NSTREAMS_BRDF
      DO_SHADOW_EFFECT       = VBRDF_Sup_In%BS_DO_SHADOW_EFFECT
      DO_EXACTONLY           = VBRDF_Sup_In%BS_DO_EXACTONLY

!  Copy linearized BRDF inputs

      DO_KERNEL_FACTOR_WFS   = VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS
      DO_KERNEL_PARAMS_WFS   = VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS
      DO_KPARAMS_DERIVS      = VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS
      N_SURFACE_WFS          = VBRDF_LinSup_In%BS_N_SURFACE_WFS
      N_KERNEL_FACTOR_WFS    = VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS
      N_KERNEL_PARAMS_WFS    = VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS

!  Copy MSR inputs

      DO_MSRCORR             = VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR
      DO_MSRCORR_EXACTONLY   = VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_EXACTONLY 
      MSRCORR_ORDER          = VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER
      n_muquad               = VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD
      n_Phiquad              = VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD

!  Main code
!  ---------

!  Set up Quadrature streams

      CALL BRDF_GAULEG ( 0.0d0, 1.0d0, QUAD_STREAMS, QUAD_WEIGHTS, NSTREAMS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = DSQRT(1.0d0-QUAD_STREAMS(I)*QUAD_STREAMS(I))
        QUAD_STRMWTS(I) = QUAD_STREAMS(I) * QUAD_WEIGHTS(I)
      enddo

!  Number of Stokes components squared
!    ** Bookkeeping for surface kernel Cox-Munk types
!    ** Only the Giss CoxMunk kernel is vectorized (as of 19 January 200

!   Additional code for complex RI Giss Cox-Munk, 15 march 2010.
!     3 parameters are PARS(1) = sigma_sq
!                      PARS(2) = Real (RI)
!                      PARS(3) = Imag (RI)

      NSTOKESSQ  = 1
      DO K = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(K) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(K) .EQ. 'GissCoxMnk' ) THEN
            N_BRDF_PARAMETERS(K) = 3
            IF ( DO_SHADOW_EFFECT ) THEN
              BRDF_PARAMETERS(K,3) = ONE
            ELSE
              BRDF_PARAMETERS(K,3) = ZERO
            ENDIF
         ELSE IF ( BRDF_NAMES(K) .EQ. 'GCMcomplex' ) THEN
            N_BRDF_PARAMETERS(K) = 3
         ENDIF
          IF ( BRDF_NAMES(K) .EQ. 'GissCoxMnk'  .OR. &
               BRDF_NAMES(K) .EQ. 'GCMcomplex' ) THEN
            NSTOKESSQ = NSTOKES * NSTOKES
         ENDIF
      ENDDO

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams
!    Warning, this shoudl be the BOA angle. OK for the non-refractive ca
!
      DO IB = 1, NBEAMS
        MUX =  DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
        SZASURCOS(IB) = MUX
        SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
      ENDDO

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = DCOS(USER_ANGLES_INPUT(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = DSQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
        COSPHI(IA) = DCOS(PHIANG(IA))
        SINPHI(IA) = DSIN(PHIANG(IA))
      ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian &
           ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF, &
             X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
             BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid &
           ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF, &
             X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
             BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ENDIF

!  Number of weighting functions, and offset

      W = 0
      WOFFSET(1) = 0
      DO K = 1, N_BRDF_KERNELS
        IF ( DO_KERNEL_FACTOR_WFS(K) ) W = W + 1
        DO P = 1, N_BRDF_PARAMETERS(K)
          IF ( DO_KERNEL_PARAMS_WFS(K,P) ) W = W + 1
        ENDDO
        IF ( K.LT.N_BRDF_KERNELS ) WOFFSET(K+1) = W
      ENDDO

      N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS
      IF ( W .ne. N_SURFACE_WFS ) stop'bookkeeping wrong'

!  Set up the MSR points
!  ---------------------

!  Air to water, Polar quadrature

      if ( DO_MSRCORR  ) THEN
         CALL brdf_gauleg ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
         DO I = 1, N_MUQUAD
            XM = X_MUQUAD(I)
            SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
            WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
         ENDDO
      endif

!  Azimuth quadrature

      if ( DO_MSRCORR  ) THEN
         N_phiquad_HALF = N_PHIQUAD / 2
         CALL brdf_gauleg ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
         DO I = 1, N_PHIQUAD_HALF
           I1 = I + N_PHIQUAD_HALF
           X_PHIQUAD(I1) = - X_PHIQUAD(I)
           W_PHIQUAD(I1) =   W_PHIQUAD(I)
         ENDDO
         DO I = 1, N_PHIQUAD
            X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
         ENDDO
      ENDIF

!  Initialise ALL outputs
!  ----------------------

!  Zero the BRDF Fourier components

      VBRDF_Sup_Out%BS_BRDF_F_0 = ZERO
      VBRDF_Sup_Out%BS_BRDF_F   = ZERO
      VBRDF_Sup_Out%BS_USER_BRDF_F_0 = ZERO
      VBRDF_Sup_Out%BS_USER_BRDF_F   = ZERO

!  Zero Exact Direct Beam BRDF

      VBRDF_Sup_Out%BS_EXACTDB_BRDFUNC = ZERO

!  Initialize surface emissivity

      VBRDF_Sup_Out%BS_EMISSIVITY      = ONE
      VBRDF_Sup_Out%BS_USER_EMISSIVITY = ONE

!  initialize linearized quantities

      VBRDF_LinSup_Out%BS_LS_BRDF_F_0      = ZERO
      VBRDF_LinSup_Out%BS_LS_BRDF_F        = ZERO
      VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0 = ZERO
      VBRDF_LinSup_Out%BS_LS_USER_BRDF_F   = ZERO

      VBRDF_LinSup_Out%BS_LS_EXACTDB_BRDFUNC = ZERO

      VBRDF_LinSup_Out%BS_LS_USER_EMISSIVITY = ONE
      VBRDF_LinSup_Out%BS_LS_EMISSIVITY      = ONE

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Local variables

        BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, MAX_BRDF_PARAMETERS
          PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, MAX_BRDF_PARAMETERS
            DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ENDIF

!  Local MSRCORR flag

        LOCAL_MSR = .false.
        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX     .or. &
             WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX .or. &
             WHICH_BRDF(K) .EQ. GISSCOXMUNK_CRI_IDX ) THEN
           LOCAL_MSR = DO_MSRCORR
        ENDIF

!  Lambertian kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. LAMBERTIAN_IDX ) THEN
          CALL VBRDF_MAKER &
             ( LAMBERTIAN_VFUNCTION, LAMBERTIAN_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL VBRDF_MAKER &
             ( ROSSTHIN_VFUNCTION, ROSSTHIN_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL VBRDF_MAKER &
             ( ROSSTHICK_VFUNCTION, ROSSTHICK_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
             ( LISPARSE_VFUNCTION_PLUS, LISPARSE_VFUNCTION_PLUS, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,  X_BRDF, &
               CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
               D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
               D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
               D_USER_EBRDFUNC )
          ELSE
            CALL VBRDF_MAKER &
             ( LISPARSE_VFUNCTION, LISPARSE_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
             ( LIDENSE_VFUNCTION_PLUS, LIDENSE_VFUNCTION_PLUS, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION,  n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
               CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
               D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
               D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
               D_USER_EBRDFUNC )
          ELSE
            CALL VBRDF_MAKER &
             ( LIDENSE_VFUNCTION, LIDENSE_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
             ( HAPKE_VFUNCTION_PLUS, HAPKE_VFUNCTION_PLUS, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION,  n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
               CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
               D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
               D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
               D_USER_EBRDFUNC )
          ELSE
            CALL VBRDF_MAKER &
             ( HAPKE_VFUNCTION, HAPKE_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL VBRDF_MAKER &
             ( ROUJEAN_VFUNCTION, ROUJEAN_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  Rahman kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
             ( RAHMAN_VFUNCTION_PLUS, RAHMAN_VFUNCTION_PLUS, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION,  n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
               CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
               D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
               D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
               D_USER_EBRDFUNC )
          ELSE
            CALL VBRDF_MAKER &
             ( RAHMAN_VFUNCTION, RAHMAN_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  Scalar-only original Cox-Munk kernel: (2 free parameters, Shadow = Th
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
               ( COXMUNK_VFUNCTION_PLUS, COXMUNK_VFUNCTION_DB_PLUS, &
                 DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
                 DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
                 NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
                 NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
                 QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
                 SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
                 CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
                 DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
                 BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
                 D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
                 D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
                 D_USER_EBRDFUNC )
         ELSE
           CALL VBRDF_MAKER &
               ( COXMUNK_VFUNCTION, COXMUNK_VFUNCTION_DB, &
                 DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
                 DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
                 NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
                 NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
                 QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
                 SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
                 X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
                 DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
                 BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  GISS Vector Cox-Munk kernel: (2 free parameters, Shadow = Third). Rea
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VBRDF_LIN_MAKER &
               ( GISSCOXMUNK_VFUNCTION_PLUS, GISSCOXMUNK_VFUNCTION_DB_PLUS, &
                 DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
                 DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
                 NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
                 NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
                 QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
                 SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, &
                 CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS, &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
                 DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
                 BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
                 D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, &
                 D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, &
                 D_USER_EBRDFUNC )
         ELSE
           CALL VBRDF_MAKER &
               ( GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION_DB, &
                 DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
                 DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
                 NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
                 NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
                 QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
                 SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
                 X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
                 DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
                 BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
          ENDIF
        ENDIF

!  New code . Giss Cox_munk with Complex RI : (3 Free parameters). Shadow
!   NO LINEARIZATION with this Kernel

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_CRI_IDX ) THEN
          CALL VBRDF_GCMCRI_MAKER &
             ( DO_USER_STREAMS, DO_SURFACE_EMISSION, DO_SHADOW_EFFECT, & 
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, n_muquad, n_phiquad, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  BPDF 2009 kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2009_IDX ) THEN
          CALL VBRDF_MAKER &
             ( BPDF2009_VFUNCTION, BPDF2009_VFUNCTION, &
               DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY, MSRCORR_ORDER, &
               DO_USER_STREAMS, DO_SURFACE_EMISSION, n_muquad, n_phiquad, &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS, &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES, &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS, &
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, &
               X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, &
               BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )
        ENDIF

!  Exact BRDFUNC
!  ------------

!  factor

        FF = BRDF_FACTORS(K)

!  Compute Exact Direct Beam BRDF

        DO O1 = 1, NSTOKESSQ
          DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                VBRDF_Sup_Out%BS_EXACTDB_BRDFUNC(O1,UM,IA,IB) = &
                VBRDF_Sup_Out%BS_EXACTDB_BRDFUNC(O1,UM,IA,IB) &
                + FF * DBKERNEL_BRDFUNC(O1,UM,IA,IB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Linearization w.r.t Kernel Factor

        W  = WOFFSET(K)
        IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
          W = W + 1
          DO O1 = 1, NSTOKESSQ
            DO IA = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  VBRDF_LinSup_Out%BS_LS_EXACTDB_BRDFUNC(W,O1,UM,IA,IB) = &
                      DBKERNEL_BRDFUNC(O1,UM,IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Linearization w.r.t Kernel parameters

        DO P = 1, BRDF_NPARS
          IF ( DERIVS(P) ) THEN
            W = W + 1
            DO O1 = 1, NSTOKESSQ
              DO IA = 1, N_USER_RELAZMS
                DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   VBRDF_LinSup_Out%BS_LS_EXACTDB_BRDFUNC(W,O1,UM,IA,IB) = &
                    FF * D_DBKERNEL_BRDFUNC(P,O1,UM,IA,IB)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) .OR. &
                          (LAMBERTIAN_KERNEL_FLAG(K) .AND. M .EQ. 0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = ONE
            DO I = 1, NSTREAMS_BRDF
              BRDF_COSAZMFAC(I) = A_BRDF(I)
              BRDF_SINAZMFAC(I) = ZERO
            ENDDO
          ELSE
            DELFAC   = TWO
            DO I = 1, NSTREAMS_BRDF
              ARGUMENT = DBLE(M) * X_BRDF(I)
              BRDF_COSAZMFAC(I) = A_BRDF(I) * DCOS ( ARGUMENT )
              BRDF_SINAZMFAC(I) = A_BRDF(I) * DSIN ( ARGUMENT )
            ENDDO
          ENDIF

!  Call

          CALL VBRDF_FOURIER &
             ( DO_USER_STREAMS, DO_SURFACE_EMISSION, &
               LAMBERTIAN_KERNEL_FLAG(K), M, NSTOKES, NSTOKESSQ, NBEAMS, &
               NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, &
               DELFAC, BRDF_FACTORS(K), BRDF_COSAZMFAC, BRDF_SINAZMFAC, &
               A_BRDF, BAX_BRDF, BRDFUNC, USER_BRDFUNC, BRDFUNC_0, &
               USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, &
               LOCAL_BRDF_F, LOCAL_BRDF_F_0, LOCAL_USER_BRDF_F, &
               LOCAL_USER_BRDF_F_0, LOCAL_EMISSIVITY, &
               LOCAL_USER_EMISSIVITY )

!  Linear call

          IF ( BRDF_NPARS .GT. 0 ) THEN
            CALL VBRDF_LIN_FOURIER &
               ( DO_USER_STREAMS, DO_SURFACE_EMISSION, &
                 LAMBERTIAN_KERNEL_FLAG(K), M, NSTOKES, NSTOKESSQ, &
                 NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, &
                 NBRDF_HALF, BRDF_NPARS, DERIVS, DELFAC, &
                 BRDF_FACTORS(K), BRDF_COSAZMFAC, BRDF_SINAZMFAC, &
                 A_BRDF, BAX_BRDF, D_BRDFUNC, D_USER_BRDFUNC, &
                 D_BRDFUNC_0, D_USER_BRDFUNC_0, &
                 D_EBRDFUNC, D_USER_EBRDFUNC, &
                 D_LOCAL_BRDF_F, D_LOCAL_BRDF_F_0, &
                 D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0, &
                 D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )
          ENDIF

!  Spherical albedo. debug only

          IF ( M .EQ. 0 ) THEN
            Q = 1
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, NSTREAMS
               SUM = ZERO
               DO J = 1, NSTREAMS
                 SUM = SUM + LOCAL_BRDF_F(Q,I,J) * QUAD_STRMWTS(J)
               ENDDO
               HELP_A = HELP_A + SUM * QUAD_STRMWTS(I)
              ENDDO
              SPHERICAL_ALBEDO(K) = HELP_A * FOUR
             ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)
!  ------------------------------------------------

!  factor

            FF = BRDF_FACTORS(K)

!  Basic Kernel sum

            DO Q = 1, NSTOKESSQ
              DO I = 1, NSTREAMS
                DO IB = 1, NBEAMS
                  VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,I,IB) = VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,I,IB) &
                  + BRDF_FACTORS(K) * LOCAL_BRDF_F_0(Q,I,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  VBRDF_Sup_Out%BS_BRDF_F(M,Q,I,J) = VBRDF_Sup_Out%BS_BRDF_F(M,Q,I,J) &
                  + BRDF_FACTORS(K) * LOCAL_BRDF_F(Q,I,J)
                ENDDO
              ENDDO
            ENDDO

!  Linearization w.r.t Kernel Factor

            W  = WOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              W = W + 1
              DO Q = 1, NSTOKESSQ
                DO I = 1, NSTREAMS
                  DO IB = 1, NBEAMS
                    VBRDF_LinSup_Out%BS_LS_BRDF_F_0(W,M,Q,I,IB) = LOCAL_BRDF_F_0(Q,I,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    VBRDF_LinSup_Out%BS_LS_BRDF_F(W,M,Q,I,J) =  LOCAL_BRDF_F(Q,I,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, BRDF_NPARS
              IF ( DERIVS(P) ) THEN
                W = W + 1
                DO Q = 1, NSTOKESSQ
                  DO I = 1, NSTREAMS
                    DO IB = 1, NBEAMS
                      VBRDF_LinSup_Out%BS_LS_BRDF_F_0(W,M,Q,I,IB) = FF*D_LOCAL_BRDF_F_0(P,Q,I,IB)
                    ENDDO
                    DO J = 1, NSTREAMS
                      VBRDF_LinSup_Out%BS_LS_BRDF_F(W,M,Q,I,J) = FF*D_LOCAL_BRDF_F(P,Q,I,J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Kernel combinations (for user-stream reflectance)
!  -------------------------------------------------

!  Basic kernel summation

            IF ( DO_USER_STREAMS ) THEN
              DO Q = 1, NSTOKESSQ
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,UM,IB) = VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,UM,IB) &
                    + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F_0(Q,UM,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,UM,J) = VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,UM,J) &
                    + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F(Q,UM,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel Factor

            W  = WOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              W = W + 1
              DO Q = 1, NSTOKESSQ
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0(W,M,Q,UM,IB) = &
                        LOCAL_USER_BRDF_F_0(Q,UM,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    VBRDF_LinSup_Out%BS_LS_USER_BRDF_F(W,M,Q,UM,J) = &
                        LOCAL_USER_BRDF_F(Q,UM,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, BRDF_NPARS
              IF ( DERIVS(P) ) THEN
                W = W + 1
                DO Q = 1, NSTOKESSQ
                  DO UM = 1, N_USER_STREAMS
                    DO IB = 1, NBEAMS
                      VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0(W,M,Q,UM,IB) = &
                          FF * D_LOCAL_USER_BRDF_F_0(P,Q,UM,IB)
                    ENDDO
                    DO J = 1, NSTREAMS
                      VBRDF_LinSup_Out%BS_LS_USER_BRDF_F(W,M,Q,UM,J) = &
                          FF * D_LOCAL_USER_BRDF_F(P,Q,UM,J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Total emissivities
!  ------------------

!  only if flagged

            IF ( DO_SURFACE_EMISSION .AND.  M .EQ. 0 ) THEN

!  Basci kernel contributions

              DO Q = 1, NSTOKES
                DO I = 1, NSTREAMS
                  VBRDF_Sup_Out%BS_EMISSIVITY(Q,I) = &
                  VBRDF_Sup_Out%BS_EMISSIVITY(Q,I) - LOCAL_EMISSIVITY(Q,I)
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  DO UI = 1, N_USER_STREAMS
                    VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,UI) = &
                    VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,UI) - LOCAL_USER_EMISSIVITY(Q,UI)
                  ENDDO
                ENDIF
              ENDDO

!  Linearization w.r.t Kernel Factor

              W  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                W = W + 1
                DO Q = 1, NSTOKES
                  DO I = 1, NSTREAMS
                    VBRDF_LinSup_Out%BS_LS_EMISSIVITY(W,Q,I) = - LOCAL_EMISSIVITY(Q,I) / FF
                  ENDDO
                  IF ( DO_USER_STREAMS ) THEN
                    DO UI = 1, N_USER_STREAMS
                      VBRDF_LinSup_Out%BS_LS_USER_EMISSIVITY(W,Q,UI) = &
                      - LOCAL_USER_EMISSIVITY(Q,UI) / FF
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, BRDF_NPARS
                IF ( DERIVS(P) ) THEN
                  W = W + 1
                  DO Q = 1, NSTOKES
                    DO I = 1, NSTREAMS
                      VBRDF_LinSup_Out%BS_LS_EMISSIVITY(W,Q,I) = - D_LOCAL_EMISSIVITY(P,Q,I)
                    ENDDO
                    IF ( DO_USER_STREAMS ) THEN
                      DO UI = 1, N_USER_STREAMS
                        VBRDF_LinSup_Out%BS_LS_USER_EMISSIVITY(W,Q,UI) = &
                        - D_LOCAL_USER_EMISSIVITY(P,Q,UI)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VBRDF_LIN_MAINMASTER

!

      SUBROUTINE VBRDF_LIN_INPUTMASTER ( &
        FILNAM, VBRDF_Sup_In, VBRDF_LinSup_In, &
        VBRDF_Sup_InputStatus )

!  Input routine for BRDF program

      USE VLIDORT_PARS
      USE VLIDORT_FINDPAR_M

      USE vbrdf_sup_inputs_def
      USE vbrdf_sup_outputs_def

      USE vbrdf_linsup_inputs_def

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(VBRDF_Sup_Inputs)   , INTENT(OUT) :: VBRDF_Sup_In
      TYPE(VBRDF_LinSup_Inputs), INTENT(OUT) :: VBRDF_LinSup_In

      TYPE(VBRDF_Input_Exception_Handling), INTENT(OUT) :: &
        VBRDF_Sup_InputStatus

!  Local variables
!  ---------------

!  Stream angle flag

      LOGICAL ::          DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL ::          DO_BRDF_SURFACE

!  Surface emission

      LOGICAL ::          DO_SURFACE_EMISSION

!  Number of Stokes components

      INTEGER ::          NSTOKES

!  Number and index-list and names of bidirectional functions

      INTEGER ::            N_BRDF_KERNELS
      INTEGER ::            WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER ::          N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      DOUBLE PRECISION :: BRDF_PARAMETERS &
          ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL ::          LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      DOUBLE PRECISION :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER ::          NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL ::          DO_SHADOW_EFFECT

!  Exact only flag (no Fourier term calculations)

      LOGICAL ::          DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      INTEGER ::          MSRCORR_ORDER
      LOGICAL ::          DO_MSRCORR_EXACTONLY
      INTEGER ::          MSRCORR_NMUQUAD
      INTEGER ::          MSRCORR_NPHIQUAD

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL ::          DO_KERNEL_FACTOR_WFS  ( MAX_BRDF_KERNELS )
      LOGICAL ::          DO_KERNEL_PARAMS_WFS  ( MAX_BRDF_KERNELS, &
                                                  MAX_BRDF_PARAMETERS )

!  Derived quantity (tells you when to do BRDF derivatives)

      LOGICAL ::          DO_KPARAMS_DERIVS  ( MAX_BRDF_KERNELS )

!  Number of surface weighting functions

      INTEGER ::          N_SURFACE_WFS
      INTEGER ::          N_KERNEL_FACTOR_WFS
      INTEGER ::          N_KERNEL_PARAMS_WFS

!  Number of discrete ordinate streams

      INTEGER ::          NSTREAMS

!  Local angle control

      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Angles

      DOUBLE PRECISION :: BEAM_SZAS (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Exception handling. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=9), PARAMETER :: PREFIX = 'BRDFSUP -'

      INTEGER ::            DUM_INDEX, DUM_NPARS
      CHARACTER (LEN=10) :: DUM_NAME
      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, J, K, L, FILUNIT, NM

!  Check list of Kernel names

      CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )

      BRDF_CHECK_NAMES = (/ &
                           'Lambertian', &
                           'Ross-thin ', &
                           'Ross-thick', &
                           'Li-sparse ', &
                           'Li-dense  ', &
                           'Hapke     ', &
                           'Roujean   ', &
                           'Rahman    ', &
                           'Cox-Munk  ', &
                           'GissCoxMnk', &
                           'GCMcomplex', &
                           'BPDF2009  '/)

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES_INPUT(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

      NSTOKES = 0

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      DO_SHADOW_EFFECT    = .FALSE.
      DO_EXACTONLY        = .FALSE.
      DO_SURFACE_EMISSION = .FALSE.

      DO_MSRCORR           = .FALSE.
      MSRCORR_ORDER        = 0
      DO_MSRCORR_EXACTONLY = .FALSE.
      MSRCORR_NMUQUAD      = 0
      MSRCORR_NPHIQUAD     = 0

      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

      N_SURFACE_WFS  = 0
      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0
      DO K = 1, MAX_BRDF_KERNELS
        DO_KPARAMS_DERIVS(K) = .false.
        DO_KERNEL_FACTOR_WFS(K) = .FALSE.
        DO L = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(K,L) = .FALSE.
        ENDDO
      ENDDO

!  number of Stokes components
!  ===========================

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Read Angle stuff
!  ================

!  user-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                          ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Solar beams
!  ===========

!  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = &
        'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Azimuth angles
!  ==============

!  Number of angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
         'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

!  Angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  User defined stream angles (should be positive)
!  ==========================

      IF ( DO_USER_STREAMS ) THEN

!  Number of angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
          'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

      ENDIF

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Surface emission flag

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  BRDF inputs

      IF ( DO_BRDF_SURFACE ) THEN

!  number of kernels, check this value

        PAR_STR = 'Number of BRDF kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
            'Number of BRDF Kernels > maximum dimension (=3)'
          ACTIONS(NM)  = &
            'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of BRDF azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  = &
        'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Main kernel input

        PAR_STR = &
           'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
                BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
               N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
          ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag

        PAR_STR = 'Do Exact-only (no Fourier) kernels?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple reflectance correction (for Cox-Munk types); General flag

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do multiple reflectance for All glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Specific MSRCORR inputs

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Do multiple reflectance for Exact-only glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR_EXACTONLY
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR scattering order

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering order for glitter kernels'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_ORDER
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR quadrature orders

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Polar quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NMUQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Azimuth quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NPHIQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Check dimensions

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_NMUQUAD .gt. max_msrs_muquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR polar quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_muquad in vlidort_pars'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
           IF ( MSRCORR_NPHIQUAD .gt. max_msrs_phiquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR azimuth quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_phiquad in vlidort_pars'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check on MSRCORR order

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_ORDER .EQ.0 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order = 0'
              ACTIONS(NM)  = 'Turn off MSRCORR flags and proceed with warning'
              DO_MSRCORR = .false. ; DO_MSRCORR_EXACTONLY = .false.
              STATUS = VLIDORT_WARNING
              NMESSAGES = NM
           ENDIF
        ENDIF

!  TEMPORARY LIMITATION TO MSR ORDER = 1
!        IF ( DO_MSRCORR ) THEN
!           IF ( MSRCORR_ORDER .GT.1 ) then
!              NM = NM + 1
!              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order  > 1'
!              ACTIONS(NM)  = 'Temporary limitation; abort for now'
!              STATUS = VLIDORT_SERIOUS
!              NMESSAGES = NM
!              GO TO 764
!           ENDIF
!        ENDIF

!  Linearized input
!  ----------------

!  Not allowed linearized inputs with GCMCRI

        PAR_STR = &
        'Kernels, indices, # pars, Factor Jacobian flag, Par Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998) &
               DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I), &
               (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)

            IF ( DUM_NAME .EQ. 'GCMcomplex' ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'GCMcomplex BRDF Kernel cannot be linearized yet'
              ACTIONS(NM)  = 'Use GISSCoxMunk kernel instead'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF

            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Kernel name not same as earlier ' // &
              'list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel ' // &
              'name'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF

            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Index name not same as earlier ' // &
              'list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel ' // &
              'Index'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF

            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input Number of BRDF parameters not same ' // &
              'as earlier list'
              ACTIONS(NM)  = 'Check second occurence of ' // &
              'N_BRDF_PARAMETERS'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF

!  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
              N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
            ENDIF
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
                N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1
              ENDIF
            ENDDO
            DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)

          ENDDO
          N_SURFACE_WFS = N_KERNEL_FACTOR_WFS+N_KERNEL_PARAMS_WFS
        ENDIF

        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 )

!  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Surface WFs > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SURFACEWFS ' // &
          'dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Check Kernel indices are within bounds. Check BRDF name is on accepted list

        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
            NM = NM + 1
            MESSAGES(NM) = &
            'Bad input: BRDF Index not on list of indices'
            ACTIONS(NM)  = &
            'Re-set input value: Look in VLIDORT_PARS for correct index'
            STATUS = VLIDORT_SERIOUS
            NMESSAGES = NM
            GO TO 764
          ELSE
            IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
              NM = NM + 1
              MESSAGES(NM) = &
              'Bad input: BRDF kernel name not one of accepted list'
              ACTIONS(NM)  = &
              'Re-set input value: Look in VLIDORT_PARS for correct name'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
          ENDIF
        ENDDO

!  End BRDF clause

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy Control inputs

      VBRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      VBRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

!  Copy Geometry results

      VBRDF_Sup_In%BS_NSTOKES           = NSTOKES
      VBRDF_Sup_In%BS_NSTREAMS          = NSTREAMS
      VBRDF_Sup_In%BS_NBEAMS            = NBEAMS
      VBRDF_Sup_In%BS_BEAM_SZAS         = BEAM_SZAS
      VBRDF_Sup_In%BS_N_USER_RELAZMS    = N_USER_RELAZMS
      VBRDF_Sup_In%BS_USER_RELAZMS      = USER_RELAZMS
      VBRDF_Sup_In%BS_N_USER_STREAMS    = N_USER_STREAMS
      VBRDF_Sup_In%BS_USER_ANGLES_INPUT = USER_ANGLES_INPUT

!  Copy BRDF inputs

      VBRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS
      VBRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      VBRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      VBRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      VBRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS
      VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG
      VBRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      VBRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF

      VBRDF_Sup_In%BS_DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT
      VBRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY

      VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR           = DO_MSRCORR
      VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_EXACTONLY = DO_MSRCORR_EXACTONLY
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER        = MSRCORR_ORDER
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD      = MSRCORR_NMUQUAD
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD     = MSRCORR_NPHIQUAD

!  Copy linearized BRDF inputs

      VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS   = DO_KERNEL_FACTOR_WFS
      VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS   = DO_KERNEL_PARAMS_WFS
      VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS      = DO_KPARAMS_DERIVS
      VBRDF_LinSup_In%BS_N_SURFACE_WFS          = N_SURFACE_WFS
      VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = N_KERNEL_FACTOR_WFS
      VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = N_KERNEL_PARAMS_WFS

!  Exception handling

      VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VBRDF_LIN_INPUTMASTER

      END MODULE vbrdf_LinSup_masters_m

