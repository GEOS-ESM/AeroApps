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
C #            LIDORT_MISCSETUPS (master), calling--            #
C #              LIDORT_PERFORMANCE_SETUP                       #
C #              LIDORT_DELTAMSCALE                             #
C #              LIDORT_SSALBINIT                               #
C #              LIDORT_QSPREP                                  #
C #              LIDORT_PREPTRANS                               #
C #                                                             #
C #            LIDORT_DIRECTBEAM   (called in LIDORT_FOURIER)   #
C #                                                             #
C #            LIDORT_CHAPMAN   (called by LIDORT_MASTER)       #
C #              BEAM_GEOMETRY_PREPARE                          #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_MISCSETUPS

C  miscellaneous setup operations, Master routine

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Code
C  ====

c  Performance set-up

      CALL LIDORT_PERFORMANCE_SETUP

C  Delta-m scaling of input quantities
C    ---- Trivial if no scattering

      CALL LIDORT_DELTAMSCALE

C  initialise single scatter albedo terms

      CALL LIDORT_SSALBINIT

C  Prepare quasi-spherical attenuation
C    -----Only if the solar source terms are present
C    -----routine will exit otherwise

      CALL LIDORT_QSPREP

C  Transmittances and Transmittance factors

      CALL LIDORT_PREPTRANS

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_DELTAMSCALE

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of set up variables (output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  local variables

      DOUBLE PRECISION FDEL, FAC2, DNL1, FDNL1
      DOUBLE PRECISION TAU, TAUFAC, DNM1, DELS
      INTEGER          N, N1, L, UT, UTA, NM1, K, IB

C  DELTAM SCALING
C  ==============

      IF ( DO_DELTAM_SCALING ) THEN

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DBLE(2*NM1+1)

C  Scaling for layer input
C  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

C  overall truncation factor

          FDEL = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

C  Scale Greek Matrix entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            PHASMOMS_TOTAL(L,N) =
     &            ( PHASMOMS_TOTAL_INPUT(L,N) - FDNL1 ) / FAC2
          ENDDO

C  Maintain phase function normalization

          PHASMOMS_TOTAL(0,N) = ONE

C  scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

C  end layer loop

        ENDDO

C  Scaling for user-defined off-grid optical depths
C     (on-grid values have already been scaled)

        IF ( DO_USER_TAUS ) THEN
          UT = 0
          DO UTA = 1, N_OUT_USERTAUS
            IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
              UT = UT + 1
              N = OFFGRID_UTAU_LAYERIDX(UT)
              TAU = USER_TAUS_INPUT(UTA) - TAUGRID_INPUT(N-1)
              TAUFAC = TAU * FAC1(N)
              OFFGRID_UTAU_VALUES(UT) = TAUFAC
            ENDIF
          ENDDO
        ENDIF

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          TAUGRID(N)     = TAUGRID_INPUT(N)
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N)
        ENDDO

        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
          DO N = 1, NLAYERS
            OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
            DO L = 0, NMOMENTS
              PHASMOMS_TOTAL(L,N) = PHASMOMS_TOTAL_INPUT(L,N)
            ENDDO
          ENDDO
        ENDIF

        IF ( DO_USER_TAUS ) THEN
          UT = 0
          DO UTA = 1, N_OUT_USERTAUS
            IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
              UT = UT + 1
              N = OFFGRID_UTAU_LAYERIDX(UT)
              TAU = USER_TAUS_INPUT(UTA) - TAUGRID(N-1)
              OFFGRID_UTAU_VALUES(UT) = TAU
            ENDIF
          ENDDO
        ENDIF

      ENDIF

C  If no solar terms, finish

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Slant optical thickness values + scaling
C  ----------------------------------------

C  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT_INPUT(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

C  Scale layer path thickness values, or not (just copy input

      IF ( DO_DELTAM_SCALING ) THEN
        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) =
     *              DELTAU_SLANT_INPUT(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT_INPUT(N,K,IB)
            ENDDO
          ENDDO
        ENDDO
       ENDIF

C  Finish module

      RETURN
      END

C
    
      SUBROUTINE LIDORT_SSALBINIT

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of set up variables (output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  local variables

      INTEGER             N, L

C  phase moment-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          OMEGA_MOMS(N,L) = OMEGA_TOTAL(N)*PHASMOMS_TOTAL(L,N)
        ENDDO
      ENDDO

C  Finish

      END

C

      SUBROUTINE LIDORT_QSPREP

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of set up variables (local output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Local variables
C  ---------------

      INTEGER          N, K, IB
      DOUBLE PRECISION S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)

C  Nothing to do if no solar sources
C  ---------------------------------

      IF ( .NOT.DO_SOLAR_SOURCES ) RETURN

C  plane-parallel case
C  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       S_T_0 = ONE
       DO IB = 1, NBEAMS
        SEC0 = ONE / X0(IB)
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_SZA(N,IB)      = X0(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
            LOCAL_SZA(N,IB)      = ZERO
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

C  pseudo-spherical case
C  ---------------------

       DO IB = 1, NBEAMS

C  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = ZERO
        DO N = 1, NLAYERS
          TAU = ZERO
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

C  set up the average secant formulation

        S_T_0 = ONE
        S_T_1 = ZERO
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) 
     &                               / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            LOCAL_SZA(N,IB)      = X0(IB)
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

C  Set the Local solar zenith angles
C  Distinguish between the refractive and non-refractive cases.

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              LOCAL_SZA(N,IB) = SUN_SZA_COSINES(N,IB)
            ELSE
              LOCAL_SZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              LOCAL_SZA(N,IB) = X0(IB)
            ELSE
              LOCAL_SZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ENDIF

       ENDDO
      ENDIF

C  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        IF ( TAU_SOLAR(IB) .GT. MAX_TAU_SPATH ) THEN
          SOLAR_BEAM_OPDEP(IB) = ZERO
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          SOLAR_BEAM_OPDEP(IB) = DEXP( - TAU_SOLAR(IB) )
          DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

C  debug

c      k = 97
c      if ( do_fdtest ) k = 98
c      do n = 1, nlayers
c       write(k,'(1p9e15.7)')(average_secant(n,ib),ib=1,nbeams)
c       write(k,'(1p9e15.7)')(initial_trans(n,ib),ib=1,nbeams)
c      enddo
c      if ( do_fdtest ) pause
 
C  finish

      RETURN
      END

C

      SUBROUTINE LIDORT_PREPTRANS

C  Prepare transmittances and transmittance factors

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup variables (output stored here)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of Multiplier variables (output stored here)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  local variables
C  ---------------

      INTEGER          N, UT, UM, IB, I
      DOUBLE PRECISION XT, SPHER, HELP

C  Transmittance factors for discrete ordinate streams
C  ===================================================

C  New code by R. Spurr, RT Solutions, 12 April 2005
C  Off-grid optical depths, RT Solutions, 30 August 2005

C  Only required for the solution saving option
C    (automatic for BVP telescoping)

      IF ( DO_SOLUTION_SAVING ) THEN

C  whole layers

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS 
            SPHER = DELTAU_VERT(N) / X(I)
            IF ( SPHER .GT. MAX_TAU_QPATH ) THEN
              T_DELT_DISORDS(I,N) = ZERO
            ELSE
              T_DELT_DISORDS(I,N) = DEXP ( - SPHER )
            ENDIF
          ENDDO
        ENDDO

C  Atmosphere Partial layers

        DO UT = 1, N_OFFGRID_USERTAUS
          N = OFFGRID_UTAU_LAYERIDX(UT)
          XT = OFFGRID_UTAU_VALUES(UT)
          DO I = 1, NSTREAMS
            HELP =  XT / X(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTDN(I,UT) = ZERO
            ELSE
              T_DISORDS_UTDN(I,UT) = DEXP(-HELP)
            ENDIF
            HELP = ( DELTAU_VERT(N) - XT ) / X(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTUP(I,UT) = ZERO
            ELSE
              T_DISORDS_UTUP(I,UT) = DEXP(-HELP)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

C  Transmittance factors for average secant stream
C  ===============================================

C   Only if solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  start solar loop

       DO IB = 1, NBEAMS

C  Whole layer Transmittance factors
C  ---------------------------------

C  layer transmittance 

        DO N = 1, NLAYERS
         IF ( N. GT. LAYER_PIS_CUTOFF(IB) ) THEN
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

C  Partial layer transmittance factors (for off-grid optical depths)
C  -----------------------------------------------------------------

        DO UT = 1, N_OFFGRID_USERTAUS
         N = OFFGRID_UTAU_LAYERIDX(UT)
         IF ( N. GT. LAYER_PIS_CUTOFF(IB) ) THEN
           T_UTDN_MUBAR(UT,IB) = ZERO
           T_UTUP_MUBAR(UT,IB) = ZERO
         ELSE
           XT = OFFGRID_UTAU_VALUES(UT)
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

C  end solar beam loop and solar sources clause

       ENDDO
      ENDIF

C  Transmittances for User Streams
C  ===============================

C  return if not flagged

      IF ( .NOT. DO_USER_STREAMS ) RETURN

C  Initial transmittances divided by user streams
C    ---- Only for solar sources

      IF ( DO_SOLAR_SOURCES ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         DO UM = 1, N_USER_STREAMS
          ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
         ENDDO
        ENDDO
       ENDDO
      ENDIF

C  Whole Layer transmittances

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

C  Partial Layer transmittances for off-grid optical depths

      DO UT = 1, N_OFFGRID_USERTAUS
        N  = OFFGRID_UTAU_LAYERIDX(UT)
        XT = OFFGRID_UTAU_VALUES(UT)
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

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_DIRECTBEAM
     I     ( FOURIER_COMPONENT,
     I       DELTA_FACTOR )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include file of reflection variables (output goes in here)

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  module arguments
C  ----------------

      DOUBLE PRECISION DELTA_FACTOR
      INTEGER          FOURIER_COMPONENT

C  Local variables
C  ---------------

      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN, REFL_ATTN
      INTEGER          I, UI, KL, IB

C  Initialize
C  ----------

C   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DIRECT_BEAM(I,IB) = ZERO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            USER_DIRECT_BEAM(UI,IB) = ZERO
          ENDDO
        ENDIF
      ENDDO

C  Attenuation of solar beam
C  -------------------------

C  New code to deal with refractive geometry case
C   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
        ELSE
         X0_BOA = X0(IB)
        ENDIF

C  There should be no flux factor here.
C    Bug fixed 18 November 2005. Earlier Italian Job!!
C    Flux Factor put back, 1 March 2007. Using 1 / pi4

        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4

        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

C  Start loop over albedo kernels
C    - reflected along the quadrature directions
C    - reflected along the user-defined directions

        DO KL = 1, N_BRDF_KERNELS

          REFL_ATTN = BRDF_FACTORS(KL) * ATTN

C  For the Lambertian case
C  -----------------------

C  All terms are the same; zero for Fourier not zero

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

C  zero it!

            DO I = 1, NSTREAMS
              A_DIRECT_BEAM(I,IB,KL) = ZERO
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                A_USER_DIRECT_BEAM(UI,IB,KL) = ZERO
              ENDDO
            ENDIF

C  Set value for Fourier = 0

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              DO I = 1, NSTREAMS
                A_DIRECT_BEAM(I,IB,KL) = REFL_ATTN
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  A_USER_DIRECT_BEAM(UI,IB,KL) = REFL_ATTN
                ENDDO
              ENDIF
            ENDIF

C  For the Bidirectional case
C  --------------------------

          ELSE

            DO I = 1, NSTREAMS
              A_DIRECT_BEAM(I,IB,KL) = 
     &           REFL_ATTN * BIREFLEC_0(KL,I,IB)
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                A_USER_DIRECT_BEAM(UI,IB,KL) = 
     *           REFL_ATTN * USER_BIREFLEC_0(KL,UI,IB)
              ENDDO
            ENDIF

          ENDIF

C  Finish Kernel loop

        ENDDO

C  Total contributions
C  -------------------

        DO KL = 1, N_BRDF_KERNELS
          DO I = 1, NSTREAMS
            DIRECT_BEAM(I,IB) = DIRECT_BEAM(I,IB) +
     &         A_DIRECT_BEAM(I,IB,KL)
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_DIRECT_BEAM(UI,IB) = 
     &                USER_DIRECT_BEAM(UI,IB) + 
     &               A_USER_DIRECT_BEAM(UI,IB,KL)
              ENDDO
          ENDIF
        ENDDO

C  end direct beam calculation

       ENDIF
      ENDDO

C  finish

      RETURN
      END

C

      SUBROUTINE LIDORT_PERFORMANCE_SETUP

C  Performance setup of flags for avoiding solutions
C  -------------------------------------------------

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of set up variables (output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Local variables
C  ---------------

      INTEGER             M, N, Q, QC, L
      LOGICAL             LOOP

C  Set performance flags
C  ---------------------

C  New section for Version 3.0.

C  Set the layer scattering flags
C  initialise (M = Fourier index) to normal  mode.

      DO M = 0, NMOMENTS
        BVP_REGULAR_FLAG(M) = .TRUE.
        DO N = 1, NLAYERS
          DO_LAYER_SCATTERING(M,N) = .TRUE.
        ENDDO
      ENDDO

C  Special clause for transmittance-only calculation
C   New flag, 31 July 2007. Return with zero SSA

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO_SOLUTION_SAVING = .TRUE.
        DO N = 1, NLAYERS
          OMEGA_TOTAL_INPUT(N) = ZERO
          DO_LAYER_SCATTERING(0,N) = .FALSE.
        ENDDO
        RETURN
      ENDIF
      
C  for M > 2 terms, examine Phase function moments -
C  no scattering if they are all zero for L > M - 1.

      IF ( DO_SOLAR_SOURCES ) THEN
       DO M = 3, NMOMENTS
        QC = NMOMENTS - M + 1
        DO N = 1, NLAYERS
          Q = 0
          DO L = M, NMOMENTS
            IF(PHASMOMS_TOTAL_INPUT(L,N).EQ.ZERO)Q=Q+1
          ENDDO
          DO_LAYER_SCATTERING(M,N) = (Q.LT.QC)
        ENDDO
       ENDDO
      ENDIF

C  Re-set solution saving flag internally
C   Do not do this.....

c      IF ( DO_SOLAR_SOURCES ) THEN
c       IF ( .NOT.DO_SOLUTION_SAVING ) THEN
c        DO M = 3, NMOMENTS
c          Q = 0
c          DO N = 1, NLAYERS
c            IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
c          ENDDO
c          IF ( Q.GT.1) DO_SOLUTION_SAVING = .TRUE.
c        ENDDO
c       ENDIF
c      ENDIF
      
C  BVP telescoping (only if do_solution_saving is set)

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        IF ( DO_BVP_TELESCOPING ) THEN
          DO M = 3, NMOMENTS
            Q = 0
            DO N = 1, NLAYERS
              IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
            ENDDO
            IF ( Q.GT.1) BVP_REGULAR_FLAG(M) = .FALSE.
          ENDDO
        ENDIF
       ENDIF
      ENDIF
      
C  Single scattering (atmosphere) usable moments in each layer
C   Bug. 2 March 2007. L must be < nmoments_input in DO WHILE
C      ( Same bug in VLIDORT, corrected 22 November 2006 )

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        IF ( DO_RAYLEIGH_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 2
          ENDDO
        ELSE IF ( DO_ISOTROPIC_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 0
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            L = 2
            LOOP = .TRUE.
            DO WHILE (LOOP.AND.L.LT.NMOMENTS_INPUT)
              L = L + 1
              LOOP= (PHASMOMS_TOTAL_INPUT(L,N).NE.ZERO)
            ENDDO
            LAYER_MAXMOMENTS(N) = L - 1
          ENDDO
        ENDIF
       ENDIF
      ENDIF
      
C  Finish

      RETURN
      END

c

      SUBROUTINE LIDORT_CHAPMAN
     O   ( FAIL, MESSAGE, TRACE )

C  This is the internal Chapman function calculation of the slant path
C  optical thickness array DELTAU_SLANT_INPUT. 

C  This module calculates DELTAU_SLANT_INPUT internally inside LIDORT
C  saving you the job of doing it yourself.

C  The following options apply:

C   1. If the plane-parallel flag is on, no further inputs are required

C   2. If the plane-parallel flag is off, then Pseudo-spherical:
C       (a) Straight line geometry, must specify
C               Earth_radius, height grid
C       (b) Refractive geometry, must specify
C               Earth_radius, height grid
C               pressure grid, temperature grid

C  The logic will be checked before the module is called.

C  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

C    Based round a call to a pure geometry module which returns slant
C    path distances which was adapted for use in the Radiant model by
C    R. Spurr during an OCO L2 intensive April 24-29, 2005.

C  Include files
C  -------------

C  parameter file for dimensioning

      INCLUDE '../includes/LIDORT.PARS'

C  Input file

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Derived Geophysical variables (output stored here)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  arguments
C  ---------

C  output status

      LOGICAL          FAIL
      CHARACTER*(*)    MESSAGE, TRACE

C  Local variables
C  ---------------

C  number of iterations (refractive case only)
C      This is debug output

      INTEGER          ITERSAVE(MAXLAYERS)

C  other local variables

      INTEGER          IB
      DOUBLE PRECISION SUN0

C  get spherical optical depths
C  ----------------------------

C  start beam loop

      DO IB = 1, NBEAMS

        SUN0 = BEAM_SZAS(IB)

        CALL BEAM_GEOMETRY_PREPARE
     I     ( MAXBEAMS, MAXLAYERS, IB, NLAYERS, FINEGRID,
     I       SUN0, EARTH_RADIUS, RFINDEX_PARAMETER,
     I       DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,
     I       HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,
     O       CHAPMAN_FACTORS, SZA_LOCAL_INPUT,
     O       ITERSAVE, FAIL, MESSAGE )

C  return if failed

        IF ( FAIL ) THEN
          TRACE = 'Geometry failure in LIDORT_CHAPMAN'
          RETURN
        ENDIF

C  end beam loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE BEAM_GEOMETRY_PREPARE
     I     ( MAXBEAMS, MAXLAYERS, IBEAM, NLAYERS, FINEGRID,
     I       SZA_GEOM_TRUE, REARTH, RFINDEX_PARAMETER,
     I       DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,
     I       HEIGHTS, PRESSURES, TEMPERATURES,
     O       CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT,
     O       ITERSAVE, FAIL, MESSAGE )

      IMPLICIT NONE
        
C  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
C  for a curved ray-traced beam through a multilayer atmosphere.

C  Coarse layering is input to the module. Values of Z, P, T  are
C  given at the layer boundaries, with the first value (index 0) at TOA.

C  The refractive geometry is assumed to start at the TOA level.

C  We also require the earth radius and the refractive index parameter
C   (For the Born-Wolf approximation)

C  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
C  In this case we do not require pressure and temperature information.
C  The calculation will then be for geometric rays.

C  The plane parallel Flag and the refractive geometry flag should not
C  both be true - this should be checked outside.

C  In the refracting case, fine-gridding of pressure and temperature is
C  done internally, temperature is interpolated linearly with height,
C  and pressure log-linearly. The refraction uses Snell's law rule.
C  Finelayer gridding assumes equidistant heights within coarse layers
C  but the number of fine layers can be varied

C  Output is specified at coarse layer boundaries

C  Module is stand-alone.

C  Reprogrammed for the OCO L2 algorithm
C   R. Spurr, RT Solutions, Inc.   April 27, 2005

C  Intended use in LIDORT and Radiant RT models.

C  Input arguments
C  ===============

C  input dimensioning

      INTEGER          MAXLAYERS, MAXBEAMS

C  Beam index

      INTEGER          IBEAM

C  number of coarse layers

      INTEGER          NLAYERS

C  number of fine layers within coarse layers

      INTEGER          FINEGRID(MAXLAYERS)

C  True solar zenith angle (degrees)

      DOUBLE PRECISION SZA_GEOM_TRUE

C  Earth radius (km)

      DOUBLE PRECISION REARTH
        
C  Refractive index parametaer (Born-Wolf approximation)

      DOUBLE PRECISION RFINDEX_PARAMETER

C  flag for plane parallel case

      LOGICAL          DO_PLANE_PARALLEL
        
C  flag for refractive geometry

      LOGICAL          DO_REFRACTIVE_GEOMETRY

C  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION HEIGHTS     (0:MAXLAYERS)
      DOUBLE PRECISION PRESSURES   (0:MAXLAYERS)
      DOUBLE PRECISION TEMPERATURES(0:MAXLAYERS)

C  Output arguments
C  ================

C  Path segments distances (km)

      DOUBLE PRECISION CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS,MAXBEAMS)

C  solar zenith angles at nadir

      DOUBLE PRECISION SZA_LEVEL_OUTPUT(0:MAXLAYERS,MAXBEAMS)

C  number of iterations (refractive case only)
C   This is debug output

      INTEGER          ITERSAVE(MAXLAYERS)

C  output status

      LOGICAL          FAIL
      CHARACTER*(*)    MESSAGE

C  Local variables
C  ===============

C  local dimensioning

      INTEGER    LOCAL_MAXLAYERS
      INTEGER    LOCAL_MAXFINELAYERS
      PARAMETER  ( LOCAL_MAXLAYERS     = 200 )
      PARAMETER  ( LOCAL_MAXFINELAYERS = 10 )
      
C  fine layer gridding for refraction

      DOUBLE PRECISION ZRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION PRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION TRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)

C  local height arrays

      DOUBLE PRECISION H(0:LOCAL_MAXLAYERS)
      DOUBLE PRECISION DELZ(LOCAL_MAXLAYERS)

C  help variables

      INTEGER          N, J, NRFINE, K, ITER, MAXF,IB
      LOGICAL          LOOP
      DOUBLE PRECISION GM_TOA, TH_TOA, MU_TOA, MU_NEXT
      DOUBLE PRECISION Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q
      DOUBLE PRECISION FU, FL, DEG_TO_RAD

        DOUBLE PRECISION LAYER_DIST,
     &      MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK,
     &      PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO,
     &      RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL

C  Standard temperature (K) and pressure (mbar).

      DOUBLE PRECISION     T_STANDARD
      PARAMETER            ( T_STANDARD = 273.16D0 )

      DOUBLE PRECISION     P_STANDARD
      PARAMETER            ( P_STANDARD = 1013.25D0 )

      DOUBLE PRECISION     STP_RATIO
      PARAMETER            ( STP_RATIO = T_STANDARD / P_STANDARD )

C  Loschmidt's number (particles/cm2/km).

      DOUBLE PRECISION  RHO_STANDARD
      PARAMETER       ( RHO_STANDARD = 2.68675D+24 )

C  Some setup operations
C  =====================

C  initialise output

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

C  check local dimensioning

      IF ( LOCAL_MAXLAYERS .LT. NLAYERS ) THEN
        MESSAGE = 'local coarse layer dimensioning insufficient'
        FAIL = .TRUE.
        RETURN
      ENDIF

C  Check fine layers do not exceed local dimensions assigned

      MAXF = 0
      DO N = 1, NLAYERS
        MAXF = MAX(MAXF,FINEGRID(N))
      ENDDO
      IF ( LOCAL_MAXFINELAYERS .LT. MAXF ) THEN
        MESSAGE = 'local fine layer dimensioning insufficient'
        FAIL = .TRUE.
        RETURN
      ENDIF
        
C  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

C  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )

C  initialize

      STH2D  = 0.0D0
        
C  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

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
c              write(*,'(i3,11F10.4)')N,(PRFINE(N,J),J=0,NRFINE)
          Z0 = Z1
          P0 = P1
          T0 = T1
          Q0 = Q1
        ENDDO
      ENDIF

C  plane-parallel case
C  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF
      
C  Refractive Geometry case
C  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        
C  Start value of SZA cosine

        MU_PREV = MU_TOA

C  start layer loop

        DO N = 1, NLAYERS

C  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1 = DASIN(SINTH1)
          PHI_0 = TH_TOA - STH1
          NRFINE = FINEGRID(N)
            
C  iteration loop

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

C  examine convergence

            DELPHI = PHI_0 - PHI_CUM
            LOOP = (DABS(DELPHI/PHI_CUM).GT.1.0D-4)

C  Fudge factors to speed up the iteration

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

C  failure

          IF ( LOOP ) THEN
            MESSAGE = 'refractive iteration not converged'
            FAIL = .TRUE.
            RETURN
          ENDIF
            
C  Update and save angle output

          MU_NEXT = DCOS(STH2D)
          MU_PREV = MU_NEXT
          SZA_LEVEL_OUTPUT(N,IB) = DACOS(MU_NEXT) / DEG_TO_RAD
          ITERSAVE(N) = ITER

        ENDDO

C  Straight line geometry
C  ======================

      ELSE
 
        DO N = 1, NLAYERS

C  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1   = DASIN(SINTH1)
          RE_UPPER = H(0)

C  solar zenith angles are all the same = input value

          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

C loop over layers K from 1 to layer N

          DO K = 1, N

C  sine-rule; PHI = earth-centered angle

            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = DASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = DSIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

C  re-set

            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2

          ENDDO

C  finish main layer loop

        ENDDO

C  Finish

      ENDIF

C  end of routine

      RETURN
      END

