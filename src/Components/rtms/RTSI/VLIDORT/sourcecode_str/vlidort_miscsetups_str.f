C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            VLIDORT_MISCSETUPS (master, calling:)            #
C #              VLIDORT_DELTAMSCALE                            #
C #              VLIDORT_SSALBINIT                              #
C #              VLIDORT_QSPREP                                 #
C #              VLIDORT_PREPTRANS                              #
C #                                                             #
C #              LAMBERTIAN_SURFACE_DIRECTBEAM                  #
C #              VLIDORT_LAMBERTIAN_DBCORRECTION                #
C #                                                             #
C #              VLIDORT_CHAPMAN                                #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_MISCSETUPS

C  miscellaneous setup operations, Master routine

c  Performance set-up is in VLIDORT_DERIVE_INPUTS (vlidort_inputs.f)
c      CALL VLIDORT_PERFORMANCE_SETUP

C  Delta-m scaling of input quantities

      CALL VLIDORT_DELTAMSCALE

C  initialise single scatter albedo terms
C    GREEKMAT x OMEGA

      CALL VLIDORT_SSALBINIT

C  Prepare quasi-spherical attenuation

      CALL VLIDORT_QSPREP

C  Transmittances and Transmittance factors

      CALL VLIDORT_PREPTRANS

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_DELTAMSCALE

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of set up variables (output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  local variables

      DOUBLE PRECISION FDEL, FAC2, DNL1, FDNL1
      DOUBLE PRECISION DNM1, DELS, DT, XTD
      INTEGER          K, K1, K2, N, N1, L, UT, UTA, NM1, IB

C  Indexing

      INTEGER          KTYPE1(4), KTYPE2(4)
      DATA             KTYPE1 / 1, 6, 11, 16 /
      DATA             KTYPE2 / 2, 5, 12, 15 /

C  First do the slant input
C  ------------------------

C  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

C  debug

c      do n = 1, nlayers
c       dels = height_grid(n-1)-height_grid(n)
c       write(*,*)n,dels,height_grid(n)
c       write(44,'(i4,101f10.5)')n,
c     &  (CHAPMAN_FACTORS(N,K,1)*(height_grid(k-1)-height_grid(k))
c     &           ,k=1,nlayers)
c      enddo
c      pause

C  DELTAM SCALING
C  ==============

      IF ( DO_DELTAM_SCALING ) THEN

C  New section by R. Spurr, RT Solutions Inc.
C   Based in part on the code in VDISORT.

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DFLOAT(2*NM1+1)

C  Scaling for layer input
C  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

C  overall truncation factor

          FDEL = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

C  Scale Greek Matrix entries

          DO L = 0, NMOMENTS
            DNL1  = DFLOAT(2*L + 1 )
            FDNL1 = FDEL * DNL1
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) =
     &            ( GREEKMAT_TOTAL_INPUT(L,N,K1) - FDNL1 ) / FAC2
                GREEKMAT_TOTAL(L,N,K2) =
     &              GREEKMAT_TOTAL_INPUT(L,N,K2) / FAC2
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) =
     &            ( GREEKMAT_TOTAL_INPUT(L,N,1) - FDNL1 ) / FAC2
            ENDIF
          ENDDO

C  Maintain phase function normalization

          GREEKMAT_TOTAL(0,N,1) = ONE

C  scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

C  end layer loop

        ENDDO

C  Scaling for user-defined off-grid optical depths
C  ------------------------------------------------

C  Scaling for user-defined off-grid optical depths
C     (on-grid values have already been scaled)

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

C  Scale layer path thickness values

        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

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
              GREEKMAT_TOTAL(L,N,K1) =
     &              GREEKMAT_TOTAL_INPUT(L,N,K1)
              GREEKMAT_TOTAL(L,N,K2) =
     &              GREEKMAT_TOTAL_INPUT(L,N,K2)
            ENDDO
          ENDDO
        ENDDO

C  Scaling for user-defined off-grid optical depths
C     (on-grid values have already been scaled)

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

C  Finish module

      RETURN
      END

C
    
      SUBROUTINE VLIDORT_SSALBINIT

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of set up variables (output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  local variables

      INTEGER          N, L, O1, O2, OM
      DOUBLE PRECISION SUM

C  phase matrix-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          DO O1 = 1, NSTOKES
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              SUM = OMEGA_TOTAL(N)*GREEKMAT_TOTAL(L,N,OM)
              !WRITE(92, *) N,L,OM,OMEGA_TOTAL(N),GREEKMAT_TOTAL(L,N,OM)
              OMEGA_GREEK(L,N,O1,O2)   = SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_QSPREP

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of set up variables (local output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Local variables
C  ---------------

      INTEGER          N, K, IB
      DOUBLE PRECISION S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)

C  Specialist code
C  ---------------

C  Cutting off the solar source arbitrarily

      IF ( DO_SPECIALIST_OPTION_3 ) THEN
        DO IB = 1, NBEAMS
          LAYER_PIS_CUTOFF(IB) = NLAYERS_CUTOFF
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          LAYER_PIS_CUTOFF(IB) = NLAYERS
        ENDDO
      ENDIF

C  TOA  cosines

      DO IB = 1, NBEAMS
        LOCAL_CSZA(0,IB) = COS_SZANGLES(IB)
      ENDDO
    
C  plane-parallel case
C  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       DO IB = 1, NBEAMS
        SEC0 = ONE / COS_SZANGLES(IB)
c        LAYER_PIS_CUTOFF(IB) = NLAYERS
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

        S_T_1 = ZERO
        S_T_0 = ONE
c        LAYER_PIS_CUTOFF(IB) = NLAYERS
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
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

C  Set the Local solar zenith cosines
C  Distinguish between the refractive and non-refractive cases.

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

C  Set Direct Beam Flag and solar beam total attenuation to surface

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

c debug
c      do n = 1, nlayers
c        write(*,*)AVERAGE_SECANT(N,1),  INITIAL_TRANS(N,1)
c      enddo
c      pause

C  finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_PREPTRANS

C  Prepare transmittances and transmittance factors

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables (output stored here)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  local variables
C  ---------------

      INTEGER          N, UT, UM, IB, I
      DOUBLE PRECISION XT, SPHER, HELP

C  Transmittance factors for discrete ordinate streams
C  ===================================================

C  New code by R. Spurr, RT Solutions, 12 April 2005
C    Required for the solution saving option
C  Partial layer code added 30 December 2005 by R. Spurr

      IF ( DO_SOLUTION_SAVING ) THEN

C  whole layers

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

C  Atmosphere Partial layers

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

C  Transmittance factors for average secant stream
C  ===============================================

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

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         IF ( N. GT. LAYER_PIS_CUTOFF(IB) ) THEN
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

C  end solar beam loop

      ENDDO

C  Transmittances for User Streams
C  ===============================

C  return if not flagged

      IF ( .NOT. DO_USER_STREAMS ) RETURN

C  Initial transmittances divided by user streams
C  ----------------------------------------------

      DO IB = 1, NBEAMS
       DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
         ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
        ENDDO
       ENDDO
      ENDDO

C  Whole Layer transmittances
C  --------------------------

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

C  Finish

      RETURN
      END

C

      SUBROUTINE LAMBERTIAN_SURFACE_DIRECTBEAM
     I     ( FOURIER_COMPONENT,
     I       DELTA_FACTOR )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include file of Direct beam reflectance variables (output goes here)

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  module arguments
C  ----------------

      DOUBLE PRECISION DELTA_FACTOR
      INTEGER          FOURIER_COMPONENT

C  Local variables
C  ---------------

      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN, REFLEC
      INTEGER          I, UI, O1, IB

C  Initialize
C  ----------

C   Safety first!  Return if there is no reflection.

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

C  Attenuation of solar beam
C  -------------------------

C  New code to deal with refractive geometry case
C   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
        ELSE
         X0_BOA = COS_SZANGLES(IB)
        ENDIF

C  There should be no flux factor here.
C    Bug fixed 18 November 2005. Earlier Italian Job!!
C       X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR / DELTA_FACTOR

        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

C  Set value for Fourier = 0
C  Natural light only, so only the first Stokes component is nonzero

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

C  end direct beam calculation

       ENDIF
      ENDDO

C  finish

      RETURN
      END

c

      SUBROUTINE VLIDORT_CHAPMAN
     O   ( FAIL, MESSAGE, TRACE )

C  This is the internal Chapman function calculation of the slant path
C  Chapman Factors required for slant-path optical thickness values. 

C  This module calculates CHAPMAN_FACTORS internally inside VLIDORT,
C  saving you the job of doing it yourself, though you can input these
C  quantities, if the DO_CHAPMAN_FACTORS flag is not set.

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

      INCLUDE '../includes/VLIDORT.PARS'

C  Input file

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Derived Geophysical variables (output stored here)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

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

      DO IB = 1, N_SZANGLES

        SUN0 = SZANGLES(IB)

        CALL BEAM_GEOMETRY_PREPARE
     I     ( MAX_SZANGLES, MAXLAYERS, IB, NLAYERS, FINEGRID,
     I       SUN0, EARTH_RADIUS, RFINDEX_PARAMETER,
     I       DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,
     I       HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,
     O       CHAPMAN_FACTORS, SZA_LOCAL_INPUT,
     O       ITERSAVE, FAIL, MESSAGE )

C  return if failed

        IF ( FAIL ) THEN
          TRACE = 'Geometry failure in VLIDORT_CHAPMAN'
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
      PARAMETER  ( LOCAL_MAXLAYERS     = 101 )
      PARAMETER  ( LOCAL_MAXFINELAYERS = 20 )
      
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

