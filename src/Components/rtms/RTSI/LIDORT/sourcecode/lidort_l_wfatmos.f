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

C ##########################################################
C #                                                        #
C # Subroutines in this Module                             #
C #                                                        #
C #     Top level routines--------------                   #
C #            UPUSER_ATMOSWF                              #
C #            DNUSER_ATMOSWF                              #
C #            MIFLUX_ATMOSWF                              #
C #                                                        #
C #     Output at quadrature angles ---------              #
C #            QUADATMOSWF_LEVEL_UP                        #
C #            QUADATMOSWF_LEVEL_DN                        #
C #            QUADATMOSWF_OFFGRID_UP                      #
C #            QUADATMOSWF_OFFGRID_DN                      #
C #                                                        #
C #     Post-processing at user angles --------            #
C #            GET_L_TOASOURCE                             #
C #            GET_L_BOASOURCE                             #
C #            L_WHOLELAYER_STERM_UP                       #
C #            L_WHOLELAYER_STERM_DN                       #
C #            QUADATMOSWF_OFFGRID_DN                      #
C #            L_PARTLAYER_STERM_UP                        #
C #            L_PARTLAYER_STERM_DN                        #
C #                                                        #
C ##########################################################

      SUBROUTINE UPUSER_ATMOSWF
     I   ( DO_INCLUDE_SURFACE,
     I     DO_INCLUDE_SURFEMISS,
     I     DO_INCLUDE_THERMEMISS,
     I     DO_INCLUDE_MVOUTPUT,
     I     DO_INCLUDE_DIRECTBEAM,
     I     SURFACE_FACTOR,
     I     FLUX_MULTIPLIER, 
     I     FOURIER_COMPONENT, IBEAM,
     I     VARIATION_INDEX, K_PARAMETERS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Input Fourier number and beam index
C  surface factor (2 for m = 0, 1 otherwise)
C  Flux multiplier = F/4.pi

      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UTA, UM, IUM, K, Q, NC, UT, I, IQD, IB

      DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION 
     &     L_LAYER_SOURCE      (MAX_USER_STREAMS,MAX_ATMOSWFS),
     &     L_MSCAT_LAYERSOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_FINAL_SOURCE

C  index

      K  = VARIATION_INDEX
      IB = IBEAM

C  Zero all Fourier component output here (safety)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_OUT_USERTAUS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              ATMOSWF_F(Q,K,UTA,IUM,IB,UPIDX) = ZERO
             ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Get the linearized BOA source terms (diffuse and direct)

      CALL GET_L_BOASOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      SURFACE_FACTOR,
     I      FOURIER_COMPONENT, IBEAM,
     I      K, K_PARAMETERS,
     O      L_BOA_MSSOURCE,
     O      L_BOA_DBSOURCE,
     O      L_BOA_THTONLY_SOURCE )

C  start the recursion

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + 
     &                             L_BOA_DBSOURCE(UM,Q)
          ENDDO
        ENDDO

C  Save the BOA source terms if present
C  otherwise zero these values (very important - bug fixed 04/20/05)

c        IF ( SAVE_LAYER_MSST ) THEN
c          IF ( DO_INCLUDE_SURFACE ) THEN
c            DO UM = LOCAL_UM_START, N_USER_STREAMS
c              DO Q = 1, K_PARAMETERS
c                ATMOSWF_MSCATBOA_STERM_F(Q,K,UM,IB) = 
c     &                     FLUX_MULTIPLIER * L_BOA_MSSOURCE(UM,Q)
c                ATMOSWF_DIRECTBOA_STERM_F(Q,K,UM,IB) = 
c     &                     FLUX_MULTIPLIER * L_BOA_DBSOURCE(UM,Q)
c              ENDDO
c            ENDDO
c          ELSE
c            DO UM = LOCAL_UM_START, N_USER_STREAMS
c              DO Q = 1, K_PARAMETERS
c                ATMOSWF_MSCATBOA_STERM_F (Q,K,UM,IB) = ZERO
c                ATMOSWF_DIRECTBOA_STERM_F(Q,K,UM,IB) = ZERO
c              ENDDO
c            ENDDO
c          ENDIF
c        ENDIF

C  end this section

      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CALL L_WHOLELAYER_STERM_UP
     I       ( IB, N, K, K_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) +
     &                  L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
                ENDDO
              ENDDO
            ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                ENDDO
              ENDDO
            ENDIF
c            IF ( SAVE_LAYER_MSST ) THEN
c              DO UM = LOCAL_UM_START, N_USER_STREAMS
c                DO Q = 1, K_PARAMETERS
c                  ATMOSWF_MSCATSTERM_F(Q,K,N,UM,IB,UPIDX) =
c     &                  L_MSCAT_LAYERSOURCE(UM,Q) * FLUX_MULTIPLIER
c                ENDDO
c              ENDDO
c            ENDIF
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADATMOSWF_OFFGRID_UP
     &       ( IB, UTA, UT, N, K, K_PARAMETERS, 
     &         DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER,
     &         L_BOA_THTONLY_SOURCE )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              DO Q = 1, K_PARAMETERS
                 ATMOSWF_F(Q,K,UTA,IQD,IB,UPIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX)
              ENDDO
            ENDDO
          ENDIF

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_UP
     I        ( IB, UT, N, K, K_PARAMETERS, 
     O          DO_INCLUDE_THERMEMISS, L_LAYER_SOURCE )
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) +
     &                L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
                  ATMOSWF_F(Q,K,UTA,IUM,IB,UPIDX) =
     &                   FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q)
                  ATMOSWF_F(Q,K,UTA,IUM,IB,UPIDX) =
     &                  FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF

          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADATMOSWF_LEVEL_UP
     *        ( IB, UTA, NLEVEL, K, K_PARAMETERS, 
     *          FLUX_MULTIPLIER, L_BOA_THTONLY_SOURCE )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              DO Q = 1, K_PARAMETERS
                 ATMOSWF_F(Q,K,UTA,IQD,IB,UPIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX)
              ENDDO
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              DO Q = 1, K_PARAMETERS
                ATMOSWF_F(Q,K,UTA,IUM,IB,UPIDX) =
     &                     FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
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

c      write(89,*)K
c      do i = 1, nstreams
c        write(89,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     & (QUAD_ATMOSWF(1,K,UTA,I,IB,upIDX),UTA=1,N_OUT_USERTAUS)
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE DNUSER_ATMOSWF
     I    ( DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_MVOUTPUT,
     I      FLUX_MULTIPLIER,
     I      FOURIER_COMPONENT, IBEAM,
     I      VARIATION_INDEX, K_PARAMETERS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local MV output control, thermal emission control

      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Input Fourier number and beam index
C  Flux multiplier = F/4.pi

      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UTA, UM, IUM, K, Q, NC, UT, I, IQD, IB
      DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_TOA_SOURCE  (MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION
     &        L_MSCAT_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_FINAL_SOURCE

C  Initialise

      K  = VARIATION_INDEX
      IB = IBEAM

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_OUT_USERTAUS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, LOCAL_UM_START 
              IUM = USEROUTPUT_INDEX(UM)
              ATMOSWF_F(Q,K,UTA,IUM,IB,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_L_TOASOURCE ( L_TOA_SOURCE, K_PARAMETERS )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = L_TOA_SOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            CALL L_WHOLELAYER_STERM_DN
     &       ( IB, N, K, K_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) +
     &                  L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
                ENDDO
              ENDDO
            ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                ENDDO
              ENDDO
            ENDIF      
c            IF ( SAVE_LAYER_MSST ) THEN
c              DO UM = LOCAL_UM_START, N_USER_STREAMS
c                DO Q = 1, K_PARAMETERS
c                  ATMOSWF_MSCATSTERM_F(Q,K,N,UM,IB,DNIDX) =
c     &                L_MSCAT_LAYERSOURCE(UM,Q) * FLUX_MULTIPLIER
c                ENDDO
c              ENDDO
c            ENDIF
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADATMOSWF_OFFGRID_DN
     &       ( IB, UTA, UT, N, K, K_PARAMETERS, 
     &         DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              DO Q = 1, K_PARAMETERS
                 ATMOSWF_F(Q,K,UTA,IQD,IB,DNIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX)
              ENDDO
            ENDDO
          ENDIF

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_DN
     I        ( IB, UT, N, K, K_PARAMETERS,
     O          DO_INCLUDE_THERMEMISS, L_LAYER_SOURCE )
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q) +
     &                L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
                ATMOSWF_F(Q,K,UTA,IUM,IB,DNIDX) =
     &                FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q)
                ATMOSWF_F(Q,K,UTA,IUM,IB,DNIDX) =
     &                FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADATMOSWF_LEVEL_DN
     *        ( IB, UTA, NLEVEL, K, K_PARAMETERS, 
     *          FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              DO Q = 1, K_PARAMETERS
                 ATMOSWF_F(Q,K,UTA,IQD,IB,DNIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX)
              ENDDO
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              DO Q = 1, K_PARAMETERS
                ATMOSWF_F(Q,K,UTA,IUM,IB,DNIDX) =
     &                       FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
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

c      write(56,*)K
c      do i = 1, nstreams
c        write(56,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     & (QUAD_ATMOSWF(1,K,UTA,I,IB,dnIDX),UTA=1,N_OUT_USERTAUS)
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE MIFLUX_ATMOSWF 
     &        ( DO_INCLUDE_DIRECTBEAM, IB, K, K_PARAMETERS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include files of linearized setupr variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine arguments
C  --------------------

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      INTEGER          IB, K, K_PARAMETERS

C  local variables
C  ----------------

      INTEGER          I, IDIR, WDIR, UTA, UT, Q, N
      DOUBLE PRECISION SUM_MI, SUM_FX, FMU0
      DOUBLE PRECISION L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

C  mean intensity and flux
C  -----------------------

C  direction loop

      DO IDIR = 1, N_DIRECTIONS

        WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

        DO UTA = 1, N_OUT_USERTAUS

C  loop over all parameters in the varying layer

          DO Q = 1, K_PARAMETERS

C  integrations

            SUM_MI = ZERO
            SUM_FX = ZERO
            DO I = 1, NSTREAMS
              SUM_MI = SUM_MI + A(I)  * QUAD_ATMOSWF(Q,K,UTA,I,IB,WDIR)
              SUM_FX = SUM_FX + AX(I) * QUAD_ATMOSWF(Q,K,UTA,I,IB,WDIR)
            ENDDO
            MINT_ATMOSWF(Q,K,UTA,IB,WDIR) = SUM_MI * HALF
            FLUX_ATMOSWF(Q,K,UTA,IB,WDIR) = SUM_FX * PI2

C  end loops

          ENDDO
        ENDDO

C  nothing to do if no solar sources

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

C  For the downward direction, add the direct beam contributions

        IF ( WDIR .EQ. DNIDX ) THEN

C  loop over all the output optical depths

          DO UTA = 1, N_OUT_USERTAUS

C  For the offgrid values

            IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
              UT = OFFGRID_UTAU_OUTINDEX(UTA)
              N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Only contributions for layers above the PI cutoff
C    L_INITIAL_TRANS is a logarithmic derivative

              IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                FMU0 = LOCAL_SZA(N,IB) * FLUX_FACTOR
                IF ( K.LE.N.OR.K.EQ.0 ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_TRANS = L_T_UTDN_MUBAR(UT,K,IB,Q) +
     &                 L_INITIAL_TRANS(N,K,IB,Q) * T_UTDN_MUBAR(UT,IB)
                    L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                    L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                    L_DIRECT_FLUX  = FMU0 * L_TRANS
                    MINT_ATMOSWF(Q,K,UTA,IB,WDIR) =
     *                  MINT_ATMOSWF(Q,K,UTA,IB,WDIR) + L_DIRECT_MEANI
                    FLUX_ATMOSWF(Q,K,UTA,IB,WDIR) =
     *                  FLUX_ATMOSWF(Q,K,UTA,IB,WDIR) + L_DIRECT_FLUX
                  ENDDO
                ENDIF
              ENDIF

C  For the on-grid balues

            ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                IF ( N.GT.0 ) THEN
                  FMU0 = LOCAL_SZA(N,IB) * FLUX_FACTOR
                  IF ( K.LE.N.OR.K.EQ.0 ) THEN
                    DO Q = 1, K_PARAMETERS
                      L_TRANS = L_T_DELT_MUBAR(N,K,IB,Q) +
     &                 L_INITIAL_TRANS(N,K,IB,Q) * T_DELT_MUBAR(N,IB)
                      L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                      L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                      L_DIRECT_FLUX  = FMU0 * L_TRANS
                      MINT_ATMOSWF(Q,K,UTA,IB,WDIR) =
     *                  MINT_ATMOSWF(Q,K,UTA,IB,WDIR) + L_DIRECT_MEANI
                      FLUX_ATMOSWF(Q,K,UTA,IB,WDIR) =
     *                  FLUX_ATMOSWF(Q,K,UTA,IB,WDIR) + L_DIRECT_FLUX
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

C  Close loops

          ENDDO
        ENDIF   

C  Continuation point for avoiding direct beam calculation

 455    CONTINUE

C  end direction loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_LEVEL_UP
     I    ( IB, UTA, NL, K, K_PARAMETERS, FLUX_MULTIPLIER,
     I      L_BOA_THTONLY_SOURCE )

C  Upwelling weighting function Fourier components at level boundary NL
C  Quadrature angles only

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IB, UTA, NL

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  linearization control

      INTEGER          K, K_PARAMETERS

C  Special case thermal transmittance - linearization

      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          N, I, I1, AA, Q, LAY
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5
      DOUBLE PRECISION FM, THELP, L_THELP

C  homogeneous and particular solution contributions SHOM and SPAR

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

      N = NL + 1
      FM = FLUX_MULTIPLIER

C  For the lowest level

      IF ( NL .EQ. NLAYERS  ) THEN

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * L_THELP
            ENDDO
          ENDDO
        ELSE

C  If this is also the layer that is varying, extra contributions

         IF ( K .EQ. NL .OR. K. EQ. 0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
                HOM2 = LCON_XVEC(I1,AA,NL) * L_T_DELT_EIGEN(AA,NL,Q)
                HOM3 =
     &              LCON(AA,NL)*T_DELT_EIGEN(AA,NL)*L_XPOS(I1,AA,NL,Q)
                HOM4 = PCON_XVEC(I1,AA,NL,Q)
                HOM5 = MCON(AA,NL) * L_XNEG(I1,AA,NL,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I1,NL,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

C  non-varying lowest layer

         ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
                HOM2 = PCON_XVEC(I1,AA,NL,Q)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WLOWER(I1,NL,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

         ENDIF

C  End scattering clause

        ENDIF

C  For other levels in the atmosphere
C  ----------------------------------

      ELSE

C  Thermal transmittance only

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
              THELP   =   BOA_THTONLY_SOURCE(I)
              DO LAY = NLAYERS, N, -1
                L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                IF ( LAY.EQ.K .OR. K.EQ.0 ) THEN
                  L_THELP = L_THELP +  L_T_WUPPER(I1,LAY,Q) / X(I)
     &                      + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
                ENDIF
                THELP = THELP * T_DELT_DISORDS(I,LAY)
     &                    + T_WUPPER(I1,LAY) / X(I)
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * L_THELP
            ENDDO
          ENDDO
        ELSE

C  If this is also the layer that is varying, extra contributions

         IF ( K .EQ. N  .OR. K. EQ. 0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q) 
                HOM2 = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
                HOM3 = MCON_XVEC(I1,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM4 = MCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XNEG(I1,AA,N,Q)
                HOM5 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WUPPER(I1,N,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

C  non-varying layer lower than the varying layer

         ELSE IF ( K.LT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q)
                HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WUPPER(I1,N,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

C  non-varying layer higher than the varying layer

         ELSE IF ( K.GT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q)
                HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FM * SHOM
            ENDDO
          ENDDO

         ENDIF

C  End clauses

        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_LEVEL_DN
     &     ( IB, UTA, NL, K, K_PARAMETERS, FLUX_MULTIPLIER )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          UTA, NL, IB

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  linearization control

      INTEGER          K, K_PARAMETERS

C  local variables
C  ---------------

      INTEGER          N, I, AA, Q, LAY
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5
      DOUBLE PRECISION THELP, L_THELP, FM

C  homogeneous and particular solution contributions SHOM and SPAR

      N = NL
      FM = FLUX_MULTIPLIER

C  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
             QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = ZERO
          ENDDO
        ENDDO

C  For other levels in the atmosphere
C  ----------------------------------

      ELSE

C  Thermal transmittance solution, build from TOA downwards
C  Scattering solution, use the Discrete Ordinate solution

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_THELP = ZERO
              THELP   = ZERO
              DO LAY = 1, N
                L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
                IF ( LAY.EQ.K .OR. K.EQ.0 ) THEN
                  L_THELP = L_THELP +  L_T_WLOWER(I,LAY,Q) / X(I)
     &                      + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
                ENDIF
                THELP = THELP * T_DELT_DISORDS(I,LAY)
     &                    + T_WLOWER(I,LAY) / X(I)
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FM * L_THELP
            ENDDO
          ENDDO
        ELSE

C  If this is also the layer that is varying, extra contributions

         IF ( K .EQ. N .OR. K.EQ.0  ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I,N,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

C  non-varying layer lower than the varying layer

         ELSE IF ( K.LT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q) 
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WLOWER(I,N,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

C  non-varying layer higher than the varying layer

         ELSE IF ( K.GT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q) 
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FM * SHOM
            ENDDO
          ENDDO

         ENDIF

C  End final clauses

        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_OFFGRID_UP
     &    ( IB, UTA, UT, N, K, K_PARAMETERS, 
     &      DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER,
     &      L_BOA_THTONLY_SOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Include files of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IB, UTA, UT, N, K, K_PARAMETERS

C  Thermal emission flag

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Special case thermal transmittance - linearization

      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          I, I1, AA, Q, LAY
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
      DOUBLE PRECISION H1, H2, H3, H4, H5, H6, FMULT, THELP, L_THELP

C  short hand

      FMULT = FLUX_MULTIPLIER

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_THELP = ZERO
            THELP     =   BOA_THTONLY_SOURCE(I)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N+1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.K .OR. K.EQ.0 ) THEN
                L_THELP = L_THELP +  L_T_WUPPER(I1,LAY,Q) / X(I)
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              THELP = THELP * T_DELT_DISORDS(I,LAY)
     &                   + T_WUPPER(I1,Q) / X(I)
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              L_THELP = L_THELP +  L_UT_T_PARTIC(I1,UT,Q) / X(I)
     &                  + THELP * L_T_DISORDS_UTUP(I,UT,Q) 
            ENDIF
            QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  ###########################################
      
C  Homogeneous
C  -----------

C  solution for N being the varying layer K

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_XVEC(I1,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
              H2 = NCON(AA,N,Q) *   XPOS(I1,AA,N)
              H3 = LCON(AA,N)   * L_XPOS(I1,AA,N,Q)
              H4 = MCON_XVEC(I1,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
              H5 = PCON(AA,N,Q) *   XNEG(I1,AA,N)
              H6 = MCON(AA,N)   * L_XNEG(I1,AA,N,Q)
              HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
              HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FMULT * SHOM
          ENDDO
        ENDDO

C  Solution for N not equal to K

      ELSE IF ( N. NE. K .AND. K.NE.0 ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
              H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
              HOM1 = T_UTDN_EIGEN(AA,UT) * H2
              HOM2 = T_UTUP_EIGEN(AA,UT) * H5
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = FMULT * SHOM
          ENDDO
        ENDDO

      ENDIF

C  Add the linearized thermal solution  (if flagged)
C    ---Only present if N = K, or K = 0 (column linearization)
C   THIS is the solution with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_UT_T_PARTIC(I1,UT,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = 
     &        QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add Linearized Classical beam solution
C  --------------------------------------

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solutions exist only for N = K or N > K

        IF ( N.GE.K .OR. K.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
                SPAR = L_WUPPER(I1,N,Q) *   T_UTDN_MUBAR(UT,IB) +
     &                   WUPPER(I1,N)   * L_T_UTDN_MUBAR(UT,K,IB,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = 
     &        QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF

C  Add Linearized Green function beam solution
C  -------------------------------------------

      ELSE

C  get the linearized Green's function multipliers

        CALL L_QUAD_GFUNCMULT
     I      ( IB, UT, N, K, K_PARAMETERS )

C  solution for N being the varying layer K

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                PAR1 = L_XPOS(I,AA,N,Q)  *   UT_GMULT_UP(AA,UT) +
     &                   XPOS(I,AA,N)    * L_UT_GMULT_UP(AA,UT,K,Q)
                PAR2 = L_XPOS(I1,AA,N,Q) *   UT_GMULT_DN(AA,UT) +
     &                   XPOS(I1,AA,N)   * L_UT_GMULT_DN(AA,UT,K,Q)
                SPAR = SPAR + PAR1 + PAR2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = 
     &           QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
            ENDDO
          ENDDO

C  Solution for N > K (N below the layer that varies)

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                PAR1 = XPOS(I,AA,N)  * L_UT_GMULT_UP(AA,UT,K,Q)
                PAR2 = XPOS(I1,AA,N) * L_UT_GMULT_DN(AA,UT,K,Q)
                SPAR = SPAR + PAR1 + PAR2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) = 
     &           QUAD_ATMOSWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
            ENDDO
          ENDDO

        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_OFFGRID_DN
     &   ( IB, UTA, UT, N, K, K_PARAMETERS, 
     &     DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Include files of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          IB, UTA, UT, N, K, K_PARAMETERS

C  Thermal emission flag

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          I, I1, AA, Q, LAY
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
      DOUBLE PRECISION H1, H2, H3, H4, H5, H6, FMULT, THELP, L_THELP

C  Short hand

      FMULT = FLUX_MULTIPLIER

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, N-1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.K .OR. K.EQ.0 ) THEN
                L_THELP = L_THELP +  L_T_WLOWER(I,LAY,Q) / X(I)
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              THELP = THELP * T_DELT_DISORDS(I,LAY)
     &                      + T_WLOWER(I,LAY) / X(I)
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              L_THELP = L_THELP +  L_UT_T_PARTIC(I,UT,Q) / X(I)
     &                  + THELP * L_T_DISORDS_UTDN(I,UT,Q) 
            ENDIF
            QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  ###########################################
      
C  Homogeneous
C  -----------

C  solution for N being the varying layer K, or K = 0 (column linearization)

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_XVEC(I,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
              H2 = NCON(AA,N,Q) *   XPOS(I,AA,N)
              H3 = LCON(AA,N)   * L_XPOS(I,AA,N,Q)
              H4 = MCON_XVEC(I,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
              H5 = PCON(AA,N,Q) *   XNEG(I,AA,N)
              H6 = MCON(AA,N)   * L_XNEG(I,AA,N,Q)
              HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
              HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FMULT * SHOM
          ENDDO
        ENDDO

C  Solution for N not equal to K

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
      
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
              H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
              HOM1 = T_UTDN_EIGEN(AA,UT) * H2
              HOM2 = T_UTUP_EIGEN(AA,UT) * H5
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = FMULT * SHOM
          ENDDO
        ENDDO

      ENDIF

C  Add the linearized thermal solution  (if flagged)
C    ---Only present if N = K, or K = 0
C   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_UT_T_PARTIC(I,UT,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = 
     &        QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add Classical beam solution
C  ---------------------------

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solutions exist only for N = K or N > K. or K = 0 (column Jacobians)

        IF ( N.GE.K .OR. K.EQ.0 ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_WUPPER(I,N,Q) *   T_UTDN_MUBAR(UT,IB) +
     &                 WUPPER(I,N)   * L_T_UTDN_MUBAR(UT,K,IB,Q)
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = 
     &           QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
            ENDDO
          ENDDO

        ENDIF


C  Add Linearized Green function beam solution
C  -------------------------------------------

      ELSE

C  get the linearized Green's function multipliers

        CALL L_QUAD_GFUNCMULT
     I      ( IB, UT, N, K, K_PARAMETERS )

C  solution for N being the varying layer K, or K = 0

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                PAR1 = L_XPOS(I1,AA,N,Q)  *   UT_GMULT_UP(AA,UT) +
     &                   XPOS(I1,AA,N)    * L_UT_GMULT_UP(AA,UT,K,Q)
                PAR2 = L_XPOS(I,AA,N,Q)   *   UT_GMULT_DN(AA,UT) +
     &                   XPOS(I,AA,N)     * L_UT_GMULT_DN(AA,UT,K,Q)
                SPAR = SPAR + PAR1 + PAR2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = 
     &           QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
            ENDDO
          ENDDO

C  Solution for N > K (N below the layer that varies)

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                PAR1 = XPOS(I1,AA,N)  * L_UT_GMULT_UP(AA,UT,K,Q)
                PAR2 = XPOS(I,AA,N)   * L_UT_GMULT_DN(AA,UT,K,Q)
                SPAR = SPAR + PAR1 + PAR2
              ENDDO
              QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) = 
     &           QUAD_ATMOSWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
            ENDDO
          ENDDO

        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE GET_L_TOASOURCE ( L_TOA_SOURCE, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Subroutine arguments
C  --------------------

      INTEGER          K_PARAMETERS
      DOUBLE PRECISION L_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          UM, Q

C  initialise TOA source function
C  ------------------------------

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_TOA_SOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

C  Finish

      END

C

        SUBROUTINE GET_L_BOASOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      SURFACE_FACTOR,
     I      FOURIER_COMPONENT, IBEAM,
     I      K, K_PARAMETERS,
     O      L_BOA_MSSOURCE,
     O      L_BOA_DBSOURCE,
     O      L_BOA_THTONLY_SOURCE )

C  Linearized Bottom of the atmosphere source term

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  albedo and Fourier/beam indices

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM

C  linearization control

      INTEGER          K, K_PARAMETERS

C  output arguments
C  ----------------

      DOUBLE PRECISION L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          M, N, J, I, UM, AA, Q, BRDF, IB, K1, LAY
      DOUBLE PRECISION L_DOWN(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION HELP  (MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION DOWN  (MAXSTREAMS)
      DOUBLE PRECISION REFLEC, L_BEAM, FAC, FACTOR_BRDF
      DOUBLE PRECISION SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

C  Starting section
C  ----------------

C  Fourier number, layer number

      M  = FOURIER_COMPONENT
      N  = NLAYERS
      IB = IBEAM

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_MSSOURCE(UM,Q) = ZERO
            L_BOA_DBSOURCE(UM,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_THTONLY_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  reflectance from surface
C  ------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN

C  1. Thermal Transmittance only
C     %%%%%%%%%%%%%%%%%%%%%%%%%%

C  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN

C  Initialise

          DO I = 1, NSTREAMS
            DOWN(I) = ZERO
            DO Q = 1, K_PARAMETERS
              L_DOWN(I,Q) = ZERO
            ENDDO
          ENDDO

C  Build

          DO LAY = 1, NLAYERS
            IF ( K.EQ.0. OR. K.EQ.LAY ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, K_PARAMETERS
                  L_DOWN(I,Q) = L_DOWN(I,Q) *   T_DELT_DISORDS(I,LAY) 
     &                          + DOWN(I)   * L_T_DELT_DISORDS(I,LAY,Q)
     &                          + L_T_WLOWER(I,LAY,Q)
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                DO Q = 1, K_PARAMETERS
                  L_DOWN(I,Q) = L_DOWN(I,Q) * T_DELT_DISORDS(I,LAY) 
                ENDDO
              ENDDO
            ENDIF
            DO I = 1, NSTREAMS
              DOWN(I) = DOWN(I)*T_DELT_DISORDS(I,LAY) + T_WLOWER(I,LAY)
            ENDDO
          ENDDO

C  Finalize

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              HELP(I,Q) = A(I) * L_DOWN(I,Q)
            ENDDO
          ENDDO

C  Continuation point for avoiding the next section

          GO TO 5432

        ENDIF

C  2. Scattering solutions
C     %%%%%%%%%%%%%%%%%%%%

C  Two cases:
C  If  K = N, this is also the layer that is varying --> Extras!
C  If  N > K with variations in layer K above N

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              L_DOWN(I,Q) = SHOM
            ENDDO
          ENDDO
        ELSE IF (K.LT.N.AND.K.NE.0) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q) 
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              L_DOWN(I,Q) = SHOM
            ENDDO
          ENDDO
        ENDIF

C  Particular integral linearization
C    Does not exist if thermal only and N > K, K = 0

        IF ( DO_SOLAR_SOURCES .OR.
     &      (DO_INCLUDE_THERMEMISS.AND.N.EQ.K) .OR.
     &      (DO_INCLUDE_THERMEMISS.AND.K.EQ.0) ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_DOWN(I,Q) = L_DOWN(I,Q) + L_WLOWER(I,N,Q)
            ENDDO
          ENDDO
        ENDIF         

C  reflectance integrand  a(j).x(j).L_DOWN(-j)

        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            HELP(I,Q) = L_DOWN(I,Q) * AX(I)
          ENDDO
        ENDDO

C  Continuation point

 5432   CONTINUE

C  reflected multiple scatter intensity at user defined-angles
C  -----------------------------------------------------------

        DO BRDF = 1, N_BRDF_KERNELS

C  .. integrate reflectance, same for all user-streams in Lambertian case

          IF ( LAMBERTIAN_KERNEL_FLAG(BRDF) ) THEN

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(BRDF)
              DO Q = 1, K_PARAMETERS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + HELP(J,Q)
                ENDDO
                REFLEC = FACTOR_BRDF * REFLEC
                IF ( DO_USER_STREAMS ) THEN
                 DO UM = LOCAL_UM_START, N_USER_STREAMS
                  L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
                 ENDDO
                ENDIF
                IF ( DO_QTHTONLY ) THEN
                 DO I = 1, NSTREAMS
                  L_BOA_THTONLY_SOURCE(I,Q) = 
     &               L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
                 ENDDO
                ENDIF
              ENDDO
            ENDIF

C  .. integrate with BRDF reflectance function at user angles 

          ELSE

            FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(BRDF)
            DO Q = 1, K_PARAMETERS
             IF ( DO_USER_STREAMS ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + HELP(J,Q) * USER_BIREFLEC(BRDF,UM,J)
                ENDDO
                REFLEC = FACTOR_BRDF * REFLEC
                L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
              ENDDO
             ENDIF
             IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + HELP(J,Q) * BIREFLEC(BRDF,I,J)
                ENDDO
                REFLEC = FACTOR_BRDF * REFLEC
                L_BOA_THTONLY_SOURCE(I,Q) =
     &               L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
             ENDIF
            ENDDO

          ENDIF

C  end loop over kernels

        ENDDO
        
C  Add direct beam if flagged
C    For K > 0, profile weighting functions
C    For K = 0, Bulk (column) weighting functions

        IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
          IF ( K .NE. 0 ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              FAC = - USER_DIRECT_BEAM(UM,IB) * DELTAU_SLANT(N,K,IB) 
              DO Q = 1, K_PARAMETERS
                L_BEAM = L_DELTAU_VERT(Q,K) * FAC
                L_BOA_DBSOURCE(UM,Q) = L_BEAM
              ENDDO
            ENDDO
          ELSE
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              FAC = - USER_DIRECT_BEAM(UM,IB) 
              DO Q = 1, K_PARAMETERS
                L_BEAM = ZERO
                DO K1 = 1, N
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K1)
     &                               * DELTAU_SLANT(N,K1,IB)
                ENDDO
                L_BOA_DBSOURCE(UM,Q) = L_BEAM * FAC
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

C  Surface emission has no variation here

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_STERM_UP
     I       ( IBEAM, GIVEN_LAYER,
     I         VARIATION_INDEX, K_PARAMETERS,
     I         DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE,
     O         L_MSCAT_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Include thermal flag

      LOGICAL          DO_INCLUDE_THERMEMISS
      
C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION
     &        L_LAYERSOURCE      (MAX_USER_STREAMS,MAX_ATMOSWFS),
     &        L_MSCAT_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          N, K, AA, UM, Q, IB
      DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, TM

C  local indices

      N  = GIVEN_LAYER
      K  = VARIATION_INDEX
      IB = IBEAM

C  Important to zero the output first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time. These quantities are always required.

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

C  Homogeneous solutions
C  =====================

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
              H3 = LCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
              H4 = MCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
              H6 = MCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
              H5 = PCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution
C     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) 
     &               + L_LAYER_TSUP_UP(UM,N,Q)*TM
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular and single scatter contributions - classical solution
C  ================================================================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WPOS(UM,N,K,Q) *   EMULT_UP(UM,N,IB) +
     &                 U_WPOS2(UM,N)    * L_EMULT_UP(UM,N,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORPOS(UM,N,Q) *   EMULT_UP(UM,N,IB) +
     &                   U_WPOS1(UM,N)     * L_EMULT_UP(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  ----------------------

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

C  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WPOS(UM,N,K,Q) *   EMULT_UP(UM,N,IB) +
     &                 U_WPOS2(UM,N)    * L_EMULT_UP(UM,N,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WPOS1(UM,N) * L_EMULT_UP(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Particular and single scatter contributions - Green's function solution
C  =======================================================================

      ELSE

C  Get the Green's function solution multipliers 
C         Only exist for N greater than or equal to K

        IF ( N.GE.K .OR. K.EQ.0 ) THEN
          CALL L_WHOLELAYER_GMULT_UP ( IB, N, K, K_PARAMETERS )
        ENDIF
      
C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N .EQ. K .OR. K.EQ.0  ) THEN

C  Add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XPOS(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                H2 = L_U_XPOS(UM,AA,N,Q) *   SGMULT_UD(AA,UM,N)
                H3 =   U_XNEG(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                H4 = L_U_XNEG(UM,AA,N,Q) *   SGMULT_UU(AA,UM,N)
                SPAR = SPAR + H1 + H2 + H3 + H4
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORPOS(UM,N,Q) *   EMULT_UP(UM,N,IB) +
     &                    U_WPOS1(UM,N)    * L_EMULT_UP(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  ----------------------

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

C  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 = U_XPOS(UM,AA,N) * L_SGMULT_DN(AA,UM,Q)
                H3 = U_XNEG(UM,AA,N) * L_SGMULT_UP(AA,UM,Q)
                SPAR = SPAR + H1 + H3
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WPOS1(UM,N) * L_EMULT_UP(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
             ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF

C  If operating in Ms-mode only, copy multiple scatter term for MSST
C  -----------------------------------------------------------------

c      IF ( DO_MSMODE_LIDORT ) THEN
c        IF ( SAVE_LAYER_MSST ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            DO Q = 1, K_PARAMETERS
c              L_MSCAT_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)
c            ENDDO
c          ENDDO
c        ENDIF
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_STERM_DN
     I       ( IBEAM, GIVEN_LAYER,
     I         VARIATION_INDEX, K_PARAMETERS,
     I         DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE,
     O         L_MSCAT_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Include thermal flag

      LOGICAL          DO_INCLUDE_THERMEMISS
      
C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION
     &        L_LAYERSOURCE      (MAX_USER_STREAMS,MAX_ATMOSWFS),
     &        L_MSCAT_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          N, K, AA, UM, Q, IB
      DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, TM

C  local indices

      N  = GIVEN_LAYER
      K  = VARIATION_INDEX
      IB = IBEAM

C  Important to zero the output first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time. These quantities are always required.

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

C  Homogeneous solutions
C  =====================

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
              H3 = LCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
              H4 = MCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
              H6 = MCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
            SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

      ELSE IF ( K.NE.0 .AND. N.NE.K ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
              H5 = PCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution
C     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) 
     &               + L_LAYER_TSUP_DN(UM,N,Q)*TM
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular and single scatter contributions - classical solution
C  ================================================================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WNEG (UM,N,K,Q) *   EMULT_DN(UM,N,IB) +
     &                 U_WNEG2(UM,N)     * L_EMULT_DN(UM,N,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORNEG(UM,N,Q) *   EMULT_DN(UM,N,IB) +
     &                     U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  ----------------------

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

C  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WNEG(UM,N,K,Q) *   EMULT_DN(UM,N,IB) +
     &                 U_WNEG2(UM,N)    * L_EMULT_DN(UM,N,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WNEG1(UM,N) * L_EMULT_DN(UM,N,K,IB,Q)
                 L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Particular and single scatter contributions - Green's function solution
C  =======================================================================

      ELSE

C  Get the Green's function solution multipliers 
C  Only exist for N greater than or equal to K, or K = 0 (Bulk)

        IF ( N.GE.K .OR. K.EQ.0 ) THEN
          CALL L_WHOLELAYER_GMULT_DN ( IB, N, K, K_PARAMETERS )
        ENDIF

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  Add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                H2 = L_U_XNEG(UM,AA,N,Q) *   SGMULT_DD(AA,UM,N)
                H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                H4 = L_U_XPOS(UM,AA,N,Q) *   SGMULT_DU(AA,UM,N)
                SPAR = SPAR + H1 + H2 + H3 + H4
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORNEG(UM,N,Q) *   EMULT_DN(UM,N,IB) +
     &                     U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  ----------------------

C  add particular solution

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                SPAR = SPAR + H1 + H3
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WNEG1(UM,N)  * L_EMULT_DN(UM,N,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF

C  If operating in Ms-mode only, copy multiple scatter term for MSST
C  -----------------------------------------------------------------

c      IF ( DO_MSMODE_LIDORT ) THEN
c        IF ( SAVE_LAYER_MSST ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            DO Q = 1, K_PARAMETERS
c              L_MSCAT_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)
c            ENDDO
c          ENDDO
c        ENDIF
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_STERM_UP
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER,
     I         VARIATION_INDEX, K_PARAMETERS,
     I         DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  layer and beam index

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Include thermal flag

      LOGICAL          DO_INCLUDE_THERMEMISS
      
C  Subroutine output arguments
C  ---------------------------

C  output linearized layer source term

      DOUBLE PRECISION L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          N, K, AA, UM, Q, UT, IB
      DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, TM

C  local indices

      N  = GIVEN_LAYER
      K  = VARIATION_INDEX
      UT = OFFGRID_INDEX
      IB = IBEAM

C  Important to zero the output first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time. These quantities are always required.

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_UD(AA,UM,UT,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_UD(AA,UM,UT)
              H3 = LCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
              H3 = H3 * UT_HMULT_UD(AA,UM,UT)
              H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_UU(AA,UM,UT,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_UU(AA,UM,UT)
              H6 = MCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
              H6 = H6 * UT_HMULT_UU(AA,UM,UT)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_UD(AA,UM,UT)
              H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_UU(AA,UM,UT)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution
C     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) 
     &               + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN 

C  Partial layer source function ( SS/Particular, Classical solution )
C  ===================================================================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WPOS (UM,N,N,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &                 U_WPOS2(UM,N)     * L_UT_EMULT_UP(UM,UT,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORPOS(UM,N,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &                     U_WPOS1(UM,N)   * L_UT_EMULT_UP(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  -----------------------

        ELSE  IF ( N .GT. K .AND. K.NE.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WPOS (UM,N,K,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &                 U_WPOS2(UM,N)     * L_UT_EMULT_UP(UM,UT,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WPOS1(UM,N) * L_UT_EMULT_UP(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Partial layer source function ( SS/Particular, Green's solution )
C  =================================================================

      ELSE

C  Get the Green's function solution multipliers
C  Only exist for N greater than or equal to K, or K = 0

        IF ( N.GE.K .OR. K.EQ.0 ) THEN
          CALL L_PARTLAYER_GMULT_UP ( IB, UT, N, K, K_PARAMETERS )
        ENDIF

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XPOS(UM,AA,N)   *  L_SGMULT_DN(AA,UM,Q)
                H2 = L_U_XPOS(UM,AA,N,Q) * UT_SGMULT_UD(AA,UM,UT)
                H3 =   U_XNEG(UM,AA,N)   *  L_SGMULT_UP(AA,UM,Q)
                H4 = L_U_XNEG(UM,AA,N,Q) * UT_SGMULT_UU(AA,UM,UT)
                SPAR = SPAR + H1 + H2 + H3 + H4
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORPOS(UM,N,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &                     U_WPOS1(UM,N)   * L_UT_EMULT_UP(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  ----------------------

        ELSE IF ( N .GT. K .AND. K. NE.0) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                  H1 = U_XPOS(UM,AA,N) * L_SGMULT_DN(AA,UM,Q)
                  H3 = U_XNEG(UM,AA,N) * L_SGMULT_UP(AA,UM,Q)
                SPAR = SPAR + H1 + H3
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WPOS1(UM,N)  * L_UT_EMULT_UP(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_STERM_DN
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER, 
     I         VARIATION_INDEX, K_PARAMETERS,
     I         DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  layer and beam indices

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Include thermal flag

      LOGICAL          DO_INCLUDE_THERMEMISS
      
C  Subroutine output arguments
C  ---------------------------

C  output linearized layer source term

      DOUBLE PRECISION L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          N, K, AA, UM, Q, UT, IB
      DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, TM

C  local indices

      N  = GIVEN_LAYER
      K  = VARIATION_INDEX
      UT = OFFGRID_INDEX
      IB = IBEAM

C  Important to zero the output first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time. These quantities are always required.

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_DD(AA,UM,UT,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_DD(AA,UM,UT)
              H3 = LCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
              H3 = H3 * UT_HMULT_DD(AA,UM,UT)
              H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_DU(AA,UM,UT,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_DU(AA,UM,UT)
              H6 = MCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
              H6 = H6 * UT_HMULT_DU(AA,UM,UT)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_DD(AA,UM,UT)
              H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_DU(AA,UM,UT)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution
C     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) 
     &               + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Partial layer source function ( SS/Particular, Classical solution )
C  ===================================================================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WNEG (UM,N,N,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &                 U_WNEG2(UM,N)     * L_UT_EMULT_DN(UM,UT,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORNEG(UM,N,Q) * UT_EMULT_DN(UM,UT,IB) +
     &                     U_WNEG1(UM,N)   * L_UT_EMULT_DN(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  -----------------------

        ELSE IF ( N .GT. K .AND. K.NE.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = L_U_WNEG (UM,N,K,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &                 U_WNEG2(UM,N)     * L_UT_EMULT_DN(UM,UT,K,IB,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WNEG1(UM,N) * L_UT_EMULT_DN(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Partial layer source function ( SS/Particular, Green's solution )
C  =================================================================

      ELSE

C  Get the Green's function solution multipliers
C  Only exist for N greater than or equal to K, or K = 0 (bulk)

        IF ( N.GE.K .OR. K.EQ.0 ) THEN
          CALL L_PARTLAYER_GMULT_DN ( IB, UT, N, K, K_PARAMETERS )
        ENDIF

C  Special case when N = K, or K = 0 (bulk)
C  ----------------------------------------

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XNEG(UM,AA,N)   *  L_SGMULT_DN(AA,UM,Q)
                H2 = L_U_XNEG(UM,AA,N,Q) * UT_SGMULT_DD(AA,UM,UT)
                H3 =   U_XPOS(UM,AA,N)   *  L_SGMULT_UP(AA,UM,Q)
                H4 = L_U_XPOS(UM,AA,N,Q) * UT_SGMULT_DU(AA,UM,UT)
                SPAR = SPAR + H1 + H2 + H3 + H4
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = L_U_WFORNEG(UM,N,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &                     U_WNEG1(UM,N)   * L_UT_EMULT_DN(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

C  Other cases when N > K
C  -----------------------

        ELSE IF ( N .GT. K .AND. K.NE.0 ) THEN

C  particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SPAR = ZERO
              DO AA = 1, NSTREAMS
                H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                SPAR = SPAR + H1 + H3
              ENDDO
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
            ENDDO
          ENDDO

C  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR = U_WNEG1(UM,N)  * L_UT_EMULT_DN(UM,UT,K,IB,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF

C  Finish

      RETURN
      END

