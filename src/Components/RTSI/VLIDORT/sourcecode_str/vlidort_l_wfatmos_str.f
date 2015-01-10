C ###############################################################
C #                                                             #
C #                    THE VLIDORT  MODEL                       #
C #                                                             #
C #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
C #  -          --         -        -        -         -        #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, Inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R                  #
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

C ##########################################################
C #                                                        #
C # Subroutines in this Module                             #
C #                                                        #
C #     Top level routines--------------                   #
C #                                                        #
C #            VLIDORT_UPUSER_ATMOSWF                      #
C #            VLIDORT_DNUSER_ATMOSWF                      #
C #                                                        #
C #     Post-processing at user angles --------            #
C #                                                        #
C #            GET_L_TOASOURCE                             #
C #            L_BOA_LAMBERTIAN_SOURCE                     #
C #                                                        #
C #            L_WHOLELAYER_STERM_UP                       #
C #            L_WHOLELAYER_STERM_DN                       #
C #            L_PARTLAYER_STERM_UP                        #
C #            L_PARTLAYER_STERM_DN                        #
C #                                                        #
C ##########################################################

      SUBROUTINE VLIDORT_UPUSER_ATMOSWF
     I   ( DO_INCLUDE_SURFACE,
     I     DO_INCLUDE_THERMEMISS,
     I     DO_INCLUDE_MVOUTPUT,
     I     DO_INCLUDE_DIRECTBEAM,
     I     SURFACE_FACTOR,    FLUX_MULTIPLIER,
     I     FOURIER_COMPONENT, IBEAM,
     I     LAYER_TO_VARY, NV_PARAMETERS,
     I     L_GET_BOA_SOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'
      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  External function

      EXTERNAL         L_GET_BOA_SOURCE

C  Input Fourier number and beam index
C  surface factor (2 for m = 0, 1 otherwise)
C  Flux multiplier = F/4.pi

      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      LOGICAL          SFLAG
      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, NV, Q, NC, UT, IB

      DOUBLE PRECISION L_CUMUL_SOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_BOA_MSSOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION L_BOA_DBSOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_LAYER_SOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_FINAL_SOURCE

C  index

      NV = LAYER_TO_VARY
      IB = IBEAM

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            DO UM = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

      IF ( DO_USER_STREAMS ) THEN

C  Get the linearized BOA source terms (diffuse and direct)

        CALL L_GET_BOA_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      SURFACE_FACTOR,
     I      FOURIER_COMPONENT, IBEAM,
     I      NV, NV_PARAMETERS,
     O      L_BOA_MSSOURCE,
     O      L_BOA_DBSOURCE )

C  Set the cumuluative source term equal to the BOA sum

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
            L_CUMUL_SOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + 
     &                                L_BOA_DBSOURCE(UM,O1,Q)
          ENDDO
         ENDDO
        ENDDO

C  Debug

        IF ( DO_DEBUG_WRITE ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
            IF ( FOURIER_COMPONENT.EQ.0.and.q.eq.2.and.nv.eq.2)then
              if (um.eq.1.and.o1.eq.1)write(83,*)UM,O1,
     &         L_BOA_MSSOURCE(UM,O1,Q), L_BOA_DBSOURCE(UM,O1,Q)
            endif
            ENDDO
           ENDDO
          ENDDO
        ENDIF

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

      DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = NLAYERS + 1 - N
            CALL L_WHOLELAYER_STERM_UP
     I       ( FOURIER_COMPONENT, IB, N, SFLAG,
     I         NV, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYER_SOURCE )
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q)
     &               +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
     &               + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_UP(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q)
     &               +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_UP
     I   ( IB, UT, N, SFLAG, NV, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O     L_LAYER_SOURCE )
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q)
     &              +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
     &              + L_T_UTUP_USERM(UT,UM,Q) *   CUMSOURCE_UP(UM,O1,NC)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) =
     &                   FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q)
     &              +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) =
     &                  FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  L_FINAL_SOURCE = 
     &                   FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
!                  if (DABS(L_FINAL_SOURCE).GT.1.0d-12 ) then
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = L_FINAL_SOURCE
!                  endif
!      if(fourier_component.eq.3.and.q.eq.3)
!     &            write(*,*)uta,um,q,L_CUMUL_SOURCE(UM,1,Q)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  debug

c        IF ( DO_DEBUG_WRITE ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            DO O1 = 1, NSTOKES
c              DO Q = 1, NV_PARAMETERS
c                if(q.eq.1.and.nv.lt.10) write(70,'(5i3,1p2e20.10)')
c     &         fourier_component,uta,nv,um,o1,       
c     &           STOKES_F(UTA,UM,IBEAM,O1,UPIDX),
c     &           ATMOSWF_F(Q,NV,UTA,UM,IBEAM,O1,UPIDX)
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF
c         if (nv.eq.10)pause

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_DNUSER_ATMOSWF
     I    ( DO_INCLUDE_THERMEMISS,
     I      FLUX_MULTIPLIER,
     I      FOURIER_COMPONENT, IBEAM,
     I      LAYER_TO_VARY, NV_PARAMETERS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  thermal emission control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Input Fourier number and beam index
C  Flux multiplier = F/4.pi

      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  local variables
C  ---------------

      LOGICAL          SFLAG
      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, NV, Q, NC, UT, IB, M

      DOUBLE PRECISION L_CUMUL_SOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_TOA_SOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_LAYER_SOURCE
     &        ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION L_FINAL_SOURCE

C  Initialise

      NV = LAYER_TO_VARY
      IB = IBEAM
      M  = FOURIER_COMPONENT

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            DO UM = 1, LOCAL_UM_START 
              DO O1 = 1, NSTOKES
                ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_L_TOASOURCE ( L_TOA_SOURCE, NV_PARAMETERS )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(UM,O1,Q) = L_TOA_SOURCE(UM,O1,Q)
            ENDDO
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

      DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = N
            CALL L_WHOLELAYER_STERM_DN
     &       ( IB, N, SFLAG, NV, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYER_SOURCE )

            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q)
     &               +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
     &               + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_DN(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q)
     &               +   T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF      
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_DN
     I    ( IB, UT, N, SFLAG, NV, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O      L_LAYER_SOURCE )
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q)
     &              +   T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
     &              + L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,O1,NC)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) =
     &                FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q)
     &              +   T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) =
     &                FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) =
     &                    FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
                ENDDO
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

C  Finish

      RETURN
      END

C

      SUBROUTINE GET_L_TOASOURCE ( L_TOA_SOURCE, NV_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Subroutine arguments
C  --------------------

      INTEGER          NV_PARAMETERS
      DOUBLE PRECISION L_TOA_SOURCE
     &           (MAX_USER_STREAMS,MAXSTOKES,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          UM, Q, O1

C  initialise TOA source function
C  ------------------------------

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, NV_PARAMETERS
          DO O1 = 1, NSTOKES
            L_TOA_SOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BOA_LAMBERTIAN_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      SURFACE_FACTOR,
     I      FOURIER_COMPONENT, IBEAM,
     I      NV, NV_PARAMETERS,
     O      L_BOA_MSSOURCE,
     O      L_BOA_DBSOURCE )

C  Linearized Bottom of the atmosphere source term

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  thermal emission control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  albedo and Fourier/beam indices

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM

C  linearization control

      INTEGER          NV, NV_PARAMETERS

C  output arguments
C  ----------------

      DOUBLE PRECISION L_BOA_MSSOURCE
     &          (MAX_USER_STREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_DBSOURCE
     &          (MAX_USER_STREAMS,MAXSTOKES,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          M, N, NN, J, I, UM, Q, IB, O1
      INTEGER          K, KO1, K0, K1, K2

      LOGICAL          DO_QTHTONLY
      DOUBLE PRECISION DOWN   (MAXSTREAMS)
      DOUBLE PRECISION L_DOWN (MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION REFLEC, L_BEAM, FAC, KMULT
      DOUBLE PRECISION SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION LXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION LXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION LXR2, NXR2, LLXR2

C  Starting section
C  ----------------

C  Fourier number, layer number

      M   = FOURIER_COMPONENT
      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_BOA_MSSOURCE(UM,O1,Q) = ZERO
              L_BOA_DBSOURCE(UM,O1,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_BOA_THTONLY_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  reflectance from surface
C  ------------------------

C  Only want the (1,1) component

      O1 = 1

C  reflectance integrand  a(j).x(j).L(downwelling)(-j)

      IF ( DO_INCLUDE_SURFACE ) THEN

C  1. Thermal Transmittance only
C     %%%%%%%%%%%%%%%%%%%%%%%%%%

C  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN

C  Initialise

          DO I = 1, NSTREAMS
            DOWN(I) = ZERO
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,Q) = ZERO
            ENDDO
          ENDDO

C  Build

          DO NN = 1, NLAYERS
            IF ( NV.EQ.0. OR. NV.EQ.NN ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,Q) = L_DOWN(I,Q) *   T_DELT_DISORDS(I,NN) 
     &                          + DOWN(I)   * L_T_DELT_DISORDS(I,NN,Q)
     &                          + L_T_WLOWER(I,NN,Q)
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,Q) = L_DOWN(I,Q) * T_DELT_DISORDS(I,NN) 
                ENDDO
              ENDDO
            ENDIF
            DO I = 1, NSTREAMS
              DOWN(I) = DOWN(I)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
            ENDDO
          ENDDO

C  Finalize

          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,Q) = QUAD_WEIGHTS(I) * L_DOWN(I,Q)
            ENDDO
          ENDDO

C  Continuation point for avoiding the next section

          GO TO 5432

        ENDIF

C  2. Scattering solutions
C     %%%%%%%%%%%%%%%%%%%%

C  Two cases:
C  (a) If  NV = N, this is also the layer that is varying --> Extras!
C  (b) If  N > NV with variations in layer NV above N

        IF ( NV .EQ. N .OR. NV.EQ.0 ) THEN

C  stream and stokes loops

          DO I = 1, NSTREAMS

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

             SHOM_R = ZERO
             DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
              LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
              MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
              NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N)
     &                     +  LXR   * L_T_DELT_EIGEN(K,N,Q)
              HOM2 = PXR + MLXR
              SHOM_R = SHOM_R + HOM1 + HOM2
             ENDDO

C  Complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N)
     &                - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
              LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N)
     &                + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)
              LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
     &                 - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
              LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
     &                 + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
              MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q)
     &                 - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
              HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N)
     &                 - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
              HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q)
     &                           - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
              HOM3CR = PXR1 + MLXR1
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR 
             ENDDO

C  real part 

             SHOM = SHOM_R + SHOM_CR
             L_DOWN(I,Q) = SHOM

C  End loops over Q and I

            ENDDO
          ENDDO

C  otherwise the varying layer is above the boundary layer

        ELSE IF ( NV.LT.N .AND.NV.NE.0 ) THEN

C  stream and stokes loops

          DO I = 1, NSTREAMS

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

             SHOM_R = ZERO
             DO K = 1, K_REAL(N)
              NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              HOM1 = NXR * T_DELT_EIGEN(K,N)
              HOM2 = PXR
              SHOM_R = SHOM_R + HOM1 + HOM2 
             ENDDO

C  Complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N)
     &                 - NXR2 * T_DELT_EIGEN(K2,N)
              HOM2CR = PXR1
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
             ENDDO

C  real part homogeneous solution

             SHOM = SHOM_R + SHOM_CR
             L_DOWN(I,Q) = SHOM

C  End loops over Q and I

            ENDDO
          ENDDO

C  End cases

        ENDIF

C  Particular integral linearization
C    Does not exist if thermal only and N > K, K = 0

        IF ( DO_SOLAR_SOURCES .OR.
     &      (DO_INCLUDE_THERMEMISS.AND.N.EQ.NV) .OR.
     &      (DO_INCLUDE_THERMEMISS.AND.NV.EQ.0) ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,Q) = L_DOWN(I,Q) + L_WLOWER(I,O1,N,Q)
            ENDDO
          ENDDO
        ENDIF        
 
C  reflectance integrand  a(j).x(j).L_DOWN(-j)

        DO Q = 1, NV_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_STRMWTS(I)
          ENDDO
        ENDDO

C  Continuation point

 5432   CONTINUE

C  ###### Lambertian reflectance (same for all user-streams)

        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
          O1 = 1
          DO Q = 1, NV_PARAMETERS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + L_DOWN(J,Q)
            ENDDO
            REFLEC = KMULT * REFLEC
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q)
     &                             + REFLEC
            ENDDO
            IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
                 L_BOA_THTONLY_SOURCE(I,Q) = 
     &               L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Add direct beam if flagged
C    For NV > 0, profile weighting functions
C    For NV = 0, Bulk (column) weighting functions

        IF ( DO_INCLUDE_DIRECTBEAM.AND..NOT.DO_DBCORRECTION ) THEN
          IF ( NV .NE. 0 ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                FAC = - USER_DIRECT_BEAM(UM,IB,O1)*DELTAU_SLANT(N,NV,IB)
                DO Q = 1, NV_PARAMETERS
                  L_BEAM = L_DELTAU_VERT(Q,NV) * FAC
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BEAM
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                FAC = - USER_DIRECT_BEAM(UM,IB,O1) 
                DO Q = 1, NV_PARAMETERS
                  L_BEAM = ZERO
                  DO K1 = 1, N
                    L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K1)
     &                               * DELTAU_SLANT(N,K1,IB)
                   ENDDO
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BEAM * FAC
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

C  End inclusion of surface terms

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_STERM_UP
     I       ( M, IBEAM, GIVEN_LAYER, SOURCETERM_FLAG,
     I         LAYER_TO_VARY, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Fourier term (debug only)

      INTEGER          M

C  Flag

      LOGICAL          SOURCETERM_FLAG

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  Thermal control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Subroutine output argument
C  --------------------------

      DOUBLE PRECISION L_LAYERSOURCE 
     &          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

C  local variables
C  ---------------

      INTEGER          N, NV, UM, O1, Q, IB
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

C   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)  = ZERO
          ENDDO
       ENDDO
      ENDDO

!  Need to go on if thermal transmittances only
!      IF ( .NOT. SOURCETERM_FLAG ) RETURN

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Homogeneous solutions
C  =====================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N .EQ. NV .OR. NV.EQ.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = ( NUXR + LLUXR ) *   HMULT_2(K,UM,N)
     &                   +  LUXR    * L_HMULT_2(K,UM,N,Q)
              H2 = ( PUXR + MLUXR ) *   HMULT_1(K,UM,N)
     &                   +  MUXR    * L_HMULT_1(K,UM,N,Q)
              SHOM_R = SHOM_R + H1 + H2 
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) -
     &                 LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) +
     &                 LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) -
     &                 MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) +
     &                 MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) -
     &                 LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) +
     &                 LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) -
     &                 MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) +
     &                 MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * HMULT_2(K1,UM,N)
     &              - ( NUXR2 + LLUXR2 ) * HMULT_2(K2,UM,N)
     &                    +  LUXR1       * L_HMULT_2(K1,UM,N,Q)
     &                    -  LUXR2       * L_HMULT_2(K2,UM,N,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * HMULT_1(K1,UM,N)
     &              - ( PUXR2 + MLUXR2 ) * HMULT_1(K2,UM,N)
     &                    +  MUXR1       * L_HMULT_1(K1,UM,N,Q)
     &                    -  MUXR2       * L_HMULT_1(K2,UM,N,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR
           ENDDO
          ENDDO
        ENDDO

C  Other cases when N not equal to NV (only variation of Integ-Cons)
C  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = NUXR *   HMULT_2(K,UM,N)
              H2 = PUXR *   HMULT_1(K,UM,N)
              SHOM_R = SHOM_R + H1 + H2 
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    NUXR1 * HMULT_2(K1,UM,N)
     &              - NUXR2 * HMULT_2(K2,UM,N)
              H2 =    PUXR1 * HMULT_1(K1,UM,N)
     &              - PUXR2 * HMULT_1(K2,UM,N)

              SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     -----Modulus 1 if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) 
     &               + L_LAYER_TSUP_UP(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular and single scatter contributions
C  ===========================================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) +
     &              UPAR_UP_2(UM,O1,N)      * L_EMULT_UP(UM,N,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR = L_UPAR_UP_1(UM,O1,N,Q) *   EMULT_UP(UM,N,IB) +
     &                UPAR_UP_1(UM,O1,N)   * L_EMULT_UP(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  Other cases when N > NV
C  -----------------------

      ELSE IF ( N .GT. NV .AND. NV.NE.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) +
     &              UPAR_UP_2(UM,O1,N)      * L_EMULT_UP(UM,N,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR =  UPAR_UP_1(UM,O1,N) * L_EMULT_UP(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End layer clause

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_STERM_DN
     I       ( IBEAM, GIVEN_LAYER, SOURCETERM_FLAG,
     I         LAYER_TO_VARY, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flag

      LOGICAL          SOURCETERM_FLAG

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  Thermal control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Subroutine output argument
C  --------------------------

      DOUBLE PRECISION L_LAYERSOURCE 
     &          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

C  local variables
C  ---------------

      INTEGER          N, NV, UM, O1, Q, IB
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

C   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Need to go on if thermal transmittances only
!      IF ( .NOT. SOURCETERM_FLAG ) RETURN

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Homogeneous solutions
C  =====================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = ( NUXR + LLUXR ) *   HMULT_1(K,UM,N)
     &                   +  LUXR    * L_HMULT_1(K,UM,N,Q)
              H2 = ( PUXR + MLUXR ) *   HMULT_2(K,UM,N)
     &                   +  MUXR    * L_HMULT_2(K,UM,N,Q)
              SHOM_R = SHOM_R + H1 + H2 
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) -
     &                 LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) +
     &                 LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) -
     &                 MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) +
     &                 MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) -
     &                 LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) +
     &                 LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) -
     &                 MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) +
     &                 MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * HMULT_1(K1,UM,N)
     &              - ( NUXR2 + LLUXR2 ) * HMULT_1(K2,UM,N)
     &                    +  LUXR1       * L_HMULT_1(K1,UM,N,Q)
     &                    -  LUXR2       * L_HMULT_1(K2,UM,N,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * HMULT_2(K1,UM,N)
     &              - ( PUXR2 + MLUXR2 ) * HMULT_2(K2,UM,N)
     &                    +  MUXR1       * L_HMULT_2(K1,UM,N,Q)
     &                    -  MUXR2       * L_HMULT_2(K2,UM,N,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  Other cases when N not equal to NV (only variation of Integ-Cons)
C  ----------------------------------

      ELSE IF ( NV.NE.0 .AND. N.NE.NV ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = NUXR * HMULT_1(K,UM,N)
              H2 = PUXR * HMULT_2(K,UM,N)
              SHOM_R = SHOM_R + H1 + H2 
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =     NUXR1 * HMULT_1(K1,UM,N)
     &              -  NUXR2 * HMULT_1(K2,UM,N)
              H2 =     PUXR1 * HMULT_2(K1,UM,N)
     &              -  PUXR2 * HMULT_2(K2,UM,N)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  End layer clause

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- Modulus 1 if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) 
     &               + L_LAYER_TSUP_DN(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular and single scatter contributions
C  ===========================================

C  Special case when N = NV, or NV = 0 (bulk)
C  -----------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) +
     &              UPAR_DN_2(UM,O1,N)      * L_EMULT_DN(UM,N,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR = L_UPAR_DN_1(UM,O1,N,Q) *   EMULT_DN(UM,N,IB) +
     &                UPAR_DN_1(UM,O1,N)   * L_EMULT_DN(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  Other cases when N > NV
C  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
          SPAR = L_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) +
     &             UPAR_DN_2(UM,O1,N)      * L_EMULT_DN(UM,N,NV,IB,Q)
          L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR =  UPAR_DN_1(UM,O1,N) * L_EMULT_DN(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End layer clause

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_STERM_UP
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER, SOURCETERM_FLAG,
     I         LAYER_TO_VARY, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flag

      LOGICAL          SOURCETERM_FLAG

C  layer and beam index

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  Thermal control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Subroutine output argument
C  --------------------------

      DOUBLE PRECISION L_LAYERSOURCE 
     &          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

C  local variables
C  ---------------

      INTEGER          N, NV, UM, O1, Q, IB, UT
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

C   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Need to go on if thermal transmittances only
!      IF ( .NOT. SOURCETERM_FLAG ) RETURN

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB = IBEAM

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)

              LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)

              H1 = ( NUXR + LLUXR ) *   UT_HMULT_UD(K,UM,UT)
     &                   +  LUXR    * L_UT_HMULT_UD(K,UM,UT,Q)
              H2 = ( PUXR + MLUXR ) *   UT_HMULT_UU(K,UM,UT)
     &                   +  MUXR    * L_UT_HMULT_UU(K,UM,UT,Q)
              SHOM_R = SHOM_R + H1 + H2 

            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) -
     &                 LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) +
     &                 LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) -
     &                 MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) +
     &                 MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) -
     &                 LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) +
     &                 LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) -
     &                 MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) +
     &                 MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_UD(K1,UM,UT)
     &              - ( NUXR2 + LLUXR2 ) * UT_HMULT_UD(K2,UM,UT)
     &                    +  LUXR1       * L_UT_HMULT_UD(K1,UM,UT,Q)
     &                    -  LUXR2       * L_UT_HMULT_UD(K2,UM,UT,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_UU(K1,UM,UT)
     &              - ( PUXR2 + MLUXR2 ) * UT_HMULT_UU(K2,UM,UT)
     &                    +  MUXR1       * L_UT_HMULT_UU(K1,UM,UT,Q)
     &                    -  MUXR2       * L_UT_HMULT_UU(K2,UM,UT,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  Other cases when N not equal to NV (only variation of Integ-Cons)
C  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = NUXR *   UT_HMULT_UD(K,UM,UT)
              H2 = PUXR *   UT_HMULT_UU(K,UM,UT)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    NUXR1 * UT_HMULT_UD(K1,UM,UT)
     &              - NUXR2 * UT_HMULT_UD(K2,UM,UT)
              H2 =    PUXR1 * UT_HMULT_UU(K1,UM,UT)
     &              - PUXR2 * UT_HMULT_UU(K2,UM,UT)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  End layer clause

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) 
     &               + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN 

C  Particular and single scatter contributions
C  ===========================================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &            UPAR_UP_2(UM,O1,N)      * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP = L_UPAR_UP_1(UM,O1,N,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &              UPAR_UP_1(UM,O1,N)   * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  Other cases when N > NV
C  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) +
     &            UPAR_UP_2(UM,O1,N)      * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP =  UPAR_UP_1(UM,O1,N) * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End layer clause

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_STERM_DN
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER, SOURCETERM_FLAG,
     I         LAYER_TO_VARY, NV_PARAMETERS, DO_INCLUDE_THERMEMISS,
     O         L_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flag

      LOGICAL          SOURCETERM_FLAG

C  layer and beam indices

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  Thermal control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Subroutine output argument
C  --------------------------

      DOUBLE PRECISION L_LAYERSOURCE 
     &          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

C  local variables
C  ---------------

      INTEGER          N, NV, UM, O1, Q, IB, UT
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

C   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Need to go on if thermal transmittances only
!      IF ( .NOT. SOURCETERM_FLAG ) RETURN

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB  = IBEAM

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)

              LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)

              H1 = ( NUXR + LLUXR ) *   UT_HMULT_DD(K,UM,UT)
     &                   +  LUXR    * L_UT_HMULT_DD(K,UM,UT,Q)
              H2 = ( PUXR + MLUXR ) *   UT_HMULT_DU(K,UM,UT)
     &                   +  MUXR    * L_UT_HMULT_DU(K,UM,UT,Q)
              SHOM_R = SHOM_R + H1 + H2 

            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) -
     &                 LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) +
     &                 LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) -
     &                 MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) +
     &                 MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) -
     &                 LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) +
     &                 LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) -
     &                 MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) +
     &                 MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_DD(K1,UM,UT)
     &              - ( NUXR2 + LLUXR2 ) * UT_HMULT_DD(K2,UM,UT)
     &                    +  LUXR1       * L_UT_HMULT_DD(K1,UM,UT,Q)
     &                    -  LUXR2       * L_UT_HMULT_DD(K2,UM,UT,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_DU(K1,UM,UT)
     &              - ( PUXR2 + MLUXR2 ) * UT_HMULT_DU(K2,UM,UT)
     &                    +  MUXR1       * L_UT_HMULT_DU(K1,UM,UT,Q)
     &                    -  MUXR2       * L_UT_HMULT_DU(K2,UM,UT,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  Other cases when N not equal to NV (only variation of Integ-Cons)
C  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

C  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

C  parameter loop

           DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = NUXR *   UT_HMULT_DD(K,UM,UT)
              H2 = PUXR *   UT_HMULT_DU(K,UM,UT)
              SHOM_R = SHOM_R + H1 + H2 
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) -
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) +
     &                 NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) -
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) +
     &                 PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    NUXR1 * UT_HMULT_DD(K1,UM,UT)
     &              - NUXR2 * UT_HMULT_DD(K2,UM,UT)
              H2 =    PUXR1 * UT_HMULT_DU(K1,UM,UT)
     &              - PUXR2 * UT_HMULT_DU(K2,UM,UT)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

C  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

C  End layer clause

      ENDIF

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
C     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) 
     &               + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular and single scatter contributions
C  ===========================================

C  Special case when N = NV, or NV = 0 (bulk)
C  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &            UPAR_DN_2(UM,O1,N)      * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
!          if (ut.eq.1)write(*,*)um,q,
!     &        L_UPAR_DN_2(UM,O1,N,NV,Q),L_UT_EMULT_DN(UM,UT,NV,IB,Q)
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP = L_UPAR_DN_1(UM,O1,N,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &              UPAR_DN_1(UM,O1,N)   * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
c            write(*,*)ib,n,nv,ut,um,q,'par2'
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  Other cases when N > NV
C  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0 ) THEN

C  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) +
     &            UPAR_DN_2(UM,O1,N)      * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

C  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP =  UPAR_DN_1(UM,O1,N) * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End layer clause

      ENDIF

C  Finish

      RETURN
      END

