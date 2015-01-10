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

C ###########################################################
C #                                                         #
C #   Contains the following Master subroutines             #
C #                                                         #
C #          UPUSER_INTENSITY (master)                      #
C #          DNUSER_INTENSITY (master)                      #
C #          MIFLUX_INTENSITY (master)                      #
C #                                                         #
C #   --For the discrete ordinate field                     #
C #                                                         #
C #          QUADINTENS_LEVEL_UP                            #
C #          QUADINTENS_LEVEL_DN                            #
C #          QUADINTENS_OFFGRID_UP                          #
C #          QUADINTENS_OFFGRID_DN                          #
C #                                                         #
C #   --For the post-processed field                        #
C #                                                         #
C #          GET_TOASOURCE                                  #
C #          GET_BOASOURCE                                  #
C #          WHOLELAYER_STERM_UP                            #
C #          WHOLELAYER_STERM_DN                            #
C #          PARTLAYER_STERM_UP                             #
C #          PARTLAYER_STERM_DN                             #
C #                                                         #
C ###########################################################

      SUBROUTINE UPUSER_INTENSITY
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_MVOUTPUT,
     I      DO_INCLUDE_DIRECTBEAM,
     I      FOURIER_COMPONENT,  IPARTIC,
     I      SURFACE_FACTOR,
     I      FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Fourier component, beam index

      INTEGER          FOURIER_COMPONENT
      INTEGER          IPARTIC

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION SURFACE_FACTOR

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER          UT, UTA, UM, IUM, I, IQD
      DOUBLE PRECISION LAYER_SOURCE       ( MAX_USER_STREAMS )
      DOUBLE PRECISION MSCAT_LAYERSOURCE  ( MAX_USER_STREAMS )
      DOUBLE PRECISION FINAL_SOURCE

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_OUT_USERTAUS
          DO UM = 1, N_USER_STREAMS
            IUM = USEROUTPUT_INDEX(UM)
            INTENSITY_F(UTA,IUM,IPARTIC,UPIDX) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  Initialize

      NUT = 0
      NC  = 0

C  Initialize post-processing recursion
C  ====================================

C  Get the BOA source terms (diffuse + direct)

      CALL GET_BOASOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      FOURIER_COMPONENT, IPARTIC, 
     I      SURFACE_FACTOR )

C  Set the cumulative source term equal to BOA values

      IF ( DO_USER_STREAMS ) THEN
        NC = 0
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) +
     &                          DIRECT_BOA_SOURCE(UM)
        ENDDO

C  Debug

c        IF ( DO_DEBUG_WRITE ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            write(98,*)UM,BOA_SOURCE(UM), DIRECT_BOA_SOURCE(UM)
c            if ( um.eq.1) then
c            IF ( FOURIER_COMPONENT.EQ.0.AND.DO_FDTEST)THEN
c       write(2,*) UM,BOA_SOURCE(UM),DIRECT_BOA_SOURCE(UM)
c            ELSE IF ( FOURIER_COMPONENT.EQ.0.AND..NOT.DO_FDTEST)THEN
c       write(3,*) UM,BOA_SOURCE(UM),DIRECT_BOA_SOURCE(UM)
c            ENDIF
c            endif
c          ENDDO
c        ENDIF

C  Save the BOA source terms if present
C  otherwise zero these values (very important - bug fixed 04/20/05)
c        IF ( SAVE_LAYER_MSST ) THEN
c          IF ( DO_INCLUDE_SURFACE ) THEN
c            DO UM = LOCAL_UM_START, N_USER_STREAMS
c              INTENSITY_MSST_BOA_F(UM,IPARTIC) =
c     &                  FLUX_MULTIPLIER * BOA_SOURCE(UM)
c              INTENSITY_DBST_BOA_F(UM,IPARTIC) =
c     &                  FLUX_MULTIPLIER * DIRECT_BOA_SOURCE(UM)
c            ENDDO
c          ELSE
c            DO UM = LOCAL_UM_START, N_USER_STREAMS
c              INTENSITY_MSST_BOA_F(UM,IPARTIC) = ZERO
c              INTENSITY_DBST_BOA_F(UM,IPARTIC) = ZERO
c            ENDDO
c          ENDIF
c        ENDIF

      ENDIF

C  Recursion Loop in Source function integration
C  =============================================

C  initialise cumulative source term loop

      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CALL WHOLELAYER_STERM_UP
     I        ( IPARTIC, N, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),
     O          DO_INCLUDE_THERMEMISS, LAYER_SOURCE, MSCAT_LAYERSOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) +
     &                   T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
            ENDDO
c            IF ( SAVE_LAYER_MSST ) THEN
c              DO UM = LOCAL_UM_START, N_USER_STREAMS
c                  INTENSITY_MSST_LAYER_F(N,UM,IPARTIC,UPIDX) =
c     &                      FLUX_MULTIPLIER * MSCAT_LAYERSOURCE(UM)
c              ENDDO
c            ENDIF
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature intensity calculation at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_OFFGRID_UP
     &   ( IPARTIC, UTA, UT, N, DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
               INTENSITY_F(UTA,IQD,IPARTIC,UPIDX) =
     &                      QUADINTENS(UTA,I,IPARTIC,UPIDX)
            ENDDO
          ENDIF

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_UP
     &       ( IPARTIC, UT, N, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),
     &         DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              FINAL_SOURCE = LAYER_SOURCE(UM) +
     &                       T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,NC)
              INTENSITY_F(UTA,IUM,IPARTIC,UPIDX) = 
     &                       FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output at layer boundaries
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_UP
     &           ( IPARTIC, UTA, NLEVEL, FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              INTENSITY_F(UTA,IQD,IPARTIC,UPIDX) =
     &                   QUADINTENS(UTA,I,IPARTIC,UPIDX)
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              INTENSITY_F(UTA,IUM,IPARTIC,UPIDX) =
     &                  FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
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

c          do i = 1, nstreams
c             write(79,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     &          (QUADINTENS(UM,I,IPARTIC,UPIDX),UM=1,N_OUT_USERTAUS)
c          enddo

c      um = 97
c      if ( do_fdtest) um = 98
c      if ( fourier_component.lt.4) then
c        write(um,*)fourier_component,ibeam,
c     &        INTENSITY_F(2,1,IBeam,UPIDX) 
c      endif
c      if (do_fdtest.and.fourier_component.eq.4.and.ibeam.eq.2)pause

C  Finish

      RETURN
      END

C

      SUBROUTINE DNUSER_INTENSITY
     I    ( DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_MVOUTPUT,
     I      FLUX_MULTIPLIER,
     I      FOURIER_COMPONENT,
     I      IPARTIC )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local inclusion flags

      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  beam index

      INTEGER          IPARTIC

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, IUM, NC, I, IQD
      DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION LAYER_SOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION FINAL_SOURCE

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_OUT_USERTAUS
          DO UM = 1, N_USER_STREAMS
            IUM = USEROUTPUT_INDEX(UM)
            INTENSITY_F(UTA,IUM,IPARTIC,DNIDX) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_TOASOURCE ( TOA_SOURCE )
        NC = 0
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          CUMSOURCE_DN(UM,NC) = TOA_SOURCE(UM)
        ENDDO
      ENDIF

C  Initialize

      NUT = 0
      NC  = 0

C  initialise cumulative source term loop

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
            CALL WHOLELAYER_STERM_DN
     &       ( IPARTIC, N, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),
     &         DO_INCLUDE_THERMEMISS, LAYER_SOURCE, MSCAT_LAYERSOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) +
     &                    T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
            ENDDO
c            IF ( SAVE_LAYER_MSST ) THEN
c              DO UM = LOCAL_UM_START, N_USER_STREAMS
c                INTENSITY_MSST_LAYER_F(N,UM,IPARTIC,DNIDX) =
c     &                    FLUX_MULTIPLIER * MSCAT_LAYERSOURCE(UM)
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
            CALL QUADINTENS_OFFGRID_DN
     &   ( IPARTIC, UTA, UT, N, DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              INTENSITY_F(UTA,IQD,IPARTIC,DNIDX) =
     &                    QUADINTENS(UTA,I,IPARTIC,DNIDX)
            ENDDO
          ENDIF

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_DN
     &       ( IPARTIC, UT, N, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),
     &         DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              FINAL_SOURCE = LAYER_SOURCE(UM) +
     &                       T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,NC)
              INTENSITY_F(UTA,IUM,IPARTIC,DNIDX) =
     &                       FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  Quadrature output at layer boundaries
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_DN
     &           ( IPARTIC, UTA, NLEVEL, FLUX_MULTIPLIER )
          ENDIF

C  Copy into storage if quad output required

          IF ( DO_QUAD_OUTPUT ) THEN
            DO I = 1, NSTREAMS
              IQD = QUADOUTPUT_INDEX(I)
              INTENSITY_F(UTA,IQD,IPARTIC,DNIDX) = 
     &                QUADINTENS(UTA,I,IPARTIC,DNIDX)
            ENDDO
          ENDIF

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              IUM = USEROUTPUT_INDEX(UM)
              INTENSITY_F(UTA,IUM,IPARTIC,DNIDX) =
     &                 FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
            ENDDO
          ENDIF

        ENDIF

C  debug

c        if ( uta .eq. n_out_usertaus ) then
c          do i = 1, nstreams
c             write(46,'(i4,f10.5,1p4e21.12)')i,xang(i),
c     &          (QUADINTENS(UM,I,IPARTIC,DNIDX),UM=1,N_OUT_USERTAUS)
c          enddo
c        endif

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

      SUBROUTINE QUADINTENS_LEVEL_UP
     &    ( IPARTIC, UTA, NLEVEL, FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          NLEVEL, UTA, IPARTIC

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          N, I, I1, AA, K
      DOUBLE PRECISION FM, THELP, SPAR, SHOM, HOM1, HOM2

C  For those optical depths at layer boundaries
C  --------------------------------------------

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the intensity at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling intensity
C  at the bottom of the atmosphere (treated separately).

      N = NLEVEL + 1
      FM = FLUX_MULTIPLIER

C  homogeneous and particular solution contributions SHOM and SPAR

C  For the lowest level

      IF ( NLEVEL .EQ. NLAYERS ) THEN

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = BOA_THTONLY_SOURCE(I)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,NLEVEL) * T_DELT_EIGEN(AA,NLEVEL)
              HOM2 = MCON_XVEC(I1,AA,NLEVEL)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WLOWER(I1,NLEVEL)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

C  For other levels in the atmosphere

      ELSE

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            THELP = BOA_THTONLY_SOURCE(I)
            DO K = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K)/X(I)
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,N)
              HOM2 = MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WUPPER(I1,N)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_LEVEL_DN
     &    ( IPARTIC, UTA, NLEVEL, FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          NLEVEL, UTA, IPARTIC

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          N, I, AA, K
      DOUBLE PRECISION FM, THELP, SPAR, SHOM, HOM1, HOM2

C  For those optical depths at layer boundaries
C  --------------------------------------------

      N = NLEVEL
      FM = FLUX_MULTIPLIER

C  Downwelling radiation at TOA ( or N = 0 ) is zero

      IF ( NLEVEL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = ZERO
        ENDDO

C  Other levels

      ELSE

C  Thermal transmittance solution, build from TOA downwards
C  Scattering solution, use the Discrete Ordinate solution

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = ZERO
            DO K = 1, N
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)/X(I)
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            SPAR = WLOWER(I,N)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
              HOM2 = MCON_XVEC(I,AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_OFFGRID_UP
     &  ( IPARTIC, UTA, UT, N, DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          N, UTA, UT, IPARTIC

C  Thermal emission flag

      LOGICAL          DO_INCLUDE_THERMEMISS

C  Flux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          I, I1, AA, K
      DOUBLE PRECISION THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I)
          DO K = NLAYERS, N+1, -1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K)/X(I)
          ENDDO
          THELP = THELP*T_DISORDS_UTUP(I,UT) + UT_T_PARTIC(I1,UT)/X(I)
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Homogeneous

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

C  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add the solar particular solutions

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = WUPPER(I1,N)*T_UTDN_MUBAR(UT,IPARTIC)
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ELSE
        CALL QUAD_GFUNCMULT ( IPARTIC, UT, N )
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = XPOS(I,AA,N)  * UT_GMULT_UP(AA,UT)
            PAR2 = XPOS(I1,AA,N) * UT_GMULT_DN(AA,UT)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_OFFGRID_DN
     &  ( IPARTIC, UTA, UT, N, DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  indices

      INTEGER          N, UTA, UT, IPARTIC

C  Thermal emission flag

      LOGICAL          DO_INCLUDE_THERMEMISS

C  FLux

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          I, I1, AA, K
      DOUBLE PRECISION THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N-1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K) / X(I)
          ENDDO
          THELP = THELP*T_DISORDS_UTDN(I,UT) + UT_T_PARTIC(I,UT) / X(I)
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Homogeneous

      DO I = 1, NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

C  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add the solar particular solutions

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO I = 1, NSTREAMS
          SPAR = WUPPER(I,N) * T_UTDN_MUBAR(UT,IPARTIC)
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ELSE
        CALL QUAD_GFUNCMULT ( IPARTIC, UT, N )
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = XPOS(I1,AA,N) * UT_GMULT_UP(AA,UT)
            PAR2 = XPOS(I,AA,N)  * UT_GMULT_DN(AA,UT)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = 
     &      QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE MIFLUX_INTENSITY 
     I       ( DO_INCLUDE_DIRECTBEAM, IPARTIC )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  input arguments

      INTEGER          IPARTIC
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  local variables

      INTEGER          I, IDIR, WDIR, UTA, UT, N
      DOUBLE PRECISION SUM_MI, SUM_FX, FMU0
      DOUBLE PRECISION DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

C  mean intensity and flux
C  -----------------------

C  direction loop

      DO IDIR = 1, N_DIRECTIONS

        WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

        DO UTA = 1, N_OUT_USERTAUS

          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            SUM_MI = SUM_MI + A(I)  * QUADINTENS(UTA,I,IPARTIC,WDIR)
            SUM_FX = SUM_FX + AX(I) * QUADINTENS(UTA,I,IPARTIC,WDIR)
          ENDDO
          MEAN_INTENSITY(UTA,IPARTIC,WDIR) = SUM_MI * HALF
          FLUX_INTEGRAL (UTA,IPARTIC,WDIR) = SUM_FX * PI2

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

              IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
                DIRECT_TRANS = INITIAL_TRANS(N,IPARTIC) *
     &                              T_UTDN_MUBAR(UT,IPARTIC)
                DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
                FMU0 = LOCAL_SZA(N,IPARTIC) * FLUX_FACTOR
                DIRECT_FLUX  = FMU0 * DIRECT_TRANS
                MEAN_INTENSITY(UTA,IPARTIC,WDIR) =
     *                  MEAN_INTENSITY(UTA,IPARTIC,WDIR) + DIRECT_MEANI
                FLUX_INTEGRAL(UTA,IPARTIC,WDIR) =
     &                  FLUX_INTEGRAL (UTA,IPARTIC,WDIR) + DIRECT_FLUX
              ENDIF

C  For the on-grid values

            ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
                IF ( N .EQ. 0 ) THEN
                  DIRECT_TRANS = ONE
                  FMU0 = LOCAL_SZA(1,IPARTIC) * FLUX_FACTOR
                ELSE
                  DIRECT_TRANS = 
     &                INITIAL_TRANS(N,IPARTIC)*T_DELT_MUBAR(N,IPARTIC)
                  FMU0 = LOCAL_SZA(N,IPARTIC) * FLUX_FACTOR
                ENDIF
                DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
                DIRECT_FLUX  = FMU0 * DIRECT_TRANS
                MEAN_INTENSITY(UTA,IPARTIC,WDIR) =
     *                  MEAN_INTENSITY(UTA,IPARTIC,WDIR) + DIRECT_MEANI
                FLUX_INTEGRAL(UTA,IPARTIC,WDIR) =
     &                  FLUX_INTEGRAL (UTA,IPARTIC,WDIR) + DIRECT_FLUX
              ENDIF
            ENDIF

C  End loop over optical depth output values

          ENDDO

C  Finish direct beam stuff

        ENDIF   

C  Continuation point for avoiding direct beam calculation

 455    CONTINUE

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE GET_TOASOURCE ( TOA_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Subroutine arguments

      DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS)

C  local variables

      INTEGER             UM

C  initialise TOA source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        TOA_SOURCE(UM) = ZERO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE GET_BOASOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      FOURIER_COMPONENT, IPARTIC,
     I      SURFACE_FACTOR )

C  Bottom of the atmosphere source term

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution, setup & reflectance variables (input, output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  surface multiplier and Fourier index

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT
      INTEGER          IPARTIC

C  local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          M, N, J, I, UM, AA, KL, K
      DOUBLE PRECISION PAR, HOM, REFLEC, KMULT, THELP(MAXSTREAMS)

C  Fourier and layer numbers

      M = FOURIER_COMPONENT
      N = NLAYERS

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  initialise boa source function

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          BOA_SOURCE(UM)        = ZERO
          DIRECT_BOA_SOURCE(UM) = ZERO
        ENDDO
      ENDIF

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

C  reflectance from surface
C  ------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN

C  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP(I) = ZERO
          ENDDO
          DO K = 1, NLAYERS
            DO I = 1, NSTREAMS
              THELP(I) = THELP(I)*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)
            ENDDO
          ENDDO
          DO I = 1, NSTREAMS
            IDOWNSURF(I) = A(I) * THELP(I)
          ENDDO
        ENDIF

C  Full solution: Downward intensity at computational angles (beam/homog)
C     --> Develop reflectance integrand  a(j).x(j).I(-j)

        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            PAR = WLOWER(I,N)
            HOM = ZERO
            DO AA = 1, NSTREAMS
              HOM = HOM + LCON_XVEC(I,AA,N)*T_DELT_EIGEN(AA,N) +
     &                    MCON_XVEC(I,AA,N)
            ENDDO
            IDOWNSURF(I) = AX(I) * ( PAR + HOM )
          ENDDO
        ENDIF

C  reflected multiple scatter intensity at user defined-angles
C  -----------------------------------------------------------

C  start loop over BRDF kernels

        DO KL = 1, N_BRDF_KERNELS

          KMULT = SURFACE_FACTOR * BRDF_FACTORS(KL)

C  ###### Lambertian reflectance (same for all user-streams)

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + IDOWNSURF(J)
              ENDDO
              REFLEC = KMULT * REFLEC
              IF ( DO_USER_STREAMS ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  ALBEDO_BOA_SOURCE(UM,KL) = REFLEC
                ENDDO
              ENDIF
              IF ( DO_QTHTONLY ) THEN
                DO I = 1, NSTREAMS
                  ALBEDO_BOA_THTONLY_SOURCE(I,KL) = REFLEC
                ENDDO
              ENDIF
            ELSE
              IF ( DO_USER_STREAMS ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  ALBEDO_BOA_SOURCE(UM,KL) = ZERO
                ENDDO
              ENDIF
            ENDIF

C  ###### bidirectional reflectance

          ELSE

            IF ( DO_USER_STREAMS ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                 REFLEC = REFLEC + IDOWNSURF(J) * USER_BIREFLEC(KL,UM,J)
                ENDDO
                ALBEDO_BOA_SOURCE(UM,KL) = KMULT * REFLEC
              ENDDO
            ENDIF

            IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + IDOWNSURF(J) * BIREFLEC(KL,I,J)
                ENDDO
                ALBEDO_BOA_THTONLY_SOURCE(I,KL) = KMULT * REFLEC
              ENDDO
            ENDIF
        
          ENDIF

C  End loop over albedo kernels

        ENDDO

C  Total BOA source term
C  ---------------------

        IF ( DO_USER_STREAMS ) THEN
          DO KL = 1, N_BRDF_KERNELS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              BOA_SOURCE(UM) = BOA_SOURCE(UM) +
     &           ALBEDO_BOA_SOURCE(UM,KL)
            ENDDO
          ENDDO
        ENDIF

C  Thermal tranmsittance only (additional quadrature terms if flagged)

        IF ( DO_QTHTONLY ) THEN
          DO KL = 1, N_BRDF_KERNELS
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) +
     &           ALBEDO_BOA_THTONLY_SOURCE(I,KL)
            ENDDO
          ENDDO
        ENDIF

C  Add direct beam if flagged

c        IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS
     &       .AND. .NOT. DO_DBCORRECTION ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DIRECT_BOA_SOURCE(UM) =
     &               DIRECT_BOA_SOURCE(UM)     +
     &               USER_DIRECT_BEAM(UM,IPARTIC)
          ENDDO
        ENDIF

C  End inclusion of surface terms

      ENDIF

C  debug

c      um = 97
c      if (   do_fdtest ) um = 98
c      if ( fourier_component.eq.0.and.ibeam.eq.1) then
c         write(um,'(1p4e17.9)')BOA_SOURCE(1)+DIRECT_BOA_SOURCE(1)
c      ENDIF

C  Add surface emission term if flagged

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_USER_STREAMS ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            BOA_SOURCE(UM) = BOA_SOURCE(UM) + SURFBB*USER_EMISSIVITY(UM)
          ENDDO
        ENDIF
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
            BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) +
     &           SURFBB * EMISSIVITY(I)
          ENDDO
        ENDIF
      ENDIF


C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_STERM_UP
     I    ( IPARTIC, N, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O      LAYERSOURCE, MSCAT_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  input

      INTEGER          N, IPARTIC
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION LAYERSOURCE       ( MAX_USER_STREAMS)
      DOUBLE PRECISION MSCAT_LAYERSOURCE ( MAX_USER_STREAMS)

C  local variables
C  ---------------

      INTEGER          AA, UM
      DOUBLE PRECISION SPAR, SHOM, SFOR1, SFOR2, TM

C  No layer source term if no scattering in the layer
C   Very important to zero both output terms (bug solved 04/20/05)

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
        MSCAT_LAYERSOURCE(UM) = ZERO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
        ENDDO
      ENDDO

C  Homogeneous first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) +
     &                  MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
        ENDDO
        LAYERSOURCE(UM) = SHOM
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution
C     Modulus 4.pi if solar sources are included (taken care of earlier)

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO UM = LOCAL_UM_START, N_USER_STREAMS
c          write(*,*)um,n,layersource(um),LAYER_TSUP_UP(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
        ENDDO
       ENDIF
      ENDIF

C  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular solar beam contributions
C     ( Classical vs. Green's function solutions )

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR2 =  EMULT_UP(UM,N,IPARTIC) * U_WPOS2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2 
        ENDDO
      ELSE
        CALL WHOLELAYER_GMULT_UP ( N, IPARTIC )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*SGMULT_UD(AA,UM,N) +
     &                    U_XNEG(UM,AA,N)*SGMULT_UU(AA,UM,N)
          ENDDO
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ENDIF

C  Options for adding the single scatter part of the beam solution
C  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
c        IF ( SAVE_LAYER_MSST ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
c          ENDDO
c        ENDIF
      ELSE
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N,IPARTIC)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_STERM_DN
     I    ( IPARTIC, N, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O      LAYERSOURCE, MSCAT_LAYERSOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  input

      INTEGER          N, IPARTIC
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)
      DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)

C  local variables
C  ---------------

      INTEGER          AA, UM
      DOUBLE PRECISION SPAR, SHOM, SFOR1, SFOR2, TM

C  No layer source term if no scattering in the layer

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
        MSCAT_LAYERSOURCE(UM) = ZERO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
        ENDDO
      ENDDO

C  Homogeneous first

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) +
     &                  MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
        ENDDO
        LAYERSOURCE(UM) = SHOM
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     ----- only with Green's function solution

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
        ENDDO
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular solar beam contributions
C     ( Classical vs. Green's function solutions )

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR2 =  EMULT_DN(UM,N,IPARTIC) * U_WNEG2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2
        ENDDO
      ELSE
        CALL WHOLELAYER_GMULT_DN ( N, IPARTIC )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*SGMULT_DD(AA,UM,N) +
     &                    U_XPOS(UM,AA,N)*SGMULT_DU(AA,UM,N)
          ENDDO
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ENDIF

C  Options
C  .. If operating in Ms-mode only, copy multiple scatter term
C  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
c        IF ( SAVE_LAYER_MSST ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c            MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
c          ENDDO
c        ENDIF
      ELSE
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N,IPARTIC)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_STERM_UP
     I      ( IPARTIC, UT, N, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O        LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  input

      INTEGER          N, UT, IPARTIC
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)

C  local variables
C  ---------------

      INTEGER          AA, UM
      DOUBLE PRECISION SPAR, SHOM, SFOR, TM

C  No layer source term if no scattering in the layer

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
        ENDDO
      ENDDO

c  homogeneous solutions

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UD(AA,UM,UT) +
     &                  MCON_UXVEC(UM,AA)*UT_HMULT_UU(AA,UM,UT)
        ENDDO
        LAYERSOURCE(UM) = SHOM
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal term (direct and diffuse)
C     ----- only with Green's function solution

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         LAYERSOURCE(UM) = LAYERSOURCE(UM)+LAYER_TSUP_UTUP(UM,UT)*TM
        ENDDO
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular solar beam contributions
C     ( Classical vs. Green's function solutions )

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR =  UT_EMULT_UP(UM,UT,IPARTIC) * U_WPOS2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ELSE
        CALL PARTLAYER_GMULT_UP ( N, UT, IPARTIC )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*UT_SGMULT_UD(AA,UM,UT) +
     &                    U_XNEG(UM,AA,N)*UT_SGMULT_UU(AA,UM,UT)
          ENDDO
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR = U_WPOS1(UM,N) * UT_EMULT_UP(UM,UT,IPARTIC)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
        ENDDO
      ENDIF
          
C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_STERM_DN
     I      ( IPARTIC, UT, N, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O        LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  input

      INTEGER          N, UT, IPARTIC
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)

C  local variables
C  ---------------

      INTEGER          AA, UM
      DOUBLE PRECISION SPAR, SHOM, SFOR, TM

C  No layer source term if no scattering in the layer

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Save some calculation time

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
        ENDDO
      ENDDO

c  homogeneous solutions

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DD(AA,UM,UT) +
     &                  MCON_UXVEC(UM,AA)*UT_HMULT_DU(AA,UM,UT)
        ENDDO
        LAYERSOURCE(UM) = SHOM
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal term (direct and diffuse)
C     ----- only with Green's function solution

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         LAYERSOURCE(UM) = LAYERSOURCE(UM)+LAYER_TSUP_UTDN(UM,UT)*TM
        ENDDO
       ENDIF
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular solar beam contributions
C     ( Classical vs. Green's function solutions )

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR = UT_EMULT_DN(UM,UT,IPARTIC) * U_WNEG2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ELSE
        CALL PARTLAYER_GMULT_DN ( N, UT, IPARTIC )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*UT_SGMULT_DD(AA,UM,UT) +
     &                    U_XPOS(UM,AA,N)*UT_SGMULT_DU(AA,UM,UT)
          ENDDO
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
        ENDDO
      ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SFOR = U_WNEG1(UM,N) * UT_EMULT_DN(UM,UT,IPARTIC)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
        ENDDO
      ENDIF

C  Finish

      RETURN
      END
