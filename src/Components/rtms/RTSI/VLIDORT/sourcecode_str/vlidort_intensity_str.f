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
C #              VLIDORT_UPUSER_INTENSITY                       #
C #              VLIDORT_DNUSER_INTENSITY                       #
C #                                                             #
C #              GET_TOASOURCE                                  #
C #              BOA_LAMBERTIAN_SOURCE                          #
C #                                                             #
C #              WHOLELAYER_STERM_UP                            #
C #              WHOLELAYER_STERM_DN                            #
C #              PARTLAYER_STERM_UP                             #
C #              PARTLAYER_STERM_DN                             #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_UPUSER_INTENSITY
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_THERMEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      FOURIER_COMPONENT,  IBEAM,
     I      SURFACE_FACTOR,
     I      GET_BOA_SOURCE,
     I      FLUX_MULTIPLIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Fourier component, beam index

      INTEGER          FOURIER_COMPONENT
      INTEGER          IBEAM

C  External function

      EXTERNAL         GET_BOA_SOURCE

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION SURFACE_FACTOR

C  Local include flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  local variables
C  ---------------

      LOGICAL          SFLAG
      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, NC, O1
      INTEGER          UT, UTA, UM
      DOUBLE PRECISION LAYER_SOURCE(MAX_USER_STREAMS,MAXSTOKES)
      DOUBLE PRECISION FINAL_SOURCE

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

      IF ( DO_USER_STREAMS ) THEN

C  Get the BOA source terms (diffuse + direct)
C    External Call, introduced 27 March 2007.

        CALL GET_BOA_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      FOURIER_COMPONENT, IBEAM, 
     I      SURFACE_FACTOR )

C  Set the cumulative source term equal to BOA values

        NC = 0
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            CUMSOURCE_UP(UM,O1,NC) = BOA_DIFFUSE_SOURCE(UM,O1) +
     &                               BOA_DIRECT_SOURCE(UM,O1)
!          write(*,*)um,BOA_DIFFUSE_SOURCE(UM,O1),
!     &                               BOA_DIRECT_SOURCE(UM,O1)
          ENDDO
        ENDDO
!        pause

C  Debug write

c        IF ( DO_DEBUG_WRITE ) THEN
c          DO UM = LOCAL_UM_START, N_USER_STREAMS
c           DO O1 = 1, NSTOKES
c            if ( um.eq.1.and.o1.eq.1) then
c              IF ( FOURIER_COMPONENT.EQ.0.AND.DO_FDTEST)THEN
c                write(81,*) UM,O1,BOA_DIFFUSE_SOURCE(UM,O1),
C     &                            BOA_DIRECT_SOURCE(UM,O1)
c              ELSE IF ( FOURIER_COMPONENT.EQ.0.AND..NOT.DO_FDTEST)THEN
c                write(82,*) UM,O1,BOA_DIFFUSE_SOURCE(UM,O1),
C     &                            BOA_DIRECT_SOURCE(UM,O1)
c              ENDIF
c            endif
c           ENDDO
c          ENDDO
c        ENDIF

      ENDIF

C  Recursion Loop in Source function integration
C  =============================================

C  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = NLAYERS + 1 - N
            CALL WHOLELAYER_STERM_UP
     &        ( N, FOURIER_COMPONENT, IBEAM, 
     &          SFLAG, DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                CUMSOURCE_UP(UM,O1,NC) = LAYER_SOURCE(UM,O1) +
     &                   T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_UP
     &         ( N, UT, IBEAM, SFLAG,
     &           DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = LAYER_SOURCE(UM,O1) +
     &                       T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,O1,NC)
                STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = 
     &                       FLUX_MULTIPLIER * FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term
C   Commented out code helps with FD testing

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,O1,NC)
                STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = FINAL_SOURCE
!                if (DABS(FINAL_SOURCE).GT.1.0d-10 ) then           ! FDHELP
!                  STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = FINAL_SOURCE   ! FDHELP
!                endif                                              ! FDHELP
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  debug write

c        IF ( DO_DEBUG_WRITE ) THEN
c         DO UM = LOCAL_UM_START, N_USER_STREAMS
c          DO O1 = 1, NSTOKES
c            if (do_fdtest) then
c                write(69,'(5i3,1pe20.10)')
c     &         fourier_component,uta,ibeam,um,o1,       
c     &           STOKES_F(UTA,UM,IBEAM,O1,UPIDX)
c            endif
c          ENDDO
c         ENDDO
c        ENDIF
         
C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

c      if ( fourier_component.eq.6) pause 'user intensity'

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_DNUSER_INTENSITY
     I    ( DO_INCLUDE_THERMEMISS,
     I      FLUX_MULTIPLIER,
     I      FOURIER_COMPONENT, IBEAM )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Thermal emissioin control

      LOGICAL          DO_INCLUDE_THERMEMISS

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  beam index

      INTEGER          IBEAM

C  local variables
C  ---------------

      LOGICAL          SFLAG
      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, NC, O1
      DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS,MAXSTOKES)
      DOUBLE PRECISION LAYER_SOURCE(MAX_USER_STREAMS,MAXSTOKES)
      DOUBLE PRECISION FINAL_SOURCE

C  Zero all Fourier components close to zenith

      IF ( DO_USER_STREAMS ) THEN
        IF ( LOCAL_UM_START .GT. 1 ) THEN
          DO UTA = 1, N_USER_LEVELS
            DO UM = 1, LOCAL_UM_START - 1
              DO O1 = 1, NSTOKES
                STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_TOASOURCE ( TOA_SOURCE )
        NC = 0
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            CUMSOURCE_DN(UM,O1,NC) = TOA_SOURCE(UM,O1)
          ENDDO
        ENDDO
      ENDIF

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
            CALL WHOLELAYER_STERM_DN
     &        ( N, IBEAM, SFLAG, DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                CUMSOURCE_DN(UM,O1,NC) = LAYER_SOURCE(UM,O1) +
     &                    T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_DN
     &    ( N, UT, IBEAM, SFLAG, DO_INCLUDE_THERMEMISS, LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = LAYER_SOURCE(UM,O1) +
     &                       T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,O1,NC)
                STOKES_F(UTA,UM,IBEAM,O1,DNIDX) =
     &                       FLUX_MULTIPLIER * FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                STOKES_F(UTA,UM,IBEAM,O1,DNIDX) =
     &                 FLUX_MULTIPLIER * CUMSOURCE_DN(UM,O1,NC)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

c        IF ( DO_DEBUG_WRITE ) THEN
c         DO UM = LOCAL_UM_START, N_USER_STREAMS
c          DO O1 = 1, NSTOKES
c            if (do_fdtest) then
c                if (uta.gt.1)write(69,'(5i3,1pe20.10)')
c     &         fourier_component,uta-1,ibeam,um,o1,       
c     &           STOKES_F(UTA,UM,IBEAM,O1,DNIDX)
c            endif
c          ENDDO
c         ENDDO
c        ENDIF

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

        SUBROUTINE GET_TOASOURCE ( TOA_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Subroutine arguments

      DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS,MAXSTOKES)

C  local variables

      INTEGER          UM, O1

C  initialise TOA source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          TOA_SOURCE(UM,O1) = ZERO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE BOA_LAMBERTIAN_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_MVOUTPUT,
     I      FOURIER_COMPONENT, IBEAM,
     I      SURFACE_FACTOR )

C  Bottom of the atmosphere source term, Lambertian

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution, setup & reflectance variables (input, output)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

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
      INTEGER          IBEAM

C  local variables
C  ---------------

      INTEGER          N, J, I, IR, IROW, UM, O1, O11
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION REFLEC, STOTAL, KMULT, LOCAL_EMISS, FP_SBB
      DOUBLE PRECISION SPAR, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR

      LOGICAL          DO_QTHTONLY
      DOUBLE PRECISION THELP(MAXSTREAMS)

C  initialise boa source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          BOA_DIFFUSE_SOURCE(UM,O1) = ZERO
          BOA_DIRECT_SOURCE(UM,O1)  = ZERO
        ENDDO
      ENDDO

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  Layer index and offset

      N = NLAYERS
      KO1 = K_REAL(N) + 1

C  1-1 element index

      O11 = 1

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
            STOKES_DOWNSURF(I,O11) = QUAD_WEIGHTS(I) * THELP(I)
          ENDDO
        ENDIF

C  Full solution: Downward intensity at computational angles (beam/homog)
C     --> Develop reflectance integrand  a(j).x(j).I(-j)

        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

C  Downward Stokes vector at surface at computational angles
C    reflectance integrand  a(j).x(j).I(-j)

         DO I = 1, NSTREAMS
           IR = ( I - 1 ) * NSTOKES
           DO O1 = 1, NSTOKES
             IROW = IR + O1
             SPAR = WLOWER(I,O1,N)
             SHOM_R = ZERO
             DO K = 1, K_REAL(N)
               LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
               MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
               HOM1 = LXR * T_DELT_EIGEN(K,N)
               HOM2 = MXR
               SHOM_R = SHOM_R + HOM1 + HOM2
             ENDDO
             SHOM_CR = ZERO
             DO K = 1, K_COMPLEX(N)
               K0 = 2*K-2
               K1 = KO1 + K0
               K2 = K1  + 1
               LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) -
     &                   LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
               LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) +
     &                   LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
               MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) -
     &                   MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
               HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N)
     &                 -LXR_CI*T_DELT_EIGEN(K2,N)
               HOM2CR = MXR_CR
               SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
             ENDDO
             STOTAL = SPAR + SHOM_R + SHOM_CR
             STOKES_DOWNSURF(I,O1) =  QUAD_STRMWTS(I) * STOTAL
           ENDDO
         ENDDO
        ENDIF

C  reflected multiple scatter intensity at user defined-angles

        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
          O1 = 1
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
          ENDDO
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = REFLEC
            ENDDO
          ENDIF
          REFLEC = KMULT * REFLEC
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            BOA_DIFFUSE_SOURCE(UM,O1) = REFLEC
          ENDDO
        ENDIF

C  Add direct beam if flagged
C   Addition of the DBCORRECTION clause is now required

        IF ( DO_INCLUDE_DIRECTBEAM.AND..NOT.DO_DBCORRECTION ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              BOA_DIRECT_SOURCE(UM,O1) = USER_DIRECT_BEAM(UM,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF

C  End inclusion of surface terms

      ENDIF

C  Add surface emission term if flagged

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        FP_SBB = SURFBB
        IF ( DO_SOLAR_SOURCES ) FP_SBB = PI4 * SURFBB
        if ( do_lambertian_surface ) then
          o1 = 1
          local_emiss = FP_SBB * ( one - lambertian_albedo )
          DO UM = LOCAL_UM_START, N_USER_STREAMS
              BOA_DIFFUSE_SOURCE(UM,O1) =
     &      BOA_DIFFUSE_SOURCE(UM,O1) + LOCAL_EMISS
          ENDDO
        else
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              BOA_DIFFUSE_SOURCE(UM,O1) =
     &        BOA_DIFFUSE_SOURCE(UM,O1) + FP_SBB*USER_EMISSIVITY(UM,O1)
            ENDDO
          ENDDO
        endif
        IF ( DO_QTHTONLY ) THEN
          if ( do_lambertian_surface ) then
            o1 = 1
            local_emiss = FP_SBB * ( one - lambertian_albedo )
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I)+LOCAL_EMISS
            ENDDO
          else
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) +
     &                 FP_SBB * EMISSIVITY(I,O1)
              ENDDO
            ENDDO
          endif
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_STERM_UP
     I   ( N, M, IBEAM, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O     LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  input control

      INTEGER          N, M, IBEAM
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Output

      DOUBLE PRECISION LAYERSOURCE      (MAX_USER_STREAMS,MAXSTOKES)

C  local variables

      LOGICAL          L1
      INTEGER          UM, O1, K, KO1, K0, K1, K2, M1
      DOUBLE PRECISION H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

C  Fourier number (debug only)

      M1 = M
      L1 = SOURCETERM_FLAG
      
C   Very important to zero output term (bug solved 04/20/05)

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          LAYERSOURCE(UM,O1)       = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Offset

      KO1 = K_REAL(N) + 1

C  Whole layer source function Homogeneous only

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
            MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
            H_R = H_R + LUXR * HMULT_2(K,UM,N) +
     &                  MUXR * HMULT_1(K,UM,N)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) -
     &               LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) +
     &               LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) -
     &               MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) +
     &               MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) -
     &                    LUX_CI * HMULT_2(K2,UM,N) +
     &                    MUX_CR * HMULT_1(K1,UM,N) -
     &                    MUX_CI * HMULT_1(K2,UM,N)
          ENDDO
          LAYERSOURCE(UM,O1) = H_R + H_CR
        ENDDO
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)
C     Modulus 4.pi if solar sources are included (taken care of earlier)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        O1 = 1
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + 
     &              LAYER_TSUP_UP(UM,N)*TM
        ENDDO
      ENDIF

C  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add particular integral contribution (solar term)

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR2 =  EMULT_UP(UM,N,IBEAM) * UPAR_UP_2(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
          ENDDO
        ENDDO
      ELSE
C  PLACEHOLDER FOR GREENS FUNCTION
      ENDIF

C  If operating in Ms-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

C  Full radiance mode, add single scatter part

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          SFOR1 = UPAR_UP_1(UM,O1,N) * EMULT_UP(UM,N,IBEAM)
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
        ENDDO
      ENDDO

C  debug. NSTOKES = 1 checks out.

c      DO UM = LOCAL_UM_START, N_USER_STREAMS
c        write(*,*)N,um, LAYERSOURCE(UM,1), EMULT_UP(UM,N,IBEAM)
c      ENDDO
c      if (n.eq.1)pause

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_STERM_DN
     I    ( N, IBEAM, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O      LAYERSOURCE ) 

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  Control input

      INTEGER          N, IBEAM
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS,MAXSTOKES)

C  local variables

      LOGICAL          L1
      INTEGER          UM, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

C  Zeroing is very important here

      L1 = SOURCETERM_FLAG
      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          LAYERSOURCE(UM,O1)       = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Offest

      KO1 = K_REAL(N) + 1

C  Whole layer source function, homoegeneous solution contribution

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
            MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
            H_R = H_R + LUXR * HMULT_1(K,UM,N) +
     &                  MUXR * HMULT_2(K,UM,N)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) -
     &               LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) +
     &               LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) -
     &               MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) +
     &               MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H_CR = H_CR + LUX_CR * HMULT_1(K1,UM,N) -
     &                    LUX_CI * HMULT_1(K2,UM,N) +
     &                    MUX_CR * HMULT_2(K1,UM,N) -
     &                    MUX_CI * HMULT_2(K2,UM,N)
          ENDDO
          LAYERSOURCE(UM,O1) = H_R + H_CR
        ENDDO
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal emission term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        O1 = 1
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1)
     &                           + LAYER_TSUP_DN(UM,N)*TM
        ENDDO
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR2 =  EMULT_DN(UM,N,IBEAM) * UPAR_DN_2(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
          ENDDO
        ENDDO
      ELSE
C  PLACEHOLDER GREENS FUNCTION
      ENDIF

C  If operating in Ms-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

C  Full radiance mode, add single scatter part

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          SFOR1 = UPAR_DN_1(UM,O1,N) * EMULT_DN(UM,N,IBEAM)
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_STERM_UP
     I  ( N, UT, IBEAM, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O    LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  Control input

      INTEGER          N, IBEAM
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS,MAXSTOKES)

C  local variables

      LOGICAL          L1
      INTEGER          UM, UT, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

C  Zeroing is very important here

      SFOR2 = ZERO
      L1 = SOURCETERM_FLAG
      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          LAYERSOURCE(UM,O1)       = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Offset

      KO1 = K_REAL(N) + 1

C  Whole layer source function, homogeneous contribution

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
            MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
            H_R = H_R + LUXR * UT_HMULT_UD(K,UM,UT) +
     &                  MUXR * UT_HMULT_UU(K,UM,UT)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) -
     &               LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) +
     &               LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) -
     &               MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) +
     &               MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H_CR = H_CR + LUX_CR * UT_HMULT_UD(K1,UM,UT) -
     &                    LUX_CI * UT_HMULT_UD(K2,UM,UT) +
     &                    MUX_CR * UT_HMULT_UU(K1,UM,UT) -
     &                    MUX_CI * UT_HMULT_UU(K2,UM,UT)
          ENDDO
          LAYERSOURCE(UM,O1) = H_R + H_CR
        ENDDO
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        O1 = 1
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + 
     &                   LAYER_TSUP_UTUP(UM,UT)*TM
        ENDDO
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR2 =  UT_EMULT_UP(UM,UT,IBEAM) * UPAR_UP_2(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
          ENDDO
        ENDDO
      ELSE
C  PLACEHOLDER GREENS FUNCTION
      ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR1 = UT_EMULT_UP(UM,UT,IBEAM) * UPAR_UP_1(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
          ENDDO
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_STERM_DN
     I  ( N, UT, IBEAM, SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS,
     O    LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main solution and multiplier variables (input/output)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  arguments
C  ---------

C  Control input

      INTEGER          N, IBEAM
      LOGICAL          SOURCETERM_FLAG
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Output

      DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS,MAXSTOKES)

C  local variables

      LOGICAL          L1
      INTEGER          UM, UT, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

C  Zeroing is very important here

      SFOR2 = ZERO
      L1 = SOURCETERM_FLAG
      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          LAYERSOURCE(UM,O1)       = ZERO
        ENDDO
      ENDDO

C  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

C  Offset

      KO1 = K_REAL(N) + 1

C  Whole layer source function, homogeneous contribution

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
            MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
            H_R = H_R + LUXR * UT_HMULT_DD(K,UM,UT) +
     &                  MUXR * UT_HMULT_DU(K,UM,UT)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) -
     &               LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) +
     &               LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) -
     &               MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) +
     &               MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H_CR = H_CR + LUX_CR * UT_HMULT_DD(K1,UM,UT) -
     &                    LUX_CI * UT_HMULT_DD(K2,UM,UT) +
     &                    MUX_CR * UT_HMULT_DU(K1,UM,UT) -
     &                    MUX_CI * UT_HMULT_DU(K2,UM,UT)
          ENDDO
          LAYERSOURCE(UM,O1) = SFOR2 + H_R + H_CR
        ENDDO
      ENDDO

C  Continuation point

 6789 continue

C  Add thermal term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        O1 = 1
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + 
     &                   LAYER_TSUP_UTDN(UM,UT)*TM
        ENDDO
      ENDIF

C  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR2 =  UT_EMULT_DN(UM,UT,IBEAM) * UPAR_DN_2(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
       ENDDO
        ENDDO
      ELSE
C  PLACEHOLDER GREENS FUNCTION
      ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR1 = UT_EMULT_DN(UM,UT,IBEAM) * UPAR_DN_1(UM,O1,N)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
          ENDDO
        ENDDO
      ENDIF

C  Finish

      RETURN
      END
