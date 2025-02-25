! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
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
! #              VLIDORT_UPUSER_INTENSITY (master)              #
! #              VLIDORT_DNUSER_INTENSITY (master)              #
! #                                                             #
! #       ...Calling....                                        #
! #                                                             #
! #              GET_TOA_SOURCE                                 #
! #              GET_BOA_SOURCE                                 #
! #                                                             #
! #              WHOLELAYER_STERM_UP                            #
! #              WHOLELAYER_STERM_DN                            #
! #              PARTLAYER_STERM_UP                             #
! #              PARTLAYER_STERM_DN                             #
! #                                                             #
! #              VLIDORT_INTEGRATED_OUTPUT (master)             #
! #                                                             #
! #       ...Calling....                                        #
! #                                                             #
! #              QUADINTENS_LEVEL_UP                            #
! #              QUADINTENS_LEVEL_DN                            #
! #              QUADINTENS_OFFGRID_UP                          #
! #              QUADINTENS_OFFGRID_DN                          #
! #                                                             #
! #              VLIDORT_CONVERGE  (master)                     #
! #              VLIDORT_CONVERGE_OBSGEO  (master)              #
! #                                                             #
! ###############################################################

      MODULE vlidort_intensity

      PRIVATE
      PUBLIC :: VLIDORT_UPUSER_INTENSITY,&
                VLIDORT_DNUSER_INTENSITY,&
                VLIDORT_INTEGRATED_OUTPUT,&
                VLIDORT_CONVERGE,&
                VLIDORT_CONVERGE_OBSGEO

      CONTAINS

      SUBROUTINE VLIDORT_UPUSER_INTENSITY ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, &
        DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT, IBEAM, &
        FLUX_MULTIPLIER, SURFACE_FACTOR, DO_DEBUG_WRITE, DO_FDTEST, &
        NSTOKES, NLAYERS, N_USER_LEVELS, DO_TOA_CONTRIBS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, &
        UTAU_LEVEL_MASK_UP, T_DELT_USERM, T_UTUP_USERM, &
        CUMTRANS, &
        DO_QUAD_OUTPUT, NSTREAMS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, USER_BRDF_F, EMISSIVITY, USER_EMISSIVITY, &
        SURFBB, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, QUAD_WEIGHTS, &
        QUAD_STRMWTS, MUELLER_INDEX, K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, WLOWER, LCON, MCON, &
        USER_DIRECT_BEAM, &
        T_DELT_DISORDS, T_DELT_EIGEN, T_WLOWER, &
        DO_SOLAR_SOURCES, DO_CLASSICAL_SOLUTION, &
        DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, &
        UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
        HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP, &
        UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, LAYER_TSUP_UTUP, BRDF_F, &
        STOKES_DOWNSURF, &
        CUMSOURCE_UP, STOKES_F, MS_CONTRIBS_F, BOA_THTONLY_SOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           IBEAM
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      LOGICAL, INTENT (IN) ::           DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::           DO_FDTEST
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::           NSTREAMS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           DO_MSMODE_THERMAL
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DOWNSURF &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: CUMSOURCE_UP &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NC, O1
      INTEGER ::          UT, UTA, UM, LUM
      DOUBLE PRECISION :: LAYER_SOURCE (MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: FINAL_SOURCE

      DOUBLE PRECISION :: BOA_DIFFUSE_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: BOA_DIRECT_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  Local user index

     LUM = 1

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Get the BOA source terms (diffuse + direct)

        CALL GET_BOA_SOURCE ( DO_MSMODE_THERMAL, &
          DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, &
          DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT, &
          DO_OBSERVATION_GEOMETRY, &
          FOURIER_COMPONENT, IBEAM, &
          SURFACE_FACTOR, &
          DO_QUAD_OUTPUT, NSTOKES, &
          NSTREAMS, NLAYERS, &
          DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
          USER_BRDF_F, &
          EMISSIVITY, USER_EMISSIVITY, &
          SURFBB, DO_THERMAL_TRANSONLY, &
          DO_DBCORRECTION, QUAD_WEIGHTS, &
          QUAD_STRMWTS, N_USER_STREAMS, &
          MUELLER_INDEX, LOCAL_UM_START, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          WLOWER, LCON, MCON, &
          USER_DIRECT_BEAM, &
          T_DELT_DISORDS, T_DELT_EIGEN, &
          T_WLOWER, BRDF_F, &
          STOKES_DOWNSURF, BOA_DIFFUSE_SOURCE, &
          BOA_DIRECT_SOURCE, BOA_THTONLY_SOURCE )

!  Set the cumulative source term equal to BOA values

        NC = 0
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              CUMSOURCE_UP(UM,O1,NC) = BOA_DIFFUSE_SOURCE(UM,O1) + &
                                       BOA_DIRECT_SOURCE(UM,O1)
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            CUMSOURCE_UP(IBEAM,O1,NC) = BOA_DIFFUSE_SOURCE(IBEAM,O1) + &
                                        BOA_DIRECT_SOURCE(IBEAM,O1)
          ENDDO
        ENDIF

!  Debug write

        IF ( DO_DEBUG_WRITE ) THEN
          write(*,*)
          write(*,*) 'at int up: boa dbg'
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO O1 = 1, NSTOKES
              if ( um.eq.ibeam .and. o1.eq.1 ) then
                IF ( FOURIER_COMPONENT.EQ.0.AND.DO_FDTEST)THEN
                  write(*,*) UM,O1,BOA_DIFFUSE_SOURCE(UM,O1),&
                                   BOA_DIRECT_SOURCE(UM,O1)
                ELSE IF ( FOURIER_COMPONENT.EQ.0.AND..NOT.DO_FDTEST)THEN
                  write(*,*) UM,O1,BOA_DIFFUSE_SOURCE(UM,O1),&
                                   BOA_DIRECT_SOURCE(UM,O1)
                ENDIF
              endif
             ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
             if ( o1.eq.1 ) then
               IF ( FOURIER_COMPONENT.EQ.0.AND.DO_FDTEST)THEN
                 write(*,*) IBEAM,O1,BOA_DIFFUSE_SOURCE(IBEAM,O1),&
                                     BOA_DIRECT_SOURCE(IBEAM,O1)
               ELSE IF ( FOURIER_COMPONENT.EQ.0.AND..NOT.DO_FDTEST)THEN
                 write(*,*) IBEAM,O1,BOA_DIFFUSE_SOURCE(IBEAM,O1),&
                                     BOA_DIRECT_SOURCE(IBEAM,O1)
               ENDIF
             endif
            ENDDO
          ENDIF
        ENDIF

      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = NLAYERS + 1 - N

            CALL WHOLELAYER_STERM_UP ( &
              N, FOURIER_COMPONENT, &
              IBEAM, SFLAG, &
              DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
              DO_OBSERVATION_GEOMETRY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_UPDN, UHOM_UPUP, &
              UPAR_UP_1, UPAR_UP_2, &
              HMULT_1, HMULT_2, &
              EMULT_UP, &
              LAYER_TSUP_UP, &
              LAYER_SOURCE )

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  IF ( DO_TOA_CONTRIBS ) THEN
                    MS_CONTRIBS_F(UM,IBEAM,O1,N) = FLUX_MULTIPLIER * &
                        CUMTRANS(N,UM) * LAYER_SOURCE(UM,O1)
                  ENDIF
                  CUMSOURCE_UP(UM,O1,NC) = LAYER_SOURCE(UM,O1) + &
                      T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,O1,NC-1)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                IF ( DO_TOA_CONTRIBS ) THEN
                  MS_CONTRIBS_F(LUM,IBEAM,O1,N) = FLUX_MULTIPLIER * &
                      CUMTRANS(N,IBEAM) * LAYER_SOURCE(IBEAM,O1)
                ENDIF
                CUMSOURCE_UP(IBEAM,O1,NC) = LAYER_SOURCE(IBEAM,O1) + &
                    T_DELT_USERM(N,IBEAM)*CUMSOURCE_UP(IBEAM,O1,NC-1)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_UP ( &
              N, UT, IBEAM, &
              SFLAG, DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
              DO_OBSERVATION_GEOMETRY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_UPDN, UHOM_UPUP, &
              UPAR_UP_1, UPAR_UP_2, &
              UT_HMULT_UU, UT_HMULT_UD, &
              UT_EMULT_UP, &
              LAYER_TSUP_UTUP, &
              LAYER_SOURCE )

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  FINAL_SOURCE = LAYER_SOURCE(UM,O1) + &
                               T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,O1,NC)
                  STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = &
                               FLUX_MULTIPLIER * FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = LAYER_SOURCE(IBEAM,O1) + &
                             T_UTUP_USERM(UT,IBEAM)*CUMSOURCE_UP(IBEAM,O1,NC)
                STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = &
                             FLUX_MULTIPLIER * FINAL_SOURCE
              ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!   Commented out code helps with FD testing

          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  FINAL_SOURCE = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,O1,NC)
                  STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = FINAL_SOURCE

!                  if (DABS(FINAL_SOURCE).GT.1.0d-10 ) then           ! FD
!                    STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = FINAL_SOURCE   ! FD
!                  endif                                              ! FD
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = FLUX_MULTIPLIER * CUMSOURCE_UP(IBEAM,O1,NC)
                STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = FINAL_SOURCE

!                if (DABS(FINAL_SOURCE).GT.1.0d-10 ) then           ! FD
!                  STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = FINAL_SOURCE  ! FD
!                endif                                              ! FD
              ENDDO
            ENDIF
          ENDIF

        ENDIF

!  debug write

        IF ( DO_DEBUG_WRITE ) THEN
          write(*,*)
          write(*,*) 'at int up: stokes_f dbg'
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                !if (do_fdtest) then
                 !write(*,*) 'LOCAL_UM_START = ',LOCAL_UM_START
                 !stop
                 if (um.eq.ibeam .and. o1.eq.1) then
                  write(*,'(4i3,1pe20.10)') &
                    fourier_component,uta,ibeam,o1, &
                    STOKES_F(UTA,UM,IBEAM,O1,UPIDX)
                 endif
                !endif
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              !if (do_fdtest) then
              if (o1.eq.1) then
                write(*,'(4i3,1pe20.10)') &
                  fourier_component,uta,ibeam,o1, &
                  STOKES_F(UTA,LUM,IBEAM,O1,UPIDX)
              endif
              !endif
            ENDDO
          ENDIF
          !stop
        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UPUSER_INTENSITY

!

      SUBROUTINE VLIDORT_DNUSER_INTENSITY ( &
        DO_INCLUDE_THERMEMISS, FOURIER_COMPONENT, &
        IBEAM, FLUX_MULTIPLIER, &
        DO_DEBUG_WRITE, NSTOKES, &
        N_USER_LEVELS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
        LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTDN_USERM, &
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        HMULT_1, HMULT_2, &
        EMULT_DN, LAYER_TSUP_DN, &
        UT_HMULT_DU, UT_HMULT_DD, &
        UT_EMULT_DN, LAYER_TSUP_UTDN, &
        CUMSOURCE_DN, STOKES_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::          IBEAM
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: CUMSOURCE_DN &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::          UT, UTA, UM, NC, O1, LUM
      DOUBLE PRECISION :: TOA_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: FINAL_SOURCE

!  Local user index

      LUM = 1

!  Zero all Fourier components close to zenith

      IF ( DO_USER_STREAMS ) THEN
        IF ( LOCAL_UM_START .GT. 1 ) THEN
          DO UTA = 1, N_USER_LEVELS
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, LOCAL_UM_START - 1
                DO O1 = 1, NSTOKES
                  STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_TOA_SOURCE ( &
          DO_OBSERVATION_GEOMETRY, NSTOKES, &
          N_USER_STREAMS, LOCAL_UM_START, IBEAM, &
          TOA_SOURCE )

        NC = 0
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              CUMSOURCE_DN(UM,O1,NC) = TOA_SOURCE(UM,O1)
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            CUMSOURCE_DN(IBEAM,O1,NC) = TOA_SOURCE(IBEAM,O1)
          ENDDO
        ENDIF
      ENDIF

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = N

            CALL WHOLELAYER_STERM_DN ( &
              N, IBEAM, SFLAG, &
              DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
              DO_OBSERVATION_GEOMETRY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_DNDN, UHOM_DNUP, &
              UPAR_DN_1, UPAR_DN_2, &
              HMULT_1, HMULT_2, &
              EMULT_DN, &
              LAYER_TSUP_DN, &
              LAYER_SOURCE )

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  CUMSOURCE_DN(UM,O1,NC) = LAYER_SOURCE(UM,O1) + &
                            T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,O1,NC-1)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                CUMSOURCE_DN(IBEAM,O1,NC) = LAYER_SOURCE(IBEAM,O1) + &
                          T_DELT_USERM(N,IBEAM)*CUMSOURCE_DN(IBEAM,O1,NC-1)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_DN ( &
              N, UT, IBEAM, SFLAG, &
              DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
              DO_OBSERVATION_GEOMETRY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_DNDN, UHOM_DNUP, &
              UPAR_DN_1, UPAR_DN_2, &
              UT_HMULT_DU, UT_HMULT_DD, &
              UT_EMULT_DN, &
              LAYER_TSUP_UTDN, &
              LAYER_SOURCE )

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  FINAL_SOURCE = LAYER_SOURCE(UM,O1) + &
                               T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,O1,NC)
                  STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = &
                               FLUX_MULTIPLIER * FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                FINAL_SOURCE = LAYER_SOURCE(IBEAM,O1) + &
                             T_UTDN_USERM(UT,IBEAM)*CUMSOURCE_DN(IBEAM,O1,NC)
                STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = &
                             FLUX_MULTIPLIER * FINAL_SOURCE
              ENDDO
            ENDIF
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = &
                         FLUX_MULTIPLIER * CUMSOURCE_DN(UM,O1,NC)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = &
                       FLUX_MULTIPLIER * CUMSOURCE_DN(IBEAM,O1,NC)
              ENDDO
            ENDIF
          ENDIF

        ENDIF

!        IF ( DO_DEBUG_WRITE ) THEN
!          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!            DO UM = LOCAL_UM_START, N_USER_STREAMS
!              DO O1 = 1, NSTOKES
!                if (do_fdtest) then
!                  if (uta.gt.1)write(*,'(5i3,1pe20.10)')
!     &              fourier_component,uta-1,ibeam,um,o1,
!     &              STOKES_F(UTA,UM,IBEAM,O1,DNIDX)
!                endif
!              ENDDO
!            ENDDO
!          ELSE
!            DO O1 = 1, NSTOKES
!              if (do_fdtest) then
!                if (uta.gt.1)write(*,'(5i3,1pe20.10)')
!     &            fourier_component,uta-1,ibeam,lum,o1,
!     &            STOKES_F(UTA,LUM,IBEAM,O1,DNIDX)
!              endif
!            ENDDO
!          ENDIF
!        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_DNUSER_INTENSITY

!

      SUBROUTINE GET_TOA_SOURCE ( &
        DO_OBSERVATION_GEOMETRY, NSTOKES, &
        N_USER_STREAMS, LOCAL_UM_START, IBEAM, &
        TOA_SOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (OUT) :: TOA_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      INTEGER :: UM, O1

!  initialise TOA source function

      if ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            TOA_SOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO
      else
        DO O1 = 1, NSTOKES
          TOA_SOURCE(IBEAM,O1) = ZERO
        ENDDO
      endif

!  Finish

      RETURN
      END SUBROUTINE GET_TOA_SOURCE

!

      SUBROUTINE GET_BOA_SOURCE ( DO_MSMODE_THERMAL, &
        DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT, &
        DO_OBSERVATION_GEOMETRY, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        DO_QUAD_OUTPUT, NSTOKES, &
        NSTREAMS, NLAYERS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        USER_BRDF_F, &
        EMISSIVITY, USER_EMISSIVITY, &
        SURFBB, DO_THERMAL_TRANSONLY, &
        DO_DBCORRECTION, QUAD_WEIGHTS, &
        QUAD_STRMWTS, N_USER_STREAMS, &
        MUELLER_INDEX, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WLOWER, LCON, MCON, &
        USER_DIRECT_BEAM, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        T_WLOWER, BRDF_F, &
        STOKES_DOWNSURF, BOA_DIFFUSE_SOURCE, &
        BOA_DIRECT_SOURCE, BOA_THTONLY_SOURCE )

!  Bottom of the atmosphere source term, Lambertian or BRDF

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DOWNSURF &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_DIFFUSE_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_DIRECT_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, J, I, IR, IROW, UM, O1, O2, O11, OM
      INTEGER ::          K, KO1, K0, K1, K2, M, IB, LUM
      DOUBLE PRECISION :: REFLEC, S_REFLEC, STOTAL, KMULT
      DOUBLE PRECISION :: SPAR, HOM1, HOM2, SHOM_R, LOCAL_EMISS, FP_SBB
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR
      LOGICAL          :: DO_QTHTONLY
      DOUBLE PRECISION :: THELP ( MAXSTREAMS )

!  Local index

      IB = IBEAM
      LUM = 1

!  initialise boa source function

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            BOA_DIFFUSE_SOURCE(UM,O1) = ZERO
            BOA_DIRECT_SOURCE(UM,O1)  = ZERO
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          BOA_DIFFUSE_SOURCE(IB,O1) = ZERO
          BOA_DIRECT_SOURCE(IB,O1)  = ZERO
        ENDDO
      ENDIF

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. &
            ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

!  Layer index and offset

      N = NLAYERS

!mick thermal fix - added if condition
      !KO1 = K_REAL(N) + 1
      IF ( .NOT. DO_THERMAL_TRANSONLY ) KO1 = K_REAL(N) + 1

!  fourier component

      M = FOURIER_COMPONENT

!  1-1 element index

      O11 = 1

!  reflectance from surface
!  ------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN

!  Thermal transmittance solution, build from TOA downwards

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

!  Full solution: Downward intensity at computational angles (beam/homog
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

!  Downward Stokes vector at surface at computational angles
!    reflectance integrand  a(j).x(j).I(-j)

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
               LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - &
                         LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
               LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + &
                         LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
               MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - &
                         MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
               HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N) &
                       -LXR_CI*T_DELT_EIGEN(K2,N)
               HOM2CR = MXR_CR
               SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
             ENDDO
             STOTAL = SPAR + SHOM_R + SHOM_CR
             STOKES_DOWNSURF(I,O1) =  QUAD_STRMWTS(I) * STOTAL
           ENDDO
         ENDDO
        ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)
!   @@@ Rob Fix, 2/9/11, Position of DO_QTHTONLY clause is wrong

        IF ( DO_LAMBERTIAN_SURFACE ) THEN

          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
            O1 = 1
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
            ENDDO
!            IF ( DO_QTHTONLY ) THEN
!              DO I = 1, NSTREAMS
!                BOA_THTONLY_SOURCE(I,O1) = REFLEC
!              ENDDO
!            ENDIF
            REFLEC = KMULT * REFLEC
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                BOA_DIFFUSE_SOURCE(UM,O1) = REFLEC
              ENDDO
            ELSE
              BOA_DIFFUSE_SOURCE(IB,O1) = REFLEC
            ENDIF
            IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
                BOA_THTONLY_SOURCE(I,O1) = REFLEC
              ENDDO
            ENDIF
          ENDIF

!  ###### bidirectional reflectance
!   @@@ Rob Fix, 2/9/11, DO_QTHTONLY clause is wrong, and wrong position

        ELSE

          KMULT = SURFACE_FACTOR

          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * &
                                          USER_BRDF_F(M,OM,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
!                  IF ( DO_QTHTONLY ) THEN
!                    DO I = 1, NSTREAMS
!                      BOA_THTONLY_SOURCE(I,O1) = REFLEC
!                    ENDDO
!                  ENDIF
               ENDDO
               BOA_DIFFUSE_SOURCE(UM,O1) = KMULT * REFLEC
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * &
                                        USER_BRDF_F(M,OM,IB,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
!                IF ( DO_QTHTONLY ) THEN
!                  DO I = 1, NSTREAMS
!                    BOA_THTONLY_SOURCE(I,O1) = REFLEC
!                  ENDDO
!                ENDIF
             ENDDO
             BOA_DIFFUSE_SOURCE(IB,O1) = KMULT * REFLEC
            ENDDO
          ENDIF

          IF ( DO_QTHTONLY ) THEN
            O11 = 1
            DO I = 1, NSTREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + STOKES_DOWNSURF(J,O11) * &
                                   BRDF_F(M,O11,I,J)
              ENDDO
              BOA_THTONLY_SOURCE(I,O11) = KMULT * REFLEC
            ENDDO
          ENDIF

        ENDIF

!  Add direct beam if flagged
!   Addition of the DBCORRECTION clause is now required

        IF ( DO_INCLUDE_DIRECTBEAM .AND. .NOT.DO_DBCORRECTION ) THEN
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                BOA_DIRECT_SOURCE(UM,O1) = BOA_DIRECT_SOURCE(UM,O1) + &
                    USER_DIRECT_BEAM(UM,IB,O1)
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              BOA_DIRECT_SOURCE(IB,O1) = BOA_DIRECT_SOURCE(IB,O1) + &
                  USER_DIRECT_BEAM(LUM,IB,O1)
            ENDDO
          ENDIF
        ENDIF

!  End inclusion of surface terms

      ENDIF

!  Add surface emission term if flagged

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        FP_SBB = SURFBB ; if ( DO_MSMODE_THERMAL ) go to 23

        if ( do_lambertian_surface ) then
          o1 = 1
          local_emiss = FP_SBB * ( one - lambertian_albedo )
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) &
                                          + LOCAL_EMISS
            ENDDO
          ELSE
            BOA_DIFFUSE_SOURCE(IB,O1) = BOA_DIFFUSE_SOURCE(IB,O1) &
                                           + LOCAL_EMISS
          ENDIF
        else
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
!  @@@ Rob fix, ordering of USER_EMISSIVITY indices wrong
!               BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) &
!                                         + FP_SBB*USER_EMISSIVITY(UM,O1)
                BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) &
                                          + FP_SBB*USER_EMISSIVITY(O1,UM)
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              BOA_DIFFUSE_SOURCE(IB,O1) = BOA_DIFFUSE_SOURCE(IB,O1) &
                                           + FP_SBB*USER_EMISSIVITY(O1,IB)
            ENDDO
          ENDIF
        endif

23      IF ( DO_QTHTONLY ) THEN
          if ( do_lambertian_surface ) then
            o1 = 1
            local_emiss = FP_SBB * ( one - lambertian_albedo )
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I,O1) = &
                 BOA_THTONLY_SOURCE(I,O1)+LOCAL_EMISS
            ENDDO
          else
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
!  @@@ Rob fix, ordering of EMISSIVITY indices wrong
!               BOA_THTONLY_SOURCE(I,O1) = BOA_THTONLY_SOURCE(I,O1) + &
!                      FP_SBB * EMISSIVITY(I,O1)
                BOA_THTONLY_SOURCE(I,O1) = BOA_THTONLY_SOURCE(I,O1) + &
                       FP_SBB * EMISSIVITY(O1,I)
              ENDDO
            ENDDO
          endif
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_BOA_SOURCE

!

      SUBROUTINE WHOLELAYER_STERM_UP ( &
        N, M, &
        IBEAM, SOURCETERM_FLAG, &
        DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
        DO_OBSERVATION_GEOMETRY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_UPDN, UHOM_UPUP, &
        UPAR_UP_1, UPAR_UP_2, &
        HMULT_1, HMULT_2, &
        EMULT_UP, &
        LAYER_TSUP_UP, &
        LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           N, M, IBEAM
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::           NSTOKES
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, M1, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  Fourier number (debug only)

      M1 = M
      L1 = SOURCETERM_FLAG

!  Local user indices

      IB  = IBEAM
      LUM = 1

!   Very important to zero output term (bug solved 04/20/05)

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          LAYERSOURCE(IB,O1) = ZERO
        ENDDO
      ENDIF

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  If no source term return
!  @@@ Moved here, after the go to 6789

      IF ( .not. L1 ) RETURN

!  Offset

      KO1 = K_REAL(N) + 1

!  Whole layer source function Homogeneous only

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            H_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
              MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
              H_R = H_R + LUXR * HMULT_2(K,UM,N) + &
                          MUXR * HMULT_1(K,UM,N)
            ENDDO
            H_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - &
                       LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
              LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + &
                       LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
              MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - &
                       MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
              MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + &
                       MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
              H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) - &
                            LUX_CI * HMULT_2(K2,UM,N) + &
                            MUX_CR * HMULT_1(K1,UM,N) - &
                            MUX_CI * HMULT_1(K2,UM,N)
            ENDDO
            LAYERSOURCE(UM,O1) = H_R + H_CR
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_UPDN(IB,O1,K,N)
            MUXR = MCON(K,N) * UHOM_UPUP(IB,O1,K,N)
            H_R = H_R + LUXR * HMULT_2(K,IB,N) + &
                        MUXR * HMULT_1(K,IB,N)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_UPDN(IB,O1,K1,N) - &
                     LCON(K2,N) * UHOM_UPDN(IB,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_UPDN(IB,O1,K2,N) + &
                     LCON(K2,N) * UHOM_UPDN(IB,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_UPUP(IB,O1,K1,N) - &
                     MCON(K2,N) * UHOM_UPUP(IB,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_UPUP(IB,O1,K2,N) + &
                     MCON(K2,N) * UHOM_UPUP(IB,O1,K1,N)
            H_CR = H_CR + LUX_CR * HMULT_2(K1,IB,N) - &
                          LUX_CI * HMULT_2(K2,IB,N) + &
                          MUX_CR * HMULT_1(K1,IB,N) - &
                          MUX_CI * HMULT_1(K2,IB,N)
          ENDDO
          LAYERSOURCE(IB,O1) = H_R + H_CR
        ENDDO
      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        O1 = 1
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) &
                               + LAYER_TSUP_UP(UM,N)*TM
          ENDDO
        ELSE
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) &
                             + LAYER_TSUP_UP(IB,N)*TM
        ENDIF
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add particular integral contribution (solar term)

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR2 =  EMULT_UP(UM,N,IB) * UPAR_UP_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR2 =  EMULT_UP(LUM,N,IB) * UPAR_UP_2(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR2
          ENDDO
        ENDIF
      ELSE
!  PLACEHOLDER FOR GREENS FUNCT
      ENDIF

!  If operating in Ms-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

!  Full radiance mode, add single scatter part

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR1 = UPAR_UP_1(UM,O1,N) * EMULT_UP(UM,N,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          SFOR1 = UPAR_UP_1(IB,O1,N) * EMULT_UP(LUM,N,IB)
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR1
        ENDDO
      ENDIF

!  debug. NSTOKES = 1 checks out.

!      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!        !write(*,*) 'at dbg 1a'
!        if ( ib.eq.1 .and. o1.eq.1 .and. um.eq.1 ) then
!          DO UM = LOCAL_UM_START, N_USER_STREAMS
!            write(*,*)N,um, LAYERSOURCE(UM,O1), EMULT_UP(UM,N,IB)
!          ENDDO
!       endif
!      ELSE
!        if ( ib.eq.1 .and. o1.eq.1) then
!         write(*,*)N,ib, LAYERSOURCE(IB,O1), EMULT_UP(LUM,N,IB)
!        endif
!      ENDIF
!      if (n.eq.1)pause

!  Finish

      RETURN
      END SUBROUTINE WHOLELAYER_STERM_UP

!

      SUBROUTINE WHOLELAYER_STERM_DN ( &
        N, IBEAM, SOURCETERM_FLAG, &
        DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
        DO_OBSERVATION_GEOMETRY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        HMULT_1, HMULT_2, &
        EMULT_DN, &
        LAYER_TSUP_DN, &
        LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N, IBEAM
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS

      LOGICAL, INTENT (IN) ::             DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  Local user indices

      IB  = IBEAM
      LUM = 1

!  Zeroing is very important here

      L1 = SOURCETERM_FLAG
      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          LAYERSOURCE(IB,O1) = ZERO
        ENDDO
      ENDIF

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  If no source term return
!  @@@ Moved here, after the go to 6789

      IF ( .not. L1 ) RETURN

!  Offset

      KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous solution contribution

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            H_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
              MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
              H_R = H_R + LUXR * HMULT_1(K,UM,N) + &
                          MUXR * HMULT_2(K,UM,N)
            ENDDO
            H_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - &
                       LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
              LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + &
                       LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
              MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - &
                       MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
              MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + &
                       MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
              H_CR = H_CR + LUX_CR * HMULT_1(K1,UM,N) - &
                            LUX_CI * HMULT_1(K2,UM,N) + &
                            MUX_CR * HMULT_2(K1,UM,N) - &
                            MUX_CI * HMULT_2(K2,UM,N)
            ENDDO
            LAYERSOURCE(UM,O1) = H_R + H_CR
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_DNDN(IB,O1,K,N)
            MUXR = MCON(K,N) * UHOM_DNUP(IB,O1,K,N)
            H_R = H_R + LUXR * HMULT_1(K,IB,N) + &
                        MUXR * HMULT_2(K,IB,N)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_DNDN(IB,O1,K1,N) - &
                     LCON(K2,N) * UHOM_DNDN(IB,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_DNDN(IB,O1,K2,N) + &
                     LCON(K2,N) * UHOM_DNDN(IB,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_DNUP(IB,O1,K1,N) - &
                     MCON(K2,N) * UHOM_DNUP(IB,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_DNUP(IB,O1,K2,N) + &
                     MCON(K2,N) * UHOM_DNUP(IB,O1,K1,N)
            H_CR = H_CR + LUX_CR * HMULT_1(K1,IB,N) - &
                          LUX_CI * HMULT_1(K2,IB,N) + &
                          MUX_CR * HMULT_2(K1,IB,N) - &
                          MUX_CI * HMULT_2(K2,IB,N)
          ENDDO
          LAYERSOURCE(IB,O1) = H_R + H_CR
        ENDDO
      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        O1 = 1
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) &
                                 + LAYER_TSUP_DN(UM,N)*TM
          ENDDO
        ELSE
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) &
                               + LAYER_TSUP_DN(IB,N)*TM
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR2 =  EMULT_DN(UM,N,IB) * UPAR_DN_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR2 =  EMULT_DN(LUM,N,IB) * UPAR_DN_2(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR2
          ENDDO
        ENDIF
      ELSE
!  PLACEHOLDER GREENS FUNCT
      ENDIF

!  If operating in Ms-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

!  Full radiance mode, add single scatter part

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            SFOR1 = UPAR_DN_1(UM,O1,N) * EMULT_DN(UM,N,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          SFOR1 = UPAR_DN_1(IB,O1,N) * EMULT_DN(LUM,N,IB)
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR1
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE WHOLELAYER_STERM_DN

!

      SUBROUTINE PARTLAYER_STERM_UP ( &
        N, UT, IBEAM, &
        SOURCETERM_FLAG, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
        DO_OBSERVATION_GEOMETRY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_UPDN, UHOM_UPUP, &
        UPAR_UP_1, UPAR_UP_2, &
        UT_HMULT_UU, UT_HMULT_UD, &
        UT_EMULT_UP, &
        LAYER_TSUP_UTUP, &
        LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N, UT, IBEAM
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  Local user indices

      IB = IBEAM
      LUM = 1

!  Zeroing is very important here

      SFOR2 = ZERO
      L1 = SOURCETERM_FLAG

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          LAYERSOURCE(IB,O1) = ZERO
        ENDDO
      ENDIF

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  If no source term return
!  @@@ Moved here, after the goto 6789

      IF ( .not. L1 ) RETURN

!  Offset

      KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous contribution

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            H_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
              MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
              H_R = H_R + LUXR * UT_HMULT_UD(K,UM,UT) + &
                          MUXR * UT_HMULT_UU(K,UM,UT)
            ENDDO
            H_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - &
                       LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
              LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + &
                       LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
              MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - &
                       MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
              MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + &
                       MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
              H_CR = H_CR + LUX_CR * UT_HMULT_UD(K1,UM,UT) - &
                            LUX_CI * UT_HMULT_UD(K2,UM,UT) + &
                            MUX_CR * UT_HMULT_UU(K1,UM,UT) - &
                            MUX_CI * UT_HMULT_UU(K2,UM,UT)
            ENDDO
            LAYERSOURCE(UM,O1) = H_R + H_CR
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_UPDN(IB,O1,K,N)
            MUXR = MCON(K,N) * UHOM_UPUP(IB,O1,K,N)
            H_R = H_R + LUXR * UT_HMULT_UD(K,IB,UT) + &
                        MUXR * UT_HMULT_UU(K,IB,UT)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_UPDN(IB,O1,K1,N) - &
                     LCON(K2,N) * UHOM_UPDN(IB,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_UPDN(IB,O1,K2,N) + &
                     LCON(K2,N) * UHOM_UPDN(IB,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_UPUP(IB,O1,K1,N) - &
                     MCON(K2,N) * UHOM_UPUP(IB,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_UPUP(IB,O1,K2,N) + &
                     MCON(K2,N) * UHOM_UPUP(IB,O1,K1,N)
            H_CR = H_CR + LUX_CR * UT_HMULT_UD(K1,IB,UT) - &
                          LUX_CI * UT_HMULT_UD(K2,IB,UT) + &
                          MUX_CR * UT_HMULT_UU(K1,IB,UT) - &
                          MUX_CI * UT_HMULT_UU(K2,IB,UT)
          ENDDO
          LAYERSOURCE(IB,O1) = H_R + H_CR
        ENDDO
      ENDIF

!  Continuation point

 6789 continue

!  Add thermal term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        O1 = 1
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + &
                           LAYER_TSUP_UTUP(UM,UT)*TM
          ENDDO
        ELSE
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + &
                         LAYER_TSUP_UTUP(IB,UT)*TM
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR2 =  UT_EMULT_UP(UM,UT,IBEAM) * UPAR_UP_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR2 =  UT_EMULT_UP(LUM,UT,IB) * UPAR_UP_2(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR2
          ENDDO
        ENDIF
      ELSE
!  PLACEHOLDER GREENS FUNCT
      ENDIF

!  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR1 = UT_EMULT_UP(UM,UT,IBEAM) * UPAR_UP_1(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR1 = UT_EMULT_UP(LUM,UT,IB) * UPAR_UP_1(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR1
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_STERM_UP

!

      SUBROUTINE PARTLAYER_STERM_DN ( &
        N, UT, IBEAM, SOURCETERM_FLAG, &
        DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
        DO_OBSERVATION_GEOMETRY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        UT_HMULT_DU, UT_HMULT_DD, &
        UT_EMULT_DN, &
        LAYER_TSUP_UTDN, &
        LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N, UT, IBEAM
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  Local user indices

      IB  = IBEAM
      LUM = 1

!  Zeroing is very important here

      SFOR2 = ZERO
      L1 = SOURCETERM_FLAG
      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          LAYERSOURCE(IB,O1) = ZERO
        ENDDO
      ENDIF

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  If no source term return
!  @@@ Moved here, after the goto 6789

      IF ( .not. L1 ) RETURN

!  Offset

      KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous contribution

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            H_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
              MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
              H_R = H_R + LUXR * UT_HMULT_DD(K,UM,UT) + &
                          MUXR * UT_HMULT_DU(K,UM,UT)
            ENDDO
            H_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - &
                       LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
              LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + &
                       LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
              MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - &
                       MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
              MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + &
                       MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
              H_CR = H_CR + LUX_CR * UT_HMULT_DD(K1,UM,UT) - &
                            LUX_CI * UT_HMULT_DD(K2,UM,UT) + &
                            MUX_CR * UT_HMULT_DU(K1,UM,UT) - &
                            MUX_CI * UT_HMULT_DU(K2,UM,UT)
            ENDDO
            LAYERSOURCE(UM,O1) = H_R + H_CR
          ENDDO
        ENDDO
      ELSE
        DO O1 = 1, NSTOKES
          H_R = ZERO
          DO K = 1, K_REAL(N)
            LUXR = LCON(K,N) * UHOM_DNDN(IB,O1,K,N)
            MUXR = MCON(K,N) * UHOM_DNUP(IB,O1,K,N)
            H_R = H_R + LUXR * UT_HMULT_DD(K,IB,UT) + &
                        MUXR * UT_HMULT_DU(K,IB,UT)
          ENDDO
          H_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LUX_CR = LCON(K1,N) * UHOM_DNDN(IB,O1,K1,N) - &
                     LCON(K2,N) * UHOM_DNDN(IB,O1,K2,N)
            LUX_CI = LCON(K1,N) * UHOM_DNDN(IB,O1,K2,N) + &
                     LCON(K2,N) * UHOM_DNDN(IB,O1,K1,N)
            MUX_CR = MCON(K1,N) * UHOM_DNUP(IB,O1,K1,N) - &
                     MCON(K2,N) * UHOM_DNUP(IB,O1,K2,N)
            MUX_CI = MCON(K1,N) * UHOM_DNUP(IB,O1,K2,N) + &
                     MCON(K2,N) * UHOM_DNUP(IB,O1,K1,N)
            H_CR = H_CR + LUX_CR * UT_HMULT_DD(K1,IB,UT) - &
                          LUX_CI * UT_HMULT_DD(K2,IB,UT) + &
                          MUX_CR * UT_HMULT_DU(K1,IB,UT) - &
                          MUX_CI * UT_HMULT_DU(K2,IB,UT)
          ENDDO
          LAYERSOURCE(IB,O1) = H_R + H_CR
        ENDDO
      ENDIF

!  Continuation point

 6789 continue

!  Add thermal term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        O1 = 1
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + &
                           LAYER_TSUP_UTDN(UM,UT)*TM
          ENDDO
        ELSE
          LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + &
                         LAYER_TSUP_UTDN(IB,UT)*TM
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular integral contribution

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR2 =  UT_EMULT_DN(UM,UT,IBEAM) * UPAR_DN_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR2 =  UT_EMULT_DN(LUM,UT,IB) * UPAR_DN_2(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR2
          ENDDO
        ENDIF
      ELSE
!  PLACEHOLDER GREENS FUNCT
      ENDIF

!  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SFOR1 = UT_EMULT_DN(UM,UT,IBEAM) * UPAR_DN_1(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            SFOR1 = UT_EMULT_DN(LUM,UT,IB) * UPAR_DN_1(IB,O1,N)
            LAYERSOURCE(IB,O1) = LAYERSOURCE(IB,O1) + SFOR1
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_STERM_DN

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument (3rd line)

      SUBROUTINE VLIDORT_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_THERMEMISS, &
        IBEAM, FLUX_MULTIPLIER, FLUX_FACTOR, &         !  added here line 3
        NSTOKES, NSTREAMS, &
        N_USER_LEVELS, &
        FLUXVEC, LAYER_PIS_CUTOFF, &
        QUAD_WEIGHTS, QUAD_STRMWTS, &
        N_DIRECTIONS, WHICH_DIRECTIONS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, &
        T_DELT_MUBAR, T_UTDN_MUBAR, &
        INITIAL_TRANS, LOCAL_CSZA, &
        DO_SOLAR_SOURCES, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, LCON, MCON, &
        T_WUPPER, UT_T_PARTIC, &
        BOA_THTONLY_SOURCE, &
        T_DELT_EIGEN, WLOWER, &
        T_DISORDS_UTDN, T_WLOWER, &
        MEAN_STOKES, FLUX_STOKES, &
        MEAN_DIRECT, FLUX_DIRECT )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument 
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: MEAN_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QSTOKES_F &
          ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  help variables

      INTEGER          :: I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1
      DOUBLE PRECISION :: SUM_MI, SUM_FX, FTRANS
      DOUBLE PRECISION :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Stokes output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADINTENS_OFFGRID_UP ( &
                N, UTA, UT, IBEAM, &
                FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                T_UTDN_MUBAR, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, &
                T_WUPPER, UT_T_PARTIC, &
                BOA_THTONLY_SOURCE, QSTOKES_F )

            ELSE
              CALL QUADINTENS_LEVEL_UP ( &
                NLEVEL, UTA, FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, &
                NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                QUAD_STREAMS, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, WLOWER, &
                LCON, MCON, &
                T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                QSTOKES_F )

            ENDIF
          ENDDO
        ENDIF

!  Downwelling Stokes output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADINTENS_OFFGRID_DN ( &
                N, UTA, UT, IBEAM, &
                FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, &
                DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTDN, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                T_UTDN_MUBAR, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, &
                LCON, MCON, &
                T_WUPPER, T_WLOWER, UT_T_PARTIC, &
                QSTOKES_F )

            ELSE
              CALL QUADINTENS_LEVEL_DN ( &
                NLEVEL, UTA, FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, &
                DO_THERMAL_TRANSONLY, &
                QUAD_STREAMS, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WLOWER, &
                LCON, MCON, &
                T_WLOWER, &
                QSTOKES_F )

            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux  output
!  -------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  Diffuse term output

          DO UTA = 1, N_USER_LEVELS
            DO O1 = 1, NSTOKES
              SUM_MI = ZERO
              SUM_FX = ZERO
              DO I = 1, NSTREAMS
                SUM_MI = SUM_MI + QUAD_WEIGHTS(I) * QSTOKES_F(UTA,I,O1)
                SUM_FX = SUM_FX + QUAD_STRMWTS(I) * QSTOKES_F(UTA,I,O1)
              ENDDO
              MEAN_STOKES(UTA,IBEAM,O1,WDIR) = SUM_MI * HALF
              FLUX_STOKES(UTA,IBEAM,O1,WDIR) = SUM_FX * PI2
           ENDDO
          ENDDO

!  nothing to do if no solar sources

          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

!  For the downward direction, add the direct beam contributions

          IF ( WDIR .EQ. DNIDX ) THEN

!  loop over all the output optical depths

           DO UTA = 1, N_USER_LEVELS

!  For the offgrid values.......
!     .....Only contributions for layers above the PI cutoff

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                DIRECT_TRANS = &
                   INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)
                DO O1 = 1, NSTOKES

!  @@@ Rob fix (forgot the Flux factor)
!                 FTRANS = FLUXVEC(O1) * DIRECT_TRANS
!                 FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!                  FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * DIRECT_TRANS
! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS

                  DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                  DIRECT_MEANI = FTRANS / PI4
                  MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
                  FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) = &
                        MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) = &
                        FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
                ENDDO
              ENDIF

!  For the on-grid values

             ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                IF ( N .EQ. 0 ) THEN
                  DIRECT_TRANS = ONE
                ELSE
                  DIRECT_TRANS = &
                     INITIAL_TRANS(N,IBEAM)*T_DELT_MUBAR(N,IBEAM)
                ENDIF
                DO O1 = 1, NSTOKES

!  @@@ Rob fix (forgot the Flux factor)
!                 FTRANS = FLUXVEC(O1) * DIRECT_TRANS
!                 FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!                  FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * DIRECT_TRANS
! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS

                  DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                  DIRECT_MEANI = FTRANS / PI4
                  MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
                  FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) = &
                        MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) = &
                        FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
                ENDDO
              ENDIF
             ENDIF
           ENDDO

!  Finish downwelling direct contribution

          ENDIF

!  Finish MV output

        ENDIF

!  Continuation point for avoiding direct beam calculation

 455    CONTINUE

!  Finish directional loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_INTEGRATED_OUTPUT

!

      SUBROUTINE QUADINTENS_LEVEL_UP ( &
        NLEVEL, UTA, FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, &
        NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        QUAD_STREAMS, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, WLOWER, &
        LCON, MCON, &
        T_WUPPER, &
        BOA_THTONLY_SOURCE, &
        QSTOKES_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NLEVEL, UTA
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F &
          ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, NL, I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      NL = NLEVEL
      N  = NL + 1

!  Lowest level contributions
!  ==========================

!  Lowest level, thermal transmittance only

      IF ( NL .EQ. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For the lowest level, scattering solution
!  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  Start main loops

        KO1 = K_REAL(NL) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              LXR = LCON(K,NL)*SOLA_XPOS(I1,O1,K,NL)
              MXR = MCON(K,NL)*SOLB_XNEG(I1,O1,K,NL)
              HOM1 = LXR * T_DELT_EIGEN(K,NL)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - &
                        LCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
              LXR_CI =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + &
                        LCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
              MXR_CR =  MCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - &
                        MCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
              HOM1CR =  LXR_CR*T_DELT_EIGEN(K1,NL) &
                      - LXR_CI*T_DELT_EIGEN(K2,NL)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I1,O1,NL)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO

!  End lowest level clause

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          DO K = NLAYERS, N, -1
            TPROP = T_WUPPER(I1,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For other levels, scattering solution
!  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

!  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
              MXR = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
              HOM1 = LXR
              HOM2 = MXR * T_DELT_EIGEN(K,N)
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - &
                        LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - &
                        MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
              MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + &
                        MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
              HOM1CR = LXR_CR
              HOM2CR =  MXR_CR * T_DELT_EIGEN(K1,N) &
                      - MXR_CI * T_DELT_EIGEN(K2,N)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_LEVEL_UP

!

      SUBROUTINE QUADINTENS_LEVEL_DN ( &
        NLEVEL, UTA, FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, &
        DO_THERMAL_TRANSONLY, &
        QUAD_STREAMS, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WLOWER, &
        LCON, MCON, &
        T_WLOWER, &
        QSTOKES_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NLEVEL, UTA
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F &
          ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL

!  Downwelling radiation at TOA ( or N = 0 ) is zero
!      A check has been performed on this.

      IF ( NLEVEL .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            QSTOKES_F(UTA,I,O1) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Other levels, Thermal transmittance-only solution

      IF ( NLEVEL .NE. 0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N
            TPROP = T_WLOWER(I,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Other levels, scattering solution
!  ---------------------------------

      IF ( NLEVEL .NE. 0 ) THEN

!  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
              MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
              HOM1 = LXR * T_DELT_EIGEN(K,N)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - &
                        LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
              LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + &
                        LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - &
                        MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N) &
                      -LXR_CI*T_DELT_EIGEN(K2,N)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_LEVEL_DN

!

      SUBROUTINE QUADINTENS_OFFGRID_UP ( &
        N, UTA, UT, IBEAM, &
        FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_UTDN_MUBAR, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, LCON, MCON, &
        T_WUPPER, UT_T_PARTIC, &
        BOA_THTONLY_SOURCE, QSTOKES_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N, UTA, UT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F &
          ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          DO K = NLAYERS, N + 1, -1
            TPROP = T_WUPPER(I1,K) / QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP =  UT_T_PARTIC(I1,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTUP(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

!  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - &
                      LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
            LXR_CI =  LCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + &
                      LCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
            MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - &
                      MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + &
                      MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT) &
                      - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT) &
                      - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

!  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

!  Finish streams/stokes loops

        ENDDO
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
           QSTOKES_F(UTA,I,O1) = &
             QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add the solar particular solutions
!     Green function solution is still a Placeholder....!

! @@@@@@@@@Rob Fix, 2/9/11, Replacement of the following code
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            SPAR = WUPPER(I1,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
!            QSTOKES_F(UTA,I,O1) = &
!              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!          ENDDO
!        ENDDO
!      ELSE
!C                P L A C E H O L D E R
!      ENDIF

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              SPAR = WUPPER(I1,O1,N)
              IF ( O1.eq.1) SPAR = SPAR - T_WUPPER(I1,N)
              SPAR = SPAR * T_UTDN_MUBAR(UT,IBEAM)
              QSTOKES_F(UTA,I,O1) = &
              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              SPAR = WUPPER(I1,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
              QSTOKES_F(UTA,I,O1) = &
              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF
! @@@@@@@@@Rob Fix, 2/9/11, END of Replacement

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_OFFGRID_UP

!

      SUBROUTINE QUADINTENS_OFFGRID_DN ( &
        N, UTA, UT, IBEAM, &
        FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTDN, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_UTDN_MUBAR, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, &
        LCON, MCON, &
        T_WUPPER, T_WLOWER, UT_T_PARTIC, &
        QSTOKES_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N, UTA, UT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F &
          ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N - 1
            TPROP = T_WLOWER(I,K) / QUAD_STREAMS(I)
            THELP = THELP*T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP = UT_T_PARTIC(I,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTDN(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

!  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - &
                      LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
            LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + &
                      LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
            MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - &
                      MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + &
                      MCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT) &
                      - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT) &
                      - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

!  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

!  Finish streams/stokes loops

        ENDDO
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
           QSTOKES_F(UTA,I,O1) = &
             QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add the solar particular solutions
!     Green function solution is still a Placeholder....!

! @@@@@@@@@Rob Fix, 2/9/11, Replacement of the following code
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        DO I = 1, NSTREAMS
!          DO O1 = 1, NSTOKES
!            SPAR = WUPPER(I,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
!            QSTOKES_F(UTA,I,O1) = &
!              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!          ENDDO
!        ENDDO
!      ELSE
!C                P L A C E H O L D E R
!      ENDIF

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              SPAR = WUPPER(I,O1,N)
              IF ( O1.eq.1) SPAR = SPAR - T_WUPPER(I,N)
              SPAR = SPAR * T_UTDN_MUBAR(UT,IBEAM)
              QSTOKES_F(UTA,I,O1) = &
              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              SPAR = WUPPER(I,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
              QSTOKES_F(UTA,I,O1) = &
              QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF
! @@@@@@@@@Rob Fix, 2/9/11, END of Replacement

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_OFFGRID_DN

!

      SUBROUTINE VLIDORT_CONVERGE ( &
        AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SS_EXTERNAL, & ! New 15 March 2012
        DO_SSFULL, DO_DOUBLE_CONVTEST, &
        DO_RAYLEIGH_ONLY, DO_UPWELLING, &
        NSTOKES, NSTREAMS, &
        NLAYERS, VLIDORT_ACCURACY, &
        N_USER_RELAZMS, N_USER_LEVELS, &
        DO_TOA_CONTRIBS, &
        DO_ALL_FOURIER, DO_DBCORRECTION, &
        DO_NO_AZIMUTH, N_CONVTESTS, &
        LOCAL_UM_START, N_DIRECTIONS, &
        WHICH_DIRECTIONS, &
        N_OUT_STREAMS, VZA_OFFSETS, &
        STOKES_SS, STOKES_DB, SS_CONTRIBS, &
        STOKES_F, MS_CONTRIBS_F, IB, &
        STOKES, FOURIER_SAVED, CONTRIBS, &
        TESTCONV, LOCAL_ITERATION )

!  convergence test on the Stokes vector intensity

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           LOCAL_N_USERAZM
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC &
          ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::           DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::           DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL, INTENT (IN) ::           DO_SS_EXTERNAL
      LOGICAL, INTENT (IN) ::           DO_SSFULL
      LOGICAL, INTENT (IN) ::           DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::           N_CONVTESTS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::           N_OUT_STREAMS
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      INTEGER, INTENT (IN) ::           IB

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      INTEGER, INTENT (INOUT) ::          FOURIER_SAVED ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (INOUT) :: CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          TESTCONV
      LOGICAL, INTENT (INOUT) ::          LOCAL_ITERATION

!  local variables
!  ---------------

      INTEGER ::          COUNT, COUNT_A, V, N
      INTEGER ::          I, IDIR, UT, UA, W, O
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, &
                          TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W) = STOKES_F(UT,I,IB,O,W)
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W) = ZERO
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_SSFULL ) THEN
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = MS_CONTRIBS_F(I,IB,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
         DO IDIR = 1, N_DIRECTIONS
          W = WHICH_DIRECTIONS(IDIR)
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IB,I) + UA
                DO O = 1, NSTOKES
!         if (v.eq.1)write(*,*)UT,STOKES(UT,V,O,W),STOKES_SS(UT,V,1,W)
                  STOKES(UT,V,O,W) = &
                    STOKES(UT,V,O,W) + STOKES_SS(UT,V,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDIF
!        pause

!    CONTRIBS. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = &
                      CONTRIBS(V,O,N) + SS_CONTRIBS(V,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!        pause

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT
!    Full single scatter option added Version 2.3

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  STOKES(UT,V,O,UPIDX) = &
                    STOKES(UT,V,O,UPIDX) + STOKES_DB(UT,V,O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    TOLD(O) = STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IB,UA,O)*STOKES_F(UT,I,IB,O,W)
                    STOKES(UT,V,O,W) = TOLD(O) + TAZM(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    TOLD(O) = STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IB,UA,O)*STOKES_F(UT,I,IB,O,W)
                    TNEW(O) = TOLD(O) + TAZM(O)
                  ENDDO
                  IF ( TAZM(1) .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM(1)/TNEW(1))
                    IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W) = TNEW(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IB) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IB) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO UA = 1, LOCAL_N_USERAZM
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              DO N = 1, NLAYERS
                DO O = 1, NSTOKES
                  TOLD(O) = CONTRIBS(V,O,N)
                  TAZM(O) = AZMFAC(I,IB,UA,O)*MS_CONTRIBS_F(I,IB,O,N)
                  CONTRIBS(V,O,N) = TOLD(O) + TAZM(O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IB) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE

!

      SUBROUTINE VLIDORT_CONVERGE_OBSGEO ( &
        AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SS_EXTERNAL, & ! New 15 March 2012
        DO_SSFULL, DO_DOUBLE_CONVTEST, &
        DO_RAYLEIGH_ONLY, DO_UPWELLING, &
        NSTOKES, NSTREAMS, &
        NLAYERS, VLIDORT_ACCURACY, &
        N_USER_RELAZMS, N_USER_LEVELS, &
        DO_TOA_CONTRIBS, &
        DO_ALL_FOURIER, DO_DBCORRECTION, &
        DO_NO_AZIMUTH, N_CONVTESTS, &
        LOCAL_UM_START, N_DIRECTIONS, &
        WHICH_DIRECTIONS, VZA_OFFSETS, &
        STOKES_SS, STOKES_DB, SS_CONTRIBS, &
        STOKES_F, MS_CONTRIBS_F, IB, &
        STOKES, FOURIER_SAVED, CONTRIBS, &
        TESTCONV, LOCAL_ITERATION )

!  convergence test on the Stokes vector intensity

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           LOCAL_N_USERAZM
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC &
          ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::           DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::           DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL, INTENT (IN) ::           DO_SS_EXTERNAL
      LOGICAL, INTENT (IN) ::           DO_SSFULL
      LOGICAL, INTENT (IN) ::           DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::           N_CONVTESTS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      INTEGER, INTENT (IN) ::           IB

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      INTEGER, INTENT (INOUT) ::          FOURIER_SAVED ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (INOUT) :: CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          TESTCONV
      LOGICAL, INTENT (INOUT) ::          LOCAL_ITERATION

!  local variables
!  ---------------

      INTEGER ::          COUNT, N
      INTEGER ::          I, IDIR, UT, UA, W, O, LUM, LUA
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, &
                          TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                STOKES(UT,IB,O,W) = STOKES_F(UT,LUM,IB,O,W)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                STOKES(UT,IB,O,W) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_SSFULL ) THEN
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                CONTRIBS(IB,O,N) = MS_CONTRIBS_F(LUM,IB,O,N)
              ENDDO
            ENDDO
          ELSE
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                CONTRIBS(IB,O,N) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
!          if (ib.eq.1)write(*,*)UT,STOKES(UT,IB,O,W),STOKES_SS(UT,IB,O,W)
                STOKES(UT,IB,O,W) = STOKES(UT,IB,O,W) + STOKES_SS(UT,IB,O,W)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!        pause

!    CONTRIBS. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                CONTRIBS(IB,O,N) = CONTRIBS(IB,O,N) + SS_CONTRIBS(IB,O,N)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!        pause

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT
!    Full single scatter option added Version 2.3

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO O = 1, NSTOKES
              STOKES(UT,IB,O,UPIDX) = STOKES(UT,IB,O,UPIDX) + STOKES_DB(UT,IB,O)
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each geometry, add Fourier component
!     - for direction, user optical depth

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO O = 1, NSTOKES
                  TOLD(O) = STOKES(UT,IB,O,W)
                  TAZM(O) = AZMFAC(LUM,IB,LUA,O)*STOKES_F(UT,LUM,IB,O,W)
                  STOKES(UT,IB,O,W) = TOLD(O) + TAZM(O)
                ENDDO
              ENDDO
            ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                TOLD(O) = STOKES(UT,IB,O,W)
                TAZM(O) = AZMFAC(LUM,IB,LUA,O)*STOKES_F(UT,LUM,IB,O,W)
                TNEW(O) = TOLD(O) + TAZM(O)
              ENDDO
              IF ( TAZM(1) .NE. ZERO ) THEN
                ACCUR = DABS(TAZM(1)/TNEW(1))
                IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                  COUNT = COUNT + 1
                ENDIF
              ELSE
                COUNT = COUNT + 1
              ENDIF
              DO O = 1, NSTOKES
                STOKES(UT,IB,O,W) = TNEW(O)
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
              LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IB) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IB) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO N = 1, NLAYERS
            DO O = 1, NSTOKES
              TOLD(O) = CONTRIBS(IB,O,N)
              TAZM(O) = AZMFAC(LUM,IB,LUA,O)*MS_CONTRIBS_F(LUM,IB,O,N)
              CONTRIBS(IB,O,N) = TOLD(O) + TAZM(O)
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IB) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE_OBSGEO

      END MODULE vlidort_intensity
