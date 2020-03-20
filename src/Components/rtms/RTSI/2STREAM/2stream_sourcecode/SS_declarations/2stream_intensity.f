C ###########################################################
C #                                                         #
C #             THE TWOSTREAM LIDORT MODEL                  #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Authors :      Robert. J. D. Spurr (1)                 #
C #                 Vijay Natraj        (2)                 #
C #                                                         #
C #  Address (1) :     RT Solutions, Inc.                   #
C #                    9 Channing Street                    #
C #                    Cambridge, MA 02138, USA             #
C #  Tel:             (617) 492 1183                        #
C #  Email :           rtsolutions@verizon.net              #
C #                                                         #
C #  Address (2) :     CalTech                              #
C #                    Department of Planetary Sciences     #
C #                    1200 East California Boulevard       #
C #                    Pasadena, CA 91125                   #
C #  Tel:             (626) 395 6962                        #
C #  Email :           vijay@gps.caltech.edu                #
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
C #          TWOSTREAM_UPUSER_INTENSITY   (master)          #
C #          TWOSTREAM_DNUSER_INTENSITY   (master)          #
C #          TWOSTREAM_CONVERGE (master)                    #
C #                                                         #
C ###########################################################

      SUBROUTINE TWOSTREAM_UPUSER_INTENSITY
     I      ( DO_INCLUDE_SURFACE, DO_MSMODE_LIDORT,
     I        DO_INCLUDE_DIRECTBEAM, DO_DB_CORRECTION,
     I        FOURIER_COMPONENT, NLAYERS, NBEAMS, N_USER_STREAMS,
     I        IPARTIC,  FLUX_MULTIPLIER,
     I        SURFACE_FACTOR, SURFTYPE, ALBEDO, USER_DIRECT_BEAM,
     I        T_DELT_EIGEN, T_DELT_USERM, STREAM_VALUE,
     I        XPOS, XNEG, WUPPER, WLOWER,
     I        LCON, LCON_XVEC, MCON, MCON_XVEC,
     I        U_XPOS, U_XNEG, U_WPOS1, U_WPOS2,
     I        HMULT_1, HMULT_2, EMULT_UP,
     I        BOA_SOURCE, DIRECT_BOA_SOURCE, IDOWNSURF,
     O        INTENSITY_F_UP, CUMSOURCE_UP )

C  Subroutine input arguments
C  --------------------------

C  Fourier component, beam index

      INTEGER          FOURIER_COMPONENT
      INTEGER          IPARTIC

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Surface stuff

      INTEGER          SURFTYPE
      DOUBLE PRECISION SURFACE_FACTOR, ALBEDO
      DOUBLE PRECISION USER_DIRECT_BEAM ( N_USER_STREAMS, NBEAMS )

C  local control flags

      LOGICAL          DO_MSMODE_LIDORT
      LOGICAL          DO_DB_CORRECTION
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN(NLAYERS)

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Solution constants of integration

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Solution constants of integration multiplied by homogeneous solutions

      DOUBLE PRECISION LCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION MCON_XVEC(2,NLAYERS)

C  General beam solutions at the Upper/Lower boundary

      DOUBLE PRECISION WUPPER(2,NLAYERS)
      DOUBLE PRECISION WLOWER(2,NLAYERS)

C  Eigenvectors defined at user-defined stream angles
C     EP for the positive KEIGEN values, EM for -ve KEIGEN

      DOUBLE PRECISION
     U        U_XPOS(N_USER_STREAMS,NLAYERS),
     U        U_XNEG(N_USER_STREAMS,NLAYERS)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION
     U        U_WPOS1(N_USER_STREAMS,NLAYERS),
     U        U_WPOS2(N_USER_STREAMS,NLAYERS)

C  solution multipliers 

      DOUBLE PRECISION 
     &      HMULT_1(N_USER_STREAMS,NLAYERS),
     &      HMULT_2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  BOA source terms

      DOUBLE PRECISION BOA_SOURCE        ( N_USER_STREAMS )
      DOUBLE PRECISION DIRECT_BOA_SOURCE ( N_USER_STREAMS )

C  Reflectance integrand  a(j).x(j).I(-j)

      DOUBLE PRECISION IDOWNSURF

C  Outputs
C  -------

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_UP
     D   (N_USER_STREAMS,NBEAMS)

C  Cumulative source terms

      DOUBLE PRECISION
     U    CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)

C  local variables
C  ---------------

C  help variables

      INTEGER          UM, N, NC
      DOUBLE PRECISION LAYERSOURCE       ( N_USER_STREAMS )
      DOUBLE PRECISION SHOM, SFOR1, SFOR2, PAR, HOM, KMULT

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_UP(UM,IPARTIC) = 0.0d0
      ENDDO

C  BOA source terms
C  ----------------

C  initialise boa source terms

      DO UM = 1, N_USER_STREAMS
        BOA_SOURCE(UM)        = 0.0d0
        DIRECT_BOA_SOURCE(UM) = 0.0d0
      ENDDO

C  Calculation

      IF ( DO_INCLUDE_SURFACE )THEN

C  Full solution: Downward intensity at computational angles (beam/homog)
C     --> Develop reflectance integrand  a(j).x(j).I(-j)

        N = NLAYERS
        PAR = WLOWER(1,N)
        HOM = LCON_XVEC(1,N)*T_DELT_EIGEN(N) + MCON_XVEC(1,N)
        IDOWNSURF = ( PAR + HOM ) * STREAM_VALUE

C  reflected multiple scatter intensity at user defined-angles

        IF ( SURFTYPE .EQ. 1 ) THEN
          KMULT = SURFACE_FACTOR * ALBEDO
          DO UM = 1, N_USER_STREAMS
            BOA_SOURCE(UM) = KMULT * IDOWNSURF
          ENDDO
        ENDIF

C  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM .AND. .NOT.DO_DB_CORRECTION ) THEN
          DO UM = 1, N_USER_STREAMS
            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
          ENDDO
        ENDIF

      ELSE
        DO UM = 1, N_USER_STREAMS
          BOA_SOURCE(UM)        = 0.0d0
          DIRECT_BOA_SOURCE(UM) = 0.0d0
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Set the cumulative source term equal to BOA values

      NC = 0
      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) +
     &                        DIRECT_BOA_SOURCE(UM)
      ENDDO

C  Recursion Loop in Source function integration
C  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N

        DO UM = 1, N_USER_STREAMS
          SHOM = LCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N) +
     &           MCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N)
          LAYERSOURCE(UM) = SHOM
        ENDDO

        DO UM = 1, N_USER_STREAMS
          SFOR2 =  EMULT_UP(UM,N,IPARTIC) * U_WPOS2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2 
       ENDDO

        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
          DO UM = 1, N_USER_STREAMS
            SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N,IPARTIC)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
          ENDDO
        ENDIF

        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_UP(UM,NC) = LAYERSOURCE(UM) +
     &           T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
        ENDDO

      ENDDO

C  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_UP(UM,IPARTIC) =
     &             FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_DNUSER_INTENSITY
     I      ( FOURIER_COMPONENT, DO_MSMODE_LIDORT,
     I        IPARTIC, FLUX_MULTIPLIER,
     I        NLAYERS, NBEAMS, N_USER_STREAMS,
     I        T_DELT_EIGEN, T_DELT_USERM, XPOS, XNEG,
     I        LCON, LCON_XVEC, MCON, MCON_XVEC,
     I        U_XPOS, U_XNEG, U_WNEG1, U_WNEG2,
     I        HMULT_1, HMULT_2, EMULT_DN,
     O        INTENSITY_F_DN,  CUMSOURCE_DN )

C  Subroutine input arguments
C  --------------------------

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  Control

      LOGICAL          DO_MSMODE_LIDORT

C  beam index

      INTEGER          IPARTIC

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN(NLAYERS)

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Solution constants of integration

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Solution constants of integration multiplied by homogeneous solutions

      DOUBLE PRECISION LCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION MCON_XVEC(2,NLAYERS)

C  Eigenvectors defined at user-defined stream angles
C     EP for the positive KEIGEN values, EM for -ve KEIGEN

      DOUBLE PRECISION
     U        U_XPOS(N_USER_STREAMS,NLAYERS),
     U        U_XNEG(N_USER_STREAMS,NLAYERS)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION
     U        U_WNEG1(N_USER_STREAMS,NLAYERS),
     U        U_WNEG2(N_USER_STREAMS,NLAYERS)

C  solution multipliers 

      DOUBLE PRECISION 
     &      HMULT_1(N_USER_STREAMS,NLAYERS),
     &      HMULT_2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  Outputs
C  -------

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_DN
     D   (N_USER_STREAMS,NBEAMS)

C  Cumulative source terms

      DOUBLE PRECISION
     U    CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

C  local variables
C  ---------------

C  Help variables

      INTEGER          UM, NC, N
      DOUBLE PRECISION LAYERSOURCE(N_USER_STREAMS)
      DOUBLE PRECISION SHOM, SFOR1, SFOR2

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_DN(UM,IPARTIC) = 0.0d0
      ENDDO

C  Initialize recursion for user-defined stream angles only

      NC = 0
      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_DN(UM,NC) = 0.0d0
      ENDDO

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term

      DO N = 1, NLAYERS
        NC = N

        DO UM = 1, N_USER_STREAMS
          SHOM = LCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N) +
     &           MCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N)
          LAYERSOURCE(UM) = SHOM
        ENDDO

        DO UM = 1, N_USER_STREAMS
          SFOR2 =  EMULT_DN(UM,N,IPARTIC) * U_WNEG2(UM,N)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2 
        ENDDO

        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
          DO UM = 1, N_USER_STREAMS
            SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N,IPARTIC)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
          ENDDO
        ENDIF

        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_DN(UM,NC) = LAYERSOURCE(UM) +
     &                    T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
        ENDDO

      ENDDO

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_DN(UM,IPARTIC) =
     &                 FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_CONVERGE
     I    ( DO_UPWELLING, DO_DNWELLING, 
     I      N_GEOMETRIES, NBEAMS, NTHREADS,
     I      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,
     I      THREAD, IBEAM, FOURIER_COMPONENT,
     I      INTENSITY_F_UP,  INTENSITY_F_DN, 
     O      INTENSITY_TOA, INTENSITY_BOA )

C  input variables
C  ---------------

C  Control

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  SS control, not required in this streamlined version
C      LOGICAL          DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

C  Numbers

      INTEGER          N_GEOMETRIES, NBEAMS, NTHREADS
      INTEGER          N_USER_STREAMS, N_USER_RELAZMS

C  FOurier component and thread, beam

      INTEGER          FOURIER_COMPONENT, THREAD, IBEAM

C  Local  azimuth factors

      INTEGER          UMOFF ( NBEAMS, N_USER_STREAMS )
      DOUBLE PRECISION AZMFAC
     &     (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_UP
     D   (N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION INTENSITY_F_DN
     D   (N_USER_STREAMS,NBEAMS)

C  Single scatter solutions, Not required here
C      DOUBLE PRECISION INTENSITY_SS_UP(N_GEOMETRIES)
C      DOUBLE PRECISION INTENSITY_SS_DN(N_GEOMETRIES)

C  Output
C  ------

      DOUBLE PRECISION INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      DOUBLE PRECISION INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

C  local variables
C  ---------------

      INTEGER          I, UA, V
      DOUBLE PRECISION TOLD, TAZM

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED do not get calculated)

c        IF ( .not. DO_SSFULL ) THEN

C  COde only for the NON-SSFULL case

          DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = INTENSITY_F_UP(I,IBEAM)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = INTENSITY_F_DN(I,IBEAM)
              ENDIF
            ENDDO
          ENDDO

C  Commented out in the streamlined version
c        IF ( DO_SSFULL ) THEN
c           DO I = 1, N_USER_STREAMS
c            DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) THEN
c                INTENSITY_TOA(V,THREAD) = 0.0d0
c              ENDIF
c              IF ( DO_DNWELLING ) THEN
c                INTENSITY_BOA(V,THREAD) = 0.0d0
c              ENDIF
c            ENDDO
c          ENDDO
c        ENDIF

C    Add the single scatter component if flagged
C  Commented out in the streamlined version
c        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
c           DO I = 1, N_USER_STREAMS
c            DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) THEN
c                INTENSITY_TOA(V,THREAD) = 
c     &             INTENSITY_TOA(V,THREAD) + INTENSITY_SS_UP(V)
c              ENDIF
c              IF ( DO_DNWELLING ) THEN
c                INTENSITY_BOA(V,THREAD) = 
c     &             INTENSITY_BOA(V,THREAD) + INTENSITY_SS_DN(V)
c              ENDIF
c            ENDDO
c          ENDDO
c       ENDIF

C  ######################
C  Fourier component = 1
C  ######################

      ELSE

C  No examination of convergence
C  -----------------------------

        DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
            V = UMOFF(IBEAM,I) + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_UP(I,IBEAM)
              INTENSITY_TOA(V,THREAD) = TOLD + TAZM
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_DN(I,IBEAM)
              INTENSITY_BOA(V,THREAD) = TOLD + TAZM
            ENDIF
          ENDDO
        ENDDO

C  Finish Fourier
   
      ENDIF

C  Finish

      RETURN
      END
