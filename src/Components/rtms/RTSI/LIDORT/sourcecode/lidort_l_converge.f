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
C #            LIDORT_L_CONVERGE (master)                       #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        IBEAM )

C  Just upgrades the weighting function Fourier cosine series

C  Include files
C  -------------

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include files of already calculated results

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  include file with geometrical indexing offsets (input)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Include file of final weighting function results (OUTPUT)

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER          FOURIER_COMPONENT, LOCAL_N_USERAZM      
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      INTEGER          IBEAM

C  local variables

      INTEGER          I, IDIR, UT, UA, Q, N, K, W, Z, V, LVARY

C   For Single scatter corrections, the quadrature output flag is
C    turned off, so that N_OUT_STREAMS = N_USER_STREAMS, and the
C    offsetting is consistent here with the usage elsewhere.

      K = 0

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED do not get calculated)

C  Full single scatter calculation is initialized to zero here (Version 2.3)

C  Profile atmospheric weighting functions
C  ---------------------------------------

       IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Diffuse field at all output angles

        IF ( .not. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = 1, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                PROFILEWF(Q,N,UT,V,W) = ATMOSWF_F(Q,N,UT,I,IBEAM,W)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = 1, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                PROFILEWF(Q,N,UT,V,W) = ZERO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C    Add the single scatter component if flagged
C     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
C     Version 3.2.   Added outgoing correction flag to this.....
C     Version 3.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                PROFILEWF(Q,N,UT,V,W) =
     &           PROFILEWF(Q,N,UT,V,W) + PROFILEWF_SS(Q,N,UT,V,W)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C   (New code, 6 May 2005)
C    Add the  component if flagged (upwelling only)
C       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
C    Full single scatter option added Version 3.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO UT = 1, N_OUT_USERTAUS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = UMOFF(IBEAM,I) + UA
               PROFILEWF(Q,N,UT,V,UPIDX) =
     &           PROFILEWF(Q,N,UT,V,UPIDX) + PROFILEWF_DB(Q,N,UT,V)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C------------------------------------------------- SECTION REMOVED 04/02/07
C  M-SCAT source term weighting functions
C  --------------------------------------
c        IF ( SAVE_LAYER_MSST ) THEN
C  downwelling output
C     Layer terms only, copy Fourier component 0
c         IF ( DO_DNWELLING ) THEN
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO K = 1, N_LAYERSOURCE_DN
c              DO I = 1, N_USER_STREAMS
c               DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_MSCATSTERM(Q,N,K,V,DNIDX) =
c     $               ATMOSWF_MSCATSTERM_F(Q,N,K,I,IBEAM,DNIDX)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
c         ENDIF
C  Upwelling output
c         IF ( DO_UPWELLING ) THEN
C  Direct beam output at BOA
C    Use exact value if flagged, otherwise copy Fourier m = 0
c          IF ( DO_DBCORRECTION ) THEN
c           DO N = 1, NLAYERS
c            IF ( LAYER_VARY_FLAG(N) ) THEN
c             DO Q = 1, LAYER_VARY_NUMBER(N)
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c               DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_DIRECTBOA_STERM(Q,N,V) =
c     &                  L_EXACTDB_SOURCE(Q,N,V)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDIF
c           ENDDO
c          ELSE
c           DO N = 1, NLAYERS
c            IF ( LAYER_VARY_FLAG(N) ) THEN
c             DO Q = 1, LAYER_VARY_NUMBER(N)
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c               DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_DIRECTBOA_STERM(Q,N,V) =
c     $                  ATMOSWF_DIRECTBOA_STERM_F(Q,N,I,IBEAM)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDIF
c          ENDDO
c          ENDIF
C  Diffuse source term at BOA (copy Fourier m = 0 )
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO I = LOCAL_UM_START, N_USER_STREAMS
c              DO UA = 1, LOCAL_N_USERAZM
c               V = UMOFF(IBEAM,I) + UA
c               ATMOSWF_MSCATBOA_STERM(Q,N,V) =
c     $                  ATMOSWF_MSCATBOA_STERM_F(Q,N,I,IBEAM)
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
C  Layer terms, copy Fourier 0 component
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO K = N_LAYERSOURCE_UP, NLAYERS
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c               DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_MSCATSTERM(Q,N,K,V,UPIDX) =
c     $             ATMOSWF_MSCATSTERM_F(Q,N,K,I,IBEAM,UPIDX)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
C  end do upwelling
c         ENDIF
c  End save layering clause
c        ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

C  end profile atmospheric WF clause

       ENDIF

C  Bulk/column atmospheric weighting functions (Version 3.3)
C  ---------------------------------------------------------

C  This section newly written, 26 September 2006, installed 14 May 2007.

       IF ( DO_COLUMN_LINEARIZATION ) THEN

C  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
          LVARY = 0
          DO Q = 1, N_TOTALCOLUMN_WFS
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_OUT_USERTAUS
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = UMOFF(IBEAM,I) + UA
               COLUMNWF(Q,UT,V,W) = ATMOSWF_F(Q,LVARY,UT,I,IBEAM,W)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_TOTALCOLUMN_WFS
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_OUT_USERTAUS
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = UMOFF(IBEAM,I) + UA
               COLUMNWF(Q,UT,V,W) = ZERO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C    Add the single scatter component if flagged
C     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
C     Version 3.2.   Added outgoing correction flag to this.....
C     Version 3.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                COLUMNWF(Q,UT,V,W) =
     &              COLUMNWF(Q,UT,V,W) + COLUMNWF_SS(Q,UT,V,W)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

C   (New code, 6 May 2005)
C    Add the  component if flagged (upwelling only)
C       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
C    Full single scatter option added Version 3.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UT = 1, N_OUT_USERTAUS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = UMOFF(IBEAM,I) + UA
               COLUMNWF(Q,UT,V,UPIDX) =
     &           COLUMNWF(Q,UT,V,UPIDX) + COLUMNWF_DB(Q,UT,V)
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

C  end bulk/column atmospheric WF clause

       ENDIF

C  Reflectance weighting functions
C  -------------------------------

       IF ( DO_SURFACE_LINEARIZATION ) THEN

C  DIffuse field at all output angles
C  Alternative - zero the output (single scatter case)

        IF ( .not. DO_SSFULL ) THEN
         DO Z = 1, N_TOTALBRDF_WFS
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_OUT_USERTAUS
            DO I = LOCAL_UM_START, N_USER_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = UMOFF(IBEAM,I) + UA
              SURFACEWF(Z,UT,V,W) = SURFACEWF_F(Z,UT,I,IBEAM,W)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ELSE
         DO Z = 1, N_TOTALBRDF_WFS
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_OUT_USERTAUS
            DO I = LOCAL_UM_START, N_USER_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = UMOFF(IBEAM,I) + UA
              SURFACEWF(Z,UT,V,W) = ZERO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

C   (New code, 6 May 2005)
C    Add the  component if flagged (upwelling only)
C       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
C    Full single scatter option added Version 2.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
         DO Z = 1, N_TOTALBRDF_WFS
          DO UT = 1, N_OUT_USERTAUS
           DO I = LOCAL_UM_START, N_OUT_STREAMS
            DO UA = 1, LOCAL_N_USERAZM
             V = UMOFF(IBEAM,I) + UA
             SURFACEWF(Z,UT,V,UPIDX) =
     &          SURFACEWF(Z,UT,V,UPIDX) + SURFACEWF_DB(Z,UT,V)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

C------------------------------------------------- SECTION REMOVED 04/02/07
C  M-SCAT source term weighting functions
C  --------------------------------------
c        IF ( SAVE_LAYER_MSST ) THEN
C  downwelling output
C     Layer terms only, copy Fourier component 0
c         IF ( DO_DNWELLING ) THEN
c          DO Z = 1, N_TOTALBRDF_WFS
c           DO K = 1, N_LAYERSOURCE_DN
c            DO I = 1, N_USER_STREAMS
c             DO UA = 1, LOCAL_N_USERAZM
c              V = UMOFF(IBEAM,I) + UA
c              SURFACEWF_MSCATSTERM(Z,K,V,DNIDX) =
c     $               SURFACEWF_MSCATSTERM_F(Z,K,I,IBEAM,DNIDX)
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDDO
c         ENDIF
C  Upwelling output
c         IF ( DO_UPWELLING ) THEN
C  Direct beam output at BOA
C    Use exact value if flagged, otherwise copy Fourier m = 0
c          IF ( DO_DBCORRECTION ) THEN
c           DO Z = 1, N_TOTALBRDF_WFS
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c             DO UA = 1, LOCAL_N_USERAZM
c              V = UMOFF(IBEAM,I) + UA
c              SURFACEWF_DIRECTBOA_STERM(Z,V) =
c     &                          LS_EXACTDB_SOURCE(Z,V)
c             ENDDO
c            ENDDO
c           ENDDO
c          ELSE
c           DO Z = 1, N_TOTALBRDF_WFS
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c             DO UA = 1, LOCAL_N_USERAZM
c              V = UMOFF(IBEAM,I) + UA
c              SURFACEWF_DIRECTBOA_STERM(Z,V) =
c     $              SURFACEWF_DIRECTBOA_STERM_F(Z,I,IBEAM)
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDIF
C  Diffuse source term at BOA (copy Fourier m = 0 )
c          DO Z = 1, N_TOTALBRDF_WFS
c           DO I = LOCAL_UM_START, N_USER_STREAMS
c            DO UA = 1, LOCAL_N_USERAZM
c             V = UMOFF(IBEAM,I) + UA
c             SURFACEWF_MSCATBOA_STERM(Z,V) =
c     $               SURFACEWF_MSCATBOA_STERM_F(Z,I,IBEAM)
c            ENDDO
c           ENDDO
c          ENDDO
C  Layer terms, copy Fourier 0 component
c          DO Z = 1, N_TOTALBRDF_WFS
c           DO K = N_LAYERSOURCE_UP, NLAYERS
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c             DO UA = 1, LOCAL_N_USERAZM
c              V = UMOFF(IBEAM,I) + UA
c              SURFACEWF_MSCATSTERM(Z,K,V,UPIDX) =
c     $             SURFACEWF_MSCATSTERM_F(Z,K,I,IBEAM,UPIDX)
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDDO
C  end do upwelling
c         ENDIF
C  End save layering clause
c        ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

C  End albedo WF clause

       ENDIF

C  If no_azimuth, then exit
C  ------------------------

       IF ( DO_NO_AZIMUTH ) RETURN

C  ######################
C  Fourier components > 0
C  ######################

      ELSE

C  Profile atmospheric weighting functions
C  ---------------------------------------

       IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = 1, N_OUT_STREAMS
               V = UMOFF(IBEAM,I) + UA
               PROFILEWF(Q,N,UT,V,W) = PROFILEWF(Q,N,UT,V,W) + 
     $            ATMOSWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

C------------------------------------------------- SECTION REMOVED 04/02/07
C  M_SCAT source term atmospheric weighting functions
c        IF ( SAVE_LAYER_MSST ) THEN
C  Downwelling layer source terms only
c         IF ( DO_DNWELLING ) THEN
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO UA = 1, LOCAL_N_USERAZM
c              DO UT = 1, N_OUT_USERTAUS
c               DO I = 1, N_OUT_STREAMS
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_MSCATSTERM(Q,N,K,V,DNIDX) =
c     $             ATMOSWF_MSCATSTERM(Q,N,K,V,DNIDX) + 
c     $            ATMOSWF_MSCATSTERM_F(Q,N,K,I,IBEAM,DNIDX)*AZMFAC(UA)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
c         ENDIF
C  upwelling
c         IF ( DO_UPWELLING ) THEN
c          W = UPIDX
C  ..... Upwelling, do Diffuse BOA source terms
C        (contributions are zero if no surface reflection)
C        ( only do direct beam if no correction)
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO UA = 1, LOCAL_N_USERAZM
c              DO K = N_LAYERSOURCE_UP, NLAYERS
c               DO I = 1, N_USER_STREAMS
c                V = UMOFF(IBEAM,I) + UA
c                ATMOSWF_MSCATBOA_STERM(Q,N,V) =
c     $            ATMOSWF_MSCATBOA_STERM(Q,N,V) + 
c     $              ATMOSWF_MSCATBOA_STERM_F(Q,N,I,IBEAM)*AZMFAC(UA)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
C  ..... Upwelling, do Diffuse BOA source terms
C        (contributions are zero if no surface reflection)
C        ( only do direct beam if no correction)
c          IF ( .NOT. DO_DBCORRECTION ) THEN
c           DO N = 1, NLAYERS
c            IF ( LAYER_VARY_FLAG(N) ) THEN
c             DO Q = 1, LAYER_VARY_NUMBER(N)
c              DO UA = 1, LOCAL_N_USERAZM
c               DO K = N_LAYERSOURCE_UP, NLAYERS
c                DO I = 1, N_USER_STREAMS
c                 V = UMOFF(IBEAM,I) + UA
c                 ATMOSWF_DIRECTBOA_STERM(Q,N,V) =
c     $            ATMOSWF_DIRECTBOA_STERM(Q,N,V) + 
c     $                ATMOSWF_DIRECTBOA_STERM_F(Q,N,I,IBEAM)*AZMFAC(UA)
c                ENDDO
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDIF
c           ENDDO
c          ENDIF
C  Do layer source terms
c          DO N = 1, NLAYERS
c           IF ( LAYER_VARY_FLAG(N) ) THEN
c            DO Q = 1, LAYER_VARY_NUMBER(N)
c             DO UA = 1, LOCAL_N_USERAZM
c              DO K = N_LAYERSOURCE_UP, NLAYERS
c               DO I = 1, N_USER_STREAMS
c                V = UMOFF(IBEAM,I) + UA
c                 ATMOSWF_MSCATSTERM(Q,N,K,V,UPIDX) =
c     $              ATMOSWF_MSCATSTERM(Q,N,K,V,UPIDX) + 
c     $            ATMOSWF_MSCATSTERM_F(Q,N,K,I,IBEAM,UPIDX)*AZMFAC(UA)
c               ENDDO
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
C  End upwelling
c         ENDIF
C  End MSST clause
c        ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

C  End profile atmospheric WF clause

       ENDIF

C  Bulk/column atmospheric weighting functions
C  -------------------------------------------

       IF ( DO_COLUMN_LINEARIZATION ) THEN

C  Variational index

         LVARY = 0

C  Add next Fourier component to output

         DO Q = 1, N_TOTALCOLUMN_WFS
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_OUT_USERTAUS
              DO I = 1, N_OUT_STREAMS
               V = UMOFF(IBEAM,I) + UA
               COLUMNWF(Q,UT,V,W) = COLUMNWF(Q,UT,V,W) + 
     $            ATMOSWF_F(Q,LVARY,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDDO

C  End bulk/column atmospheric WF clause

       ENDIF

C  albedo weighting functions (non-Lambertian only)
C  ------------------------------------------------

       IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

C  Copy full output

         DO Z = 1, N_TOTALBRDF_WFS
          DO UA = 1, LOCAL_N_USERAZM
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_OUT_USERTAUS
             DO I = 1, N_OUT_STREAMS
              V = UMOFF(IBEAM,I) + UA
              SURFACEWF(Z,UT,V,W) = SURFACEWF(Z,UT,V,W) +
     &          SURFACEWF_F(Z,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

C------------------------------------------------- SECTION REMOVED 04/02/07
C  M-SCAT source term atmospheric weighting functions
c         IF ( SAVE_LAYER_MSST ) THEN
C  Downwelling layer source terms only
c          IF ( DO_DNWELLING ) THEN
c           DO Z = 1, N_TOTALBRDF_WFS
c            DO UA = 1, LOCAL_N_USERAZM
c             DO K = 1, N_LAYERSOURCE_DN
c              DO I = 1, N_USER_STREAMS
c               V = UMOFF(IBEAM,I) + UA
c               SURFACEWF_MSCATSTERM(Z,K,V,DNIDX) =
c     $            SURFACEWF_MSCATSTERM(Z,K,V,DNIDX) + 
c     $            SURFACEWF_MSCATSTERM_F(Z,K,I,IBEAM,DNIDX)*AZMFAC(UA)
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDIF
C  Upwelling
c          IF ( DO_UPWELLING ) THEN
C  Direct beam output at BOA (only if no DB correction applied)
c           IF ( .NOT.DO_DBCORRECTION ) THEN
c            DO Z = 1, N_TOTALBRDF_WFS
c             DO I = LOCAL_UM_START, N_USER_STREAMS
c              DO UA = 1, LOCAL_N_USERAZM
c               V = UMOFF(IBEAM,I) + UA
c               SURFACEWF_DIRECTBOA_STERM(Z,V) =
c     $            SURFACEWF_DIRECTBOA_STERM(Z,V) +
c     $            SURFACEWF_DIRECTBOA_STERM_F(Z,I,IBEAM)*AZMFAC(UA)
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDIF
C  Diffuse source term at BOA
c           DO Z = 1, N_TOTALBRDF_WFS
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c             DO UA = 1, LOCAL_N_USERAZM
c              V = UMOFF(IBEAM,I) + UA
c              SURFACEWF_MSCATBOA_STERM(Z,V) =
c     $          SURFACEWF_MSCATBOA_STERM(Z,V) +
c     $          SURFACEWF_MSCATBOA_STERM_F(Z,I,IBEAM)*AZMFAC(UA)
c             ENDDO
c            ENDDO
c           ENDDO
C  Layer terms
c           DO Z = 1, N_TOTALBRDF_WFS
c            DO K = N_LAYERSOURCE_UP, NLAYERS
c             DO I = LOCAL_UM_START, N_USER_STREAMS
c              DO UA = 1, LOCAL_N_USERAZM
c               V = UMOFF(IBEAM,I) + UA
c               SURFACEWF_MSCATSTERM(Z,K,V,UPIDX) =
c     $            SURFACEWF_MSCATSTERM(Z,K,V,UPIDX) +
c     $            SURFACEWF_MSCATSTERM_F(Z,K,I,IBEAM,UPIDX)
c              ENDDO
c             ENDDO
c            ENDDO
c           ENDDO
C  end do upwelling
c          ENDIF
C  End save layering clause
c         ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

C  End albedo WF clauses

        ENDIF
       ENDIF

C  end Fourier clause

      ENDIF

C  Finish

      RETURN
      END
