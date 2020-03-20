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
   
      SUBROUTINE VLIDORT_L_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        IBEAM )

C  Just upgrades the weighting function Fourier cosine series

C  Include files
C  -------------

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Include files of already calculated results

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  Include file of final weighting function results

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER          FOURIER_COMPONENT, LOCAL_N_USERAZM      
      INTEGER          IBEAM
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS,MAXSTOKES)

C  local variables

      INTEGER          I, IDIR, UT, UA, Q, N, K, W, Z, V, O1, LVARY

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
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O1 = 1, NSTOKES
                 PROFILEWF(Q,N,UT,V,O1,W)
     &                     = ATMOSWF_F(Q,N,UT,I,IBEAM,O1,W)
                ENDDO
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
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O1 = 1, NSTOKES
                 PROFILEWF(Q,N,UT,V,O1,W) = ZERO
                ENDDO
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
C     Version 2.1.   Added outgoing correction flag to this.....
C     Version 2.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O1 = 1, NSTOKES
                 PROFILEWF(Q,N,UT,V,O1,W) =
     &           PROFILEWF(Q,N,UT,V,O1,W) + PROFILEWF_SS(Q,N,UT,V,O1,W)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C    Add the  component if flagged (upwelling only)
C       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
C    Full single scatter option added Version 2.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO UT = 1, N_USER_LEVELS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 PROFILEWF(Q,N,UT,V,O1,UPIDX) =
     &          PROFILEWF(Q,N,UT,V,O1,UPIDX) + PROFILEWF_DB(Q,N,UT,V,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C  end atmospheric WF clause

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
            DO UT = 1, N_USER_LEVELS
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,W) =
     &                    ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_TOTALCOLUMN_WFS
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,W) = ZERO
               ENDDO
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
             DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O1 = 1, NSTOKES
                  COLUMNWF(Q,UT,V,O1,W) =
     &              COLUMNWF(Q,UT,V,O1,W) + COLUMNWF_SS(Q,UT,V,O1,W)
                ENDDO
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
            DO UT = 1, N_USER_LEVELS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,UPIDX) =
     &           COLUMNWF(Q,UT,V,O1,UPIDX) + COLUMNWF_DB(Q,UT,V,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

C  end bulk/column atmospheric WF clause

       ENDIF

C  Reflectance weighting functions
C  -------------------------------

C  DIffuse field at all output angles
C  Alternative - zero the output (single scatter case)

       IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( .not. DO_SSFULL ) THEN
         DO Z = 1, N_TOTALBRDF_WFS
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_USER_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO O1 = 1, NSTOKES
               SURFACEWF(Z,UT,V,O1,W) = SURFACEWF_F(Z,UT,I,IBEAM,O1,W)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ELSE
         DO Z = 1, N_TOTALBRDF_WFS
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_USER_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO O1 = 1, NSTOKES
               SURFACEWF(Z,UT,V,O1,W) = ZERO
              ENDDO
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
          DO UT = 1, N_USER_LEVELS
           DO I = LOCAL_UM_START, N_OUT_STREAMS
            DO UA = 1, LOCAL_N_USERAZM
             V = VZA_OFFSETS(IBEAM,I) + UA
             DO O1 = 1, NSTOKES
              SURFACEWF(Z,UT,V,O1,UPIDX) =
     &          SURFACEWF(Z,UT,V,O1,UPIDX) + SURFACEWF_DB(Z,UT,V,O1)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

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
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                PROFILEWF(Q,N,UT,V,O1,W) = PROFILEWF(Q,N,UT,V,O1,W) + 
     $          ATMOSWF_F(Q,N,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

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
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,W) = COLUMNWF(Q,UT,V,O1,W) + 
     $        ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
!       if(q.eq.3.and.ut.ne.5)write(*,*)
!     &         fourier_component,ut,v,ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W)
               ENDDO
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
            DO UT = 1, N_USER_LEVELS
             DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO O1 = 1, NSTOKES
               SURFACEWF(Z,UT,V,O1,W) = SURFACEWF(Z,UT,V,O1,W) +
     &          SURFACEWF_F(Z,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

C  End albedo WF clause

        ENDIF
       ENDIF

C  end Fourier clause

      ENDIF

C  Finish

      RETURN
      END
