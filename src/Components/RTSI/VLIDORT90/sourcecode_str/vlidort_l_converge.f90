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
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################


      MODULE vlidort_l_converge_module

      PRIVATE
      PUBLIC :: VLIDORT_L_CONVERGE

      CONTAINS

      SUBROUTINE VLIDORT_L_CONVERGE ( &
        FOURIER_COMPONENT, LOCAL_N_USERAZM, &
        IBEAM, AZMFAC, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SS_EXTERNAL, & ! New 15 March 2012
        DO_SSFULL, DO_UPWELLING, &
        NSTOKES, NLAYERS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, &
        DO_DBCORRECTION, DO_NO_AZIMUTH, &
        N_USER_STREAMS, LOCAL_UM_START, &
        N_DIRECTIONS, WHICH_DIRECTIONS, &
        N_OUT_STREAMS, &
        VZA_OFFSETS, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, DO_SURFACE_LINEARIZATION, &
        N_SURFACE_WFS, &
        PROFILEWF_SS, COLUMNWF_SS, &
        PROFILEWF_DB, COLUMNWF_DB, SURFACEWF_DB, &
        SURFACEWF_F, ATMOSWF_F, &
        SURFACEWF, PROFILEWF, COLUMNWF )

!  Just upgrades the weighting function Fourier cosine series

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::          LOCAL_N_USERAZM
      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: AZMFAC &
          ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL, INTENT (IN) ::          DO_SS_EXTERNAL
      LOGICAL, INTENT (IN) ::          DO_SSFULL
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          N_OUT_STREAMS
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: ATMOSWF_F &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables

      INTEGER :: I, IDIR, UT, UA, Q, N, K, W, Z, V, O1, LVARY

!   For Single scatter corrections, the quadrature output flag is
!    turned off, so that N_OUT_STREAMS = N_USER_STREAMS, and the
!    offsetting is consistent here with the usage elsewhere.

      K = 0

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version

!  Profile atmospheric weighting functions
!  ---------------------------------------

       IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Diffuse field at all output angles

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
                 PROFILEWF(Q,N,UT,V,O1,W) &
                           = ATMOSWF_F(Q,N,UT,I,IBEAM,O1,W)
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

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
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
                 PROFILEWF(Q,N,UT,V,O1,W) = &
                 PROFILEWF(Q,N,UT,V,O1,W) + PROFILEWF_SS(Q,N,UT,V,O1,W)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO UT = 1, N_USER_LEVELS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 PROFILEWF(Q,N,UT,V,O1,UPIDX) = &
                PROFILEWF(Q,N,UT,V,O1,UPIDX) + PROFILEWF_DB(Q,N,UT,V,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  end atmospheric WF clause

       ENDIF

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

       IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Diffuse field at all output angles

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
                 COLUMNWF(Q,UT,V,O1,W) = &
                          ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W)
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

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O1 = 1, NSTOKES
                  COLUMNWF(Q,UT,V,O1,W) = &
                    COLUMNWF(Q,UT,V,O1,W) + COLUMNWF_SS(Q,UT,V,O1,W)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

!   (New code, 6 May 2005)
!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 3.3

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UT = 1, N_USER_LEVELS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,UPIDX) = &
                 COLUMNWF(Q,UT,V,O1,UPIDX) + COLUMNWF_DB(Q,UT,V,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end bulk/column atmospheric WF clause

       ENDIF

!  Reflectance weighting functions
!  -------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)

       IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( .not. DO_SSFULL ) THEN
         DO Z = 1, N_SURFACE_WFS
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
         DO Z = 1, N_SURFACE_WFS
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

!   (New code, 6 May 2005)
!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
         DO Z = 1, N_SURFACE_WFS
          DO UT = 1, N_USER_LEVELS
           DO I = LOCAL_UM_START, N_OUT_STREAMS
            DO UA = 1, LOCAL_N_USERAZM
             V = VZA_OFFSETS(IBEAM,I) + UA
             DO O1 = 1, NSTOKES
              SURFACEWF(Z,UT,V,O1,UPIDX) = &
                SURFACEWF(Z,UT,V,O1,UPIDX) + SURFACEWF_DB(Z,UT,V,O1)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

!  End albedo WF clause

       ENDIF

!  If no_azimuth, then exit
!  ------------------------

       IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

       IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

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
                PROFILEWF(Q,N,UT,V,O1,W) = PROFILEWF(Q,N,UT,V,O1,W) + &
                ATMOSWF_F(Q,N,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  End profile atmospheric WF clause

       ENDIF

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

       IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Variational index

         LVARY = 0

!  Add next Fourier component to output

         DO Q = 1, N_TOTALCOLUMN_WFS
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
               V = VZA_OFFSETS(IBEAM,I) + UA
               DO O1 = 1, NSTOKES
                 COLUMNWF(Q,UT,V,O1,W) = COLUMNWF(Q,UT,V,O1,W) + &
                 ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
!       if(q.eq.3.and.ut.ne.5)write(*,*)
!     &         fourier_component,ut,v,ATMOSWF_F(Q,LVARY,UT,I,IBEAM,O1,W
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDDO

!  End bulk/column atmospheric WF clause

       ENDIF

!  albedo weighting functions (non-Lambertian only)
!  ------------------------------------------------

       IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

!  Copy full output

         DO Z = 1, N_SURFACE_WFS
          DO UA = 1, LOCAL_N_USERAZM
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
             DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO O1 = 1, NSTOKES
               SURFACEWF(Z,UT,V,O1,W) = SURFACEWF(Z,UT,V,O1,W) + &
                SURFACEWF_F(Z,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

!  End albedo WF clause

        ENDIF
       ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_CONVERGE

      END MODULE vlidort_l_converge_module

