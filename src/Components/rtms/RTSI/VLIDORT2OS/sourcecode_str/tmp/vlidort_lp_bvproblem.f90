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
! #  For linearizations involving atmospheric parameters        #
! #                                                             #
! #            LP_BVP_SOLUTION_MASTER                           #
! #            LP_BVP_COLUMN_SETUP (re-named, Version 2.4)      #
! #            LP_BEAMSOLUTION_NEQK (re-named)                  #
! #            LP_BEAMSOLUTION_NNEK (re-named)                  #
! #                                                             #
! #  For linearizations with telescoped boundary value problem  #
! #                                                             #
! #            LP_BVPTEL_SOLUTION_MASTER                        #
! #            LP_BVPTEL_COLUMN_SETUP                           #
! #                                                             #
! ###############################################################


      MODULE vlidort_lp_bvproblem

      PRIVATE
      PUBLIC :: LP_BVP_SOLUTION_MASTER, &
                LP_BVPTEL_SOLUTION_MASTER

      CONTAINS

      SUBROUTINE LP_BVP_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_THERMEMISS, &
        VARIATION_INDEX, N_WEIGHTFUNCS, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, NLAYERS, NSTREAMS_2, &
        NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DELTAU_SLANT, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, DIRECT_BEAM, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_T_WUPPER, L_T_WLOWER, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
        LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_BVEC, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
        LAYER_TO_VARY, N_LAYER_WFS, &
        DO_PLANE_PARALLEL, &
        N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        L_WLOWER, L_WUPPER, NCON, PCON, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS
      USE VLIDORT_LPC_BVPROBLEM

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          VARIATION_INDEX
      INTEGER, INTENT (IN) ::          N_WEIGHTFUNCS
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION, INTENT (INOUT) ::  NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  boundary condition flags

      LOGICAL :: MODIFIED_BCL3, MODIFIED_BCL4

!  error tracing variables

      INTEGER :: STATUS_SUB

!  helper variables

      DOUBLE PRECISION :: COL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for profile linearizations
!  ----------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
      MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

      CALL LP_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_THERMEMISS, &
        MODIFIED_BCL3, MODIFIED_BCL4, &
        VARIATION_INDEX, N_WEIGHTFUNCS, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DELTAU_SLANT, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, DIRECT_BEAM, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_T_WUPPER, L_T_WLOWER, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        DO_PLANE_PARALLEL, LP_INITIAL_TRANS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )

!  Back-substitution

      CALL L_BVP_BACKSUB ( &
        VARIATION_INDEX, N_WEIGHTFUNCS, &
        NLAYERS, NTOTAL, N_SUBDIAG, &
        N_SUPDIAG, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        LCON, MCON, &
        COL2_WF, SCOL2_WF, &
        NCON, PCON, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_BVP_SOLUTION_MASTER

!

      SUBROUTINE LP_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_THERMEMISS, &
        MODIFIED_BCL3, MODIFIED_BCL4, &
        LAYER_TO_VARY, N_LAYER_WFS, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DELTAU_SLANT, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, DIRECT_BEAM, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_T_WUPPER, L_T_WLOWER, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        DO_PLANE_PARALLEL, LP_INITIAL_TRANS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )

      USE VLIDORT_PARS
      USE VLIDORT_LPC_BVPROBLEM

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          MODIFIED_BCL3, MODIFIED_BCL4
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) :: COL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          Q, N, N1, I, I1, IR, IROW, CM, C0, O1, NV
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION :: T1, T2, T1R, T1I, T2R, T2I

      DOUBLE PRECISION :: R2_L_BEAM &
          ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      LOGICAL ::          REGULAR_BCL3, REGULAR_BCL4

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  Layer to vary

      NV = LAYER_TO_VARY

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS
            DO N = 1, NLAYERS
              L_WUPPER(I,O1,N,Q) = ZERO
              L_WLOWER(I,O1,N,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,O1,NV,Q) = L_T_WUPPER(I,NV,Q)
            L_WLOWER(I,O1,NV,Q) = L_T_WLOWER(I,NV,Q)
          ENDDO
        ENDDO
      ENDIF

!  Get the linearized beam solution for the varying layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LP_BEAMSOLUTION_NEQK ( &
          FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS, &
          NSTOKES, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          NSTREAMS_2, DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, &
          BVEC, &
          LP_T_DELT_MUBAR, LP_BVEC, &
          L_WUPPER, L_WLOWER )
      ENDIF

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT. MODIFIED_BCL3
      REGULAR_BCL4 = .NOT. MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

!  Start  loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS

!  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  No variation case (BCL1)

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS
            COL2_WF(IROW,Q) = ZERO
          ENDDO
         ENDDO
        ENDDO

      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------------

!  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN

!  .. nothing varying in these layers

        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0 = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_LAYER_WFS
                COL2_WF(CM,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

!  offsets

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contribution

            L_BEAM  = + L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  End BCL3 condition

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

!  offsets

        N = LAYER_TO_VARY
        N1  = N + 1
        C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LP_BEAMSOLUTION_NNEK ( &
            FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS, &
            DO_PLANE_PARALLEL, NSTOKES, &
            DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
            NSTREAMS_2, DO_LAYER_SCATTERING, &
            T_DELT_MUBAR, INITIAL_TRANS, &
            BVEC, LP_INITIAL_TRANS, &
            LP_T_DELT_MUBAR, LP_BVEC, &
            L_WUPPER, L_WLOWER )
        ENDIF

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  End BCL4 condition

      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final
!  ---------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

!  offsets

          N1  = N + 1
          C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL LP_BEAMSOLUTION_NNEK ( &
              FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS, &
              DO_PLANE_PARALLEL, NSTOKES, &
              DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
              NSTREAMS_2, DO_LAYER_SCATTERING, &
              T_DELT_MUBAR, INITIAL_TRANS, &
              BVEC, LP_INITIAL_TRANS, &
              LP_T_DELT_MUBAR, LP_BVEC, &
              L_WUPPER, L_WLOWER )
          ENDIF

!  .. contributions from beam solution (direct assign). No homog. variat

          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)
                COL2_WF(CM,Q) = L_BEAM
              ENDDO
            ENDDO
          ENDDO

!  end layer loop

        ENDDO

!  end BCL5 boundary conditions

      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------------

      N = NLAYERS

!  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN

!  get the linearized downward-reflected term

        CALL L_BVP_SURFACE_SETUP ( &
          DO_INCLUDE_SURFACE, MODIFIED_BCL4, &
          IBEAM, FOURIER_COMPONENT, &
          SURFACE_FACTOR, N_LAYER_WFS, &
          NSTOKES, NSTREAMS, &
          NLAYERS, DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, BRDF_F, &
          QUAD_STRMWTS, NSTKS_NSTRMS, &
          MUELLER_INDEX, T_DELT_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, L_T_DELT_EIGEN, &
          L_SOLA_XPOS, L_SOLB_XNEG, &
          L_WLOWER, &
          R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  ordinary BCL6 Component loop

      ELSE

!  get the linearized downward-reflected term

        CALL L_BVP_SURFACE_SETUP ( &
          DO_INCLUDE_SURFACE, MODIFIED_BCL4, &
          IBEAM, FOURIER_COMPONENT, &
          SURFACE_FACTOR, N_LAYER_WFS, &
          NSTOKES, NSTREAMS, &
          NLAYERS, DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, BRDF_F, &
          QUAD_STRMWTS, NSTKS_NSTRMS, &
          MUELLER_INDEX, T_DELT_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, L_T_DELT_EIGEN, &
          L_SOLA_XPOS, L_SOLB_XNEG, &
          L_WLOWER, &
          R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops : beam contributions only

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS
              L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
              COL2_WF(CM,Q) = - L_BEAM
            ENDDO
          ENDDO
        ENDDO

!  End lowest layer boundary value linearization

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * &
                      DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM)
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  debug

!      if ( do_debug_write ) then
!        DO N = 1, NTOTAL
!          write(95,'(4i4,1p4e17.9)')
!     &      FOURIER_COMPONENT,IBEAM,LAYER_TO_VARY,N,
!     &                 COL2_WF(N,1)
!        ENDDO
!      ENDIF
!       pause

!  finish

      RETURN
      END SUBROUTINE LP_BVP_COLUMN_SETUP

!

      SUBROUTINE LP_BEAMSOLUTION_NEQK ( &
        FOURIER, IB, N, N_LAYER_WFS, &
        NSTOKES, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        NSTREAMS_2, DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, &
        BVEC, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        L_WUPPER, L_WLOWER )

!  Linearization of beam particular integral in one layer.

!  In this module, this is the Layer that contains the variation.
!                   N = LAYER_TO_VARY

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IB
      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER          :: I, O1, Q
      DOUBLE PRECISION :: CONST, WDEL, VAR1
      DOUBLE PRECISION :: LBSOL(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS,2)

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive
!    Exit (Solutions have already been zeroed).

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) &
             .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        RETURN
      ENDIF

!  Classical solution
!  ==================

!  Very simple, same code for all situations

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        DO Q = 1, N_LAYER_WFS
          VAR1 = LP_T_DELT_MUBAR(N,N,IB,Q) * CONST
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
              LBSOL(I,O1,Q,1) = CONST * LP_BVEC(I,O1,N,N,Q)
              LBSOL(I,O1,Q,2) = WDEL  * LBSOL(I,O1,Q,1) &
                                    + VAR1 * BVEC(I,O1,N)
            ENDDO
          ENDDO
        ENDDO
      ELSE
!  GREENS FUNCTION PLACEHOLDER
      ENDIF

! Add to existing solution

      DO Q = 1, N_LAYER_WFS
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) + LBSOL(I,O1,Q,1)
            L_WLOWER(I,O1,N,Q) = L_WLOWER(I,O1,N,Q) + LBSOL(I,O1,Q,2)
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_BEAMSOLUTION_NEQK

!

      SUBROUTINE LP_BEAMSOLUTION_NNEK ( &
        FOURIER, IB, N, NV, NV_PARAMETERS, &
        DO_PLANE_PARALLEL, NSTOKES, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        NSTREAMS_2, DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, &
        BVEC, LP_INITIAL_TRANS, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        L_WUPPER, L_WLOWER )

!  Linearization of beam particular integral in one layer only.
!  This is for a layer N not equal to varying layer NV.

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IB
      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NV
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER ::          I, O1, Q
      DOUBLE PRECISION :: CONST, WDEL, TRANS2, VAR1, VAR2, VAR_U

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) &
             .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        RETURN
      ENDIF

!  Classical solution
!  ==================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

!  initialise

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        TRANS2  = CONST * WDEL

!  Distinguish two cases
!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
!  ..(b) plane-parallel and qs for n = 1

        IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
          DO Q = 1, NV_PARAMETERS
            VAR1 = LP_T_DELT_MUBAR(N,NV,IB,Q) * CONST
            VAR2 = LP_INITIAL_TRANS(N,NV,IB,Q)
            DO I = 1, NSTREAMS_2
              DO O1 = 1, NSTOKES
                VAR_U = VAR2 * BVEC(I,O1,N) + LP_BVEC(I,O1,N,NV,Q)
                L_WUPPER(I,O1,N,Q) = CONST  * VAR_U
                L_WLOWER(I,O1,N,Q) = TRANS2 * VAR_U &
                                       + VAR1 * BVEC(I,O1,N)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, NV_PARAMETERS
            VAR1 = LP_INITIAL_TRANS(N,NV,IB,Q) * CONST
            DO I = 1, NSTREAMS_2
              DO O1 = 1, NSTOKES
                L_WUPPER(I,O1,N,Q) = VAR1 * BVEC(I,O1,N)
                L_WLOWER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) * WDEL
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_BEAMSOLUTION_NNEK

!

      SUBROUTINE LP_BVPTEL_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        LAYER_TO_VARY, N_LAYER_WFS, &
        FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, NLAYERS, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        DO_PROFILE_LINEARIZATION, &
        L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        DO_SPECIALIST_OPTION_2, NSTREAMS_2, &
        DELTAU_SLANT, &
        L_DELTAU_VERT, DIRECT_BEAM, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, &
        T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
        LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_BVEC, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
        BANDTELMAT2, IPIVOTTEL, &
        SMAT2, SIPIVOT, &
        L_WLOWER, L_WUPPER, NCON, PCON, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS
      USE LAPACK_TOOLS
      USE VLIDORT_LPC_BVPROBLEM

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOTTEL ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) ::  NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  error tracing variables

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Other local help variables

      INTEGER ::          I, I1, N, N1, NAC, O1, IC, ICOW, Q
      INTEGER ::          K, KO1, K0, K1, K2, C0, NS
      INTEGER ::          IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION :: SPAR, SHOM, L_HOM1, L_HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION :: LXR, MXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, MXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, MXR2, NXR2, PXR2, LLXR2, MLXR2

      DOUBLE PRECISION :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  Initialise

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

      Q = 0

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

      CALL LP_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        LAYER_TO_VARY, N_LAYER_WFS, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        DELTAU_SLANT, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        DIRECT_BEAM, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, &
        BVEC, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        DO_PLANE_PARALLEL, LP_INITIAL_TRANS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )

!  Solve linearized BVP: Several Active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, &
              N_BVTELMATRIX_SUPDIAG, N_LAYER_WFS, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              COLTEL2_WF, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = 'DGBTRS call (NLAYERS>1)in LP_BVPTEL_SOLUTION_MASTER'
         STATUS  = VLIDORT_SERIOUS
         RETURN
        ENDIF

!  Set linearized integration constants, active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTKS_NSTRMS_2

!  set real constants from the solution vector

         DO K = 1, K_REAL(N)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K,N,Q) = COLTEL2_WF(C0+IROW,Q)
            PCON(K,N,Q) = COLTEL2_WF(C0+IROW1,Q)
          ENDDO
         ENDDO

!  set complex constants from the solution vector

         KO1 = K_REAL(N) + 1
         DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K1,N,Q) = COLTEL2_WF(C0+IROW,    Q)
            NCON(K2,N,Q) = COLTEL2_WF(C0+IROW_S,  Q)
            PCON(K1,N,Q) = COLTEL2_WF(C0+IROW1,   Q)
            PCON(K2,N,Q) = COLTEL2_WF(C0+IROW1_S, Q)
          ENDDO
         ENDDO

!  End number of telescoped layers

        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAC = ACTIVE_LAYERS(1)

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, N_LAYER_WFS, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (NLAYERS=1)in LP_BVPTEL_SOLUTION_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set real constants from the solution vector

        DO K = 1, K_REAL(NAC)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K,NAC,Q) = SCOL2_WF(IROW,Q)
            PCON(K,NAC,Q) = SCOL2_WF(IROW1,Q)
          ENDDO
        ENDDO

!  set complex constants from the solution vector

        KO1 = K_REAL(NAC) + 1
        DO K = 1, K_COMPLEX(NAC)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(NAC)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(NAC)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K1,NAC,Q) = SCOL2_WF(IROW,    Q)
            NCON(K2,NAC,Q) = SCOL2_WF(IROW_S,  Q)
            PCON(K1,NAC,Q) = SCOL2_WF(IROW1,   Q)
            PCON(K2,NAC,Q) = SCOL2_WF(IROW1_S, Q)
          ENDDO
        ENDDO

      ENDIF

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Profile linearization
!  *********************

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required if the first active
!       layer is the varying layer

      NAC = ACTIVE_LAYERS(1)
      IF ( NAC .GT. 1 ) THEN

        N1 = NAC - 1

!  Case 1. If active layer is also the varying layer
!  -------------------------------------------------

        IF ( LAYER_TO_VARY.EQ.NAC ) THEN

!  start stream, stokes loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1

!  start parameter loop

           DO Q = 1, N_LAYER_WFS

!  real homogeneous solutions
!    Needs checking------------------------ 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I1,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I1,O1,K,NAC)
              MXR  = MCON(K,NAC)   *   SOLB_XNEG(I1,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I1,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I1,O1,K,NAC,Q)
              L_HOM1 = NXR + LLXR
              L_HOM2 =   T_DELT_EIGEN(K,NAC)   * ( PXR + MLXR ) + &
                       L_T_DELT_EIGEN(K,NAC,Q) *   MXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

!  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I1,O1,K1,NAC) &
                      - NCON(K2,NAC,Q) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC) &
                      - PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC) &
                      + PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)

              MXR1  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K1,NAC) &
                      - MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
              MXR2  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K2,NAC) &
                      + MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I1,O1,K1,NAC,Q) &
                       - LCON(K2,NAC) * L_SOLA_XPOS(I1,O1,K2,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q) &
                       - MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q)
              MLXR2  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q) &
                       + MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q)

              L_HOM1CR = NXR1 + LLXR1
              L_HOM2CR = + T_DELT_EIGEN(K1,NAC)     * ( PXR1 + MLXR1 ) &
                         - T_DELT_EIGEN(K2,NAC)     * ( PXR2 + MLXR2 ) &
                         + L_T_DELT_EIGEN(K1,NAC,Q) *   MXR1 &
                         - L_T_DELT_EIGEN(K2,NAC,Q) *   MXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WUPPER(I1,O1,NAC,Q)
             PCON(ICOW,N1,Q) = SPAR + SHOM
             NCON(ICOW,N1,Q) = ZERO

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  End case 1

        ENDIF

!  Case 2. Active layer is below varying layer
!  -------------------------------------------

        IF ( LAYER_TO_VARY.LT.NAC) THEN

!  start stream, stokes loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1

!  start parameter loop

           DO Q = 1, N_LAYER_WFS

!  real homogeneous solutions
!      Needs Checking ---------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR = NCON(K,NAC,Q)*SOLA_XPOS(I1,O1,K,NAC)
              PXR = PCON(K,NAC,Q)*SOLB_XNEG(I1,O1,K,NAC)
              L_HOM1 = NXR
              L_HOM2 = PXR * T_DELT_EIGEN(K,NAC)
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

!  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I1,O1,K1,NAC) &
                      - NCON(K2,NAC,Q) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC) &
                      - PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC) &
                      + PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)

              L_HOM1CR = NXR1
              L_HOM2CR =  + T_DELT_EIGEN(K1,NAC) * PXR1 &
                          - T_DELT_EIGEN(K2,NAC) * PXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
             ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WUPPER(I1,O1,NAC,Q)
             PCON(ICOW,N1,Q) = SPAR + SHOM
             NCON(ICOW,N1,Q) = ZERO

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  End case 2

        ENDIF

!  Case 3. Active layer is above varying layer
!  -------------------------------------------

        IF ( LAYER_TO_VARY.GT.NAC) THEN

!  No solutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1
           DO Q = 1, N_LAYER_WFS
             PCON(ICOW,N1,Q) = ZERO
             NCON(ICOW,N1,Q) = ZERO
           ENDDO
          ENDDO
         ENDDO

!  End case 3

        ENDIF

      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards
!   Additional linearizations if you are passing through the varying lay

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              NCON(IROW,N,Q) = ZERO
              PCON(IROW,N,Q) = &
                  T_DELT_DISORDS(I,N1) * PCON(IROW,N1,Q)
            ENDDO
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              PCON(IROW,N,Q) = PCON(IROW,N,Q) &
                  + L_T_DELT_DISORDS(I,N1,Q) * MCON(IROW,N1)
            ENDDO
          ENDDO
         ENDDO
        ENDIF
      ENDDO

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!   -- PCON values are zero (no upwelling radiation)
!   -- NCON values propagated downwards from bottom of last active layer

!  layer immediately below Last active layer
!    .... Require linearized solutions at bottom of last active layer
!    .... Additional linearizations is the last active layer is also
!         the varying layer.

      NAC = ACTIVE_LAYERS (NLAYERS_TEL)
      IF ( NAC .LT. NLAYERS ) THEN
        N1 = NAC + 1

!  Case 1. If active layer equals the varying layer
!  ------------------------------------------------

        IF ( LAYER_TO_VARY .EQ. NAC ) THEN

!  start stream, stokes loops

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1

!  start parameter loop

           DO Q = 1, N_LAYER_WFS

!  real homogeneous solutions
!   Needs checking------------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I,O1,K,NAC)
              LXR  = LCON(K,NAC)   *   SOLA_XPOS(I,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I,O1,K,NAC,Q)
              L_HOM2 = PXR + MLXR
              L_HOM1 =   T_DELT_EIGEN(K,NAC)   * ( NXR + LLXR ) + &
                       L_T_DELT_EIGEN(K,NAC,Q) *   LXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

!  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC) &
                      - NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC) &
                      + NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I,O1,K1,NAC) &
                      - PCON(K2,NAC,Q) *   SOLB_XNEG(I,O1,K2,NAC)

              LXR1  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K1,NAC) &
                      - LCON(K2,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
              LXR2  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K2,NAC) &
                      + LCON(K2,NAC) *   SOLA_XPOS(I,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q) &
                       - LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q)
              LLXR2  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q) &
                       + LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I,O1,K1,NAC,Q) &
                       - MCON(K2,NAC) * L_SOLB_XNEG(I,O1,K2,NAC,Q)

              L_HOM2CR = PXR1 + MLXR1
              L_HOM1CR = + T_DELT_EIGEN(K1,NAC)     * ( NXR1 + LLXR1 ) &
                         - T_DELT_EIGEN(K2,NAC)     * ( NXR2 + LLXR2 ) &
                         + L_T_DELT_EIGEN(K1,NAC,Q) *   LXR1 &
                         - L_T_DELT_EIGEN(K2,NAC,Q) *   LXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WLOWER(I,O1,NAC,Q)
!             PCON(IROW,N1,Q) = SPAR + SHOM
!             NCON(IROW,N1,Q) = ZERO
             NCON(IROW,N1,Q) = SPAR + SHOM
             PCON(IROW,N1,Q) = ZERO

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  End case 1

        ENDIF

!  Case 2. Active layer is below varying layer
!  -------------------------------------------

        IF ( LAYER_TO_VARY.LT.NAC) THEN

!  start stream, stokes loops

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1

!  start parameter loop

           DO Q = 1, N_LAYER_WFS

!  real homogeneous solutions
!   Needs checking------------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR = NCON(K,NAC,Q) * SOLA_XPOS(I,O1,K,NAC)
              PXR = PCON(K,NAC,Q) * SOLB_XNEG(I,O1,K,NAC)
              L_HOM2 = PXR
              L_HOM1 = NXR * T_DELT_EIGEN(K,NAC)
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

!  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC) &
                      - NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC) &
                      + NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I,O1,K1,NAC) &
                      - PCON(K2,NAC,Q) *   SOLB_XNEG(I,O1,K2,NAC)

              L_HOM2CR = PXR1
              L_HOM1CR =  + T_DELT_EIGEN(K1,NAC) * NXR1 &
                          - T_DELT_EIGEN(K2,NAC) * NXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
             ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WLOWER(I,O1,NAC,Q)
             NCON(IROW,N1,Q) = SPAR + SHOM
             PCON(IROW,N1,Q) = ZERO

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  End case 2

        ENDIF

!  Case 3. Active layer is above varying layer
!  -------------------------------------------

        IF ( LAYER_TO_VARY.GT.NAC) THEN

!  No solutions

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           DO Q = 1, N_LAYER_WFS
             PCON(IROW,N1,Q) = ZERO
             NCON(IROW,N1,Q) = ZERO
           ENDDO
          ENDDO
         ENDDO

!  End case 3

        ENDIF

      ENDIF

!  other layers to bottom of medium: propagate downwards.
!   Additional variation if you are passing through the varying layer.

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              PCON(IROW,N,Q) = ZERO
              NCON(IROW,N,Q) = &
                  T_DELT_DISORDS(I,N1) * NCON(IROW,N1,Q)
            ENDDO
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              NCON(IROW,N,Q) = NCON(IROW,N,Q) &
                  + L_T_DELT_DISORDS(I,N1,Q) * LCON(IROW,N1)
            ENDDO
          ENDDO
         ENDDO
        ENDIF
      ENDDO

!  debug

!      if ( layer_to_vary .eq. 1 ) then
!       if ( fourier_component .EQ.3 .and.ibeam.eq.1 ) then
!        do n = 1, nlayers
!          write(99,'(a,i3,1p2e18.10)')'hey',n,ncon(3,n,1),pcon(3,n,1)
!        enddo
!       endif
!      endif

!  debug

!      q = 96
!      if ( fourier_component.eq.0.and.ibeam.eq.1
!     &          .and.layer_to_vary.eq.12) then
!        DO N = 1, NLAYERS
!          DO I = 1, NSTREAMS*2
!          write(q,'(3i4,1p4e17.9)')IBEAM,N,I, NCON(I,N,1), PCON(I,N,1)
!         ENDDO
!        ENDDO
!      ENDIF

!  finish

      RETURN
      END SUBROUTINE LP_BVPTEL_SOLUTION_MASTER

!

      SUBROUTINE LP_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        LAYER_TO_VARY, N_LAYER_WFS, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        DELTAU_SLANT, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        DIRECT_BEAM, &
        DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
        DO_LAYER_SCATTERING, &
        T_DELT_MUBAR, INITIAL_TRANS, &
        BVEC, &
        LP_T_DELT_MUBAR, LP_BVEC, &
        DO_PLANE_PARALLEL, LP_INITIAL_TRANS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )

      USE VLIDORT_PARS
      USE VLIDORT_LPC_BVPROBLEM

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT(OUT) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT(OUT) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          STATUS

      INTEGER ::          Q, N, I, I1, IR, CM, C0, IROW, O1, N1, NS
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: CPOS,CNEG,L_HOM_R,L_HOM_CR,L_BEAM
      DOUBLE PRECISION :: T1,T2,T1R,T1I,T2R,T2I,FAC

      DOUBLE PRECISION :: R2_L_BEAM &
          ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Try this safety-first zeroing

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,LAYER_TO_VARY,Q) = ZERO
            L_WLOWER(I,O1,LAYER_TO_VARY,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  status

      status = 0

!  Get the linearized solutions for the layer that is varying
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        IF ( N.EQ.LAYER_TO_VARY ) THEN
          CALL LP_BEAMSOLUTION_NEQK ( &
            FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS, &
            NSTOKES, &
            DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
            NSTREAMS_2, DO_LAYER_SCATTERING, &
            T_DELT_MUBAR, INITIAL_TRANS, &
            BVEC, &
            LP_T_DELT_MUBAR, LP_BVEC, &
            L_WUPPER, L_WLOWER )
        ENDIF
      ENDDO

!  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 3456

!  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        DO Q = 1, MAX_ATMOSWFS
          COLTEL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  top of first active layer, first boundary condition
!  ---------------------------------------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)
      C0 = 0

!  If this active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

!  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contribution

            L_BEAM  = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  otherwise if varying layer is above first active layer, there are bea
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL LP_BEAMSOLUTION_NNEK ( &
          FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS, &
          DO_PLANE_PARALLEL, NSTOKES, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          NSTREAMS_2, DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, &
          BVEC, LP_INITIAL_TRANS, &
          LP_T_DELT_MUBAR, LP_BVEC, &
          L_WUPPER, L_WLOWER )

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(CM,Q) = - L_WUPPER(I,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

!  Intermediate boundaries between active layers
!  ---------------------------------------------

      DO NS = 1, NLAYERS_TEL - 1

!  offsets

       N  = ACTIVE_LAYERS(NS)
       N1 = N + 1
       C0 = NS*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  if N is the varying layer, immediately above boundary

       IF ( N .EQ. LAYER_TO_VARY ) THEN

!  Get the linearized beam solution for the next layer N1

        CALL LP_BEAMSOLUTION_NNEK ( &
          FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS, &
          DO_PLANE_PARALLEL, NSTOKES, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          NSTREAMS_2, DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, &
          BVEC, LP_INITIAL_TRANS, &
          LP_T_DELT_MUBAR, LP_BVEC, &
          L_WUPPER, L_WLOWER )

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  If N1 is the varying layer, immediately below boundary
!    Only require contributions from this layer

       ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

!  Beam contribution

            L_BEAM  = + L_WUPPER(I,O1,N1,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                   L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                   L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                   - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                   + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  non-zero variations if LAYER_TO_VARY is an active layer above N
!    Get the linearized beam solution for the next layer
!  .. contributions from beam solutions on both sides.

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

!  Get the linearized beam solution for the next layer

        CALL LP_BEAMSOLUTION_NNEK ( &
          FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS, &
          DO_PLANE_PARALLEL, NSTOKES, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          NSTREAMS_2, DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, &
          BVEC, LP_INITIAL_TRANS, &
          LP_T_DELT_MUBAR, LP_BVEC, &
          L_WUPPER, L_WLOWER )

!  .. contributions from beam solution (direct assign). No homog. variat

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS
              L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)
              COLTEL2_WF(CM,Q) = L_BEAM
            ENDDO
          ENDDO
        ENDDO

!  Finish different types of boundary condition linearizations

       ENDIF

!  End loop over intermediate active layer boundaries

      ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)
      C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  If this is the surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
             .AND.  N.EQ.NLAYERS ) THEN

!  Modified BCL4M Component loop

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP ( &
           DO_INCLUDE_SURFACE, .TRUE., &
           IBEAM, FOURIER_COMPONENT, &
           SURFACE_FACTOR, N_LAYER_WFS, &
           NSTOKES, NSTREAMS, &
           NLAYERS, DO_LAMBERTIAN_SURFACE, &
           LAMBERTIAN_ALBEDO, BRDF_F, &
           QUAD_STRMWTS, NSTKS_NSTRMS, &
           MUELLER_INDEX, T_DELT_EIGEN, &
           K_REAL, K_COMPLEX, &
           SOLA_XPOS, L_T_DELT_EIGEN, &
           L_SOLA_XPOS, L_SOLB_XNEG, &
           L_WLOWER, &
           R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  start loops

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  If the varying layer is above the surface layer then

        ELSE

!  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP ( &
           DO_INCLUDE_SURFACE, .FALSE., &
           IBEAM, FOURIER_COMPONENT, &
           SURFACE_FACTOR, N_LAYER_WFS, &
           NSTOKES, NSTREAMS, &
           NLAYERS, DO_LAMBERTIAN_SURFACE, &
           LAMBERTIAN_ALBEDO, BRDF_F, &
           QUAD_STRMWTS, NSTKS_NSTRMS, &
           MUELLER_INDEX, T_DELT_EIGEN, &
           K_REAL, K_COMPLEX, &
           SOLA_XPOS, L_T_DELT_EIGEN, &
           L_SOLA_XPOS, L_SOLB_XNEG, &
           L_WLOWER, &
           R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  start loops : beam contributions only

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
             COLTEL2_WF(CM,Q) = - L_BEAM
           ENDDO
          ENDDO
         ENDDO

        ENDIF

!  Otherwise, there's no surface to consider.............

      ELSE

!  If active layer is varying layer, need full calculation.

       IF ( N .EQ. LAYER_TO_VARY ) THEN

!  start loops

        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS

!  Beam contributions

           L_BEAM = - L_WLOWER(I1,O1,N,Q)

!  Linearized Real homogeneous solution contributions

           L_HOM_R  = ZERO
           DO K = 1, K_REAL(N)
             CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
             CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + &
                  L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
             T1 = LCON(K,N) * CPOS
             T2 = MCON(K,N) * CNEG
             L_HOM_R = L_HOM_R + T1 + T2
           ENDDO

!  Linearized Complex homogeneous solution contributions

           L_HOM_CR  = ZERO
           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
             K0 = 2*K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) - &
                  L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
             T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
             T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
             T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
             L_HOM_CR = L_HOM_CR + T1 + T2
           ENDDO

!  Final contribution

           COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  otherwise use beam solution linearizations propagated downwards
!  from the layer that is varying (already computed for this layer)

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

!  start loops

        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(CM,Q) = - L_WLOWER(I1,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

!  Finally end the surface condition

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Only for the Specialist option # 2, and only if the surface
!   is included (Fourier = 0) and the lowest active layer is
!   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * &
                      DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM)
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

!  Continuation point for single layer stuffs
!  ==========================================

 3456 CONTINUE

!  zero column vector. Extremely important.

      DO I = 1, NSTKS_NSTRMS_2
        DO Q = 1, MAX_ATMOSWFS
          SCOL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  top of active layer
!  -------------------

!  If active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

!  Start  loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS

!  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            SCOL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  otherwise if varying layer is above active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL LP_BEAMSOLUTION_NNEK ( &
          FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS, &
          DO_PLANE_PARALLEL, NSTOKES, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          NSTREAMS_2, DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, &
          BVEC, LP_INITIAL_TRANS, &
          LP_T_DELT_MUBAR, LP_BVEC, &
          L_WUPPER, L_WLOWER )

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(IROW,Q) = - L_WUPPER(I,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

!  Bottom of active layer
!  ----------------------

      C0 = NSTKS_NSTRMS

!  If this is the  surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
             .AND.  N.EQ.NLAYERS ) THEN

!  Modified BCL4M Component loop

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP ( &
           DO_INCLUDE_SURFACE, .TRUE., &
           IBEAM, FOURIER_COMPONENT, &
           SURFACE_FACTOR, N_LAYER_WFS, &
           NSTOKES, NSTREAMS, &
           NLAYERS, DO_LAMBERTIAN_SURFACE, &
           LAMBERTIAN_ALBEDO, BRDF_F, &
           QUAD_STRMWTS, NSTKS_NSTRMS, &
           MUELLER_INDEX, T_DELT_EIGEN, &
           K_REAL, K_COMPLEX, &
           SOLA_XPOS, L_T_DELT_EIGEN, &
           L_SOLA_XPOS, L_SOLB_XNEG, &
           L_WLOWER, &
           R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  start loops

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  If the varying layer is above the surface layer then

        ELSE

!  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP ( &
           DO_INCLUDE_SURFACE, .FALSE., &
           IBEAM, FOURIER_COMPONENT, &
           SURFACE_FACTOR, N_LAYER_WFS, &
           NSTOKES, NSTREAMS, &
           NLAYERS, DO_LAMBERTIAN_SURFACE, &
           LAMBERTIAN_ALBEDO, BRDF_F, &
           QUAD_STRMWTS, NSTKS_NSTRMS, &
           MUELLER_INDEX, T_DELT_EIGEN, &
           K_REAL, K_COMPLEX, &
           SOLA_XPOS, L_T_DELT_EIGEN, &
           L_SOLA_XPOS, L_SOLB_XNEG, &
           L_WLOWER, &
           R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

!  start loops : beam contributions only

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
             SCOL2_WF(CM,Q) = - L_BEAM
           ENDDO
          ENDDO
         ENDDO

        ENDIF

!  Otherwise, there's no surface to consider.............

      ELSE

!  If active layer is varying layer, need full calculation.

        IF ( N .EQ. LAYER_TO_VARY ) THEN

!  start loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

!  Beam contributions

            L_BEAM = - L_WLOWER(I1,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            SCOL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

           ENDDO
          ENDDO
         ENDDO

!  otherwise use beam solution linearizations propagated downwards
!  from the layer that is varying (already computed for this layer)

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

!  start loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             SCOL2_WF(CM,Q) = - L_WLOWER(I1,O1,N,Q)
           ENDDO
          ENDDO
         ENDDO

        ENDIF

!  End the surface condition

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Only for the Specialist option # 2, and only if the surface
!   is included (Fourier = 0) and the lowest active layer is
!   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * &
                      DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM)
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                SCOL2_WF(CM,Q) = SCOL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

!  finish

      RETURN
      END SUBROUTINE LP_BVPTEL_COLUMN_SETUP

      END MODULE vlidort_lp_bvproblem

