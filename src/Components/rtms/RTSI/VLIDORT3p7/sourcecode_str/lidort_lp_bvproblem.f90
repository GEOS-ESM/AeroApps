! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #  For linearizations involving atmospheric parameters        #
! #            LP_BVP_SOLUTION_MASTER                           #
! #            LP_BVP_COLUMN_SETUP                              #
! #            LP_BVP_SURFACE_SETUP                             #
! #            LP_BEAMSOLUTION_NEQK                             #
! #            LP_BEAMSOLUTION_NNEK                             #
! #                                                             #
! #  For linearizations (atmospheric) using telescoped BVP      #
! #            LP_BVPTEL_SOLUTION_MASTER                        #
! #            LP_BVPTEL_COLUMN_SETUP                           #
! #    Placeholder, Version 3.3, modify telescoped problem      #
! #                                                             #
! ###############################################################

module lidort_lp_bvproblem

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  LIDORT Use dependencies

   USE lidort_aux,      only : DGBTRS, DGETRS
   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

!private
public

contains

SUBROUTINE LP_BVP_SOLUTION_MASTER &
      ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL,   & ! Input
        DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
        NLAYERS, NSTREAMS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,  & ! Input
        VARIATION_INDEX, N_WEIGHTFUNCS,                               & ! Input
        DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
        SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
        T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
        GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
        L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous solution 
        L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + Thermal
        NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output - Linearized Constants + PI
        STATUS, MESSAGE, TRACE )                                        ! Output - Exception handling

!  Linearization of the Boundary Problem Solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Flag

      LOGICAL  , intent(in)  ::  DO_PLANE_PARALLEL

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)   :: TAYLOR_ORDER

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER_COMPONENT, IBEAM

!  Number of streams

      INTEGER  , intent(in)   :: NSTREAMS, NSTREAMS_2

!  Number of layers

      INTEGER  , intent(in)   :: NLAYERS

!  BVProblem Band matrix control

      INTEGER  , intent(in)   :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Linearization control

      INTEGER  , intent(in)   :: VARIATION_INDEX, N_WEIGHTFUNCS

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  surface factor = 1+delta(m,0). Albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(in)  ::  LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  ::  MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Homogeneous and thermal solution variables
!  ------------------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------


!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  boundary condition flags

      LOGICAL    :: MODIFIED_BCL3, MODIFIED_BCL4

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI, CN

!  Other local help variables 

      INTEGER    :: I, I1, Q, N, AA

!  Initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Linearization of the regular BVP case
!  =====================================

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
      MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

!  Column vector setup

      CALL LP_BVP_COLUMN_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL,   & ! Input, Flags
           DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
           NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL,                        & ! Input, Numbers
           MODIFIED_BCL3, MODIFIED_BCL4, VARIATION_INDEX, N_WEIGHTFUNCS, & ! Input
           DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
           SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
           T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
           GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
           LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
           L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous
           L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
           L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT,  COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) VARIATION_INDEX
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for layer '//CN//'DGBTRS call in LP_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for 1-layer: DGETRS call in LP_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(I,N,Q) = SCOL2_WF(I,Q)
            PCON(I,N,Q) = SCOL2_WF(I1,Q)
          ENDDO
        ENDDO

      ENDIF

!  linearized BVP results
!  ======================

!  Associated quantities

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO AA = 1, NSTREAMS
             DO Q = 1, N_WEIGHTFUNCS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE LP_BVP_SOLUTION_MASTER

!

SUBROUTINE LP_BVP_COLUMN_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL,   & ! Input, Flags
           DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IPARTIC,          & ! Input, Flags and order
           NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL,                        & ! Input, Numbers
           MODIFIED_BCL3, MODIFIED_BCL4, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input
           DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
           SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
           T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
           GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
           LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
           L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous
           L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
           L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  Linearized column vector setup (profile weighting functions)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  ::  DO_PLANE_PARALLEL

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  ::  TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  boundary condition flags

      LOGICAL  , intent(in)  :: MODIFIED_BCL3, MODIFIED_BCL4

!  Linearization control

      INTEGER  , intent(in)  :: LAYER_TO_VARY, N_LAYER_WFS

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IPARTIC

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  Factor, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous
!  -----------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  help variables

      INTEGER    :: Q, AA, N, N1, I, I1, CM, C0, K, M
      REAL(fpk)  :: CPOS, CNEG, L_HOM, L_PARTIC, L_BEAM, FAC

!  Local linearized reflectance arrays

      REAL(fpk)  :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Boundary options

      LOGICAL    :: REGULAR_BCL3, REGULAR_BCL4

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  Proxies (Layer to vary, Fourier)

      M = FOURIER
      K = LAYER_TO_VARY

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          DO N = 1, NLAYERS
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
            L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
          ENDDO
        ENDDO
      ENDIF

!  Get the linearized beam solution for the first layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LP_BEAMSOLUTION_NEQK  &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, K, N_LAYER_WFS, M, IPARTIC,     & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
      ENDIF

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2_WF(I,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

!  No variation case (BCL1)

      ELSE

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            COL2_WF(I,Q) = ZERO
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
          C0  = N1*NSTREAMS_2 - NSTREAMS
          DO I = 1, NSTREAMS_2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTREAMS_2 - NSTREAMS

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PARTIC  = + L_WUPPER(I,N,Q)
            L_HOM = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

        N  = LAYER_TO_VARY
        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, K, N_LAYER_WFS, M, IPARTIC, & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF

!  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_WF(CM,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final layer
!  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

          N1 = N + 1
          C0  = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, K, N_LAYER_WFS, M, IPARTIC, & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
          ENDIF

!  .. contributions from beam solution (direct assign). No homog. variation

          DO I = 1, NSTREAMS_2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
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

        CALL LP_BVP_SURFACE_SETUP                          &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, M,                     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Compute the solution

        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
              CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_WF(CM,Q) = - L_PARTIC - L_HOM
         ENDDO
        ENDDO

!  ordinary BCL6 Component loop
 
      ELSE

!  get the linearized downward-reflected term

        CALL LP_BVP_SURFACE_SETUP                          &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, M,                     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Compute the solution

        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)              
            COL2_WF(CM,Q) = - L_PARTIC
          ENDDO
        ENDDO
        
      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          FAC = - DIRECT_BEAM(I,IPARTIC) * DELTAU_SLANT(N,LAYER_TO_VARY,IPARTIC) 
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
          ENDDO
        ENDDO
      ENDIF

!  debug

!      if ( layer_to_vary.eq.18 ) then
!        write(67,*)layer_to_vary
!        do n = 1, ntotal
!           write(67,*)n,COL2_WF(n,1),COL2_WF(n,2)
!        enddo
!      endif
!      if ( layer_to_vary.eq.18 ) pause

!  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE LP_BVP_COLUMN_SETUP

!

SUBROUTINE LP_BVP_SURFACE_SETUP                            &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, FOURIER_COMPONENT,     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Linearized surface reflectance terms

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  input arguments
!  ---------------

!  BRDF flag

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Number of weighting functions

      INTEGER  , intent(in)  :: N_LAYER_WFS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Flag for type of boundary condition

      LOGICAL  , intent(in)  :: MODIFIED_BCL4

!  overall surface flag and surface factor, albedo

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(FPK), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HSP_U, HSM_U, REFL_P, REFL_B, REFL_M
      REAL(fpk)  :: FACTOR
      INTEGER    :: AA, J, Q, N, I, M

!  Initial section
!  ---------------

!  Always zero the result to start

      DO I = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          R2_L_PARTIC(I,Q)    = ZERO
          DO AA = 1, NSTREAMS
            R2_L_HOMP(I,AA,Q) = ZERO
            R2_L_HOMM(I,AA,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  FOurier component

      M = FOURIER_COMPONENT

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

      N = NLAYERS

!  Particular integral parts

      DO J = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          PV_W(J,Q) = L_WLOWER(J,N,Q) * QUAD_STRMWTS(J)
        ENDDO
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BCL4 ) THEN
        DO J = 1, NSTREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              HSP_U = L_XPOS(J,AA,N,Q) *   T_DELT_EIGEN(AA,N) + &
                        XPOS(J,AA,N)   * L_T_DELT_EIGEN(AA,N,Q)
              HSM_U = L_XNEG(J,AA,N,Q)
              HV_P(J,AA,Q) = QUAD_STRMWTS(J)*HSP_U
              HV_M(J,AA,Q) = QUAD_STRMWTS(J)*HSM_U
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Integrated Downward reflection (Calculation, Lambertian case)
!  -------------------------------------------------------------

     IF ( .not. DO_BRDF_SURFACE ) THEN

!  amplitude

        FACTOR = SURFACE_FACTOR * ALBEDO

!  Only a solution for the isotropic part.

        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, N_LAYER_WFS

!  Particular solution

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + PV_W(J,Q)
            ENDDO
            REFL_B = REFL_B * FACTOR
            DO I = 1, NSTREAMS
              R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B
            ENDDO

!  Homogeneous solutions (only for modified BC)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
                  REFL_P = REFL_P + HV_P(J,AA,Q)
                  REFL_M = REFL_M + HV_M(J,AA,Q)
                ENDDO
                REFL_P = REFL_P * FACTOR
                REFL_M = REFL_M * FACTOR
                DO I = 1, NSTREAMS
                  R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                  R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
                ENDDO
              ENDDO
            ENDIF

!  end parameter loop

          ENDDO
        ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!  ----------------------------------------------------------------

      ELSE

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS

!  particular solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

            REFL_B = ZERO
            DO J = 1, NSTREAMS
!              REFL_B = REFL_B + PV_W(J,Q) * BRDF_F(M,J,I)
              REFL_B = REFL_B + PV_W(J,Q) * BRDF_F(M,I,J)
            ENDDO
            REFL_B = REFL_B * SURFACE_FACTOR
            R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B

!  homogeneous solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
!                  REFL_P = REFL_P + HV_P(J,AA,Q) * BRDF_F(M,J,I)
!                  REFL_M = REFL_M + HV_M(J,AA,Q) * BRDF_F(M,J,I)
                  REFL_P = REFL_P + HV_P(J,AA,Q) * BRDF_F(M,I,J)
                  REFL_M = REFL_M + HV_M(J,AA,Q) * BRDF_F(M,I,J)
                ENDDO
                REFL_P = REFL_P * SURFACE_FACTOR
                REFL_M = REFL_M * SURFACE_FACTOR
                R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
              ENDDO
            ENDIF

!  end parameter and stream loops

          ENDDO
        ENDDO

!  End BRDF clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_BVP_SURFACE_SETUP

!

SUBROUTINE LP_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, N_WEIGHTFUNCS, M, IB,        & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N
!   Profile variation also in this layer (N = K)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  ::  M, IB, N

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

! mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Local linearized Green's functioN multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT
      REAL(fpk)  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)

!  No linearized particular solution beyond the cutoff layer. ALSO -
!  Nothing if layer is inactive (Does not depend on solution saving)
!    Exit (Solutions have not necessarily already been zeroed).

!      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N)) &
!             .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
      IF (.NOT.DO_LAYER_SCATTERING(M,N)  .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_WEIGHTFUNCS
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Green's function solution
!  =========================

!  Set up linearizations of GAMMA constants
!  ----------------------------------------

!  Distinguish two cases:
!  ..(a) quasi-spherical for n > 1
!  ..(b) plane-parallel or QS for n=1

!  Linearizations of optical depth integrations
!  Linearized Green function multipliers

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

      DO AA = 1, NSTREAMS

         ZDEL  = T_DELT_EIGEN(AA,N)
         ZWDEL = ZDEL * WDEL
         AST   = CONST * ATERM_SAVE(AA,N) 
         BST   = CONST * BTERM_SAVE(AA,N) 

!  Downwelling, Make allowances for Taylor series

        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS   = GAMMA_M(AA,N)
           DELTA = DELTAU_VERT(N)
           LAM   = AVERAGE_SECANT(N,IB)
           DO Q = 1, N_WEIGHTFUNCS
              L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_KEG  = L_KEIGEN(AA,N,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
           ENDDO
        ELSE
           MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
           DO Q = 1, N_WEIGHTFUNCS
              L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
              L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(AA,N)
           ENDDO
        ENDIF
        DO Q = 1, N_WEIGHTFUNCS
           L_AST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_ATERM_SAVE(AA,N,Q)
           L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
        ENDDO

!  Upwelling

        MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
        DO Q = 1, N_WEIGHTFUNCS
           L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q)   ; L_KEG  = L_KEIGEN(AA,N,Q)
           L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
           L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(AA,N)
        ENDDO
        DO Q = 1, N_WEIGHTFUNCS
           L_BST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_BTERM_SAVE(AA,N,Q)
           L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
        ENDDO

!  End discrete ordinate loop

      ENDDO

!  Set linearized form of particular integral at boundaries
!  --------------------------------------------------------

      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          S_P_U = ZERO
          S_P_L = ZERO
          S_M_U = ZERO
          S_M_L = ZERO
          DO AA = 1, NSTREAMS
            S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N) + &
                              GFUNC_UP(AA,N) * L_XPOS(I1,AA,N,Q)
            S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N) +  &
                              GFUNC_UP(AA,N) * L_XPOS(I,AA,N,Q)
            S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N) +  &
                              GFUNC_DN(AA,N) * L_XPOS(I,AA,N,Q)
            S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N) + &
                              GFUNC_DN(AA,N) * L_XPOS(I1,AA,N,Q)
          ENDDO
!mick fix 2/17/12 - switched for additional code below
          LBSOL(I,Q,1)  = S_P_U
          LBSOL(I1,Q,1) = S_M_U
          LBSOL(I1,Q,2) = S_M_L
          LBSOL(I,Q,2)  = S_P_L
        ENDDO
      ENDDO

!mick fix 2/17/12 - code added to include both solar & thermal sources.
!                   like LC_BEAMSOLUTION_NEQK now
! Add to existing solution

      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS_2
          L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
          L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BEAMSOLUTION_NEQK

!

SUBROUTINE LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, K, K_PARAMETERS, M, IB,      & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N
!   Profile variation Not in this layer (N =/ K)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  :: M, IB, N

!  Varying lyaer, number of varying parameters (input)

      INTEGER  , intent(in)  :: K, K_PARAMETERS

!  Input optical properties after delta-M scaling. 

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Local linearized Green's function multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT, T1
      REAL(fpk)  :: L_WDEL, L_LAM, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

!  Green's function solution
!  =========================

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if layer is inactive (Does not depend on solution saving)
!    Zero the solutions and exit.

      IF (.NOT.DO_LAYER_SCATTERING(M,N) .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS_2
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Set up linearizations of GAMMA constants

!  Mick note 9/25/2013 - LP_GAMMA_P & LP_GAMMA_M were previously
!  normalized by GAMMA_P & GAMMA_M, respectively.  With the new
!  definitions of GAMMA_P & GAMMA_M (now straight sums and differences
!  instead of their corresponding reciprocals), we define LP_GAMMA_M as
!  follows:

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

        DO AA = 1, NSTREAMS

          ZDEL  = T_DELT_EIGEN(AA,N)
          ZWDEL = ZDEL * WDEL
          AST   = CONST * ATERM_SAVE(AA,N) 
          BST   = CONST * BTERM_SAVE(AA,N) 

!mick fix 9/4/2013 - small numbers analysis added

!  Downwelling, Make allowances for Taylor series

          IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
            EPS   = GAMMA_M(AA,N)
            DELTA = DELTAU_VERT(N)
            LAM   = AVERAGE_SECANT(N,IB)
            DO Q = 1, K_PARAMETERS
              L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, ZERO, ZERO, L_LAM, WDEL, LAM, L_MULT(Q) )
            ENDDO
          ELSE
            MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
            DO Q = 1, K_PARAMETERS
              L_WDEL =  LP_T_DELT_MUBAR (N,K,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
              L_MULT(Q) = ( - L_WDEL - MULT * L_LAM ) / GAMMA_M(AA,N)
            ENDDO
          ENDIF
          DO Q = 1, K_PARAMETERS
            L_AST =  LP_INITIAL_TRANS(N,K,IB,Q) 
            L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
          ENDDO

!  Upwelling

          MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
          DO Q = 1, K_PARAMETERS
            L_WDEL =  LP_T_DELT_MUBAR (N,K,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
            L_MULT(Q) = - ( L_WDEL*ZDEL + MULT*L_LAM ) / GAMMA_P(AA,N)
          ENDDO
          DO Q = 1, K_PARAMETERS
            L_BST =  LP_INITIAL_TRANS(N,K, IB,Q)
            L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
          ENDDO

!  End eigenloop

        ENDDO

!  ..(b) plane-parallel and qs for n = 1

      ELSE

        DO AA = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            T1 = LP_INITIAL_TRANS(N,K,IB,Q)
            L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * T1
            L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * T1
          ENDDO
        ENDDO

      ENDIF

!  Set linearized form of particular integral

      DO Q = 1, K_PARAMETERS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          S_P_U = ZERO
          S_P_L = ZERO
          S_M_U = ZERO
          S_M_L = ZERO
          DO AA = 1, NSTREAMS
            S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N)
            S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N)
            S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N)
            S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N)
          ENDDO
!mick fix 2/17/12 - switched for additional code below
          LBSOL(I,Q,1)  = S_P_U
          LBSOL(I1,Q,1) = S_M_U
          LBSOL(I1,Q,2) = S_M_L
          LBSOL(I,Q,2)  = S_P_L
        ENDDO
      ENDDO

!mick fix 2/17/12 - code added to include both solar & thermal sources.
! Add to existing solution

      DO Q = 1, K_PARAMETERS
        DO I = 1, NSTREAMS_2
          L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
          L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BEAMSOLUTION_NNEK

!

SUBROUTINE LP_BVPTEL_SOLUTION_MASTER                              &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, FOURIER,        & ! Input, Flags, Order
             IBEAM, NLAYERS, NSTREAMS, NSTREAMS_2, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input, Numbers
             N_SUPDIAG, N_SUBDIAG, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE, & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                  & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam Quantities
             BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,                        & ! Input, BVProblem
             T_DELT_DISORDS, L_T_DELT_DISORDS,  LCON_XVEC, MCON_XVEC,       & ! Input, Extinction, non active layers 
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                          & ! Input, Homogeneous + constants
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE,                                    & ! Input, Linearized Greens
             L_WUPPER, L_WLOWER, NCON, PCON, NCON_XVEC, PCON_XVEC,          & ! Output
             STATUS, MESSAGE, TRACE )                                         ! Output

!  Linearization of the Telescoped Boundary Problem Solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                              ZERO, LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: layer varying and number of varying parameters (input)

      INTEGER  , intent(in)  ::  LAYER_TO_VARY, N_LAYER_WFS

!  BVProblem Band matrix control

      INTEGER  , intent(in)  ::  N_SUBDIAG, N_SUPDIAG

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2       (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDTELMAT2 (MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT    (MAXSTREAMS_2)

!  discrete ordinate factors (BVP telescoping, solutions saving)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER    :: I, I1, K, K1, Q, M
      INTEGER    :: NS, N, N1, NAF, NAL, AA, C0, CM, CP
      REAL(fpk)  :: SHOM, L_HOM1, L_HOM2

!  Initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Wf count

      M = FOURIER
      Q = 0

!  Set up linearized BVP, column vector  B' where AX = B'
!  ======================================================

!  Bulk: Compute the main column B' where AX = B'

      CALL LP_BVPTEL_COLUMN_SETUP &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, NLAYERS, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,    & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,                & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                  & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                          & ! Input, Homogeneous solution
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input, Linearized Homogeneous 
             L_ATERM_SAVE, L_BTERM_SAVE,                                    & ! Input, Linearized Greens 
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                       ! Output

!  Solve linearized BVP, several active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, N_LAYER_WFS, &
               BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2_WF, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = 'DGBTRS call in L_BVPTEL_BEAMSOLUTION_MASTER'
         STATUS  = LIDORT_SERIOUS
         RETURN
        ENDIF

!  Set linearized integration constants, active layers

        C0 = - NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTREAMS_2
         DO I = 1, NSTREAMS
          CM = C0 + I
          CP = CM + NSTREAMS
          DO Q = 1, N_LAYER_WFS
           NCON(I,N,Q) = COLTEL2_WF(CM,Q)
           PCON(I,N,Q) = COLTEL2_WF(CP,Q)
          ENDDO
         ENDDO
        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ('N', NSTREAMS_2, N_LAYER_WFS,  & 
                    SMAT2, MAXSTREAMS_2, SIPIVOT, SCOL2_WF, MAXSTREAMS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in LP_BVPTEL_BEAMSOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
        ENDIF

!  Set linearized integration constants for active layer

        N = ACTIVE_LAYERS(1)
        DO K = 1, NSTREAMS
         K1 = K + NSTREAMS 
         DO Q = 1, N_LAYER_WFS
           NCON(K,N,Q) = SCOL2_WF(K,Q)
           PCON(K,N,Q) = SCOL2_WF(K1,Q)
         ENDDO
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              NCON_XVEC(I,K,N,Q) = NCON(K,N,Q) * XPOS(I,K,N)
              PCON_XVEC(I,K,N,Q) = PCON(K,N,Q) * XNEG(I,K,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required if the first active
!       layer is the varying layer

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        IF ( LAYER_TO_VARY.EQ.NAF ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM1 = NCON_XVEC(I1,AA,NAF,Q) + LCON(AA,NAF)*L_XPOS(I1,AA,NAF,Q)
                L_HOM2 = T_DELT_EIGEN(AA,NAF) *          &
                  ( PCON_XVEC(I1,AA,NAF,Q)             + & 
                    MCON(AA,NAF)*L_XNEG(I1,AA,NAF,Q) ) + &
               L_T_DELT_EIGEN(AA,NAF,Q) * MCON_XVEC(I1,AA,NAF)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY.LT.NAF) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM1 = NCON_XVEC(I1,AA,NAF,Q)
                L_HOM2 = T_DELT_EIGEN(AA,NAF) * PCON_XVEC(I1,AA,NAF,Q)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N1,Q) = ZERO
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  For remaining non-active atmospheri! layers to TOA, propagate upwards.
!   Additional linearizations if you are passing through the varying layer.

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            NCON(I,N,Q) = ZERO
            PCON(I,N,Q) = T_DELT_DISORDS(I,N1) * PCON(I,N1,Q)
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N,Q) = PCON(I,N,Q) + L_T_DELT_DISORDS(I,N1,Q) * MCON(I,N1)
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
!    .... Additional linearizations if the last active layer is also
!         the varying layer.

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        IF ( LAYER_TO_VARY .EQ. NAL ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM2 = PCON_XVEC(I,AA,NAL,Q) + MCON(AA,NAL)*L_XNEG(I,AA,NAL,Q)
                L_HOM1 = T_DELT_EIGEN(AA,NAL) *          &
                   ( NCON_XVEC(I,AA,NAL,Q)             + &
                     LCON(AA,NAL)*L_XPOS(I,AA,NAL,Q) ) + &
                L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I,AA,NAL)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
              PCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY .LT. NAL ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM2 = PCON_XVEC(I,AA,NAL,Q) 
                L_HOM1 = T_DELT_EIGEN(AA,NAL) * NCON_XVEC(I,AA,NAL,Q)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
              PCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N1,Q) = ZERO
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  other layers to bottom of medium: propagate downwards.
!   Additional variation if you are passing through the varying layer.

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            PCON(I,N,Q) = ZERO
            NCON(I,N,Q) = T_DELT_DISORDS(I,N1) * NCON(I,N1,Q)
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              NCON(I,N,Q) = NCON(I,N,Q) + L_T_DELT_DISORDS(I,N1,Q) * LCON(I,N1)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Associated quantities for inactive layers
!  -----------------------------------------

!  atmosphere layers with no scattering

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
             DO Q = 1, N_LAYER_WFS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
             ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE LP_BVPTEL_SOLUTION_MASTER

!

SUBROUTINE LP_BVPTEL_COLUMN_SETUP                                 &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, NLAYERS, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,    & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,                & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                  & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                          & ! Input, Homogeneous solution
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input, Linearized Homogeneous 
             L_ATERM_SAVE, L_BTERM_SAVE,                                    & ! Input, Linearized Greens 
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                       ! Output

!  Column setup for the linearized telescoped BVP

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: layer and number of varying parameters (input)

      INTEGER  , intent(in)  ::  LAYER_TO_VARY, N_LAYER_WFS

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  M, IBEAM

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER    :: Q,AA,N,N1,NS,I,I1,CM,C0,NAF
      REAL(fpk)  :: CPOS, CNEG, L_HOM, L_BEAM

!  Try this safety-first zeroing

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          DO N = 1, NLAYERS
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Get the linearized solutions for all active layers
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        IF ( N.EQ.LAYER_TO_VARY ) THEN
          CALL LP_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, N_LAYER_WFS, M, IBEAM,       & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF
      ENDDO

!  special case for only 1 active layer

      IF ( NLAYERS_TEL .GT. 1 ) THEN

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
       NAF = N

!  If this active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

       IF ( LAYER_TO_VARY .EQ. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COLTEL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

!  otherwise if varying layer is above first active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,  & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
             COLTEL2_WF(I,Q) = - L_WUPPER(I,N,Q)
          ENDDO
        ENDDO

       ENDIF

!  Intermediate boundaries between active layers
!  ---------------------------------------------

       DO NS = 1, NLAYERS_TEL - 1

!  offsets

        N  = ACTIVE_LAYERS(NS)
        N1 = N + 1
        C0 = NS*NSTREAMS_2 - NSTREAMS

!  if N is the varying layer, immediately above boundary
!  Get the linearized beam solution for the next layer N1

        IF ( N .EQ. LAYER_TO_VARY ) THEN

         CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM, & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

         DO I = 1, NSTREAMS_2
           CM = C0 + I
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
             L_HOM  = ZERO
             DO AA = 1, NSTREAMS
               CNEG = L_XNEG(I,AA,N,Q)
               CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + &
                    L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
               L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
             ENDDO
             COLTEL2_WF(CM,Q) = L_BEAM - L_HOM
           ENDDO
         ENDDO

!  If N1 is the varying layer, immediately below boundary
!    Only require contributions from this layer

        ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

         DO I = 1, NSTREAMS_2
           CM = C0 + I
           DO Q = 1, N_LAYER_WFS
             L_BEAM  = + L_WUPPER(I,N1,Q)
             L_HOM = ZERO
             DO AA = 1, NSTREAMS
               CPOS = L_XPOS(I,AA,N1,Q)
               CNEG = T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + &
                    L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
               L_HOM = L_HOM + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
             ENDDO
             COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
           ENDDO
         ENDDO

!  non-zero variations if LAYER_TO_VARY is an active layer above N
!    Get the linearized beam solution for the next layer
!  .. contributions from beam solutions on both sides.

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

         CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM, & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

         DO I = 1, NSTREAMS_2
           CM = C0 + I
           DO Q = 1, N_LAYER_WFS
             COLTEL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
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
       C0 = (NS-1)*NSTREAMS_2 + NSTREAMS

!  If last active layer is varying

       IF ( N .EQ. LAYER_TO_VARY ) THEN

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CNEG = L_XNEG(I1,AA,N,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

!  If varying layer is above layer and is active, beam contribution

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            COLTEL2_WF(CM,Q) = - L_BEAM
          ENDDO
        ENDDO

       ENDIF

!  Special case, single-active-layer case

      ELSE

!  zero column vector
!   Probably not needed

       DO I = 1, NSTREAMS_2
         DO Q = 1, MAX_ATMOSWFS
           SCOL2_WF(I,Q) = ZERO
         ENDDO
       ENDDO

!  top of active layer
!  -------------------

       NS = 1
       N = ACTIVE_LAYERS(NS)

!  If active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

       IF ( LAYER_TO_VARY .EQ. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
             CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            SCOL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

!  otherwise if varying layer is above active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,  & ! Input, Numbers
             DELTAU_VERT, LAYER_PIS_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = - L_WUPPER(I,N,Q)
          ENDDO
        ENDDO

       ENDIF

!  Bottom of active layer
!  ----------------------

!  If active layer is varying layer, need full calculation.

       IF ( N .EQ. LAYER_TO_VARY ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CNEG = L_XNEG(I1,AA,N,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

!  otherwise use beam solution linearizations propagated downwards
!  from the layer that is varying (already computed for this layer)

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            I1 = I + NSTREAMS
            SCOL2_WF(I1,Q) = - L_WLOWER(I1,N,Q)
          ENDDO
        ENDDO

       ENDIF

!  End layer clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE LP_BVPTEL_COLUMN_SETUP

!  End

end module lidort_lp_bvproblem
