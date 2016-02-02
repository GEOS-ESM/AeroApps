! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
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

! ###########################################################
! #                                                         #
! #  2OS History :                                          #
! #                                                         #
! #     Mark 1: 2007 for OCO-Mk1, publication               #
! #     Mark 2: 2009, with BRDFs                            #
! #     Mark 3: 2013, Re-engineered Model:                  #
! #              * Including Multiple Geometry              #
! #              * Merging with Stand-alone FO code         #
! #              * New linearization for Bulk properties    #
! #     Mark 4: 2014, integration with VLIDORT              #
! #                                                         #
! #                                                         #
! ###########################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_2OSCORR_LPS_MASTER (master)              #
! #                                                             #
! ###############################################################


      MODULE vlidort_2OScorr_lps_master_m

      USE VLIDORT_PARS
      USE vlidort_2OScorr_utilities
      USE vlidort_2OScorr_routines
      USE vlidort_2OScorr_la_routines
      USE vlidort_2OScorr_lps_routines

      PRIVATE
      PUBLIC :: VLIDORT_2OSCORR_LPS_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_2OSCORR_LPS_MASTER ( &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup, &
        VLIDORT_Out, &
        VLIDORT_LinFixIn, &
        VLIDORT_LinModIn, &
        VLIDORT_LinSup, &
        VLIDORT_LinOut )

      USE VLIDORT_PARS

      USE VLIDORT_Inputs_def
      USE VLIDORT_Sup_InOut_def
      USE VLIDORT_Outputs_def
      USE VLIDORT_LinInputs_def
      USE VLIDORT_LinSup_InOut_def
      USE VLIDORT_LinOutputs_def

      IMPLICIT NONE

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN) :: VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)    :: VLIDORT_Sup

!  VLIDORT output structure. This must be Intent(InOut)

      TYPE(VLIDORT_Outputs), INTENT (INOUT)      :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs), INTENT (IN)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (IN) :: VLIDORT_LinModIn

!  VLIDORT linearized supplements structure

      TYPE(VLIDORT_LinSup_InOut), INTENT (INOUT)    :: VLIDORT_LinSup

!  VLIDORT linearized output structure. This must be Intent(InOut)

      TYPE(VLIDORT_LinOutputs), INTENT (INOUT)      :: VLIDORT_LinOut

!  *******************
!  LOCAL 2OS Variables
!  *******************

!  Inputs/Outputs
!  ==============

!  MaxGeoms parameter
!    This could be set to MAX_GEOMETRIES or MAXBEAMS from VLIDORT_PARS
!    IF you set it to MAX_Geometries, then there could be redundancy.
!    To be set independently here, then checked against MAXBEAMS

   integer, parameter :: MaxGeoms = 5

!  Boolean flags

   logical :: do_deltam_scaling
   logical :: do_plane_parallel
   logical :: do_Lambertian

!  Linearization control
!   Note that we distinguish betwen the VLIDORT and 2OS profile WF indexing

   logical :: do_LP_Jacobians
   logical :: do_LS_Jacobians
   integer :: npars_vli(MAXLAYERS),npars(MAXLAYERS), nspars

!  Control integers

   integer :: nlayers, nstokes, nstreams, ngeoms, nfoumax 
   integer :: nmoments, ncoeffs(MAXLAYERS) 

!  Quadrature

   real(fpk) :: qstreams(MAXSTREAMS), qweights(MAXSTREAMS)

!  Fourier convergence

   real(fpk) :: convergence_eps

!  Geometries stored in array Geoms(i,j,*)
!    i = 1/2/3 = SZA/VZA/AZM, j = 1/2/3/4 = Degrees/Radians/Cosines/sines

   real(fpk) :: geoms(3,4,MaxGeoms)

!  Chapman factors

   real(fpk) :: SunChapman (MAXLAYERS,MAXLAYERS,MaxGeoms)

!  optical properties (these will be deltam-scaled)
!     omegas_1 = omegas * ( 1.0 - truncfac )
!     omegas_2 = omegas * truncfac

   real(fpk) :: opdeps(MAXLAYERS)
   real(fpk) :: omegas(MAXLAYERS)
   real(fpk) :: truncfac(MAXLAYERS)
   real(fpk) :: omegas_1(MAXLAYERS)
   real(fpk) :: omegas_2(MAXLAYERS)
   real(fpk) :: coeffs(0:MAXMOMENTS,MAXLAYERS,6)

!  Linearized optical properties 

   real(fpk) :: LA_opdeps(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: LA_omegas(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: LA_truncfac(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: LA_omegas_1(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: LA_omegas_2(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: LA_coeffs(0:MAXMOMENTS,MAXLAYERS,6,MAX_ATMOSWFS)

!  Surface stuff. Either Lambertian or Pre-prepared BRDF from VLIDORT_Sup

   real(fpk) :: albedo

!  First-order results, only for Deltam-scaling

   real(fpk) :: FO_Zmatrix(MAXLAYERS,4,MaxGeoms)
   real(fpk) :: FO_R1saved(0:MAXLAYERS,4,MaxGeoms)

!  R2 reflectance and Icorr. The Main results!

   real(fpk) :: R2(4,MaxGeoms)
   real(fpk) :: Icorr(MaxGeoms)

   real(fpk) :: LP_R2   (4,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk) :: LS_R2   (4,MAX_SURFACEWFS,MaxGeoms)

   real(fpk) :: LP_Icorr(MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk) :: LS_Icorr(MAX_SURFACEWFS,MaxGeoms)

!  Local 2OS variables
!  ===================

!  Average secant output

   real(fpk) :: avg_secants   (MAXLAYERS,MaxGeoms)
   real(fpk) :: LP_avg_secants(MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  Local output from subroutines "Calculate_Multipliers" and "Calculate_Multipliers_LPPLUS"
!  ----------------------------------------------------------------------------------------

!  First-order Transmittances and multipliers. Only defined for Quadrature directions.

   real(fpk)  :: O1_QVTrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk)  :: O1_QATrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)

   real(fpk)  :: O1_QVMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk)  :: O1_QAMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  Second-order Transmittances and multipliers

   real(fpk)  :: O2_AVTrans(MAXLAYERS,MaxGeoms)
   real(fpk)  :: O2_AVMult(MAXLAYERS,MaxGeoms)

   real(fpk)  :: O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk)  :: O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk)  :: O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk)  :: O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  Linearized First-order Transmittances and multipliers

   real(fpk)  :: LP_O1_QVTrans(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O1_QVMult (MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk)  :: LP_O1_QATrans(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O1_QAMult (MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  Linearized Second-order Transmittances and multipliers

   real(fpk)  :: LP_O2_AVTrans(MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O2_AVMult (MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk)  :: LP_O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk)  :: LP_O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  Local output from "Calculate_Zmatrix" and "Calculate_Zmatrix_LAPLUS"
!  -------------------------------------

   real(fpk)  :: Ptc(2,Maxstreams,4,4)
   real(fpk)  :: Pts(2,Maxstreams,4,4)
   real(fpk)  :: Prc(2,Maxstreams,4,4)
   real(fpk)  :: Prs(2,Maxstreams,4,4)

!  Linearized source phase matrices

   real(fpk)  :: LP_Ptc(2,Maxstreams,4,4,MAX_ATMOSWFS)
   real(fpk)  :: LP_Pts(2,Maxstreams,4,4,MAX_ATMOSWFS)
   real(fpk)  :: LP_Prc(2,Maxstreams,4,4,MAX_ATMOSWFS)
   real(fpk)  :: LP_Prs(2,Maxstreams,4,4,MAX_ATMOSWFS)

!  Local Reflectance matrices
!  --------------------------

!  First-order Fourier component variables

   real(fpk)  :: R1cscal(2,MAXSTREAMS)
   real(fpk)  :: R1c    (2,MAXSTREAMS,4,4)
   real(fpk)  :: R1s    (2,MAXSTREAMS,4,4)

   real(fpk)  :: LP_R1cscal(2,MAXSTREAMS,    MAXLAYERS,MAX_ATMOSWFS)
   real(fpk)  :: LP_R1c    (2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk)  :: LP_R1s    (2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)

   real(fpk)  :: LS_R1cscal(2,MAXSTREAMS,    MAX_SURFACEWFS)
   real(fpk)  :: LS_R1c    (2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk)  :: LS_R1s    (2,MAXSTREAMS,4,4,MAX_SURFACEWFS)

!  Second-order Fourier component variables

   real(fpk)  :: R2cscal, LP_R2cscal  ( MAXLAYERS, MAX_ATMOSWFS ), LS_R2cscal  ( MAX_SURFACEWFS )
   real(fpk)  :: R2c(4),  LP_R2c   ( 4, MAXLAYERS, MAX_ATMOSWFS ), LS_R2c   ( 4, MAX_SURFACEWFS )
   real(fpk)  :: R2s(4),  LP_R2s   ( 4, MAXLAYERS, MAX_ATMOSWFS ), LS_R2s   ( 4, MAX_SURFACEWFS )

!  Local help variables
!  ********************

!  Legendre setups (formerly saved routines)
!     Setups for the Zmatrix thread-safe alternative

   integer   :: ll, msq
   real(fpk) :: bin2mm,bin2m2,binf2,help,xmu_m1,xmu_p1, urootm
   integer   :: Tns1, Tns2
   real(fpk) :: Trl(0:2*MAXMOMENTS)
   real(fpk) :: Tr2lp1(0:MAXMOMENTS)
   real(fpk) :: Trlsq(0:MAXMOMENTS)
   real(fpk) :: Tsqlm(0:MAXMOMENTS),Tsql4(0:MAXMOMENTS)
   real(fpk) :: Tqroot6, Ttwom, Trmsq, Tbinfac
   real(fpk) :: Txmu(MAXSTREAMS_p2) 
   real(fpk) :: Txsi(MAXSTREAMS_p2)
   real(fpk) :: Txmu_m1_sq(MAXSTREAMS_p2)
   real(fpk) :: Txmu_p1_sq(MAXSTREAMS_p2)
   real(fpk) :: TPlm_2all(MAXSTREAMS_p2,5)
   real(fpk) :: TPlm_mlp(MAXSTREAMS_p2)
   real(fpk) :: TPlm_mlm(MAXSTREAMS_p2)

!  help

   logical   :: almost
   logical   :: nextm
   integer   :: m, ndcoeffs, N, N1, L, i, p, p1, k, O1
   real(fpk) :: xa, xv, xs, muv, mus, phi, W1, W2, LP_xa(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: piconv, fluxfac, conversion

!  Debug testing 23-24 October 2014. Validate against 2013 results.

!   logical, parameter :: do_debug_23oct14 = .true.
   logical, parameter :: do_debug_23oct14 = .false.
   logical            :: regular_ps,enhanced_ps
   integer            :: nd, ld, kd
   real(fpk)          :: dtr, dn, z, aa, ray1, ray2
   real(fpk)          :: molabs(MaxLayers) 
   real(fpk)          :: molsca(MaxLayers) 
   real(fpk)          :: aerext(MaxLayers) 
   real(fpk)          :: opdeps_r(MaxLayers) 
   real(fpk)          :: omegas_r(MaxLayers) 

!  verbosity parameter

!   logical, parameter :: verbo = .true.
   logical, parameter :: verbo = .false.

!  Exception handling

   integer       :: nmessages
   character*100 :: messages, actions

!  Input Check section
!  ===================

!mick fix 4/21/2015 - added status initializations
   VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK  = VLIDORT_SUCCESS
   VLIDORT_Out%Status%TS_2OSCORR_NCHECKMESSAGES = 0
   VLIDORT_Out%Status%TS_2OSCORR_CHECKMESSAGES(0:0) = ' '
   VLIDORT_Out%Status%TS_2OSCORR_ACTIONS(0:0)       = ' '

   VLIDORT_Out%Status%TS_2OSCORR_STATUS_CALCULATION = VLIDORT_SUCCESS
   VLIDORT_Out%Status%TS_2OSCORR_MESSAGE = ' '
   VLIDORT_Out%Status%TS_2OSCORR_TRACE_1 = ' '

   nmessages = 0

!  check parameter

   if ( MaxGeoms.ne.MAX_SZANGLES ) then
      messages = 'MaxGeoms parameter in 2OS correction not equal to MAX_SZANGLES parameter in VLIDORT'
      actions  = 'Reset MaxGeoms parameter in 2OS correction'
      nmessages = nmessages + 1
      VLIDORT_Out%Status%TS_2OSCORR_CHECKMESSAGES(nmessages)     = Adjustl(Trim(messages))
      VLIDORT_Out%Status%TS_2OSCORR_ACTIONS      (nmessages)     = Adjustl(Trim(actions))
   endif

!  Check observational Geometry flag is turned on

   if ( .not. VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) then
      messages = '2OS correction must work with Observation Geometry flag set'
      actions  = 'Reset DO_OBSERVATION_GEOMETRY in VLIDORT inputs'
      nmessages = nmessages + 1
      VLIDORT_Out%Status%TS_2OSCORR_CHECKMESSAGES(nmessages)     = Adjustl(Trim(messages))
      VLIDORT_Out%Status%TS_2OSCORR_ACTIONS      (nmessages)     = Adjustl(Trim(actions))
   endif

!  Check upwelling 

   if ( .not. VLIDORT_FixIn%Bool%TS_DO_UPWELLING ) THEN
      messages = '2OS correction only works for upwelling radiation'
      actions  = 'Reset DO_UPWELLING flag to TRUE in VLIDORT inputs'
      nmessages = nmessages + 1
      VLIDORT_Out%Status%TS_2OSCORR_CHECKMESSAGES(nmessages)     = Adjustl(Trim(messages))
      VLIDORT_Out%Status%TS_2OSCORR_ACTIONS      (nmessages)     = Adjustl(Trim(actions))
   endif

!  Check TOA

   if ( VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) .ne. zero ) then
      messages = '2OS correction only works for upwelling radiation at TOA'
      actions  = 'Reset USERLEVELS(1) = 0.0d0 in VLIDORT inputs'
      nmessages = nmessages + 1
      VLIDORT_Out%Status%TS_2OSCORR_CHECKMESSAGES(nmessages)     = Adjustl(Trim(messages))
      VLIDORT_Out%Status%TS_2OSCORR_ACTIONS      (nmessages)     = Adjustl(Trim(actions))
   endif

!  Exit if any failure

   if ( nmessages .ne. 0 ) then
      VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK = VLIDORT_SERIOUS
      VLIDORT_Out%Status%TS_2OSCORR_NCHECKMESSAGES    = nmessages
      return
   endif

!  Copy (from VLIDORT Type structures) or derive relevant inputs
!  =============================================================

!  flags

   do_plane_parallel = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL
   do_Lambertian     = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE
   do_deltam_scaling = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING

!  debug case

   if ( do_debug_23oct14 ) goto 455

!  Numbers

   nlayers  = VLIDORT_FixIn%Cont%TS_NLAYERS
   nstokes  = VLIDORT_FixIn%Cont%TS_NSTOKES
   nstreams = VLIDORT_FixIn%Cont%TS_NSTREAMS
   ngeoms   = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS

!  Linearization control

   do_LS_Jacobians = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
   do_LP_Jacobians = VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
   do n = 1, nlayers
      n1 = nlayers - n  + 1
      npars_vli(n)  = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n)
      npars(n1) = npars_vli(n)
   enddo
   nspars = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

!  some settings here can be set simply (default)
!        or by examining GreekMat(1) non-zero entries

   nfoumax  = 2*nstreams - 1
   nmoments = 2*nstreams - 1
   ncoeffs  = nmoments

!  Gaussian qudrature is re-done

   call Gaussquad_2os (MAXSTREAMS,nstreams,zero,one,qstreams,qweights)

!  Accuracy and albedo

   convergence_eps = VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY
   albedo          = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO

!  Geometries

   do L = 1, ngeoms
      do i = 1, 3
         geoms(i,1,L) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(L,i)
         geoms(i,2,L) = geoms(i,1,L) * deg_to_rad
         geoms(i,3,L) = cos(geoms(i,2,L))
         geoms(i,4,L) = sin(geoms(i,2,L))
      enddo
   enddo

!  Chapman factors. Not reversed (top down)

   SunChapman(:,:,1:ngeoms) = VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS(:,:,1:ngeoms)

!  Bulk optical properties. Reversed (bottom up)

   opdeps = zero ; omegas = zero ; truncfac = zero ; omegas_1 = zero ; omegas_2 = zero
   if ( do_deltam_scaling ) then
     do n = 1, nlayers
       n1 = nlayers + 1 - n
       opdeps(n1)   = VLIDORT_Out%Diag%TS_DELTAU_VERT(n)
       omegas(n1)   = VLIDORT_Out%Diag%TS_OMEGA_TOTAL(n)
       truncfac(n1) = VLIDORT_Out%Diag%TS_TRUNC_FACTOR(n)
       omegas_1(n1) = omegas(n1)
       omegas_2(n1) = omegas(n1) * truncfac(n1)       ! Not going to be required, I think.
     enddo
   else
     do n = 1, nlayers
       n1 = nlayers + 1 - n
       opdeps(n1)   = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)
       omegas(n1)   = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
       omegas_1(n1) = omegas(n1) 
     enddo
   endif

!  Linearized bulk optical properties
!     Single normalized quantities for 2OS

   LA_opdeps = zero ; LA_omegas = zero ; LA_truncfac = zero ; LA_omegas_1 = zero ; LA_omegas_2 = zero
   if ( do_LP_Jacobians ) then
     if ( do_deltam_scaling ) then
       do n = 1, nlayers
         n1 = nlayers + 1 - n
         do p = 1, npars_vli(n)
            p1 = npars(n1)
            LA_opdeps(n1,p1)   = opdeps(n1) * VLIDORT_LinOut%Diag%TS_L_DELTAU_VERT(p,n)
            LA_omegas(n1,p1)   = omegas(n1) * VLIDORT_LinOut%Diag%TS_L_OMEGA_TOTAL(p,n)
            LA_truncfac(n1,p1) = VLIDORT_LinOut%Diag%TS_L_TRUNC_FACTOR(p,n) ! Already single normalized
            LA_omegas_1(n1,p1) = LA_omegas(n1,p1) 
            LA_omegas_2(n1,p1) = LA_omegas(n1,p1) * truncfac(n1) + omegas(n1) * LA_truncfac(n1,p1)  ! Not needed, probably.
         enddo
       enddo
     else
       do n = 1, nlayers
         n1 = nlayers + 1 - n
         do p = 1, npars_vli(n)
            p1 = npars(n1)
            LA_opdeps(n1,p1)   = opdeps(n1) * VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(p,n)
            LA_omegas(n1,p1)   = omegas(n1) * VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(p,n)
            LA_omegas_1(n1,p1) = LA_omegas(n1,1) 
         enddo
       enddo
     endif
   endif

!  Scattering matrix optical properties
!   TAKE CAREFUL NOTE of the INDEXING SYSTEM for COEFFS, and the sign conventions

   coeffs = zero
   if ( do_deltam_scaling ) then
     do n = 1, nlayers
       n1 = nlayers + 1 - n
       coeffs(0:maxmoments,n1,1) =   VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,1)
       coeffs(0:maxmoments,n1,5) = - VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,2)
       coeffs(0:maxmoments,n1,2) =   VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,6)
       coeffs(0:maxmoments,n1,3) =   VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,11)
       coeffs(0:maxmoments,n1,6) = - VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,15)
       coeffs(0:maxmoments,n1,4) =   VLIDORT_Out%Diag%TS_GREEKMAT_TOTAL(0:MAXMOMENTS,n,16)
     enddo
   else
     do n = 1, nlayers
       n1 = nlayers + 1 - n
       coeffs(0:maxmoments,n1,1) =   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,1)
       coeffs(0:maxmoments,n1,5) = - VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,2)
       coeffs(0:maxmoments,n1,2) =   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,6)
       coeffs(0:maxmoments,n1,3) =   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,11)
       coeffs(0:maxmoments,n1,6) = - VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,15)
       coeffs(0:maxmoments,n1,4) =   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,n,16)
     enddo
   endif

!  Linearized optical properties - Single normalized quantities for 2OS
!   TAKE CAREFUL NOTE of the INDEXING SYSTEM for COEFFS

   LA_coeffs = zero
   if ( do_LP_Jacobians ) then
     if ( do_deltam_scaling ) then
       do n = 1, nlayers
         n1 = nlayers + 1 - n
         do p = 1, npars_vli(n)
           p1 = npars(n1)
           do m = 0, MAXMOMENTS
             LA_coeffs(m,n1,1,p1) = coeffs(m,n1,1) * VLIDORT_LinOut%Diag%TS_L_GREEKMAT_TOTAL(p,m,n,1)
             LA_coeffs(m,n1,5,p1) = coeffs(m,n1,5) * VLIDORT_LinOut%Diag%TS_L_GREEKMAT_TOTAL(p,m,n,2)
             LA_coeffs(m,n1,2,p1) = coeffs(m,n1,2) * VLIDORT_LinOut%diag%TS_L_GREEKMAT_TOTAL(p,m,n,6)
             LA_coeffs(m,n1,3,p1) = coeffs(m,n1,3) * VLIDORT_LinOut%Diag%TS_L_GREEKMAT_TOTAL(p,m,n,11)
             LA_coeffs(m,n1,6,p1) = coeffs(m,n1,6) * VLIDORT_LinOut%Diag%TS_L_GREEKMAT_TOTAL(p,m,n,15)
             LA_coeffs(m,n1,4,p1) = coeffs(m,n1,4) * VLIDORT_LinOut%Diag%TS_L_GREEKMAT_TOTAL(p,m,n,16)
           enddo
         enddo
       enddo
     else
       do n = 1, nlayers
         n1 = nlayers + 1 - n
         do p = 1, npars_vli(n)
           p1 = npars(n1)
           do m = 0, MAXMOMENTS
             LA_coeffs(m,n1,1,p1) = coeffs(m,n1,1) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,1)
             LA_coeffs(m,n1,5,p1) = coeffs(m,n1,5) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,2)
             LA_coeffs(m,n1,2,p1) = coeffs(m,n1,2) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,6)
             LA_coeffs(m,n1,3,p1) = coeffs(m,n1,3) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,11)
             LA_coeffs(m,n1,6,p1) = coeffs(m,n1,6) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,15)
             LA_coeffs(m,n1,4,p1) = coeffs(m,n1,4) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(p,m,n,16)
           enddo
         enddo
       enddo
     endif
   endif

!  Surface stuff: albedo. [BRDF stuff is set within the Fourier Loop]

   albedo          = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO

!  First-order saved fields. Only required for deltam scaling - Do not now think this.
!   Status 10/24/14. Think we will not need this.

   FO_Zmatrix = zero ; FO_R1saved = zero
   if ( do_deltam_scaling ) then
      do L = 1, ngeoms
         FO_R1saved(0,1:nstokes,L) = VLIDORT_Sup%SS%TS_FO_R1SAVED(L,0,1:nstokes)
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            FO_Zmatrix(n1,1:nstokes,L) = VLIDORT_Sup%SS%TS_FO_ZMATRIX(L,n,1:nstokes,1)
            FO_R1saved(n1,1:nstokes,L) = VLIDORT_Sup%SS%TS_FO_R1SAVED(L,n,1:nstokes)
         enddo
      enddo
   endif

!  Skip this section

   go to 456

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        START OF DEBUG SECTION INPUT 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Debug section
!  =============

!  Continuation point

455 continue

   do_LP_Jacobians = .true.  ; npars  = 1
   do_LS_Jacobians = .false. ; nspars = 1

!  Get the Dumped data input

   if ( do_deltam_scaling ) then
      open(87,file='debugfiles/fort.87_DM',status='old')
   else
      open(87,file='debugfiles/fort.87_Rg',status='old')
   endif
   read(87,*)nlayers,nstreams,nstokes,nmoments,nfoumax,regular_ps,enhanced_ps
   read(87,*)convergence_eps,phi,xs,xv
   do k = 1, nstreams
      read(87,*)kd,qstreams(k),qweights(k)
   enddo
   do n = 1, nlayers
      read(87,*)omegas_r(n),omegas_2(n),opdeps_r(n),ncoeffs(n)
      if ( n.eq.31 ) omegas_r(n) = omegas_r(n-1)
      if ( n.eq.12 ) ray1 = omegas_r(n)*opdeps_r(n)
      if ( n.eq.32 ) ray2 = omegas_r(n)*opdeps_r(n)
      do l = 0, 2*nstreams-1
         read(87,*)nd,ld,(coeffs(l,n,k),k=1,6)
      enddo
   enddo
   close(87)

!  gas fiddle

   aa = 0.95
   do n = 1, nlayers
     if ( ncoeffs(n).gt.2.and.n.gt.12.and.n.lt.32 ) then
        dn = dble(n) ; molsca(n)  = ( (dn-12.0)*ray2 + (32.0-dn)*ray1 ) / 20.0
        z = opdeps_r(n) * omegas_r(n) - molsca(n) 
        aerext(n) = z / aa
        molabs(n) = opdeps_r(n) - molsca(n)  - aerext(n) 
      else
        molabs(n) = (one-omegas_r(n))*opdeps_r(n)
        molsca(n) = omegas_r(n)*opdeps_R(n)
        aerext(n) = zero
      endif
   enddo

!  Optical properties and linearization for the Jacobian

   do n = 1, nlayers
      opdeps(n)      =  molabs(n) + molsca(n) + aerext(n)
      omegas_1(n)    =  ( molsca(n) + aa * aerext(n) ) / opdeps(n)
      LA_opdeps(n,1)   = molabs(n)
      LA_omegas_1(n,1) = - omegas_1(n) * LA_opdeps(n,1) / opdeps(n)
   enddo
   LA_omegas_2 = zero ;  LA_coeffs   = zero

!   Get the dumped chapman factos

   SunChapman = zero
   open(88,file='debugfiles/fort.88',status='old')
       do n = 1, nlayers
         do k = 1, n
           read(88,*)SunChapman(n,k,1)
         enddo
       enddo
   close(88)

!  Get the saved FO results dimp (only for testing the deltam-code)

   open(89,file='debugfiles/fort.89_DM',status='old')
      read(89,*)(FO_R1saved(0,k,1),k=1,nstokes)
      do n = 1, nlayers
         read(89,*)(FO_Zmatrix(n,k,1),k=1,nstokes)
         read(89,*)(FO_R1saved(n,k,1),k=1,nstokes)
      enddo
   close(89)

!  Set albedo

   Albedo = 0.10_fpk

!  Set geometry

   ngeoms = 1  ; dtr = acos(-1.0d0)/180.0d0
   geoms(1,1,1) = acos(xs)/dtr 
   geoms(1,2,1) = acos(xs) 
   geoms(1,3,1) = xs
   geoms(1,4,1) = sqrt(one-xs*xs)
   geoms(2,1,1) = acos(xv)/dtr 
   geoms(2,2,1) = acos(xv) 
   geoms(2,3,1) = xv
   geoms(2,4,1) = sqrt(one-xv*xv)
   geoms(3,1,1) = phi 
   geoms(3,2,1) = phi*dtr 
   geoms(3,3,1) = cos(phi*dtr)
   geoms(3,4,1) = sin(phi*dtr)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        END OF DEBUG SECTION INPUT 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point for avoiding debug input

456 continue

!  Start Calculation
!  =================

!  Initialize major output

   R2    = zero ; Icorr    = zero
   LP_R2 = zero ; LP_Icorr = zero
   LS_R2 = zero ; LS_Icorr = zero

!  Calculate average secants for the solar transmittances

   if ( do_LP_Jacobians ) then
      call Set_avsecant_LPPlus &
         ( MaxGeoms, do_plane_parallel,do_LP_jacobians,nlayers,ngeoms,npars,geoms,&
           SunChapman,opdeps,LA_opdeps,avg_secants,LP_avg_secants )
   else
      call Set_avsecant &
        ( MaxGeoms, do_plane_parallel, nlayers, ngeoms, &
          geoms, SunChapman, opdeps, avg_secants )
   endif

!  debug
!   do n = 26,60
!      write(*,*)avg_secants(n,1), LA_opdeps(n,1),LP_avg_secants(n,44,1,1)
!   enddo
!   pause

!  calculate Multipliers for First- and Second-order calculations
!  --------------------------------------------------------------

   if ( do_LP_Jacobians ) then
      call Calculate_Multipliers_LPPlus &
        ( MaxGeoms, nstreams, nlayers, ngeoms, npars, do_LP_Jacobians, & ! Inputs
          qstreams, avg_secants, geoms, opdeps, omegas_1,              & ! Inputs
          LP_avg_secants, LA_opdeps, LA_omegas_1,                      & ! Inputs
          O1_QVTrans, O1_QVMult, O1_QATrans, O1_QAMult,                & ! First-order output
          LP_O1_QVTrans, LP_O1_QVMult, LP_O1_QATrans, LP_O1_QAMult,    & ! First-order output
          O2_AVTrans, O2_AVMult, O2_QVMult_d, O2_QAMult_d,             & ! Second-order output
          O2_QAVMult_du, O2_QAVMult_dd,                                & ! Second-order output
          LP_O2_AVTrans, LP_O2_AVMult, LP_O2_QVMult_d, LP_O2_QAMult_d, & ! Second-order output
          LP_O2_QAVMult_du, LP_O2_QAVMult_dd  )                          ! Second-order output
   else
      call Calculate_Multipliers &
         (MaxGeoms, nstreams, nlayers, ngeoms,            & ! Inputs
          qstreams, avg_secants, geoms, opdeps, omegas_1, & ! Inputs
          O1_QVTrans, O1_QVMult, O1_QATrans, O1_QAMult,   & ! First-order output
          O2_AVTrans, O2_AVMult, O2_QVMult_d,O2_QAMult_d, & ! Second-order output
          O2_QAVMult_du, O2_QAVMult_dd  )                   ! Second-order output
   endif

!  Debug

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   do n = 1, nlayers
!       do L = 1, nstreams
!          write(67,*)O1_QVTrans(n,L,1),O1_QVMult(n,L,1), O1_QATrans(n,L,1), O1_QAMult(n,L,1)
!      enddo
!   enddo
!   do n = 1, nlayers
!       do L = 1, nstreams
!          write(68,*)LP_O1_QVTrans(n,L,1,1),LP_O1_QVMult(n,L,1,1), LP_O1_QATrans(n,L,23,1,1), LP_O1_QAMult(n,L,23,1,1)
!      enddo
!   enddo
!   do n = 1, nlayers
!       do L = 1, nstreams
!          write(67,*)O2_AVTrans(n,1),O2_QVMult_d(n,L,1), &
!           O2_QAMult_d(n,L,1),O2_QAVMult_du(n,L,1),O2_QAVMult_dd(n,L,1)
!      enddo
!   enddo
!   do n = 1, nlayers
!       do L = 1, nstreams
!          write(68,*)LP_O2_AVTrans(n,23,1,1),LP_O2_QVMult_d(n,L,23,1,1),&
!           LP_O2_QAMult_d(n,L,23,1,1),LP_O2_QAVMult_du(n,L,23,1,1),LP_O2_QAVMult_dd(n,L,23,1,1)
!      enddo
!   enddo
!   pause'Multiplier output'
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!     Setups for the thread-safe alternative

    Tqroot6 = -quarter*sqrt(six)
    Trl(0) = zero ; Trlsq(0) = zero ; Tr2lp1(0) = one
    do l = 1, MAXMOMENTS
       ll = l + MAXMOMENTS
       Trl(l)  = real(l,fpk) ; Trlsq(l) = Trl(l) * Trl(l)
       Trl(ll) = real(ll,fpk)
       Tr2lp1(l) = two*Trl(l) + one
    enddo
    Tns1 = nstreams + 1 ; Tns2 = nstreams + 2
    Txmu(1:nstreams) = qstreams(1:nstreams)

!  Start geometry loop
!  -------------------

   do L = 1, ngeoms

      if (verbo) print *,' adding: treat Geometry ', L

!  Perform Fourier loop
!  --------------------

!  Initialize Fourier loop

      m = -1
      nextm = .true.

!  Cosines and secants

      mus = geoms(1,3,L) ; xs  = one / mus
      muv = geoms(2,3,L) ; xv  = one / muv
      phi = geoms(3,1,L)

!  All cosines, sines, Initial Spherical functions
!     Setups for the thread-safe alternative

      Txmu(Tns1) = geoms(1,3,L)
      Txmu(Tns2) = geoms(2,3,L)
      do i = 1, Tns2
        help = one - Txmu(i) * Txmu(i)
        Txsi(i) = sqrt(abs(help))
        xmu_m1 = one - Txmu(i) ; Txmu_m1_sq(i) = xmu_m1 * xmu_m1
        xmu_p1 = one + Txmu(i) ; Txmu_p1_sq(i) = xmu_p1 * xmu_p1
        TPlm_2all(i,1) =  Tqroot6 * Txsi(i) * Txsi(i)
        TPlm_2all(i,2) = - half * Txsi(i)* xmu_m1
        TPlm_2all(i,3) =   half * Txsi(i)* xmu_p1
        TPlm_2all(i,4) = - quarter * Txmu_m1_sq(i)
        TPlm_2all(i,5) = - quarter * Txmu_p1_sq(i)
      enddo

!  Open loop

      do while (nextm)

!  increase Fourier number by 1

         m = m+1

!  initialise the Fourier components for R1 and R2 (first/second order). Very Important

         R2c    = zero ; R2s    = zero ; R2cscal    = zero
         LP_R2c = zero ; LP_R2s = zero ; LP_R2cscal = zero
         LS_R2c = zero ; LS_R2s = zero ; LS_R2cscal = zero
         R1c    = zero ; R1s    = zero ; R1cscal    = zero
         LP_R1c = zero ; LP_R1s = zero ; LP_R1cscal = zero
         LS_R1c = zero ; LS_R1s = zero ; LS_R1cscal = zero

!  setup the Fourier comonent source terms at BOA for the R1 (first order)
!    These are pre-initialized. Only the cosine terms for m = 0 survive for Lambertian
!    BRDF case: Copy and re-order the BRDF supplement inputs to VLIDORT
!     Old code: Glint surface (Cox-Munk, Giss). Surface Type # 2
!     Old code: Land surface BRDFs.             Surface Type # 3/4

         if (do_Lambertian) then
            if ( m.eq.0 ) then
               R1cscal(1:2,1:nstreams)  = Albedo
               R1c(1:2,1:nstreams,1,1)  = Albedo
               if ( do_LS_Jacobians ) then
                  Ls_R1cscal(1:2,1:nstreams,1:nspars)     = one
                  Ls_R1c    (1:2,1:nstreams,1,1,1:nspars) = one
               endif
            endif
         else
            R1cscal(1,1:nstreams) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,1,1:nstreams,L)
            R1cscal(2,1:nstreams) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F(m,1,1:nstreams,L)
            R1c(1,1:nstreams,1,1) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,1,1:nstreams,L)
            R1c(1,1:nstreams,1,2) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,2,1:nstreams,L)
            R1c(1,1:nstreams,2,1) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,5,1:nstreams,L)
            R1c(1,1:nstreams,2,2) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,6,1:nstreams,L)
            R1s(1,1:nstreams,1,3) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,3,1:nstreams,L)
            R1s(1,1:nstreams,2,3) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,7,1:nstreams,L)
            R1s(1,1:nstreams,3,1) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,9,1:nstreams,L)
            R1s(1,1:nstreams,3,2) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,10,1:nstreams,L)
            R1s(1,1:nstreams,3,3) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,11,1:nstreams,L)
            R1s(1,1:nstreams,4,4) = VLIDORT_Sup%BRDF%TS_BRDF_F_0   (m,16,1:nstreams,L)
            R1c(2,1:nstreams,1,1) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,1,1:nstreams,L)
            R1c(2,1:nstreams,1,2) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,2,1:nstreams,L)
            R1c(2,1:nstreams,2,1) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,5,1:nstreams,L)
            R1c(2,1:nstreams,2,2) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,6,1:nstreams,L)
            R1s(2,1:nstreams,1,3) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,3,1:nstreams,L)
            R1s(2,1:nstreams,2,3) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,7,1:nstreams,L)
            R1s(2,1:nstreams,3,1) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,9,1:nstreams,L)
            R1s(2,1:nstreams,3,2) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,10,1:nstreams,L)
            R1s(2,1:nstreams,3,3) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,11,1:nstreams,L)
            R1s(2,1:nstreams,4,4) = VLIDORT_Sup%BRDF%TS_USER_BRDF_F   (m,16,1:nstreams,L)

            if ( do_LS_Jacobians ) then
              LS_R1cscal(1,1:nstreams,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,1,1:nstreams,L)
              LS_R1cscal(2,1:nstreams,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F (1:nspars,m,1,1:nstreams,L)
              LS_R1c(1,1:nstreams,1,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,1,1:nstreams,L)
              LS_R1c(1,1:nstreams,1,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,2,1:nstreams,L)
              LS_R1c(1,1:nstreams,2,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,5,1:nstreams,L)
              LS_R1c(1,1:nstreams,2,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,6,1:nstreams,L)
              LS_R1s(1,1:nstreams,1,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,3,1:nstreams,L)
              LS_R1s(1,1:nstreams,2,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,7,1:nstreams,L)
              LS_R1s(1,1:nstreams,3,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,9,1:nstreams,L)
              LS_R1s(1,1:nstreams,3,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,10,1:nstreams,L)
              LS_R1s(1,1:nstreams,4,4,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0    (1:nspars,m,16,1:nstreams,L)
              LS_R1s(1,1:nstreams,3,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,11,1:nstreams,L)
              LS_R1c(2,1:nstreams,1,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,1,1:nstreams,L)
              LS_R1c(2,1:nstreams,1,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,2,1:nstreams,L)
              LS_R1c(2,1:nstreams,2,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,5,1:nstreams,L)
              LS_R1c(2,1:nstreams,2,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,6,1:nstreams,L)
              LS_R1s(2,1:nstreams,1,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,3,1:nstreams,L)
              LS_R1s(2,1:nstreams,2,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,7,1:nstreams,L)
              LS_R1s(2,1:nstreams,3,1,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,9,1:nstreams,L)
              LS_R1s(2,1:nstreams,3,2,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,10,1:nstreams,L)
              LS_R1s(2,1:nstreams,3,3,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,11,1:nstreams,L)
              LS_R1s(2,1:nstreams,4,4,1:nspars) = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:nspars,m,16,1:nstreams,L)
            endif

         endif

!  Legendre setups for the Zmatrix thread-safe alternative

         if (m.eq.0) Ttwom = one ; if (m.gt.0) Ttwom = Ttwom * half
         msq = m * m ; Trmsq = real(msq,fpk)
         bin2mm = one
         do n = 1, m
            bin2mm = bin2mm *Trl(n+m)/Trl(n)
         enddo
         Tbinfac = Ttwom * sqrt(bin2mm) ; binf2 = zero ; bin2m2 = zero
         if (m.ge.2) then
            bin2m2 = bin2mm * Trl(m) * Trl(m-1) / ( Trl(m+1)*Trl(m+2) )
            binf2  = -Ttwom * sqrt(bin2m2)
         endif
         if ( m.gt.2 ) then
            do i = 1, Tns2
               urootm = Txsi(i)**(m-2)
               TPlm_mlm(i) = binf2 * urootm * Txmu_m1_sq(i)
               TPlm_mlp(i) = binf2 * urootm * Txmu_p1_sq(i)
            enddo
         endif
         Tsqlm(m) = zero ; Tsqlm(m+1:nmoments) = sqrt(Trlsq(m+1:nmoments)-Trmsq)
         Tsql4(2) = zero ; Tsql4(3:nmoments)   = sqrt(Trlsq(3:nmoments)  -four)

!  here is the Recursion loop
!  --------------------------

         do n = 1, nlayers

            if (verbo) print *,' adding: treat Fourier term m = ',m
            if (verbo) print *,' adding: treat layer ', N

!  Average secants

            xa = avg_secants(n,L)
            if ( do_LP_Jacobians ) then
               do k = 1, nlayers
                  LP_xa(k,1:npars(k)) = LP_avg_secants(n,k,1:npars(k),L)
               enddo
            endif

!  Other shorthand

            W1 = omegas_1(n) ; W2 = omegas_2(n)

!  local layer phase matrices for single Fourier component

!            ndcoeffs = min(nmoments,ncoeffs(n))
!            if ( do_LP_Jacobians ) then
!               call Calculate_Zmatrix_LAPlus &
!                  ( MaxGeoms,m,N,L,do_LP_Jacobians,            & ! Inputs
!                    nstreams,nstokes,ndcoeffs,npars(n),        & ! Inputs
!                    qstreams, geoms, coeffs(0:MAXMOMENTS,N,:), & ! Inputs
!                    LA_coeffs(0:MAXMOMENTS,N,:,:),             & ! Inputs
!                    Ptc,Pts,Prc,Prs,LP_Ptc,LP_Pts,LP_Prc,LP_Prs )
!            else
!               call Calculate_Zmatrix &
!                  ( MaxGeoms,m,N,L,nstokes,nstreams,ndcoeffs,   & !I
!                    qstreams, geoms, coeffs(0:MAXMOMENTS,N,:),  & !I
!                    Ptc,Pts,Prc,Prs)
!            endif

!  local layer phase matrices for single Fourier component.
!         ALTERNATIVE Thread-Safe Treatment

            ndcoeffs = min(nmoments,ncoeffs(n))
            if ( do_LP_Jacobians ) then
               call Calculate_Zmatrix_LAPlus_ALT &
                ( do_LP_Jacobians, m, nstokes, nstreams, Tns1, Tns2,       & ! Inputs
                  ndcoeffs, npars(n), Trl, Tr2Lp1, Tsql4, Tsqlm, Tbinfac,     & ! Inputs
                  coeffs(0:MAXMOMENTS,N,:), LA_coeffs(0:MAXMOMENTS,N,:,:), & ! Inputs
                  Txmu, Txsi, TPlm_2all, TPlm_mlp, TPlm_mlm,               & ! Inputs
                  Ptc, Pts, Prc, Prs, LP_Ptc, LP_Pts, LP_Prc, LP_Prs )
            else
               call Calculate_Zmatrix_ALT &
                ( m, nstokes, nstreams, Tns1, Tns2, ndcoeffs, Trl, Tr2Lp1, & ! Inputs
                  Tsql4, Tsqlm, Tbinfac, coeffs(0:MAXMOMENTS,N,:),         & ! Inputs
                  Txmu, Txsi, TPlm_2all, TPlm_mlp, TPlm_mlm,               & ! Inputs
                  Ptc,Pts,Prc,Prs )
            endif

!             if ( m.eq.0.and.n.eq.30 ) then
!             do k = 1, nstreams
!               write(93,*)k,Ptc(1,k,1,1),Ptc(1,k,2,1),Ptc(1,k,1,3),Ptc(1,k,2,3)
!             enddo
!             do k = 1, nstreams
!               write(93,*)k,Prc(2,k,1,1),Prc(2,k,2,1),Ptc(2,k,3,1),Ptc(2,k,4,1)
!             enddo
!             pause 'zmat 1'
!             endif

!  Second order calculation

            if ( do_LP_Jacobians .or. do_LS_Jacobians ) then
               call Calculate_SecondOrder_LPSPlus &
                ( do_LP_Jacobians,do_LS_Jacobians,m,n,nstreams,nstokes,nlayers,npars,nspars,  & ! Input Control
                  qweights,omegas_1(n),xv,xa,   LA_omegas_1(n,:), LP_xa,         & ! Input Basics
                  O2_AVTrans(n,L),              LP_O2_AVTrans(n,:,:,L),          & ! Input R2 Multiplier
                  O2_QVMult_d(n,:,L),           LP_O2_QVMult_d(n,:,:,:,L),       & ! Input R2 Multipliers 
                  O2_QAMult_d(n,:,L),           LP_O2_QAMult_d(n,:,:,:,L),       & ! Input R2 Multipliers
                  O2_QAVMult_du(n,:,L),         LP_O2_QAVMult_du(n,:,:,:,L),     & ! Input R2 Multipliers 
                  O2_QAVMult_dd(n,:,L),         LP_O2_QAVMult_dd(n,:,:,:,L),     & ! Input R2 Multipliers 
                  Ptc,Pts,Prc,Prs,              LP_Ptc,LP_Pts,LP_Prc,LP_Prs,     & ! Input scatt matrices
                  R1c,R1s,R1cscal,LP_R1c,LP_R1s,LP_R1cscal,LS_R1c,LS_R1s,LS_R1cscal,  & ! Input R1 sources
                  R2c,R2s,R2cscal,LP_R2c,LP_R2s,LP_R2cscal,LS_R2c,LS_R2s,LS_R2cscal )   ! Output R2 reflectances
            else
               call Calculate_SecondOrder &
                ( m, n, nstreams, nstokes,                    & ! Input Control
                  qweights, omegas_1(n), xv, xa,              & ! Input local
                  Ptc,Pts,Prc,Prs,                            & ! Input Local scattering matrices
                  O2_AVTrans(n,L),                            & ! Input R2 Multiplier
                  O2_QVMult_d  (n,:,L), O2_QAMult_d  (n,:,L), & ! Input R2 Multipliers
                  O2_QAVMult_du(n,:,L), O2_QAVMult_dd(n,:,L), & ! Input R2 Multipliers 
                  R1c,R1s,R1cscal,                            & ! Input R1 source field
                  R2c,R2s,R2cscal )                             ! Output R2 reflectance
            endif

!  Add the direct term from delta-M scaling
!   Only add this for Fourier = 0, and then only to R2cscal and R2c
!   Debugged, 18 December 2012

            if ( do_deltam_scaling .and. m.eq.0 ) then
               if ( do_LP_Jacobians ) then
!  PLACEHOLDER
               else
                  call Scale_SecondOrder &
                    (MaxGeoms,nstokes,n,L,xv,xa,O2_AVTrans(n,L),  & ! I
                     omegas_1(n),omegas_2(n), opdeps(n),          & ! I
                     FO_Zmatrix, FO_R1saved,                      & ! I
                     R2c,R2cscal )                                  ! Modified
               endif
            endif

!  debug 3/4/15
!       if(m.eq.0)write(*,*)m,n,R1cscal(1:2,7),LS_R1cscal(1:2,7,1)
!       if(m.eq.0)write(*,*)m,n,R2c(1:2),LS_R2c(1:2,1)
!       if(m.eq.0)write(*,*)m,n,R1c(1:2,7,1,1),LS_R1c(1:2,7,1,1,1)
!       if(m.eq.0)write(*,*)m,n,LP_R1cscal(1:2,7,25,1)
!       if(m.eq.0)write(*,*)m,layer,R2c(1:2),LS_R2c(1:2,1)
!       if (m.eq.0.and.n.eq.nlayers)pause'gronk'

!  Get the first order of scattering (using previous layer)

            if ( do_LP_Jacobians .or. do_LS_Jacobians ) then
               call Calculate_FirstOrder_LPSPlus &
                ( do_LP_Jacobians, do_LS_Jacobians,m,n,nstreams,nlayers,npars,nspars, & ! Input Control
                  Prc,Prs,            LP_Prc,LP_Prs,                 & ! Scattering input
                  O1_QVTrans(n,:,L),  LP_O1_QVTrans(n,:,:,L),        & ! First-order multiplier input
                  O1_QVMult(n,:,L),   LP_O1_QVMult(n,:,:,L),         & ! First-order multiplier input
                  O1_QATrans(n,:,L),  LP_O1_QATrans(n,:,:,:,L),      & ! First-order multiplier input
                  O1_QAMult(n,:,L),   LP_O1_QAMult(n,:,:,:,L),       & ! First-order multiplier input
                  R1c,R1s,R1cscal,LP_R1c,LP_R1s,LP_R1cscal,LS_R1c,LS_R1s,LS_R1cscal ) ! First-order output

            else
               call Calculate_FirstOrder &
                ( m, n, nstreams, Prc, Prs,            & ! Control + scattering   input
                  O1_QVTrans(n,:,L), O1_QVMult(n,:,L), & ! First-order multiplier input
                  O1_QATrans(n,:,L), O1_QAMult(n,:,L), & ! First-order multiplier input
                  R1c, R1s, R1cscal )                    ! First-order output
            endif

!        if(m.eq.0.and.do_LP_Jacobians)      write(34,*)n,LP_R2cscal(17,1),R2cscal
!        if(m.eq.0.and..not.do_LP_Jacobians) write(35,*)n,R2cscal
!       if(m.eq.0)write(*,*)m,layer,R2c(1:2),LS_R2c(1:2,1)
!       if (m.eq.0.and.n.eq.nlayers)pause'gronk'

!  end recursion

         enddo

!  Add fourier component

         if ( do_LP_Jacobians ) then
            call Add_Fourier_Component_LPSPlus &
             ( m, nstokes, phi, do_lp_jacobians, do_ls_jacobians, nlayers, npars, nspars, & ! Inputs
               R2c,R2s,R2cscal,LP_R2c,LP_R2s,LP_R2cscal,LS_R2c,LS_R2s,LS_R2cscal,         & ! Inputs
               R2(:,L), Icorr(L), LP_R2(:,:,:,L), LP_Icorr(:,:,L), LS_R2(:,:,L), LS_Icorr(:,L) )
         else
            call Add_Fourier_Component &
             ( m, nstokes, phi, R2c, R2s, R2cscal, & ! Inputs      
               R2(:,L), Icorr(L) )
         endif

!  debug 3/4/15
!         write(*,*)M,LS_R2(1:2,1,L), LS_Icorr(1,L) 

!  Test Fourier convergence. (Use VZA and SZA cosines = muv, mus)

!mick fix 4/21/2015 - added "almost" to argument list
         call Test_Fourier_Convergence &
             (verbo,m,nfoumax,nstokes,convergence_eps,muv,mus,almost,nextm,R2s,R2c)

!  End fourier loop

      enddo

!  End answers (multiply by SZA cosine)

      R2(1:nstokes,L) = R2(1:nstokes,L) * mus
      Icorr(L)        = Icorr(L)        * mus
      if ( do_LP_Jacobians ) then
         do k = 1, nlayers 
            do p = 1, npars(k)
               LP_R2(1:nstokes,k,p,L) = LP_R2(1:nstokes,k,p,L) * mus
               LP_Icorr(k,p,L)        = LP_Icorr(k,p,L)        * mus
            enddo
         enddo
      endif
      if ( do_LS_Jacobians ) then
         LS_R2(1:nstokes,1:nspars,L) = LS_R2(1:nstokes,1:nspars,L) * mus
         LS_Icorr(1:nspars,L)        = LS_Icorr(1:nspars,L)        * mus
      endif

!  End Geometry Index loop

   enddo

!  Interpret results to VLIDORT output
!  ===================================

!  Copy to vlidort output. Note reverse ordering of Profile WFs
!    VLIDORT's 2OS output is the smae indexing as the MS and SS ordering

   do L = 1, ngeoms
      VLIDORT_Out%Main%TS_2OSCORR_R2(L,1:nstokes) = R2(1:nstokes,L)
      VLIDORT_Out%Main%TS_2OSCORR_ICORR(L)        = Icorr(L)
      if ( do_LP_Jacobians ) then
         do n = 1, nlayers
            n1 = nlayers - n + 1
            do p = 1, npars_vli(n)
               p1 = npars(n1)
               VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(p,n,L,1:nstokes) = LP_R2(1:nstokes,n1,p1,L)
               VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(p,n,L)        = LP_Icorr(n1,p1,L)
            enddo
         enddo
      endif
      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(p,L,1:nstokes) = LS_R2(1:nstokes,p,L)
            VLIDORT_LinOut%surf%TS_2OSCORR_LS_ICORR(p,L)        = LS_Icorr(L,p)
         enddo
      endif
   enddo

!  Apply 2OS correction to VLIDORT Radiance/Jacobian output
!  ========================================================

!   2OS results are MS-only, must be added to FO Radiance/Jacobian results
!   R2 and Icorr must be multiplied by FLUXFAC * PI/mu0 for Radiance/Jacobian output 

   piconv = acos( - 1.0d0 ) ; fluxfac = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Correct the Stokes Intensity output from VLIDORT (which includes first-order).
!  Compute the Q and U components from first and second order of scattering results.
!  Same thing for the Profile and Surface Jacobians, if flagged
!       Jacobian Indexing is now completely that of VLIDORT.

   do L = 1, ngeoms

!  multiplying factor. Should be mu0/pi, Vijay 11 May 2015

!      Conversion = fluxfac * piconv / geoms(1,3,L)
      Conversion = fluxfac * geoms(1,3,L) / piconv

!  Stokes output

      VLIDORT_Out%Main%TS_STOKES(1,L,1,upidx) = VLIDORT_Out%Main%TS_STOKES(1,L,1,upidx) &
                                              + VLIDORT_Out%Main%TS_2OSCORR_ICORR(L) * Conversion
      DO O1 = 2, nstokes
          VLIDORT_Out%Main%TS_STOKES(1,L,O1,upidx) = VLIDORT_Sup%SS%TS_STOKES_SS(1,L,o1,upidx) &
                                                   + VLIDORT_Sup%SS%TS_STOKES_DB(1,L,o1)       &
                                                   + VLIDORT_Out%Main%TS_2OSCORR_R2(L,O1) * Conversion
      enddo

!   Repeat process for Profile Jacobians

      if ( do_LP_Jacobians ) then
         do n = 1, nlayers
            do p = 1, npars_vli(n)
               VLIDORT_LinOut%Prof%TS_PROFILEWF(p,n,1,L,1,upidx) = &
                            VLIDORT_LinOut%Prof%TS_PROFILEWF(p,n,1,L,1,upidx) &
                          + VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(p,n,L) * Conversion
               do o1 = 2, nstokes
                 VLIDORT_LinOut%Prof%TS_PROFILEWF(p,n,1,L,o1,upidx) = &
                            VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(p,n,1,L,O1,upidx) &
                          + VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB(p,n,1,L,O1)       &
                          + VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(p,n,L,O1) * Conversion
               enddo
            enddo
         enddo
      endif

!  Repeat process for Surface Jacobians

      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            VLIDORT_LinOut%Surf%TS_SURFACEWF(p,1,L,1,upidx) = VLIDORT_LinOut%Surf%TS_SURFACEWF(p,1,L,1,upidx) &
                                                            + VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(p,L) * Conversion
            do o1 = 2, nstokes
              VLIDORT_LinOut%Surf%TS_SURFACEWF(p,1,L,O1,upidx) = VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(p,1,L,O1) &
                                                               + VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(p,L,O1) * Conversion
            enddo
         enddo
      endif

!  End Geometry loop

   enddo

!  finish

end subroutine VLIDORT_2OSCORR_LPS_MASTER

!  End module

END MODULE vlidort_2OScorr_lps_master_m

