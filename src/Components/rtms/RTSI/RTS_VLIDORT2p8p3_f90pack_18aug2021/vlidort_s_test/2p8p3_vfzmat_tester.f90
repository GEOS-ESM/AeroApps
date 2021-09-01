
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

      program Vfzmat_Tester

!  This is a new Version 2.8 driver. Created 9/19/16
!  This is designed to test the new VFZMAT supplement
!     R. Spurr. RT Solutions Inc.

!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  supplement Modules

      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

      USE VLIDORT_AUX_m,    Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m, Only : VLIDORT_INPUT_MASTER, VLIDORT_Sup_Init
      USE VLIDORT_MASTERS_m
 
      IMPLICIT NONE

!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT debug input control

      LOGICAL :: DO_DEBUG_INPUT

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  Help variables
!   Rob 3/22/15. renamed thread --> task

      INTEGER ::          I,IB,K,L,LDUM,N,N6,N61,NDUM,NTASKS,O1,T,TASK,V,UM,UTA
      DOUBLE PRECISION :: KD, WAER, TAER, PARCEL, RAYWT, AERWT, DEPOL, F11, F12, F33, F34
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT, MU
      LOGICAL ::          DO_STD_VFZMAT_TEST

      INTEGER ::          OutUnit

!  Proxies for the VFZMAT supplement
!  ---------------------------------

!  Flags (observational geometry, Sunlight), Geometry numbers
 
      LOGICAL ::          DO_OBSGEOMS, DO_SUNLIGHT, DO_UPWELLING, DO_DNWELLING
      INTEGER ::          N_GEOMS, N_SZAS, N_VZAS, N_AZMS
      INTEGER ::          OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Angles. Convention as for  VLIDORT

      DOUBLE PRECISION :: SZAS (MAX_SZANGLES)
      DOUBLE PRECISION :: VZAS (MAX_USER_VZANGLES)
      DOUBLE PRECISION :: AZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: OBSGEOMS (MAX_GEOMETRIES,3)

!  Input data from file
!  --------------------

!  The Rayleigh (molecular atmopshere)

      DOUBLE PRECISION :: MOLOMG ( MAXLAYERS ), MOLEXT ( MAXLAYERS )
      DOUBLE PRECISION :: RAYMOMS ( 0:2,MAXLAYERS )

!  The aerosol data:  FMatrix data for interpolation

      Integer, parameter :: Max_InAngles = 361!181
      integer            :: N_InAngles
      DOUBLE PRECISION   :: InAngles    ( Max_InAngles )
      DOUBLE PRECISION   :: InFmatrices ( MAXLAYERS, Max_InAngles, 6 )
      Logical            :: Exist_InFmatrices ( MAXLAYERS )

!  The VFZMAT supplemental variables for Rayleigh
!  Output Fmatrices (Calculated from Coefficients), Zmatrices, rayleigh coefficients
!     --- The Rayleigh coefficients are "PROBLEM_RAY"

      DOUBLE PRECISION :: RayFmatrices_up   ( MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION :: RayFmatrices_dn   ( MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION :: RayZmatrices_up   ( MAX_GEOMETRIES, 4, 4 )
      DOUBLE PRECISION :: RayZmatrices_dn   ( MAX_GEOMETRIES, 4, 4 )
      DOUBLE PRECISION :: RayCoeffs(0:2,6)

!  The VFZMAT supplemental variables for Aerosols
!  Output Fmatrices (Interpolated), Zmatrices, Fmatrix coefficients

      DOUBLE PRECISION :: OutFmatrices_up   ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION :: OutFmatrices_dn   ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION :: Zmatrices_up      ( MAX_GEOMETRIES, MAXLAYERS, 4, 4 )
      DOUBLE PRECISION :: Zmatrices_dn      ( MAX_GEOMETRIES, MAXLAYERS, 4, 4 )
      DOUBLE PRECISION :: FmatCoeffs        ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 6 )

!  VLIDORT Proxies
!  ---------------

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
      LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      INTEGER ::          NLAYERS, NSTOKES, NFINELAYERS

!  Height proxy, # moments

      INTEGER ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Optical proxies. Fmatrix proxies are new for Version 2.8

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 
      DOUBLE PRECISION :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 

!  Proxies (Map no longer required for Version 2.8)

!      INTEGER ::          UTAMAP(MAX_USER_LEVELS,10)
      INTEGER ::          N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  VLIDORT standard output preparation

      INTEGER ::          N_GEOMETRIES
      INTEGER ::          N_SZANGLES

!  Saved results
!  -------------

      INTEGER, PARAMETER :: MAX_TASKS = 10

      DOUBLE PRECISION :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAX_TASKS )

!  #####################################
!      START PROGRAM
!  #####################################

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_s_test/2p8p3_VLIDORT_ReadInput.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

!  really should use the F-matrix capability

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT = .true.

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( '2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

!  Initialize some input variables not handled by VLIDORT_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )

!  1/31/21. Version 2.8.3, New variables must be set by hand

      VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-20

!  1/31/21. Version 2.8.3, Set the DO_MSSTS flag to generate output for the MS sphericity corrections

      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .false.

!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.

!  1/31/21. Version 2.8.3, Set the NSTOKES2 FOURIER0 flag by hand

      VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2 = .true.

!  Set some proxies for use with the new supplement
!  ------------------------------------------------

!  General

      do_Sunlight  = .true.
      do_upwelling = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      do_dnwelling = VLIDORT_FixIn%Bool%TS_DO_DNWELLING

      nstokes      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers      = VLIDORT_FixIn%Cont%TS_nlayers
      n_szangles   = VLIDORT_ModIn%MSunrays%TS_n_szangles

      do_ObsGeoms  = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY

      offsets = 0 ; n_szas  = 0 ; n_vzas  = 0 ; n_azms  = 0 ;
      szas = zero ; vzas = zero ; azms = zero ; obsgeoms = zero

!  Angles for VFZMAT
 
      if ( do_ObsGeoms ) then
        n_geoms = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
        do v = 1, n_geoms
           obsgeoms(v,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(v,1:3)
        enddo
      else
        n_szas  = VLIDORT_ModIn%MSunrays%TS_n_szangles
        szas(1:n_szas) = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:n_szas)
        n_vzas  = VLIDORT_ModIn%MUserVal%TS_n_user_vzangles
        vzas(1:n_vzas) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:n_vzas)
        n_azms  = VLIDORT_ModIn%MUserVal%TS_n_user_relazms
        azms(1:n_azms) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:n_azms)
        n_geoms = n_szas * n_vzas * n_azms
        DO IB = 1, n_szas
          DO UM = 1, n_vzas
            offsets(ib,um) = n_azms * ( n_vzas * ( ib - 1 ) + um - 1 ) 
          enddo
        enddo
      endif

!  Set Output-level Proxies (saved values)

      n_user_levels                = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels(1:n_user_levels) = VLIDORT_ModIn%MUserVal%TS_user_levels(1:n_user_levels)

!  Version 2.7. Map is necessary because FO output is ONLY AT LEVEL BOUNDARIES
!  Version 2.8. Map no longer required, as FO code now has partials
!      utamap = 0
!      do uta = 1, n_user_levels
!         utamap(uta,1:4) = uta
!      enddo
!      utamap(1,5:8) = 1 ; utamap(2,5:8) = 2 ; utamap(5,5:8) = 3

!  Rob 3/22/15. Add wavelength Diagnostic (just a dummy here)
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Get the pre-prepared Molecular atmosphere
!  -----------------------------------------

      height_grid = zero
      open(45,file='vlidort_s_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n),molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)
      depol = ( 1.0d0- 2.0 * raymoms(2,1) ) / ( 1.0d0 + raymoms(2,1) )

!  Initialise optical proxies

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero
      Fmatrix_up           = zero
      Fmatrix_dn           = zero

!  Surface

      lambertian_albedo    = 0.05d0

!  Add Aerosols bottom 6 layers, spread evenly, with total AOD 0.5, SSAlbedo = 0.95

      n6 = nlayers - 6 ; n61 = n6 + 1
      waer = 0.95d0 ; taer = 0.5d0

!  F-matrix Inputs at 181 angles @ 1-degree separation

      InFmatrices = zero
      Exist_InFmatrices = .false.
      Exist_InFmatrices(n61:nlayers) = .true.

!  Get F-matrix Input from Mie. Same in all 6 lowest layers

      do_std_vfzmat_test = .true.
      !do_std_vfzmat_test = .false.
      if ( do_std_vfzmat_test ) then
        !Use original F-matrix file (with 181 angles @ 1-degree separation)
        N_InAngles = 181
        open(1,file='vlidort_s_test/RTSMie_InScatMat.Res', status='old')
        ngreek_moments_input = 80
      else
        !Use Deirmendjian cumulus cloud C.1 F-matrix file to compare with solar test results
        !(with 361 angles @ 0.5-degree separation)
        N_InAngles = 361
        open(1,file='vlidort_s_test/RTSMie_vfzmat_test.dat', status='old')
        !Skip header lines & other unneeded lines
        do i=1,389
          read(1,*)
        enddo
        ngreek_moments_input = 365
      endif

      do k = 1, N_InAngles
        read(1,*) InAngles(k),f11,f33,f12,f34
        !write(*,*) InAngles(k),f11,f33,f12,f34
        mu = cos(InAngles(k)*deg_to_rad) 
        InFmatrices(n61:nlayers,k,1) =  f11
        InFMatrices(n61:nlayers,k,2) =  f11
        InFMatrices(n61:nlayers,k,3) =  f33
        InFMatrices(n61:nlayers,k,4) =  f33
        InFMatrices(n61:nlayers,k,5) =  f12
        InFMatrices(n61:nlayers,k,6) = -f34
      enddo
      close(1)

!  Call to the supplement Rayleigh-master

      Call vfzmat_Rayleigh &
       ( max_geometries, max_szangles, max_user_vzangles, max_user_relazms, & ! Input  Dimensions (VLIDORT)
         do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,              & ! Input  Flags
         nstokes, n_geoms, n_szas, n_vzas, n_azms,                          & ! Input  Numbers
         offsets, deg_to_rad, szas, vzas, azms, obsgeoms, Depol,            & ! Input  Geometries + Depol
         RayFmatrices_up, RayFmatrices_dn, RayZmatrices_up, RayZmatrices_dn, RayCoeffs )

!  Call to the supplement aerosol master

      Call vfzmat_Master &
       ( maxmoments_input, max_geometries, max_szangles, max_user_vzangles,        & ! Input  Dimensions (VLIDORT)
         max_user_relazms, maxlayers,                                              & ! Input  Dimensions (VLIDORT)
         max_InAngles, N_InAngles, InAngles, InFmatrices, Exist_InFmatrices,       & ! Input  Fmatrices
         do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
         ngreek_moments_input, nlayers, nstokes, n_geoms, n_szas, n_vzas, n_azms,  & ! Input  Numbers
         offsets, deg_to_rad, szas, vzas, azms, obsgeoms,                          & ! Input  Geometries
         OutFmatrices_up, OutFmatrices_dn, Zmatrices_up, Zmatrices_dn, FMatCoeffs )

!  Set optical properties for the atmosphere:

!  Fix the non-aerosol layers

      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        greekmat_total_input(0:2,n,1) = RayCoeffs(0:2,1)
        fmatrix_up(n,1:n_geoms,1) = RayFmatrices_up(1:n_geoms,1)
        fmatrix_dn(n,1:n_geoms,1) = RayFmatrices_dn(1:n_geoms,1)
      enddo

!  Fix layers with aerosol

      parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
      do n = n61, nlayers
        aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
        aersca = aerext * waer      ; molsca = molomg(n) * molext(n)
        totext = molext(n) + aerext ; totsca = molsca    + aersca
        raywt  = molsca / totsca    ; aerwt  = aersca / totsca
        deltau_vert_input(n) = totext
        omega_total_input(n) = totsca / totext
        greekmat_total_input(0:2,n,1) = raywt * RayCoeffs(0:2,1) + aerwt * FmatCoeffs(n,0:2,1)
        greekmat_total_input(3:ngreek_moments_input,n,1) =  aerwt * FmatCoeffs(n,3:ngreek_moments_input,1)
        fmatrix_up(n,1:n_geoms,1) = raywt * RayFmatrices_up(1:n_geoms,1) + aerwt * OutFmatrices_up(1:n_geoms,n,1)
        fmatrix_dn(n,1:n_geoms,1) = raywt * RayFmatrices_dn(1:n_geoms,1) + aerwt * OutFmatrices_dn(1:n_geoms,n,1)
      enddo

!  Copy to type-structure inputs

      VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP = fmatrix_up
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN = fmatrix_dn

!  Number of tasks
!   --- Version 2.7 (3/22/15) increase to 8 tasks (rename from tasks)
!   --- Version 2.8 (9/18/16) reduced  to 6 tasks

      ntasks = 6
!      ntasks = 1

!  Start task loop

      do task = 1, ntasks
!      do task = 4, 6, 2

        t = task

!  2p8. The 6 tasks are now the following-----
!    Task 1: No Single-scatter (FO) correction, no delta-M scaling
!    Task 2: No Single-scatter (FO) correction, With delta-M scaling
!    Task 3: FO Regular  PS, Single-scatter correction, With delta-M scaling (SUN in curved atmosphere)
!    Task 4: FO Enhanced PS, Single-scatter correction, With delta-M scaling (SUN+LOS in curved atmosphere)
!    Task 5: As task 4, but with Solution Saving
!    Task 6: As task 4, but with Solution Saving and BVP Telescoping

!  2p7. These two tasks were included in Version 2.7, Now removed for version 2.8
!     Task 3: Ingoing-only Single-scatter correction, With delta-M scaling (SUN in curved atmosphere)
!     Task 4: In/Outgoing  Single-scatter correction, With delta-M scaling (SUN+LOS in curved atmosphere)

!  2p8. Now has partials in the FO code.

        if ( task .eq. 1 ) then
           DO_FOCORR          = .false.
           DO_FOCORR_NADIR    = .false.
           DO_FOCORR_OUTGOING = .false.
           DO_DELTAM_SCALING  = .false.
           DO_SOLUTION_SAVING = .false.
           DO_BVP_TELESCOPING = .false.
           NFINELAYERS        = 0
       else if ( task .eq. 2 ) then
           DO_FOCORR          = .false.
           DO_FOCORR_NADIR    = .false.
           DO_FOCORR_OUTGOING = .false.
           DO_DELTAM_SCALING  = .true.
           DO_SOLUTION_SAVING = .false.
           DO_BVP_TELESCOPING = .false.
           NFINELAYERS        = 0
       else if ( task .eq. 3 ) then
           DO_FOCORR          = .true.
           DO_FOCORR_NADIR    = .true.
           DO_FOCORR_OUTGOING = .false.
           DO_DELTAM_SCALING  = .true.
           DO_SOLUTION_SAVING = .false.
           DO_BVP_TELESCOPING = .false.
           NFINELAYERS        = 0
        else if ( task .eq. 4 ) then
           DO_FOCORR          = .true.
           DO_FOCORR_NADIR    = .false.
           DO_FOCORR_OUTGOING = .true.
           DO_DELTAM_SCALING  = .true.
           DO_SOLUTION_SAVING = .false.
           DO_BVP_TELESCOPING = .false.
           NFINELAYERS        = 4            ! Non zero here
        else if ( task .eq. 5 ) then
           DO_FOCORR          = .true.
           DO_FOCORR_NADIR    = .true.
           DO_FOCORR_OUTGOING = .false.
           DO_DELTAM_SCALING  = .true.
           DO_SOLUTION_SAVING = .true.
           DO_BVP_TELESCOPING = .false.
           NFINELAYERS        = 4           ! Non zero here
        else if ( task .eq. 6 ) then
           DO_FOCORR          = .true.
           DO_FOCORR_NADIR    = .false.
           DO_FOCORR_OUTGOING = .true.
           DO_DELTAM_SCALING  = .true.
           DO_SOLUTION_SAVING = .true.
           DO_BVP_TELESCOPING = .true.
           NFINELAYERS        = 4            ! Non zero here
        endif
        write(*,*)'Doing task # ', task

!  Copy to type-structure inputs.
!    Rob 3/22/15/. Add FO_CALC setting, now named DO_FOCORR(V2.8)

        VLIDORT_ModIn%MBool%TS_DO_FOCORR          = DO_FOCORR
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = DO_FOCORR_NADIR
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = DO_FOCORR_OUTGOING
        VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = DO_DELTAM_SCALING
        VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = DO_SOLUTION_SAVING
        VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = DO_BVP_TELESCOPING
        VLIDORT_FixIn%Cont%TS_NFINELAYERS         = NFINELAYERS

!  VLIDORT call

        do_debug_input = .false.
        CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling, write-up (optional)

        CALL VLIDORT_WRITE_STATUS ( &
          '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

        n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Save results. Version 2.8, don't need the utamap.

        do uta = 1, n_user_levels
          do o1 = 1, nstokes
            do v = 1, n_geometries
              stokes_PT(uta,v,o1,upidx,t) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
              stokes_PT(uta,v,o1,dnidx,t) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            enddo
            do v = 1, n_szangles
              MEANST_PT(uta,v,o1,upidx,t) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
              MEANST_PT(uta,v,o1,dnidx,t) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                          + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
              FLUX_PT(uta,v,o1,upidx,t)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
              FLUX_PT(uta,v,o1,dnidx,t)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                          + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)
            enddo
          enddo
        enddo

!  End multi-task loop

      enddo

!  Set some local variables

      n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Headings for output file
!    -- Versioni 2.7 Rob 3/22/15. increase to 8 tasks (rename from tasks)
!    -- Versioni 2.8 Rob 9/18/16. reduced  to 6 tasks now

      OutUnit = 36
      OPEN(36,file = 'vlidort_s_test/results_vfzmat_tester.all', status = 'unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes &
         ( OutUnit, 1, 'I', max_tasks, n_geometries, n_szangles,  n_user_levels, 6, &
           user_levels, STOKES_PT, MEANST_PT, FLUX_PT )

!  Close file

      close(OutUnit)

!  Finish

      write(*,*)
      write(*,*)'Main program finished successfully'

      stop

      end program Vfzmat_Tester

!******************************************************************************
      subroutine write_stokes &
        ( un, o1, a1, max_tasks, n_geometries, n_szangles,  n_user_levels, nt, &
          user_levels, STOKES_PT, MEANST_PT, FLUX_PT )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, &
                                 maxstokes, max_directions, upidx, dnidx

      implicit none

!  input variables

      integer, intent(in) ::           max_tasks, n_geometries, n_szangles, n_user_levels
      integer, intent(in) ::           un, nt, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in)  :: user_levels ( max_user_levels )

      DOUBLE PRECISION, INTENT(IN)  :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN)  :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN)  :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAX_TASKS )

!  local variables

      integer :: n, v, t

!  Write stokes component

      write(un,'(/T32,a/T32,a/)') &
        '  STOKES-'//a1//' , tasks 1-6', &
        '  ===================='
      write(un,'(a,T32,6(a16,2x)/)')'Geometry    Level/Output', &
                                    'No FOCorr No DM','No FOCorr + DM ', &
                                    'FO_ReglrPS + DM','FO_EnhncPs + DM', &
                                    '#4 + Solsaving ','#5 + BVPTelscpe'

      do v = 1, n_geometries
        do n = 1, n_user_levels
          write(un,366)v,'Upwelling @',user_levels(n), (stokes_PT(n,v,o1,upidx,t),t=1,nt)
        enddo
        do n = 1, n_user_levels
          write(un,366)v,'Dnwelling @',user_levels(n), (stokes_PT(n,v,o1,dnidx,t),t=1,nt)
        enddo
        write(un,*)' '
      enddo

!  Write fluxes

      write(un,'(/T32,a/T32,a/)') &
        'Stokes-'//a1//' ACTINIC + REGULAR FLUXES, tasks 1-2', &
        '============================================'
      write(un,'(a,T31,4(a16,2x)/)')' Sun SZA    Level/Output', &
                                    ' Actinic No DM ',' Actinic + DM  ', &
                                    ' Regular No DM ',' Regular + DM  '

      do v = 1, n_szangles
        do n = 1, n_user_levels
          write(un,367)v,'Upwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,upidx,t),t=1,2), &
                (FLUX_PT(n,v,o1,upidx,t),t=1,2)
        enddo
        do n = 1, n_user_levels
          write(un,367)v,'Dnwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,dnidx,t),t=1,2), &
                (FLUX_PT(n,v,o1,dnidx,t),t=1,2)
        enddo
        write(un,*)' '
      enddo

366   format(i5,T11,a,f6.2,2x,8(1x,1pe16.7,1x))
367   format(i5,T11,a,f6.2,2x,4(1x,1pe16.7,1x))

      end subroutine write_stokes
