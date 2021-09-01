
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

      program Thermal_Tester

!  This is the Version 2.8 driver. Created 9/18/16 from 2.7 Driver
!     R. Spurr. RT Solutions Inc.

!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VBRDF_SUP_AUX_m, Only : VBRDF_READ_ERROR
      USE VBRDF_SUP_MOD_m

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

      USE VLIDORT_GETPLANCK_m
      USE VLIDORT_AUX_m,    Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m, Only : VLIDORT_INPUT_MASTER, VLIDORT_Sup_Init, VLIDORT_BRDF_Sup_Init
      USE VLIDORT_MASTERS_m

      USE VLIDORT_VBRDF_SUP_ACCESSORIES_m

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

!  VLIDORT BRDF supplement input structure

      TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

!  VLIDORT BRDF supplement output structure

      TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out

!  VLIDORT VBRDF supplement output status structure

      TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus

!  VBRDF supplement / VLIDORT VBRDF-related inputs consistency check status. Rob new 3/22/15

      TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_VBRDFCheck_Status

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  Saved variables

      INTEGER ::          N_SZANGLES_SAVE, N_USER_RELAZMS_SAVE

!  Help variables
!   Rob 3/22/15. renamed thread --> task

      INTEGER ::          NGREEKMAT_ENTRIES, L, LDUM, LCHECK
      INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)

      INTEGER ::          M, N,N6,UTA,UTAM,NDUM,V,T,NTASKS,TASK,O1,I
      DOUBLE PRECISION :: KD, GAER, WAER, TAER, PARCEL, RAYWT, AERWT
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT
      DOUBLE PRECISION :: MOLOMG ( MAXLAYERS ), MOLEXT ( MAXLAYERS )
      DOUBLE PRECISION :: AMOMS(6), AERMOMS(6,0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), DEPOL, BETA2
      DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM

      INTEGER ::          OutUnit

!  Thermal emission setups

      LOGICAL ::            THERMFAIL
      CHARACTER (LEN=70) :: THERMALMESSAGE
      DOUBLE PRECISION ::   SURFTEMP, AIRTEMPS(0:MAXLAYERS)
      DOUBLE PRECISION ::   WNUMLO, WNUMHI
      INTEGER ::            NTEMPS, SMALLV

!  BRDF supplement variables

      LOGICAL ::          BS_DO_BRDF_SURFACE
      LOGICAL ::          BS_DO_USER_STREAMS
      LOGICAL ::          BS_DO_SURFACE_EMISSION

      INTEGER ::          BS_NSTOKES

      INTEGER ::          BS_NSTREAMS
      INTEGER ::          BS_NMOMENTS_INPUT

      INTEGER ::          BS_NBEAMS
      DOUBLE PRECISION :: BS_BEAM_SZAS ( MAXBEAMS )

      INTEGER ::          BS_N_USER_RELAZMS
      DOUBLE PRECISION :: BS_USER_RELAZMS (MAX_USER_RELAZMS)

      INTEGER ::          BS_N_USER_STREAMS
      DOUBLE PRECISION :: BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

      LOGICAL ::          DO_DEBUG_RESTORATION

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_LAMBERTIAN_SURFACE
      LOGICAL ::          DO_THERMAL_EMISSION, DO_SURFACE_EMISSION
      LOGICAL ::          DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY
      LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
      LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      INTEGER ::          NLAYERS, NSTOKES, NFINELAYERS, N_THERMAL_COEFFS

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

      DOUBLE PRECISION :: SURFBB, THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: FLUX_FACTOR

      INTEGER ::          NSTREAMS
      INTEGER ::          N_SZANGLES
      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )
      LOGICAL ::          DO_USER_VZANGLES
      INTEGER ::          N_USER_VZANGLES
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER ::          N_USER_RELAZMS
      DOUBLE PRECISION :: USER_RELAZMS ( MAX_USER_RELAZMS )
      LOGICAL ::          DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY

!  Proxies (Map no longer required for Version 2.8)

!      INTEGER ::          UTAMAP(MAX_USER_LEVELS,10)
      INTEGER ::          N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  VLIDORT standard output preparation

      INTEGER ::          N_GEOMETRIES

!  Saved results

      INTEGER, PARAMETER :: MAX_TASKS = 10

      DOUBLE PRECISION :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS  )
      DOUBLE PRECISION :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS  )
      DOUBLE PRECISION :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAX_TASKS  )

!  #####################################
!      START PROGRAM
!  #####################################

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Initialise error file output flag

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_ReadInput.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

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

!  Initialize some input variables not handled by VLIDORT_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )

!  Set Output-level Proxies (saved values)

      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels(1:n_user_levels) = VLIDORT_ModIn%MUserVal%TS_user_levels(1:n_user_levels)

!  Version 2.7. Map is necessary because FO output is ONLY AT LEVEL BOUNDARIES
!  Version 2.8. Map no longer required, as FO code now has partials
!      utamap = 0
!      do uta = 1, n_user_levels
!         utamap(uta,1:4) = uta
!      enddo
!      utamap(1,5:8) = 1 ; utamap(2,5:8) = 2 ; utamap(5,5:8) = 3

!  Rob 3/22/15. Add wavelength Diagnostic (just a dummy here)
!   *** Important for checking against BRDF wavelength.
!   Check is only done for the "New Cox-Munk case" (not relevant here)
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 3.333d0

!  Copy to some local variables

      NSTOKES                 = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS                = VLIDORT_FixIn%Cont%TS_NSTREAMS
      N_SZANGLES              = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      SZANGLES                = VLIDORT_ModIn%MSunrays%TS_SZANGLES
      DO_USER_VZANGLES        = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      N_USER_VZANGLES         = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      USER_VZANGLES           = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT
      N_USER_RELAZMS          = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS            = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS
      DO_OBSERVATION_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      DO_DOUBLET_GEOMETRY     = VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Save some values

      N_SZANGLES_SAVE     = N_SZANGLES
      N_USER_RELAZMS_SAVE = N_USER_RELAZMS

!  Get the pre-prepared atmosphere
!  -------------------------------

      height_grid = 0.0d0
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      open(1,file='vlidort_v_test/ProblemIII.Moms',status = 'unknown')
      do L = 1, 5
        read(1,*)
      enddo
      LCHECK = 0
      L      = -1
      cutoff = 1.0d-06
      DO while (L.lt.106 .and. LCHECK.lt.2)
         read(1,*)LDUM,(amoms(k),k=1,6)
         if ( amoms(1).gt.cutoff ) then
           L = L + 1
           do k = 1, 6
             aermoms(k,L) = amoms(k)
           enddo
           if ( Lcheck.eq.1 ) Lcheck = Lcheck + 1
         else
           Lcheck = Lcheck + 1
           if ( Lcheck .eq. 1 ) then
             L = L + 1
             do k = 1, 6
               aermoms(k,L) = amoms(k)
             enddo
           endif
         endif
      enddo
      close(1)
      ngreek_moments_input = L

!  Rayleigh scattering law

!  @@@ Rob Fix 6/13/13. Formula for depol was wrong, Index was wrong.
!                       However, this only affects the circular-polarization
!                       Also the (3,1) instead of (4,1) will not affect result

!  Old code
!      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / &
!                ( 1.0d0 + 2.0d0*raymoms(2,NLAYERS) )
      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / &
                ( 1.0d0 + raymoms(2,NLAYERS) )

      beta2 = raymoms(2,NLAYERS)
      PROBLEM_RAY = 0.0d0
      PROBLEM_RAY(1,0) =  1.D0

!  Old code, index was wrong
!      PROBLEM_RAY(3,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(4,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - DSQRT(6.D0) * beta2
      PROBLEM_RAY(2,2) =  6.D0 * beta2

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Thermal stuff
!  -------------

!   Array of temperatures should be 24 levels

      open(1,file= 'vlidort_v_test/input_temp23.dat', status='old')
      read(1,*)SURFTEMP
      read(1,*)NTEMPS
      if ( NTEMPS .ne. NLAYERS ) then
         write(*,*)'Number of layers in temperature data set not equal to input NLAYERS'
         stop 'Test program V2p8p3_thermal_tester stopped'
      endif
      do n = 0, NTEMPS
         read(1,*)ndum, AIRTEMPS(NTEMPS-N)
      enddo
      close(1)

!  Do the BB inputs

      wnumhi = 2500.0d0 + 1.0d0
      wnumlo = 2500.0d0 - 1.0d0
      call get_planckfunction &
         ( wnumlo,wnumhi, SURFTEMP, &
           SURFBB, SMALLV, THERMFAIL, THERMALMESSAGE )
      do n = 0, nlayers
         call get_planckfunction &
         ( wnumlo,wnumhi, AIRTEMPS(N), &
           THERMAL_BB_INPUT(N), SMALLV, THERMFAIL, THERMALMESSAGE )
      enddo
      if ( THERMFAIL ) THEN
         write(*,'(A)')Trim(THERMALMESSAGE)
         stop 'Failure from call to get_planckfunction in program V2p8p3_thermal_tester'
      endif

!  Set some VLIDORT variables

      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = SURFBB
      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS) = THERMAL_BB_INPUT(0:NLAYERS)

!     pause'thermal check'

!  Main section
!  ------------

!  Number of tasks

!     ntasks = 4        !  thermal, 4 Lambertian
      ntasks = 6        !  thermal, 4 Lambertian + 2 BRDF

!  Task 1: Thermal only, with scattering, + delta-M, Lambertian
!  Task 2: Thermal transmittance only
!  Task 3: Crossover FOCORR_NADIR    FO_Single-scatter correction + delta-M, Lamb
!  Task 4: Crossover FOCORR_OUTGOING FO Single-scatter correction + delta-M, Lamb
!  Task 5: Crossover FOCORR_OUTGOING FO Single-scatter correction + delta-M, BRDF1
!  Task 6: Crossover FOCORR_OUTGOING FO Single-scatter correction + delta-M, BRDF3

!  Set some VLIDORT variables

      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS    = 2
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION = .TRUE.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION = .TRUE.

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .TRUE.
      !VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .TRUE.

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = 1.0d-03

!  Set some local variables

      DO_SURFACE_EMISSION = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION

!  Initialise output arrays

      STOKES_PT = ZERO
      MEANST_PT = ZERO
      FLUX_PT   = ZERO

      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT    = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS  = .false.

!  Start task loop

       do task = 1, ntasks
       !do task = 5, 5
       !do task = 6, 6

!  Initialise

        deltau_vert_input = zero
        omega_total_input = zero
        greekmat_total_input = zero
        Fmatrix_up           = zero
        Fmatrix_dn           = zero

!  Create optical properties, every task

        n6 = nlayers - 6
        do n = 1, n6
          deltau_vert_input(n) = molext(n)
          omega_total_input(n) = molomg(n)
          DO K = 1, ngreekmat_entries
            GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
            DO L = 0, 2
! Old         GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(ck,L)
              GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
            ENDDO
            GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
          ENDDO
        enddo

!  Fix layers with aerosol

        waer = 0.99999d0 ; taer = 0.5d0
        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          aersca = aerext * waer
          molsca = molomg(n) * molext(n)
          totext = molext(n) + aerext
          totsca = molsca    + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          deltau_vert_input(n) = totext
          omega_total_input(n) = totsca / totext
          DO K = 1, ngreekmat_entries
            GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
            DO L = 0, 2
! Old         MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
              MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*SK*AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,N,GK) = MOM
            ENDDO
            DO L = 3, NGREEK_MOMENTS_INPUT
              MOM = AERWT * SK * AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,n,gk) = MOM
            ENDDO
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
        enddo

!  Copy to type-structure inputs

        VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ngreek_moments_input
        VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid

        VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
        VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
        VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

        VLIDORT_FixIn%Optical%TS_FMATRIX_UP           = fmatrix_up ! zero here
        VLIDORT_FixIn%Optical%TS_FMATRIX_DN           = fmatrix_dn ! zero here

!  Task assignation

        t = task
        if ( task .eq. 1 ) then
           DO_SOLAR_SOURCES      = .false.
           DO_THERMAL_TRANSONLY  = .false.
           DO_FOCORR             = .false.
           DO_FOCORR_NADIR       = .false.
           DO_FOCORR_OUTGOING    = .false.
           DO_DELTAM_SCALING     = .true.
           NFINELAYERS           = 0
           DO_LAMBERTIAN_SURFACE = .true.
           LAMBERTIAN_ALBEDO     = 0.05d0
           N_SZANGLES     = 1
           N_USER_RELAZMS = 1
        else if ( task .eq. 2 ) then
           DO_SOLAR_SOURCES      = .false.
           DO_THERMAL_TRANSONLY  = .true.
           DO_FOCORR             = .false.
           DO_FOCORR_NADIR       = .false.
           DO_FOCORR_OUTGOING    = .false.
           DO_DELTAM_SCALING     = .false.
           NFINELAYERS           = 0
           DO_LAMBERTIAN_SURFACE = .true.
           LAMBERTIAN_ALBEDO     = 0.05d0
           N_SZANGLES     = 1
           N_USER_RELAZMS = 1
        else if ( task .eq. 3 ) then
           DO_SOLAR_SOURCES      = .true.
           DO_THERMAL_TRANSONLY  = .false.
           DO_FOCORR             = .true.
           DO_FOCORR_NADIR       = .true.
           DO_FOCORR_OUTGOING    = .false.
           DO_DELTAM_SCALING     = .true.
           NFINELAYERS           = 0
           DO_LAMBERTIAN_SURFACE = .true.
           LAMBERTIAN_ALBEDO     = 0.05d0
           N_SZANGLES     = N_SZANGLES_SAVE
           N_USER_RELAZMS = N_USER_RELAZMS_SAVE
        else if ( task .eq. 4 ) then
           DO_SOLAR_SOURCES      = .true.
           DO_THERMAL_TRANSONLY  = .false.
           DO_FOCORR             = .true.
           DO_FOCORR_NADIR       = .false.
           DO_FOCORR_OUTGOING    = .true.
           DO_DELTAM_SCALING     = .true.
           NFINELAYERS           = 4            ! Non zero here
           DO_LAMBERTIAN_SURFACE = .true.
           LAMBERTIAN_ALBEDO     = 0.05d0
           N_SZANGLES     = N_SZANGLES_SAVE
           N_USER_RELAZMS = N_USER_RELAZMS_SAVE
        else if ( task .eq. 5 ) then
           DO_SOLAR_SOURCES      = .true.
           DO_THERMAL_TRANSONLY  = .false.
           DO_FOCORR             = .true.
           DO_FOCORR_NADIR       = .false.
           DO_FOCORR_OUTGOING    = .true.
           DO_DELTAM_SCALING     = .true.
           NFINELAYERS           = 4
           DO_LAMBERTIAN_SURFACE = .false.      ! BRDF surface here
           LAMBERTIAN_ALBEDO     = 0.0d0
           N_SZANGLES     = N_SZANGLES_SAVE
           N_USER_RELAZMS = N_USER_RELAZMS_SAVE
        else if ( task .eq. 6 ) then
           DO_SOLAR_SOURCES      = .true.
           DO_THERMAL_TRANSONLY  = .false.
           DO_FOCORR             = .true.
           DO_FOCORR_NADIR       = .false.
           DO_FOCORR_OUTGOING    = .true.
           DO_DELTAM_SCALING     = .true.
           NFINELAYERS           = 4            ! Non zero here
           DO_LAMBERTIAN_SURFACE = .false.      ! BRDF surface here
           LAMBERTIAN_ALBEDO     = 0.0d0
           N_SZANGLES     = N_SZANGLES_SAVE
           N_USER_RELAZMS = N_USER_RELAZMS_SAVE
        endif

        write(*,*)'Doing task # ', task

!  Copy to type-structure inputs

        VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES     = DO_SOLAR_SOURCES
        VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY = DO_THERMAL_TRANSONLY

        VLIDORT_ModIn%MBool%TS_DO_FOCORR            = DO_FOCORR          ! New 3/22/15
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR      = DO_FOCORR_NADIR
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING   = DO_FOCORR_OUTGOING

        VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = DO_DELTAM_SCALING
        VLIDORT_FixIn%Cont%TS_NFINELAYERS           = NFINELAYERS
        VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
        VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO  = LAMBERTIAN_ALBEDO
        VLIDORT_ModIn%MSunrays%TS_N_SZANGLES        = N_SZANGLES
        VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS    = N_USER_RELAZMS

!  For Lambertian cases, set BRDF inputs to zero

        if ( task .lt. 5 ) then
          CALL VLIDORT_BRDF_Sup_Init ( VLIDORT_Sup )
          VLIDORT_Sup%BRDF%TS_EMISSIVITY(1,:)      = 1.0d0 - LAMBERTIAN_ALBEDO
          VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1,:) = 1.0d0 - LAMBERTIAN_ALBEDO
        endif

!  BRDF supplement call if task is 5 
!    Hard-wire the supplement inputs. gcmcri

        if ( task .eq. 5 ) then

!  Set BS inputs (hard-wired, copy VLIDORT Main inputs)

          BS_NMOMENTS_INPUT    = 2*NSTREAMS - 1
          DO_DEBUG_RESTORATION = .false.

          VBRDF_Sup_In%BS_DO_BRDF_SURFACE     = .true.
          VBRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_VZANGLES
          VBRDF_Sup_In%BS_DO_SOLAR_SOURCES    = DO_SOLAR_SOURCES
          VBRDF_Sup_In%BS_DO_USER_OBSGEOMS    = DO_OBSERVATION_GEOMETRY
          VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY
          VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

          VBRDF_Sup_In%BS_NSTOKES  = NSTOKES
          VBRDF_Sup_In%BS_NSTREAMS = NSTREAMS

          VBRDF_Sup_In%BS_NBEAMS    = N_SZANGLES
          VBRDF_Sup_In%BS_BEAM_SZAS = ZERO
          DO I = 1, N_SZANGLES
            VBRDF_Sup_In%BS_BEAM_SZAS(I) = SZANGLES(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_STREAMS    = N_USER_VZANGLES
          VBRDF_Sup_In%BS_USER_ANGLES_INPUT = ZERO
          DO I = 1, N_USER_VZANGLES
            VBRDF_Sup_In%BS_USER_ANGLES_INPUT(I) = USER_VZANGLES(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_RELAZMS = N_USER_RELAZMS
          VBRDF_Sup_In%BS_USER_RELAZMS   = ZERO
          DO I = 1,N_USER_RELAZMS
            VBRDF_Sup_In%BS_USER_RELAZMS(I) = USER_RELAZMS(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_OBSGEOMS = 0
          VBRDF_Sup_In%BS_USER_OBSGEOMS   = ZERO

!  Set BS inputs (supplement only, by hand)
!    9/29/14 Change names DBONLY, remove DO_EXACT flag

!  Misc control

          VBRDF_Sup_In%BS_DO_SHADOW_EFFECT     = .true.
!          VBRDF_Sup_In%BS_DO_EXACT             = .true.
          VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY = .false.

          VBRDF_Sup_In%BS_NSTREAMS_BRDF     = 100
          VBRDF_Sup_In%BS_N_BRDF_KERNELS    = 1

          VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .false.

!  Multiple-scattering Glitter options

          VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = .false.
          VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = .false.
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = 0
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = 40
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = 100

!  Initialize

          VBRDF_Sup_In%BS_WHICH_BRDF        = 0
          VBRDF_Sup_In%BS_BRDF_NAMES        = '          '
          VBRDF_Sup_In%BS_BRDF_FACTORS      = 0.0d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS   = 0.0d0

!  Specific BRDF inputs

!  CoxMunk. parameters are sig-sq and RFindex-sq
          !VBRDF_Sup_In%BS_WHICH_BRDF(1)        = 9
          !VBRDF_Sup_In%BS_BRDF_NAMES(1)        = 'Cox-Munk  '
          !VBRDF_Sup_In%BS_BRDF_FACTORS(1)      = 0.1d0
          !VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1) = 2
          !VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1) = 0.0798d0
          !VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2) = 1.779556d0

!  GCMComplex, parameters are sig-sq, Real(RF), Imag(RF)
          VBRDF_Sup_In%BS_WHICH_BRDF(1)        = 11
          VBRDF_Sup_In%BS_BRDF_NAMES(1)        = 'GCMcomplex'
          VBRDF_Sup_In%BS_BRDF_FACTORS(1)      = 1.0d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1) = 3
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1) = 0.0798d0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2) = 1.32d0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3) = 0.0001d0

!  Block of New 2.7 variables (scaling & NCM - all dummies here)

!  Albedo scaling
!   9/29/14. Add WSABSA OUTPUT

          VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .false.
          VBRDF_Sup_In%BS_DO_WSA_SCALING   = .false.
          VBRDF_Sup_In%BS_DO_BSA_SCALING   = .false.
          VBRDF_Sup_In%BS_WSA_VALUE        = 0.0d0
          VBRDF_Sup_In%BS_BSA_VALUE        = 0.0d0

!  New Cox-Munk

          VBRDF_Sup_In%BS_DO_NewCMGLINT    = .false.
          VBRDF_Sup_In%BS_DO_NewGCMGLINT   = .false.
          VBRDF_Sup_In%BS_SALINITY         = 0.0d0
          VBRDF_Sup_In%BS_WAVELENGTH       = 3333.0d0   !Diagnostic
          VBRDF_Sup_In%BS_WINDSPEED        = 0.0d0
          VBRDF_Sup_In%BS_WINDDIR          = 0.0d0
          VBRDF_Sup_In%BS_DO_GlintShadow   = .false.
          VBRDF_Sup_In%BS_DO_FoamOption    = .false.
          VBRDF_Sup_In%BS_DO_FacetIsotropy = .false.

        ENDIF

!  BRDF supplement call if task is 6 
!    Hard-wire the supplement inputs. 3-kernel example.

        if ( task .eq. 6  ) then

!  Set BS inputs (hard-wired, copy VLIDORT Main inputs)

          BS_NMOMENTS_INPUT    = 2*NSTREAMS - 1
          DO_DEBUG_RESTORATION = .false.

          VBRDF_Sup_In%BS_DO_BRDF_SURFACE     = .true.
          VBRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_VZANGLES
          VBRDF_Sup_In%BS_DO_SOLAR_SOURCES    = DO_SOLAR_SOURCES
          VBRDF_Sup_In%BS_DO_USER_OBSGEOMS    = DO_OBSERVATION_GEOMETRY
          VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY
          VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

          VBRDF_Sup_In%BS_NSTOKES  = NSTOKES
          VBRDF_Sup_In%BS_NSTREAMS = NSTREAMS

          VBRDF_Sup_In%BS_NBEAMS = N_SZANGLES
          VBRDF_Sup_In%BS_BEAM_SZAS = ZERO
          DO I = 1, N_SZANGLES
            VBRDF_Sup_In%BS_BEAM_SZAS(I) = SZANGLES(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_STREAMS = N_USER_VZANGLES
          VBRDF_Sup_In%BS_USER_ANGLES_INPUT = ZERO
          DO I = 1, N_USER_VZANGLES
            VBRDF_Sup_In%BS_USER_ANGLES_INPUT(I) = USER_VZANGLES(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_RELAZMS = N_USER_RELAZMS
          VBRDF_Sup_In%BS_USER_RELAZMS   = ZERO
          DO I = 1,N_USER_RELAZMS
            VBRDF_Sup_In%BS_USER_RELAZMS(I) = USER_RELAZMS(I)
          ENDDO

          VBRDF_Sup_In%BS_N_USER_OBSGEOMS = 0
          VBRDF_Sup_In%BS_USER_OBSGEOMS   = ZERO

!  Set BS inputs (supplement only, by hand)
!    9/29/14 Change names DBONLY, remove DO_EXACT flag

!  Misc control

          VBRDF_Sup_In%BS_DO_SHADOW_EFFECT     = .true.
!          VBRDF_Sup_In%BS_DO_EXACT             = .true.
          VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY = .false.

          VBRDF_Sup_In%BS_NSTREAMS_BRDF    = 100
          VBRDF_Sup_In%BS_N_BRDF_KERNELS   = 3

          VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .false.

!  Multiple-scattering Glitter options

          VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = .false.
          VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = .false.
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = 0
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = 40
          VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = 100

!  Initialize

          VBRDF_Sup_In%BS_WHICH_BRDF        = 0
          VBRDF_Sup_In%BS_BRDF_NAMES        = '          '
          VBRDF_Sup_In%BS_BRDF_FACTORS      = 0.0d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS   = 0.0d0

!  Specific BRDF inputs: 3-kernel

          VBRDF_Sup_In%BS_WHICH_BRDF(1)   = 2
          VBRDF_Sup_In%BS_BRDF_NAMES(1)   = 'Ross-thin '
          VBRDF_Sup_In%BS_BRDF_FACTORS(1) = 0.3d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1) = 0

          VBRDF_Sup_In%BS_WHICH_BRDF(2)   = 5
          VBRDF_Sup_In%BS_BRDF_NAMES(2)   = 'Li-dense  '
          VBRDF_Sup_In%BS_BRDF_FACTORS(2) = 0.1d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(2) = 2
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,1) = 2.0d0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,2) = 1.0d0

          VBRDF_Sup_In%BS_WHICH_BRDF(3)   = 9
          VBRDF_Sup_In%BS_BRDF_NAMES(3)   = 'Cox-Munk  '
          VBRDF_Sup_In%BS_BRDF_FACTORS(3) = 0.1d0
          VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(3) = 2
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(3,1) = 0.0798d0
          VBRDF_Sup_In%BS_BRDF_PARAMETERS(3,2) = 1.779556d0

!  Block of New 2.7 variables (scaling & NCM - all dummies here)

!  Albedo scaling
!   9/29/14. Add WSABSA OUTPUT

          VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .false.
          VBRDF_Sup_In%BS_DO_WSA_SCALING   = .false.
          VBRDF_Sup_In%BS_DO_BSA_SCALING   = .false.
          VBRDF_Sup_In%BS_WSA_VALUE        = 0.0d0
          VBRDF_Sup_In%BS_BSA_VALUE        = 0.0d0

!  New Cox-Munk

          VBRDF_Sup_In%BS_DO_NewCMGLINT    = .false.
          VBRDF_Sup_In%BS_DO_NewGCMGLINT   = .false.
          VBRDF_Sup_In%BS_SALINITY         = 0.0d0
          VBRDF_Sup_In%BS_WAVELENGTH       = 3333.0d0  !Diagnostic
          VBRDF_Sup_In%BS_WINDSPEED        = 0.0d0
          VBRDF_Sup_In%BS_WINDDIR          = 0.0d0
          VBRDF_Sup_In%BS_DO_GlintShadow   = .false.
          VBRDF_Sup_In%BS_DO_FoamOption    = .false.
          VBRDF_Sup_In%BS_DO_FacetIsotropy = .false.

        ENDIF

        IF ( TASK .GT. 4 ) THEN

!  Accessory check. Added 3/22/15
!  This is the checking routine

          CALL VLIDORT_VBRDF_INPUT_CHECK ( &
            VBRDF_Sup_In,             & ! Inputs
            VLIDORT_FixIn,            & ! Inputs
            VLIDORT_ModIn,            & ! Inputs
            VLIDORT_VBRDFCheck_Status ) ! Outputs

          IF ( VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK .ne. VLIDORT_SUCCESS ) &
            CALL VLIDORT_VBRDF_INPUT_CHECK_ERROR ( &
              'V2p8p3_VLIDORT_VBRDFcheck.log', VLIDORT_VBRDFCheck_Status )

!  Perform supplement calculation for task 5 and 6
!    The output will now be used for the Main calculation

          CALL VBRDF_MAINMASTER ( &
            DO_DEBUG_RESTORATION,   &  ! Inputs
            BS_NMOMENTS_INPUT,      &  ! Inputs
            VBRDF_Sup_In,           &  ! Inputs
            VBRDF_Sup_Out,          &  ! Outputs
            VBRDF_Sup_OutputStatus )   ! Output Status

          IF ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. VLIDORT_SUCCESS ) then
            DO M = 1, VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
               write(*,'(a,I3,A,A)')'BRDF call Failure: Message # ',M,' = ',&
                           trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(M))
            enddo
            Stop 'VBRDF supplement master failed, in program V2p8p3_thermal_tester'
          ENDIF     

!  Copy VBRDF Sup outputs to VLIDORT's VBRDF Sup inputs (std only)

          CALL SET_VLIDORT_VBRDF_INPUTS ( &
            VBRDF_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn, & !Inputs
            VLIDORT_Sup )                                  !Outputs

        ENDIF

!  VLIDORT call

        do_debug_input = .false.
        CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling, write-up (optional)

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Write final file for all LOCAL tasks
!  ====================================

!  Copy some variables

      n_szangles    = VLIDORT_ModIn%MSunrays%TS_n_szangles
      n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Title of output file

      OutUnit = 36
      if ( nstokes .eq. 1 ) then
        OPEN(OutUnit,file='vlidort_v_test/results_thermal_tester_I000.all',status='unknown')
      else if ( nstokes .eq. 3 ) then
        OPEN(OutUnit,file='vlidort_v_test/results_thermal_tester_IQU0.all',status='unknown')
      else if ( nstokes .eq. 4 ) then
        OPEN(OutUnit,file='vlidort_v_test/results_thermal_tester_IQUV.all',status='unknown')
      endif

!  Write Stokes Intensity and Mean values

      call write_stokes_thermal &
         ( OutUnit, 1, 'I', max_tasks, n_geometries, n_szangles,  n_user_levels, ntasks, &
           user_levels, Stokes_PT, MEANST_PT, FLUX_PT )

!  Write Stokes Q and U and Mean values

      if ( NSTOKES .ge. 3 ) then
        call write_stokes_thermal &
         ( OutUnit, 2, 'Q', max_tasks, n_geometries, n_szangles,  n_user_levels, ntasks, &
           user_levels, Stokes_PT, MEANST_PT, FLUX_PT )

        call write_stokes_thermal &
         ( OutUnit, 3, 'U', max_tasks, n_geometries, n_szangles,  n_user_levels, ntasks, &
           user_levels, Stokes_PT, MEANST_PT, FLUX_PT )
      endif

!  Write Stokes V and Mean values

      if ( NSTOKES .eq. 4 ) then
        call write_stokes_thermal &
         ( OutUnit, 4, 'V', max_tasks, n_geometries, n_szangles,  n_user_levels, ntasks, &
           user_levels, Stokes_PT, MEANST_PT, FLUX_PT )
      endif

!  Close file

      close(OutUnit)

!  Close error file if it has been opened

      write(*,*)
      IF ( OPENFILEFLAG ) then
        !close( LIDORT_ERRUNIT )   not active yet
        write(*,*)'Main program executed with internal defaults,'// &
                  ' warnings in "V2p8p3_VLIDORT_Execution.log"'
      ELSE
        write(*,*)'Main program finished successfully'
      ENDIF

!  Finish

      stop
      end program Thermal_Tester

!**********************************************************************
      subroutine write_stokes_thermal &
        ( un, o1, a1, max_tasks, n_geometries, n_szangles,  n_user_levels, nt, &
          user_levels, stokes_pt, MEANST_PT, FLUX_PT )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, &
                                 maxstokes, max_directions, upidx, dnidx

      implicit none

!  input variables

      integer, intent(in) ::           max_tasks, n_geometries, n_szangles, n_user_levels
      integer, intent(in) ::           un, nt, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in)  :: user_levels ( max_user_levels )

      double precision, intent(in)  :: stokes_pt &
          ( max_user_levels, max_geometries, maxstokes, max_directions, max_tasks )
      double precision, intent(in)  :: MEANST_PT &
          ( max_user_levels, max_szangles, maxstokes, max_directions, max_tasks )
      double precision, intent(in)  :: FLUX_PT &
          ( max_user_levels, max_szangles, maxstokes ,max_directions, max_tasks )

!  local variables

      integer :: n, v, t

!  Write stokes component

      write(un,'(/T32,a/T32,a/)') &
        '  STOKES-'//a1//' , tasks 1-6', &
        '  ===================='
      write(un,'(a,T32,6(a16,2x)/)')'Geometry    Level/Output', &
                                    'ThermalOnly +DM','ThermTransOnly ', &
                                    'Xo+FORg+DM+Lamb','Xo+FOEh+DM+Lamb', &
                                    'Xo+FOEh+DM+CxMk','Xo+FOEh+DM+3ker'

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
        'Stokes-'//a1//' ACTINIC + REGULAR FLUXES, tasks 1-6', &
        '============================================'
      write(un,'(a,T31,2(a21,87x)/)')' Sun SZA    Level/Output', &
                                      ' Actinic 1 through 6 ', &
                                      ' Regular 1 through 6 '

      do v = 1, n_szangles
        do n = 1, n_user_levels
          write(un,367)v,'Upwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,upidx,t),t=1,nt), &
                (FLUX_PT(n,v,o1,upidx,t),t=1,nt)
        enddo
        do n = 1, n_user_levels
          write(un,367)v,'Dnwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,dnidx,t),t=1,nt), &
                (FLUX_PT(n,v,o1,dnidx,t),t=1,nt)
        enddo
        write(un,*)' '
      enddo

366   format(i5,T11,a,f6.2,2x,8(1x,1pe16.7,1x))
367   format(i5,T11,a,f6.2,2x,16(1x,1pe16.7,1x))

      end subroutine write_stokes_thermal
