
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

      program Solar_Tester

!  This is the Version 2.8 driver. Created 9/18/16 from 2.7 Driver
!     R. Spurr. RT Solutions Inc.

!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

      USE VLIDORT_AUX_m,    Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m, Only : VLIDORT_INPUT_MASTER, VLIDORT_Sup_Init
      USE VLIDORT_MASTERS_m

!  Implicit none

      IMPLICIT NONE

!  Dimensioning integer parameters

      INTEGER, PARAMETER :: MAXWVN = 10!50!100!1000

!  Number of tasks. Rob 3/22/15. Increase to 8 (2 additional FO tests)

      INTEGER, PARAMETER :: MAXTASKS = 6

!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT debug input control

      LOGICAL :: DO_DEBUG_INPUT

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn, VLIDORT_FixIn_SAVE
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn, VLIDORT_ModIn_SAVE

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup, VLIDORT_Sup_SAVE

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  Help variables

      INTEGER ::          NGREEKMAT_ENTRIES, L, LDUM, LCHECK
      INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)

      INTEGER ::          N,N6,UTA,UTAM,NDUM,V,T,NTASKS,TASK,O1
      DOUBLE PRECISION :: KD, GAER, WAER, TAER, PARCEL, RAYWT, AERWT
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT
      DOUBLE PRECISION :: MOLOMG ( MAXLAYERS ), MOLEXT ( MAXLAYERS )
      DOUBLE PRECISION :: AMOMS(6), AERMOMS(6,0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), DEPOL, BETA2
      DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM

      CHARACTER(LEN=4) :: C2

      INTEGER ::          OutUnit

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
      !DOUBLE PRECISION :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 
      !DOUBLE PRECISION :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 

!  Proxies (Map no longer required for Version 2.8)

!      INTEGER ::          UTAMAP(MAX_USER_LEVELS,10)
      INTEGER ::          N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  VLIDORT standard output preparation

      INTEGER ::          N_GEOMETRIES
      INTEGER ::          N_SZANGLES

!  Saved results

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: OMP_N_GEOMETRIES
      DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: STOKES_PT, MEANST_PT, FLUX_PT

!  Local shared array

      INTEGER, DIMENSION ( MAXTASKS, MAXWVN ) :: VLIDORT_TID_SAVE

!  OpenMP tests (general)

      INTEGER       :: TEST, N_OMP_TESTS, OMP_MAXTHREADS
      INTEGER       :: TID, OMP_NTHREADS
      INTEGER       :: wvn, numwvn, ompt

      CHARACTER (LEN=1)  :: charTID
      CHARACTER (LEN=1)  :: CH1 
      CHARACTER (LEN=5)  :: CH2,CH3
      CHARACTER (LEN=30) :: TAIL

!  OpenMP tests (timing)

      INTEGER       :: n_core, time_divider
      REAL          :: omp_e1, omp_e2

!  OpenMP functions

      INTEGER :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!  Test file output logical

      LOGICAL :: DO_STD_OUTPUT_FILE = .TRUE.!.FALSE.

!  Start of code
!  =============

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Set test parameters

      N_CORE = 2!4
      N_OMP_TESTS = 2!3

      numwvn = maxwvn
      ntasks = maxtasks

!  Begin test loop

      DO TEST = 1, N_OMP_TESTS
      !DO TEST = 3, N_OMP_TESTS
      !DO TEST = 1, 1
      !DO TEST = 2, 2

      write(*,*)
      write(*,*) '******************************************'
      write(*,'(1x,a,i1)') 'Doing test ',TEST

!  Set total number of threads (optional)
!      Note: if not set, OpenMP will set the number of threads
!            to the number of cores on the system
!      Note: if not set, OpenMP will default to dividing up the
!            wvn loop below equally among the available threads

      if ( test == 1 ) then
        OMP_MAXTHREADS = 1
      elseif ( test == 2 ) then
        OMP_MAXTHREADS = 2
      elseif ( test == 3 ) then
        OMP_MAXTHREADS = 4
      endif

      !write(*,*)
      !write(*,*) 'OMP_MAXTHREADS = ',OMP_MAXTHREADS

!  Set total number of threads

      !write(*,*)
      !write(*,*) 'Setting total number of threads'
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Start timing

      call cpu_time(omp_e1)

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_ReadInput.cfg', & ! Input
        VLIDORT_FixIn_SAVE, & ! Outputs
        VLIDORT_ModIn_SAVE, & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

!  1/31/21. Version 2.8.3, New variables must be set by hand

      VLIDORT_FixIn_SAVE%Cont%TS_ASYMTX_TOLERANCE = 1.0d-20

!  1/31/21. Version 2.8.3, Set the DO_MSSTS flag to generate output for the MS sphericity corrections

      VLIDORT_FixIn_SAVE%Bool%TS_DO_MSSTS = .false.

!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.

!  1/31/21. Version 2.8.3, Set the NSTOKES2 FOURIER0 flag by hand

      VLIDORT_FixIn_SAVE%Bool%TS_DO_FOURIER0_NSTOKES2 = .true.

!  Initialize some input variables not handled by VLIDORT_INPUT_MASTER
!     These are the same for all Local tasks

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup_SAVE )

!  Set Output-level Proxies (saved values)

      nstokes       = VLIDORT_FixIn_SAVE%Cont%TS_nstokes
      n_user_levels = VLIDORT_FixIn_SAVE%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn_SAVE%MUserVal%TS_user_levels
      n_szangles    = VLIDORT_ModIn_SAVE%MSunrays%TS_n_szangles

!  Version 2.7. Map is necessary because FO output is ONLY AT LEVEL BOUNDARIES
!  Version 2.8. Map no longer required, as FO code now has partials
!      utamap = 0
!      do uta = 1, n_user_levels
!         utamap(uta,1:4) = uta
!      enddo
!      utamap(1,5:8) = 1 ; utamap(2,5:8) = 2 ; utamap(5,5:8) = 3

!  Set output file designation for later use

      if ( nstokes .eq. 1 ) then
        c2 = 'I000'
      else if ( nstokes .eq. 3 ) then
        c2 = 'IQU0'
      else if ( nstokes .eq. 4 ) then
        c2 = 'IQUV'
      endif

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Get the pre-prepared atmosphere
!  -------------------------------

      height_grid = 0.0d0
      nlayers = VLIDORT_FixIn_SAVE%Cont%TS_nlayers
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

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Allocate output memory

      ALLOCATE ( & 
        OMP_N_GEOMETRIES ( MAXTASKS, MAXWVN ), &
        STOKES_PT        ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, &
                           MAX_DIRECTIONS, MAXTASKS, MAXWVN ), &
        MEANST_PT        ( MAX_USER_LEVELS, MAX_SZANGLES,   MAXSTOKES, &
                           MAX_DIRECTIONS, MAXTASKS, MAXWVN ), &
        FLUX_PT          ( MAX_USER_LEVELS, MAX_SZANGLES,   MAXSTOKES, &
                           MAX_DIRECTIONS, MAXTASKS, MAXWVN ) )

!     Initialize output

        OMP_N_GEOMETRIES = 0
        STOKES_PT = ZERO
        MEANST_PT = ZERO
        FLUX_PT   = ZERO

!     Begin parallel region
!   Rob 3/22/15. Add Output Level variables to Shared Space

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     CMASK, GMASK, RMASK, SMASK, N_USER_LEVELS, USER_LEVELS, &
!$OMP     NSTOKES, N_SZANGLES, NUMWVN, NTASKS, NLAYERS, NGREEKMAT_ENTRIES, &
!$OMP     HEIGHT_GRID, MOLEXT, MOLOMG, PROBLEM_RAY, &
!$OMP     AERMOMS, NGREEK_MOMENTS_INPUT, &
!$OMP     VLIDORT_FixIn_SAVE, VLIDORT_ModIn_SAVE, VLIDORT_Sup_SAVE, &
!$OMP     OMP_N_GEOMETRIES, STOKES_PT, MEANST_PT, FLUX_PT, &
!$OMP     VLIDORT_TID_SAVE ) &
!$OMP   FIRSTPRIVATE (OPENFILEFLAG)

!     Obtain thread number

      TID = OMP_GET_THREAD_NUM()
      !TID = 0
      
!     Obtain and display total number of threads and local thread

      IF (TID == 0) THEN
        OMP_NTHREADS = OMP_GET_NUM_THREADS()
        !OMP_NTHREADS = 1
        !write(*,*)
        !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
        write(*,*)
        write(*,'(1x,a,i1,a)') 'Running vlidort driver with ',OMP_NTHREADS,' thread(s)'
        write(*,'(1x,a,i1/)')  ' *** Tests for NSTOKES = ',nstokes
      END IF

      write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      do wvn = 1, numwvn

      !write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn

!  Initialize some inputs

      VLIDORT_FixIn = VLIDORT_FixIn_SAVE
      VLIDORT_ModIn = VLIDORT_ModIn_SAVE
      VLIDORT_Sup   = VLIDORT_Sup_SAVE

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero
      !Fmatrix_up           = zero ! 2p8
      !Fmatrix_dn           = zero ! 2p8

!  Rob 3/22/15. Add wavelength Diagnostic (just a dummy here)
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Fix the non-aerosol layers

      n6 = nlayers - 6
      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
          DO L = 0, 2
! Old       GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(ck,L)
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
          GK = gmask(k); sk = dble(smask(k)); ck = cmask(k) ; rk = rmask(k)
          DO L = 0, 2
! Old       MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
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

!  Surface

      lambertian_albedo = 0.05d0

!  Copy local control integers

      VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo

      !VLIDORT_FixIn%Optical%TS_FMATRIX_UP = fmatrix_up ! zero here
      !VLIDORT_FixIn%Optical%TS_FMATRIX_DN = fmatrix_dn ! zero here

!  Start task loop.
!   --- Version 2.7 (3/22/15) increase to 8 tasks (rename from tasks)
!   --- Version 2.8 (9/18/16) reduced  to 6 tasks

      do task = 1, ntasks

        t = task

        !write(*,*) 'doing LID task  = ',task

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
           DO_FOCORR_NADIR    = .false.
           DO_FOCORR_OUTGOING = .true.
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
        !write(*,*)'Doing task # ', task

!    Rob 3/22/15/. Add FO_CALC setting, now named DO_FOCORR (V2.8)

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

        write(charTID,'(i1)') TID
        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution_TID_' // charTID // '.log', &
          VLIDORT_ERRUNIT+TID, OPENFILEFLAG, VLIDORT_Out%Status )

!  Close error file if it has been opened

        IF ( OPENFILEFLAG ) then
          CLOSE(VLIDORT_ERRUNIT+TID)
          write(*,*) ; write(*,*) 'For OMP thread ',TID
          write(*,*)'V2p8p3_solar_OMP_tester.exe: program executed with internal '// &
           'defaults, warnings in "V2p8p3_VLIDORT_Execution_TID_' // charTID // '.log"'
        ENDIF

!  Set some variables

        n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Save results for the current wavenumber
!    -- Version 2.8, don't need the utamap.

        do uta = 1, n_user_levels
          do o1 = 1, nstokes
            do v = 1, n_geometries
               stokes_PT(uta,v,o1,upidx,t,wvn) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
               stokes_PT(uta,v,o1,dnidx,t,wvn) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            enddo
            do v = 1, n_szangles
               MEANST_PT(uta,v,o1,upidx,t,wvn) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
               MEANST_PT(uta,v,o1,dnidx,t,wvn) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                               + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
               FLUX_PT(uta,v,o1,upidx,t,wvn)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
               FLUX_PT(uta,v,o1,dnidx,t,wvn)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                               + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)
            enddo
          enddo
        enddo

!  Save number of geometries & OMP thread information

        OMP_N_GEOMETRIES(task,wvn) = n_geometries
        VLIDORT_TID_SAVE(task,wvn) = TID + 1

!  End task loop

      enddo

!  End "wavenumber" loop

      enddo

!$OMP END DO

      !IF (TID == 0) call cpu_time(omp_e2)

!  End parallel region

!$OMP END PARALLEL


!  Timing tests (OpenMP)
!  =====================

      call cpu_time(omp_e2)

      if (omp_nthreads <= n_core) then
        time_divider = omp_nthreads
      else
        time_divider = n_core
      endif

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        numwvn, omp_nthreads, (omp_e2 - omp_e1)/real(time_divider)

!  Write Final file for all threads
!  ================================

!  Set some variables

      n_geometries  = OMP_N_GEOMETRIES(1,1)

!  Prepare output file

      OutUnit = 36

      if ( do_std_output_file ) then
        !place all output in one file
        write(ch1,'(i1)')   omp_nthreads
        write(ch2,'(i5.5)') numwvn 
        tail = '_nt' // ch1 // '_nwn' // ch2

        OPEN(OutUnit,file='vlidort_v_test/results_solar_tester_' &
             // c2 // '.all' // trim(tail),status='replace')
      endif

!  Write output. (Just a sampling of Results)

      !do wvn = 1, numwvn
      do wvn = 1, numwvn, numwvn/omp_nthreads

        if ( .not. do_std_output_file ) then
          !place one set of outputs per file
          write(ch1,'(i1)')   omp_nthreads
          write(ch2,'(i5.5)') numwvn 
          write(ch3,'(i5.5)') wvn
          tail = '_nt' // ch1 // '_nwn' // ch2 // '_wvn' // ch3

          OPEN(OutUnit,file='vlidort_v_test/results_solar_tester_' &
               // c2 // '.all' // trim(tail),status='replace')
        endif

!  Write Stokes Intensity and Mean values

        call write_stokes &
           ( OutUnit, 1, 'I', wvn, VLIDORT_TID_SAVE, DO_STD_OUTPUT_FILE, &
             maxtasks, n_geometries, n_szangles, n_user_levels, maxwvn, &
             user_levels, Stokes_PT(:,:,:,:,:,wvn), &
             MEANST_PT(:,:,:,:,:,wvn), FLUX_PT(:,:,:,:,:,wvn) )

!  Write Stokes Q and U and Mean values

        if ( nstokes .ge. 3 ) then
          call write_stokes &
           ( OutUnit, 2, 'Q', wvn, VLIDORT_TID_SAVE, DO_STD_OUTPUT_FILE, &
             maxtasks, n_geometries, n_szangles, n_user_levels, maxwvn, &
             user_levels, Stokes_PT(:,:,:,:,:,wvn), &
             MEANST_PT(:,:,:,:,:,wvn), FLUX_PT(:,:,:,:,:,wvn) )

          call write_stokes &
           ( OutUnit, 3, 'U', wvn, VLIDORT_TID_SAVE, DO_STD_OUTPUT_FILE, &
             maxtasks, n_geometries, n_szangles, n_user_levels, maxwvn, &
             user_levels, Stokes_PT(:,:,:,:,:,wvn), &
             MEANST_PT(:,:,:,:,:,wvn), FLUX_PT(:,:,:,:,:,wvn) )
        endif

!  Write Stokes V and Mean values

        if ( nstokes .eq. 4 ) then
          call write_stokes &
           ( OutUnit, 4, 'V', wvn, VLIDORT_TID_SAVE, DO_STD_OUTPUT_FILE, &
             maxtasks, n_geometries, n_szangles, n_user_levels, maxwvn, &
             user_levels, Stokes_PT(:,:,:,:,:,wvn), &
             MEANST_PT(:,:,:,:,:,wvn), FLUX_PT(:,:,:,:,:,wvn) )
        endif

        if ( .not. do_std_output_file ) close(OutUnit)

!  End file write loop

      enddo

      if ( do_std_output_file ) close(OutUnit)

!     Deallocate output memory

      DEALLOCATE ( OMP_N_GEOMETRIES, STOKES_PT, MEANST_PT, FLUX_PT )

!  End test loop

      enddo

!  Finish

      write(*,*) '******************************************'
      write(*,*)
      write(*,*)'Main program finished successfully'

      stop

      end program Solar_Tester

!  ****************************************************************************
!  ****************************************************************************
      subroutine write_stokes &
        ( un, o1, a1, wvn, vlidort_tid_save, do_std_output_file, &
          maxtasks, n_geometries, n_szangles, n_user_levels, maxwvn, &
          user_levels, STOKES_PT, MEANST_PT, FLUX_PT )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, &
                                 maxstokes, max_directions, upidx, dnidx

      implicit none

!  input variables

      integer, intent(in) ::           maxtasks, n_geometries, n_szangles, n_user_levels, &
                                       maxwvn
      integer, intent(in) ::           un, o1, wvn
      character (len=1), intent(in) :: a1
      double precision, intent(in)  :: user_levels ( max_user_levels )

      integer, intent(in)           :: vlidort_tid_save ( maxtasks, maxwvn )
      logical, intent(in)           :: do_std_output_file

      DOUBLE PRECISION, INTENT(IN)  :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAXTASKS )
      DOUBLE PRECISION, INTENT(IN)  :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAXTASKS )
      DOUBLE PRECISION, INTENT(IN)  :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAXTASKS )

!  local variables

      integer :: n, v, t, ntasks

      ntasks = maxtasks

!  Write stokes component

      if ( do_std_output_file ) then
        write(un,'(/T32,a,I4/T32,a/)') &
          '  STOKES-'//a1//' , tasks 1-6, WVN Point # ',wvn, &
          '  ======================================'
        write(un,'(a,T32,6(8x,I2,8x)//a,T32,6(a16,2x)/)')&
          'OPENMP_Thread # = ',(vlidort_tid_save(t,wvn),t=1,ntasks), &
                                    'Geometry    Level/Output', &
                                    'No FOCorr No DM','No FOCorr + DM ', &
                                    'FO_ReglrPS + DM','FO_EnhncPs + DM', &
                                    '#4 + Solsaving ','#5 + BVPTelscpe'
! 2p7          'No SSCorr No DM','No SSCorr + DM ','SSNadir  + DM  ',&
! 2p7         'SSOutgoing + DM','FO Reg PS + DM ','FO Enh PS + DM ','#5 + Solsaving ','#6 + BVPTelscpe'
      else
        write(un,'(/T32,a/T32,a/)') &
          '  STOKES-'//a1//' , tasks 1-6', &
          '  ======================'
        write(un,'(a,T32,6(a16,2x)/)')&
                                    'Geometry    Level/Output', &
                                    'No FOCorr No DM','No FOCorr + DM ', &
                                    'FO_ReglrPS + DM','FO_EnhncPs + DM', &
                                    '#4 + Solsaving ','#5 + BVPTelscpe'
! 2p7          'No SSCorr No DM','No SSCorr + DM ', 'SSNadir  + DM  ',&
! 2p7         'SSOutgoing + DM','FO Reg PS + DM ','FO Enh PS + DM ','#5 + Solsaving ','#6 + BVPTelscpe'
      endif

      do v = 1, n_geometries
        do n = 1, n_user_levels
          write(un,366)v,'Upwelling @',user_levels(n), (stokes_PT(n,v,o1,upidx,t),t=1,ntasks)
        enddo
        do n = 1, n_user_levels
          write(un,366)v,'Dnwelling @',user_levels(n), (stokes_PT(n,v,o1,dnidx,t),t=1,ntasks)
        enddo
        write(un,*)' '
      enddo

!  Write fluxes

      if ( do_std_output_file ) then
        write(un,'(/T32,a,I4/T32,a/)') &
          'Stokes-'//a1//' ACTINIC + REGULAR FLUXES, tasks 1-2, WVN Point # ',wvn, &
          '============================================================='
        write(un,'(a,T31,4(8x,I2,8x)//a,T32,4(a16,2x)/)')&
          'OPENMP_Thread # = ',(vlidort_tid_save(t,wvn),t=1,2),(vlidort_tid_save(t,wvn),t=1,2),&
          ' Sun SZA    Level/Output', &
          ' Actinic No DM ',' Actinic + DM  ',&
          ' Regular No DM ',' Regular + DM  '
      else
        write(un,'(/T32,a/T32,a/)') &
          'Stokes-'//a1//' ACTINIC + REGULAR FLUXES, tasks 1-2', &
          '=============================================='
        write(un,'(a,T31,4(a16,2x)/)')' Sun SZA    Level/Output', &
                                      ' Actinic No DM ',' Actinic + DM  ', &
                                      ' Regular No DM ',' Regular + DM  '
      endif

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
