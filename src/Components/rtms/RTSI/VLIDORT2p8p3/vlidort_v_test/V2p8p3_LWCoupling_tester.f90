
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
                                                   
program LWCoupling_tester

!  This 2.8.1 VLIDORT Driver tests the Water-leaving adjustment facility
!    PROGRAMMED BY ROBERT SPURR, 4/10/19, 5/10/19
  
!  Used Module files

   USE VSLEAVE_SUP_AUX_m, Only : VSLEAVE_READ_ERROR
   USE VSLEAVE_SUP_MOD_m

   USE VLIDORT_PARS_m
   USE VLIDORT_IO_DEFS_m

   USE VLIDORT_AUX_m     , only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
   USE VLIDORT_INPUTS_m  , only : VLIDORT_INPUT_MASTER

   USE VLIDORT_MASTERS_m

   USE VLIDORT_VSLEAVE_SUP_ACCESSORIES_m

!  Implicit none

   IMPLICIT NONE

!  VLIDORT file inputs status structure

   TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT debug input control

      LOGICAL :: DO_DEBUG_INPUT

!  VLIDORT input structures

   TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

!  LIDORT supplements i/o structure

   TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

!  VLIDORT output structure

   TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VSLEAVE supplement file inputs status structure

   TYPE(VSLEAVE_Input_Exception_Handling) :: VSLEAVE_Sup_InputStatus

!  VSLEAVE supplement input structures

   TYPE(VSLEAVE_Sup_Inputs)               :: VSLEAVE_Sup_In

!  VSLEAVE supplement output structure

   TYPE(VSLEAVE_Sup_Outputs)              :: VSLEAVE_Sup_Out
   TYPE(VSLEAVE_Output_Exception_Handling):: VSLEAVE_Sup_OutputStatus

!  VSLEAVE supplement / VLIDORT VSLEAVE-related inputs consistency check status

   TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_SLEAVECheck_Status

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  Help variables
!   Rob 3/22/15. renamed thread --> task

      INTEGER ::          NGREEKMAT_ENTRIES, LDUM, LCHECK
      INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)
      INTEGER ::          N6,NDUM,T,NTASKS,TASK,O1
      integer ::          nbeams, nvzas, nazms
      integer ::          i, l, n, v, vd, w, uta, nmi, ib, um, ia

      DOUBLE PRECISION :: AMOLSCA, TOTSCA, TOTEXT,  RAYWT, AERWT
      DOUBLE PRECISION :: MOLOMG ( MAXLAYERS ), MOLEXT ( MAXLAYERS )
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), DEPOL, BETA2
      DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM, KD, LAMBDA

      INTEGER ::          OutUnit

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_OBSERVATION_GEOMETRY
      LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      INTEGER ::          NLAYERS, NSTOKES, NFINELAYERS
      CHARACTER*90 ::     TRACE
      
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

      INTEGER ::          N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  VLIDORT standard output preparation

      INTEGER ::          N_GEOMETRIES
      INTEGER ::          N_SZANGLES

!  Saved results

      DOUBLE PRECISION :: I_base_TOA(2,max_geometries)
      DOUBLE PRECISION :: I_base_BOA(2,max_geometries)

!  Timing test

      INTEGER :: irun
      REAL    :: e1,e2

!  Start of code
!  =============

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

   GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
   RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!   CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
   CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
   SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Wavelength

   lambda = 440.0d0

   !  Define a proxy

   nmi = 2

!  Initialize error file output flag

   OPENFILEFLAG = .false.

!  START PARTS LOOP
!  ================

!   w = 1 : PART I.  Do Full Iterative LW calculation and Output WLADJUSTED surface-leaving
!   w = 2 : PART II. Set External_Wleave flag and use External data

   do w = 1, 2

!  Read SLEAVE inputs

      CALL VSLEAVE_INPUTMASTER ( &
        'vlidort_v_test/V2p8p3_VSLEAVE_LWCoupling.cfg', & ! Input
         VSLEAVE_Sup_In,         & ! Outputs
         VSLEAVE_Sup_InputStatus ) ! Outputs

      IF ( VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VSLEAVE_READ_ERROR ( 'V2p8p3_VSLEAVE_LWCoupling.log', VSLEAVE_Sup_InputStatus )

!  Define some inputs not handled by VSLEAVE_LIN_INPUTMASTER

      VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = 'vlidort_v_test/data'

!  Redefine SLEAVE input wavelength (override SLEAVE config file)

      VSLEAVE_Sup_In%SL_FL_WAVELENGTH = lambda !in nm

!  Save some inputs to local variables

      DO_OBSERVATION_GEOMETRY = VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS

!  Call and exception handling
      
      CALL VSLEAVE_MAINMASTER ( &
        VSLEAVE_Sup_In,          & ! Inputs
        VSLEAVE_Sup_Out,         & ! Outputs
        VSLEAVE_Sup_OutputStatus ) ! Output Status

      if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.vlidort_success ) THEN
         TRACE = 'VSLEAVE_MAINMASTER failed' ; go to 678
      endif

!  Read LIDORT Main inputs
! 
      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_LWCoupling.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_LWCoupling.log', VLIDORT_InputStatus )

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

!  Redefine some VLIDORT input variables (override VLIDORT config file)

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = lambda/1000.0d0 !in um

!  Proxies

      nbeams = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      nvzas  = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      nazms  = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      nstokes       = VLIDORT_FixIn%Cont%TS_nstokes
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels(1:n_user_levels) = VLIDORT_ModIn%MUserVal%TS_user_levels(1:n_user_levels)

!  IMPORTANT - Check compatibility of VSLEAVE and VLIDORT Main inputs

      CALL VLIDORT_VSLEAVE_INPUT_CHECK ( &
        VSLEAVE_Sup_In,             & ! Inputs
        VLIDORT_FixIn,              & ! Inputs
        VLIDORT_ModIn,              & ! Inputs
        VLIDORT_SLEAVECheck_Status )  ! Outputs

!  Exception handling

      IF ( VLIDORT_SLEAVECheck_Status%TS_STATUS_INPUTCHECK .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_VSLEAVE_INPUT_CHECK_ERROR ( &
           'V2p8p3_VLIDORT_VSLEAVE_checkcall.log', VLIDORT_SLEAVECheck_Status )

!  Copy VSLEAVE Sup outputs to VLIDORT's VSLEAVE Sup inputs (std only)

      CALL SET_VLIDORT_VSLEAVE_INPUTS ( &
        VSLEAVE_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn,  & !Inputs
        VLIDORT_Sup )                                     !Outputs

!  Set the two calls
      
      if ( w.eq.1 ) then
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE  = .false.
         VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = .true.
         VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION      = .true.
      else if ( w.eq.2 ) then
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE  = .true.
         VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = .false.
         VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION      = .false.
      endif

!  External setting

      if ( w.eq.2 ) then
         write(*,*)
         write(*,*)'PART II - Reading external field (PART I WLAdjusted data)'
         open(1,file='vlidort_v_test/Sleave_Isotropic_WLAdjusted.dat',status='old')
         do v = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
            read(1,*)vd, VLIDORT_Sup%Sleave%TS_SLTERM_ISOTROPIC(1,v)
         enddo
         close(1)
      else
         write(*,*)'PART I  - Iterative WL calculation'
      endif

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Get the pre-prepared atmosphere
!  -------------------------------

      height_grid = 0.0d0
      nlayers = 23 ; VLIDORT_FixIn%Cont%TS_nlayers = nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

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

!  Initialise

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero
      Fmatrix_up           = zero ! 2p8
      Fmatrix_dn           = zero ! 2p8

!  Fix optical properties all layers

      do n = 1, nlayers
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
          DO L = 0, 2
            GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
       ENDDO
      enddo

      lambertian_albedo = 0.2d0

!  Copy to type-structure inputs

      VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = 2
      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP = fmatrix_up ! zero here
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN = fmatrix_dn ! zero here

!   VLIDORT call

      call cpu_time(e1)
      do_debug_input = .false.
      do irun = 1, 10
         if (mod(irun,10).eq.0) write(*,*)'doing run # ', irun
         CALL VLIDORT_MASTER ( do_debug_input,  &
             VLIDORT_FixIn, &
             VLIDORT_ModIn, &
             VLIDORT_Sup,   &
             VLIDORT_Out )
      enddo
      call cpu_time(e2)
      write(*,*)' CPU time ', e2-e1
      
!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS (  &
        'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Save baseline results

      do v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
        I_base_TOA(w,v) = VLIDORT_Out%Main%TS_Stokes(1,v,1,upidx)
        I_base_BOA(w,v) = VLIDORT_Out%Main%TS_Stokes(2,v,1,dnidx)
      enddo

!  write Sleave results
!  write Sleave results ( PART I )
!    - Write Unadjusted Sleave output (Direct from SLEAVE supplement)
!    - Write WLadjusted Sleave output (from LIDORT Transflux Adjustment)

      if ( w.eq.1 ) then
         open(1,file='vlidort_v_test/Sleave_Isotropic_UNadjusted.dat',status='unknown')
         open(2,file='vlidort_v_test/Sleave_Isotropic_WLAdjusted.dat',status='unknown')
         do v = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
            write(1,300) v, VLIDORT_Sup%Sleave%TS_SLTERM_ISOTROPIC(1,v)
            write(2,300) v, VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(1,v)
         enddo
         close(1) ; close(2) ; close(3)
      endif

300   format(1x,i2.2,1x,1pe16.7)

!  End two-part loop
      
   enddo
   
!  write results
   
   open(1,file = 'vlidort_v_test/results_LWCoupling_TOAUp_BOADn_Radiances.all',status = 'unknown')
   write(1,'(/A/A//A/A/)')     '                   Radiances with Isotropic  Water-leaving scenario',&
                               '                   ================================================',&
                               '                                    TOA Upwelling                BOA Downwelling',&
                               'Geom     SZA/VZA/AZM            Adjusted      External        Adjusted      External'
   do ib = 1, nbeams
      do um = 1, nvzas
         do ia = 1, nazms
            v = nazms*nvzas*(ib-1) + nazms*(um-1) + ia
            write(1,67) v, VLIDORT_ModIn%MSunrays%TS_SZANGLES(ib),&
                           VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(um),&
                           VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(ia), I_base_TOA(1:2,v), I_base_BOA(1:2,v)
   enddo ; enddo ; enddo
67 format(i2,2x,3f7.2,1x,2(2x,1p2e14.6))
   close(1)
   
!  Finish

   write(*,*)
   write(*,*) 'Main program finished successfully'
   stop

!mick fix 3/22/2017 - added error finish section
!  Error finish

678   continue

      write(*,*)
      write(*,'(1x,a)') trim(TRACE)
      write(*,*)'Number of error messages from SLEAVE calculation = ',VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
      do i = 1, VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
         write(*,'(a,i2,a,a)')' Message # ', i, ': ',adjustl(trim(VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)))
      enddo

   stop
end program LWCoupling_Tester

