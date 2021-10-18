
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

     PROGRAM Siewert2000_validation

!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

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

!  local dimensions

      INTEGER, PARAMETER :: MAXMOM = 100, MAXLAY = 20

!  help variables

      INTEGER ::          UM,UA,L,N,V,UT,O1,NSTOKES_SQ
      DOUBLE PRECISION :: DIFF, TAUTOTAL, TAUDIFF, TOPHEIGHT

!  VLIDORT standard input preparation

      INTEGER ::          NLAYERS
      INTEGER ::          NSTOKES
      INTEGER ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION :: HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Displaying results

      INTEGER ::          N_USER_VZANGLES
      INTEGER ::          N_USER_LEVELS
      INTEGER ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION :: USER_STREAMS ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: STOKES ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  aerosol data problem 2a

      DOUBLE PRECISION :: PROBLEM_IIA(6,0:11)

      DATA (PROBLEM_IIA(1,L),L = 0,11)/ &
         0.000000D0,  0.000000D0,  3.726079D0,  2.202868D0, &
         1.190694D0,  0.391203D0,  0.105556D0,  0.020484D0, &
         0.003097D0,  0.000366D0,  0.000035D0,  0.000003D0 /
      DATA (PROBLEM_IIA(2,L),L = 0,11)/ &
         1.000000D0,  2.104031D0,  2.095158D0,  1.414939D0, &
         0.703593D0,  0.235001D0,  0.064039D0,  0.012837D0, &
         0.002010D0,  0.000246D0,  0.000024D0,  0.000002D0 /
      DATA (PROBLEM_IIA(3,L),L = 0,11)/ &
         0.000000D0,  0.000000D0, -0.116688D0, -0.209370D0, &
        -0.227137D0, -0.144524D0, -0.052640D0, -0.012400D0, &
        -0.002093D0, -0.000267D0, -0.000027D0, -0.000002D0 /
      DATA (PROBLEM_IIA(4,L),L = 0,11)/ &
         0.915207D0,  2.095727D0,  2.008624D0,  1.436545D0, &
         0.706244D0,  0.238475D0,  0.056448D0,  0.009703D0, &
         0.001267D0,  0.000130D0,  0.000011D0,  0.000001D0 /
      DATA (PROBLEM_IIA(5,L),L = 0,11)/ &
         0.000000D0,  0.000000D0,  0.065456D0,  0.221658D0, &
         0.097752D0,  0.052458D0,  0.009239D0,  0.001411D0, &
         0.000133D0,  0.000011D0,  0.000001D0,  0.000000D0 /
      DATA (PROBLEM_IIA(6,L),L = 0,11)/ &
         0.000000D0,  0.000000D0,  3.615946D0,  2.240516D0, &
         1.139473D0,  0.365605D0,  0.082779D0,  0.013649D0, &
         0.001721D0,  0.000172D0,  0.000014D0,  0.000001D0 /

!  Test 2
!  ======

!  description
!     Problem IIA Slab, Siewert (2000)
!     Single layer of optical thickness 1, Lambertian albedo 0.0
!     Aerosol scattering problem IIA, SSA = 0.973257
!     plane parallel, classical solution
!     1 solar angles, 11 viewing angles, 3 azimuths (33 geometries)
!     Upwelling/downwelling radiance at 7 Tau values
!     WITH POLARIZATION

      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_Siewert2000_validation.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

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

!  Siewert value

      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1) = DACOS(0.6d0)/DEG_TO_RAD

!  Define some local variables

      NSTOKES    = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTOKES_SQ = NSTOKES * NSTOKES

!  Initialize Greek matrix

      DO N = 1, MAXLAYERS
        DO L = 0, MAXMOMENTS
          DO O1 = 1, NSTOKES_SQ
            GREEKMAT_TOTAL_INPUT(L,N,O1) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Set up of optical properties

      deltau_vert_input = zero
      omega_total_input = zero
      greekmat_total_input = zero

      height_grid = zero
      nlayers              = VLIDORT_FixIn%Cont%TS_NLAYERS
      ngreek_moments_input = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

      topheight = 10.0d0
      tautotal  = 1.0d0
      diff = dble(nlayers)
      height_grid(0) = topheight
      taudiff = tautotal / diff
      do n = 1, nlayers
       height_grid(n) = topheight - dble(n)*topheight/diff
       deltau_vert_input(n) = taudiff
       omega_total_input(n) = 0.973527d0
       do l = 0, ngreek_moments_input
        greekmat_total_input(l,n,1)  = problem_iia(2,l)
        if ( l.eq.0) greekmat_total_input(l,n,1) = 1.0d0
        if ( nstokes .gt. 1 ) then
          greekmat_total_input(l,n,2)  =   problem_iia(3,l)
          greekmat_total_input(l,n,5)  =   problem_iia(3,l)
          greekmat_total_input(l,n,6)  =   problem_iia(1,l)
          greekmat_total_input(l,n,11) =   problem_iia(6,l)
          greekmat_total_input(l,n,12) = - problem_iia(5,l)
          greekmat_total_input(l,n,15) =   problem_iia(5,l)
          greekmat_total_input(l,n,16) =   problem_iia(4,l)
        endif
       enddo
      enddo

!  Copy to type-structure inputs

      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

!  Call to VLIDORT_MASTER

      write(*,'(/a)')'Siewert2000 validation status:'
      do_debug_input = .false.
      CALL VLIDORT_MASTER ( do_debug_input, &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup,   &
        VLIDORT_Out )

!  Exception handling, write-up (optional)

      OPENFILEFLAG = .false.
      CALL VLIDORT_WRITE_STATUS ( &
        'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Write results

      N_USER_VZANGLES = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      USER_STREAMS    = COS(VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT * DEG_TO_RAD)
      N_USER_LEVELS   = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
      VZA_OFFSETS     = VLIDORT_Out%Main%TS_VZA_OFFSETS
      STOKES          = VLIDORT_Out%Main%TS_STOKES

      OPEN(1,FILE='vlidort_v_test/results_Siewert2000_validation.all', STATUS='unknown')

      UA = 1

      write(1,*)'Table 2'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      write(1,*)'Table 3'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      UA = 2

      write(1,*)'Table 4'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      write(1,*)'Table 5'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      write(1,*)'Table 6'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,3,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,3,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      write(1,*)'Table 7'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,4,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,4,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      UA = 3

      write(1,*)'Table 8'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      write(1,*)'Table 9'
      DO UM = 1, N_USER_VZANGLES
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') -USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
      ENDDO
      DO UM = N_USER_VZANGLES, 1, -1
        V = VZA_OFFSETS(1,UM) + UA
        write(1,'(f8.2,1p7e15.5)') USER_STREAMS(UM), &
            (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
      ENDDO

      CLOSE(1)

      write(*,*)
      write(*,*)'Main program finished successfully'

      stop

      END PROGRAM Siewert2000_validation
