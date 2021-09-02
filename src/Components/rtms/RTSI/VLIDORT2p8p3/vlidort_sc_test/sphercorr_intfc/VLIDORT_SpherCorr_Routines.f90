
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

MODULE VLIDORT_SPHERCORR_ROUTINES_m

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3-point and multi-point Corrections

!  3-point correction history --
!    V1: 03/01/20. New 3pt Correction

!  Multi-point correction history --
!    V1: 03/01/20. New  Correction

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Package
!     R. Spurr. RT Solutions Inc.

!  HERE IS THE VLIDORT IMPLEMENTATION
!  ==================================

!  Module files for VLIDORT.

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_MASTERS_m
 
!  Module for the geometry conversions

      USE SPHERCORR_GEOMETRY_CONVERSIONS_m

!  All three routines are public

public  :: MS2pt_SPHERICITY_CORRECTION_V3, &
           MS3pt_SPHERICITY_CORRECTION_V1, &
           MSMpt_SPHERICITY_CORRECTION_V1

contains

subroutine MS2pt_SPHERICITY_CORRECTION_V3 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(OUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(OUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector with 2-point corr
      DOUBLE PRECISION, INTENT(OUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector for BOA geo.     4 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector with MS2pt corr. 5 = 1 + 3

!  Twilight flag and exception handling

      LOGICAL      , intent(out) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(out) :: Local_Message, Local_Action

!  Local
!  -----

!  Variables for doing the sphericity calculation

      DOUBLE PRECISION :: MU, MU0, MU1, D01, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      DOUBLE PRECISION :: EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA(3), OBSGEOMS_TOA(3), COSSCAT
      DOUBLE PRECISION :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NLAYERS, NS, O1, DIRIDX

!  Start Code
!  ==========

!  Initialize

      FO_STOKES_BOAGEOM   = zero
      MS_STOKES_BOAGEOM   = zero
      MS_STOKES_SPHERCORR = zero
      STOKES_BOAGEOM      = zero
      STOKES_SPHERCORR    = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Local_Message = ' '
      Local_Action  = ' '

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 2-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 2
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]

!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set second geometry through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Conversion. 12/18/19. Add Direction index.

      Call BOATOA_2point_conversion ( DIRIDX, &
         EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3), & ! input TOA/BOA heights, BOA Geometry
         cosscat, OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3) )                   ! output TOA geometry, scattering angle

!  Debug
!write(*,*)DIRIDX,OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3)
!write(*,*)DIRIDX,OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3)
!write(*,*)DIRIDX,ACOS(COSSCAT)/DEG_TO_RAD
!stop

!  Reset the VLIDORT second geometrical input (TOA)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_TOA(3)

!  re-set geometry numbers. 2 Geometries

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 2
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = 2
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 2
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = 2

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_TOA(1).gt.90.0 ) then
         Twilight_flag = .true. ; return
      endif

!  debug - Check this against the FO geometry value.
!WRITE(*,*)VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3) 
!WRITE(*,*)VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) 

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling (simply done here)

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS2pt_SPHERICITY_CORRECTION_V3'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store Stokes-vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS) = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS) = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  Stokes-vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,2,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,2,1:NS) 
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,2,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,2,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,2,1:NS,DNIDX)
      endif

!write(*,*)FO_STOKES_BOAGEOM(1),FO_STOKES_TOAGEOM(1), VLIDORT_Out%Main%TS_STOKES(1,1:2,1,DNIDX) ; stop

!  perform 2-POINT SPHERICITY CORRECTION, Radiative Transfer
!  =========================================================

!  set the surface MS source term

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  get the Layer MS source terms by mid-point interpolation with Cos(SZA)
!    -- Works, except for the twilight and azimuth-flip conditions
!      Mu0 = COS(OBSGEOMS_BOA(1)*DEG_TO_RAD)
!      Mu1 = COS(OBSGEOMS_TOA(1)*DEG_TO_RAD)
!      D01 = ONE / ( Mu1 - Mu0 )
!      do N = 1, nlayers
!         Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N-1)) )
!         F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
!         DO O1 = 1, NS
!            AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,O1,N)
!write(33,*)N,AMSST(N,O1),VLIDORT_Out%Main%TS_PATHGEOMS(1,N)/DEG_TO_RAD
!         ENDDO
!      ENDDO

!  Alternative: get the Layer MS source terms by end-point interpolation with Cos(VZA)
!     -- Always works, however may not be accurate around azimuth-flip conditions

      Mu0 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
      Mu1 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
      D01 = ONE / ( Mu1 - Mu0 )
      do N = 1, nlayers
         Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
         F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
         DO O1 = 1, NS
           AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,O1,N)
         ENDDO
!write(*,*)N,F0,VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1,N),VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1,N)
      ENDDO

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final Radiance answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MS2pt_SPHERICITY_CORRECTION_V3

!

subroutine MS3pt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations

      REAL(fpk), INTENT(OUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(OUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(OUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector vector with 3-point corr
      REAL(fpk), INTENT(OUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector vector for BOA geo.     4 = 1 + 2.
      REAL(fpk), INTENT(OUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector vector with MS3pt corr. 5 = 1 + 3.

!  Twilight flag and exception handling

      LOGICAL      , intent(out) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(out) :: Local_Message, Local_Action

!  Local
!  -----

!  Variables for doing the sphericity calculation

      REAL(fpk) :: MU, MU0, MU1, MU2, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      REAL(fpk) :: EARTH_RADIUS, TOA_HEIGHT, MID_HEIGHT, OBSGEOMS_BOA(3), OBSGEOMS_MID(3), OBSGEOMS_TOA(3), COSSCAT
      REAL(fpk) :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NLAYERS, NS, DIRIDX

!  Start Code
!  ==========

!  Initialize output.

      FO_STOKES_BOAGEOM   = zero
      MS_STOKES_BOAGEOM   = zero
      MS_STOKES_SPHERCORR = zero
      STOKES_BOAGEOM      = zero
      STOKES_SPHERCORR    = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Local_Message = ' '
      Local_Action  = ' '

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 3-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set other geometries through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Geometrical Conversion. 3/1/20. Add MID_HEIGHT value

      Call BOATOA_3point_conversion ( DIRIDX, &
         EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA,         & ! input TOA height, BOA Geometry
         MID_HEIGHT, cosscat, OBSGEOMS_TOA, OBSGEOMS_MID ) ! output TOA/MID geometries, scattering angle, mid-height

!  Reset the VLIDORT second (mid) and third (TOA) geometrical input

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_MID(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_MID(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_MID(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_MID(3)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(3,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(3)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(3)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(3)            = OBSGEOMS_TOA(3)

!  re-set geometry numbers. 3 Geometries

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 3
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = 3
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 3
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = 3

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_TOA(1).gt.90.0 ) then
         Twilight_flag = .true. ; return
      endif

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS3pt_SPHERICITY_CORRECTION_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store STOKES Vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  STOKES Vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,3,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,3,1:NS) 
         STOKES_TOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,3,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,3,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,3,1:NS,DNIDX)
      endif

!  perform 3-point SPHERICITY CORRECTION, Recurrence Relation
!  ==========================================================

!  set the surface MS source term, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS) = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  Alternative: get the Layer MS source terms by 3-point interpolation with Cos(VZA) using exponential relaxation
!     -- Always works, however may not be accurate around azimuth-flip conditions
!     -- THIS DOES NOT WORK GOOD AND WAS ABANDONED................................
!      Mu2 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
!      Mu1 = COS(OBSGEOMS_MID(2)*DEG_TO_RAD) ! Average value
!      Mu0 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
!     write(*,*)'mu2,mu1,mu0',mu2, mu1, mu0
!     write(397,'(I2,f10.6)')0,COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,0))
!      D21 = ONE / (Mu2 - Mu1)
!      do N = 1, nlayers
!         Mu = COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N))
!            ZTOP = ( VLIDORT_Out%Main%TS_LAYER_MSSTS(3,1:NS,N) - VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) )
!            ZBOT = ( VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) - VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N) )
!            Z = ZTOP/ZBOT ; KAY = Log(Z) * D21
!            BCON = ZBOT * EXP ( KAY * MU1 )  / ( Z - ONE )
!            ACON = VLIDORT_Out%Main%TS_LAYER_MSSTS(2,N) - BCON * EXP ( - KAY * MU1 )
!            AMSST(N,1:NS) = ACON + BCON * EXP ( - KAY * MU )
!        write(397,'(i2,f10.6,1p3e16.6,2x,1pe16.6)')&
!N,COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)),VLIDORT_Out%Main%TS_LAYER_MSSTS(1:3,1:NS,N), AMSST(N,1:NS)
!        ENDDO

!  Alternative: get the Layer MS source terms by 3-point Linear interpolation with Cos(VZA) 
!     -- Always works, however may not be accurate around azimuth-flip conditions
!     -- in any layer, Mu is the layer average cosine, and we interpolate to this value.

      Mu2 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
      Mu1 = COS(OBSGEOMS_MID(2)*DEG_TO_RAD) ! Average value of MU0 and MU2
      Mu0 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
      do N = 1, nlayers
        Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
        IF ( Mu.gt.mu1 ) then
          F0 = ( Mu1 - Mu ) / ( Mu1 - Mu0 ) ; F1 = ONE - F0
          AMSST(N,1:NS) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(3,1:NS,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N)
        ELSE
          F0 = ( Mu2 - Mu ) / ( Mu2 - Mu1 ) ; F1 = ONE - F0
          AMSST(N,1:NS) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N)
        ENDIF
      enddo

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final Radiance answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MS3pt_SPHERICITY_CORRECTION_V1

!

subroutine MSMpt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations

      REAL(fpk), INTENT(OUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(OUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(OUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector vector with 3-point corr
      REAL(fpk), INTENT(OUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector vector for BOA geo.          4 = 1 + 2.
      REAL(fpk), INTENT(OUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector vector with Multipoint corr. 5 = 1 + 3.

!  Twilight flag and exception handling

      LOGICAL      , intent(out) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(out) :: Local_Message, Local_Action

!  Local
!  -----

      INTEGER, parameter :: maxmults = MAXLAYERS + 1

!  Variables for doing the sphericity calculation

      REAL(fpk) :: TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      REAL(fpk) :: EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_ALL(maxmults,3),  OBSGEOMS_BOA(3), COSSCAT
      REAL(fpk) :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NMULT, NG, NG1, NS, NLAYERS, DIRIDX

!  Start Code
!  ==========

!  Initialize

      FO_STOKES_BOAGEOM   = zero
      MS_STOKES_BOAGEOM   = zero
      MS_STOKES_SPHERCORR = zero
      STOKES_BOAGEOM      = zero
      STOKES_SPHERCORR    = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Local_Message = ' '
      Local_Action  = ' '

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 3-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set other geometries through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Geometrical Conversion. 3/1/20. Multipoint calculation

      NMULT = NLAYERS + 1
      CALL BOATOA_Mpoint_conversion ( &
          maxmults, maxlayers, nlayers, nmult, DirIdx,                     & ! Input numbers
          EARTH_RADIUS, VLIDORT_FixIn%Chapman%TS_height_grid, OBSGEOMS_BOA, & ! input Heights, BOA Geometry
          cosscat, OBSGEOMS_ALL )                                            ! output all geometries, Cosine scattering angle

!  Multi-points: Reset the VLIDORT geometrical input at all layer boundaries from BOA to TOA
!    NMULT is TOA geometry

      DO N = 2, NMULT
        VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(N,1:3) = OBSGEOMS_ALL(N,1:3)
        VLIDORT_ModIn%MSunrays%TS_SZANGLES(N)                = OBSGEOMS_ALL(N,1)
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(N)     = OBSGEOMS_ALL(N,2)
        VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(N)            = OBSGEOMS_ALL(N,3)
      ENDDO

!  re-set geometry numbers. NMULT Geometries in all

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = NMULT
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = NMULT
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = NMULT
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = NMULT

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_ALL(NMULT,1).gt.90.0 ) then
         Twilight_flag = .true. ; return
      endif

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MSMpt_SPHERICITY_CORRECTION_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store STOKES Vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS)    = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)       = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS)    = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)       = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  DEBUG: Check geometry
!      do n = 0, nlayers
!         NG = NMULT - n
!  write(*,*)VLIDORT_Out%Main%TS_PATHGEOMS(1,N)/DEG_TO_RAD,VLIDORT_Out%Main%TS_PATHGEOMS(2,N)/DEG_TO_RAD,OBSGEOMS_ALL(NG,1:2)
!      enddo
!stop

!  STOKES Vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,NMULT,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,NMULT,1:NS) 
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,NMULT,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,NMULT,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,NMULT,1:NS,DNIDX)
      endif

!  perform multi-point SPHERICITY CORRECTION, Recurrence Relation
!  ==============================================================

!  set the surface MS source term, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  Alternative: get the Layer MS source terms by layer averaging the multi-point values
!     -- Always works, however may not be accurate around azimuth-flip conditions

      do N = 1, nlayers
         NG = NMULT - N + 1 ; NG1 = NG - 1
         AMSST(N,1:NS) = HALF * ( VLIDORT_Out%Main%TS_LAYER_MSSTS(NG,1:NS,N) + VLIDORT_Out%Main%TS_LAYER_MSSTS(NG1,1:NS,N))
!            write(398,'(i2,2f10.6,2x,1p2e16.6,2x,1pe16.6)')N,&
!               COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N)),COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)),&
!               VLIDORT_Out%Main%TS_LAYER_MSSTS(NMULT,1:NS,N), VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N), AMSST(N,1:NS)
      ENDDO

!  Verbose debug output

      if ( verbose ) then
        do NG = 1, NMULT
          write(399,'(i2,2f8.5,1x,1p23e11.4)')NG,&
               COS(OBSGEOMS_ALL(NG,1)*DEG_TO_RAD),COS(OBSGEOMS_ALL(NG,2)*DEG_TO_RAD),&
               VLIDORT_Out%Main%TS_LAYER_MSSTS(NG,1:NS,1:NLAYERS)
        ENDDO
        stop
      endif

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
!mick fix 1/5/2021 - trim dim passed from SMSST
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final STOKES Vector answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MSMpt_SPHERICITY_CORRECTION_V1

!  End module

END MODULE VLIDORT_SPHERCORR_ROUTINES_m

