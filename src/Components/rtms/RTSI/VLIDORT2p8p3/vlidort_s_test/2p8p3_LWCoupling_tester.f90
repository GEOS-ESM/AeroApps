
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
                                                   
program VSLEAVE_LWCoupling_Tester

!  This 2.8.1 VLIDORT Driver tests the Water-leaving adjustment facility
!    PROGRAMMED BY ROBERT SPURR, 4/10/19, 5/9/19
  
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

   integer, parameter :: maxlambdas = 900
   integer, parameter :: maxlevels  = maxlayers + 1

   integer :: nlambdas, data_nheights, nlayers, nmoments_input
   integer :: i, l, n, v, vd, w, uta, nmi
   character(len=99) :: trace

   logical :: openfileflag
   logical :: do_observation_geometry

   real(fpk) :: aerwt, raywt, raysca, aersca, aertau, tottau, totsca
   real(fpk) :: q24, g24, qgr, ggr, ang, aerang, aerssa, aerasy
   real(fpk) :: cumaer1, cumaer2, aerprv1, aerprv2, aers_550, aert_550
   real(fpk) :: rho_1, rho_2, rho_a, diff, temp, psurf, ray2
   real(fpk) :: aerphas(maxmoments_input),aer550(maxlayers)

   real(fpk) :: lambdas(maxlambdas), lambda_start
   real(fpk) :: rayleigh_xsec(maxlambdas),rayleigh_depol(maxlambdas)

   real(fpk) :: height(maxlevels),tempr(maxlevels)
   real(fpk) :: data_pressures(maxlevels),data_heights(maxlevels)
   real(fpk) :: data_temperatures(maxlevels)

   real(fpk) :: height_grid(0:maxlayers),layer_temperatures(maxlayers)
   real(fpk) :: layer_aircolumns(maxlayers)

   real(fpk) :: deltau_vert_input(maxlayers)
   real(fpk) :: omega_total_input(maxlayers)
   real(fpk) :: phasmoms_total_input(0:maxmoments_input,maxlayers)

   real(fpk) :: lambertian_albedo

!  Saved results

   real(fpk) :: I_base_TOA(maxlambdas,max_geometries)
   real(fpk) :: I_base_BOA(maxlambdas,max_geometries)

!  Loschmidt's number (particles/cm3).

   real(fpk), parameter :: RHO_STANDARD = 2.68675D+19

!  STP values.

   real(fpk), parameter :: PZERO  = 1013.25D0, TZERO  = 273.15D0, PTZERO = PZERO / TZERO
   real(fpk), parameter :: RHO_ZERO = RHO_STANDARD * TZERO / PZERO, CONSTANT = 1.0D+05 * RHO_ZERO

!  CO2 PPMV mixing ratio (for the Rayleigh stuff)

   real(fpk), parameter :: CO2_PPMV_MIXRATIO = 390.0d0

!  Problem setup
!  -------------

!  Wavelengths

   lambda_start = 440.0d0
   nlambdas     = 1
   do w = 1, nlambdas
      lambdas(w) = lambda_start + 1.0d0*dble(w-1)
   enddo

!  Get the pressures, temperatures and heights from FILE.

   OPEN(1,FILE='vlidort_s_test/data/temp_psurf.prf_35',STATUS='OLD')
   READ(1,*)DATA_NHEIGHTS,psurf
   if ( Data_nheights.gt.maxlayers+1) Stop 'need more dimensioning'
   DO I = 1, DATA_NHEIGHTS
       READ(1,*)HEIGHT(I), TEMPR(I)
   ENDDO
   CLOSE(1)

!  Set the pressures

   DATA_PRESSURES(DATA_NHEIGHTS)=psurf
   DO I = DATA_NHEIGHTS-1, 1 ,-1
      DATA_PRESSURES(I) = DATA_PRESSURES(&
       I+1)*exp(-9.81*28.9/8314.0*(height(DATA_NHEIGHTS-I+1)-&
       height(DATA_NHEIGHTS-I))*1.0E3/2.0*(1.0/tempr(DATA_NHEIGHTS-I+1)&
       +1.0/tempr(DATA_NHEIGHTS-I)))
   ENDDO

!  Set the data heights

   DO I = 1, DATA_NHEIGHTS
      DATA_HEIGHTS(I)=HEIGHT(DATA_NHEIGHTS-I+1)
      DATA_TEMPERATURES(I)=TEMPR(DATA_NHEIGHTS-I+1)
   ENDDO

!  Set the height grid, number of layers, Layer temperatures and Aircolumns

   HEIGHT_GRID = 0.0d0
   NLAYERS = DATA_NHEIGHTS - 1
   HEIGHT_GRID(0) = DATA_HEIGHTS(1)
   DO I = 1, NLAYERS
      RHO_1 = DATA_PRESSURES(I) / DATA_TEMPERATURES(I)
      RHO_2 = DATA_PRESSURES(I+1) / DATA_TEMPERATURES(I+1)
      TEMP  = 0.5D0*(DATA_TEMPERATURES(I)+DATA_TEMPERATURES(I+1))
      RHO_A = 0.5D0 * CONSTANT * ( RHO_1 + RHO_2 )
      DIFF  = DATA_HEIGHTS(I) - DATA_HEIGHTS(I+1)
      HEIGHT_GRID(I) = DATA_HEIGHTS(I+1)
      LAYER_TEMPERATURES(I) = TEMP
      LAYER_AIRCOLUMNS(I)   = DIFF * RHO_A
   ENDDO

!  Rayleigh function

   call Rayleigh_function &
    ( MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
      NLAMBDAS,   LAMBDAS, &
      RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Loading areas

   aers_550 = 0.01d0
   q24 = 2.0d0 * aers_550 / (height_grid(0)-height_grid(24))
   g24 = q24 /  (height_grid(0)-height_grid(24))

   aert_550 = 0.15d0
   qgr = 2.0d0 * aert_550 / (height_grid(24)-height_grid(nlayers))
   ggr = qgr / (height_grid(24)-height_grid(nlayers))

!  Basic layering for aerosol

   aerprv1 = 0.0d0
   aerprv2 = 0.0d0
   do n = 1, nlayers
      if ( n.le.24 ) then
         cumaer1 = g24 * (height_grid(0)-height_grid(n))
         aer550(n) = 0.5*(cumaer1+aerprv1)*(height_grid(n-1)-height_grid(n))
         aerprv1 = cumaer1
      else
         cumaer2 = ggr * (height_grid(24)-height_grid(n))
         aer550(n) = 0.5*(cumaer2+aerprv2)*(height_grid(n-1)-height_grid(n))
         aerprv2 = cumaer2
      endif
   enddo

!  Other aerosol inputs

   aerAng = 1.01d0
   aerssa = 0.99d0
   aerasy = 0.75d0
   nmoments_input = 50
   do l = 1, nmoments_input
      aerphas(L) = dble(2*L+1) * aerasy ** dble(l)
   enddo

   !  Define a proxy

   nmi = nmoments_input

!  Initialize error file

   OPENFILEFLAG = .false.

!  Start wavelength loop
!  ---------------------

   do w = 1, nlambdas
      write(*,*) 'w = ',w

!  Read SLEAVE inputs

      CALL VSLEAVE_INPUTMASTER ( &
        'vlidort_s_test/2p8p3_VSLEAVE_LWCoupling.cfg', & ! Input
         VSLEAVE_Sup_In,         & ! Outputs
         VSLEAVE_Sup_InputStatus ) ! Outputs

      IF ( VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VSLEAVE_READ_ERROR ( '2p8p3_VSLEAVE_LWCoupling.log', VSLEAVE_Sup_InputStatus )

!  Define some inputs not handled by VSLEAVE_LIN_INPUTMASTER

      VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = 'vlidort_s_test/data'

!  Redefine SLEAVE input wavelength (override SLEAVE config file)

      VSLEAVE_Sup_In%SL_FL_WAVELENGTH = lambdas(w) !in nm

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
        'vlidort_s_test/2p8p3_VLIDORT_LWCoupling.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( '2p8p3_VLIDORT_LWCoupling.log', VLIDORT_InputStatus )

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

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = lambdas(w)/1000.0_fpk !in um

!  IMPORTANT - Check compatibility of VSLEAVE and VLIDORT Main inputs

      CALL VLIDORT_VSLEAVE_INPUT_CHECK ( &
        VSLEAVE_Sup_In,             & ! Inputs
        VLIDORT_FixIn,              & ! Inputs
        VLIDORT_ModIn,              & ! Inputs
        VLIDORT_SLEAVECheck_Status )  ! Outputs

!  Exception handling

      IF ( VLIDORT_SLEAVECheck_Status%TS_STATUS_INPUTCHECK .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_VSLEAVE_INPUT_CHECK_ERROR ( &
           '2p8p3_VLIDORT_VSLEAVE_checkcall.log', VLIDORT_SLEAVECheck_Status )

!  Copy VSLEAVE Sup outputs to VLIDORT's VSLEAVE Sup inputs (std only)

      CALL SET_VLIDORT_VSLEAVE_INPUTS ( &
        VSLEAVE_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn,  & !Inputs
        VLIDORT_Sup )                                     !Outputs

!  External setting

      IF ( VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE ) then
         write(*,*)'Reading external field'
         open(1,file='External_Wleave.dat',status='old')
         do v = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
            read(1,*)vd, VLIDORT_Sup%Sleave%TS_SLTERM_ISOTROPIC(1,v)
         enddo
         close(1)
      endif

!  Angstrom and ray2mom

      ray2 = ( 1.0d0 - rayleigh_depol(w) ) / ( 2.0d0 + rayleigh_depol(w) )
      ang  = ( lambdas(w) / 550.0d0 ) ** (-aerAng)

!      write(*,*)w,lambdas(w),ang,rayleigh_depol(w)

!  Make VLIDORT optical properties

      do n = 1, nlayers
         aertau = aer550(n) * ang
         aersca = aertau * aerssa
         raysca = LAYER_AIRCOLUMNS(n) * rayleigh_xsec(w)
         tottau = raysca + aertau
         totsca = raysca + aersca
         deltau_vert_input(n) = tottau
         omega_total_input(n) = totsca / tottau
         raywt = raysca / totsca ; aerwt = 1.0d0 - raywt
         phasmoms_total_input(0,n) = 1.0d0
         phasmoms_total_input(1,n) = aerwt * aerphas(1)
         phasmoms_total_input(2,n) = aerwt * aerphas(2) + raywt * ray2
         do L = 3, nmoments_input
            phasmoms_total_input(L,n) = aerwt * aerphas(L)
         enddo
      enddo

      lambertian_albedo = 0.02d0

!  Debug

!      if ( w.eq.1) then
!        do n = 1, nlayers
!          write(*,*)w,n,deltau_vert_input(n),omega_total_input(n),&
!                    (phasmoms_total_input(L,n), L = 0,20)
!        enddo
!      endif

!  Set VLIDORT Type structures

      VLIDORT_FixIn%Cont%TS_nstokes                   = 1      
      VLIDORT_FixIn%Cont%TS_nlayers                   = nlayers      
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input     = nmoments_input
      VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = height_grid(0:nlayers)

      VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)          = deltau_vert_input(1:nlayers)
      VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)         = omega_total_input(1:nlayers)
      VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1) = phasmoms_total_input(0:nmi,1:nlayers)
      VLIDORT_FixIn%Optical%TS_lambertian_albedo                     = lambertian_albedo

!   VLIDORT call

      do_debug_input = .false.
      CALL VLIDORT_MASTER ( do_debug_input,  &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup,   &
        VLIDORT_Out )

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS (  &
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Save baseline results

      do v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
        I_base_TOA(w,v) = VLIDORT_Out%Main%TS_Stokes(1,v,1,upidx)
        I_base_BOA(w,v) = VLIDORT_Out%Main%TS_Stokes(2,v,1,dnidx)
      enddo

!  write Sleave results

      if ( w.eq.1 ) then
         open(1,file='vlidort_s_test/Sleave_Isotropic_UNadjusted.dat',status='unknown')
         do v = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
            write(1,300) v, VLIDORT_Sup%Sleave%TS_SLTERM_ISOTROPIC(1,v)
         enddo
         close(1)
      endif
      
!  End wavelength loop
      
   enddo
   
!  write results
   
   open(1,file = 'vlidort_s_test/results_LWCoupling_TOAUp_BOADn_Radiances.all',&
        status = 'unknown')
   do v = 1,  VLIDORT_Out%Main%TS_N_GEOMETRIES
      write(1,310) v, I_base_TOA(1:nlambdas,v), I_base_BOA(1:nlambdas,v)
   enddo
   close(1)

   IF ( VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT ) then
      open(1,file='vlidort_s_test/Sleave_Isotropic_WLAdjusted.dat',status='unknown')
      do v = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
         write(1,300) v, VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(1,v)
      enddo
      close(1)
   endif

300   format(1x,i2.2,1x,1pe16.7)
310   format(1x,i2.2,2(1x,1pe16.7))

!  Finish

   write(*,*)
   write(*,*) 'Main program finished successfully'
   stop

!mick fix 3/22/2017 - added error finish section
!  Error finish

678   continue

   write(*,*)
   write(*,'(1x,a)') trim(TRACE)
   write(*,*)'Number of error messages from SLEAVE calculation = ',&
              VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
   do i = 1, VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
      write(*,'(a,i2,a,a)')' Message # ', i, ': ',&
         adjustl(trim(VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)))
   enddo

   stop
end program VSLEAVE_LWCoupling_Tester

!

SUBROUTINE RAYLEIGH_FUNCTION &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS, &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.

      USE VLIDORT_PARS_m, only: fpk

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  wavelength

      INTEGER          FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      real(fpk) FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )

!  CO2 mixing ratio

      real(fpk) CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      real(fpk) RAYLEIGH_XSEC  ( FORWARD_MAXLAMBDAS )
      real(fpk) RAYLEIGH_DEPOL ( FORWARD_MAXLAMBDAS )

!  Local variables
!  ---------------

      INTEGER          W
      real(fpk) MASS_DRYAIR
      real(fpk) NMOL, PI, CONS
      real(fpk) MO2,MN2,MARG,MCO2,MAIR
      real(fpk) FO2,FN2,FARG,FCO2,FAIR
      real(fpk) LAMBDA_C,LAMBDA_M,LPM2,LP2
      real(fpk) N300M1,NCO2M1,NCO2
      real(fpk) NCO2SQ, NSQM1,NSQP2,TERM
      real(fpk) S0_A, S0_B
      real(fpk) S1_A, S1_B, S1_C, S1_D, S1_E
      real(fpk) S2_A
      real(fpk) S3_A, S3_B, S3_C, S3_D, S3_E

!  data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      PARAMETER        ( S0_A = 15.0556D0 )
      PARAMETER        ( S0_B = 28.9595D0 )

      PARAMETER        ( S1_A = 8060.51D0 )
      PARAMETER        ( S1_B = 2.48099D+06 )
      PARAMETER        ( S1_C = 132.274D0 )
      PARAMETER        ( S1_D = 1.74557D+04 )
      PARAMETER        ( S1_E = 39.32957D0 )

      PARAMETER        ( S2_A = 0.54D0 )

      PARAMETER        ( S3_A = 1.034D0 )
      PARAMETER        ( S3_B = 3.17D-04 )
      PARAMETER        ( S3_C = 1.096D0 )
      PARAMETER        ( S3_D = 1.385D-03 )
      PARAMETER        ( S3_E = 1.448D-04 )

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop

      DO W = 1, FORWARD_NLAMBDAS

!  wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  end loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION
