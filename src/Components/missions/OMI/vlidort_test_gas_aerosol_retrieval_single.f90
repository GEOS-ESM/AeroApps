#  define I_AM_MAIN
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"

      program vlidort_test_gas_retrieval

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use Chem_MieMod
!  Module files for Rayleigh+aerosol ops related
      USE create_atm_prfl_single
!  Module file for pre-processed data readin
 
      USE gravity, only: gnot, gee
      USE rayleigh
      USE rd_fp_single
!  Module files for VLIDORT

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE, intrinsic :: iso_fortran_env , only : REAL64, REAL32
      USE VLIDORT_AUX_m
      USE VLIDORT_INPUTS_m, only : VLIDORT_INPUT_MASTER
      USE VLIDORT_MASTERS_m

!  Module files for VBRDF: everything in here

      USE VBRDF_SUP_MOD_m

!  Module files for SLeave: everything in here

      USE VSLEAVE_SUP_MOD_m

!  Module for the VLIDORT/VBRDF checker routine
      USE VLIDORT_VBRDF_SUP_ACCESSORIES_m

!  Module for the VLIDORT/VSLEAVE checker routine
      USE VLIDORT_VSLEAVE_SUP_ACCESSORIES_m

     !     USE VLIDORT_SUP_ACCESSORIES

      IMPLICIT NONE
     include "mpif.h"

! ESMF Objects
!----------------
  type(ESMF_Config)       :: cf
  type(ESMF_VM)           :: vm 


!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  for non-Lambertian surface

!  VBRDF supplement file inputs status structure

      TYPE(VBRDF_Input_Exception_Handling)   :: VBRDF_Sup_InputStatus

!  VLIDORT VBRDF supplement input structure

      TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

!  VLIDORT VBRDF supplement output structure

      TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out

!  VLIDORT VBRDF supplement output status structure

      TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus

!  VBRDF supplement / VLIDORT VBRDF-related inputs consistency check status

   TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_VBRDFCheck_Status

! === start water-leaving contribution defs

!  VSLEAVE supplement file inputs status structure

   TYPE(VSLEAVE_Input_Exception_Handling)   :: VSLEAVE_Sup_InputStatus

!  VLIDORT VSLEAVE supplement input structure

   TYPE(VSLEAVE_Sup_Inputs)                 :: VSLEAVE_Sup_In

!  VLIDORT VSLEAVE supplement output structure

   TYPE(VSLEAVE_Sup_Outputs)                :: VSLEAVE_Sup_Out
   TYPE(VSLEAVE_Output_Exception_Handling):: VSLEAVE_Sup_OutputStatus

!  VSLEAVE supplement / VLIDORT VSLEAVE-related inputs consistency check status

   TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_VSLEAVECheck_Status

!  Additional checking for the supplements

   TYPE(VLIDORT_Exception_Handling)       :: VBRDF_VSLEAVECheck_Status

   type (aer_data)       :: aer_data_rec

! === end water-leaving related defs

!  pre-processed omi pxl based input data

!  MERRA-2 aerosol profile for omi pxl 

!  VBRDF supplement variables

      INTEGER       :: BS_NMOMENTS_INPUT
      LOGICAL       :: DO_DEBUG_RESTORATION
      LOGICAL       :: exist

!  water leaving supplement variables

      DOUBLE PRECISION :: CHLORCONC

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  save atm profiles for debugging
      LOGICAL ::          debug
      LOGICAL ::          do_arz_inp

!  use oz xsec at exact wvs (for PACE simulation)
      LOGICAL ::          exactozxsc

!     DOUBLE PRECISION :: btmopd, btmssa
      real, allocatable  :: wbrc(:,:) ! two bracket wvs for oz_abs interpolation

      
      integer(kind=4), parameter :: FORWARD_MAXLAMBDAS = 10000
      real(kind=8) :: e_mobyrad=2.57365E-02, LW_0, LW_1, dy, e_heights(0:23)
      real(kind=8) :: K_Isotropic, K_Direct, K_SOLARVERT, Lw_External, simrad,simrad_0, simjac, trans_atmos
      real(kind=8), parameter :: conv      = 1.0d-08
      integer  , parameter :: maxiter = 1

      real(kind=4), dimension(FORWARD_MAXLAMBDAS) :: FORWARD_LAMBDAS
      real(kind=8), dimension(FORWARD_MAXLAMBDAS) :: NOUTM1
      integer, parameter :: sp = real32
      real(kind=sp) :: g_actual



!  omi geom related variables
      integer ::       io, ndts, wave_ind, lay, j, w, wv_ind
      integer ::       ix,nx,npx ! # of x-track & along-track pixels of omi
!     integer ::       n_geometries ! temp var
      integer ::       n_obsgeoms ! # obs
      real, allocatable :: svgeoms(:),wvs(:) 
      real ::          obsgeoms(3) ! sza, vza and raa
      real ::          w1,w2,ws ! start,end and step of wvs
!     real ::          wvs(nwv), lwf ! wv & criterion for ocean pxl
      real ::          lndf,fillval,kd, tot_tau ! land frac, fill value, help value
      real ::          elv,pelv,sfw,wdir ! pxl avg DEM & P, sfc windspeed & direction
      real ::          wl_rad ! =iso term or user term dep on do_isotropic setup
!     real ::          SL_TRANS ! =Atmospheric Transmittance at ocean surface for a given sza and nstoke (from pre-calculated tables)
      real ::          toaR,boaR,boaF,boaT,tbrf,ubrf ! toa and boa rads, true & unadj sfc brf
      real ::          flux_factor   ! solar flux for sfc brf calculation
!  pxls related info
      real, allocatable, dimension(:) :: lats, lons,pelvs,sfws,wdirs,chlors
      
      integer(kind=4) :: FORWARD_NLAMBDAS
      real(kind=8) :: CO2_PPMV_MIXRATIO, MCO2, MAIR

      INTEGER     ::   error ! Error flag

!  define some constants first
      real, parameter :: deltH=500.0  ! H btw TOA & last level
      real, parameter :: T0=273.13 ! K zero
      real, parameter :: fo2=0.21 ! o2 fraction in dry air

!  omi pixel-based atmos profile related variables

!     integer, allocatable ::       subp(:) ! subset index in 72 layers
      real, allocatable, dimension (:)  :: ai,bi ! GMI coeffs to create P grid
      real, allocatable, dimension (:)  :: P_grid, T_grid ! level P & T

      real(8), allocatable :: rayleigh_xsec(:),rayleigh_depol(:)
      real(8), allocatable :: oz_xsec(:,:,:) ! [2,nwv,3]
      real(8), allocatable :: layer_aircol(:), layer_ozcol(:,:), layer_o4col(:)
      real(4), allocatable :: accoz(:), oz_abs(:,:)
      real(4), allocatable :: o4_xsec(:,:)
      real(8), allocatable :: ozprf(:,:)


      real ::          lat,lon  ! omi pxl location
      real ::          t1,t2  ! level Ts for a layer

!  control file read in variables





      character (len=4), allocatable :: cwv(:)
      integer (kind = 4):: ip,i,iw,iwv,ilay,err,npxl,ip0,ip1, wv0, wv1,  wv0_out,wv1_out, out_wv_ind1,out_wv_ind2, lev

      real ::          start, finish ! for cpu time monitoring
      real ::          t_brdf, t_wl,t_ler,t_tot ! tot cpu time for each module

!  Help variables

      INTEGER ::          N,L,NDUM,LDUM,V,V1,K
      DOUBLE PRECISION :: beta2,depol
      real(8), allocatable :: bt3(:)
      character*1 ::      cpx,coz ! pixel #, oz prf from coz
      CHARACTER*90 ::     TRACE

!  VLIDORT standard input preparation; proxy variables

      LOGICAL ::          DO_RAYLEIGH_ONLY
      LOGICAL ::          do_o3, do_o2o2,do_tocscaling,do_dblaer, do_leaving, trawl  , set_Lw
      INTEGER ::          nlevels,NSTOKES
      INTEGER ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION :: height_grid (0:MAXLAYERS)
      DOUBLE PRECISION, allocatable :: omega_total_input (:, :)
      DOUBLE PRECISION, allocatable :: deltau_vert_input (:, :)
      DOUBLE PRECISION, allocatable :: greekmat_ray_input (:,:,:)
      DOUBLE PRECISION, allocatable :: greekmat_aer_input (:,:,:, :)
      DOUBLE PRECISION, allocatable :: greekmat_total_input (:,:,:, :)
      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO
      real :: tocstn ! ref toc for scaling
      real (kind = 4)  :: p0  = 1013.25   ! sfc P
      real (kind = 4), dimension (72)   :: h_mid,delp, t_mid



      INTEGER :: AOS_MODEL
      LOGICAL :: ABTOT_CONTROL


!  Aerosol
   LOGICAL            :: DO_AEROSOLS
   integer            :: aer_nmom  !  number of moments       
   integer, allocatable :: p2a_ind(:)  !  p_ij ->a_k matchup index
!  LOGICAL,      DIMENSION( MAXLAYERS )            :: AEROSOL_LAYERFLAGS
   real(fpk), allocatable, DIMENSION(:, :) :: deltau_ray, deltau_ozone, deltau_o4

   real(fpk), allocatable, DIMENSION(:, :) :: aerosol_deltau, aerosol_ssalbs
   real(fpk), allocatable, DIMENSION(:, :, :) :: aerosol_deltau_all, aerosol_ssalbs_all
   real(fpk), allocatable :: AEROSOL_SCATMOMS(:,:,:,:)

   real      ::  lgwvs(4) ! available wvs for legendre coefs
   character (len=3) :: lgcwv(4)
    LOGICAL :: DO_DEBUG_INPUT

    DOUBLE PRECISION :: I_base_TOA(2)
    DOUBLE PRECISION :: EPS, EPSFAC, LW_BASIC_SAVED

    character*255 :: fp_filename, out_filename

!  detector dep varibles
    NAMELIST /VSLEAVE/ VSLEAVE_Sup_In
    NAMELIST /VBRDF/ VBRDF_Sup_In
    NAMELIST /VLIDORT/ VLIDORT_FixIn, VLIDORT_ModIn
    !NAMELIST /MERRA2/ aer_data_rec

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status                            ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)                 ! MPI status
  integer                               :: myid, npet, CoresPerNode                    ! MPI dimensions and processor id
  logical                               :: amOnFirstNode

! System tracking variables
! -----------------------------
  character(len=*), parameter           :: Iam = 'omi_vlidort'

! Initialize MPI with ESMF
! ------------------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet)
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'


    call GET_COMMAND_ARGUMENT(1,fp_filename)
    call GET_COMMAND_ARGUMENT(2,out_filename)


    open(14,file=out_filename)
    
    !write(14,*) 'Year, Month, Day, Orbit, Line, Row, AOD, VZA, Wavelength, Lw_Basic_0,Lw_Basic, Lw_Adjusted, &
    !     Lw_Retrieved, Kdir, Kiso, Trans_Atmos, Jacobian, Radiance_0, Radiance, ES'
    
    write(14,*) 'Wave, Lw_Bas_0,Lw_Bas, Lw_Adj, Kdir, Ksolar, Kiso, TAtm, Jac, Rad_0, Rad, ES'

    AOS_MODEL = 99

   !  3/3/22. ABTOT_CONTROL MUST BE FALSE to compare with Version 2.8.3 SLEAVE !!!!!!!!!!!!

    !ABTOT_CONTROL = .true.
    ABTOT_CONTROL = .false.
   
!  epsfac

   eps = 0.01_fpk
   epsfac = one + eps

!  Initialize error file
   OPENFILEFLAG = .false.

!  save the atm profile inputs
   debug = .true.

   do_leaving = .true.
   do_o3 = .true.
   do_o2o2 = .true.
   do_aerosols = .true. 


   fillval=-999.0




   
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Step 1. VBRDF/VSLEAVE inputs + check consistency
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Step 1a. VBRDF control Read input, abort if failed
!  =================================================


      CALL VBRDF_INPUTMASTER ( &
        'vlidort_test/VBRDF_ReadInput.cfg' , & ! Input
        VBRDF_Sup_In,         & ! Outputs
        VBRDF_Sup_InputStatus ) ! Outputs
      


!  Exception handling. Use Type structure directly

!  Step 1b. VSLEAVE control Read input, abort if failed
!  =================================================

      CALL VSLEAVE_INPUTMASTER ( &
        'vlidort_test/VSLEAVE_ReadInput.cfg' , & ! Input
        VSLEAVE_Sup_In,         & ! Outputs
        VSLEAVE_Sup_InputStatus ) ! Outputs


!  set Vsleave_DataPath. Input setting, 12/28/15. New variable from Rob






!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Step 2. VLIDORT control Read input, abort if failed
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_test/VLIDORT_ReadInput.cfg', & ! Input for wleaving consideration
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs




      VBRDF_Sup_In%BS_DO_NewCMGLINT  = .false.
      VBRDF_Sup_In%BS_DO_NewGCMGLINT = .true.
      VBRDF_Sup_In%BS_WHICH_BRDF(1) = 16


      !VBRDF_Sup_In%BS_DO_NewCMGLINT  = .true.
      !VBRDF_Sup_In%BS_DO_NewGCMGLINT = .false.
      !VBRDF_Sup_In%BS_WHICH_BRDF(1) = 15

      !flux_factor=VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR
      !flux_factor = 3.1415926535897403 * flux_factor
      !VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor


      VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-12
      VLIDORT_FixIn%Bool%TS_DO_MSSTS                       = .False.
      VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2 = .True.

      VLIDORT_FIXIN%BOOL%TS_DO_WLADJUSTED_OUTPUT=.true.
      VLIDORT_FIXIN%BOOL%TS_DO_TOA_CONTRIBS = .true.
      
      VLIDORT_FIXIN%BOOL%TS_DO_DNWELLING=.true.

      !  Want these flags fixed -- external off, Coupling On for both parts.
      



! =======need to reset wv in vbrdf and vsleaving readins


!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Step 3. Check VBRDF/VSLEAVE inputs against VLIDORT inputs
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Step 3a. Check VBRDF inputs against VLIDORT inputs, abort if failed
!  ==================================================================

!  Routine in vlidort_sup2p7_accessories.f90

!  if not doing WL, manually input wv to satify VLIDORT_VSLEAVE_INPUT_CHECKER

!     CALL VLIDORT_VBRDF_INPUT_CHECKER ( &



      CALL VLIDORT_VBRDF_INPUT_CHECK ( &
        VBRDF_Sup_In,               & ! Inputs
        VLIDORT_FixIn,              & ! Inputs
        VLIDORT_ModIn,              & ! Inputs
        VLIDORT_VBRDFCheck_Status )   ! Outputs



!  Step 3b. Check VSLEAVE against VLIDORT inputs, abort if failed
!  ==============================================================

!  Routine in vlidort_sup_accessories.f90

!  if not doing WL, manually input wv to satify VLIDORT_VSLEAVE_INPUT_CHECKER




!     CALL VLIDORT_VSLEAVE_INPUT_CHECKER ( &
      CALL VLIDORT_VSLEAVE_INPUT_CHECK ( &  ! v2.8.1b
        VSLEAVE_Sup_In,             & ! Inputs
        VLIDORT_FixIn,              & ! Inputs
        VLIDORT_ModIn,              & ! Inputs
        VLIDORT_VSLEAVECheck_Status ) ! Outputs



!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Step 5. Prepare Optical inputs (5a), then call VLIDORT (5b)
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Step 5a. Prepare optical property inputs for the atmosphere
!  move 5a out of the user orbit/ler loop because atm not changing w/ orbit 

!  Define some input variables not handled by VLIDORT_INPUT_MASTER

      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS    = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF  = 0

      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID    = ZERO
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID = ZERO
      VLIDORT_FixIn%Chapman%TS_FINEGRID         = ZERO

      VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = ZERO ! new
      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT  = ZERO ! new

!  Define default input variables for vbrdf & vsleave
      lambertian_albedo = VLIDORT_FixIn%Optical%TS_lambertian_albedo
      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F_0         = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F           = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0    = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F      = ZERO
      VLIDORT_Sup%BRDF%TS_EMISSIVITY       = ZERO
      VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY  = ZERO

      !VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      !VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      !VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
      !VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!  Get the pre-prepared atmosphere: heights (km) at layer boundaries (from TOA),Rayleigh x-sec etc


      nlevels = VLIDORT_FixIn%Cont%TS_nlayers

      CALL rd_aerpfl(fp_filename,fillval,aer_data_rec)
      
      
      
      if (.not. allocated(bt3)) then
         allocate(bt3(aer_data_rec%nwv))
         allocate(rayleigh_xsec(aer_data_rec%nwv),rayleigh_depol(aer_data_rec%nwv))
         allocate(greekmat_ray_input(0:2,aer_data_rec%nwv,MAXSTOKES_SQ))
         allocate(oz_xsec(2,aer_data_rec%nwv,3),wbrc(2,aer_data_rec%nwv))

      endif



      if (.not. allocated(omega_total_input)) then 
         allocate(oz_abs(aer_data_rec%nwv,nlevels))
         allocate(o4_xsec(aer_data_rec%nwv,nlevels))
         allocate(omega_total_input(aer_data_rec%nwv,MAXLAYERS), deltau_vert_input(aer_data_rec%nwv,MAXLAYERS))
         allocate(deltau_ray(aer_data_rec%nwv,MAXLAYERS))
         allocate(deltau_ozone(aer_data_rec%nwv,MAXLAYERS),deltau_o4(aer_data_rec%nwv,MAXLAYERS)) 
         allocate(greekmat_aer_input(0:MAXMOMENTS_INPUT,aer_data_rec%nwv,nlevels,MAXSTOKES_SQ)) 
         allocate(greekmat_total_input(0:MAXMOMENTS_INPUT,aer_data_rec%nwv,MAXLAYERS,MAXSTOKES_SQ)) 
         allocate(P_grid(nlevels),T_grid(nlevels)) ! from toa -> boa
         allocate(layer_aircol(nlevels),layer_o4col(nlevels))
         allocate(aerosol_deltau(nlevels,aer_data_rec%nwv), aerosol_ssalbs(nlevels,aer_data_rec%nwv))
         
         
      endif
      ! some will be initialize in prep_profile_iops again
      deltau_vert_input = ZERO
      omega_total_input = ZERO
      greekmat_ray_input = ZERO
      greekmat_aer_input = ZERO
      greekmat_total_input = ZERO
      
      

      
      write(*,'(A,i4,A,i2,A,i2,A,i5,A,i4,A,i2)')  'Processing pixel for Year: ',aer_data_rec%year,' Month: ', &
           aer_data_rec%month,' Day: ',aer_data_rec%day,' Orbit: ', aer_data_rec%orbit,' Line: ', aer_data_rec%line,& 
           ' Row: ',aer_data_rec%row

      !Converting ozone profiles from mmr to column DU
      aer_data_rec%oz_layer =  (1.0e6)*aer_data_rec%oz_layer*(28.95949/47.998)*0.0078914*aer_data_rec%delp(:)*0.001
      
         
         
      if (do_aerosols) then 
         VLIDORT_ModIn%MBool%TS_do_deltam_scaling     =.true.
         VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY   = .false.
         
         do iw = 1,aer_data_rec%nwv
            call rd_phsfunc (trim(fp_filename),fillval,aer_data_rec,aer_data_rec%wvs(iw))
            if (iw .eq. 1) then 
               if (.not. allocated(p2a_ind)) then
                  allocate (p2a_ind(aer_data_rec%nPol))
               endif
               !  setup p2a_index matchup index (P_ij -> a_k in Mie F_matrix in vlidort)
               p2a_ind = (/ 1, 2, 5, 3, 4, 6 /)
               
               aer_nmom=aer_data_rec%nmom
               if (.not. allocated(AEROSOL_SCATMOMS)) then
                  allocate(AEROSOL_SCATMOMS(aer_data_rec%nPol,0:aer_nmom-1,aer_data_rec%nwv,nlevels))
               endif
               ! initialize
               AEROSOL_SCATMOMS = ZERO
            endif ! iw=1
            ! need to apply matchup index btw a_k needed by vlidort and p_ij provided in
               ! phase function file
            AEROSOL_SCATMOMS(:,0:aer_nmom-1,iw,1:nlevels)= &
                 aer_data_rec%pmom(p2a_ind,1:aer_nmom,1:nlevels)
            !
            aerosol_deltau(1:nlevels,iw)=aer_data_rec%dtau(iw,:)
            
            
            aerosol_ssalbs(1:nlevels,iw)=aer_data_rec%ssa(iw,:)
         enddo
         where(aerosol_ssalbs .gt. 0.999999) aerosol_ssalbs=0.999999
         where(AEROSOL_SCATMOMS .lt. 1.0d-6) AEROSOL_SCATMOMS=1.0d-6
         
      else
         !deltam scaling used for aerosol only
         VLIDORT_ModIn%MBool%TS_do_deltam_scaling     =.false.
         VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY   = .true.
         
      endif
      

      P_grid(1)=0.01*(1.0+aer_data_rec%delp(1)) ! 1st level, Pa->hPa, TOA=1 Pa
      do n=2,nlevels
         P_grid(n)=aer_data_rec%delp(n-1)+0.01*aer_data_rec%delp(n) ! top -> btm, Pa->hPa
      enddo
      
      
      nstokes = VLIDORT_FixIn%Cont%TS_nstokes
      flux_factor=VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR
      
     
      
      
      call prep_profile_iops (aer_data_rec%wvs,aer_data_rec%nwv,nstokes,nlevels, & ! input
           do_aerosols,aer_nmom-1,AEROSOL_SCATMOMS, & ! input
           greekmat_ray_input,greekmat_aer_input,ngreek_moments_input,&
           rayleigh_xsec,rayleigh_depol) 


      !  pre-process Rayleigh atm for tau and omega
      
      
      n_obsgeoms=VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      
      if ( .not. VBRDF_Sup_In%BS_DO_BRDF_SURFACE ) then
         write(*,*) 'do Lambertian sfc'
         write(*,*) 'sfc alb=', lambertian_albedo
      else ! do brdf sfc
         lambertian_albedo = 0.00d0 ! initial for non-Lambertian surface
      endif
      
      !  Do not want debug restoration. Normal number of moments
      
      DO_DEBUG_RESTORATION = .FALSE.
      BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1
      
      !  initialise cpu time counts
      t_brdf = 0.0d0; t_wl = 0.0d0 ; t_ler = 0.0d0
      


      VSLEAVE_Sup_In%SL_CHLORCONC = aer_data_rec%chl
      VSLEAVE_Sup_In%SL_WINDSPEED = aer_data_rec%wsp
      

      elv=0.001d0*aer_data_rec%h_mid(nlevels) ! m -> km
      
      
      !  Initialise
      height_grid = 0.0d0

      ! calculate level h from h_mid
      do n=1,nlevels-1
         height_grid(n)=0.5*(aer_data_rec%h_mid(n)+aer_data_rec%h_mid(n+1))
      enddo


      if (aer_data_rec%h_mid(nlevels) .ge. 1000.0d0*elv) then
         height_grid(nlevels)=1000.0d0*elv ! use Terrain H as boa elv; km -> m
      else ! h_mid is below elv
         height_grid(nlevels)=0.5d0*aer_data_rec%h_mid(nlevels)
         elv=0.001d0*height_grid(nlevels) ! reset boa elv; m -> km
      endif
      height_grid(0)=height_grid(1)+deltH ! TOA H=deltH + first level H      
      
        
        

      ! for Rayleight atm
      height_grid=0.001d0*height_grid ! -> km
      
      
      CO2_PPMV_MIXRATIO = 360.0d0
      FORWARD_NLAMBDAS = aer_data_rec%nwv
      
      DO i=1, FORWARD_NLAMBDAS
         FORWARD_LAMBDAS(i) = aer_data_rec%wvs(i)
      ENDDO
      
      
      CALL RAYLEIGH_FUNCTION( &
           FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
           FORWARD_NLAMBDAS,   FORWARD_LAMBDAS, &
             rayleigh_xsec, rayleigh_depol, NOUTM1 )
      
      MCO2 = 1.0D-6 * CO2_PPMV_MIXRATIO
      MAIR = 15.0556D0 * MCO2 + MAIR_NOCO2
      
      lat = aer_data_rec%lat
      
      do i = 1,aer_data_rec%nwv
         do lay = 1, 72
            !g_actual = GEE(GNOT(lat), aer_data_rec%h_mid(lay), lat)
            g_actual = GEE(GNOT(lat),0.0, lat)
            deltau_ray(i,lay) = rayleigh_xsec(i)*(aer_data_rec%delp(lay))*10.*NA/MAIR/g_actual
         enddo
      enddo
      
       !deltau_ray = 0.14350495127854637/72.
      
      layer_o4col=0.0d0
      do n = nlevels,1,-1 ! boa -> toa
         t1=aer_data_rec%T(n) ! at lower bound
         !  calc layer_aircol_density & layer_o4col
         
         if (n .gt. 1) then ! for a layer
            
            t2=0.5*(aer_data_rec%T(n)+aer_data_rec%T(n-1)) ! layer T -> level T
            kd=height_grid(n-1)-height_grid(n)
            beta2=0.5*(P_grid(n)/t1+P_grid(n-1)/t2)
            if (do_o2o2) layer_o4col(n)=fo2*fo2*beta2*beta2*kd
         endif ! n>1
      enddo ! n=nlayers -> 1


      call ozone_xsec(aer_data_rec%wvs,aer_data_rec%nwv,oz_xsec,wbrc)
      if (do_o2o2) call o2o2_xsec_fspecm (aer_data_rec%wvs,aer_data_rec%nwv,aer_data_rec%delp,nlevels,o4_xsec)
      !aer_data_rec%o4_xsec = o4_xsec

      oz_abs=0.0d0
      do n = nlevels,1,-1 ! boa -> toa
          t1=aer_data_rec%T(n) ! at lower bound
          if (n .eq. 1) then ! the 1st layer (toa)
             kd=aer_data_rec%T(1)-T0 ! F -> C
          else
             kd=aer_data_rec%T(n)-T0 ! assume T_grid is layer mean T
          endif
          bt3=oz_xsec(1,:,1)+oz_xsec(1,:,2)*kd+ &
                  oz_xsec(1,:,3)*kd*kd ! oz_abs at wbrc(1)
          oz_abs(:,n)=bt3 ! keep oz_abs at wbrc(1)
       enddo
       
       !do lay = 1, 72
       !   print *, 'Layer: ',lay,' Layer Ozone: ',aer_data_rec%oz_layer(lay),' Ozone Abs: ', oz_abs(1,lay)
       !enddo
       !stop
       call generate_vlidort_iops (deltau_ray,aer_data_rec%wvs,aer_data_rec%nwv,nlevels,1, & ! input
           aer_data_rec%oz_layer,oz_abs,o4_xsec,layer_o4col, & ! input 
           greekmat_ray_input, greekmat_aer_input, & ! input
           do_aerosols, aerosol_deltau, aerosol_ssalbs, & !  Inputs
            deltau_vert_input,omega_total_input,greekmat_total_input) ! output
      
      


      do i = 1,aer_data_rec%nwv
         do lay = 1, 72
            deltau_ozone(i,lay) = oz_abs(i,lay)*aer_data_rec%oz_layer(lay)
            deltau_o4(i,lay) = o4_xsec(i,lay)*layer_o4col(lay)
         enddo
      enddo
      
      
      
      !  Copy to type-structure inputs
      
      
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid
      
      
        

        
      !  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     Step 4. Call VBRDF/VSLEAVE supplements, copy results to VLIDORT inputs
      !  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      !  Step 4a. Call the VBRDF supplement
      !  ----------------------------------
      call cpu_time(start) ! start time
      
      do iw=1,aer_data_rec%nwv ! for testing detector dep geoms

         print *, 'TOTAL TAU is ', sum(deltau_vert_input(iw,:))
         PRINT *, 'Ozone Tau: ',SUM(deltau_ozone(iw,:))
         PRINT *, 'O4 Tau: ',SUM(deltau_o4(iw,:))
         PRINT *, 'Rayleigh Tau: ',SUM(deltau_ray(iw,:))
         print *, 'Aerosol Tau: ',sum(aerosol_deltau(:,iw))
         
           
           !      write(*,*) 'iw=',iw
           
           
           
           obsgeoms(1)=aer_data_rec%sza
           obsgeoms(2)=aer_data_rec%vza
           obsgeoms(3)=aer_data_rec%raa

           VSLEAVE_Sup_In%SL_USER_OBSGEOMS(1,:)=obsgeoms(:)
           VSLEAVE_Sup_In%SL_BEAM_SZAS(1)=obsgeoms(1)
           VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT(1)=obsgeoms(2)
           VSLEAVE_Sup_In%SL_USER_RELAZMS(1)=obsgeoms(3)
           VSLEAVE_Sup_In%SL_WAVELENGTH=aer_data_rec%wvs(iw)*0.001d0 ! nm -> micron
           
           
           !  re-assign obsgeoms_input w/ read-in data
           VBRDF_Sup_In%BS_USER_OBSGEOMS(1,:)=obsgeoms(:)
           VBRDF_Sup_In%BS_BEAM_SZAS(1)=obsgeoms(1)
           VBRDF_Sup_In%BS_USER_ANGLES_INPUT(1)=obsgeoms(2)
           VBRDF_Sup_In%BS_USER_RELAZMS(1)=obsgeoms(3)
           VBRDF_Sup_In%BS_WINDSPEED = aer_data_rec%wsp

           
           
           !  VLIDORT BRDF call
           VBRDF_Sup_In%BS_WAVELENGTH=aer_data_rec%wvs(iw)*0.001d0 ! nm -> micron
           
           trawl = .true. ; j = 0
           do while (trawl .and. j.lt.maxiter )
              set_Lw = .false. ; j = j + 1 ; if (j.eq.1)  Set_Lw = .true.



              CALL VBRDF_MAINMASTER ( &
                DO_DEBUG_RESTORATION,    & ! Inputs
                BS_NMOMENTS_INPUT,       & ! Inputs
                VBRDF_Sup_In,            & ! Inputs
                VBRDF_Sup_Out,           & ! Outputs
                VBRDF_Sup_OutputStatus )   ! Output Status

           !  Exception handling (turn this off for operational )
        
        

           !  Step 4b. Copy VBRDF results to VLIDORT Tyope structures
           !  -------------------------------------------------------
           
           VLIDORT_Sup%BRDF%TS_BRDF_F_0        = VBRDF_Sup_Out%BS_BRDF_F_0
           VLIDORT_Sup%BRDF%TS_BRDF_F          = VBRDF_Sup_Out%BS_BRDF_F
           VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = VBRDF_Sup_Out%BS_USER_BRDF_F_0
           VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = VBRDF_Sup_Out%BS_USER_BRDF_F
           VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC ! new but ...
           VLIDORT_Sup%BRDF%TS_EMISSIVITY       = ZERO
           VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY  = ZERO
           

           VSLEAVE_Sup_In%SL_WAVELENGTH=aer_data_rec%wvs(iw)*0.001d0 !
           

           if (set_Lw) then 

              CALL VSLEAVE_MAINMASTER ( aos_model, ABTOT_CONTROL, &
                   VSLEAVE_Sup_In,          & ! Inputs
                   VSLEAVE_Sup_Out,         & ! Outputs
                   VSLEAVE_Sup_OutputStatus ) ! Output Status
              if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
           
                 TRACE = 'VSLEAVE_MAINMASTER failed'
                 write(*,'(1x,a)') trim(TRACE)
                 
                 write(*,*)'Number of error messages from VSLEAVE calculation = ',&
                      VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
                 write(*,*) VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)
                 do i = 1, VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
                    write(*,'(a,i2,a,a)')' Message # ', i, ': ',&
                         VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)

                 enddo
                 
              endif

              
              CALL SET_VLIDORT_VSLEAVE_INPUTS ( &
                   VSLEAVE_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn,  & !Inputs
                   VLIDORT_Sup )
 


 

              LW_EXTERNAL  = VSLEAVE_Sup_Out%SL_LW_BASIC (1,1)
              K_DIRECT     = VSLEAVE_Sup_Out%SL_KSCALE_USERANGLES(1,1,1,1)
              K_SOLARVERT     = VSLEAVE_Sup_Out%SL_KSCALE_SOLARVERT(1,1)
              K_ISOTROPIC  = VSLEAVE_Sup_Out%SL_KSCALE_ISOTROPIC (1,1)
              print *, 'K_direct', K_DIRECT, 'K_isotropic', K_ISOTROPIC, 'K_solar',K_SOLARVERT

           endif

           IF ( .not. set_Lw ) then
              VLIDORT_Sup%Sleave%TS_LW_BASIC(1,1)  = LW_External
           endif

   
           
           VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,:)=obsgeoms(:)
           VLIDORT_ModIn%MSunRays%TS_SZANGLES(1)=obsgeoms(1)
           VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1)=obsgeoms(2)
           VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1)=obsgeoms(3)
           
           if (VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS .eq. 2) then 
              VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(2)=VLIDORT_FixIn%Cont%TS_nlayers
           endif
           

           VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT=height_grid(nlevels)

           
           
           
           VLIDORT_FixIn%Optical%TS_atmos_wavelength     = 0.001d0*aer_data_rec%wvs(iw) ! (mu)
           VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
           VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input(iw,:)
           VLIDORT_FixIn%Optical%TS_greekmat_total_input = &
                greekmat_total_input(:,iw,:,:)
           VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input(iw,:)



           do w = 1, 2

              !  Set the two calls.
              
              if ( w.eq.1 ) then
                 !         write(*,*)'PART I  - Iterative WL calculation, Baseline'
                 LW_BASIC_SAVED = VLIDORT_Sup%Sleave%TS_LW_BASIC(1,1)
              else if ( w.eq.2 ) then
                 !         write(*,*)'PART II - Iterative WL calculation, Use perturbed LW basic'
                 VLIDORT_Sup%Sleave%TS_LW_BASIC(1,1) = VLIDORT_Sup%Sleave%TS_LW_BASIC(1,1) * epsfac
              endif

        
              DO_DEBUG_INPUT = .true.
              

              open(87,file="VSLEAVE_input.nml")
              write(87,nml=VSLEAVE)

              open(88,file="VBRDF_input.nml")
              write(88,nml=VBRDF)

              open(89,file="VLIDORT_input.nml")
              write(89,nml=VLIDORT)

              !open(81,file="MERRA2_input.nml")
              !write(81,nml=MERRA2)



              CALL VLIDORT_MASTER ( DO_DEBUG_INPUT, & 
                   VLIDORT_FixIn, &
                   VLIDORT_ModIn, &
                   VLIDORT_Sup,   &
                   VLIDORT_Out )
              
              CALL VLIDORT_WRITE_STATUS ( &
                   'VLIDORT_Execution.log', &
                   95, OPENFILEFLAG, &
                   VLIDORT_Out%Status )
              
              trans_atmos = VLIDORT_Out%Main%TS_TRANS_ATMOS_FINAL(1)
              I_base_TOA(w)=VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
              !print *, v
              boaT=VLIDORT_Out%Main%TS_flux_diffuse(2,1,1,dnidx) + &! diffuse dwelling
                   VLIDORT_Out%Main%TS_DNFLUX_DIRECT(2,1,1) !,dnidx) ! direct dwelling
              wl_rad=VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(1,1,1,1)

           simrad_0 = I_base_TOA(1)
           simrad = I_base_TOA(2)
           simjac = ( I_base_TOA(2) - I_base_TOA(1) ) / eps / LW_BASIC_SAVED

           dy = e_mobyrad - simrad
           Lw_1 = LW_External + ( dy / simjac)
           if ( abs((LW_External/Lw_1) - 1.0d0) .lt. conv ) trawl = .false.


        enddo

        print *, 'VZA ',aer_data_rec%vza
        write(14,'(F5.1,11(",",1pe12.5))')  aer_data_rec%wvs(iw), LW_BASIC_SAVED,LW_BASIC_SAVED*epsfac, wl_rad,  &
              K_DIRECT, K_SOLARVERT,K_ISOTROPIC, trans_atmos, simjac, simrad_0,simrad, boaT

        !write(14,'(i4,2(",",i2),",",i5,",",i4,",",i2,",",F12.7,",",F12.7,",",F12.7,9(",",1pe16.8))') aer_data_rec%year, & 
        !     aer_data_rec%month,aer_data_rec%day,aer_data_rec%orbit, aer_data_rec%line, aer_data_rec%row, &
        !     sum(aer_data_rec%dtau(iw,:)), aer_data_rec%vza,aer_data_rec%wvs(iw), LW_BASIC_SAVED,Lw_External, wl_rad,  &
        !     Lw_1*trans_atmos, K_DIRECT, K_ISOTROPIC, trans_atmos, simjac, simrad_0,simrad, boaT

        
        write(*,'(/A,i3)')' Iteration # ',j
        write(*,'(A,1pe16.8)')'   -- Atmospheric Transmittance ',trans_atmos
        write(*,'(A,1pe16.8)')'   -- Lw Adjusted value = ',wl_rad
        write(*,'(A,1pe16.8)')'   -- Lw Basic value (used in retrieval) = ',Lw_External
        write(*,'(A,1pe16.8)')'   -- Lw Basic value * trans_atmos (used in retrieval) = ',Lw_External*trans_atmos
        write(*,'(A,1pe16.8)')'   -- Jacobian = ',simjac
        write(*,'(A,1pe16.8)')'   -- Imeas - Isim: ',dy
        write(*,'(A,1pe16.8)')'   -- [Imeas - Isim]/dy ',dy/simjac
        write(*,'(A,1pe16.8)')'   -- TOA Radiance for this iteration   = ',simrad
        write(*,'(A,1pe16.8)')'   -- Es for this iteration   = ',boaT
        write(*,'(A,1pe16.8)')'   -- Apply inversion, new Lw_dir Basic value = ',Lw_1*trans_atmos


        LW_External = Lw_1



        enddo ! nwv

     enddo
  
     close(14)
     stop 'VLIDORT run was successful'
  call ESMF_Finalize(__RC__)
end program vlidort_test_gas_retrieval
