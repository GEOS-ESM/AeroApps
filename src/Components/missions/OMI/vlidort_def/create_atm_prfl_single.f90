MODULE create_atm_prfl_single

! Description: generate opd and ssa profile for Rayleigh+aerosol atm
!              for each OMI pixel based on P & T profiles at pxl location
! Revision History:
!     initial version: 2017-March -01
!     1st revised ver: 2017-Nov-29 (added background aerosols)
!     2nd revised ver: 2019-May-20 (added pixel-based aerosols)
!     3rd revised ver: 2019-May-29 (added hyperspectral wvs)
!     complete history: svn log file://omi/cm/svn/app/

! Author: Wenhan Qin (SSAI)
!*******************************************************************************

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            create_atmos_iops (master)                       #
! #                                                             #
! ###############################################################

  
  PRIVATE

  public :: prep_profile_iops
  public :: ozone_xsec
  public :: generate_vlidort_iops
  
  public :: o2o2_xsec_fspecm

! public :: rd_GMI_T

contains

!*******************************************************************************
subroutine o2o2_xsec_fspecm (lambdas, nlambdas,T,nlev, & ! input
           o4_xsec)  ! output

  integer, INTENT(IN)  ::  nlambdas, nlev
  real, INTENT(IN) :: T(nlev),lambdas(nlambdas)
  real(4), INTENT(OUT) :: o4_xsec(nlambdas,nlev)
  
  !  some constants
  real(8), parameter :: const=1.0D-10 ! cm^5/mole^2 -> m^5/mole^2
  integer, parameter :: nt=5 ! # of Ts for o2o2_xsec
  !  real, parameter :: eps=0.001
  
  !  Other help variables
  integer         :: N, L, W
  real, parameter :: fval=-9.99900E+10
  REAL*8          :: xsec(nlambdas,nt),xsec_inp(nt),d
  REAL            :: kd,Tx(nt) !, xsec_wv(nlambdas)
  character(len=100) :: junk

  data Tx /203.0, 233.0, 253.0, 273.0, 293.0/
  
  !  open(10,file='pace_o2o2_Thalman_Volkamer_273K.dat',status='old') ! [cm^5 molecule^-2]
  write(*,*) 'reading o2o2 xsec file: o4_Thalman_Volkamer_Tall_5nm.xs'
  open(10,file='o4_Thalman_Volkamer_Tall_5nm.xs',status='old') ! [cm^5 molecule^-2]
  read(10,*) junk
  read(10,*) junk
  read(10,*) junk
  do w=1,nlambdas
     read(10,*) kd,(xsec(w,l),l=1, nt) ! make sure wvs in xsec match lambdas
  enddo
  close(10)
  where(abs(xsec) .ge. 0.9*abs(fval)) xsec=0.0d0 ! zero out fillval: don't need this if fval<0
  
  do w = 1,nlambdas
     if (xsec(w,1) .gt. 0.0d0) then ! no zero xsec
        xsec_inp=xsec(w,:) ! for all NTs
     else
        o4_xsec(w,:)=0.0d0
        cycle ! zero o2o2 xsec, next w
     endif
     do n=1,nlev
        kd=T(n)
        L = minloc(abs(kd-Tx),dim=1)
        if (kd .lt. Tx(L)) then ! below Tx(L)
           if (L .eq. 1) then
              o4_xsec(w,n)=xsec_inp(L) ! below Tx(1)
              cycle
           else
              d=(xsec_inp(L)-xsec_inp(L-1))/(Tx(L)-Tx(L-1))
              o4_xsec(w,n)=d*(kd-Tx(L-1))+xsec_inp(L-1)
           endif
        else if (kd .gt. Tx(L)) then ! above Tx(L)
           if (L .eq. nt) then
              o4_xsec(w,n)=xsec_inp(L) ! above Tx(nt)
              cycle
           else
              d=(xsec_inp(L+1)-xsec_inp(L))/(Tx(L+1)-Tx(L))
              o4_xsec(w,n)=d*(kd-Tx(L))+xsec_inp(L)
           endif
        else ! equal Tx(L)
           o4_xsec(w,n)=xsec_inp(L)
        endif
     enddo ! done nlev
  enddo ! done nlambdas
  o4_xsec=const*o4_xsec ! -> m^5/mole^2
  return
end subroutine o2o2_xsec_fspecm

  subroutine ozone_xsec (lambdas, nlambdas,oz_xsec,wvbrc)  ! output
    integer, INTENT(IN)  ::  nlambdas
    real             :: lambdas(nlambdas), wvbrc(2,nlambdas),wx
    integer, parameter :: nloc=21
    real, parameter :: p0=1013.25
    integer :: nwlvac
    double precision, INTENT(OUT) :: oz_xsec(2,nlambdas,3)
    real(8), allocatable :: c(:,:)
    real,allocatable :: wlvac(:) ! wl=A
    character(len=4) :: loc_id(nloc)
    
    !  Other help variables
    integer         :: N, L, W
    REAL            :: kd


    nwlvac=43
    wvbrc=0.0 ! not needed
    allocate(c(nwlvac,3),wlvac(nwlvac))

    open(10,file='o3abs_brion_195_660_vacfinal.dat_atmcm_5nm',status='old')
    write(*,*) 'reading oz xsec file: o3abs_brion_195_660_vacfinal.dat_atmcm_5nm'

    do w=1,nwlvac
     read(10,*) wlvac(w),c(w,1),c(w,2),c(w,3)
    enddo

    ! find the lambdas start and end wv posts in wlvac
    L=minloc(abs(lambdas(1)-wlvac),dim=1)
    n=minloc(abs(lambdas(nlambdas)-wlvac),dim=1)

    do w = 1,nlambdas
       oz_xsec(1,w,:) = c(L,:)
       oz_xsec(2,w,:) = c(L,:)
       L=L+1
    enddo
   return
 end subroutine ozone_xsec
    
 





  subroutine prep_profile_iops (lambdas, nlambdas,NSTOKES, nlay,& ! input
           do_aerosols,aerosol_nscatmoms,AEROSOL_SCATMOMS, & ! input
           greekmat_ray_input,greekmat_aer_input,ngreek_moments_input,&
           rayleigh_xsec,rayleigh_depol) ! output

!  purposes -- determine parameters independent of layer # 

   USE VLIDORT_PARS_m, only: MAXMOMENTS_INPUT, MAXSTOKES_SQ, fpk

   real,      dimension(nlambdas)  :: lambdas
   real(fpk), dimension(nlambdas)  :: rayleigh_xsec,rayleigh_depol

!  CO2 PPMV mixing ratio (for the Rayleigh stuff)
   DOUBLE PRECISION  CO2_PPMV_MIXRATIO
   PARAMETER        ( CO2_PPMV_MIXRATIO = 385.0d0 ) ! by Eun-Su's suggestion

   INTEGER ::          W, K, L, nlambdas,nlay
   INTEGER ::          NSTOKES, NGREEKMAT_ENTRIES,NGREEK_MOMENTS_INPUT

   INTEGER ::          RK, CK, GK, SK, RMASK(8), GMASK(8),smask(8),CMASK(8)
   DOUBLE PRECISION :: DEPOL, BETA2
   DOUBLE PRECISION :: PROBLEM_RAY(6,0:2)
   DOUBLE PRECISION :: greekmat_ray_input ( 0:,:, : )
   DOUBLE PRECISION :: greekmat_aer_input ( 0:,:,:, : )

!  Aerosol control

   logical, INTENT(IN)  :: do_aerosols
   integer, INTENT(IN)  :: aerosol_nscatmoms
   REAL(fpk), INTENT(IN)  :: AEROSOL_SCATMOMS(:,0:,:,:)

!  Problem setup
!  -------------

!  Initialize IOPS variables

      RAYLEIGH_XSEC = 0.0d0
      RAYLEIGH_DEPOL = 0.0d0
      GREEKMAT_RAY_INPUT  = 0.0d0
      GREEKMAT_AER_INPUT  = 0.0d0
      PROBLEM_RAY = 0.0d0

!  Variables controlling placements and signs in the 4x4 Greek-matrix entries

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /) ! for Rayleigh
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /) ! for Aerosols from v_tester code

!  Set masking limits
      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8


!  Set some variables
      if ( do_aerosols ) then
        ngreek_moments_input = aerosol_nscatmoms
      else
!  Rayleigh scattering only, Number of expansion coefficients is always 2
        ngreek_moments_input = 2
      endif
!     write(*,*) '# of legendre terms=',aerosol_nscatmoms
!     write(*,*) 'NGREEK_MOMENTS_INPUT=',NGREEK_MOMENTS_INPUT

!  Rayleigh function

   call Rayleigh_function &
    ( nLAMBDAS, CO2_PPMV_MIXRATIO, &
!     NLAMBDAS,   LAMBDAS, &
      nlambdas,   lambdas, &
      RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

   rayleigh_xsec = 1.0D-4*rayleigh_xsec ! -> m^2

   !Set to 0 for testbed comparison
   !rayleigh_depol = 0.   

   do w = 1, nlambdas
     depol = rayleigh_depol(w)
     !depol  = 2.8622142884383864E-002 
     !print *, depol

!    write(*,*) 'calculated depol=',depol
     beta2 = ( 1.D0 - DEPOL ) / (2.D0 + DEPOL )
     PROBLEM_RAY(1,0) =  1.D0
     PROBLEM_RAY(4,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
     PROBLEM_RAY(1,2) =  beta2
     PROBLEM_RAY(5,2) =  - DSQRT(6.D0) * beta2
     PROBLEM_RAY(2,2) =  6.D0 * beta2

     DO K = 1, ngreekmat_entries
       GK = gmask(k); rk = rmask(k); sk = dble(smask(k)); ck = cmask(k) 
       do L = 0, ngreek_moments_input
         if (L .le. 2) greekmat_ray_input(L,w,gk) = PROBLEM_RAY(rk,L)
         if (do_aerosols) GREEKMAT_AER_INPUT(L,w,1:nlay,gk) = &
            sk*AEROSOL_SCATMOMS(ck,L,w,1:nlay)
       ENDDO
     enddo
     
   enddo ! all wvs

   return
end subroutine prep_profile_iops

!*******************************************************************************

!*******************************************************************************
subroutine generate_vlidort_iops (deltau_rayleigh ,lambdas, nlambdas,nlayers,ngases, & ! input
             layer_gascolumns, gas_xsec, o4_xsec,layer_o4columns,& ! input
             greekmat_ray_input,greekmat_aer_input,& ! input
             do_aerosols, aerosol_deltau, aerosol_ssalbs, & ! Inputs
   deltau_vert_input,omega_total_input,greekmat_total_input ) ! output

   USE VLIDORT_PARS_m, only: MAXMOMENTS_INPUT, MAXSTOKES_SQ, fpk

   integer, INTENT(IN)  ::  nlambdas, nlayers, ngases
   real             :: lambdas(nlambdas) !, lambda_start
!  trace gas absorbers
   REAL, INTENT(IN)  :: gas_xsec(:,:),layer_gascolumns(:), o4_xsec(:,:)
!  O2-O2 xsec
   REAL(fpk), INTENT(IN) :: layer_o4columns(:)
!  Air density and Rayleigh xsec

   REAL(fpk), intent(IN) :: GREEKMAT_RAY_INPUT (0:,:, :)
   REAL(fpk), intent(IN) :: deltau_rayleigh(:,:)
!  Aerosol control
   logical, INTENT(IN)  :: do_aerosols
!  Aerosol optical properties
!  logical, INTENT(IN)  :: AEROSOL_LAYERFLAGS(nlayers)
   REAL(fpk), INTENT(IN)  :: AEROSOL_DELTAU(:,:)
   REAL(fpk), INTENT(IN)  :: AEROSOL_SSALBS(:,:)
   REAL(fpk), intent(IN) :: GREEKMAT_AER_INPUT (0:,:,:, :)
!  output as input to vlidort
   REAL(fpk), intent(out) :: deltau_vert_input(:,:)
   REAL(fpk), intent(out) :: omega_total_input(:,:)
   REAL(fpk), intent(out) :: GREEKMAT_TOTAL_INPUT ( 0:,:,:, : )

!  Local variables
!  ===============

!  Limiting value for OMEGA (total single scattering albedo)
!  If there is no absorption, this is maximum value short of the prefect
!  scattering case(OMEGA = 1.0)

   real, parameter    :: omega_lim = 0.99999
   !real, parameter    :: omega_lim = 0.99999995000000000


!  Other help variables

   integer         :: N, L, G, W
   REAL(kind=8)    :: RAY_SCDEP, GAS_ABS, GAS_OPDEP, MOL_OPDEP, O4_OPDEP
   REAL(kind=8)    :: TOT_OPDEP, TOT_SCDEP, OMEGA

!  aerosol related variables
   logical      :: AERFLAG_N
   REAL(fpk)    :: AER_OPDEP, AER_SCDEP
   REAL(fpk)    :: RAY_WT, AER_WT, AER_SSALB

!  start of code

!  Initialize VLIDORT IOPS

   DELTAU_VERT_INPUT    = 0.0d0
   OMEGA_TOTAL_INPUT    = 0.0d0
   GREEKMAT_TOTAL_INPUT = 0.0d0

!  assuming aerosol is in every layer of nlayers
   AERFLAG_N     = DO_AEROSOLS ! .AND. AEROSOL_LAYERFLAGS(nlayers-N+1)

   do w = 1, nlambdas

!  Start layer loop

   DO N = 1, NLAYERS

!  Rayeigh scattering optical depth
      RAY_SCDEP = deltau_rayleigh(w,N)

!  O2-O2 optical depth

      O4_OPDEP = o4_xsec(W,N)*layer_o4columns(N)

!  Trace gases optical depths. GAS_OPDEP = Sum(GAS_ABS)

      GAS_OPDEP = 0.0D0
!  turn off ngases loop; so only works for ngases=1
      if (ngases .gt. 1) stop 'need to turn on ngases loop first !'
 
!     DO G = 1, NGASES
!        GAS_ABS      = LAYER_GASCOLUMNS(N,G) * gas_xsec(w,N,g)
         GAS_ABS      = LAYER_GASCOLUMNS(N) * gas_xsec(w,N)
         GAS_OPDEP    = GAS_OPDEP + GAS_ABS
!     ENDDO
!  write(*,*) 'GAS_OPDEP=',n,GAS_OPDEP
!  write(*,*) 'GAS col & xsec=',n,LAYER_GASCOLUMNS(N), gas_xsec(w,N)

!  Total molecular optical depth for the layer

      MOL_OPDEP =  RAY_SCDEP + GAS_OPDEP + O4_OPDEP

!  Total optical extinction and scattering depths for layer
!  These are just copies of the above, for the Rayleigh-only case

      TOT_OPDEP = MOL_OPDEP
      TOT_SCDEP = RAY_SCDEP

!  Layer flags

!     AERFLAG_N     = DO_AEROSOLS  .AND. AEROSOL_LAYERFLAGS(nlayers-N+1)

!  Add aerosols if present in this layer

      AER_OPDEP = 0.0d0
      AER_SCDEP = 0.0d0
      IF ( AERFLAG_N ) THEN
         AER_OPDEP = AEROSOL_DELTAU(N,W) ! dtau toa -> boa
         AER_SSALB = AEROSOL_SSALBS(N,W)
         AER_SCDEP = AER_OPDEP * AER_SSALB
         TOT_OPDEP = AER_OPDEP + TOT_OPDEP
         TOT_SCDEP = AER_SCDEP + TOT_SCDEP
      ENDIF

!  Scatter weighting

      RAY_WT                  = RAY_SCDEP / TOT_SCDEP
      IF ( AERFLAG_N ) AER_WT = AER_SCDEP / TOT_SCDEP

!  VLIDORT Bulk properties. Use limiting (toggle) value for OMEGA

      deltau_vert_input(w,N) = TOT_OPDEP
      OMEGA                = TOT_SCDEP / TOT_OPDEP

      IF ( OMEGA .GT. OMEGA_LIM ) THEN
         OMEGA_TOTAL_INPUT(w,N) = OMEGA_LIM
      ELSE
         OMEGA_TOTAL_INPUT(w,N) = TOT_SCDEP / TOT_OPDEP
      ENDIF
!  check values
!     write(*,*) 'w,n,dtau,ssa',w,n,deltau_vert_input(w,N),OMEGA_TOTAL_INPUT(w,N)
!     if (n .gt. 6) stop 'check the values first'

!  Fill up the Greek-matrix of coefficients 
!  SCATTERING LAW moments, Rayleigh first

!     GREEKMAT_TOTAL_INPUT(:,w,N,:)=RAY_WT*GREEKMAT_RAY_INPUT(:,w,:)
! if define GREEKMAT_RAY_INPUT(0:2,:,:)
      GREEKMAT_TOTAL_INPUT(0:2,w,N,:)=RAY_WT*GREEKMAT_RAY_INPUT(0:2,w,:)
      !GREEKMAT_TOTAL_INPUT(0:2,w,N,:)=GREEKMAT_RAY_INPUT(0:2,w,:)

!  Add aerosols
      
      IF ( AERFLAG_N ) GREEKMAT_TOTAL_INPUT(:,w,N,:)= &
         GREEKMAT_TOTAL_INPUT(:,w,N,:)+AER_WT*GREEKMAT_AER_INPUT(:,w,N,:)

!  Phase function normalization

      GREEKMAT_TOTAL_INPUT(0,w,N,1) = 1.0d0

   ENDDO !  End layer loop

   ENDDO !  End wavelength loop

   RETURN
end subroutine generate_vlidort_iops

!*******************************************************************************
SUBROUTINE RAYLEIGH_FUNCTION &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS, &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.

      IMPLICIT NONE

!  Input arguments
!  ---------------
!  wavelength

      INTEGER          FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
!     DOUBLE PRECISION FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )
! use real to be consistent with the input from main prog
      real             FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )

!  CO2 mixing ratio

      DOUBLE PRECISION CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      DOUBLE PRECISION RAYLEIGH_XSEC  ( FORWARD_MAXLAMBDAS )
      DOUBLE PRECISION RAYLEIGH_DEPOL ( FORWARD_MAXLAMBDAS )

!  Local variables
!  ---------------

      INTEGER          W
      DOUBLE PRECISION MASS_DRYAIR
      DOUBLE PRECISION NMOL, PI, CONS
      DOUBLE PRECISION MO2,MN2,MARG,MCO2,MAIR
      DOUBLE PRECISION FO2,FN2,FARG,FCO2,FAIR
      DOUBLE PRECISION LAMBDA_C,LAMBDA_M,LPM2,LP2
      DOUBLE PRECISION N300M1,NCO2M1,NCO2
      DOUBLE PRECISION NCO2SQ, NSQM1,NSQP2,TERM
      DOUBLE PRECISION S0_A, S0_B
      DOUBLE PRECISION S1_A, S1_B, S1_C, S1_D, S1_E
      DOUBLE PRECISION S2_A
      DOUBLE PRECISION S3_A, S3_B, S3_C, S3_D, S3_E

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

!*******************************************************************************
!  End module

End  module create_atm_prfl_single

