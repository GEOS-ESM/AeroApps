MODULE generate_vlidort_iops

  
  PRIVATE


  public :: generate_vlidort_iops


contains

subroutine generate_vlidort_iops (deltau_rayleigh ,lambdas, nlambdas,nlayers,ngases, & ! input
             layer_gascolumns, gas_xsec, o4_xsec,layer_o4columns,& ! input
             greekmat_ray_input,greekmat_aer_input,& ! input
             do_aerosols, aerosol_deltau, aerosol_ssalbs, & ! Inputs
   deltau_vert_input,omega_total_input,greekmat_total_input ) ! output

   USE VLIDORT_PARS_m, only: MAXMOMENTS_INPUT, MAXSTOKES_SQ, fpk

   integer, INTENT(IN)  ::  nlambdas, nlayers, ngases
   real             :: lambdas(nlambdas) !, lambda_start
!  trace gas absorbers
   double precision :: gas_xsec(:,:),layer_gascolumns(:)
!  O2-O2 xsec
   double precision :: o4_xsec(:,:),layer_o4columns(:)
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

End  module generate_vlidort_iops      

