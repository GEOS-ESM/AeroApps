      module VLIDORT_LambMod

      implicit NONE

!     PRIVATE

      PUBLIC  VLIDORT_LambInit ! Initialize once
      PUBLIC  VLIDORT_LambRun  ! Run for each profile
      PUBLIC  VLIDORT_LambFin  ! Finalize once
      PUBLIC  VLIDORT_LER      ! Lambertian Equivalent Reflectivity
      PUBLIC  VLIDORT_AI       ! Aerosol Index
      
      include 'VLIDORT.PARS'         ! VLIDORT constants
      include 'VLIDORT_INPUTS.VARS'  ! VLIDORT input variables
      include 'VLIDORT_BOOKKEEP.VARS'  
      include 'VLIDORT_RESULTS.VARS' ! VLIDORT output variables 


      type VLIDORT_Lamb

         logical         :: initialized = .false.
         real*8          :: wavelength          ! in [nm]
         real*8          :: albedo              ! surface albedo
         real*8          :: solar_zenith        ! in degrees
         real*8          :: relative_azymuth    ! in degrees
         real*8          :: sensor_zenith       ! in degrees
         integer         :: nMom                ! number of momemts read (phase function) 
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer ::   g(:)              ! aerosol asymmetry factor
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]
         real*8, pointer:: pmom(:,:,:)          ! components of the scattering phase matrix

       end type VLIDORT_Lamb

      Contains

!.............................................................................

      subroutine VLIDORT_LambInit (self, km, rc, scalar)
    
      type(VLIDORT_Lamb), intent(inout)        :: self
      integer,              intent(in)         :: km     ! number of atmospheric layers
      integer,              intent(out)        :: rc     ! error code
      

      logical, optional,    intent(in)  :: scalar ! if true, perform scalar calculation
                                                  !  default is VECTOR calculation

!                     ---

      logical :: vector



      if ( present(scalar) ) then
         vector = .not. scalar
      else
         vector = .true.
      end if
      rc = 0



!     Initialize input variables in commons
!     -------------------------------------
      call VLIDORT_INIT_CONTROL_VARS()
      call VLIDORT_INIT_SURFACE_VARS()
      call VLIDORT_INIT_THERMAL_VARS()
      call VLIDORT_INIT_MODEL_VARS()

!     Customize those parameters that do not vary from pixel to pixel
!     ---------------------------------------------------------------


!                         Modes of Operation
!                         ------------------

      DO_FULLRAD_MODE    = .true.  ! Do full Stokes vector calculation?
      DO_SSCORR_NADIR    = .false. ! Do nadir single scatter correction?
      DO_SSCORR_OUTGOING = .true.  ! Do outgoing single scatter correction?
      DO_SSFULL          = .false. ! Do Full-up single scatter calculation?
      DO_DBCORRECTION    = .true.  ! Do direct beam correction?
      DO_DOUBLE_CONVTEST = .true.  ! Perform double convergence test?

!                            Solar Sources
!                            -------------

      DO_SOLAR_SOURCES    = .true.     ! Include solar sources?
      DO_PLANE_PARALLEL   = .false.    ! Plane-parallel treatment of direct beam?
      DO_CHAPMAN_FUNCTION = .true.     ! Perform internal Chapman function calculation?
      DO_REFRACTIVE_GEOMETRY = .false. ! Beam path with refractive atmosphere?

!                         Performance Control
!                         -------------------

      DO_RAYLEIGH_ONLY     = .false. ! Rayleigh atmosphere only?
      DO_DELTAM_SCALING    = .true.  ! Include Delta-M scaling?
      DO_SSCORR_TRUNCATION = .false. ! Additional Delta-M scaling for SS correction?
      DO_SOLUTION_SAVING   = .false. ! Solution saving mode?
      DO_BVP_TELESCOPING   = .false. ! Boundary value problem telescoping mode?

!                      User-defined output control
!                      ---------------------------

      DO_UPWELLING = .true.     ! Upwelling output?
      DO_DNWELLING = .false.    ! Downwelling output?
      DO_USER_VZANGLES = .true. ! User-defined viewing zenith angles?

      DO_ADDITIONAL_MVOUT = .false.  ! Generate mean value output additionally?
      DO_MVOUT_ONLY       = .false. ! Generate only mean value output?

!                           Write Control
!                           -------------

      DO_DEBUG_WRITE    = .false. ! Debug write?
      DO_WRITE_INPUT    = .false. ! Input control write?
      DO_WRITE_SCENARIO = .false. ! Input scenario write?
      DO_WRITE_FOURIER  = .false. ! Fourier component output write?
      DO_WRITE_RESULTS  = .false. ! Results write?

      INPUT_WRITE_FILENAME    = '/dev/null' ! filename for input write
      SCENARIO_WRITE_FILENAME = '/dev/null' ! filename for scenario write
      FOURIER_WRITE_FILENAME  = '/dev/null' ! Fourier output filename
      RESULTS_WRITE_FILENAME  = '/dev/null' ! filename for main output

!                     Stokes/streams/layers/moments
!                     -----------------------------

      if ( vector ) then
         NSTOKES = 3        ! Number of Stokes vector components
      else
         NSTOKES = 1
      end if
      NSTREAMS = 6         ! Number of half-space streams
      NLAYERS = km          ! Number of atmospheric layers
      NFINELAYERS = 2       ! Number of fine layers (outgoing sphericity correction)
      NGREEK_MOMENTS_INPUT = 500 ! Number of scattering matrix expansion coefficients
 
!                            accuracy input
!                            --------------

      VLIDORT_ACCURACY = 0.0001  ! Fourier series convergence
      FLUX_FACTOR = 1.00d0      ! Solar flux constant; =1 if no solar sources.

!                        Pseudo-spherical inputs
!                        -----------------------

      EARTH_RADIUS = 6371.0 ! Earth radius (km)
      RFINDEX_PARAMETER = 0.000288 ! Refractive index parameter
      GEOMETRY_SPECHEIGHT = 0.0 ! Input geometry specification height [km]

!                       Solar/Viewing Geometry 
!                       ----------------------

      N_SZANGLES = 1       ! Number of BOA solar zenith angles 
      N_USER_RELAZMS = 1   ! Number of user defined relative azimuth angles    
      N_USER_VZANGLES = 1  ! Number of user-defined stream angles     
      N_USER_LEVELS = 1    ! Number of user-defined vertical output levels

!                          Lambertian Only
!                          ---------------

      DO_LAMBERTIAN_SURFACE = .true. ! Lambertian surface?
      LAMBERTIAN_ALBEDO = 0.0  ! Lambertian (isotropic) input albedo, will be reset
                                     !  during run
      

!                          Thermal controls
!                          ---------------- 

      DO_THERMAL_EMISSION  = .false.  ! Do thermal emission?


      DO_THERMAL_TRANSONLY = .false.  ! Do thermal emission, transmittance only?
      N_THERMAL_COEFFS     = 2        ! Number of thermal coefficients
      DO_SURFACE_EMISSION  =  .false. ! Do Surface emission?

!                            Error Checking
!                            --------------

      if ( NSTOKES  .GT. MAXSTOKES  )                   rc = 1
      if ( NSTREAMS .GT. MAXSTREAMS )                   rc = 2
      if ( NLAYERS  .GT. MAXLAYERS  )                   rc = 3
      if ( NFINELAYERS .GT. MAXFINELAYERS )             rc = 4
      if ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT)  rc = 5
      if ( N_SZANGLES .GT. MAX_SZANGLES )               rc = 6
      if ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS )       rc = 7
      if ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES )     rc = 8
      if ( N_USER_LEVELS .GT. MAX_USER_LEVELS )         rc = 9 
      if ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS )   rc =  10
     
!     All done
!     --------
      self%initialized = .true.

      end subroutine VLIDORT_LambInit

!.............................................................................
 
      subroutine VLIDORT_LambRun (self, radiance, reflectance, AOT, rc, 
     &                            scalar, aerosol)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      
      type(VLIDORT_Lamb), intent(inout)   :: self        ! Contains most input
      real*8,             intent(out)     :: radiance    ! TOA radiance
      real*8,             intent(out)     :: reflectance ! TOA reflectance
      real*8,             intent(out)     :: AOT         ! Aerosol Optical Thickness
      integer,            intent(out)     :: rc

      logical, optional,  intent(in)      :: scalar  ! If True, do scalar calculation
      logical,            intent(in)      :: aerosol ! if False, Rayleigh only

!                           ----

      integer :: STATUS_INPUTCHECK, STATUS_CALCULATION 
      logical :: vector

!                           ----

!     local variables
!     ---------------
      integer                                          :: i, j, k, l, m, n
      integer                                          :: ncoeffs, k1, IDR, ierror
      real*8                                           :: ray_l      
      real*8                                           :: tau_l 
      real*8                                           :: ssa_l 
      real*8                                           :: g_l 
      real*8                                           :: tau_ext
      real*8                                           :: tau_scat
      real*8                                           :: ssa_tot
      real*8                                           :: raysmom2
      real*8                                           :: gammamom2
      real*8                                           :: alphamom2 
      real*8                                           :: deltamom1 
      real*8                                           :: aerswt 
      real*8                                           :: rayswt
      real*8                                           :: factor      
      real*8                                           :: COEF_DEPOL
      real*8                                           :: wmicron
      real*8                                           :: x
      real*8, dimension(NLAYERS)                       :: Ray
      real*8, dimension(0:MAXMOMENTS_INPUT)            :: aersmom 
      real*8, dimension(0:MAXMOMENTS_INPUT,NLAYERS,16) :: aervmoms 
      real*8, dimension(0:2, 16)                       :: rayvmoms
      real*8                                           :: difz
      real*8                                           :: somray

!     Volume Rayleigh scattering coefficient (depends closely on the 
!     thermodynamic conditions of the atm and varies with height, 
!     product of the molecular number density of air by the total ray 
!     scattering cross section).
      real*8, dimension(NLAYERS+1)                     :: Vol 
      real*8, dimension(NLAYERS+1)                     :: sect 
      real*8, parameter                                :: pi = 3.14159
      real*8, parameter                                :: DEPOL_RATIO = 0.030
      
      
       rc = 0
 
       if ( .not. self%initialized ) then
         rc = 1
         return
       end if
      
       if ( present(scalar) ) then
         vector = .not. scalar
      else
         vector = .true.
      end if

      if ( DO_UPWELLING ) then
         IDR = 1
      else
         IDR = 2
      end if
     
      SZANGLES(1) = self%solar_zenith
      USER_RELAZMS(1) = self%relative_azymuth
      USER_VZANGLES(1) = self%sensor_zenith
      USER_LEVELS(1) = 0.0
      LAMBERTIAN_ALBEDO = self%albedo  
      
      print *, 'ALBEDO & SZA', self%wavelength,LAMBERTIAN_ALBEDO,SZANGLES(1) 
      
  
!               Calculation of the Rayleigh-Scattering Optical Depth
!               ----------------------------------------------------
      height_grid (0) = self%ze(1) * 1.E-3  ! en km
      pressure_grid (0) = self%pe(1) * 1.E-2 ! en hPa
      temperature_grid (0) = self%te(1)
      
      do i = 1, NLAYERS 
         height_grid (i) = self%ze(i+1) * 1.E-3 ! en km
         pressure_grid (i) = self%pe(i+1) * 1.E-2 ! en hPa
         temperature_grid (i) = self%te(i+1)
      end do
            
!     Rayleigh extinction profile from Bodhaine et al., (1999) 
!     and Tomasi et al., (2005)
!     (wavelength in micrometer, pressure in hpa)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wmicron = self%wavelength * 1.E-3  ! micrometer
!      print*, 'VLIDORT 2'
      do j = 0, NLAYERS 
        
         Vol(j) = 3.69296E-18 * A(wmicron) * pressure_grid(j) 
     &          / temperature_grid(j)
         Vol(j) = Vol(j) * 1.E5 ! en km-1

        sect(j) = Vol(j)/2.546899E19 * 1013.25/pressure_grid(j) 
     &          * temperature_grid(j)/288.15

!         print*,*, "sect", j,height_grid(j), sect(j) ! test de compa with Bodhaine
           
      end do
       
!     logarithmique interpolation procedure of the Rayleigh profile 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      somray = 0.
     
      do j = 0, NLAYERS-1
         difz = height_grid(j) - height_grid(j+1)
!         ray(j+1) = difz * Vol(j) * (Vol(j+1)/Vol(j))**0.5 
         ray(j+1) = difz * (Vol(j) + Vol(j+1))/2.
         somray = somray + ray(j+1) ! Total Tau Ray
      end do

!     greek moments for Rayleigh only( only if vector )
!     ----------------------------- 
      if ( vector ) then

         gammamom2 = -SQRT(6.) * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
         alphamom2 = 6 * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
         deltamom1 = 3 * (1 - 2 * DEPOL_RATIO) / (2 + DEPOL_RATIO)
         raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 

!DEPOL_RATIO = (6 * COEF_DEPOL - 6) / (3 + 7 * COEF_DEPOL) 
      
         rayvmoms(0,1) = 1.0
         rayvmoms(1,1) = 0.0
         rayvmoms(2,1) = raysmom2
         rayvmoms(0,2) = 0.0
         rayvmoms(1,2) = 0.0
         rayvmoms(2,2) = gammamom2
         do k = 3, 4
            do l = 0, 2
               rayvmoms(l,k) = 0.0
            end do
         end do
         rayvmoms(0,5) = 0.0
         rayvmoms(1,5) = 0.0
         rayvmoms(2,5) = gammamom2 
         rayvmoms(0,6) = 0.0
         rayvmoms(1,6) = 0.0
         rayvmoms(2,6) = alphamom2 
         do k = 7, 15
            do l = 0, 2
               rayvmoms(l,k) = 0.0
            end do
         end do
         rayvmoms(0,16) = 0.0
         rayvmoms(1,16) = deltamom1
         rayvmoms(2,16) = 0.0 
     
      end if

!     Loop over the layers:
!     ---------------------
      AOT = 0.0
      do i = 1, NLAYERS  
         ray_l = ray(i)         ! indice l for  each layer

         if ( .not. aerosol ) then
            tau_l = 0.0 
            ssa_l = 0.0 
            g_l = 0.0 
         else
            tau_l = self%tau(i)
            ssa_l = self%ssa(i) 
            g_l = self%g(i) 
            AOT = AOT + tau_l
         end if

!        total optical depths for extinction and scattering 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tau_ext = ray_l + tau_l
         tau_scat = ray_l +  ssa_l * tau_l
         
!        single scattering albedo total
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ssa_tot = tau_scat / tau_ext
         if ( ssa_tot > 0.99999 ) then
            ssa_tot = 0.99999
         end if
     
         deltau_vert_input(i) = tau_ext
         omega_total_input(i) = ssa_tot 

!            Compute Henyey-Greenstein phase function, including Rayleigh
!            ------------------------------------------------------------

!     SCALAR testing Rayleigh second moment
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wmicron = self%wavelength * 1.E-3  ! micrometer
     
      raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 
      if (scalar) then    
    
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (tau_l == 0.0) then
     
            greekmat_total_input(0,i,1) = 1.0
            greekmat_total_input(1,i,1) = 0.0
            greekmat_total_input(2,i,1) = raysmom2
            do l = 3, NGREEK_MOMENTS_INPUT       
               greekmat_total_input(l,i,1) = 0.0
            end do 
         end if
     

!     SCALAR phase function moments (aerosol + Rayleigh)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (tau_l /= 0.0) then
             aerswt = ssa_l * tau_l / tau_scat  
             rayswt = ray_l / tau_scat    
             aersmom(0) = 1.0
             aersmom(1) = 3.0 * g_l
             aersmom(2) = 5.0 * g_l * aersmom(1) / 3.0

             greekmat_total_input(0,i,1) = 1.0
             greekmat_total_input(1,i,1) = aersmom(1) * aerswt
             greekmat_total_input(2,i,1) = raysmom2 * rayswt + aersmom(2) * aerswt
         
             do l = 3, NGREEK_MOMENTS_INPUT        
                factor = REAL(2*l+1) / REAL(2*l-1) 
                aersmom(l) = factor * g_l * aersmom(l-1)
                greekmat_total_input(l,i,1) = aersmom(l) * aerswt
             end do  
          end if ! end if tau_l /= 0.0
      end if  ! end scalar
      
      

!     VECTOR phase function moments 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if ( vector ) then  
      
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (tau_l == 0.0) then
   
           do k=1, 16
             do l=0, 2
              greekmat_total_input(l,i,k) = rayvmoms(l,k)
             end do
           end do

            do k = 1, 16 
              do l = 3, NGREEK_MOMENTS_INPUT
             greekmat_total_input(l,i,k) = 0.0
             end do
           end do

      end if
      
      
!     VECTOR phase function moments (aerosol + Rayleigh)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (tau_l /= 0.0) then 
             
         aerswt = ssa_l * tau_l / tau_scat  
         rayswt = ray_l / tau_scat 

!         if (g_l > 0.9.and.g_l<=1.0) then
!            self%nmom = 500
!         else if (g_l > 0.8.and.g_l<=0.9) then
!           self%nmom = 200
!         else 
!            self%nmom = 100
!         end if
!        print*, 'ok1', self%nmom
             
         if ( self%nmom -1 < NGREEK_MOMENTS_INPUT ) then 
            do l= 0, self%nmom-1
               
               aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11             
               aervmoms(l,i,6)  = self%pmom(i,l+1,1) ! P11               
               aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12                
               aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12             
               aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33           
               aervmoms(l,i,16) = self%pmom(i,l+1,3) ! P33          
               aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34            
               aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34 = P44
            end do
         
            do l = self%nmom, NGREEK_MOMENTS_INPUT
               aervmoms(l,i,1)  = 0.0 ! P11
               aervmoms(l,i,6)  = 0.0 ! P11 = P22
               aervmoms(l,i,2)  = 0.0 ! P12
               aervmoms(l,i,5)  = 0.0 ! P12
               aervmoms(l,i,11) = 0.0 ! P33
               aervmoms(l,i,16) = 0.0 ! P33
               aervmoms(l,i,12) = 0.0 ! P34
               aervmoms(l,i,15) = 0.0 ! P34 = P44
        
            end do
        
         else 
         
            do l= 0, NGREEK_MOMENTS_INPUT
               aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11
               aervmoms(l,i,6)  = self%pmom(i,l+1,1) ! P11
               aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12
               aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12
               aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33
               aervmoms(l,i,16) = self%pmom(i,l+1,3) ! P33
               aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34
               aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34 = P44
            end do
       
         end if

         do k = 1, 16
            do l = 0,2
               greekmat_total_input(l,i,k) = rayvmoms(l,k) * rayswt 
     &                                     + aervmoms(l,i,k) * aerswt
            end do
           
            do l = 3, NGREEK_MOMENTS_INPUT
               greekmat_total_input(l,i,k) = aervmoms(l,i,k) * aerswt
            end do
                  
         end do
         greekmat_total_input(0,i,1) = 1.0
      
     
      end if                    ! end if tau_l /= 0.0

      end if                    ! end if vector
      
!     end layer loop
!     ---------------
      end do
      
           
!     Call the Lambertian driver for doing the actual calculation
!     -----------------------------------------------------------
      call VLIDORT_Master_Lambertian ( STATUS_INPUTCHECK, STATUS_CALCULATION )
      if ( STATUS_INPUTCHECK /= 0 ) rc = 1
      if ( STATUS_CALCULATION /= 0 ) rc = 2
      if ( rc /= 0 ) return

!     Return TOA radiance
!     -------------------
       RADIANCE = 0.0
       REFLECTANCE = 0.0
       if ( scalar ) then
       RADIANCE = STOKES(1, 1, 1, IDR)             
       REFLECTANCE = (pi * RADIANCE) / ( cos(SZANGLES(1)*pi/180.0) * FLUX_FACTOR ) 
       end if
   
       if ( vector ) then
       RADIANCE = STOKES(1, 1, 1, IDR)  
       REFLECTANCE = (pi * RADIANCE) / ( cos(SZANGLES(1)*pi/180.0) * FLUX_FACTOR )   
       end if
       
       print *, 'RADIANCE',self%wavelength, RADIANCE   
      Contains  

!.................................................................................
! determination of the function A(wmicron)
!-----------------------------------      
    
         function A(X)
!  pour le Rayleigh, X = lambda en micron
       implicit none
       real*8 :: A
       real*8, intent(in) :: X
       real*8 :: depol_ratio
       real*8 :: XX, XX2, C
       real*8 :: depol1, depol2, depol3, depol4
       real*8 :: coef_depol
       real*8 :: RI

      XX=1./X

      XX2=XX*XX
     
!  ns-1  sans approximation
      RI=8060.77+2481070/(132.274-XX2)+17456.3/(39.32957-XX2)


! autre maniere de determiner coef_depol
        depol1 = 1.034 + 3.17E-4 * XX2
        depol2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2
        depol3 = 1.
        depol4 = 1.15
!depol5 = 1.001 if we 
        C = 0.030 ! for standard air , concentration of CO2 = 300 ppmv

        coef_depol= (78.084 * depol1 + 20.946 * depol2 + 0.934 * depol3 +  C * depol4) / (78.084 + 20.946 + 0.934 + C)

         
      A=coef_depol * XX2 * XX2 * RI * RI
      
      return

      end function A      
    
      end subroutine VLIDORT_LambRun

!.............................................................................

      subroutine VLIDORT_LambFin (self, rc)
      type(VLIDORT_Lamb), intent(inout) :: self
      integer,              intent(out) :: rc
      self%initialized = .false.
      
      rc = 0
      
      end subroutine VLIDORT_LambFin

!..............................................................................

      subroutine VLIDORT_LER (self, rad, refl, rc)
 
      ! First compute the spherical albedo and the transmission of the 
      ! atmosphere purely rayleigh and
      ! the Lambertian Equivalent Reflectivity R*
      !-------------------------------------------------------------------------------------

       type(VLIDORT_Lamb), intent(inout) :: self 
       real*8,  intent(in)               :: rad       !Radiance mesure (atm Ray + aer)
       real*8,  intent(out)              :: refl      !surface reflectivity R*           
       integer, intent(out)              :: rc
       real*8                            :: spher_alb !spherical albedo
       real*8                            :: trans     !transmission
       real*8                            :: reflectance 
       real*8                            :: beta
       real*8                            :: rad0      !Ray atmo and no surface reflection
       real*8                            :: rad1      !Ray atmo and surface albedo = 0.1
       real*8                            :: rad2      !Ray atmo and surface albedo = 0.2
       real*8                            :: rad_aer
       real*8                            :: rad_ray
       real*8                            :: delta_rad
       real*8                            :: AOT
      
       
      
       ! atmosphere Rayleigh and surface albedo = 0.0
       ! --------------------------------------------
       self%albedo = 0.0
       call VLIDORT_LambRun (self, rad0, reflectance, AOT, rc, .false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.1
       ! --------------------------------------------
       self%albedo = 0.1
       call VLIDORT_LambRun (self, rad1, reflectance, AOT, rc, .false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.2
       ! --------------------------------------------
       self%albedo = 0.2
       call VLIDORT_LambRun (self, rad2, reflectance, AOT, rc, .false., .false.)

       beta =  ( rad1 - rad0 ) / ( rad2 - rad0 )
       spher_alb = (( beta * 0.2 ) - 0.1 ) / ( 0.1 * 0.2 * ( beta - 1 )) ! Spherical albedo
       trans = ( 1 - ( spher_alb * 0.1 )) * ( rad1 - rad0 ) / 0.1        ! transmission
       
!       print*,*, 'spher_alb', spher_alb, trans
 
       ! Compute the Reflectivity R*
       ! --------------------------
       ! Compute the radiance for an atmosphere purely Ray and surface albedo = 0.0
       !----------------------------------------------------------------------------
       self%albedo = 0.0
       call VLIDORT_LambRun (self, rad_ray, reflectance, AOT, rc, .false., .false.)
!       print*,*, 'rad_ray', rad_ray      

       delta_rad =  rad - rad_ray    
        
       refl = delta_rad / ( trans + ( spher_alb * delta_rad )) ! Reflectivity R*

!       print*,*, 'refl', refl

       end subroutine VLIDORT_LER
       
!.........................................................................

      subroutine VLIDORT_AI (self, rad, refl, AI, rc)   

      ! Compute the Aerosol Index from the reflectivity calculate 
      ! from the previous subroutine
      ! ---------------------------------------------------------
       type(VLIDORT_Lamb), intent(inout)  :: self
       real*8,             intent(in)     :: rad
       real*8,             intent(in)     :: refl      !reflectivity R*
       integer,            intent(out)    :: rc
       real*8,             intent(out)    :: AI        !Aerosol Index
       real*8                             :: rad_nouv  !calcul with R* as the surface albedo param
       real*8                             :: reflectance 
       real*8                             :: alpha  
       
       real*8                             :: AOT 
       real*8                             :: g_tot
      ! Compute the radiance for an atmosphere Ray with the R* as the surface alb parameter
      !-----------------------------------------------------------------------------------

      self%albedo = refl  
      
      call VLIDORT_LambRun (self, rad_nouv, reflectance, AOT, rc, .false., .false.)

      if (rad > 0.0) then
      alpha = rad / rad_nouv
      AI = -100 * log10(alpha) ! Calcul of the Aerosol Index
      else
      AI = self%MISSING
      end if
     
     
      
      end subroutine VLIDORT_AI


!................................................................................
! function coef_depolarisation ( KING FACTOR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      function COEF_DEPOL(X)
!      implicit none
!      real*8              :: COEF_DEPOL ! King Factor
!      real*8              :: DEPOL_RATIO ! depolarisation ratio
!      real*8, intent(in)  :: X
!      real*8              :: XX, XX2, C
!      real*8              :: DEPOL1, DEPOL2, DEPOL3, DEPOL4
      
!      XX=1./X

!      XX2=XX*XX
     
!     Calculation of the King Factor
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      DEPOL1 = 1.034 + 3.17E-4 * XX2
!      DEPOL2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2
!      DEPOL3 = 1.
!      DEPOL4 = 1.15
!      C = 0.030 ! for standard air , concentration of CO2 = 300 ppmv

!      COEF_DEPOL= (78.084 * DEPOL1 + 20.946 * DEPOL2 + 0.934 * DEPOL3 +  C * DEPOL4) / (78.084 + 20.946 + 0.934 + C)
    
!      return

!      end function COEF_DEPOL

      end module VLIDORT_LambMod
