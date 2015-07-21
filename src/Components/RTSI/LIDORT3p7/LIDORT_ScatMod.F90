      module LIDORT_ScatMod
 
      USE LIDORT_PARS
      
      USE BRDF_SUP_MOD
      USE LIDORT_IO_DEFS
      
      USE LIDORT_AUX
      USE LIDORT_INPUTS
      USE LIDORT_MASTERS

      USE LIDORT_Mod
      USE LIDORT_SurfaceMod

      implicit NONE

      PUBLIC  LIDORT_Run      ! Run for each profile (pixel)
      PUBLIC  LIDORT_LER      ! Lambertian Equivalent Reflectivity
      PUBLIC  LIDORT_AI       ! Aerosol Index
          

      type LIDORT_scat

!         logical         :: initialized = .false.
         real*8          :: wavelength          ! in [nm]
         integer         :: nMom                ! number of momemts read (phase function) 
         integer         :: nPol                ! number of components of the scattering matrix
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer ::   g(:)              ! aerosol asymmetry factor
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]
         real*8, pointer:: pmom(:,:,:)          ! components of the scattering phase matrix

         type(LIDORT_Surface) :: Surface

      end type LIDORT_scat
       
      Contains
!.............................................................................
 
      subroutine LIDORT_Run (self,radiance, reflectance, ROT, aerosol,rc)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE LIDORT_PARS
         
      USE LIDORT_IO_DEFS
      USE BRDF_SUP_MOD
      
      USE LIDORT_AUX
      USE LIDORT_INPUTS
      USE LIDORT_MASTERS
     

      type(LIDORT_scat), intent(inout)   :: self        ! Contains most input
      real*8,             intent(out)     :: radiance    ! TOA radiance
      real*8,             intent(out)     :: reflectance ! TOA reflectance
      real*8,dimension(:),intent(inout)   :: ROT         ! Rayleigh Optical Thickness
      integer,            intent(out)     :: rc

      logical,            intent(in)      :: aerosol ! if False, Rayleigh only

!                           ----

      integer :: STATUS_INPUTCHECK, STATUS_CALCULATION 

!                           ----

!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: ncoeffs, k1, IDR, ierror
      integer                                            :: NLAYERS
      real*8                                             :: ray_l      
      real*8                                             :: tau_l 
      real*8                                             :: ssa_l 
      real*8                                             :: g_l 
      real*8                                             :: tau_ext
      real*8                                             :: tau_scat
      real*8                                             :: ssa_tot
      real*8                                             :: raysmom2
      real*8                                             :: gammamom2
      real*8                                             :: alphamom2 
      real*8                                             :: deltamom1 
      real*8                                             :: aerswt 
      real*8                                             :: rayswt
      real*8                                             :: factor      
      real*8                                             :: COEF_DEPOL
      real*8                                             :: wmicron
      real*8                                             :: AOT         ! Aerosol Optical Thickness


   
      real*8                                             :: x
      real*8, dimension(MAXLAYERS)                       :: Ray
      real*8, dimension(0:MAXMOMENTS_INPUT)              :: aersmom 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms 
      real*8, dimension(0:2, 16)                         :: rayvmoms
      real*8                                             :: difz
      real*8                                             :: somray
      logical                                            :: DO_BRDF_SURFACE
      real*8                                             :: LAMBERTIAN_ALBEDO
      integer                                            :: NGREEK_MOMENTS_INPUT
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS)    :: greekmat_total_input
      real*8, dimension(0:MAXLAYERS)                     :: height_grid                     
      real*8, dimension(0:MAXLAYERS)                     :: pressure_grid 
      real*8, dimension(0:MAXLAYERS)                     :: temperature_grid

      real*8, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, dimension(MAXLAYERS)                       :: omega_total_input

      real*8,dimension(MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)  :: INTENSITY
      real*8                                                           :: FLUX_FACTOR

!     Volume Rayleigh scattering coefficient (depends closely on the 
!     thermodynamic conditions of the atm and varies with height, 
!     product of the molecular number density of air by the total ray 
!     scattering cross section).
      real*8, dimension(MAXLAYERS+1)                     :: Vol 
      real*8, dimension(MAXLAYERS+1)                     :: sect 
      real*8, parameter                                  :: pi = 4.*atan(1.0)
      real*8, parameter                                  :: DEPOL_RATIO = 0.030
       
      
       rc = 0
 
      if ( .not. self%Surface%Base%initialized ) then
        rc = 1
        return
      end if
      
      
      if ( self%Surface%Base%VIO%LIDORT_FixIn%Bool%TS_DO_UPWELLING ) then
         IDR = 1
      else
         IDR = 2
      end if
     
!                          Lambertian OR BRDF surface
!                          --------------------------

      if ( self%Surface%sfc_type == 1 ) then               ! Lambertian surface ?
         DO_BRDF_SURFACE = .false.        
         LAMBERTIAN_ALBEDO = self%Surface%albedo           ! Lambertian (isotropic) input albedo
      else                                   
         DO_BRDF_SURFACE = .true. 
         LAMBERTIAN_ALBEDO = 0.0                           ! Use BRDF -> albedo set to 0
 
      end if

      self%Surface%Base%VIO%LIDORT_Fixin%Bool%TS_DO_BRDF_SURFACE = DO_BRDF_SURFACE
      self%Surface%Base%VIO%LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = LAMBERTIAN_ALBEDO
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_BRDF_F_0        = self%Surface%Base%VIO%BRDF_Sup_Out%BS_BRDF_F_0
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_BRDF_F          = self%Surface%Base%VIO%BRDF_Sup_Out%BS_BRDF_F
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = self%Surface%Base%VIO%BRDF_Sup_Out%BS_USER_BRDF_F_0
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_USER_BRDF_F     = self%Surface%Base%VIO%BRDF_Sup_Out%BS_USER_BRDF_F
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = self%Surface%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_EMISSIVITY      = self%Surface%Base%VIO%BRDF_Sup_Out%BS_EMISSIVITY
      self%Surface%Base%VIO%LIDORT_Sup%BRDF%TS_USER_EMISSIVITY = self%Surface%Base%VIO%BRDF_Sup_Out%BS_USER_EMISSIVITY

!                         Angles (SZA, viewing, relatuve azimuth), Level
!                        ------------------------------------------------      
      self%Surface%Base%VIO%LIDORT_ModIn%MSunRays%TS_BEAM_SZAS(1)           = self%Surface%solar_zenith
      self%Surface%Base%VIO%LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1)        = self%Surface%relat_azimuth
      self%Surface%Base%VIO%LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1)   = self%Surface%sensor_zenith
      self%Surface%Base%VIO%LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1)         = 0.0     ! This is TOA   

 
                      
      
      NGREEK_MOMENTS_INPUT = self%Surface%Base%VIO%LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT
      NLAYERS = self%Surface%Base%VIO%LIDORT_FixIn%Cont%TS_NLAYERS     
     
!               Calculation of the Rayleigh-Scattering Optical Depth
!               ----------------------------------------------------
   
      height_grid(0) = self%ze(1) * 1.E-3  ! en km
      pressure_grid(0) = self%pe(1) * 1.E-2 ! en hPa
      temperature_grid(0) = self%te(1)


      do i = 1, NLAYERS 
         height_grid(i) = self%ze(i+1) * 1.E-3 ! en km
         pressure_grid(i) = self%pe(i+1) * 1.E-2 ! en hPa
         temperature_grid(i) = self%te(i+1)
      end do
            
      self%Surface%Base%VIO%LIDORT_FixIn%Chapman%TS_height_grid = height_grid
      self%Surface%Base%VIO%LIDORT_FixIn%Chapman%TS_pressure_grid = pressure_grid
      self%Surface%Base%VIO%LIDORT_FixIn%Chapman%TS_temperature_grid = temperature_grid


!     Rayleigh extinction profile from Bodhaine et al., (1999) 
!     and Tomasi et al., (2005)
!     (wavelength in micrometer, pressure in hpa)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wmicron = self%wavelength * 1.E-3  ! micrometer
      
      do j = 0, NLAYERS 
        
         Vol(j) = 3.69296E-18 * A(wmicron) * pressure_grid(j) / temperature_grid(j)
         Vol(j) = Vol(j) * 1.E5 ! en km-1

        sect(j) = Vol(j)/2.546899E19 * 1013.25/pressure_grid(j) * temperature_grid(j)/288.15

          
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


!     Loop over the layers:
!     ---------------------
      AOT = 0.0
      ROT = 0.0
      do i = 1, NLAYERS  
         ray_l = ray(i)         ! indice l for  each layer
         ROT(i) = ray(i)
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

!        Compute Henyey-Greenstein phase function, including Rayleigh
!        ------------------------------------------------------------

!        SCALAR testing Rayleigh second moment
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         wmicron = self%wavelength * 1.E-3  ! micrometer
     
         raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 

!           Phase function moments (Rayleigh only)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (tau_l == 0.0) then
  
            greekmat_total_input(0,i) = 1.0
            greekmat_total_input(1,i) = 0.0
            greekmat_total_input(2,i) = raysmom2
            do l = 3, NGREEK_MOMENTS_INPUT        
               greekmat_total_input(l,i) = 0.0
            end do 
         end if
  

!           SCALAR phase function moments (aerosol + Rayleigh)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (tau_l /= 0.0) then
            aerswt = ssa_l * tau_l / tau_scat  
            rayswt = ray_l / tau_scat    
            aersmom(0) = 1.0
            aersmom(1) = 3.0 * g_l
            aersmom(2) = 5.0 * g_l * aersmom(1) / 3.0

            greekmat_total_input(0,i) = 1.0
            greekmat_total_input(1,i) = aersmom(1) * aerswt
            greekmat_total_input(2,i) = raysmom2 * rayswt + aersmom(2) * aerswt
         
            do l = 3, NGREEK_MOMENTS_INPUT         
               factor = REAL(2*l+1) / REAL(2*l-1) 
               aersmom(l) = factor * g_l * aersmom(l-1)
               greekmat_total_input(l,i) = aersmom(l) * aerswt
            end do  
         end if ! end if tau_l /= 0.0
            
!     end layer loop
!     ---------------
      end do
  
      self%Surface%Base%VIO%LIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = deltau_vert_input
      self%Surface%Base%VIO%LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = omega_total_input
      self%Surface%Base%VIO%LIDORT_FixIn%Optical%TS_PHASMOMS_TOTAL_INPUT = greekmat_total_input
      
!      print*, 'CALL MASTER', rc
 
!     Call the MASTER driver for doing the actual calculation
!     -----------------------------------------------------------
      call LIDORT_MASTER (self%Surface%Base%VIO%LIDORT_FixIn, &
           self%Surface%Base%VIO%LIDORT_ModIn, &
           self%Surface%Base%VIO%LIDORT_Sup, &
           self%Surface%Base%VIO%LIDORT_Out)
  
      if ( self%Surface%Base%VIO%LIDORT_Out%Status%TS_STATUS_INPUTCHECK /= 0 ) then
         rc = 3
         write(*,*) 'LIDORT_MASTER STATUS_INPUTCHECK RETURNED ERROR'
         write(*,*) self%Surface%Base%VIO%LIDORT_Out%Status%TS_STATUS_INPUTCHECK
         write(*,*) self%Surface%Base%VIO%LIDORT_Out%Status%TS_CHECKMESSAGES
      end if

      if ( self%Surface%Base%VIO%LIDORT_Out%Status%TS_STATUS_CALCULATION /= 0 ) then
         write(*,*) 'LIDORT_MASTER STATUS_CALCULATION RETURNED ERROR'
         write(*,*) self%Surface%Base%VIO%LIDORT_Out%Status%TS_STATUS_CALCULATION
         write(*,*) self%Surface%Base%VIO%LIDORT_Out%Status%TS_MESSAGE
         rc = 4
      end if
      if ( rc /= 0 ) return
      
      FLUX_FACTOR = self%Surface%Base%VIO%LIDORT_FixIn%SunRays%TS_FLUX_FACTOR
      INTENSITY   = self%Surface%Base%VIO%LIDORT_Out%Main%TS_INTENSITY ! output of LIDORT_MASTER subroutine
      
!     Return TOA radiance
!     -------------------
      RADIANCE = 0.0
      REFLECTANCE = 0.0
      RADIANCE = INTENSITY(1, 1, IDR)
           
      REFLECTANCE = (pi * RADIANCE) / ( cos(self%Surface%Base%VIO%LIDORT_ModIn%MSunRays%TS_BEAM_SZAS(1)*pi/180.0) * FLUX_FACTOR ) 
      
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
    
      end subroutine LIDORT_Run

!..............................................................................

      subroutine LIDORT_LER (self, rad, refl, rc)
 
      ! First compute the spherical albedo and the transmission of the 
      ! atmosphere purely rayleigh and
      ! the Lambertian Equivalent Reflectivity R*
      !-------------------------------------------------------------------------------------

       type(LIDORT_scat), intent(inout)           :: self 
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
       real*8,allocatable,dimension(:)   :: ROT
       real*8                            :: SSA
       
       allocate (ROT(self%Surface%Base%VIO%LIDORT_FixIn%Cont%TS_NLAYERS))
      
       ! atmosphere Rayleigh and surface albedo = 0.0
       ! --------------------------------------------
       self%Surface%sfc_type = 1 !Lambertian
       self%Surface%albedo = 0.0
       call LIDORT_Run (self, rad0, reflectance, ROT, .false., rc)

       ! atmosphere Rayleigh and surface albedo = 0.1
       ! --------------------------------------------
       self%Surface%sfc_type =1
       self%Surface%albedo = 0.1
       call LIDORT_Run (self, rad1, reflectance, ROT, .false., rc)

       ! atmosphere Rayleigh and surface albedo = 0.2
       ! --------------------------------------------
       self%Surface%sfc_type = 1
       self%Surface%albedo = 0.2
       call LIDORT_Run (self, rad2, reflectance, ROT, .false., rc)

       beta =  ( rad1 - rad0 ) / ( rad2 - rad0 )
       spher_alb = (( beta * 0.2 ) - 0.1 ) / ( 0.1 * 0.2 * ( beta - 1 )) ! Spherical albedo
       trans = ( 1 - ( spher_alb * 0.1 )) * ( rad1 - rad0 ) / 0.1        ! transmission
       
!       print*, 'spher_alb', spher_alb, trans
 
       ! Compute the Reflectivity R*
       ! --------------------------
       ! Compute the radiance for an atmosphere purely Ray and surface albedo = 0.0
       !----------------------------------------------------------------------------
       self%Surface%sfc_type = 1 
       self%Surface%albedo = 0.0
       call LIDORT_Run (self, rad_ray, reflectance, ROT, .false., rc)
        

       delta_rad =  rad - rad_ray    
        
       refl = delta_rad / ( trans + ( spher_alb * delta_rad )) ! Reflectivity R*

      

       end subroutine LIDORT_LER
       
!.........................................................................

      subroutine LIDORT_AI (self,rad, refl, AI, rc)   

      ! Compute the Aerosol Index from the reflectivity calculate 
      ! from the preVOus subroutine
      ! ---------------------------------------------------------
      type(LIDORT_scat),         intent(inout)   :: self
      real*8,             intent(in)     :: rad
      real*8,             intent(in)     :: refl      !reflectivity R*
      integer,            intent(out)    :: rc
      real*8,             intent(out)    :: AI        !Aerosol Index
      real*8                             :: rad_nouv  !calcul with R* as the surface albedo param
      real*8                             :: reflectance 
      real*8                             :: alpha  
      real*8                             :: SSA
      real*8,allocatable,dimension(:)    :: ROT
      real*8                             :: g_tot

      allocate (ROT(self%Surface%Base%VIO%LIDORT_FixIn%Cont%TS_NLAYERS))
      ! Compute the radiance for an atmosphere Ray with the R* as the surface alb parameter
      !-----------------------------------------------------------------------------------
      self%Surface%sfc_type = 1 !Lambertian
      self%Surface%albedo = refl  
      
      call LIDORT_Run (self, rad_nouv, reflectance, ROT, .false., rc)
      
      if (rad > 0.0) then
      alpha = rad / rad_nouv
      AI = -100 * log10(alpha) ! Calcul of the Aerosol Index
      else
      AI = self%MISSING
      end if
     
      
      end subroutine LIDORT_AI
!.........................................................................
      end module LIDORT_ScatMod
