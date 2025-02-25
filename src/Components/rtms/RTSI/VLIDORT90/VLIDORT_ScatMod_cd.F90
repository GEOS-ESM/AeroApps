      module VLIDORT_ScatMod_cd
 
      USE VLIDORT_PARS
      
      USE VBRDF_SUP_MOD
      USE VLIDORT_IO_DEFS
      
      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
      USE VLIDORT_MASTERS

      USE VLIDORT_Mod
      USE VLIDORT_SurfaceMod

      implicit NONE

      PUBLIC  VLIDORT_Run      ! Run for each profile (pixel)
      PUBLIC  VLIDORT_LER      ! Lambertian Equivalent Reflectivity
      PUBLIC  VLIDORT_AI       ! Aerosol Index
          

      type VLIDORT_scat

!         logical         :: initialized = .false.
         real*8          :: wavelength          ! in [nm]
!         real*8          :: albedo              ! surface albedo
!         real*8          :: solar_zenith        ! in degrees
!         real*8          :: relative_azymuth    ! in degrees
!         real*8          :: sensor_zenith       ! in degrees
         integer         :: nMom                ! number of momemts read (phase function aerosol) 
         integer         :: nPol                ! number of components of the scattering matrix
         integer         :: nMom_cd             ! number of momemts phase function cloud
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: cod(:)              ! cloud optical depth
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer :: ssac(:)             ! cloud ssa
         real*8, pointer ::   g(:)              ! aerosol asymmetry factor
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]
         real*8, pointer :: pmom(:,:,:)         ! components of the scattering phase matrix
         real*8, pointer :: cdsmom(:)         ! cloud phase function P11

         type(VLIDORT_Surface) :: Surface

       end type VLIDORT_scat
       
      Contains
!.............................................................................
 
      subroutine VLIDORT_Run (self,radiance, reflectance, AOT, rc,Q,U, & 
                                 scalar,aerosol)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE VLIDORT_PARS
         
      USE VLIDORT_IO_DEFS
      USE VBRDF_SUP_MOD
      
      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
      USE VLIDORT_MASTERS
     

      type(VLIDORT_scat), intent(inout)   :: self        ! Contains most input
      real*8,             intent(out)     :: radiance    ! TOA radiance
      real*8,             intent(out)     :: reflectance ! TOA reflectance
      real*8,             intent(out)     :: AOT         ! Aerosol Optical Thickness
      real*8,             intent(out)     :: U           ! U Stokes component
      real*8,             intent(out)     :: Q           ! Q Stokes component
!      real*8,             intent(out)     :: SSA         ! Single Scattering Albedo
      integer,            intent(out)     :: rc

      logical,            intent(in)      :: scalar  ! If True, do scalar calculation
      logical,            intent(in)      :: aerosol ! if False, Rayleigh only

!                           ----

      integer :: STATUS_INPUTCHECK, STATUS_CALCULATION 
      logical :: vector

!                           ----

!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: ncoeffs, k1, IDR, ierror
      integer                                            :: NLAYERS
      integer                                            :: NSTOKES
      real*8                                             :: ray_l      
      real*8                                             :: tau_l 
      real*8                                             :: cod_l 
      real*8                                             :: ssa_l 
      real*8                                             :: ssac_l
      real*8                                             :: g_l 
      real*8                                             :: tau_ext
      real*8                                             :: tau_scat
      real*8                                             :: ssa_tot
      real*8                                             :: raysmom2
      real*8                                             :: gammamom2
      real*8                                             :: alphamom2 
      real*8                                             :: deltamom1 
      real*8                                             :: cdswt
      real*8                                             :: aerswt 
      real*8                                             :: rayswt
      real*8                                             :: factor      
      real*8                                             :: COEF_DEPOL
      real*8                                             :: wmicron

   
      real*8                                             :: x, V
      real*8, dimension(MAXLAYERS)                       :: Ray
      real*8, dimension(0:MAXMOMENTS_INPUT)              :: aersmom 
      real*8, dimension(0:MAXMOMENTS_INPUT)              :: cdsmom 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms 
      real*8, dimension(0:2, 16)                         :: rayvmoms
      real*8                                             :: difz
      real*8                                             :: somray
      logical                                            :: DO_LAMBERTIAN_SURFACE
      real*8                                             :: LAMBERTIAN_ALBEDO
      integer                                            :: NGREEK_MOMENTS_INPUT
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: greekmat_total_input      
      real*8, dimension(0:MAXLAYERS)                     :: height_grid                     
      real*8, dimension(0:MAXLAYERS)                     :: pressure_grid 
      real*8, dimension(0:MAXLAYERS)                     :: temperature_grid

      real*8, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, dimension(MAXLAYERS)                       :: omega_total_input

      real*8, dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS)                       :: STOKES
      real*8                                             :: FLUX_FACTOR

!     Volume Rayleigh scattering coefficient (depends closely on the 
!     thermodynamic conditions of the atm and varies with height, 
!     product of the molecular number density of air by the total ray 
!     scattering cross section).
      real*8, dimension(MAXLAYERS+1)                     :: Vol 
      real*8, dimension(MAXLAYERS+1)                     :: sect 
      real*8, parameter                                  :: pi = 3.14159
      real*8, parameter                                  :: DEPOL_RATIO = 0.030
       
      
       rc = 0
      

      if ( .not. self%Surface%Base%initialized ) then
        rc = 1
        return
      end if
      
      if ( scalar ==.true. ) then
        vector = .not. scalar
      else
        vector = .true.
      end if

!                     Stokes/streams/layers/moments
!                     -----------------------------
      if ( vector ) then
         NSTOKES = 3                                       ! Number of Stokes vector components
      else
         NSTOKES = 1
        
      end if

      if ( NSTOKES  .GT. MAXSTOKES  )   rc = 2  

      self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NSTOKES          = NSTOKES
      
      if ( self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_UPWELLING ) then
         IDR = 1
      else
         IDR = 2
      end if
     
     
!                          Lambertian OR BRDF surface
!                          --------------------------
!      print*, 'SURFACE TYPE', self%Surface%sfc_type
      if ( self%Surface%sfc_type == 1 ) then               ! Lambertian surface ?
         DO_LAMBERTIAN_SURFACE = .true.        
         LAMBERTIAN_ALBEDO = self%Surface%albedo           ! Lambertian (isotropic) input albedo
      else                                   
         DO_LAMBERTIAN_SURFACE = .false. 
         LAMBERTIAN_ALBEDO = 0.0                           ! Use BRDF -> albedo set to 0
 
      end if

      self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = LAMBERTIAN_ALBEDO
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_BRDF_F_0        = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_BRDF_F_0
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_BRDF_F          = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_BRDF_F
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_BRDF_F_0
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_BRDF_F
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_EMISSIVITY      = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_EMISSIVITY
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_EMISSIVITY

!                         Angles (SZA, viewing, relatuve azimuth), Level
!                        ------------------------------------------------      
      self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1) = self%Surface%solar_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1) = self%Surface%relat_azimuth
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1) = self%Surface%sensor_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = 0.0     ! This is TOA              
      
      NGREEK_MOMENTS_INPUT = self%Surface%Base%VIO%VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
      NLAYERS = self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NLAYERS     
     
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
            
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_height_grid = height_grid
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_pressure_grid = pressure_grid
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_temperature_grid = temperature_grid



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

!      print*, 'somray', somray
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
!      SSA = 0.0
      do i = 1, NLAYERS  
         ray_l = ray(i)         ! indice l for  each layer
         cod_l = self%cod(i)
         ssac_l = self%ssac(i)
         if ( .not. aerosol ) then
            tau_l = 0.0 
            ssa_l = 0.0 
            g_l = 0.0 
         else
            tau_l = self%tau(i)            
            ssa_l = self%ssa(i)             
            g_l = self%g(i) 
            AOT = AOT + tau_l
!            SSA = SSA + ssa_l*tau_l
         end if

!        total optical depths for extinction and scattering 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tau_ext = ray_l + tau_l + cod_l
         tau_scat = ray_l +  ssa_l * tau_l + ssac_l * cod_l

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
!      print *,'ok scalar'
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if ((tau_l == 0.0).and.(cod_l == 0.0)) then
     
            greekmat_total_input(0,i,1) = 1.0
            greekmat_total_input(1,i,1) = 0.0
            greekmat_total_input(2,i,1) = raysmom2
            do l = 3, NGREEK_MOMENTS_INPUT        
               greekmat_total_input(l,i,1) = 0.0
            end do 
         end if
     

!     SCALAR phase function moments (aerosol + Rayleigh + cloud)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if ((tau_l /= 0.0).OR.(cod_l /= 0.0)) then
             cdswt = ssac_l * cod_l / tau_scat 
             aerswt = ssa_l * tau_l / tau_scat  
             rayswt = ray_l / tau_scat    
             aersmom(0) = 1.0
             aersmom(1) = 3.0 * g_l
             aersmom(2) = 5.0 * g_l * aersmom(1) / 3.0

          ! if phase fct cloud has less/or no moments 
          ! than NGREEK_MOMENTS_INPUT-> complete with 0
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             cdsmom = 0
             
             do l= 0,self%nMom_cd
                cdsmom(l) = self%cdsmom(l+1) * REAL(2*l+1)
             end do    


             greekmat_total_input(0,i,1) = 1.0
             greekmat_total_input(1,i,1) = aersmom(1) * aerswt + cdsmom(1)* cdswt
             greekmat_total_input(2,i,1) = raysmom2 * rayswt + aersmom(2) * aerswt + cdsmom(2)* cdswt
         
     
             do l = 3, NGREEK_MOMENTS_INPUT         
                factor = REAL(2*l+1) / REAL(2*l-1) 
                aersmom(l) = factor * g_l * aersmom(l-1)
                greekmat_total_input(l,i,1) = aersmom(l) * aerswt + cdsmom(l)* cdswt
             end do  
          end if ! end if tau_l /= 0.0
      end if  ! end scalar
      
!     VECTOR phase function moments 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if ( vector ) then  
!      print *, 'OK VECTOR'
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (tau_l == 0.0) then
!           print *, 'tau_l==0'
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
!          print*, 'ok taul /=0'
        
         aerswt = ssa_l * tau_l / tau_scat  
         rayswt = ray_l / tau_scat 
       do k = 1, 16
         do l= 0, NGREEK_MOMENTS_INPUT
            aervmoms(l,i,k) = 0.0
         end do
       end do

   
!         if ( self%nmom -1 < NGREEK_MOMENTS_INPUT ) then 
            do l= 0, self%nmom-1               
               aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11 
               aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12                
               aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12 = P21                            
               aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33 
               aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34            
               aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34
           
               if (self%nPol == 4) then
               aervmoms(l,i,6)  = self%pmom(i,l+1,1) ! P22 = P11 for spherical               
               aervmoms(l,i,16) = self%pmom(i,l+1,3) ! P44 = P33         
               else  if (self%nPol == 6) then   
               aervmoms(l,i,6)  = self%pmom(i,l+1,5) ! P22  for non spherical               
               aervmoms(l,i,16) = self%pmom(i,l+1,6) ! P44  
               else
               rc = 1
               end if
            end do
         

         do k = 1, 16
            do l = 0,2
            greekmat_total_input(l,i,k) = rayvmoms(l,k) * rayswt+ aervmoms(l,i,k) * aerswt
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
  
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = deltau_vert_input
      self%Surface%Base%VIO%VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = omega_total_input
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = greekmat_total_input
     
 
!     Call the MASTER driver for doing the actual calculation
!     -----------------------------------------------------------
      call VLIDORT_MASTER (self%Surface%Base%VIO%VLIDORT_FixIn, &
          self%Surface%Base%VIO%VLIDORT_ModIn, &
          self%Surface%Base%VIO%VLIDORT_Sup, &
          self%Surface%Base%VIO%VLIDORT_Out)

      if ( self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_INPUTCHECK /= 0 ) rc = 3
      if ( self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_CALCULATION /= 0 ) rc = 4
      if ( rc /= 0 ) return
      
      STOKES = self%Surface%Base%VIO%VLIDORT_Out%Main%TS_STOKES ! output of VLIDORT_MASTER subroutine
      FLUX_FACTOR = self%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR

!     Return TOA radiance
!     -------------------
      RADIANCE = 0.0
      REFLECTANCE = 0.0
      Q = 0
      U = 0
      V = 0
      if ( scalar ) then
       RADIANCE = STOKES(1, 1, 1, IDR)
              
       REFLECTANCE = (pi * RADIANCE) / ( cos(self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1)*pi/180.0) * FLUX_FACTOR ) 
      
      end if
  
      if ( vector ) then
       RADIANCE = STOKES(1, 1, 1, IDR)
       Q = STOKES(1, 1, 2, IDR)
       U = STOKES(1, 1, 3, IDR)
!       V = STOKES(1, 1, 4, IDR)
!      print *, 'RAD', RADIANCE, Q, U
       REFLECTANCE = (pi * RADIANCE) / ( cos(self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1)*pi/180.0) * FLUX_FACTOR )   
      end if
!     print *, 'RADIANCE',self%wavelength, RADIANCE

!      STOP

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
    
      end subroutine VLIDORT_Run

!..............................................................................

      subroutine VLIDORT_LER (self, rad, refl, rc)
 
      ! First compute the spherical albedo and the transmission of the 
      ! atmosphere purely rayleigh and
      ! the Lambertian Equivalent Reflectivity R*
      !-------------------------------------------------------------------------------------

       type(VLIDORT_scat), intent(inout)           :: self 
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
       real*8                            :: SSA
       real*8                            :: U, Q
       
      
       ! atmosphere Rayleigh and surface albedo = 0.0
       ! --------------------------------------------
       self%Surface%sfc_type = 1 !Lambertian
       self%Surface%albedo = 0.0
       call VLIDORT_Run (self, rad0, reflectance, AOT,  rc,Q,U, .false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.1
       ! --------------------------------------------
       self%Surface%sfc_type =1
       self%Surface%albedo = 0.1
       call VLIDORT_Run (self, rad1, reflectance, AOT, rc, Q,U,.false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.2
       ! --------------------------------------------
       self%Surface%sfc_type = 1
       self%Surface%albedo = 0.2
       call VLIDORT_Run (self, rad2, reflectance, AOT, rc,Q,U, .false., .false.)

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
       call VLIDORT_Run (self, rad_ray, reflectance, AOT, rc, Q,U,.false., .false.)
       print*, 'rad_RAYLEIGH', rad_ray      

       delta_rad =  rad - rad_ray    
        
       refl = delta_rad / ( trans + ( spher_alb * delta_rad )) ! Reflectivity R*

       print*, 'refl', refl

       end subroutine VLIDORT_LER
       
!.........................................................................

      subroutine VLIDORT_AI (self,rad, refl, AI, rc)   

      ! Compute the Aerosol Index from the reflectivity calculate 
      ! from the preVOus subroutine
      ! ---------------------------------------------------------
       type(VLIDORT_scat),         intent(inout)   :: self
       real*8,             intent(in)     :: rad
       real*8,             intent(in)     :: refl      !reflectivity R*
       integer,            intent(out)    :: rc
       real*8,             intent(out)    :: AI        !Aerosol Index
       real*8                             :: rad_nouv  !calcul with R* as the surface albedo param
       real*8                             :: reflectance 
       real*8                             :: alpha  
       real*8                             :: SSA
       real*8                             :: AOT 
       real*8                             :: g_tot
       real*8                             :: Q, U
      ! Compute the radiance for an atmosphere Ray with the R* as the surface alb parameter
      !-----------------------------------------------------------------------------------
      self%Surface%sfc_type = 1 !Lambertian
      self%Surface%albedo = refl  
      
      call VLIDORT_Run (self, rad_nouv, reflectance, AOT, rc,Q,U, .false., .false.)
      print*, 'calcul AI', refl
      if (rad > 0.0) then
      alpha = rad / rad_nouv
      AI = -100 * log10(alpha) ! Calcul of the Aerosol Index
      else
      AI = self%MISSING
      end if
     
      
      end subroutine VLIDORT_AI
!.........................................................................
      end module VLIDORT_ScatMod_cd
