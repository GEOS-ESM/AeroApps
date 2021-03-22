
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p1                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr                          #
! #                Matt Christi                                 #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p1                             #
! #  Release Date :   31 August 2019                            #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
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
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.1 of the VLIDORT software library.          #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2019.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of VLIDORT_2p8p1 ( Version 2.8.1 )            #
! #                                                                 #
! # VLIDORT_2p8p1 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p1 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p1  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ############################################################### Old
! # Subroutines in the CORRECTIONS Module                       # Old
! #                                                             # Old
! #      Version 2.0. Nadir single scatter correction           # Old
! #            VLIDORT_SSCORR_NADIR (master)                    # Old
! #              VLIDORTSS_FMATRICES                            # Old
! #              VLIDORTSS_FMATRICES_MULTI                      # Old
! #                                                             # Old
! #      Version 2.2. outgoing sphericity correction            # Old
! #              2.3. partial-layer integration                 # Old
! #            VLIDORT_SSCORR_OUTGOING (master)                 # Old
! #                 SSCORR_OUTGOING_ZMATRIX                     # Old
! #                 OUTGOING_INTEGRATION_UP                     # Old
! #                 OUTGOING_INTEGRATION_DN                     # Old
! #                 VLIDORT_DBCORRECTION                        # Old
! #                                                             # Old
! #      Version 2.4RTC ------- Notes -------                   # Old
! #                                                             # Old
! #            Additional correction to ZMATRIX and FMATRIX     # Old
! #            routines for Azimuth > 180. HvdM, Eqs(94/95)     # Old
! #            Implemented by V. Natraj and R. Spurr, 5/1/09    # Old
! #            Correctly done, October 2010.                    # Old
! #            Correction for Downwelling SSCORR_OUTGOING       # Old
! #            Separate Geometry routine implemented.           # Old
! #            Implemented by R. Spurr, 10/06/10                # Old
! #            Consolidated DB correction (Lambertian/BRDF)     # Old
! #                                                             # Old
! ###############################################################

! ###############################################################
! #                                                             #
! #            VFO_MASTER_INTERFACE                             #
! #                                                             #
! ###############################################################

!  this is the standard interface to the FO code - 

!  Version 2.0 - 2.7. Internal SSCORR/DBCORR routines (Old)
!  Version 2.7. VFO interface to FO code version 1.4 was written
!  Version 2.8. VFO interface to FO code version 1.5 (upgrade)
!  Version 2.8. Internal SSCORR/DBCORR routines removed.

!  Version 2.8.1. 4/9/19.
!     -- Need to add the FO Surface-leaving assignation and saved cumulative transmittance

MODULE vlidort_vfo_interface_m

      USE VFO_Master_m

      PRIVATE
      PUBLIC :: VFO_MASTER_INTERFACE

      CONTAINS

!

      SUBROUTINE VFO_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_DELTAM_SCALING,          & ! Input SS control flags
        DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_PARTLAYERS, DO_FMATRIX,     & ! input Model control flags
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,       & ! Input Surface flags
        NSTOKES, NLAYERS, NFINELAYERS, NGREEK_MOMENTS_INPUT,                                & ! Input numbers
        N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES, N_USER_RELAZMS, USER_RELAZMS, & ! Input geometry
        N_USER_LEVELS, USER_LEVEL_MASK_UP, USER_LEVEL_MASK_DN, N_PARTLAYERS,                & ! Input levels  control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
        EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER, FLUXVEC,                             & ! Input Flux/Heights
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
        DELTAU_VERT, FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR, BB_INPUT,                        & ! Inputs (Optical - Regular)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                        & ! Inputs (Optical - Surface)
        SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output Stokes vectors
        FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES, CUMTRANS, SLTERM,                        & ! Output Stokes vectors
        FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Exception handling

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                 MAXLAYERS, MAXSTOKES, MAXMOMENTS_INPUT, MAXSTOKES_SQ, MAXMOMENTS,   &
                                 MAX_USER_STREAMS, MAXBEAMS, MAX_GEOMETRIES, MAXFINELAYERS,          &
                                 MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO, PIE, DEG_TO_RAD

      IMPLICIT NONE

!  Parameter argument
!  ------------------

      INTEGER, PARAMETER  :: FFP = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

!  Flags. Fmatrix flag (7/7/16), Surface leaving inputs (8/3/16), Introduced for Version 2.8
!     Water-leaving flag introduced for Version 2.8.1. 4/9/19.
      
      LOGICAL, INTENT(IN) ::            DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::            DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_EMISSION

      LOGICAL, INTENT(IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::            DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) ::            DO_SSCORR_OUTGOING

      LOGICAL, INTENT(IN) ::            DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING
      LOGICAL, INTENT(IN) ::            DO_FMATRIX

      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_OBSERVATION_GEOMETRY

      LOGICAL, INTENT(IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::            DO_WATER_LEAVING
      LOGICAL, INTENT(IN) ::            DO_SL_ISOTROPIC

!  Control integers

      INTEGER, INTENT(IN) ::            NSTOKES
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS
      INTEGER, INTENT(IN) ::            NGREEK_MOMENTS_INPUT

!  Geometry

      INTEGER, INTENT(IN) ::            N_SZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_VZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_RELAZMS
      DOUBLE PRECISION, INTENT(IN) ::   USER_RELAZMS  ( MAX_USER_RELAZMS )

!  Levels. Now require the Level Mask inputs, 9/17/16.

      INTEGER, INTENT(IN) ::            N_USER_LEVELS
      INTEGER, INTENT (IN) ::           USER_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           USER_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!      DOUBLE PRECISION, INTENT(IN) ::   USER_LEVELS ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs, added 9/17/16

      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  other inputs

      DOUBLE PRECISION, INTENT(IN) ::   EARTH_RADIUS
      DOUBLE PRECISION, INTENT(IN) ::   HEIGHT_GRID ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   SS_FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT(IN) ::   FLUXVEC ( MAXSTOKES )

!  Optical unscaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   BB_INPUT ( 0:MAXLAYERS )

!  optical scaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   TRUNC_FACTOR ( MAXLAYERS )

!  F-matrices, introduced for Version 2.8 of VLIDORT, and Version 1.5 of FO
!mick fix 9/19/2017 - swapped layer & geo dimensions

      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  surface BRDF and EMiss.

      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT(IN) ::   SURFBB
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Surface leaving inputs, 8/3/16, for Version 2.8

      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Outputs
!  =======

!  Solar

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite
!mick mod 9/19/2017 - added FO_STOKES_ATMOS & FO_STOKES_SURF to vector FO output

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  4/9/19. Additional output for the sleave correction

      DOUBLE PRECISION, Intent(out)     :: CUMTRANS ( max_user_levels, max_geometries )

!  4/9/19. Surface leaving FO assignation

      DOUBLE PRECISION, intent(out)     :: SLTERM ( MAXSTOKES, MAX_GEOMETRIES )

!  Exception handling

      LOGICAL, INTENT (OUT)             :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1, TRACE_2

!  FO VARIABLES TAKEN DIRECTLY FROM THE VLIDORT INPUT
!  ==================================================

!  Dimensioning integers

       !INTEGER   :  MAXLAYERS, MAXMOMENTS_INPUT, MAX_USER_LEVELS
       !INTEGER   :: N_USER_LEVELS, NSTOKES, N_USER_LEVELS

!  Sources control, including thermal, deltam_scaling, directional flags

      !LOGICAL   :: DO_SOLAR_SOURCES
      !LOGICAL   :: DO_THERMAL_EMISSION
      !LOGICAL   :: DO_SURFACE_EMISSION
      !LOGICAL   :: DO_DELTAM_SCALING
      !LOGICAL   :: DO_UPWELLING, DO_DNWELLING

!  FluxVector

      !REAL(FFP) :: FLUXVEC ( MAXSTOKES )

!  surface property

      !REAL(FFP) :: SURFBB

!  Local variables
!  ===============

!  Max dimensions
!  --------------

      INTEGER   :: MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, MAXFINE

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

      INTEGER   :: NGEOMS, NSZAS, NVZAS, NAZMS, NFINE, NMOMENTS_INPUT

!  Configuration inputs
!  --------------------

!  Flags (sphericity flags should be mutually exclusive)

      LOGICAL   :: DO_PLANPAR
      LOGICAL   :: DO_ENHANCED_PS

!  Lambertian surface flag

      LOGICAL   :: DO_LAMBERTIAN

!  Vector sunlight flag

      LOGICAL   :: DO_SUNLIGHT

!  Observational geometry flag

      LOGICAl   :: do_ObsGeom

!  General inputs
!  --------------

!  DTR = degrees-to-Radians

      REAL(FFP) :: DTR

!  Critical adjustment for cloud layers

      LOGICAL   :: DoCrit
      REAL(FFP) :: Acrit

!  Earth radius + heights. Partial heights added 9/17/16

      REAL(FFP) :: ERADIUS
      REAL(FFP) :: HEIGHTS ( 0:MAXLAYERS )
      REAL(FFP) :: PARTIAL_HEIGHTS ( MAX_PARTLAYERS )

!  Output levels. Not required now, 9/17/16.
!      INTEGER   :: FO_USER_LEVELS ( MAX_USER_LEVELS )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees). Note the dimensioning now,

      REAL(FFP) :: obsgeom_boa  ( MAX_GEOMETRIES, 3 )
      REAL(FFP) :: theta_boa ( MAX_SZANGLES )        !SZA
      REAL(FFP) :: alpha_boa ( MAX_USER_VZANGLES )   !UZA
      REAL(FFP) :: phi_boa   ( MAX_USER_RELAZMS )    !RAA

!  Optical inputs
!  --------------

!  Solar flux

      REAL(FFP) :: FLUX

!  Atmosphere

      REAL(FFP) :: EXTINCTION  ( MAXLAYERS )
      REAL(FFP) :: DELTAUS     ( MAXLAYERS )
      REAL(FFP) :: OMEGA       ( MAXLAYERS )
      REAL(FFP) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, MAXSTOKES, MAXSTOKES )

!  For TMS correction

      REAL(FFP) :: TRUNCFAC    ( MAXLAYERS )

!  Thermal inputs. Use directly, don't need to copy
!      REAL(FFP) :: BB_INPUT ( 0:MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      REAL(FFP) :: REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES )

!  Surface properties - emissivity
!    - polarized surface emissivity. 12/11/17 Rob Add.

!      REAL(FFP) :: EMISS ( MAX_USER_VZANGLES )
      REAL(FFP) :: EMISS ( MAXSTOKES, MAX_USER_VZANGLES )

!  Surface leaving FO term (8/3/16, Introduced).  4/9/19/ This is now an output
!      REAL(FFP) :: SLTERM ( MAXSTOKES, MAX_GEOMETRIES )

!  Other variables

      integer   :: ns, nv, na, nv_offset(max_szangles), na_offset(max_szangles,max_user_vzangles)
      INTEGER   :: G, GK, L, N, O1, O2, UT
      INTEGER   :: LUM=1, LUA=1

!  Define VFO_MASTER inputs
!  ========================

!  Note: argument passing for variables defined elsewhere commented out here

!  Max dimensions.

      maxgeoms          = MAX_GEOMETRIES
      maxfine           = MAXFINELAYERS
      maxszas           = MAX_SZANGLES
      maxvzas           = MAX_USER_VZANGLES
      maxazms           = MAX_USER_RELAZMS

!  Mode

      do_ObsGeom  = DO_OBSERVATION_GEOMETRY

!  Dimensions

      if ( do_ObsGeom ) then
         ngeoms  = N_SZANGLES
      else
         ngeoms  = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS
      endif
      nszas   = N_SZANGLES
      nvzas   = N_USER_VZANGLES
      nazms   = N_USER_RELAZMS

      nfine           = NFINELAYERS
      nmoments_input  = NGREEK_MOMENTS_INPUT

!  Offsets

      na_offset = 0 ; nv_offset = 0
      if ( .not. do_obsgeom ) then
        do ns = 1, nszas
          nv_offset(ns) = nvzas * nazms * (ns - 1) 
          do nv = 1, nvzas
            na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
          enddo
        enddo
      endif

!  Configuration inputs. Removed "regular_ps" flag 9/17/16 (Redundant)

      if (DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (DO_SSCORR_NADIR) then
          do_enhanced_ps = .FALSE.
        else if (DO_SSCORR_OUTGOING) then
          do_enhanced_ps = .TRUE.
        end if
      end if

      do_lambertian      = DO_LAMBERTIAN_SURFACE
      do_sunlight        = .TRUE. !set for now

!  General inputs

      dtr    = DEG_TO_RAD

!      DoCrit = .FALSE. !set for now
!      Acrit  = ZERO    !set for now

!  Rob Fix 3/16/15. Default values changed 
!   [ Note for future - these should be inputs ]

      DoCrit = .TRUE.     
      Acrit  = 1.0d-12

!  Earth Radius and heights. Partials added, 9/17/16

      eradius             = EARTH_RADIUS
      heights(0:nlayers)  = HEIGHT_GRID(0:NLAYERS)
      do ut = 1, n_partlayers
        n = partlayers_layeridx(ut)
        partial_heights(ut) = heights(n-1) - partlayers_values(ut) * ( heights(n-1) - heights(n) )
      enddo

!  This code no longer required, now using same arrays as VLIDORT
!      fo_user_levels(1:n_user_levels) = INT(user_levels(1:n_user_levels))

!  Geometry inputs. all angles in Degrees.
!  These are good for Lattice and ObsGeom

      obsgeom_boa = zero
      if ( do_ObsGeom) then
         obsgeom_boa(1:ngeoms,1) = SZANGLES      (1:ngeoms)
         obsgeom_boa(1:ngeoms,2) = USER_VZANGLES (1:ngeoms)
         obsgeom_boa(1:ngeoms,3) = USER_RELAZMS  (1:ngeoms)
      endif
      theta_boa(1:N_SZANGLES)       = SZANGLES(1:N_SZANGLES)
      alpha_boa(1:N_USER_VZANGLES)  = USER_VZANGLES(1:N_USER_VZANGLES)
      phi_boa  (1:N_USER_RELAZMS)   = USER_RELAZMS (1:N_USER_RELAZMS)

!  Optical inputs

      flux     = SS_FLUX_MULTIPLIER !=FLUX_FACTOR/(4.0d0*PIE)

! @@ Rob 7/30/13
!   If deltam-scaling, must use scaled deltau !!!

      extinction = 0.0d0
      if (do_deltam_scaling) then
        extinction(1:nlayers) = DELTAU_VERT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT(1:NLAYERS)
      else
        extinction(1:nlayers) = DELTAU_VERT_INPUT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT_INPUT(1:NLAYERS)
      endif

!  Single scattering albedo is unscaled

      omega(1:nlayers)      = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Version 2.8 (FO 1.5). Optional use of Fmatrix input
!    Greekmat will not be needed in this case

      greekmat = zero 
      if ( .not. do_fmatrix ) then
        do o1=1,nstokes
          do o2=1,nstokes
            gk = 4*(o1-1) + o2
            do l=0,nmoments_input
              do n=1,nlayers
                greekmat(n,l,o1,o2) = GREEKMAT_TOTAL_INPUT(L,N,GK)
              end do
            end do
          end do
        end do
      endif

! For TMS correction only

      truncfac = zero
      if (do_deltam_scaling) then
        do n=1,nlayers
          truncfac(n) = TRUNC_FACTOR(N)
        end do
      end if

!  copy thermal BB input. Now used directly.
!      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

!  BRDF reflection now also for Lattice.

      reflec = ZERO
      if (do_lambertian) then
        do g = 1,ngeoms
          reflec(1,1,g) = LAMBERTIAN_ALBEDO
        end do
      else
        if ( do_ObsGeom ) then
           do g=1,ngeoms
             do o1=1,nstokes
               do o2=1,nstokes
                 gk = 4*(o1-1) + o2
                 reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,lum,lua,g)
               end do
             end do
           end do
        else
          do ns = 1, nszas
            do nv = 1, nvzas
              do na = 1, nazms
                g = na_offset(ns,nv) + na
                do o1=1,nstokes
                  do o2=1,nstokes
                    gk = 4*(o1-1) + o2
                    reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,nv,na,ns)
                  end do
                end do
              end do
            end do
          end do
        end if
      end if

!  Surface leaving, introduced 8/3/16. Proper Initializing, 8/3/16
! 4/9/19 ==> This is now an output.

      SLTERM = zero
      if ( do_surface_leaving ) then
        if ( do_ObsGeom ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            SLTERM(1,1:ngeoms) = SLTERM_ISOTROPIC(1,1:ngeoms)
          ELSE
            DO O1 = 1, NSTOKES
              SLTERM(O1,1:ngeoms) = SLTERM_USERANGLES(O1,LUM,LUA,1:ngeoms)
            ENDDO
          ENDIF
        else
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(1,g) = SLTERM_ISOTROPIC(1,ns)
            end do ; end do ; end do 
         ELSE
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do o1 = 1,nstokes
                  SLTERM(O1,g) = SLTERM_USERANGLES(O1,nv,na,ns)
               end do
            end do ; end do ; end do 
          ENDIF
        endif
      endif

!  Emissivity. Proper Initializing, 8/3/16
!    - polarized surface emissivity. 12/11/17 Rob Add.

      emiss = zero
      if ( do_surface_emission ) then
        if ( do_ObsGeom ) then
          do g=1,ngeoms
!            emiss(g) = USER_EMISSIVITY(1,g)
            emiss(1:nstokes,g) = USER_EMISSIVITY(1:nstokes,g)
          end do
        else
          do nv=1,nvzas
!            emiss(nv) = USER_EMISSIVITY(1,nv)
            emiss(1:nstokes,nv) = USER_EMISSIVITY(1:nstokes,nv)
          end do
        end if
      endif

!  Call VFO_MASTER
!  ===============

!  Upgraded to FO Version 1.5, 7/7/16. Surface Leaving 8/3/16. Partials 9/17/16.
!   4/9/19. CUMTRANS now an additional output. add water-leaving control.

      CALL VFO_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms,                                          & ! Input max dims
         maxlayers, max_partlayers, maxfine, maxmoments_input, max_user_levels,        & ! Input max dims
         do_solar_sources, do_sunlight, do_thermal_emission, do_surface_emission,      & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_deltam_scaling, do_fmatrix, do_obsgeom,        & ! Input flags (general)
         do_lambertian, do_surface_leaving, do_water_leaving,                          & ! Input flags (surface)
         do_Partlayers, do_planpar, do_enhanced_ps,                                    & ! Input flags (geoms)
         nstokes, ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,         & ! Inputs (control-numbers)
         n_user_levels, user_level_mask_up, user_level_mask_dn,                        & ! Inputs (control-levels)
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx,   & ! Inputs (control-partial)
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                   & ! Input general
         obsgeom_boa, alpha_boa, theta_boa, phi_boa, flux, fluxvec,                    & ! Input geometry/Flux
         extinction, deltaus, omega, greekmat, fmatrix_up, fmatrix_dn,                 & ! Input atmos optical
         truncfac, bb_input, surfbb, emiss, reflec, slterm,                            & ! Input thermal/surf optical
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                     & ! Output
         fo_stokes_atmos, fo_stokes_surf, fo_stokes, cumtrans,                         & ! Output
         fail, message, trace_1, trace_2 )                                               ! Exception handling

!  Done

      RETURN
END SUBROUTINE VFO_MASTER_INTERFACE

!  finish module

END MODULE vlidort_vfo_interface_m


