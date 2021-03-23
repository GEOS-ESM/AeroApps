
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p2                         #
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
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
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
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! #            VFO_MASTER_INTERFACE                             #
! #                                                             #
! ###############################################################

!  this is the standard interface to the FO code - 

!  Version 2.0 - 2.7. Internal FOCORR/DBCORR routines (Old)
!  Version 2.7. VFO interface to FO code version 1.4 was written
!  Version 2.8. VFO interface to FO code version 1.5 (upgrade)
!  Version 2.8. Internal FOCORR/DBCORR routines removed.

!  Version 2.8.1. 4/9/19.
!     -- Need to add the FO Surface-leaving assignation and saved cumulative transmittance

!  4/15/20. Version 2.8.2.
!      -- Argument list considerably reduced in size
!      -- No geometrical variables (all removed), now precalculated in FOGeometry (input)
!      -- No greekmoms_total_input setting, no DO_FMATRIX and no ngreek_moments_input
!      -- ssflux and na_offset are input directly (not re-set here)
!      -- Exception handling no longer required
!      -- Add DO_DOUBLET_GEOMETRY flag and Offset array (ND_OFFSET)

      MODULE vlidort_vfo_interface_m

!  User dependencies

      USE VLIDORT_Setups_def_m
      USE VFO_RTCalc_Master_m

      PRIVATE
      PUBLIC :: VFO_MASTER_INTERFACE

      CONTAINS

!

      SUBROUTINE VFO_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
        DO_PLANPAR, DO_ENHANCED_PS, DO_DELTAM_SCALING,                                      & ! Input FO control flags
        DO_UPWELLING, DO_DNWELLING, DO_OBSGEOM, DO_DOUBLET, DO_PARTLAYERS,                  & ! input Model control flags
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,       & ! Input Surface flags
        NSTOKES, NLAYERS, NGEOMS, NSZAS, NVZAS, NAZMS, NA_OFFSET, ND_OFFSET,                & ! Input Numbers/Offsets
        N_USER_LEVELS, USER_LEVEL_MASK_UP, USER_LEVEL_MASK_DN, N_PARTLAYERS,                & ! Input levels  control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, FOGeometry,           & ! Input partial/Geometry
        HEIGHT_GRID, SSFLUX, FLUXVEC, DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                 & ! Inputs (Optical)
        DELTAU_VERT, FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR,                                  & ! Inputs (Optical)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SLTERM_ISOTROPIC, SLTERM_USERANGLES,            & ! Inputs (Surface)
        BB_INPUT, SURFBB, USER_EMISSIVITY,                                                  & ! Inputs (Thermal)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output Stokes vectors
        FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES, CUMTRANS, SLTERM )                        ! Output Stokes vectors

!  4/15/20. Version 2.8.2. Drop MAXMOMENTS, MAXMOMENTS_INPUT, MAXFINELAYERS
!  4/15/20. Version 2.8.2. Argument list follows that for LIDORT 3.8.2
!  4/15/20. Version 2.8.2. Add DO_DOUBLET_GEOMETRY flag and Offset array (ND_OFFSET)

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                 MAXLAYERS, MAXSTOKES, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,     &
                                 MAX_GEOMETRIES, MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO

      IMPLICIT NONE

!  Parameter argument
!  ------------------

      INTEGER, PARAMETER  :: FFP = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

!  Flags. Fmatrix flag (7/7/16), Surface leaving inputs (8/3/16), Introduced for Version 2.8
!     ==> Water-leaving flag introduced for Version 2.8.1. 4/9/19.
!     ==> 4/15/20. Version 2.8.2. Add Do_Doublet flag

      LOGICAL, INTENT(IN) ::            DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::            DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_EMISSION

      LOGICAL, INTENT(IN) ::            DO_PLANPAR
      LOGICAL, INTENT(IN) ::            DO_ENHANCED_PS

      LOGICAL, INTENT(IN) ::            DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING

      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_OBSGEOM
      LOGICAL, INTENT(IN) ::            DO_DOUBLET

      LOGICAL, INTENT(IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::            DO_WATER_LEAVING
      LOGICAL, INTENT(IN) ::            DO_SL_ISOTROPIC

!  Control integers
!    ==> 4/15/20. Version 2.8.2. Add Do_Doublet Offset

      INTEGER  , INTENT(IN) :: NSTOKES, NLAYERS, NGEOMS, NSZAS, NVZAS, NAZMS
      INTEGER  , INTENT(IN) :: NA_OFFSET ( MAXBEAMS, MAX_USER_STREAMS )
      INTEGER  , INTENT(IN) :: ND_OFFSET ( MAXBEAMS )

!  Type structure, FOGeometry input
  
      TYPE(VLIDORT_Geometry_FO), INTENT(IN)      :: FOGeometry

!  require the Level Mask inputs

      INTEGER, INTENT (IN) :: N_USER_LEVELS
      INTEGER, INTENT (IN) :: USER_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: USER_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
!      Real(fpk), INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs.

      LOGICAL  , INTENT(IN) :: DO_PARTLAYERS
      INTEGER  , INTENT(IN) :: N_PARTLAYERS
      LOGICAL  , INTENT(IN) :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER  , INTENT(IN) :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER  , INTENT(IN) :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Height grid

      DOUBLE PRECISION, INTENT(IN) ::   HEIGHT_GRID ( 0:MAXLAYERS )

!  other inputs

      DOUBLE PRECISION, INTENT(IN) ::   SSFLUX
      DOUBLE PRECISION, INTENT(IN) ::   FLUXVEC ( MAXSTOKES )

!  Optical unscaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  optical scaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   TRUNC_FACTOR ( MAXLAYERS )

!  F-matrices, introduced for Version 2.8 of VLIDORT, and Version 1.5 of FO
!mick fix 9/19/2017 - swapped layer & geo dimensions

      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Thermal

      DOUBLE PRECISION, INTENT(IN) ::    BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::    SURFBB
      DOUBLE PRECISION, INTENT(IN) ::    USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  surface BRDF

      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

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

!  Local variables
!  ===============

!  Vector sunlight flag

      LOGICAL   :: DO_SUNLIGHT

!  Vector emissivity flag (12/11/17 Rob add.)
!      4/15/20. Version 2.8.2. Introduced as an input to the FO call

      LOGICAL   :: do_Polarized_Emissivity

!  Atmosphere

      REAL(FFP) :: EXTINCTION  ( MAXLAYERS )
      REAL(FFP) :: DELTAUS     ( MAXLAYERS )
      REAL(FFP) :: OMEGA       ( MAXLAYERS )
      REAL(FFP) :: TRUNCFAC    ( MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      REAL(FFP) :: REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES )

!  Other variables

      integer   :: ns, nv, na
      INTEGER   :: G, GK, O1, O2
      INTEGER   :: LUM=1, LUA=1

!  Define VFO_MASTER inputs
!  ========================

!  4/15/20. Version 2.8.2.
!      -- No geometrical variables (all removed), now precalculated in FOGeometry
!      -- No phasmoms setting, no DO_FMATRIX and no ngreek_moments_input
!      -- ssflux and na_offset are input directly (not re-set here)
!      -- Exception handling no longer required
!      -- Add Doublet geometry control

!  Set sunlight flag

      do_sunlight        = .TRUE. !set for now

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

! For TMS correction only

      truncfac = zero
      if (do_deltam_scaling) then
         truncfac(1:nlayers) = TRUNC_FACTOR(1:nlayers)
      end if

!  BRDF reflection now also for Lattice.
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet option for BRDF reflection.

      reflec = ZERO
      if (do_lambertian_surface) then
        reflec(1,1,1:ngeoms) = LAMBERTIAN_ALBEDO
      else
        if ( do_ObsGeom ) then
          do o1=1,nstokes ; do o2=1,nstokes
            gk = 4*(o1-1) + o2
            reflec(o1,o2,1:ngeoms) = EXACTDB_BRDFUNC(gk,lum,lua,1:ngeoms)
          end do ; end do
        else  if ( do_doublet ) then
          do ns = 1, nszas ; do nv = 1, nvzas
            g = nd_offset(ns) + nv
            do o1=1,nstokes ; do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,nv,LUA,ns)
            end do ; end do
          end do ; end do
        else
          do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
            g = na_offset(ns,nv) + na
            do o1=1,nstokes ; do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,nv,na,ns)
            end do ; end do
          enddo ; enddo ; enddo
        end if
      end if

!  Surface leaving, introduced 8/3/16. Proper Initializing, 8/3/16
!   ==> 4/9/19. This is now an output.
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet option for Surface leaving.

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
        else if ( do_Doublet ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(1,g) = SLTERM_ISOTROPIC(1,ns)
            end do ; end do 
          ELSE
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               do o1 = 1,nstokes
                  SLTERM(O1,g) = SLTERM_USERANGLES(O1,nv,LUA,ns)
               end do
            end do ; end do
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

!  Set do_Polarized_Emissivity flag. 12/11/17 Rob add.
!    Now set outside this module , using the following code
!mick fix 3/2/2020 - trimmed summing of USER_EMISSIVITY elements

      do_Polarized_Emissivity = .false.
      if ( nstokes.gt.1 .and. do_surface_emission ) then
         do nv = 1, nvzas
           !if ( SUM(USER_EMISSIVITY(2:4,nv)).ne.zero ) do_Polarized_Emissivity = .true.
           if ( SUM(USER_EMISSIVITY(2:nstokes,nv)).ne.zero ) do_Polarized_Emissivity = .true.
         enddo
      endif

!  Call VFO_RTCALC_MASTER
!  ======================

!  Upgraded to FO Version 1.5, 7/7/16. Surface Leaving 8/3/16. Partials 9/17/16.
!   4/9/19. CUMTRANS now an additional output. add water-leaving control.

!  4/15/20. Version 2.8.2.
!      -- Name of subroutine changed. Input arguments rearranged somewhat
!      -- FOGeometry structure passed instead of many geometry variables
!      -- Removed: greekmat, do_fmatrix, ngreek_moments_input.  flux --> ssflux
!      -- Dimensioning removed; now not so stand-alone!!!!
!      -- Exception handling removed
!      -- Add doublet geometry flag and offset to input list

      CALL VFO_RTCALC_MASTER &
       ( do_solar_sources, do_thermal_emission, do_surface_emission, do_Polarized_Emissivity,      & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_lambertian_surface, do_surface_leaving, do_water_leaving,  & ! Input flags (Up/dn/surface)
         do_sunlight, do_deltam_scaling, do_obsgeom, do_doublet, do_Partlayers, do_planpar,        & ! Input flags (RT Control)
         do_enhanced_ps, nstokes, ngeoms, nszas, nvzas, nazms, na_offset, nd_offset, nlayers,      & ! Inputs (Numbers/offsets)
         n_user_levels, user_level_mask_up, user_level_mask_dn,                                    & ! Inputs (control-levels)
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx, FOGeometry,   & ! Inputs (control-partial)
         ssflux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,            & ! Input atmos optical
         bb_input, surfbb, user_emissivity, reflec, slterm,                                        & ! Input thermal/surface/Geometry
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                                 & ! Output
         fo_stokes_atmos, fo_stokes_surf, fo_stokes, cumtrans )                                      ! Output

!  Done

      RETURN
END SUBROUTINE VFO_MASTER_INTERFACE

!  finish module

END MODULE vlidort_vfo_interface_m


