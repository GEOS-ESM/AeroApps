
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
! #            VFO_LCS_MASTER_INTERFACE                         #
! #                                                             #
! ###############################################################

!  Version 2.0 - 2.7. Internal SSCORR/DBCORR routines (Old)
!  Version 2.7. VFO interface to FO code version 1.4 was written
!  Version 2.8. VFO interface to FO code version 1.5 (upgrade)
!  Version 2.8. Internal SSCORR/DBCORR routines removed.

!   -- 4/9/19. Need to add the FO Surface-leaving assignation and saved cumulative transmittance
!              Also their linearizations

!  4/15/20. Version 2.8.2.
!      -- Argument list considerably reduced in size
!      -- No geometrical variables (all removed), now precalculated in FOGeometry (input)
!      -- No phasmoms_total_input setting, no do_Phasfunc and no nmoments_input
!      -- ssflux and na_offset are input directly (not re-set here)
!      -- Exception handling no longer required
!      -- Add DO_DOUBLET_GEOMETRY flag and Offset array (ND_OFFSET)

      MODULE vlidort_vfo_lcs_interface_m

      use VLIDORT_Setups_def_m
      USE VFO_RTCalc_LinMasters_m, Only : VFO_RTCALC_LCS_MASTER

      PUBLIC :: VFO_LCS_MASTER_INTERFACE

      CONTAINS

      SUBROUTINE VFO_LCS_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                      & ! Input Sources flags
        DO_PLANPAR, DO_ENHANCED_PS, DO_DELTAM_SCALING,                                   & ! Input SS control flags
        DO_UPWELLING, DO_DNWELLING, DO_OBSGEOM, DO_DOUBLET, DO_PARTLAYERS,               & ! input Model control flags
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,    & ! Input Surface flags
        DO_COLUMNWFS, DO_SURFACEWFS, DO_SLEAVEWFS,                                       & ! Input Jacobian flags
        NSTOKES, NLAYERS, NGEOMS, NSZAS, NVZAS, NAZMS, NA_OFFSET, ND_OFFSET,             & ! Input numbers
        N_USER_LEVELS, USER_LEVEL_MASK_UP, USER_LEVEL_MASK_DN, N_PARTLAYERS,             & ! Input levels  control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                    & ! Input partial
        N_COLUMNWFS, LVARYFMAT, N_SLEAVEWFS, N_REFLECWFS, N_SURFACEWFS, FOGeometry,      & ! Input Jacobians/Geometry
        HEIGHT_GRID, SSFLUX, FLUXVEC, DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, DELTAU_VERT, & ! Inputs (Optical)
        FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR, LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC,        & ! Inputs (Optical/Surface) 
        SLTERM_ISOTROPIC, SLTERM_USERANGLES,  BB_INPUT, SURFBB, USER_EMISSIVITY,         & ! Inputs (Sleave/Thermal)
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT,                         & ! Inputs (Optical - Lin Atmos)
        L_FMATRIX_UP, L_FMATRIX_DN, L_TRUNC_FACTOR, LS_USER_EMISSIVITY,                  & ! Inputs (Optical - Lin Atmos)
        LS_EXACTDB_BRDFUNC, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,               & ! Inputs (Optical - Lin Surf)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output - Stokes vectors
        FO_COLUMNWF_SS,  FO_COLUMNWF_DB, FO_COLUMNWF_DTA, FO_COLUMNWF_DTS,                  & ! Output - Column Jacobians
        FO_SURFACEWF_DB, FO_SURFACEWF_DTS,                                                  & ! Output - Surface Jacobians
        FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES,                                         & ! Output - Stokes composites
        FO_COLUMNWF_ATMOS, FO_COLUMNWF_SURF, FO_COLUMNWF, FO_SURFACEWF,                     & ! Output - Jacobian composites
        CUMTRANS, SLTERM, LC_CUMTRANS, LSSL_SLTERM )                                          ! Output - Auxiliary

!  4/15/20. Version 2.8.2. Drop MAXMOMENTS, MAXMOMENTS_INPUT, MAXFINELAYERS

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                 MAXLAYERS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS, MAX_GEOMETRIES,   &       
                                 MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS, MAXSTOKES_SQ,          &
                                 MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO, ONE, SMALLNUM

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
      LOGICAL, INTENT(IN) ::            DO_WATER_LEAVING    ! 4/9/19 Added
      LOGICAL, INTENT(IN) ::            DO_SL_ISOTROPIC

!  Linearization flags. Sleave flag 8/3/16.

      LOGICAL, INTENT(IN) ::            DO_COLUMNWFS
      LOGICAL, INTENT(IN) ::            DO_SURFACEWFS
      LOGICAL, INTENT(IN) ::            DO_SLEAVEWFS

!  Control integers
!    ==> 4/15/20. Version 2.8.2. Add Do_Doublet Offset

      INTEGER, INTENT(IN) ::            NSTOKES, NLAYERS, NGEOMS, NSZAS, NVZAS, NAZMS
      INTEGER, INTENT(IN) ::            NA_OFFSET ( MAXBEAMS, MAX_USER_STREAMS )
      INTEGER, INTENT(IN) ::            ND_OFFSET ( MAXBEAMS )

!  Type structure, FOGeometry input
  
      TYPE(VLIDORT_Geometry_FO), INTENT(IN)      :: FOGeometry

!  Linearization control. 
!    Sleave control 8/3/16. Note that N_SURFACEWFS = N_REFLECWFS + N_SLEAVEWFS

      INTEGER, INTENT(IN) ::            N_COLUMNWFS
      LOGICAL, INTENT(IN) ::            LVARYFMAT  ( MAXLAYERS, MAX_ATMOSWFS )
      INTEGER, INTENT(IN) ::            N_SURFACEWFS
      INTEGER, INTENT(IN) ::            N_REFLECWFS
      INTEGER, INTENT(IN) ::            N_SLEAVEWFS

!  require the Level Mask inputs

      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           USER_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           USER_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs, added 9/17/16

      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

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

      DOUBLE PRECISION, INTENT(IN) ::   BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   SURFBB
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  surface BRDF

      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Surface leaving inputs, 8/3/16, for Version 2.8

      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized optical properties

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

!mick fix 9/19/2017 - swapped layer & geo indices

      !DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )

!  Linearized surface properties

      DOUBLE PRECISION, INTENT(IN) ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  Surface leaving linearizations. New for Version 2.8, 8/3/16

      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  Outputs
!  =======
!mick mod 9/19/2017 - added FO_STOKES_ATMOS, FO_STOKES_SURF, FO_COLUMNWF_ATMOS, FO_COLUMNWF_SURF (new output from vector FO code)

!  Standard
!  --------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  4/9/19. Additional output for the sleave correction

      DOUBLE PRECISION, INTENT (INOUT)  :: CUMTRANS ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  4/9/19. Surface leaving FO assignation

      DOUBLE PRECISION, INTENT (INOUT)  :: SLTERM (  MAXSTOKES, MAX_GEOMETRIES )

!  Linearized
!  ----------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_SS &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DB &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DB &
                                           ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DTA &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DTS &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DTS &
                                           ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_ATMOS &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_SURF &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF &
                                           ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF &
                                           ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  4/9/19. Additional linearized output for the sleave correction

      DOUBLE PRECISION, INTENT (INOUT) :: LC_CUMTRANS ( MAX_ATMOSWFS, max_user_levels, max_geometries )
      DOUBLE PRECISION, INTENT (INOUT) :: LSSL_SLTERM ( MAXSTOKES, MAX_SLEAVEWFS, MAX_GEOMETRIES )

!  Local variables
!  ===============
!mick fix 3/2/2020 - LS_EMISS restored for now
!                  - changed 2nd dim of LS_EMISS from MAX_GEOMETRIES back to MAX_USER_VZANGLES

!  Vector sunlight flag

      LOGICAL   :: DO_SUNLIGHT

!  Vector emissivity flag (12/11/17 Rob add.)
!      4/15/20. Version 2.8.2. Introduced as an input to the FO call

      LOGICAL   :: do_Polarized_Emissivity

!  Surface properties - reflective (could be the albedo)

      REAL(FFP) :: REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES )
      REAL(FFP) :: LS_REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES, MAX_SURFACEWFS )
      REAL(FFP) :: LS_EMISS  ( MAXSTOKES, MAX_USER_VZANGLES, MAX_SURFACEWFS )

!  optical

      REAL(FFP) :: EXTINCTION  ( MAXLAYERS )
      REAL(FFP) :: DELTAUS     ( MAXLAYERS )
      REAL(FFP) :: OMEGA       ( MAXLAYERS )
      REAL(FFP) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

      REAL(FFP) :: L_EXTINCTION ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_DELTAUS    ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_OMEGA      ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_TRUNCFAC   ( MAXLAYERS, MAX_ATMOSWFS )

!  Other variables

      integer   :: n, ns, nv, na, G, par, spar, GK, O1, O2
      INTEGER   :: LUM=1, LUA=1

!  Define VFO_LCS_MASTER inputs
!  ============================

!  4/15/20. Version 2.8.2.
!      -- No geometrical variables (all removed), now precalculated in FOGeometry
!      -- No GREEKMAT setting, no do_FMATRIX and no NGREEK_MOMENTS_INPUT
!      -- ssflux and na_offset are input directly (not re-defined here)
!      -- Exception handling no longer required
!      -- Add Doublet geometry control

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

!  Version 2.8. Copying not needed
!      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

!  BRDF reflection now also for Lattice.
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet option for BRDF reflection.

      reflec = ZERO
      if (do_lambertian_surface) then
         reflec(1,1,1:ngeoms) = LAMBERTIAN_ALBEDO
      else
        if ( do_ObsGeom ) then
           do g=1,ngeoms
            do o1=1,nstokes ; do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,lum,lua,g)
            end do ; end do
           end do
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

!  Surface leaving, introduced 8/3/16. Proper Initializing
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

!  Linearized control
!  ==================

!  4/15/20. Version 2.8.2. This stuff removed.
      !Default setting: see below for modification when Greekmat moments vary
!      Lvarymoms(1:nlayers,1:max_atmoswfs) = .FALSE.

!  Linearized optical

      !Atmospheric quantities
      L_extinction = ZERO
      L_deltaus    = ZERO
      L_omega      = ZERO
      L_truncfac   = ZERO

!  Recall that VLIDORT takes FULLY     NORMALIZED linearized atmospheric inputs (x/y)*(dy/dx)
!  But FO code    only takes PARTIALLY NORMALIZED linearized atmospheric inputs x*(dy/dx)
!    -- thus, we must compensate for this difference when formulating the inputs here!

      do n=1,nlayers
        if (do_columnwfs) then
          do par=1,n_columnwfs

            if (do_deltam_scaling) then
              L_extinction(n,par) = (DELTAU_VERT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT(N) * L_DELTAU_VERT(PAR,N)
            else
              L_extinction(n,par) = (DELTAU_VERT_INPUT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT_INPUT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT_INPUT(N)*L_DELTAU_VERT_INPUT(PAR,N)
            endif

            L_omega(n,par)      = OMEGA_TOTAL_INPUT(N) * L_OMEGA_TOTAL_INPUT(PAR,N)

            if (do_deltam_scaling) L_truncfac(n,par) = L_TRUNC_FACTOR(PAR,N)

!  4/15/20. Version 2.8.2. This stuff removed.
!  Version 2.8 (FO 1.5). Optional use of Fmatrix input
!    Linearized Greekmat will not be needed in this case
!  Check for variation of Greek moments associated with
!    Jacobian wrt current atmospheric parameter
!mick fix 9/19/2017 - added .not. to the block IF condition (i.e. using traditional GMat
!                     moments, NOT the new FMat facility) 
!                   - changed linearized GMat moment IF condition
!                   - added ELSE section for proper passing of linearized FO F-matrix input
!            if ( .not. do_Fmatrix ) then
!              do o1=1,nstokes ; do o2=1,nstokes
!                  gk = 4*(o1-1) + o2
!                  do l=0,nmoments_input
!                    L_greekmat(n,l,o1,o2,par) = GREEKMAT_TOTAL_INPUT(L,N,GK) * L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK)
!                    !if ( (l==1).and.(gk==1).and.(L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK) > 1.0d-8) ) Lvarymoms(n,par) = .TRUE.
!                    if ( (gk==1).and.(ABS(L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK)) >= 1000.0d0*SMALLNUM) ) Lvarymoms(n,par) = .TRUE.
!                  end do
!               end do ; end do
!            else
!               do g=1,ngeoms
!                  FO_L_FMATRIX_UP(n,g,1:6,par) = L_FMATRIX_UP(par,n,g,1:6)
!                  FO_L_FMATRIX_DN(n,g,1:6,par) = L_FMATRIX_DN(par,n,g,1:6)
!               end do
!            endif

          end do
        end if
      end do

!  Surface quantities
!  ------------------

!  Initialize (LS Emiss not required)

      LS_reflec   = ZERO
      LSSL_slterm = ZERO

!  reflectance
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet option for Linearized reflectance.

      if (do_surfacewfs) then
        if (do_lambertian_surface) then
          do spar=1,n_reflecwfs
            LS_reflec(1,1,1:ngeoms,spar) = ONE
          end do
        else
          do spar=1,n_reflecwfs
            if ( do_ObsGeom ) then
              do g=1,ngeoms
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,LUM,LUA,G)
                end do ; end do
              end do
            else if ( do_Doublet ) then
              do ns = 1, nszas ; do nv = 1, nvzas
                g = nd_offset(ns) + nv
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,NV,LUA,NS)
                end do ; end do
              end do ; end do
            else
              do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                g = na_offset(ns,nv) + na
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,NV,NA,NS)
                end do ; end do
              end do ; end do ; end do
            endif
          end do
        end if
      endif

!mick fix 9/19/2017 - added "do_surface_emission" to IF condition
!  - polarized surface emissivity. 12/11/17 Rob Add.
!  4/15/20. Version 2.8.2. This Emissivity stuff removed
!mick fix 3/2/2020  - turned back on due to the dimension swapping here and
!                     consistency of "LS_emiss" dimensions with other similar
!                     arrays inside VFO_RTCALC_LCS_MASTER

      if ( do_surface_emission .and. do_surfacewfs ) then
        do spar=1,n_reflecwfs
          if ( do_ObsGeom ) then
            do g=1,ngeoms
              LS_emiss(1:nstokes,g,spar) = LS_USER_EMISSIVITY(spar,1:nstokes,g)
            end do
          else
            do nv=1,nvzas
             LS_emiss(1:nstokes,nv,spar) = LS_USER_EMISSIVITY(spar,1:nstokes,nv)
            end do
          endif
        end do
      end if

!  Surface leaving, New section 8/3/16 for Version 2.8.  4/9/19 ==> This is now an output.
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet option for Linearized surface leaving.

      if ( do_surface_leaving .and. do_sleavewfs ) then
         do spar = 1, n_sleavewfs
            if ( do_ObsGeom ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_SLTERM(1,1:ngeoms,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,1:ngeoms)
               ELSE
                  DO O1 = 1, NSTOKES
                     LSSL_SLTERM(O1,1:ngeoms,spar) = LSSL_SLTERM_USERANGLES(spar,O1,LUM,LUA,1:ngeoms)
                  ENDDO
               ENDIF
            else if ( do_Doublet ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas
                     g = nd_offset(ns) + nv
                     LSSL_SLTERM(1,g,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,ns)
                  end do ; end do
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas 
                     g = nd_offset(ns) + nv
                     do o1 = 1,nstokes
                        LSSL_SLTERM(O1,g,spar) = LSSL_SLTERM_USERANGLES(spar,O1,nv,LUA,ns)
                     end do
                  end do ; end do
               ENDIF
            else
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     LSSL_SLTERM(1,g,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,ns)
                  end do ; end do ; end do 
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     do o1 = 1,nstokes
                        LSSL_SLTERM(O1,g,spar) = LSSL_SLTERM_USERANGLES(spar,O1,nv,na,ns)
                     end do
                  end do ; end do ; end do 
               ENDIF
            endif
         enddo
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

!  Call VFO_LCS_MASTER
!  ===================

!  Upgraded to FO Version 1.5, 7/11/16
!   4/9/19. CUMTRANS and LC_CUMTRANS now additional output. add water-leaving control.

!  4/15/20. Version 2.8.2.
!      -- Name of subroutine changed. Input arguments rearranged somewhat
!      -- FOGeometry structure passed instead of many geometry variables
!      -- Removed: GREEKMAT, DO_FMATRIX, NGREEK_MOMENTS_INPUT.  flux --> ssflux
!      -- Dimensioning removed; now not so stand-alone!!!!
!      -- Exception handling removed
!      -- Add doublet geometry flag and offset to input list

      CALL VFO_RTCALC_LCS_MASTER &
       ( do_solar_sources, do_thermal_emission, do_surface_emission, do_Polarized_Emissivity,     & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_lambertian_surface, do_surface_leaving, do_water_leaving, & ! Input flags (dir/surface)
         do_sunlight, do_deltam_scaling, do_obsgeom, do_doublet, do_partlayers, do_planpar,       & ! Input flags (surface/geoms)
         do_enhanced_ps, do_columnwfs, do_surfacewfs, do_sleavewfs,                               & ! Input Lin flags
         nstokes, ngeoms, nszas, nvzas, nazms, na_offset, nd_offset, nlayers,                     & ! Input Numbers/Offsets
         n_user_levels, user_level_mask_up, user_level_mask_dn,                                   & ! Input control-levels
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx, FOGeometry,  & ! Input control-partial/Geometry
         n_columnwfs, LvaryFmat, n_reflecwfs, n_sleavewfs, n_surfacewfs,                          & ! Input Lin control
         ssflux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,           & ! Input atmos optical
         bb_input, surfbb, user_emissivity, reflec, slterm,                                       & ! Input thermal/surface
         L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, L_fmatrix_dn,                & ! Input Lin atmos optical
         LS_emiss, LS_reflec, LSSL_slterm,                                                        & ! Input Lin thermal/surf
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                                & ! Output - Stokes   Vectors
         fo_columnwf_ss,  fo_columnwf_db,  fo_columnwf_dta, fo_columnwf_dts,                      & ! Output - Column Jacobians
         fo_surfacewf_db,  fo_surfacewf_dts,                                                      & ! Output - Surface  Jacobians
         fo_stokes_atmos,  fo_stokes_surf, fo_stokes,                                             & ! Output - Stokes   composites
         fo_columnwf_atmos, fo_columnwf_surf, fo_columnwf, fo_surfacewf,                          & ! Output - Jacobian composites
         cumtrans, LC_cumtrans )                                                                    ! Output - ancillary

!  Done

      RETURN
END SUBROUTINE VFO_LCS_MASTER_INTERFACE

!  finish module

END MODULE vlidort_vfo_lcs_interface_m


