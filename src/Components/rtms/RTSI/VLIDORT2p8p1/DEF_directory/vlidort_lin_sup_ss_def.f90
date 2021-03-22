
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


      MODULE VLIDORT_LinSup_SS_def_m

!  This module contains the following structures,
!  with intents :

!       VLIDORT_LinSup_SS_Col    nested in VLIDORT_LinSup_SS_InOut
!      VLIDORT_LinSup_SS_Prof    nested in VLIDORT_LinSup_SS_InOut
!      VLIDORT_LinSup_SS_Surf    nested in VLIDORT_LinSup_SS_InOut
!     VLIDORT_LinSup_SS_InOut    Intent(InOut)

      USE VLIDORT_PARS_m, Only : fpk, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS, &
                                 MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS 

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Col


!  SS atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_COLUMNWF_SS

!  DB atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES ) :: TS_COLUMNWF_DB


      END TYPE VLIDORT_LinSup_SS_Col

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Prof


!  SS atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_PROFILEWF_SS

!  DB atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES ) :: TS_PROFILEWF_DB


      END TYPE VLIDORT_LinSup_SS_Prof

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS_Surf


!  SS surface weighting functions

!      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
!        MAXSTOKES, MAX_DIRECTIONS ) :: TS_SURFACEWF_SS

!  DB surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES ) :: TS_SURFACEWF_DB


      END TYPE VLIDORT_LinSup_SS_Surf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSup_SS


      TYPE(VLIDORT_LinSup_SS_Col)   :: Col
      TYPE(VLIDORT_LinSup_SS_Prof)  :: Prof
      TYPE(VLIDORT_LinSup_SS_Surf)  :: Surf


      END TYPE VLIDORT_LinSup_SS

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_LinSup_SS_Col,&
                VLIDORT_LinSup_SS_Prof,&
                VLIDORT_LinSup_SS_Surf,&
                VLIDORT_LinSup_SS

      END MODULE VLIDORT_LinSup_SS_def_m
