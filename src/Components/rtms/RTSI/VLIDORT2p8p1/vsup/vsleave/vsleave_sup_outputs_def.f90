
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

module VSLEAVE_Sup_Outputs_def_m

!  This module contains the following structures:

!  VSLEAVE_Sup_Outputs - Intent(In)  for VLIDORT,
!                        Intent(Out) for VSLEAVE_Sup

!mick mod 9/19/2017 - added type structure VSLEAVE_Output_Exception_Handling

use VLIDORT_PARS_m

implicit none

! #####################################################################
! #####################################################################

type VSLEAVE_Sup_Outputs

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), dimension ( MAXSTOKES, MAXBEAMS ) :: SL_SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: SL_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, &
        MAXBEAMS )   :: SL_SLTERM_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, &
        MAXBEAMS )   :: SL_USER_SLTERM_F_0

!   New Diagnostic, R. Spurr, 05 October 2015.

!     This is the Total Atmospheric flux (diffuse and direct) at BOA
!     for a Rayleigh atmosphere. It has a Stokes 4-vector form!
!     Flux is specified at the SLEAVE input wavelength, and is a
!     function of solar zenith angle.

!  In SLEAVE, this quantity comes either from the Gordon/Wang (1994)
!     approximate formula, or from a look-up table created offline by
!     running LIDORT for a 35-layer Rayleigh atmosphere for 15 SZAs
!     (0-88 degs.), and for wavelengths 270-900 nm @ 10 nm intervals

       REAL(fpk), dimension ( MAXBEAMS )  :: SL_TRANS_ATMOS

end type VSLEAVE_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE VSLEAVE_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: SL_STATUS_INPUTREAD
      INTEGER      :: SL_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_INPUTACTIONS

      END TYPE VSLEAVE_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VSLEAVE_Output_Exception_Handling

!  Exception handling for Output. New code, 19 Sep 2017 (mick add)
!     Message Length should be at least 120 Characters

      INTEGER      :: SL_STATUS_OUTPUT
      INTEGER      :: SL_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_OUTPUTMESSAGES

      END TYPE VSLEAVE_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VSLEAVE_Sup_Outputs, &
             VSLEAVE_Input_Exception_Handling, &
             VSLEAVE_Output_Exception_Handling

end module VSLEAVE_Sup_Outputs_def_m