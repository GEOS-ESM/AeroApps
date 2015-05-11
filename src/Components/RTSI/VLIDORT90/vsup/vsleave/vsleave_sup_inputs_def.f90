
! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

module VSLEAVE_Sup_Inputs_def

!  Version 2.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  This module contains the following structures:

!  VSLEAVE_Sup_Inputs - Intent(In) for VSLEAVE_Sup

use VLIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type VSLEAVE_Sup_Inputs

!  General control variables
!  -------------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: SL_DO_SLEAVING

!  Isotropic flag

      LOGICAL :: SL_DO_ISOTROPIC

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: SL_DO_EXACT
      LOGICAL :: SL_DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: SL_DO_FLUORESCENCE

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: SL_DO_SOLAR_SOURCES
      LOGICAL :: SL_DO_USER_OBSGEOMS

!  Geometry and integer control
!  ----------------------------

!  Number of Stokes components

      INTEGER :: SL_NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: SL_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: SL_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: SL_BEAM_SZAS

!  user-defined relative azimuths

      INTEGER                                 :: SL_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: SL_USER_RELAZMS

!  User-defined zenith angle input 

      LOGICAL                                 :: SL_DO_USER_STREAMS
      INTEGER                                 :: SL_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: SL_USER_ANGLES_INPUT

!  !@@ Observational geometry inputs

      INTEGER                                    :: SL_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: SL_USER_OBSGEOMS

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SL_SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: SL_CHLORCONC

!  Input wavelength in [Microns]

      REAL(fpk) :: SL_WAVELENGTH

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: SL_WINDSPEED, SL_WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: SL_NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: SL_DO_GlintShadow
      LOGICAL   :: SL_DO_FoamOption
      LOGICAL   :: SL_DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: SL_FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: SL_FL_Latitude, SL_FL_Longitude

!  Input Epoch

      INTEGER :: SL_FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: SL_FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: SL_FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: SL_FL_InputGAUSSIANS(3,2)

end type VSLEAVE_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VSLEAVE_Sup_Inputs

end module VSLEAVE_Sup_Inputs_def

