
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

      MODULE VLIDORT_LinInputs_def

!  This Module contains the following VLIDORT input Structures, with Intents :

!        VLIDORT_Fixed_LinControl    nested in VLIDORT_Fixed_LinInputs
!        VLIDORT_Fixed_LinOptical    nested in VLIDORT_Fixed_LinInputs
!         VLIDORT_Fixed_LinInputs    Intent(In)

!     VLIDORT_Modified_LinControl    nested in VLIDORT_Modified_LinInputs
!      VLIDORT_Modified_LinInputs    Intent(InOut)

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinControl

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, dimension(MAXLAYERS)  :: TS_LAYER_VARY_FLAG
      INTEGER, dimension(MAXLAYERS)  :: TS_LAYER_VARY_NUMBER

!  Total number of column Jacobians

      INTEGER  :: TS_N_TOTALCOLUMN_WFS

!  Total number of profile Jacobians

      INTEGER  :: TS_N_TOTALPROFILE_WFS

!  Total number of surface Jacobians

      INTEGER  :: TS_N_SURFACE_WFS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      INTEGER  :: TS_N_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Jacobian names

      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_COLUMNWF_NAMES
      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_PROFILEWF_NAMES

      END TYPE VLIDORT_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinOptical


!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: &
        TS_L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: &
        TS_L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, &
        MAXLAYERS, MAXSTOKES_SQ) :: TS_L_GREEKMAT_TOTAL_INPUT


      END TYPE VLIDORT_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinInputs


      TYPE(VLIDORT_Fixed_LinControl)    :: Cont
      TYPE(VLIDORT_Fixed_LinOptical)    :: Optical


      END TYPE VLIDORT_Fixed_LinInputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_LinControl

!  Control linearization

      LOGICAL  :: TS_DO_COLUMN_LINEARIZATION
      LOGICAL  :: TS_DO_PROFILE_LINEARIZATION
      LOGICAL  :: TS_DO_ATMOS_LINEARIZATION

      LOGICAL  :: TS_DO_SURFACE_LINEARIZATION
      LOGICAL  :: TS_DO_LINEARIZATION

!  This flag moved from Fixed Lin Control, Version 2.7

      LOGICAL  :: TS_DO_SIMULATION_ONLY

!  BlackBody Jacobian Flags, Introduced March 26th 2014, Version 2.7

      LOGICAL  :: TS_DO_ATMOS_LBBF
      LOGICAL  :: TS_DO_SURFACE_LBBF

!  These two flags have been superseded in Version 2.7
!      LOGICAL  :: TS_DO_SURFBB_LINEARIZATION  !  REPLACED BY TS_DO_SURFACE_LBBF
!      LOGICAL  :: TS_DO_LTE_LINEARIZATION     !  REPLACED BY TS_DO_ATMOS_LBBF 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      LOGICAL  :: TS_DO_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      END TYPE VLIDORT_Modified_LinControl

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_LinInputs


      TYPE(VLIDORT_Modified_LinControl)    :: MCont


      END TYPE VLIDORT_Modified_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Fixed_LinControl, &
                VLIDORT_Fixed_LinOptical, &
                VLIDORT_Fixed_LinInputs, &
                VLIDORT_Modified_LinControl, &
                VLIDORT_Modified_LinInputs

   END MODULE VLIDORT_LinInputs_def
