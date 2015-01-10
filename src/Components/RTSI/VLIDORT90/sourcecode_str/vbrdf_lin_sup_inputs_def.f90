
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
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

module VBRDF_LinSup_Inputs_def

!  This module contains the following structures:

!  VBRDF_LinSup_Inputs - Intent(In) for VBRDF_LinSup

use VLIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type VBRDF_LinSup_Inputs

!  Linearizaion material
!  ---------------------

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL, dimension ( MAX_BRDF_KERNELS ) :: &
        BS_DO_KERNEL_FACTOR_WFS
      LOGICAL, dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: &
        BS_DO_KERNEL_PARAMS_WFS

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, dimension ( MAX_BRDF_KERNELS )  :: BS_DO_KPARAMS_DERIVS

!  number of surface weighting functions

      INTEGER    :: BS_N_SURFACE_WFS
      INTEGER    :: BS_N_KERNEL_FACTOR_WFS
      INTEGER    :: BS_N_KERNEL_PARAMS_WFS

end type VBRDF_LinSup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VBRDF_LinSup_Inputs

end module VBRDF_LinSup_Inputs_def

