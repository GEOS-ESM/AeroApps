
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


module VSLEAVE_Sup_Outputs_def

!  This module contains the following structures:

!  VSLEAVE_Sup_Outputs - Intent(In)  for VLIDORT,
!                        Intent(Out) for VSLEAVE_Sup

use VLIDORT_PARS

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

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VSLEAVE_Sup_Outputs, &
             VSLEAVE_Input_Exception_Handling

end module VSLEAVE_Sup_Outputs_def
