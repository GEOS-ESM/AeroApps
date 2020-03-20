! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

      module LIDORT_Sup_BRDF_def

!  Version 3.7, Internal threading removed, 02 May 2014

!  This module contains the following Structures

!      LIDORT_Sup_BRDF      Intent(In) for LIDORT,
!                           Intent(Out) for LIDORT BRDFSup

      use LIDORT_PARS, only : fpk, MAXMOMENTS, MAXSTREAMS, MAXBEAMS, &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS

      implicit none

! #####################################################################
! #####################################################################

      type LIDORT_Sup_BRDF

!  Exact (direct bounce) BRDF

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: TS_EXACTDB_BRDFUNC

!  Fourier components of BRDF, in the following order
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )   :: TS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS ) :: TS_BRDF_F
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )   :: TS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS ) :: TS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAXSTREAMS )       :: TS_EMISSIVITY
      REAL(fpk), dimension ( MAX_USER_STREAMS ) :: TS_USER_EMISSIVITY

      end type LIDORT_Sup_BRDF

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Sup_BRDF

      end module LIDORT_Sup_BRDF_def
