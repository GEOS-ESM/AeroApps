
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

      MODULE VLIDORT_Sup_BRDF_def

!  This module contains the following structures:

!     VLIDORT_Sup_BRDF    Intent(In) for VLIDORT,
!                         Intent(Out) for VLIDORT BRDFSup

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Sup_BRDF


!  Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), dimension ( MAXSTOKES_SQ, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: TS_EXACTDB_BRDFUNC

!  Fourier components of BRDF, in the following order (same all threads)
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, &
        MAXBEAMS )   :: TS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, &
        MAXSTREAMS ) :: TS_BRDF_F
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, &
        MAXBEAMS )   :: TS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, &
        MAXSTREAMS ) :: TS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAXSTOKES, MAXSTREAMS )       :: TS_EMISSIVITY
      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS ) :: TS_USER_EMISSIVITY


      END TYPE VLIDORT_Sup_BRDF

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Sup_BRDF

      END MODULE VLIDORT_Sup_BRDF_def
