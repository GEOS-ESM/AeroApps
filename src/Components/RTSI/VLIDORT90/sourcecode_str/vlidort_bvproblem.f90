! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            BVP_MATRIXSETUP_MASTER      (master)             #
! #             BVP_SURFACE_SETUP_HOM                           #
! #             BVP_MATRIX_INIT                                 #
! #             BVP_MATRIX_SETUP                                #
! #             BVP_MATRIX_SVD                                  #
! #                                                             #
! #            BVP_SOLUTION_MASTER      (master)                #
! #             BVP_SURFACE_SETUP_BEAM                          #
! #             BVP_COLUMN_SETUP                                #
! #             BVP_BACKSUB                                     #
! #                                                             #
! # Telescped BVP: Subroutines in this Module                   #
! #                                                             #
! #            BVPTEL_MATRIXSETUP_MASTER      (master)          #
! #             BVPTEL_MATRIX_INIT                              #
! #             BVPTEL_MATRIX_INIT_OMP                          #
! #             BVPTEL_MATRIX_SETUP                             #
! #             BVPTEL_MATRIX_SVD                               #
! #                                                             #
! #            BVPTEL_SOLUTION_MASTER      (master)             #
! #             BVPTEL_COLUMN_SETUP                             #
! #             BVPTEL_BACKSUB                                  #
! #                                                             #
! #            BMAT_ROWMASK    (integer function)               #
! #            BTELMAT_ROWMASK (integer function)               #
! #                                                             #
! ###############################################################


      MODULE vlidort_bvproblem

      PRIVATE
      PUBLIC :: BVP_MATRIXSETUP_MASTER, &
                BVP_SOLUTION_MASTER, &
                BVPTEL_MATRIXSETUP_MASTER, &
                BVPTEL_SOLUTION_MASTER

      CONTAINS

      SUBROUTINE BVP_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, NSTKS_NSTRMS, &
        MUELLER_INDEX, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        NSTKS_NSTRMS_2, &
        T_DELT_EIGEN, &
        R2_HOMP, R2_HOMM, AXBID_F, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::            FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::            K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::   T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) ::  BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (OUT) ::           IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER :: STATUS_SUB

      INTEGER :: NMIN ( MAXTOTAL )
      INTEGER :: NMAX ( MAXTOTAL )
      INTEGER :: KALL

!  Initialize exception handling

      STATUS  = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, NSTKS_NSTRMS, &
        MUELLER_INDEX, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        R2_HOMP, R2_HOMM, &
        AXBID_F )

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT ( &
        NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        NSTOKES, NSTREAMS, &
        NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        T_DELT_EIGEN, &
        R2_HOMP, R2_HOMM, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  SVD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_SVD ( &
        NLAYERS, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, SMAT2, &
        IPIVOT, SIPIVOT, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVP_MATRIXSETUP_MASTER

!

      SUBROUTINE BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, NSTKS_NSTRMS, &
        MUELLER_INDEX, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        R2_HOMP, R2_HOMM, &
        AXBID_F )

!  This is the Lambertian or BRDF surface routine
!     Reflected homogeneous solutions

!  merged July 26, 2010. Using new BRDF material
!   RT SOLUTIONS Inc,

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  local variables
!  ---------------

      INTEGER ::         I,J,O1,O2,OM,K,KO1,K0,K1,K2,M,NL
      DOUBLE PRECISION :: AXJ, KMULT, FACTOR_KERNEL
      DOUBLE PRECISION :: H_1,   H_2,   H_1_S,   H_2_S
      DOUBLE PRECISION :: H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION :: H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO K = 1, NSTKS_NSTRMS
          DO O1 = 1, NSTOKES
            R2_HOMP(I,O1,K) = ZERO
            R2_HOMM(I,O1,K) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Fourier component

      M = FOURIER_COMPONENT

!  Last layer

      NL = NLAYERS

!  Offset

      KO1 = K_REAL(NL) + 1

!  Return if Fourier component not zero for Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE &
           .and. FOURIER_COMPONENT .GT. 0 ) RETURN

!  Go to total BRDF case for non-Lmabertian surface

      IF ( .not.  DO_LAMBERTIAN_SURFACE ) goto 877

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

!  Integrate Downward streams of particular solutions

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  Homogeneous real solutions

      DO K = 1, K_REAL(NLAYERS)
        DO I = 1, NSTREAMS
          H_1 = ZERO
          H_2 = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1 = H_1 + AXJ*SOLA_XPOS(J,1,K,NLAYERS)
            H_2 = H_2 + AXJ*SOLB_XNEG(J,1,K,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K) = H_1 * KMULT
          R2_HOMM(I,1,K) = H_2 * KMULT
        ENDDO
      ENDDO

!  Homogeneous complex solutions

      KO1 = K_REAL(NLAYERS) + 1
      DO K = 1, K_COMPLEX(NLAYERS)
        K0 = 2*K - 2
        K1 = KO1 + K0
        K2 = K1  + 1
        DO I = 1, NSTREAMS
          H_1_CR = ZERO
          H_2_CR = ZERO
          H_1_CI = ZERO
          H_2_CI = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1_CR = H_1_CR + AXJ*SOLA_XPOS(J,1,K1,NLAYERS)
            H_2_CR = H_2_CR + AXJ*SOLB_XNEG(J,1,K1,NLAYERS)
            H_1_CI = H_1_CI + AXJ*SOLA_XPOS(J,1,K2,NLAYERS)
            H_2_CI = H_2_CI + AXJ*SOLB_XNEG(J,1,K2,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K1) = H_1_CR * KMULT
          R2_HOMM(I,1,K1) = H_2_CR * KMULT
          R2_HOMP(I,1,K2) = H_1_CI * KMULT
          R2_HOMM(I,1,K2) = H_2_CI * KMULT
        ENDDO
      ENDDO

!  debug Lambertian

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          DO K = 1, NSTKS_NSTRMS
!           write(*,'(3i3,1p2e14.6)')
!     &             I,O1,K,R2_HOMP(I,O1,K),R2_HOMM(I,O1,K)
!          ENDDO
!        ENDDO
!      ENDDO
!      pause'r2homp'

!  return from Lambertian

      RETURN

!  Total BRDF case
!  ===============

!  Continuation point

877   continue

!  Integrate Downward streams of particular solutions

      FACTOR_KERNEL = SURFACE_FACTOR

!  help variable

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO J = 1, NSTREAMS
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
!mick fix - 2/15/2011 - reactivated 1st statement; commented out 2nd statement
              AXBID_F(I,J,OM) = QUAD_STRMWTS(J) * BRDF_F(M,OM,I,J)
!              AXBID_F(I,J,OM) = QUAD_STRMWTS(J) * BRDF_F(M,OM,J,I)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  homogeneous real solutions

      DO K = 1, K_REAL(NL)
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            H_1 = ZERO
            H_2 = ZERO
            DO J = 1, NSTREAMS
              H_1_S = ZERO
              H_2_S = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                H_1_S = H_1_S + AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K,NL)
                H_2_S = H_2_S + AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K,NL)
              ENDDO
              H_1 = H_1 + H_1_S
              H_2 = H_2 + H_2_S
            ENDDO
            R2_HOMP(I,O1,K) = FACTOR_KERNEL * H_1
            R2_HOMM(I,O1,K) = FACTOR_KERNEL * H_2
          ENDDO
        ENDDO
      ENDDO

!  homogeneous complex solutions

      DO K = 1, K_COMPLEX(NL)
        K0 = 2*K - 2
        K1 = KO1 + K0
        K2 = K1  + 1
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            H_1_CR = ZERO
            H_2_CR = ZERO
            H_1_CI = ZERO
            H_2_CI = ZERO
            DO J = 1, NSTREAMS
              H_1_S_CR = ZERO
              H_2_S_CR = ZERO
              H_1_S_CI = ZERO
              H_2_S_CI = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                H_1_S_CR = H_1_S_CR + &
                    AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K1,NL)
                H_2_S_CR = H_2_S_CR + &
                    AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K1,NL)
                H_1_S_CI = H_1_S_CI + &
                    AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K2,NL)
                H_2_S_CI = H_2_S_CI + &
                    AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K2,NL)
              ENDDO
              H_1_CR = H_1_CR + H_1_S_CR
              H_2_CR = H_2_CR + H_2_S_CR
              H_1_CI = H_1_CI + H_1_S_CI
              H_2_CI = H_2_CI + H_2_S_CI
            ENDDO
            R2_HOMP(I,O1,K1) = FACTOR_KERNEL * H_1_CR
            R2_HOMM(I,O1,K1) = FACTOR_KERNEL * H_2_CR
            R2_HOMP(I,O1,K2) = FACTOR_KERNEL * H_1_CI
            R2_HOMM(I,O1,K2) = FACTOR_KERNEL * H_2_CI
          ENDDO
        ENDDO
      ENDDO

!  debug total  BRDF case

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          DO K = 1, NSTKS_NSTRMS
!           write(*,'(3i3,1p2e14.6)')
!     &             2,1,K,R2_HOMP(2,1,K),R2_HOMM(2,1,K)
!          ENDDO
!        ENDDO
!      ENDDO
!      pause'r2homp'

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_HOM

!

      SUBROUTINE BVP_MATRIX_INIT ( &
        NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  Initialise the compressed matrix

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

      INTEGER, INTENT (OUT) ::          NMIN ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAX ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALL
      DOUBLE PRECISION, INTENT (OUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  Local variables
!  ---------------

!  compression function

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  help variables

      INTEGER :: I, J, N3, JS, JF, IS, LF, L, I1, NMJ, NXJ

!  special case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  compression row indices, NMIN, NMAX and KALL

      DO J = 1, N_SUPDIAG + 1
        NMIN(J) = 1
      ENDDO
      DO J = N_SUPDIAG + 2, NTOTAL
        NMIN(J) = J - N_SUPDIAG
      ENDDO
      DO J = 1, NTOTAL - N_SUBDIAG
        NMAX(J) = J + N_SUBDIAG
      ENDDO
      DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
        NMAX(J) = NTOTAL
      ENDDO
      KALL = N_SUBDIAG + N_SUPDIAG + 1

!  Former code to assign Array - now replaced by external function
!    Jukaa Kujanpaa, FMI, August 2005
!      DO I = 1, NTOTAL
!        DO J = 1, NTOTAL
!          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
!            BMAT_ROWMASK(I,J) = KALL + I - J
!          ENDIF
!        ENDDO
!      ENDDO

!  compression matrix zeroing

      N3 = NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      LF = NLAYERS - 2

!  upper band top

      JS = NSTKS_NSTRMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTKS_NSTRMS
        DO J = JS, JF + I
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  upper band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 - NSTKS_NSTRMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTKS_NSTRMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
        JS = IS - N3 + 1
        JF = IS - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band bottom

      JS = LF * NSTKS_NSTRMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTKS_NSTRMS
      DO I = 1, NSTKS_NSTRMS
        I1 = I + IS
        DO J = JS + I, JF
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BVP_MATRIX_INIT

!

      SUBROUTINE BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        NSTOKES, NSTREAMS, &
        NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        T_DELT_EIGEN, &
        R2_HOMP, R2_HOMM, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  Fills up the compressed matrix directly

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      INTEGER, INTENT (IN) ::           NMIN ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           NMAX ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           KALL

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  compression function
!  --------------------

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  Local variables
!  ---------------

      INTEGER ::          I, I1, N, N1, C0, CE_OFFSET, O1, M
      INTEGER ::          EP, EM, EPC, EPS, EMS
      INTEGER ::          CP, CM, CEP, CEM, CEPS, CEMS, CEM1, CEP1
      INTEGER ::          ROWEP, ROWEM, ROWEPS, ROWEMS, ROWCEP, ROWCEM
      INTEGER ::          ROWCEPS, ROWCEMS, ROWCEP1, ROWCEM1
      INTEGER ::          IR, IROW, KO11, KO1, K0, K1, K2
      DOUBLE PRECISION :: XPNET, XMNET
      DOUBLE PRECISION :: X1R, X2R, X1I, X2I, X1R_H, X1I_H

!  Initialization

      EP = 0

!  Fourier component

      M = FOURIER_COMPONENT

!  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

!mick fix 2/17/11 - turn back on initialise (fix under test!)
!  Initialise - makes no difference
      DO EP = 1, NTOTAL
        DO EM = 1, MAXBANDTOTAL
          BANDMAT2(EM,EP) = ZERO
        ENDDO
      ENDDO

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS

            ROWEP = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)

            BANDMAT2(ROWEP,EP) = SOLA_XPOS(I,O1,EP,N)
            BANDMAT2(ROWEM,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  complex solutions (REWORKED)

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1

            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS

            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

            ROWEP  = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM  = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
            ROWEPS = BMAT_ROWMASK(IROW,EPS,NMIN(EPS),NMAX(EPS),KALL)
            ROWEMS = BMAT_ROWMASK(IROW,EMS,NMIN(EMS),NMAX(EMS),KALL)

            BANDMAT2(ROWEP,EP)   =   X1R
            BANDMAT2(ROWEPS,EPS) = - X1I
            BANDMAT2(ROWEM,EM)   =   X2R
            BANDMAT2(ROWEMS,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  intermediate layer boundaries
!  -----------------------------

      C0 = - NSTKS_NSTRMS
      DO N = 2, NLAYERS
        N1 = N - 1
        KO1  = K_REAL(N) + 1
        KO11 = K_REAL(N1) + 1
        C0   = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW

!  real solutions for layer above the boundary

            DO EP = 1, K_REAL(N1)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              ROWCEP = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP)   = &
                       T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
              BANDMAT2(ROWCEM,CEM)   = &
                         SOLB_XNEG(I,O1,EP,N1)
            ENDDO

!  real solutions for layer below the boundary
!   ( Note the change of sign !!! )

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEP1 = CEP + NSTKS_NSTRMS_2
              CEM1 = CEM + NSTKS_NSTRMS_2

              ROWCEP1 = BMAT_ROWMASK(CM,CEP1,NMIN(CEP1),NMAX(CEP1),KALL)
              ROWCEM1 = BMAT_ROWMASK(CM,CEM1,NMIN(CEM1),NMAX(CEM1),KALL)

              BANDMAT2(ROWCEP1,CEP1) = &
                         - SOLA_XPOS(I,O1,EP,N)
              BANDMAT2(ROWCEM1,CEM1) = &
                         - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
            ENDDO

!  complex solutions for layer above boundary

            DO EPC = 1, K_COMPLEX(N1)
              K0 = 2*EPC - 2
              K1 = KO11 + K0
              K2 = K1  + 1

              EP  = K_REAL(N1) + EPC
              EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1) &
                   -SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
              X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1) &
                   +SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
              X2R = SOLB_XNEG(I,O1,K1,N1)
              X2I = SOLB_XNEG(I,O1,K2,N1)

              ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP ,CEP)  =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM ,CEM)  =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I
            ENDDO

!  complex solutions for layer below boundary
!   ( Note the change of sign !!! )

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   = - X1R
              BANDMAT2(ROWCEPS,CEPS) =   X1I
              BANDMAT2(ROWCEM,CEM)   = - X2R
              BANDMAT2(ROWCEMS,CEMS) =   X2I
            ENDDO

          ENDDO
        ENDDO
      ENDDO

!  bottom BC (with albedo additions if flagged)

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      C0  = C0 + NSTKS_NSTRMS_2
      CE_OFFSET = C0 - NSTKS_NSTRMS

!  With albedo
!  -----------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)

              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM,CEM)   =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I

!         iunit = 20
!        IF (do_fdtest)iunit=19
!         write(iunit,4567)I,O1,EPC,
!     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N)+
!     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N),
!     &                   R2_HOMM(I,O1,EPC)
! 4567  format   (3i4,1p2e20.10)

            ENDDO

          ENDDO
        ENDDO

!  No albedo

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions
!    ----- Deep bug found for X1I, 17 December 2005.

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM,CEM)   =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I
           ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal completion

      RETURN

!  special case for 1 layer only
!  =============================

345   CONTINUE

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
            SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  REWORKED complex solutions

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            SMAT2(IROW,EP)  =   X1R
            SMAT2(IROW,EPS) = - X1I
            SMAT2(IROW,EM)  =   X2R
            SMAT2(IROW,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  bottom BC (with albedo additions)
!  ---------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  bottom BC (with No albedo)
!  --------------------------

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  REWORKED Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal return and finish

      RETURN
      END SUBROUTINE BVP_MATRIX_SETUP

!

      SUBROUTINE BVP_MATRIX_SVD ( &
        NLAYERS, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, SMAT2, &
        IPIVOT, SIPIVOT, &
        STATUS, MESSAGE, TRACE )

!  Solves the boundary value problem.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

      INTEGER, INTENT (OUT) ::            IPIVOT ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::            SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SVD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF &
           ( NTOTAL, NTOTAL, SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call (nlayers=1) in BVP_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_MATRIX_SVD

!

      SUBROUTINE BVP_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_SURFEMISS, FOURIER_COMPONENT, &
        IBEAM, SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        WLOWER, &
        EMISSIVITY, SURFBB, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DIRECT_BEAM, &
        WUPPER, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        N_SUBDIAG, N_SUPDIAG, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        AXBID_F, &
        COL2, SCOL2, &
        R2_BEAM, LCON, MCON, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFEMISS
      INTEGER, INTENT (IN) ::            FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::            IBEAM
      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::   WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SURFBB
      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::   DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::   WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      LOGICAL, INTENT (IN) ::            DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::            DO_FDTEST
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG
      INTEGER, INTENT (IN) ::            K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::            IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::   SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::            SIPIVOT ( MAXSTRMSTKS_2 )
      DOUBLE PRECISION, INTENT (IN) ::   AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!mick fix 7/29/2014 - placed COL2 & SCOL2 in call to make VLIDORT threadsafe
      DOUBLE PRECISION, INTENT (INOUT) ::  COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) ::  SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

      DOUBLE PRECISION, INTENT (INOUT) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

      INTEGER ::          STATUS_SUB

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  --Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        WLOWER, AXBID_F, &
        R2_BEAM )

!  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVP_COLUMN_SETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE, &
        DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS, &
        LAMBERTIAN_ALBEDO, &
        IBEAM, NSTOKES, NSTREAMS, &
        NLAYERS, EMISSIVITY, SURFBB, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        R2_BEAM, DIRECT_BEAM, &
        WUPPER, WLOWER, &
        COL2, SCOL2 )

!  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVP_BACKSUB ( &
        FOURIER_COMPONENT, IBEAM, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        NLAYERS, &
        NTOTAL, N_SUBDIAG, &
        N_SUPDIAG, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        COL2, SCOL2, &
        LCON, MCON, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVP_SOLUTION_MASTER

!

      SUBROUTINE BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        WLOWER, AXBID_F, &
        R2_BEAM )

!  This is the Lambertian or BRDF surface routine
!     Reflected beam solutions

!  merged July 26, 2010. Using new BRDF material
!   RT SOLUTIONS Inc,

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (INOUT) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      DOUBLE PRECISION :: REFL_B, REFL_B_S, FACTOR_KERNEL, AXJ
      INTEGER          :: I, J, O1, O2, OM, NL

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          R2_BEAM(I,O1)    = ZERO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if not the Fourier m = 0 component (Lambertian surface)

      IF ( DO_LAMBERTIAN_SURFACE .and. &
             FOURIER_COMPONENT .GT. 0 ) RETURN

!  Skip if non-Lambertian surface

      if ( .not. DO_LAMBERTIAN_SURFACE ) GOTO 878

!  Lambertian surface
!  ==================

!  Integrate Downward streams of particular solutions

      REFL_B = ZERO
      DO J = 1, NSTREAMS
        AXJ = QUAD_STRMWTS(J)
        REFL_B = REFL_B + AXJ*WLOWER(J,1,NLAYERS)
      ENDDO

!  set reflectance

      REFL_B = REFL_B * SURFACE_FACTOR * LAMBERTIAN_ALBEDO
      DO I = 1, NSTREAMS
        R2_BEAM(I,1) = REFL_B
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
!        ENDDO
!      ENDDO
!      pause'r2beam'

!  Return for Lambertian case

      RETURN

!  BRDF surface
!  ============

!  Continuation point

878   CONTINUE

!  Last layer

      NL = NLAYERS

!  Integrate Downward streams of particular solutions

      FACTOR_KERNEL = SURFACE_FACTOR
      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          REFL_B = ZERO
          DO J = 1, NSTREAMS
            REFL_B_S = ZERO
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              REFL_B_S = REFL_B_S + &
                       AXBID_F(I,J,OM) * WLOWER(J,O2,NL)
            ENDDO
            REFL_B = REFL_B + REFL_B_S
          ENDDO
         R2_BEAM(I,O1) = FACTOR_KERNEL * REFL_B
        ENDDO
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
!        ENDDO
!      ENDDO
!      pause'r2beam'

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_BEAM

!

      SUBROUTINE BVP_COLUMN_SETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE, &
        DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS, &
        LAMBERTIAN_ALBEDO, &
        IBEAM, NSTOKES, NSTREAMS, &
        NLAYERS, EMISSIVITY, SURFBB, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        R2_BEAM, DIRECT_BEAM, &
        WUPPER, WLOWER, &
        COL2, SCOL2 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: I,I1,LAY,LAY1,C0,CM,O1,IR,IROW,NSTOKES_FIRST
      DOUBLE PRECISION :: LOCAL_EMISS, FP_SBB

!  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

!  zero column vector

      DO I = 1, NTOTAL
        COL2(I,IBEAM) = ZERO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )
!  -----------------------------

      DO LAY = 2, NLAYERS

        LAY1 = LAY - 1
        C0 = LAY1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = WUPPER(I,O1,LAY) - WLOWER(I,O1,LAY1)
          ENDDO
        ENDDO

      ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------

      LAY = NLAYERS
      C0 = (LAY-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  with non-zero albedo, include integrated downward reflectances

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO

!  no albedo, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO

      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = COL2(CM,IBEAM) + DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Add thermal emission of ground surface (only to final level)
!  ------------------------------------------------------------

!mick fix 2/14/2012 - changed treatment of emissivity in this NLAYERS > 1
!                     section to be consistent with LIDORT

!      IF ( DO_INCLUDE_SURFEMISS ) THEN
!        O1 = 1
!        LOCAL_EMISS = ONE - LAMBERTIAN_ALBEDO
!        FP_SBB = SURFBB
!        IF ( DO_LAMBERTIAN_SURFACE ) THEN
!          DO I = 1, NSTREAMS
!            CM = C0 + NSTOKES*(I-1) + O1
!            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * local_emiss
!          ENDDO
!        ELSE
!          DO I = 1, NSTREAMS
!            CM = C0 + NSTOKES*(I-1) + O1
!!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!!           COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(I,O1)
!            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(O1,I)
!          ENDDO
!        ENDIF
!      ENDIF

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          CM = C0 + NSTOKES*(I-1) + O1
!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!           COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(I,O1)
          COL2(CM,IBEAM) = COL2(CM,IBEAM) + SURFBB * EMISSIVITY(O1,I)
        ENDDO
      ENDIF

!    if ( do_include_surfemiss)then
!       write(*,*)irow,cm,col2(cm,ibeam),FP_SBB * EMISSIVITY(1,1) ; pause
!    endif

!  debug

!      do i = 1, ntotal
!        write(*,*)i,COL2(i,IBEAM)
!      enddo
!      pause'hello'

!  normal completion

      RETURN

!  special case

345   CONTINUE

!  zero column vector

      DO I = 1, NTOTAL
        SCOL2(I,IBEAM) = ZERO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------------------------------------

!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance

      C0 = NSTKS_NSTRMS
      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM)+DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Add thermal emission of ground surface (only to final level)
!  ------------------------------------------------------------

!mick fix 2/14/2012 - simplified loop structure here to make it like the
!                     NLAYERS > 1 thermal section above

!      IF ( DO_INCLUDE_SURFEMISS ) THEN
!        NSTOKES_FIRST = 1
!        DO I = 1, NSTREAMS
!          IR = NSTOKES*(I-1)
!          DO O1 = 1, NSTOKES_FIRST
!            IROW = IR + O1
!            CM = C0 + IROW
!
!!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!!           SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB*EMISSIVITY(I,O1)
!            SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB*EMISSIVITY(O1,I)
!
!          ENDDO
!        ENDDO
!      ENDIF

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          CM = C0 + NSTOKES*(I-1) + O1
!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!           COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(I,O1)
          SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB * EMISSIVITY(O1,I)
        ENDDO
      ENDIF

!
!  finish

      RETURN
      END SUBROUTINE BVP_COLUMN_SETUP

!

      SUBROUTINE BVP_BACKSUB ( &
        FOURIER, IBEAM, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        NLAYERS, &
        NTOTAL, N_SUBDIAG, &
        N_SUPDIAG, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        COL2, SCOL2, &
        LCON, MCON, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::             FOURIER
      INTEGER, INTENT (IN) ::             IBEAM
      LOGICAL, INTENT (IN) ::             DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::             DO_FDTEST
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NTOTAL
      INTEGER, INTENT (IN) ::             N_SUBDIAG
      INTEGER, INTENT (IN) ::             N_SUPDIAG
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::             IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::    SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI
      CHARACTER (LEN=2) :: CB

!  mick fix 7/23/2014. Local arrays defined
      DOUBLE PRECISION  :: LOC_COL2  ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2 ( MAXSTRMSTKS_2, 1 )

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COL2(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              LOC_COL2, MAXTOTAL, INFO )
        COL2(1:NTOTAL,IBEAM) = LOC_COL2(1:NTOTAL,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB, Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,LAY) = COL2(C0+IROW,IBEAM)
            MCON(K,LAY) = COL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,LAY) = COL2(C0+IROW,   IBEAM)
            LCON(K2,LAY) = COL2(C0+IROW_S, IBEAM)
            MCON(K1,LAY) = COL2(C0+IROW1,  IBEAM)
            MCON(K2,LAY) = COL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

!  debug-----------------------------------------------

!        if ( do_debug_write) then
!         if (fourier.eq.1 ) then
!         DO LAY = 1, NLAYERS
!          KO1 = K_REAL(LAY) + 1
!          DO K = 1, K_REAL(LAY)
!           write(50,'(4i3,1p2e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K,LAY), MCON(K,LAY)
!          ENDDO
!          DO K = 1, K_COMPLEX(LAY)
!           K0 = 2*K - 2
!           K1 = KO1 + K0
!           K2 = K1  + 1
!           write(50,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K1,LAY), MCON(K1,LAY),
!     &                LCON(K2,LAY), MCON(K2,LAY)
!          ENDDO
!         ENDDO
!        endif
!        ENDIF

!        DO LAY = 20, NLAYERS
!          DO K = 1, K_REAL(LAY)
!             write(97,'(3i5,1p2e15.7)')FOURIER,LAY,K,
!     &                LCON(K,LAY),MCON(K,LAY)
!            ENDDO
!           ENDDO
!          IF ( FOURIER.EQ.0)STOP

!        write(*,*)'reg',fourier,ibeam,lcon(1,7)

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTKS_NSTRMS_2,1) = SCOL2(1:NSTKS_NSTRMS_2,IBEAM)
        CALL DGETRS ( 'N', NSTKS_NSTRMS_2, 1, &
                      SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                      LOC_SCOL2, MAXSTRMSTKS_2, INFO )
        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = LOC_SCOL2(1:NSTKS_NSTRMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1) in BVP_BACKSUB, Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,LAY) = SCOL2(IROW,IBEAM)
          MCON(K,LAY) = SCOL2(IROW1,IBEAM)
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW = K + K_REAL(LAY)
          IROW1 = IROW + NSTKS_NSTRMS
          IROW_S = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,LAY) = SCOL2(IROW,    IBEAM)
          LCON(K2,LAY) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,LAY) = SCOL2(IROW1,   IBEAM)
          MCON(K2,LAY) = SCOL2(IROW1_S, IBEAM)
        ENDDO

!  debug-----------------------------------------------

!        IF ( DO_DEBUG_WRITE .AND. DO_FDTEST ) THEN
!         DO LAY = 1, NLAYERS
!          KO1 = K_REAL(LAY) + 1
!          DO K = 1, K_REAL(LAY)
!           WRITE(40,'(4I3,2E20.10)') FOURIER,IBEAM,K,LAY, &
!                      LCON(K,LAY), MCON(K,LAY)
!          ENDDO
!          DO K = 1, K_COMPLEX(LAY)
!           K0 = 2*K - 2
!           K1 = KO1 + K0
!           K2 = K1  + 1
!           WRITE(50,'(4I3,4E20.10)') FOURIER,IBEAM,K,LAY, &
!                      LCON(K1,LAY), MCON(K1,LAY), &
!                      LCON(K2,LAY), MCON(K2,LAY)
!          ENDDO
!         ENDDO
!        ENDIF

!      IF ( DO_DEBUG_WRITE ) THEN
!       IF ( DO_FDTEST ) THEN
!        DO N = 1, NLAYERS
!          DO I = 1, NSTKS_NSTRMS
!           write(30,'(3i3,1p2e20.10)')
!     &         FOURIER,IBEAM,I,LCON(I,N),MCON(I,N)
!          ENDDO
!        ENDDO
!       ENDIF
!      ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_BACKSUB

!

      SUBROUTINE BVPTEL_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        NSTOKES, NSTREAMS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, &
        MUELLER_INDEX, K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NSTREAMS_2, NSTKS_NSTRMS_2, T_DELT_EIGEN, &
        DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT, &
        N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, &
        N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        BANDTELMAT2, IPIVOTTEL, &
        R2_HOMP, R2_HOMM, AXBID_F, &
        SMAT2, SIPIVOT, &
        STATUS, MESSAGE, TRACE )

!  Sets up the telescoped boundary value problem.
!    Standard case: Fourier > 0. With surface reflection term

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::            FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      INTEGER, INTENT (IN) ::            NLAYERS
      LOGICAL, INTENT (IN) ::            DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::            DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::            K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::   T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      LOGICAL, INTENT (INOUT) ::         DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::         BVTEL_FOURIER_COMPONENT

      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (OUT) ::           NLAYERS_TEL
      INTEGER, INTENT (OUT) ::           ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (OUT) ::           IPIVOTTEL ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER :: STATUS_SUB

      INTEGER :: NMINTEL ( MAXTOTAL )
      INTEGER :: NMAXTEL ( MAXTOTAL )
      INTEGER :: KALLTEL

!  start code
!  ----------

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVPTEL_MATRIX_INIT_OMP ( &
        FOURIER_COMPONENT, NLAYERS, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT, &
        N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, &
        N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        NMINTEL, NMAXTEL, KALLTEL, &
        BANDTELMAT2 )

!  Do the surface setup if required

!  Necessary to have reflected solutions in the lower boundary
!  layer, even when BVP Telescoping only for Lambertian albedo case.

!  Specialist option 2 only, and only if the lowest active layer
!    is the lowest atmospheric layer. Lambertian only.

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          CALL BVP_SURFACE_SETUP_HOM  ( &
            DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
            SURFACE_FACTOR, &
            NSTOKES, NSTREAMS, &
            NLAYERS, DO_LAMBERTIAN_SURFACE, &
            LAMBERTIAN_ALBEDO, BRDF_F, &
            QUAD_STRMWTS, NSTKS_NSTRMS, &
            MUELLER_INDEX, &
            K_REAL, K_COMPLEX, &
            SOLA_XPOS, SOLB_XNEG, &
            R2_HOMP, R2_HOMM, &
            AXBID_F )
        ENDIF
      ENDIF

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVPTEL_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        R2_HOMP, R2_HOMM, &
        NMINTEL, NMAXTEL, KALLTEL, &
        BANDTELMAT2, SMAT2 )

!  SVD decomposition of compressed boundary values matrix

      CALL BVPTEL_MATRIX_SVD ( &
        NSTKS_NSTRMS_2, N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        NLAYERS_TEL, &
        BANDTELMAT2, SMAT2, &
        IPIVOTTEL, SIPIVOT, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  return

      RETURN
      END SUBROUTINE BVPTEL_MATRIXSETUP_MASTER

!

      SUBROUTINE BVPTEL_MATRIX_INIT ( &
        FOURIER_COMPONENT, NLAYERS, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        DO_BVTEL_INITIAL, &
        N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, &
        N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        NMINTEL, NMAXTEL, KALLTEL, &
        BANDTELMAT2 )

!  Initialise the compressed matrix

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )

      LOGICAL, INTENT (INOUT) ::        DO_BVTEL_INITIAL

      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (OUT) ::          NLAYERS_TEL
      INTEGER, INTENT (OUT) ::          ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (OUT) ::          NMINTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAXTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALLTEL
      DOUBLE PRECISION, INTENT (OUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )

!  Local variables
!  ---------------

      INTEGER :: I, J, NS, N, N3

      INTEGER, SAVE ::          N_BVTELMATRIX_SIZE_SAVE
      INTEGER, SAVE ::          N_BVTELMATRIX_SUPDIAG_SAVE
      INTEGER, SAVE ::          N_BVTELMATRIX_SUBDIAG_SAVE
      INTEGER, SAVE ::          NLAYERS_TEL_SAVE
      INTEGER, SAVE ::          ACTIVE_LAYERS_SAVE ( MAXLAYERS )
      INTEGER, SAVE ::          NMINTEL_SAVE ( MAXTOTAL )
      INTEGER, SAVE ::          NMAXTEL_SAVE ( MAXTOTAL )
      INTEGER, SAVE ::          KALLTEL_SAVE
      DOUBLE PRECISION, SAVE :: BANDTELMAT2_SAVE ( MAXBANDTOTAL, MAXTOTAL )

!  compression function

!      INTEGER  BTELMAT_ROWMASK
!      EXTERNAL BTELMAT_ROWMASK

!  set up
!  ------

      IF ( DO_BVTEL_INITIAL ) THEN

!  Determine active layers in atmosphere

        NS = 0
        ACTIVE_LAYERS = 0
        DO N = 1, NLAYERS
         IF ( DO_LAYER_SCATTERING(FOURIER_COMPONENT,N) ) THEN
          NS = NS + 1
          ACTIVE_LAYERS(NS) = N
         ENDIF
        ENDDO
        NLAYERS_TEL = NS

!  size of Reduced BVTEL matrix

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS * NLAYERS_TEL

!  save

        ACTIVE_LAYERS_SAVE(1:NLAYERS) = ACTIVE_LAYERS(1:NLAYERS)
        NLAYERS_TEL_SAVE              = NLAYERS_TEL
        N_BVTELMATRIX_SIZE_SAVE       = N_BVTELMATRIX_SIZE

!  Exit if only one active layer

        IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

!  Number of sub and super diagonals

        N_BVTELMATRIX_SUPDIAG = 3 * NSTKS_NSTRMS - 1
        N_BVTELMATRIX_SUBDIAG = 3 * NSTKS_NSTRMS - 1

!  save

        N_BVTELMATRIX_SUPDIAG_SAVE = N_BVTELMATRIX_SUPDIAG
        N_BVTELMATRIX_SUBDIAG_SAVE = N_BVTELMATRIX_SUBDIAG

!  compression row indices

        DO J = 1, N_BVTELMATRIX_SUPDIAG + 1
          NMINTEL(J) = 1
        ENDDO
        DO J = N_BVTELMATRIX_SUPDIAG + 2, N_BVTELMATRIX_SIZE
          NMINTEL(J) = J - N_BVTELMATRIX_SUPDIAG
        ENDDO

        DO J = 1, N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG
          NMAXTEL(J) = J + N_BVTELMATRIX_SUBDIAG
        ENDDO
        N3 = N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG + 1
        DO J = N3, N_BVTELMATRIX_SIZE
          NMAXTEL(J) = N_BVTELMATRIX_SIZE
        ENDDO

        KALLTEL = N_BVTELMATRIX_SUBDIAG + N_BVTELMATRIX_SUPDIAG + 1

!  save

        NMINTEL_SAVE(1:N_BVTELMATRIX_SIZE) = NMINTEL(1:N_BVTELMATRIX_SIZE)
        NMAXTEL_SAVE(1:N_BVTELMATRIX_SIZE) = NMAXTEL(1:N_BVTELMATRIX_SIZE)
        KALLTEL_SAVE = KALLTEL

!  Avoid fancy zeroing - adopt kludge
!   Potential Danger point

        DO I = 1, MAXBANDTOTAL
          DO J = 1, N_BVTELMATRIX_SIZE
            BANDTELMAT2(I,J) = ZERO
          ENDDO
        ENDDO

!  save

        BANDTELMAT2_SAVE(:,1:N_BVTELMATRIX_SIZE) = &
          BANDTELMAT2(:,1:N_BVTELMATRIX_SIZE)

!  Control

        GO TO 346

!  control point

 345    CONTINUE

!  single layer setting

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS

 346    CONTINUE

!  reset

        DO_BVTEL_INITIAL = .FALSE.

      ELSE

        ACTIVE_LAYERS(1:NLAYERS) = ACTIVE_LAYERS_SAVE(1:NLAYERS)
        NLAYERS_TEL              = NLAYERS_TEL_SAVE
        N_BVTELMATRIX_SIZE       = N_BVTELMATRIX_SIZE_SAVE

        IF ( NLAYERS_TEL .EQ. 1 ) GO TO 347

        N_BVTELMATRIX_SUPDIAG    = N_BVTELMATRIX_SUPDIAG_SAVE
        N_BVTELMATRIX_SUBDIAG    = N_BVTELMATRIX_SUBDIAG_SAVE
        NMINTEL(1:N_BVTELMATRIX_SIZE) = NMINTEL_SAVE(1:N_BVTELMATRIX_SIZE)
        NMAXTEL(1:N_BVTELMATRIX_SIZE) = NMAXTEL_SAVE(1:N_BVTELMATRIX_SIZE)
        KALLTEL                       = KALLTEL_SAVE

        BANDTELMAT2(:,1:N_BVTELMATRIX_SIZE) = &
          BANDTELMAT2_SAVE(:,1:N_BVTELMATRIX_SIZE)

!  Control

        GO TO 348

!  control point

 347    CONTINUE

!  single layer setting

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS

 348    CONTINUE

!  end initialization

      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_INIT

!

      SUBROUTINE BVPTEL_MATRIX_INIT_OMP ( &
        FOURIER_COMPONENT, NLAYERS, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT, &
        N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, &
        N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        NMINTEL, NMAXTEL, KALLTEL, &
        BANDTELMAT2 )

!  Initialise the compressed matrix

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

      LOGICAL, INTENT (INOUT) ::        DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::        BVTEL_FOURIER_COMPONENT

      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (OUT) ::          NLAYERS_TEL
      INTEGER, INTENT (OUT) ::          ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (OUT) ::          NMINTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAXTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALLTEL
      DOUBLE PRECISION, INTENT (OUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )

!  Local variables
!  ---------------

      INTEGER :: I, J, NS, N, N3, LOCAL_FOURIER_COMPONENT

!  compression function

!      INTEGER  BTELMAT_ROWMASK
!      EXTERNAL BTELMAT_ROWMASK

!mick fix 7/29/2014 - Changed subroutine to define output vars each time to make VLIDORT 
!                     threadsafe.  Removed the saved vars and the original if block with
!                     its else section in this version of the subroutine.  Also added
!                     if block to save the original fourier component used to perform
!                     these setup calculations.

!  Set up
!  ------

!  Save fourier component used the first time subroutine is called

      IF ( DO_BVTEL_INITIAL ) THEN
        LOCAL_FOURIER_COMPONENT = FOURIER_COMPONENT
        BVTEL_FOURIER_COMPONENT = FOURIER_COMPONENT
        DO_BVTEL_INITIAL = .FALSE.
      ELSE
        LOCAL_FOURIER_COMPONENT = BVTEL_FOURIER_COMPONENT
      ENDIF

!  Determine active layers in atmosphere

      NS = 0
      ACTIVE_LAYERS = 0
      DO N = 1, NLAYERS
       IF ( DO_LAYER_SCATTERING(LOCAL_FOURIER_COMPONENT,N) ) THEN
        NS = NS + 1
        ACTIVE_LAYERS(NS) = N
       ENDIF
      ENDDO
      NLAYERS_TEL = NS

!  Size of Reduced BVTEL matrix

      N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS * NLAYERS_TEL

!  Exit if only one active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

!  Number of sub and super diagonals

      N_BVTELMATRIX_SUPDIAG = 3 * NSTKS_NSTRMS - 1
      N_BVTELMATRIX_SUBDIAG = 3 * NSTKS_NSTRMS - 1

!  Compression row indices

      DO J = 1, N_BVTELMATRIX_SUPDIAG + 1
        NMINTEL(J) = 1
      ENDDO
      DO J = N_BVTELMATRIX_SUPDIAG + 2, N_BVTELMATRIX_SIZE
        NMINTEL(J) = J - N_BVTELMATRIX_SUPDIAG
      ENDDO

      DO J = 1, N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG
        NMAXTEL(J) = J + N_BVTELMATRIX_SUBDIAG
      ENDDO
      N3 = N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG + 1
      DO J = N3, N_BVTELMATRIX_SIZE
        NMAXTEL(J) = N_BVTELMATRIX_SIZE
      ENDDO

      KALLTEL = N_BVTELMATRIX_SUBDIAG + N_BVTELMATRIX_SUPDIAG + 1

!  Avoid fancy zeroing - adopt kludge
!   Potential Danger point

      DO I = 1, MAXBANDTOTAL
        DO J = 1, N_BVTELMATRIX_SIZE
          BANDTELMAT2(I,J) = ZERO
        ENDDO
      ENDDO

!  Control

      GO TO 346

!  Control point

 345  CONTINUE

!  Single layer setting

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS

 346  CONTINUE

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_INIT_OMP

!

      SUBROUTINE BVPTEL_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        R2_HOMP, R2_HOMM, &
        NMIN, NMAX, KALL, &
        BANDTELMAT2, SMAT2 )

!  Fills up the matrix directly (compressed or 1-layer)

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      INTEGER, INTENT (IN) ::           NMIN ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           NMAX ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           KALL

      DOUBLE PRECISION, INTENT (INOUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )

      DOUBLE PRECISION, INTENT (OUT)   :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  Local variables
!  ---------------

      INTEGER ::          I, J, I1, N, N1, O1, IR, IROW, NS
      INTEGER ::          EP, EM, EPC, EPS, EMS
      INTEGER ::          CP, CM, CEP, CEM, CEPS, CEMS, CEP1, CEM1
      INTEGER ::          ROWEP, ROWEM, ROWEPS, ROWEMS, ROWCEP, ROWCEM
      INTEGER ::          ROWCEPS, ROWCEMS, ROWCEP1, ROWCEM1
      DOUBLE PRECISION :: XPNET, XMNET, X1R, X2R, X1I, X2I, X1R_H, X1I_H
      INTEGER ::          KO1, K0, K1, K2, KO11
      INTEGER ::          CE_OFFSET, C0

!  compression function

!      INTEGER  BTELMAT_ROWMASK
!      EXTERNAL BTELMAT_ROWMASK

!  Initialization

      EP = 0

!  If Nlayers = 1, go to special case

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

!  top BC for first active layer 1: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS

            ROWEP = BTELMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM = BTELMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)

            BANDTELMAT2(ROWEP,EP)  = &
                         SOLA_XPOS(I,O1,EP,N)
            BANDTELMAT2(ROWEM,EM)  = &
                         SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  complex solutions (REWORKED)

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1

            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS

            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

            ROWEP  = BTELMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM  = BTELMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
            ROWEPS = BTELMAT_ROWMASK(IROW,EPS,NMIN(EPS),NMAX(EPS),KALL)
            ROWEMS = BTELMAT_ROWMASK(IROW,EMS,NMIN(EMS),NMAX(EMS),KALL)

            BANDTELMAT2(ROWEP,EP)   =   X1R
            BANDTELMAT2(ROWEPS,EPS) = - X1I
            BANDTELMAT2(ROWEM,EM)   =   X2R
            BANDTELMAT2(ROWEMS,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  intermediate layer boundaries

      C0 = - NSTKS_NSTRMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        KO1  = K_REAL(N) + 1
        KO11 = K_REAL(N1) + 1
        C0   = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW

!  real solutions for layer above the boundary

            DO EP = 1, K_REAL(N1)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              ROWCEP = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDTELMAT2(ROWCEP,CEP)   = &
                       T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
              BANDTELMAT2(ROWCEM,CEM)   = &
                         SOLB_XNEG(I,O1,EP,N1)
            ENDDO

!  real solutions for layer below the boundary
!   ( Note the change of sign !!! )

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEP1 = CEP + NSTKS_NSTRMS_2
              CEM1 = CEM + NSTKS_NSTRMS_2

              ROWCEP1 = BTELMAT_ROWMASK(CM,CEP1,NMIN(CEP1),NMAX(CEP1),KALL)
              ROWCEM1 = BTELMAT_ROWMASK(CM,CEM1,NMIN(CEM1),NMAX(CEM1),KALL)

              BANDTELMAT2(ROWCEP1,CEP1) = &
                         - SOLA_XPOS(I,O1,EP,N)
              BANDTELMAT2(ROWCEM1,CEM1) = &
                       - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
            ENDDO

!  complex solutions for layer above boundary

            DO EPC = 1, K_COMPLEX(N1)
              K0 = 2*EPC - 2
              K1 = KO11 + K0
              K2 = K1  + 1

              EP  = K_REAL(N1) + EPC
              EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1) &
                   -SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
              X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1) &
                   +SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
              X2R = SOLB_XNEG(I,O1,K1,N1)
              X2I = SOLB_XNEG(I,O1,K2,N1)

              ROWCEP  = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BTELMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BTELMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDTELMAT2(ROWCEP,CEP)   =   X1R
              BANDTELMAT2(ROWCEPS,CEPS) = - X1I
              BANDTELMAT2(ROWCEM,CEM)   =   X2R
              BANDTELMAT2(ROWCEMS,CEMS) = - X2I
            ENDDO

!  complex solutions for layer below boundary
!   ( Note the change of sign !!! )

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              ROWCEP  = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BTELMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BTELMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDTELMAT2(ROWCEP,CEP)   = - X1R
              BANDTELMAT2(ROWCEPS,CEPS) =   X1I
              BANDTELMAT2(ROWCEM,CEM)   = - X2R
              BANDTELMAT2(ROWCEMS,CEMS) =   X2I
            ENDDO

          ENDDO
        ENDDO
      ENDDO

!  bottom BC for Lowest active layer
!   Normally no surface additions, except Specialist # 2

      N = ACTIVE_LAYERS(NLAYERS_TEL)
      KO1 = K_REAL(N) + 1
      C0  = C0 + NSTKS_NSTRMS_2
      CE_OFFSET = C0 - NSTKS_NSTRMS

      IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)

              ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDTELMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDTELMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)

              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)

              ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BTELMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BTELMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDTELMAT2(ROWCEP,CEP)   =   X1R
              BANDTELMAT2(ROWCEPS,CEPS) = - X1I
              BANDTELMAT2(ROWCEM,CEM)   =   X2R
              BANDTELMAT2(ROWCEMS,CEMS) = - X2I
            ENDDO

!  End stokes and streams loops

          ENDDO
        ENDDO

!  No surface albedo reflections

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)

              ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDTELMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDTELMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)

              ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BTELMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BTELMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDTELMAT2(ROWCEP,CEP)   =   X1R
              BANDTELMAT2(ROWCEPS,CEPS) = - X1I
              BANDTELMAT2(ROWCEM,CEM)   =   X2R
              BANDTELMAT2(ROWCEMS,CEMS) = - X2I
           ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal completion

      RETURN

!  special case. Only 1 active layer

345   CONTINUE

!  Set up BVP matrix for the active layer
!  ======================================

!  active layer for telescoped BVP

      N = ACTIVE_LAYERS(1)
      KO1 = K_REAL(N) + 1

!  initialize using the SMAT2 matrix

      DO I = 1, NSTKS_NSTRMS_2
        DO J = 1, NSTKS_NSTRMS_2
          SMAT2(I,J) = ZERO
        ENDDO
      ENDDO

!  top of the layer (downwelling only)
!  -----------------------------------

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
            SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  REWORKED complex solutions

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            SMAT2(IROW,EP)  =   X1R
            SMAT2(IROW,EPS) = - X1I
            SMAT2(IROW,EM)  =   X2R
            SMAT2(IROW,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  bottom BC. Only albedo additions for specialist option 2
!  --------------------------------------------------------

      IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  bottom BC (with No albedo). This is the usual situation
!  -------------------------------------------------------

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  REWORKED Complex solutions
!    Note to self. The second set of solutions seems correct.

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              IF ( N.EQ.NLAYERS ) THEN
                SMAT2(CP,CEP)  =   X1R   ! apparently wrong
                SMAT2(CP,CEPS) = - X1I   ! apparently wrong
                SMAT2(CP,CEM)  =   X2R   ! apparently wrong
                SMAT2(CP,CEMS) = - X2I   ! apparently wrong
              ELSE
                SMAT2(CP,CEP)  =   X1R  ! tested 29 December 2005
                SMAT2(CP,CEPS) = - X1I  ! tested 29 December 2005
                SMAT2(CP,CEM)  =   X2R  ! tested 29 December 2005
                SMAT2(CP,CEMS) = - X2I  ! tested 29 December 2005
              ENDIF
            ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal return and finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_SETUP

!

      SUBROUTINE BVPTEL_MATRIX_SVD ( &
        NSTKS_NSTRMS_2, N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        NLAYERS_TEL, &
        BANDTELMAT2, SMAT2, &
        IPIVOTTEL, SIPIVOT, &
        STATUS, MESSAGE, TRACE )

! TELESCOPED boundary value problem SVD decomposition.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::             NLAYERS_TEL

      DOUBLE PRECISION, INTENT (INOUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

      INTEGER, INTENT (OUT) ::            IPIVOTTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::            SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER           :: INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SVD the BVPTEL matrix: With compression (multilayers)
!  ----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF &
           ( N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SIZE, &
             N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG, &
             BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF &
           (  NSTKS_NSTRMS_2, NSTKS_NSTRMS_2, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF (nlayers_tel=1)call in BVPTEL_MATRIX_SVD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_SVD

!

      SUBROUTINE BVPTEL_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        FOURIER_COMPONENT, IBEAM, &
        SURFACE_FACTOR, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        NSTOKES, NSTREAMS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        QUAD_STRMWTS, MUELLER_INDEX, WLOWER, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        DIRECT_BEAM, WUPPER, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        IPIVOT, SMAT2, SIPIVOT, &
        N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        BANDTELMAT2, IPIVOTTEL, &
        AXBID_F, &
        COLTEL2, SCOL2, &
        R2_BEAM, LCON, MCON, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           IBEAM
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::           SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SUBDIAG
      DOUBLE PRECISION, INTENT (IN) ::  BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::           IPIVOTTEL ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!mick fix 7/29/2014 - placed COLTEL2 & SCOL2 in call to make VLIDORT threadsafe
      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

      DOUBLE PRECISION, INTENT (INOUT) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

      INTEGER ::          STATUS_SUB

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  This is suitable for Lambertian surfaces only.
!  --Additional setups for the albedo layer.
!    Only required for the Specialise Option 2 case

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          CALL BVP_SURFACE_SETUP_BEAM ( &
            DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
            SURFACE_FACTOR, &
            NSTOKES, NSTREAMS, &
            NLAYERS, DO_LAMBERTIAN_SURFACE, &
            LAMBERTIAN_ALBEDO, &
            QUAD_STRMWTS, MUELLER_INDEX, &
            WLOWER, AXBID_F, &
            R2_BEAM )
        ENDIF
      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        IBEAM, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        R2_BEAM, DIRECT_BEAM, &
        WUPPER, WLOWER, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        SCOL2, COLTEL2 )

!  --Solve the boundary problem for this Fourier component (back substit

      CALL BVPTEL_BACKSUB ( &
        IBEAM, NSTOKES, &
        NSTREAMS, NLAYERS, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, WLOWER, &
        IPIVOT, &
        SMAT2, SIPIVOT, &
        N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        BANDTELMAT2, IPIVOTTEL, &
        COLTEL2, &
        SCOL2, LCON, MCON, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  return

      RETURN
      END SUBROUTINE BVPTEL_SOLUTION_MASTER

!

      SUBROUTINE BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        IBEAM, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_SPECIALIST_OPTION_2, &
        NSTREAMS_2, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        R2_BEAM, DIRECT_BEAM, &
        WUPPER, WLOWER, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, &
        ACTIVE_LAYERS, &
        SCOL2, COLTEL2 )

!  Sets up the telescoped boundary value problem, RHS vector
!    Standard case: Fourier > 0.
!         Suitable for Lambertian surfaces.
!    Surface reflection only appears in the Specialist Option 2 case.

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_2
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER ::         I, I1, N, N1, NS, C0, CM
      INTEGER ::         IR, O1, IROW

!  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

!  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        COLTEL2(I,IBEAM) = ZERO
      ENDDO

!  Upper boundary for first active layer: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COLTEL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
        ENDDO
      ENDDO

!  intermediate layer boundaries

      C0 = - NSTKS_NSTRMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        C0 = C0 + NSTKS_NSTRMS_2
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COLTEL2(CM,IBEAM) = WUPPER(I,O1,N) - WLOWER(I,O1,N1)
          ENDDO
        ENDDO
      ENDDO

!  Lower boundary for last active layer NLAYERS_TEL:
!  No albedo, as this is FOURIER > 0

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)
      C0 = C0 + NSTKS_NSTRMS_2
      IF ( DO_INCLUDE_SURFACE .AND. DO_SPECIALIST_OPTION_2 ) THEN
        IF ( N.EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COLTEL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                COLTEL2(CM,IBEAM) = COLTEL2(CM,IBEAM) &
                        + DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COLTEL2(CM,IBEAM) = - WLOWER(I1,O1,N)
          ENDDO
        ENDDO
      ENDIF

!  Normal return

      RETURN

!  Continuity point for single layer setup

 345  CONTINUE

!  active layer for telescoped BVP

      N = ACTIVE_LAYERS(1)

!  initialize using the SCOL2 matrix

      DO I = 1, NSTKS_NSTRMS_2
        SCOL2(I,IBEAM) = ZERO
      ENDDO

!  Set up BVP column for the active layer
!  ======================================

!  Upper boundary for layer (downwelling only)
!  -------------------------------------------

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
        ENDDO
      ENDDO

!  lower boundary for layer.
!    Surface only required for the Specialist  2 option

      C0 = NSTKS_NSTRMS

      IF ( DO_INCLUDE_SURFACE .AND. DO_SPECIALIST_OPTION_2 ) THEN
        IF ( N.EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) &
                        + DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N)
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_COLUMN_SETUP

!

      SUBROUTINE BVPTEL_BACKSUB ( &
        IBEAM, NSTOKES, &
        NSTREAMS, NLAYERS, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, WLOWER, &
        IPIVOT, &
        SMAT2, SIPIVOT, &
        N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
        NLAYERS_TEL, ACTIVE_LAYERS, &
        BANDTELMAT2, IPIVOTTEL, &
        COLTEL2, &
        SCOL2, LCON, MCON, &
        STATUS, MESSAGE, TRACE )

!  Solves the telescoped boundary value problem.
!    Standard case: Fourier > 0. No surface reflection term
!         Suitable for Lambertian surfaces.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::             IBEAM
      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::    T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (IN) ::             IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::    SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::             NLAYERS_TEL
      INTEGER, INTENT (IN) ::             ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::             IPIVOTTEL ( MAXTOTAL )

      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           K, KO1, K0, K1, K2
      INTEGER ::           I, I1, N, NS, N1, NAC, O1, IC, ICOW, INFO
      INTEGER ::           C0, IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION ::  SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION ::  SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION ::  LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI
      CHARACTER (LEN=3) :: CI
      CHARACTER (LEN=2) :: CB

      DOUBLE PRECISION  :: LOC_COLTEL2 ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2   ( MAXSTRMSTKS_2, 1 )

!  start code
!  ----------

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Back-substitution for multi-layer BVP TEL
!  =========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution using RHS column vector COLTEL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COLTEL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1) = COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM)
        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, &
              N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG, 1, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              LOC_COLTEL2, MAXTOTAL, INFO )
        COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM) = LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = &
              'DGBTRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set integration constants for active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTKS_NSTRMS_2
          KO1 = K_REAL(N) + 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,N) = COLTEL2(C0+IROW,IBEAM)
            MCON(K,N) = COLTEL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,N) = COLTEL2(C0+IROW,   IBEAM)
            LCON(K2,N) = COLTEL2(C0+IROW_S, IBEAM)
            MCON(K1,N) = COLTEL2(C0+IROW1,  IBEAM)
            MCON(K2,N) = COLTEL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

!        write(*,*)'tel',fourier,ibeam,lcon(1,7)

!  Solve the boundary problem: Single Layer only
!  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  Active layer

        N = ACTIVE_LAYERS(1)
        KO1 = K_REAL(N) + 1

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTKS_NSTRMS_2,1) = SCOL2(1:NSTKS_NSTRMS_2,IBEAM)
        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, 1, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              LOC_SCOL2, MAXSTRMSTKS_2, INFO )
        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = LOC_SCOL2(1:NSTKS_NSTRMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = &
              'DGETRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set real constants from the solution vector

        DO K = 1, K_REAL(N)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,N) = SCOL2(IROW,IBEAM)
          MCON(K,N) = SCOL2(IROW1,IBEAM)
        ENDDO

!  set complex constants from the solution vector

        DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,N) = SCOL2(IROW,    IBEAM)
          LCON(K2,N) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,N) = SCOL2(IROW1,   IBEAM)
          MCON(K2,N) = SCOL2(IROW1_S, IBEAM)
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Set integration constants for non-active layers
!  ===============================================

!  Transmittance layers above active layer
!  ---------------------------------------

!   -- LCON values are zero (no downwelling radiation)
!   -- MCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!    .... Require solutions at top of active layer

      NAC = ACTIVE_LAYERS(1)
      KO1 = K_REAL(NAC) + 1

      IF ( NAC .GT. 1 ) THEN

        N1 = NAC - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            ICOW = IC + O1

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
             LXR = LCON(K,NAC) * SOLA_XPOS(I1,O1,K,NAC)
             MXR = MCON(K,NAC) * SOLB_XNEG(I1,O1,K,NAC)
             HOM1 = LXR
             HOM2 = MXR * T_DELT_EIGEN(K,NAC)
             SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAC) * SOLA_XPOS(I1,O1,K1,NAC) - &
                       LCON(K2,NAC) * SOLA_XPOS(I1,O1,K2,NAC)
              MXR_CR = MCON(K1,NAC) * SOLB_XNEG(I1,O1,K1,NAC) - &
                       MCON(K2,NAC) * SOLB_XNEG(I1,O1,K2,NAC)
              MXR_CI = MCON(K1,NAC) * SOLB_XNEG(I1,O1,K2,NAC) + &
                       MCON(K2,NAC) * SOLB_XNEG(I1,O1,K1,NAC)
              HOM1CR = LXR_CR
              HOM2CR = MXR_CR*T_DELT_EIGEN(K1,NAC) &
                     - MXR_CI*T_DELT_EIGEN(K2,NAC)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,NAC)
            MCON(ICOW,N1) = SPAR + SHOM
            LCON(ICOW,N1) = ZERO

          ENDDO
        ENDDO

      ENDIF

!  other layers to top, just propagate

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            LCON(IROW,N) = ZERO
            MCON(IROW,N) = T_DELT_DISORDS(I,N1) * MCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

!  Transmittance layers below active layer
!  ---------------------------------------

!  layer immediately below active layer
!    .... Require solutions at bottom of active layer

      NAC = ACTIVE_LAYERS(NLAYERS_TEL)
      KO1 = K_REAL(NAC) + 1
      IF ( NAC .LT. NLAYERS ) THEN
        N1 = NAC + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
              LXR = LCON(K,NAC) * SOLA_XPOS(I,O1,K,NAC)
              MXR = MCON(K,NAC) * SOLB_XNEG(I,O1,K,NAC)
              HOM1 = LXR * T_DELT_EIGEN(K,NAC)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAC) * SOLA_XPOS(I,O1,K1,NAC) - &
                       LCON(K2,NAC) * SOLA_XPOS(I,O1,K2,NAC)
              LXR_CI = LCON(K1,NAC) * SOLA_XPOS(I,O1,K2,NAC) + &
                       LCON(K2,NAC) * SOLA_XPOS(I,O1,K1,NAC)
              MXR_CR = MCON(K1,NAC) * SOLB_XNEG(I,O1,K1,NAC) - &
                       MCON(K2,NAC) * SOLB_XNEG(I,O1,K2,NAC)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,NAC) &
                      -LXR_CI*T_DELT_EIGEN(K2,NAC)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part and add particular solution
!    Note to self. The  MINUS sign appears to be correct !  WHY?

!            SHOM = SHOM_R + SHOM_CR   ! apparently wrong
            SHOM = SHOM_R - SHOM_CR
            SPAR = WLOWER(I,O1,NAC)
            LCON(IROW,N1) = SPAR + SHOM
            MCON(IROW,N1) = ZERO

          ENDDO
        ENDDO

      ENDIF

!  other layers to bottom

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            MCON(IROW,N) = ZERO
            LCON(IROW,N) = T_DELT_DISORDS(I,N1) * LCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

!  debug

!        IF ( FOURIER .EQ. 3 ) THEN
!        DO N = 1, NLAYERS
!          DO K = 1, K_REAL(N)
!             WRITE(98,'(3i5,1p2e15.7)')FOURIER,N,K,
!     &        LCON(K,N), MCON(K,N)
!            ENDDO
!        ENDDO
!        ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_BACKSUB

!

      INTEGER FUNCTION BMAT_ROWMASK ( &
        I, J, NMINJ, NMAXJ, KALL )

!  Integer function for Band-matrix compression
!   Replaces Array BMAT_ROWMASK(I,J) which was too large.
!   Replacement by J. Kujanpaa, FMI, August 2005.

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: I
      INTEGER, INTENT (IN) :: J
      INTEGER, INTENT (IN) :: NMINJ
      INTEGER, INTENT (IN) :: NMAXJ
      INTEGER, INTENT (IN) :: KALL

!  Assign the function (Keep the hard-wired Stop for now)
!    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINJ) .AND. (I.LE.NMAXJ) ) THEN
         BMAT_ROWMASK = KALL + I - J
      ELSE
         BMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BMAT_ROWMASK:',I,J
         STOP
      ENDIF

!  Finish

      RETURN
      END FUNCTION BMAT_ROWMASK

!

      INTEGER FUNCTION BTELMAT_ROWMASK ( &
        I, J, NMINJ, NMAXJ, KALL )

!  Integer function for Band-matrix compression
!   Replaces Array BTELMAT_ROWMASK(I,J) which was too large.
!   Replacement by J. Kujanpaa, FMI, August 2005.

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: I
      INTEGER, INTENT (IN) :: J
      INTEGER, INTENT (IN) :: NMINJ
      INTEGER, INTENT (IN) :: NMAXJ
      INTEGER, INTENT (IN) :: KALL

!  Assign the function (Keep the hard-wired Stop for now)
!    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINJ) .AND. (I.LE.NMAXJ) ) THEN
         BTELMAT_ROWMASK = KALL + I - J
      ELSE
         BTELMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BTELMAT_ROWMASK:',I,J
         STOP
      ENDIF

!  Finish

      RETURN
      END FUNCTION BTELMAT_ROWMASK

      END MODULE vlidort_bvproblem

