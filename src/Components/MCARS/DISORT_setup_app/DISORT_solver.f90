module DISORT_solver

	implicit none
	
	private
	
	integer, parameter :: MAXIT = 1000, MXSQT = 1000
	real*8, parameter :: one= 1.0d0, two = 2.0d0, zero = 0.0d0
	logical, parameter :: DELTAM = .true., CORINT = .false.

	real*8 :: TOL
	INTEGER, parameter ::  MXCLY = 26, MXULV = 27, MXCMU = 64, MXUMU = 5, & ! MXCMU was 16 and 32
               MXPHI = 5, MI = MXCMU / 2, MI9M2 = 9*MI - 2, &
               NNLYRI = MXCMU*MXCLY      

	real, parameter :: approx_limit = 0.3, approx_limit_tan = 0.2
	
	real :: TINY, HUGE, PI, PLKAVG_SIGDPI, PLKAVG_CONC, DITHER, PIMYPI
	real :: RPD
	real :: SQT(MXSQT)
	REAL :: CMU_S(MXCMU),CWT_S(MXCMU), PHIRAD( MXPHI )
		
	integer :: NN

	
	DOUBLE PRECISION :: COSPHI_S(0:MXCMU-1,1:MXPHI)
	

	public :: DISORT_DP_GC, init_solver

contains

	subroutine init_solver
	
  		 USE RDI1MACH_f90, ONLY : R1MACH, D1MACH
		
		 real, parameter :: SIGMA = 5.67032E-8 
! number of streams
		 integer, parameter :: NSTR = 64


		 integer :: NS, IQ


		
         TINY   = D1MACH( 1 )
         HUGE   = D1MACH( 2 )
		 TOL = D1MACH(4)
		 
		 PI = acos(-1.0)
		 
         PLKAVG_SIGDPI = SIGMA / PI
         PLKAVG_CONC   = 15. / PI**4
		 
         DO NS = 1, MXSQT
            SQT( NS ) = SQRT( FLOAT( NS ) )
		 end do
		
         IF( DITHER.LT.1.E-10 ) DITHER = 10.*DITHER

         RPD  = PI / 180.0

		 
		 PIMYPI = PI * 1.0d0
		 
         
		 NN   = NSTR / 2
		CMU_S(:) = 0.0
		CWT_S(:) = 0.0

      CALL QGAUSN( NN, CMU_S, CWT_S )
!c                                  ** Downward (neg) angles and weights
      DO IQ = 1, NN
         CMU_S( IQ + NN ) = -CMU_S( IQ )
         CWT_S( IQ + NN ) = CWT_S( IQ )
      ENDDO
	
	
	end subroutine init_solver


      SUBROUTINE DISORT_DP_GC( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
                WVNMLO,WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
                        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
                        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, &
                        PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY, &
                        MAXULV, MAXUMU, MAXPHI, MAXMOM, UU,UUS)
     
		
	 
    	USE RDI1MACH_f90, ONLY : R1MACH, D1MACH

		implicit none

!c
!c               I N T E R N A L    V A R I A B L E S
!c
!c   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
!c                     of Eqs. SS(12), STWJ(8E), STWL(23f)
!c                     (used only in SOLEIG)
!c
!c   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
!c                     of Eqs. SS(12), STWJ(8E), STWL(23f)
!c                     (used only in SOLEIG)
!c
!c   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
!c                     (see each subroutine for definition)
!c
!c   B()               Right-hand side vector of Eq. SC(5) going into
!c                     SOLVE0,1;  returns as solution vector
!c                     vector  L, the constants of integration
!c
!c   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
!c                     given azimuthal component.  First index always
!c                     refers to a computational angle.  Second index:
!c                     if zero, refers to incident beam angle UMU0;
!c                     if non-zero, refers to a computational angle.
!c
!c   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
!c                     tational angles.
!c
!c   BPLANK            Intensity emitted from bottom boundary
!c
!c   CBAND()           Matrix of left-hand side of the linear system
!c                     Eq. SC(5), scaled by Eq. SC(12);  in banded
!c                     form required by LINPACK solution routines
!c
!c   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)
!c
!c   CMU(IQ)           Computational polar angles (Gaussian)
!c
!c   CWT(IQ)           Quadrature weights corresponding to CMU
!c
!c   CORINT            When set TRUE, correct intensities for
!c                     delta-scaling effects (see Nakajima and Tanaka,
!c                     1988). When FALSE, intensities are not corrected.
!c                     In general, CORINT should be set true when beam
!c                     source is present (FBEAM is not zero) and DELTAM
!c                     is TRUE in a problem including scattering.
!c                     However, execution is faster when CORINT is FALSE,
!c                     and intensities outside the aureole may still be
!c                     accurate enough.  When CORINT is TRUE, it is
!c                     important to have a sufficiently high order of
!c                     Legendre approximation of the phase function. This
!c                     is because the intensities are corrected by
!c                     calculating the single-scattered radiation, for
!c                     which an adequate representation of the phase
!c                     function is crucial.  In case of a low order
!c                     Legendre approximation of an otherwise highly
!c                     anisotropic phase function, the intensities might
!c                     actually be more accurate when CORINT is FALSE.
!c                     When only fluxes are calculated (ONLYFL is TRUE),
!c                     or there is no beam source (FBEAM=0.0), or there
!c                     is no scattering (SSALB=0.0 for all layers) CORINT
!c                     is set FALSE by the code.
!c
!c   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
!c                     is the number of the Fourier component in the
!c                     azimuth cosine expansion
!c
!c   DELTAM            TRUE,  use delta-M method ( see Wiscombe, 1977 );
!c                     FALSE, do not use delta-M method. In general, for
!c                     a given number of streams, intensities and
!c                     fluxes will be more accurate for phase functions
!c                     with a large forward peak if DELTAM is set true.
!c                     Intensities close to the forward scattering
!c                     direction are often less accurate, however, when
!c                     the delta-M method is applied. The intensity
!c                     correction of Nakajima and Tanaka is used to
!c                     improve the accuracy of the intensities.
!c
!c   DITHER            Small quantity subtracted from single-scattering
!c                     albedos of unity, in order to avoid using special
!c                     case formulas;  prevents an eigenvalue of exactly
!c                     zero from occurring, which would cause an
!c                     immediate overflow
!c
!c   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
!c                     if DELTAM = TRUE, otherwise equal to DTAUC)
!c
!c   EMU(IU)           Bottom-boundary directional emissivity at user
!c                     angles.
!c
!c   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)
!c
!c   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
!c                     SOLEIG; stored permanently in  GC
!c
!c   EXPBEA(LC)        Transmission of direct beam in delta-M optical
!c                     depth coordinates
!c
!c   FLYR(LC)          Separated fraction in delta-M method
!c
!c   GL(K,LC)          Phase function Legendre polynomial expansion
!c                     coefficients, calculated from PMOM by
!c                     including single-scattering albedo, factor
!c                     2K+1, and (if DELTAM=TRUE) the delta-M
!c                     scaling
!c
!c   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
!c                     g  in Eq. SC(1)
!c
!c   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
!c                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
!c                       G without the L factor )
!c
!c   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
!c                     routines
!c
!c   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)
!c
!c   KCONV             Counter in azimuth convergence test
!c
!c   LAYRU(LU)         Computational layer in which user output level
!c                     UTAU(LU) is located
!c
!c   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
!c                     obtained by solving scaled version of Eq. SC(5)
!c
!c   LYRCUT            TRUE, radiation is assumed zero below layer
!c                     NCUT because of almost complete absorption
!c
!c   NAZ               Number of azimuthal components considered
!c
!c   NCUT              Computational layer number in which absorption
!c                     optical depth first exceeds ABSCUT
!c
!c   OPRIM(LC)         Single scattering albedo after delta-M scaling
!c
!c   PASS1             TRUE on first entry, FALSE thereafter
!c
!c   PKAG(0:LC)        Integrated Planck function for internal emission
!c
!c   PRNTU0(L)         logical flag to trigger printing of azimuthally-
!c                     averaged intensities:
!c                       L    quantities printed
!c                      --    ------------------
!c                       1    azimuthally-averaged intensities at user
!c                               levels and computational polar angles
!c                       2    azimuthally-averaged intensities at user
!c                               levels and user polar angles
!c
!c   PSI0(IQ)          Sum just after square bracket in  Eq. SD(9)
!c
!c   PSI1(IQ)          Sum in  Eq. STWL(31d)
!c
!c   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
!c                     given azimuthal component.  First index always
!c                     refers to a user angle.  Second index:
!c                     if zero, refers to incident beam angle UMU0;
!c                     if non-zero, refers to a computational angle.
!c
!c   SQT(k)            Square root of k (used only in LEPOLY for
!c                     computing associated Legendre polynomials)
!c
!c   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)
!c
!c   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
!c                     DELTAM = TRUE, otherwise equal to TAUC)
!c
!c   TPLANK            Intensity emitted from top boundary
!c
!c   UUM(IU,LU)        Expansion coefficients when the intensity
!c                     (u-super-M) is expanded in Fourier cosine series
!c                     in azimuth angle
!c
!c   U0C(IQ,LU)        Azimuthally-averaged intensity at quadrature
!c                     angle
!c
!c   U0U(IU,LU)        If ONLYFL = FALSE, azimuthally-averaged intensity
!c                     at user angles and user levels
!c
!c                     If ONLYFL = TRUE and MAXUMU.GE.NSTR,
!c                     azimuthally-averaged intensity at computational
!c                     (Gaussian quadrature) angles and user levels;
!c                     the corresponding quadrature angle cosines are
!c                     returned in UMU.  If MAXUMU.LT.NSTR, U0U will be
!c                     zeroed, and UMU, NUMU will not be set.
!c
!c   UTAUPR(LU)        Optical depths of user output levels in delta-M
!c                     coordinates;  equal to  UTAU(LU) if no delta-M
!c
!c   WK()              scratch array
!c
!c   XR0(LC)           X-sub-zero in expansion of thermal source func-
!c                     tion preceding Eq. SS(14)(has no mu-dependence);
!c                     b-sub-zero in Eq. STWL(24d)
!c
!c   XR1(LC)           X-sub-one in expansion of thermal source func-
!c                     tion; see  Eqs. SS(14-16); b-sub-one in STWL(24d)
!c
!c   YLM0(L)           Normalized associated Legendre polynomial
!c                     of subscript L at the beam angle (not saved
!c                     as function of superscipt M)
!c
!c   YLMC(L,IQ)        Normalized associated Legendre polynomial
!c                     of subscript L at the computational angles
!c                     (not saved as function of superscipt M)
!c
!c   YLMU(L,IU)        Normalized associated Legendre polynomial
!c                     of subscript L at the user angles
!c                     (not saved as function of superscipt M)
!c
!c   Z()               scratch array used in SOLVE0, ALBTRN to solve
!c                     a linear system for the constants of integration
!c
!c   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)
!c
!c   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
!c                     angles from an equation derived from SS(16)
!c
!c   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)
!c
!c   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
!c                     angles from an equation derived from SS(16)
!c
!c   ZBEAM(IU,LC)      Particular solution for beam source
!c
!c   ZJ(IQ)            Right-hand side vector  X-sub-zero in
!c                     Eq. SS(19), also the solution vector
!c                     Z-sub-zero after solving that system
!c
!c   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ
!c
!c   ZPLK0(IQ,LC)      Permanent storage for the thermal source
!c                     vectors  Z0  obtained by solving  Eq. SS(16)
!c
!c   ZPLK1(IQ,LC)      Permanent storage for the thermal source
!c                     vectors  Z1  obtained by solving  Eq. SS(16)
!c
!c +-------------------------------------------------------------------+
!c
!c  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):
!c
!c       MXCLY  = Max no. of computational layers
!c       MXULV  = Max no. of output levels
!c       MXCMU  = Max no. of computation polar angles
!c       MXUMU  = Max no. of output polar angles
!c       MXPHI  = Max no. of output azimuthal angles
!c       MXSQT  = Max no. of square roots of integers (for LEPOLY)

!c +-------------------------------------------------------------------+

!c     .. Parameters ..

!c     ..
!c     .. Scalar Arguments ..

      CHARACTER, intent(in) :: HEADER*127
      LOGICAL, intent(in) ::   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER, intent(in) ::   IBCND, MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, NLYR, NMOM, NPHI, NSTR, NTAU, NUMU
      REAL, intent(in) ::      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,UMU0, WVNMHI, WVNMLO
     
 !c     ..
!c     .. Array Arguments ..

      LOGICAL, intent(inout) ::   PRNT( 5 )
      REAL, intent(inout) ::       DTAUC( : ), &
                PHI( : ), PMOM( 0:, : ), &
                SSALB( : ), &
               TEMPER( 0: ),  &
               UMU( : ), UTAU( : ), &
               UU( :,:,: ),UUS(:,:,:)
          
!c     ..
!c     .. Local Scalars ..

      LOGICAL   LYRCUT 
      INTEGER   IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL, NCOS, NCUT, NS
	  integer :: iazm, ier, istr
      REAL      AZERR, AZTERM, BPLANK, COSPHI, DELM0, DUM, SGN, TPLANK
      REAL       ANGCOS(1)

!c     ..
!c     .. Local Arrays ..

      LOGICAL   PRNTU0( 2 )
      INTEGER   IPVT( NNLYRI ), LAYRU( MXULV )

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MXCMU, MXCMU ), &
               B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ), &
               CBAND( MI9M2, NNLYRI ), CC( MXCMU, MXCMU ), &
	           CMU( MXCMU ), CWT( MXCMU ), DTAUCP( MXCLY ), &
               EMU( MXUMU ), EVAL( MI ), EVECC( MXCMU, MXCMU ), &
               EXPBEA( 0:MXCLY ), FLDIR( MXULV ), FLDN( MXULV ), &
               FLYR( MXCLY ), GC( MXCMU, MXCMU, MXCLY ), &
               GL( 0:MXCMU, MXCLY ), GU( MXUMU, MXCMU, MXCLY ), &
               KK( MXCMU, MXCLY ), LL( MXCMU, MXCLY ), OPRIM( MXCLY ), &
               PHASA( MXCLY ), PHAST( MXCLY ), PHASM( MXCLY ), &
               PKAG( 0:MXCLY ), PSI0( MXCMU ), &
               PSI1( MXCMU ), RMU( MXUMU, 0:MI ), &
               TAUC( 0:MXCLY ), TAUCPR( 0:MXCLY ), U0C( MXCMU, MXULV ), &
               U0U( MXUMU, MXULV ), UTAUPR( MXULV ), &
               UUM( MXUMU, MXULV ), WK( MXCMU ), XR0( MXCLY ), &
               XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ), &
               YLMU( 0:MXCMU, MXUMU ), Z( NNLYRI ), Z0( MXCMU ), &
               Z0U( MXUMU, MXCLY ), Z1( MXCMU ), Z1U( MXUMU, MXCLY ), &
               ZBEAM( MXUMU, MXCLY )
      REAL      ZJ( MXCMU ), ZPLK0( MXCMU, MXCLY ),ZPLK1( MXCMU, MXCLY ), ZZ( MXCMU, MXCLY )

		real :: mazim_phirad

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ), WKD( MXCMU )
     
!c     ..


		integer :: start_time, end_time, crate, cmax

   
		TAUC = 0.0
	

      DO 30 LC = 1, NLYR

         IF( SSALB( LC ).EQ.1.0 ) SSALB( LC ) = 1.0 - DITHER
         TAUC( LC ) = TAUC( LC - 1 ) + DTAUC( LC )
!c       		print*,SSALB(LC)

   30 CONTINUE


!c                                ** Check input dimensions and variables

!c                                 ** Zero internal and output arrays


			iazm = 1
!c         DO iazm = 1,NPHI         
         	PHIRAD(iazm) = (PHI(iazm) - PHI0) * PIMYPI/180.d0
         	COSPHI_S(0,iazm) = 1.d0
         	if (PHIRAD(iazm) < approx_limit) then 
         		COSPHI_S(1,iazm) = 1. - PHIRAD(iazm)*PHIRAD(iazm)/2.
				else
	         	COSPHI_S(1,iazm) = COS(PHIRAD(iazm))
				endif
         
         DO istr = 2,NSTR-1
       		COSPHI_S(istr,iazm) = 2.d0 * COSPHI_S(istr-1,iazm) * COSPHI_S(1,iazm) - COSPHI_S(istr-2,iazm) 
         ENDDO
!c         ENDDO


!c                                 ** Perform various setup operations

      CALL SETDIS( CMU_S, DELTAM, DTAUC, DTAUCP, &
             EXPBEA, FBEAM, FLYR, &
                  GL, IBCND, LAYRU, LYRCUT, MAXMOM, MAXUMU, MXCMU, &
                  NCUT, NLYR, NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, &
                  CORINT, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU, &
                  UTAUPR, UMU, UMU0, USRTAU, USRANG )


!c                              ** Handle special case for getting albedo
!c                              ** and transmissivity of medium for many
!c                              ** beam angles at once

!c                                   ** Calculate Planck functions
      IF( .NOT.PLANK ) THEN

         BPLANK = 0.0
         TPLANK = 0.0
		  PKAG = 0.0

      ELSE

         TPLANK = 0.0
         BPLANK =       PLKAVG( WVNMLO, WVNMHI, BTEMP )

         DO 40 LEV = 0, NLYR
            PKAG( LEV ) = PLKAVG( WVNMLO, WVNMHI, TEMPER( LEV ) )
   40    CONTINUE

      END IF
		
!c ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
!c           (EQ STWJ 5, STWL 6)

      KCONV  = 0
      NAZ    = NSTR - 1
!c                                    ** Azimuth-independent case

      IF( FBEAM.EQ.0.0 .OR. ABS(1.-UMU0).LT.1.E-5 .OR. ( ABS(1.-UMU(1)).LT.1.E-5 ) ) NAZ = 0




      DO 180 MAZIM = 0, NAZ

         IF( MAZIM.EQ.0 ) DELM0  = 1.0
         IF( MAZIM.GT.0 ) DELM0  = 0.0

!c                             ** Get normalized associated Legendre
!c                             ** polynomials for
!c                             ** (a) incident beam angle cosine
!c                             ** (b) computational and user polar angle
!c                             **     cosines


         IF( FBEAM.GT.0.0 ) THEN

            NCOS   = 1
            ANGCOS = -UMU0

          CALL LEPOLY( NCOS, MAZIM, MXCMU, NSTR-1, ANGCOS, SQT, YLM0 )

			
         END IF


		 CALL LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, SQT, YLMU )

         CALL LEPOLY( NN, MAZIM, MXCMU, NSTR-1, CMU_S, SQT, YLMC )


!c                       ** Get normalized associated Legendre polys.
!c                       ** with negative arguments from those with
!c                       ** positive arguments; Dave/Armstrong Eq. (15),
!c                       ** STWL(59)


         SGN  = -1.0

         DO 70 L = MAZIM, NSTR - 1

            SGN  = -SGN

            DO 60 IQ = NN + 1, NSTR
               YLMC( L, IQ ) = SGN*YLMC( L, IQ - NN )
   60       CONTINUE

   70    CONTINUE
!c                                 ** Specify users bottom reflectivity
!c                                 ** and emissivity properties
         IF( .NOT.LYRCUT )THEN		 
            CALL SURFAC( ALBEDO, LAMBER, MI, MAZIM, MXUMU, NN,  BDR, EMU, BEM, RMU )

		ENDIF



!c ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO 80 LC = 1, NCUT

!c                      ** Solve eigenfunction problem in Eq. STWJ(8B),
!c                      ** STWL(23f); return eigenvalues and eigenvectors

            CALL SOLEIG( AMB, APB, ARRAY, CMU_S, CWT_S, GL( 0,LC ), MI, &
                        MAZIM, MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL, &
                        KK( 1,LC ), GC( 1,1,LC ), AAD, EVECCD, EVALD, &
                        WKD )
!c                                  ** Calculate particular solutions of
!c                                  ** Eq. SS(18), STWL(24a) for incident
!c                                  ** beam source
            IF( FBEAM.GT.0.0 ) &
              CALL UPBEAM( ARRAY, CC, CMU_S, DELM0, FBEAM, GL( 0,LC ), &
                            IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK, &
                            YLM0, YLMC, ZJ, ZZ( 1,LC ) )


!c                              ** Calculate particular solutions of Eq.
!c                              ** SS(15), STWL(25) for thermal emission
!c                              ** source
!c
			
            IF( PLANK .AND. MAZIM.EQ.0 ) THEN

               XR1( LC ) = 0.0

               IF( DTAUCP( LC ).GT.0.0 ) XR1( LC ) = ( PKAG( LC ) - PKAG( LC-1 ) ) / DTAUCP( LC )

               XR0( LC ) = PKAG( LC-1 ) - XR1( LC )*TAUCPR( LC-1 )

               CALL UPISOT( ARRAY, CC, CMU_S, IPVT, MXCMU, NN, NSTR, &
                          OPRIM( LC ), WK, XR0( LC ), XR1( LC ), &
                           Z0, Z1, ZPLK0( 1,LC ), ZPLK1( 1,LC ) )

			else
				Z0 = 0.
				Z1 = 0.
            END IF

			
			
			
			


            IF( .NOT.ONLYFL .AND. USRANG ) THEN

!c                                            ** Interpolate eigenvectors
!c                                            ** to user angles

            CALL TERPEV( CWT_S, EVECC, GL( 0,LC ), GU( 1,1,LC ), MAZIM, &
                           MXCMU, MXUMU, NN, NSTR, NUMU, WK, YLMC, &
                           YLMU )
!c                                            ** Interpolate source terms
!c                                            ** to user angles

            CALL TERPSO( CWT_S, DELM0, FBEAM, GL( 0,LC ), MAZIM, MXCMU, &
                           PLANK, NUMU, NSTR, OPRIM( LC ), PI, YLM0, &
                           YLMC, YLMU, PSI0, PSI1, XR0( LC ), &
                           XR1( LC ), Z0, Z1, ZJ, ZBEAM( 1,LC ), &
                           Z0U( 1,LC ), Z1U( 1,LC ) )

            END IF

   80    CONTINUE

!c ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


!c                      ** Set coefficient matrix of equations combining
!c                      ** boundary and layer interface conditions

         CALL SETMTX( BDR, CBAND, CMU_S, CWT_S, DELM0, DTAUCP, GC, KK, &
                     LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, &
                     NNLYRI, NN, NSTR, TAUCPR, WK )

!c                      ** Solve for constants of integration in homo-
!c                      ** geneous solution (general boundary conditions)

		if ( .NOT. (PLANK .and. MAZIM .eq. 0)) then 
			ZPLK0 = 0.
			ZPLK1 = 0.
		endif

       CALL SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU_S, CWT_S, EXPBEA, &
                     FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, MI, &
                     MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI, PI, &
                     TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

!c                                  ** Compute upward and downward fluxes


         UUM = 0.0

!c       ! USRANG is always true
!c                                     ** Compute azimuthal intensity
!c                                     ** components at user angles

         CALL USRINT( BPLANK, CMU_S, CWT_S, DELM0, DTAUCP, EMU, EXPBEA, &
                        FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL, &
                        LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT, NLYR, &
                        NN, NSTR, PLANK, NUMU, NTAU, PI, RMU, TAUCPR, &
                        TPLANK, UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U, &
                        ZZ, ZPLK0, ZPLK1, UUM )



         IF( MAZIM.EQ.0 ) THEN
!c                               ** Save azimuthally averaged intensities

!C    we only ever have one angle set
			IU = 1
			J = 1

            DO 130 LU = 1, NTAU

				U0U( IU, LU ) = UUM( IU, LU )
				UU( IU, LU, J ) = UUM( IU, LU )

  130       CONTINUE
!c                              ** Print azimuthally averaged intensities
!c                              ** at user angles



         ELSE
!c                                ** Increment intensity by current
!c                                ** azimuthal component (Fourier
!c                                ** cosine series);  Eq SD(2), STWL(6)
            AZERR  = 0.0

!C        We can get away with this because we ever only run one angle. 
			J=1
			IU = 1

               
			   IF(NSTR <= 16)THEN	
			   	mazim_phirad = MAZIM*PHIRAD( J )
			   	if (mazim_phirad < approx_limit) then 
			   		COSPHI = 1. - mazim_phirad*mazim_phirad/2.0
			   	else
                 	COSPHI = COS( mazim_phirad )
               endif
        		ELSE
               COSPHI = COSPHI_S(MAZIM,J)
            ENDIF               

               DO 160 LU = 1, NTAU

                     AZTERM = UUM( IU, LU )*COSPHI
                     UU( IU, LU, J ) = UU( IU, LU, J ) + AZTERM
                     AZERR  = MAX( AZERR, RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )

  160          CONTINUE

            IF( AZERR.LE.ACCUR ) KCONV  = KCONV + 1

		IF( (KCONV.GE. 2) ) then 
			GO TO 190
		  endif

         END IF

  180 CONTINUE

!c ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


  190 CONTINUE
  

!c                                    ** Apply Nakajima/Tanaka intensity
!c                                    ** corrections


!c                                          ** Print intensities

      END subroutine DISORT_DP_GC






      SUBROUTINE ASYMTX( AA, EVEC, EVAL, M, IA, IEVEC, IER, WKD, AAD, &
                        EVECD, EVALD )

!c    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!c       Solves eigenfunction problem for real asymmetric matrix
!c       for which it is known a priori that the eigenvalues are real.

!c       This is an adaptation of a subroutine EIGRF in the IMSL
!c       library to use real instead of complex arithmetic, accounting
!c       for the known fact that the eigenvalues and eigenvectors in
!c       the discrete ordinate solution are real.  Other changes include
!c       putting all the called subroutines in-line, deleting the
!c       performance index calculation, updating many DO-loops
!c       to Fortran77, and in calculating the machine precision
!c       TOL instead of specifying it in a data statement.

!c       EIGRF is based primarily on EISPACK routines.  The matrix is
!c       first balanced using the Parlett-Reinsch algorithm.  Then
!c       the Martin-Wilkinson algorithm is applied.

!c       There is a statement 'J  = WKD( I )' that converts a double
!c       precision variable to an integer variable, that seems dangerous
!c       to us in principle, but seems to work fine in practice.

!c       References:
!c          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!c             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!c             Sources and Development of Mathematical Software,
!c             Prentice-Hall, Englewood Cliffs, NJ
!c         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!c             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!c         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!c             Clarendon Press, Oxford


!c     .. Scalar Arguments ..

      INTEGER, intent(in) ::   IA, IEVEC, M
	  integer, intent(inout) :: IER
!c     ..
!c     .. Array Arguments ..

      REAL, intent(inout) ::     AA( IA, M ), EVAL( M ), EVEC( IEVEC, M )
      DOUBLE PRECISION, intent(inout) :: AAD( IA, M ), EVALD( M ), EVECD( IA, M ), &
                      WKD( * )
!c     ..
!c     .. Local Scalars ..

      LOGICAL   NOCONV, NOTLAS
      INTEGER   I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
      DOUBLE PRECISION COL, DISCRI, F, G, H, &
                       P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T, &
                       UU, VV, W, X, Y, Z, sqrt_discri
!c     ..
!c     .. External Functions ..

!c     ..
      real*8, parameter :: C1 = 0.4375D0, C2=0.5D0, C3=0.75D0 , &
               C4=0.95D0, C5=16.D0, C6=256.D0
               


      IER  = 0


!c                           ** Handle 1x1 and 2x2 special cases
      IF( M.EQ.1 ) THEN

         EVAL( 1 )   = AA( 1,1 )
         EVEC( 1,1 ) = 1.0
         RETURN

      ELSE IF( M.EQ.2 ) THEN

         DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 + 4.*AA( 1,2 )*AA( 2,1 )


         SGN  = ONE

         IF( AA( 1,1 ) .LT. AA( 2,2 ) ) SGN  = - ONE
			sqrt_discri = SQRT( DISCRI )
         EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*sqrt_discri )
         EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*sqrt_discri )
         EVEC( 1,1 ) = 1.0
         EVEC( 2,2 ) = 1.0

         IF( AA( 1,1 ) .EQ. AA( 2,2 ) .AND. &
            ( AA( 2,1 ).EQ.0.0 .OR. AA( 1,2 ).EQ.0.0 ) ) THEN

            RNORM = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) + ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
            W     = TOL * RNORM
            EVEC( 2,1 ) =   AA( 2,1 ) / W
            EVEC( 1,2 ) = - AA( 1,2 ) / W

         ELSE

            EVEC( 2,1 ) = AA( 2,1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
            EVEC( 1,2 ) = AA( 1,2 ) / ( EVAL( 2 ) - AA( 1,1 ) )

         END IF

         RETURN

      END IF

        AAD(1:M,1:M) = 1.d0 * AA(1:M,1:M)


!c                                ** Initialize output variables
      IER  = 0
      EVALD(1:M) = 0.d0
      EVECD(:,:) = 0.d0

	
      DO 40 I = 1, M

         EVECD( I, I ) = ONE

   40 CONTINUE

!c                  ** Balance the input matrix and reduce its norm by
!c                  ** diagonal similarity transformation stored in WK;
!c                  ** then search for rows isolating an eigenvalue
!c                  ** and push them down
      RNORM  = ZERO
      L  = 1
      K  = M

   50 CONTINUE
      KKK  = K

      DO 90 J = KKK, 1, -1

         ROW  = ZERO

         DO 60 I = 1, K
            IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )
   60    CONTINUE

         IF( ROW.EQ.ZERO ) THEN

            WKD( K ) = J

            IF( J.NE.K ) THEN

               DO 70 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, K )
                  AAD( I, K ) = REPL
   70          CONTINUE

               DO 80 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( K, I )
                  AAD( K, I ) = REPL
   80          CONTINUE

            END IF

            K  = K - 1
            GO TO  50

         END IF

   90 CONTINUE
!c                                ** Search for columns isolating an
!c                                ** eigenvalue and push them left
  100 CONTINUE
      LLL  = L

      DO 140 J = LLL, K

         COL  = ZERO

         DO 110 I = L, K
            IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )
  110    CONTINUE

         IF( COL.EQ.ZERO ) THEN

            WKD( L ) = J

            IF( J.NE.L ) THEN

               DO 120 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, L )
                  AAD( I, L ) = REPL
  120          CONTINUE

               DO 130 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( L, I )
                  AAD( L, I ) = REPL
  130          CONTINUE

            END IF

            L  = L + 1
            GO TO  100

         END IF

  140 CONTINUE

!c                           ** Balance the submatrix in rows L through K
  
        WKD(L:K) = ONE

  160 CONTINUE
      NOCONV = .FALSE.

      DO 220 I = L, K

         COL  = ZERO
         ROW  = ZERO

         DO 170 J = L, K

            IF( J.NE.I ) THEN
               COL  = COL + ABS( AAD( J,I ) )
               ROW  = ROW + ABS( AAD( I,J ) )
            END IF

  170    CONTINUE

         F  = ONE
         G  = ROW / C5
         H  = COL + ROW

  180    CONTINUE
         IF( COL.LT.G ) THEN

            F    = F*C5
            COL  = COL*C6
            GO TO  180

         END IF

         G  = ROW*C5

  190    CONTINUE
         IF( COL.GE.G ) THEN

            F    = F / C5
            COL  = COL / C6
            GO TO  190

         END IF
!c                                                ** Now balance
         IF( ( COL + ROW ) / F.LT.C4*H ) THEN

            WKD( I ) = WKD( I )*F
            NOCONV = .TRUE.

            DO 200 J = L, M
               AAD( I, J ) = AAD( I, J ) / F
  200       CONTINUE

  
  			   AAD(1:K,I) = AAD(1:K,I) * F

         END IF

  220 CONTINUE


      IF( NOCONV ) GO TO  160
!c                                   ** Is A already in Hessenberg form?
      IF( K-1 .LT. L+1 ) GO TO  370

!c                                   ** Transfer A to a Hessenberg form
      DO 310 N = L + 1, K - 1

         H  = ZERO
         WKD( N + M ) = ZERO
         SCALE  = ZERO
!c                                                 ** Scale column
         DO 230 I = N, K
            SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
  230    CONTINUE

         IF( SCALE.NE.ZERO ) THEN

            DO 240 I = K, N, -1
               WKD( I + M ) = AAD( I, N - 1 ) / SCALE
               H  = H + WKD( I + M )**2
  240       CONTINUE

            G    = - SIGN( SQRT( H ), WKD( N + M ) )
            H    = H - WKD( N + M )*G
            WKD( N + M ) = WKD( N + M ) - G
!c                                            ** Form (I-(U*UT)/H)*A
            DO 270 J = N, M

               F  = ZERO

               DO 250 I = K, N, -1
                  F  = F + WKD( I + M )*AAD( I, J )
  250          CONTINUE

               DO 260 I = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
  260          CONTINUE

  270       CONTINUE
!c                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 300 I = 1, K

               F  = ZERO

               DO 280 J = K, N, -1
                  F  = F + WKD( J + M )*AAD( I, J )
  280          CONTINUE

               DO 290 J = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
  290          CONTINUE

  300       CONTINUE

            WKD( N + M ) = SCALE*WKD( N + M )
            AAD( N, N - 1 ) = SCALE*G

         END IF

  310 CONTINUE


      DO 360 N = K - 2, L, -1

         N1   = N + 1
         N2   = N + 2
         F  = AAD( N + 1, N )

         IF( F.NE.ZERO ) THEN

            F  = F*WKD( N + 1 + M )

            DO 320 I = N + 2, K
               WKD( I + M ) = AAD( I, N )
  320       CONTINUE

            IF( N + 1.LE.K ) THEN

               DO 350 J = 1, M

                  G  = ZERO

                  DO 330 I = N + 1, K
                     G  = G + WKD( I + M )*EVECD( I, J )
  330             CONTINUE

                  G  = G / F

                  DO 340 I = N + 1, K
                     EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
  340             CONTINUE

  350          CONTINUE

            END IF

         END IF

  360 CONTINUE


  370 CONTINUE

      N  = 1

      DO 390 I = 1, M

         DO 380 J = N, M
            RNORM  = RNORM + ABS( AAD( I,J ) )
  380    CONTINUE

         N  = I

         IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )

  390 CONTINUE

      N  = K
      T  = ZERO

!c                                      ** Search for next eigenvalues
  400 CONTINUE
      IF( N.LT.L ) GO TO  550

      IN  = 0
      N1  = N - 1
      N2  = N - 2
!c                          ** Look for single small sub-diagonal element
  410 CONTINUE

      DO 420 I = L, N

         LB  = N + L - I

         IF( LB.EQ.L ) GO TO  430

         S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )

         IF( S.EQ.ZERO ) S  = RNORM

         IF( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430

  420 CONTINUE


  430 CONTINUE
      X  = AAD( N, N )

      IF( LB.EQ.N ) THEN
!c                                        ** One eigenvalue found
         AAD( N, N ) = X + T
         EVALD( N ) = AAD( N, N )
         N  = N1
         GO TO  400

      END IF

      Y  = AAD( N1, N1 )
      W  = AAD( N, N1 )*AAD( N1, N )

      IF( LB.EQ.N1 ) THEN
!c                                        ** Two eigenvalues found
         P  = ( Y - X )*C2
         Q  = P**2 + W
         Z  = SQRT( ABS( Q ) )
         AAD( N, N ) = X + T
         X  = AAD( N, N )
         AAD( N1, N1 ) = Y + T
!c                                        ** Real pair
         Z  = P + SIGN( Z, P )
         EVALD( N1 ) = X + Z
         EVALD( N ) = EVALD( N1 )

         IF( Z.NE.ZERO ) EVALD( N ) = X - W / Z

         X  = AAD( N, N1 )
!c                                  ** Employ scale factor in case
!c                                  ** X and Z are very small
         R  = SQRT( X*X + Z*Z )
         P  = X / R
         Q  = Z / R
!c                                             ** Row modification
         DO 440 J = N1, M
            Z  = AAD( N1, J )
            AAD( N1, J ) = Q*Z + P*AAD( N, J )
            AAD( N, J ) = Q*AAD( N, J ) - P*Z
  440    CONTINUE
!c                                             ** Column modification
         DO 450 I = 1, N
            Z  = AAD( I, N1 )
            AAD( I, N1 ) = Q*Z + P*AAD( I, N )
            AAD( I, N ) = Q*AAD( I, N ) - P*Z
  450    CONTINUE
!c                                          ** Accumulate transformations
         DO 460 I = L, K
            Z  = EVECD( I, N1 )
            EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
            EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
  460    CONTINUE

         N  = N2
         GO TO  400

      END IF


      IF( IN.EQ.30 ) THEN

!c                    ** No convergence after 30 iterations; set error
!c                    ** indicator to the index of the current eigenvalue
         IER  = N
         GO TO  700

      END IF
!c                                                  ** Form shift
      IF( IN.EQ.10 .OR. IN.EQ.20 ) THEN

         T  = T + X

         DO 470 I = L, N
            AAD( I, I ) = AAD( I, I ) - X
  470    CONTINUE

         S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
         X  = C3*S
         Y  = X
         W  = -C1*S**2

      END IF


      IN  = IN + 1

!c                ** Look for two consecutive small sub-diagonal elements

      DO 480 J = LB, N2
         I  = N2 + LB - J
         Z  = AAD( I, I )
         R  = X - Z
         S  = Y - Z
         P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
         Q  = AAD( I + 1, I + 1 ) - Z - R - S
         R  = AAD( I + 2, I + 1 )
         S  = ABS( P ) + ABS( Q ) + ABS( R )
         P  = P / S
         Q  = Q / S
         R  = R / S

         IF( I.EQ.LB ) GO TO  490

         UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
         VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) + ABS( AAD( I+1, I+1 ) ) )

         IF( UU .LE. TOL*VV ) GO TO  490

  480 CONTINUE

  490 CONTINUE
      AAD( I+2, I ) = ZERO

      DO 500 J = I + 3, N
         AAD( J, J - 2 ) = ZERO
         AAD( J, J - 3 ) = ZERO
  500 CONTINUE

!c             ** Double QR step involving rows K to N and columns M to N

      DO 540 KA = I, N1

         NOTLAS = KA.NE.N1

         IF( KA.EQ.I ) THEN

            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )

            IF( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )

         ELSE

            P  = AAD( KA, KA - 1 )
            Q  = AAD( KA + 1, KA - 1 )
            R  = ZERO

            IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )

            X  = ABS( P ) + ABS( Q ) + ABS( R )

            IF( X.EQ.ZERO ) GO TO  540

            P  = P / X
            Q  = Q / X
            R  = R / X
            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            AAD( KA, KA - 1 ) = -S*X

         END IF

         P  = P + S
         X  = P / S
         Y  = Q / S
         Z  = R / S
         Q  = Q / P
         R  = R / P
!c                                              ** Row modification
         DO 510 J = KA, M

            P  = AAD( KA, J ) + Q*AAD( KA + 1, J )

            IF( NOTLAS ) THEN

               P  = P + R*AAD( KA + 2, J )
               AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z

            END IF

            AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
            AAD( KA, J ) = AAD( KA, J ) - P*X

  510    CONTINUE
!c                                                 ** Column modification
         DO 520 II = 1, MIN( N, KA + 3 )

            P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*AAD( II, KA + 2 )
               AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R

            END IF

            AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
            AAD( II, KA ) = AAD( II, KA ) - P

  520    CONTINUE
!c                                          ** Accumulate transformations
         DO 530 II = L, K

            P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*EVECD( II, KA + 2 )
               EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R

            END IF

            EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
            EVECD( II, KA ) = EVECD( II, KA ) - P

  530    CONTINUE

  540 CONTINUE

      GO TO  410
!c                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

      IF( RNORM.NE.ZERO ) THEN

         DO 580 N = M, 1, -1

            N2   = N
            AAD( N, N ) = ONE

            DO 570 I = N - 1, 1, -1

               W  = AAD( I, I ) - EVALD( N )

               IF( W.EQ.ZERO ) W  = TOL*RNORM

               R  = AAD( I, N )

               DO 560 J = N2, N - 1
                  R  = R + AAD( I, J )*AAD( J, N )
  560          CONTINUE

               AAD( I, N ) = -R / W
               N2   = I

  570       CONTINUE

  580    CONTINUE
!c                      ** End backsubstitution vectors of isolated evals
         DO 600 I = 1, M

            IF( I.LT.L .OR. I.GT.K ) THEN

               DO 590 J = I, M
                  EVECD( I, J ) = AAD( I, J )
  590          CONTINUE

            END IF

  600    CONTINUE
!c                                   ** Multiply by transformation matrix
         IF( K.NE.0 ) THEN

            DO 630 J = M, L, -1

               DO 620 I = L, K

                  Z  = ZERO

                  DO 610 N = L, MIN( J, K )
                     Z  = Z + EVECD( I, N )*AAD( N, J )
  610             CONTINUE

                  EVECD( I, J ) = Z

  620          CONTINUE

  630       CONTINUE

         END IF

      END IF


      DO 650 I = L, K

         DO 640 J = 1, M
            EVECD( I, J ) = EVECD( I, J ) * WKD( I )
  640    CONTINUE

  650 CONTINUE

!c                           ** Interchange rows if permutations occurred
      DO 670 I = L-1, 1, -1

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 660 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  660       CONTINUE

         END IF

  670 CONTINUE


      DO 690 I = K + 1, M

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 680 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  680       CONTINUE

         END IF

  690 CONTINUE

!c                         ** Put results into output arrays
  700 CONTINUE

      DO 720 J = 1, M

         EVAL( J ) = EVALD( J )

         DO 710 K = 1, M
            EVEC( J, K ) = EVECD( J, K )
  710    CONTINUE

  720 CONTINUE

      END SUBROUTINE ASYMTX



      SUBROUTINE SETDIS( CMU, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, &
                        FLYR, GL, IBCND, LAYRU, LYRCUT, MAXMOM, MAXUMU, &
                        MXCMU, NCUT, NLYR, NTAU, NN, NSTR, PLANK, NUMU, &
                        ONLYFL, CORINT, OPRIM, PMOM, SSALB, TAUC, &
                        TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU, &
                        USRANG )

!c          Perform miscellaneous setting-up operations
!c
!c    INPUT :  all are DISORT input variables (see DOC file)
!c
!c
!c    O U T P U T     V A R I A B L E S:
!c
!c       NTAU,UTAU   if USRTAU = FALSE (defined in DISORT.doc)
!c       NUMU,UMU    if USRANG = FALSE (defined in DISORT.doc)
!c
!c       EXPBEA      transmission of direct beam
!c
!c       FLYR        separated fraction in delta-M method
!c
!c       GL          phase function Legendre coefficients multiplied
!c                   by (2L+1) and single-scatter albedo
!c
!c       LAYRU       Computational layer in which UTAU falls
!c
!c       LYRCUT      flag as to whether radiation will be zeroed
!c                   below layer NCUT
!c
!c       NCUT        computational layer where absorption
!c                   optical depth first exceeds  ABSCUT
!c
!c       NN          NSTR / 2
!c
!c       OPRIM       delta-M-scaled single-scatter albedo
!c
!c       TAUCPR      delta-M-scaled optical depth
!c
!c       UTAUPR      delta-M-scaled version of  UTAU
!c
!c   Called by- DISORT
!c   Calls- QGAUSN, ERRMSG
!c ---------------------------------------------------------------------

!c     .. Scalar Arguments ..

      LOGICAL, intent(in) ::   CORINT, DELTAM, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER, intent(in) ::   IBCND, MAXMOM, MAXUMU, MXCMU, NLYR, NN, NSTR, &
               NTAU, NUMU
      REAL, intent(in) ::      FBEAM, UMU0
	  integer, intent(inout) :: NCUT
	  logical, intent(inout) :: LYRCUT
!c     ..
!c     .. Array Arguments ..

      INTEGER, intent(inout)::   LAYRU( * )
      REAL, intent(inout)::      CMU( : ), DTAUC( * ), DTAUCP( * ), &
               EXPBEA( 0:* ), FLYR( * ), GL( 0:MXCMU, * ), OPRIM( * ), &
               PMOM( 0:MAXMOM, * ), SSALB( * ), TAUC( 0:* ), &
               TAUCPR( 0:* ), UMU( MAXUMU ), UTAU( * ), UTAUPR( * )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, IU, K, LC, LU
      REAL      ABSTAU, F, YESSCT

!c     ..
      real, parameter :: ABSCUT =  10. 


!c                        ** Apply delta-M scaling and move description
!c                        ** of computational layers to local variables
      EXPBEA( 0 ) = 1.0
      TAUCPR( 0 ) = 0.0
      ABSTAU      = 0.0
      YESSCT      = 0.0

      DO 40 LC = 1, NLYR

         YESSCT = YESSCT + SSALB( LC )

!c         PMOM( 0, LC ) = 1.0

         IF( ABSTAU.LT.ABSCUT ) NCUT  = LC

         ABSTAU = ABSTAU + ( 1. - SSALB( LC ) )*DTAUC( LC )


            F  = PMOM( NSTR, LC )
		
			if (F == 0.) then 
				OPRIM( LC )  = SSALB( LC )
				DTAUCP( LC ) = DTAUC( LC )
				TAUCPR( LC ) = TAUCPR( LC - 1 ) + DTAUC( LC )

                DO K = 0, NSTR - 1
					GL(K, LC) = 0.
					if (PMOM(K, LC) .ne. 0.) GL( K, LC ) = ( 2*K + 1 )*SSALB( LC )* PMOM( K,LC )  
				end do

			else 
				OPRIM( LC )  = SSALB( LC )*( 1. - F ) / ( 1. - F*SSALB(LC) )
				DTAUCP( LC ) = ( 1. - F*SSALB( LC ) )*DTAUC( LC )
				TAUCPR( LC ) = TAUCPR( LC - 1 ) + DTAUCP( LC )

				DO 30 K = 0, NSTR - 1
					GL(K, LC) = 0.
					if (PMOM(K, LC) .ne. 0.) GL( K, LC ) = ( 2*K + 1 )*OPRIM( LC )* &
                            ( PMOM( K,LC ) - F ) / ( 1. - F )
   30			CONTINUE
			endif


         FLYR( LC ) = F
         EXPBEA( LC ) = 0.0

         IF( FBEAM.GT.0.0 ) EXPBEA( LC ) = EXP( -TAUCPR( LC )/UMU0 )

   40 CONTINUE
!c                      ** If no thermal emission, cut off medium below
!c                      ** absorption optical depth = ABSCUT ( note that
!c                      ** delta-M transformation leaves absorption
!c                      ** optical depth invariant ).  Not worth the
!c                      ** trouble for one-layer problems, though.
      LYRCUT = .FALSE.
      IF( ABSTAU.GE.ABSCUT .AND. .NOT.PLANK ) LYRCUT = .TRUE.
      IF( .NOT.LYRCUT ) NCUT = NLYR

!c                             ** Set arrays defining location of user
!c                             ** output levels within delta-M-scaled
!c                             ** computational mesh
      DO 70 LU = 1, NTAU

         DO 50 LC = 1, NLYR

            IF( UTAU( LU ).GE.TAUC( LC-1 ) .AND. &
               UTAU( LU ).LE.TAUC( LC ) ) GO TO  60

   50    CONTINUE
         LC   = NLYR

   60    CONTINUE
         UTAUPR( LU ) = TAUCPR( LC - 1 ) + &
                                    ( 1. - SSALB( LC )*FLYR( LC ) )* &
                                    ( UTAU( LU ) - TAUC( LC-1 ) )
         LAYRU( LU ) = LC

   70 CONTINUE

      END SUBROUTINE SETDIS



      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK, &
                        LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, &
                        NNLYRI, NN, NSTR, TAUCPR, WK )

!c        Calculate coefficient matrix for the set of equations
!c        obtained from the boundary conditions and the continuity-
!c        of-intensity-at-layer-interface equations;  store in the
!c        special banded-matrix format required by LINPACK routines
!c
!c
!c
!c   BAND STORAGE
!c
!c      LINPACK requires band matrices to be input in a special
!c      form where the elements of each diagonal are moved up or
!c      down (in their column) so that each diagonal becomes a row.
!c      (The column locations of diagonal elements are unchanged.)
!c
!c      Example:  if the original matrix is
!c
!c          11 12 13  0  0  0
!c          21 22 23 24  0  0
!c           0 32 33 34 35  0
!c           0  0 43 44 45 46
!c           0  0  0 54 55 56
!c           0  0  0  0 65 66
!c
!c      then its LINPACK input form would be:
!c
!c           *  *  *  +  +  +  , * = not used
!c           *  * 13 24 35 46  , + = used for pivoting
!c           * 12 23 34 45 56
!c          11 22 33 44 55 66
!c          21 32 43 54 65  *
!c
!c      If A is a band matrix, the following program segment
!c      will convert it to the form (ABD) required by LINPACK
!c      band-matrix routines:
!c
!c               N  = (column dimension of A, ABD)
!c               ML = (band width below the diagonal)
!c               MU = (band width above the diagonal)
!c               M = ML + MU + 1
!c               DO J = 1, N
!c                  I1 = MAX(1, J-MU)
!c                  I2 = MIN(N, J+ML)
!c                  DO I = I1, I2
!c                     K = I - J + M
!c                     ABD(K,J) = A(I,J)
!c                  END DO
!c               END DO
!c
!c      This uses rows  ML+1  through  2*ML+MU+1  of ABD.
!c      The total number of rows needed in ABD is  2*ML+MU+1 .
!c      In the example above, N = 6, ML = 1, MU = 2, and the
!c      row dimension of ABD must be >= 5.
!c

!c     .. Scalar Arguments ..

      LOGICAL, intent(in) ::   LAMBER, LYRCUT
      INTEGER, intent(in) ::   MI, MI9M2, MXCMU, NCUT, NN, NNLYRI, NSTR
	  integer, intent(inout) :: NCOL
      REAL, intent(in) ::     DELM0
!c     ..
!c     .. Array Arguments ..

      REAL, intent(inout) ::     BDR( MI, 0:MI ), CBAND( MI9M2, NNLYRI ), CMU( MXCMU ), &
               CWT( MXCMU ), DTAUCP( * ), GC( MXCMU, MXCMU, * ), &
               KK( MXCMU, * ), TAUCPR( 0:* ), WK( MXCMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      REAL      EXPA, SUM


       CBAND = 0.0

      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
	  
!c                         ** Use continuity conditions of Eq. STWJ(17)
!c                         ** to form coefficient matrix in STWJ(20);
!c                         ** employ scaling transformation STWJ(22)
      DO 60 LC = 1, NCUT

         DO 10 IQ = 1, NN
            WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
   10    CONTINUE

         JCOL  = 0

         DO 30 IQ = 1, NN

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 20 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
               IROW  = IROW + 1
   20       CONTINUE

            JCOL  = JCOL + 1

   30    CONTINUE


         DO 50 IQ = NN + 1, NSTR

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL


            DO 40 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )* &
                                               WK( NSTR + 1 - IQ )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
               IROW  = IROW + 1
   40       CONTINUE

            JCOL  = JCOL + 1

   50    CONTINUE

   60 CONTINUE
!c                  ** Use top boundary condition of STWJ(20a) for
!c                  ** first layer
      JCOL  = 0

      DO 80 IQ = 1, NN

         EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
         IROW  = NSHIFT - JCOL + NN

         DO 70 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
            IROW  = IROW + 1
   70    CONTINUE

         JCOL  = JCOL + 1

   80 CONTINUE


      DO 100 IQ = NN + 1, NSTR

         IROW  = NSHIFT - JCOL + NN

         DO 90 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
            IROW  = IROW + 1
   90    CONTINUE

         JCOL  = JCOL + 1

  100 CONTINUE
!c                           ** Use bottom boundary condition of
!c                           ** STWJ(20c) for last layer

      NNCOL = NCOL - NSTR
      JCOL  = 0

      DO 130 IQ = 1, NN

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR

         DO 120 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

!c                          ** No azimuthal-dependent intensity if Lam-
!c                          ** bert surface; no intensity component if
!c                          ** truncated bottom layer

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )

            ELSE

               SUM  = 0.0

               DO 110 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )* &
                          GC( NN + 1 - K, IQ, NCUT )
  110          CONTINUE

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) - &
                                     ( 1.+ DELM0 )*SUM
            END IF

            IROW  = IROW + 1

  120    CONTINUE

         JCOL  = JCOL + 1

  130 CONTINUE


      DO 160 IQ = NN + 1, NSTR

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         EXPA   = WK( NSTR + 1 - IQ )

         DO 150 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA

            ELSE

               SUM  = 0.0

               DO 140 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )* &
                              GC( NN + 1 - K, IQ, NCUT )
  140          CONTINUE

               CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) - &
                                     ( 1.+ DELM0 )*SUM )*EXPA
            END IF

            IROW  = IROW + 1

  150    CONTINUE

         JCOL  = JCOL + 1

  160 CONTINUE

      END SUBROUTINE SETMTX

	

      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM, &
                        MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC, &
                        AAD, EVECCD, EVALD, WKD )

!c         Solves eigenvalue/vector problem necessary to construct
!c         homogeneous part of discrete ordinate solution; STWJ(8b),
!c         STWL(23f)
!c         ** NOTE ** Eigenvalue problem is degenerate when single
!c                    scattering albedo = 1;  present way of doing it
!c                    seems numerically more stable than alternative
!c                    methods that we tried
!c

!c     .. Scalar Arguments ..

      INTEGER, intent(in)::   MAZIM, MI, MXCMU, NN, NSTR
!c     ..
!c     .. Array Arguments ..

      REAL, intent(inout) ::      AMB( MI, MI ), APB( MI, MI ), ARRAY( MI, * ), &
               CC( MXCMU, MXCMU ), CMU( MXCMU ), CWT( MXCMU ), &
               EVAL( MI ), EVECC( MXCMU, MXCMU ), GC( MXCMU, MXCMU ), &
               GL( 0:MXCMU ), KK( MXCMU ), YLMC( 0:MXCMU, MXCMU )
      DOUBLE PRECISION, intent(inout) :: AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ), &
                      WKD( MXCMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IER, IQ, JQ, KQ, L
      REAL      ALPHA, BETA, GPMIGM, GPPLGM, SUM

!c                             ** Calculate quantities in Eqs. SS(5-6),
!c                             ** STWL(8b,15,23f)
      DO 40 IQ = 1, NN

         DO 20 JQ = 1, NSTR

            SUM  = 0.0
            DO 10 L = MAZIM, NSTR - 1
               SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
   10       CONTINUE

            CC( IQ, JQ ) = 0.5*SUM*CWT( JQ )

   20    CONTINUE

         DO 30 JQ = 1, NN
!c                             ** Fill remainder of array using symmetry
!c                             ** relations  C(-mui,muj) = C(mui,-muj)
!c                             ** and        C(-mui,-muj) = C(mui,muj)

            CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
            CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )

!c                                       ** Get factors of coeff. matrix
!c                                       ** of reduced eigenvalue problem

            ALPHA  = CC( IQ, JQ ) / CMU( IQ )
            BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
            AMB( IQ, JQ ) = ALPHA - BETA
            APB( IQ, JQ ) = ALPHA + BETA

   30    CONTINUE

         AMB( IQ, IQ ) = AMB( IQ, IQ ) - 1.0 / CMU( IQ )
         APB( IQ, IQ ) = APB( IQ, IQ ) - 1.0 / CMU( IQ )

   40 CONTINUE
!c                      ** Finish calculation of coefficient matrix of
!c                      ** reduced eigenvalue problem:  get matrix
!c                      ** product (alfa+beta)*(alfa-beta); SS(12),
!c                      ** STWL(23f)
      DO 70 IQ = 1, NN

         DO 60 JQ = 1, NN

            SUM  = 0.
            DO 50 KQ = 1, NN
               SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
   50       CONTINUE

            ARRAY( IQ, JQ ) = SUM

   60    CONTINUE

   70 CONTINUE
!c                      ** Find (real) eigenvalues and eigenvectors

      CALL ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WKD, AAD, &
                  EVECCD, EVALD )


      DO 80 IQ = 1, NN
         EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
         KK( IQ + NN ) = EVAL( IQ )
!c                                      ** Add negative eigenvalue
         KK( NN + 1 - IQ ) = -EVAL( IQ )
   80 CONTINUE

!c                          ** Find eigenvectors (G+) + (G-) from SS(10)
!c                          ** and store temporarily in APB array
      DO 110 JQ = 1, NN

         DO 100 IQ = 1, NN

            SUM  = 0.
            DO 90 KQ = 1, NN
               SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
   90       CONTINUE

            APB( IQ, JQ ) = SUM / EVAL( JQ )

  100    CONTINUE

  110 CONTINUE


      DO 130 JQ = 1, NN

         DO 120 IQ = 1, NN

            GPPLGM = APB( IQ, JQ )
            GPMIGM = EVECC( IQ, JQ )
!c                                ** Recover eigenvectors G+,G- from
!c                                ** their sum and difference; stack them
!c                                ** to get eigenvectors of full system
!c                                ** SS(7) (JQ = eigenvector number)

            EVECC( IQ,      JQ ) = 0.5*( GPPLGM + GPMIGM )
            EVECC( IQ + NN, JQ ) = 0.5*( GPPLGM - GPMIGM )

!c                                ** Eigenvectors corresponding to
!c                                ** negative eigenvalues (corresp. to
!c                                ** reversing sign of 'k' in SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )

  120    CONTINUE

  130 CONTINUE

      END SUBROUTINE SOLEIG

	
      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA, &
                        FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, &
                        MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI, &
                        PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )


!c     .. Scalar Arguments ..

      LOGICAL, intent(in) ::   LAMBER, LYRCUT
      INTEGER, intent(in) ::   MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL, intent(in) ::      BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0
!c     ..
!c     .. Array Arguments ..

      INTEGER, intent(inout) ::   IPVT( * )
      REAL, intent(inout) ::      B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ), &
               CBAND( MI9M2, NNLYRI ), CMU( MXCMU ), CWT( MXCMU ), &
               EXPBEA( 0:* ), LL( MXCMU, * ), TAUCPR( 0:* ), &
               Z( NNLYRI ), ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ), &
               ZZ( MXCMU, * )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IPNT, IQ, IT, JQ, LC, NCD
      REAL      RCOND, SUM
!c     ..
!c     .. External Subroutines ..

      EXTERNAL  SGBCO, SGBSL

		B(1:NNLYRI) = 0.0
      
!c                              ** Construct B,  STWJ(20a,c) for
!c                              ** parallel beam + bottom reflection +
!c                              ** thermal emission at top and/or bottom

      IF( MAZIM.GT.0 .AND. FBEAM.GT.0.0 ) THEN

!c                                         ** Azimuth-dependent case
!c                                         ** (never called if FBEAM = 0)
         IF( LYRCUT .OR. LAMBER ) THEN

!c               ** No azimuthal-dependent intensity for Lambert surface;
!c               ** no intensity component for truncated bottom layer

            DO 10 IQ = 1, NN
!c                                                  ** Top boundary
               B( IQ ) = -ZZ( NN + 1 - IQ, 1 )
!c                                                  ** Bottom boundary

               B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )

   10       CONTINUE


         ELSE

            DO 30 IQ = 1, NN

               B( IQ ) = -ZZ( NN + 1 - IQ, 1 )

               SUM  = 0.
               DO 20 JQ = 1, NN
                  SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )* &
                              ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
   20          CONTINUE

               B( NCOL - NN + IQ ) = SUM
               IF( FBEAM.GT.0.0 ) B( NCOL - NN + IQ ) = SUM + &
                  ( BDR( IQ,0 )*UMU0*FBEAM/PI - ZZ( IQ+NN,NCUT ) )* &
                  EXPBEA( NCUT )

   30       CONTINUE

         END IF
!c                             ** Continuity condition for layer
!c                             ** interfaces of Eq. STWJ(20b)
         IT  = NN

         DO 50 LC = 1, NCUT - 1

            DO 40 IQ = 1, NSTR
               IT  = IT + 1
               B( IT ) = ( ZZ( IQ,LC+1 ) - ZZ( IQ,LC ) )*EXPBEA( LC )
   40       CONTINUE

   50    CONTINUE


      ELSE
!c                                   ** Azimuth-independent case

         IF( FBEAM.EQ.0.0 ) THEN

            DO 60 IQ = 1, NN
!c                                      ** Top boundary

               B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK

   60       CONTINUE


            IF( LYRCUT ) THEN
!c                               ** No intensity component for truncated
!c                               ** bottom layer
               DO 70 IQ = 1, NN
!c                                      ** Bottom boundary

                  B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) - &
                                         ZPLK1( IQ + NN, NCUT ) * &
                                         TAUCPR( NCUT )
   70          CONTINUE


            ELSE

               DO 90 IQ = 1, NN

                  SUM  = 0.
                  DO 80 JQ = 1, NN
                     SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )* &
                             ( ZPLK0( NN+1-JQ, NCUT ) + &
                               ZPLK1( NN+1-JQ, NCUT ) *TAUCPR( NCUT ) )
   80             CONTINUE

                  B( NCOL - NN + IQ ) = 2.*SUM + BEM( IQ )*BPLANK - &
                                       ZPLK0( IQ + NN, NCUT ) - &
                                       ZPLK1( IQ + NN, NCUT ) * &
                                       TAUCPR( NCUT )
   90          CONTINUE

            END IF ! if LYRCUT
!c                             ** Continuity condition for layer
!c                             ** interfaces, STWJ(20b)
            IT  = NN
            DO 110 LC = 1, NCUT - 1

               DO 100 IQ = 1, NSTR
                  IT  = IT + 1
                  B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) + &
                           ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )* &
                           TAUCPR( LC )
  100          CONTINUE

  110       CONTINUE


         ELSE ! if FBEAM is not 0.0


            DO 120 IQ = 1, NN
               B( IQ ) = -ZZ( NN + 1 - IQ, 1 ) - &
                        ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
  120       CONTINUE

            IF( LYRCUT ) THEN

               DO 130 IQ = 1, NN
                  B( NCOL-NN+IQ ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT) &
                                   - ZPLK0(IQ+NN, NCUT) &
                                   - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  130          CONTINUE


            ELSE

               DO 150 IQ = 1, NN

                  SUM  = 0.
                  DO 140 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ) &
                               * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT) &
                                 + ZPLK0(NN+1-JQ, NCUT) &
                                 + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
  140             CONTINUE

                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI &
                                     - ZZ(IQ+NN, NCUT) ) * EXPBEA(NCUT) &
                                 + BEM(IQ) * BPLANK &
								 - ZPLK0(IQ+NN, NCUT) &
                                 - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  150          CONTINUE

            END IF ! if LYRCUT


            IT  = NN

            DO 170 LC = 1, NCUT - 1

               DO 160 IQ = 1, NSTR

                  IT  = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC) &
                         + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) + &
                         ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)
  160          CONTINUE

  170       CONTINUE

         END IF ! if FBEAM check
         
      END IF
!c                     ** Find L-U (lower/upper triangular) decomposition
!c                     ** of band matrix CBAND and test if it is nearly
!c                     ** singular (note: CBAND is destroyed)
!c                     ** (CBAND is in LINPACK packed format)
      RCOND  = 0.0
      NCD    = 3*NN - 1

      CALL SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )


!c                   ** Solve linear system with coeff matrix CBAND
!c                   ** and R.H. side(s) B after CBAND has been L-U
!c                   ** decomposed.  Solution is returned in B.

      CALL SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )

!c                   ** Zero CBAND (it may contain 'foreign'
!c                   ** elements upon returning from LINPACK);
!c                   ** necessary to prevent errors

      CBAND = 0.0

      DO 190 LC = 1, NCUT

         IPNT  = LC*NSTR - NN

         DO 180 IQ = 1, NN
            LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
            LL( IQ + NN,     LC ) = B( IQ + IPNT )
  180    CONTINUE

  190 CONTINUE

      END SUBROUTINE SOLVE0
	

      SUBROUTINE SURFAC( ALBEDO, LAMBER, MI, MAZIM, &
                        MXUMU, NN, BDR, EMU, BEM, RMU )

!c       Computes user's surface bidirectional properties, STWL(41)

!c     .. Parameters ..
!c     ..
!c     .. Scalar Arguments ..

      LOGICAL, intent(in) ::    LAMBER
      INTEGER, intent(in) ::    MAZIM, MI, MXUMU, NN
      REAL, intent(in) ::       ALBEDO
!c     ..
!c     .. Array Arguments ..
      REAL, intent(inout) ::       BDR( MI, 0:MI ), BEM( MI ),  EMU( MXUMU ), &
              RMU( MXUMU, 0:MI )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IU
      
!c We're only ever doing lambertian calculations here, so this subroutine has been 
!c summarily gutted. 

      IF( LAMBER .AND. MAZIM.EQ.0 ) THEN   !outer IF
      
      	  BEM(1:NN) 	 = 1.0 - ALBEDO
      	  BDR(1:NN,0:NN) = ALBEDO

		  IU = 1
               
		  RMU(IU,0:NN) = ALBEDO
		  EMU( IU ) = 1.0 - ALBEDO
				


	  ELSE
	  
			BDR = 0.0
			BEM=0.0

			EMU = 0.0
			RMU = 0.0


      END IF
      

      END SUBROUTINE SURFAC 


      SUBROUTINE TERPEV( CWT, EVECC, GL, GU, MAZIM, MXCMU, MXUMU, NN, &
                        NSTR, NUMU, WK, YLMC, YLMU )

!c         Interpolate eigenvectors to user angles; Eq SD(8)

!c   Called by- DISORT, ALBTRN
!c --------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      INTEGER   MAZIM, MXCMU, MXUMU, NN, NSTR, NUMU
!c     ..
!c     .. Array Arguments ..

      REAL      CWT( MXCMU ), EVECC( MXCMU, MXCMU ), GL( 0:MXCMU ), &
               GU( MXUMU, MXCMU ), WK( MXCMU ), YLMC( 0:MXCMU, MXCMU ), &
               YLMU( 0:MXCMU, MXUMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, IU, JQ, L
      REAL      SUM
!c     ..


      DO 50 IQ = 1, NSTR

         DO 20 L = MAZIM, NSTR - 1
!c                                   ** Inner sum in SD(8) times all
!c                                   ** factors in outer sum but PLM(mu)
            SUM  = 0.0
            DO 10 JQ = 1, NSTR
               SUM  = SUM + CWT( JQ )*YLMC( L, JQ )*EVECC( JQ, IQ )
   10       CONTINUE

            WK( L + 1 ) = 0.5*GL( L )*SUM

   20    CONTINUE
!c                                    ** Finish outer sum in SD(8)
!c                                    ** and store eigenvectors
			IU = 1
!c         DO 40 IU = 1, NUMU
			

            SUM  = 0.
            DO 30 L = MAZIM, NSTR - 1
               SUM  = SUM + WK( L + 1 )*YLMU( L, IU )
   30       CONTINUE

            IF( IQ.LE.NN ) GU( IU, IQ + NN )       = SUM
            IF( IQ.GT.NN ) GU( IU, NSTR + 1 - IQ ) = SUM

!c   40    CONTINUE

   50 CONTINUE

      END SUBROUTINE TERPEV


      SUBROUTINE TERPSO( CWT, DELM0, FBEAM, GL, MAZIM, MXCMU, PLANK, &
                        NUMU, NSTR, OPRIM, PI, YLM0, YLMC, YLMU, PSI0, &
                        PSI1, XR0, XR1, Z0, Z1, ZJ, ZBEAM, Z0U, Z1U )

!c         Interpolates source functions to user angles, Eq. STWL(30)
!c
!c
!c    I N P U T      V A R I A B L E S:
!c
!c       CWT    :  Weights for Gauss quadrature over angle cosine
!c
!c       DELM0  :  Kronecker delta, delta-sub-m0
!c
!c       GL     :  Delta-M scaled Legendre coefficients of phase function
!c                 (including factors 2L+1 and single-scatter albedo)
!c
!c       MAZIM  :  Order of azimuthal component
!c
!c       OPRIM  :  Single scattering albedo
!c
!c       XR0    :  Expansion of thermal source function, Eq. STWL(24d)
!c
!c       XR1    :  Expansion of thermal source function Eq. STWL(24d)
!c
!c       YLM0   :  Normalized associated Legendre polynomial
!c                 at the beam angle
!c
!c       YLMC   :  Normalized associated Legendre polynomial
!c                 at the quadrature angles
!c
!c       YLMU   :  Normalized associated Legendre polynomial
!c                 at the user angles
!c
!c       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16), STWL(26a)
!c
!c       Z1     :  Solution vectors Z-sub-one  of Eq. SS(16), STWL(26b)
!c
!c       ZJ     :  Solution vector Z-sub-zero after solving Eq. SS(19),
!c                 STWL(24b)
!c
!c       (remainder are DISORT input variables)
!c
!c
!c    O U T P U T     V A R I A B L E S:
!c
!c       ZBEAM  :  Incident-beam source function at user angles
!c
!c       Z0U,Z1U:  Components of a linear-in-optical-depth-dependent
!c                 source (approximating the Planck emission source)
!c
!c
!c   I N T E R N A L    V A R I A B L E S:
!c
!c       PSI0  :  Sum just after square bracket in  Eq. SD(9)
!c       PSI1  :  Sum in Eq. STWL(31d)
!c
!c   Called by- DISORT
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      LOGICAL, intent(in)::   PLANK
      INTEGER, intent(in)::   MAZIM, MXCMU, NSTR, NUMU
      REAL , intent(in)::     DELM0, FBEAM, OPRIM, PI, XR0, XR1
!c     ..
!c     .. Array Arguments ..

      REAL, intent(inout) ::     CWT( MXCMU ), GL( 0:MXCMU ), PSI0( MXCMU ), &
               PSI1( MXCMU ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ), &
               YLMU( 0:MXCMU, * ), Z0( MXCMU ), Z0U( * ), Z1( MXCMU ), &
               Z1U( * ), ZBEAM( * ), ZJ( MXCMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, IU, JQ
      REAL      FACT, PSUM, PSUM0, PSUM1, SUM, SUM0, SUM1
!c     ..


      IF( FBEAM.GT.0.0 ) THEN
!c                                  ** Beam source terms; Eq. SD(9)

         DO 20 IQ = MAZIM, NSTR - 1

            PSUM   = 0.
            DO 10 JQ = 1, NSTR
               PSUM  = PSUM + CWT( JQ )*YLMC( IQ, JQ )*ZJ( JQ )
   10       CONTINUE

            PSI0( IQ + 1 ) = 0.5*GL( IQ )*PSUM

   20    CONTINUE

         FACT   = ( 2. - DELM0 )*FBEAM / ( 4.0*PI )

			IU = 1
!c         DO 40 IU = 1, NUMU

            SUM    = 0.
            DO 30 IQ = MAZIM, NSTR - 1
               SUM  = SUM + YLMU( IQ, IU )* &
                         ( PSI0( IQ+1 ) + FACT*GL( IQ )*YLM0( IQ ) )
   30       CONTINUE

            ZBEAM( IU ) = SUM

!c   40    CONTINUE

      END IF


      IF( PLANK .AND. MAZIM.EQ.0 ) THEN

!c                          ** Thermal source terms, STWJ(27c), STWL(31c)
!c

         DO 60 IQ = MAZIM, NSTR - 1

           PSUM0  = 0.0
            PSUM1  = 0.0
            DO 50 JQ = 1, NSTR
               PSUM0  = PSUM0 + CWT( JQ )*YLMC( IQ, JQ )*Z0( JQ )
               PSUM1  = PSUM1 + CWT( JQ )*YLMC( IQ, JQ )*Z1( JQ )
   50       CONTINUE

            PSI0( IQ + 1 ) = 0.5*GL( IQ ) * PSUM0
            PSI1( IQ + 1 ) = 0.5*GL( IQ ) * PSUM1

   60    CONTINUE

			IU = 1
!c         DO 80 IU = 1, NUMU

            SUM0   = 0.0
            SUM1   = 0.0
            DO 70 IQ = MAZIM, NSTR - 1
               SUM0  = SUM0 + YLMU( IQ, IU ) * PSI0( IQ + 1 )
               SUM1  = SUM1 + YLMU( IQ, IU ) * PSI1( IQ + 1 )
   70       CONTINUE

            Z0U( IU ) = SUM0 + ( 1. - OPRIM ) * XR0
            Z1U( IU ) = SUM1 + ( 1. - OPRIM ) * XR1

!c   80    CONTINUE

      END IF

      END SUBROUTINE TERPSO



      SUBROUTINE UPBEAM( ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM, &
                        MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ, &
                        ZZ )

!c         Finds the incident-beam particular solution of SS(18),
!c         STWL(24a)
!c
!c   I N P U T    V A R I A B L E S:
!c
!c       CC     :  C-sub-ij in Eq. SS(5)
!c
!c       CMU    :  Abscissae for Gauss quadrature over angle cosine
!c
!c       DELM0  :  Kronecker delta, delta-sub-m0
!c
!c       GL     :  Delta-M scaled Legendre coefficients of phase function
!c                 (including factors 2L+1 and single-scatter albedo)
!c
!c       MAZIM  :  Order of azimuthal component
!c
!c       YLM0   :  Normalized associated Legendre polynomial
!c                 at the beam angle
!c
!c       YLMC   :  Normalized associated Legendre polynomial
!c                 at the quadrature angles
!c
!c       (remainder are DISORT input variables)
!c
!c
!c   O U T P U T    V A R I A B L E S:
!c
!c       ZJ     :  Right-hand side vector X-sub-zero in SS(19),STWL(24b);
!c                 also the solution vector Z-sub-zero after solving
!c                 that system
!c
!c       ZZ     :  Permanent storage for ZJ, but re-ordered
!c
!c
!c   I N T E R N A L    V A R I A B L E S:
!c
!c       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19),
!c                   STWL(24b)
!c       IPVT   :  Integer vector of pivot indices required by LINPACK
!c       WK     :  Scratch array required by LINPACK
!c
!c   Called by- DISORT
!c   Calls- SGECO, SGESL
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      INTEGER, intent(in) ::   MAZIM, MXCMU, NN, NSTR
      REAL, intent(in) ::       DELM0, FBEAM, PI, UMU0
!c     ..
!c     .. Array Arguments ..

      INTEGER , intent(inout) ::   IPVT( * )
      REAL, intent(inout) ::       ARRAY( MXCMU, MXCMU ), CC( MXCMU, MXCMU ), CMU( MXCMU ), &
               GL( 0:MXCMU ), WK( MXCMU ), YLM0( 0:MXCMU ), &
               YLMC( 0:MXCMU, * ), ZJ( MXCMU ), ZZ( MXCMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, JOB, JQ, K
      REAL      RCOND, SUM
!c     ..
!c     .. External Subroutines ..

      EXTERNAL  SGECO, SGESL
!c     ..


      DO 30 IQ = 1, NSTR

         DO 10 JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
   10    CONTINUE

         ARRAY( IQ, IQ ) = 1.+ CMU( IQ ) / UMU0 + ARRAY( IQ, IQ )

         SUM  = 0.
         DO 20 K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
   20    CONTINUE

         ZJ( IQ ) = ( 2.- DELM0 )*FBEAM*SUM / ( 4.*PI )
   30 CONTINUE


!c                  ** Find L-U (lower/upper triangular) decomposition
!c                  ** of ARRAY and see if it is nearly singular
!c                  ** (NOTE:  ARRAY is altered)
      RCOND  = 0.0

      CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )


!c                ** Solve linear system with coeff matrix ARRAY
!c                ** (assumed already L-U decomposed) and R.H. side(s)
!c                ** ZJ;  return solution(s) in ZJ
      JOB  = 0

      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB )


      DO 40 IQ = 1, NN
         ZZ( IQ + NN )     = ZJ( IQ )
         ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
   40 CONTINUE

      END SUBROUTINE UPBEAM


      SUBROUTINE UPISOT( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR, OPRIM, &
                        WK, XR0, XR1, Z0, Z1, ZPLK0, ZPLK1 )

!c       Finds the particular solution of thermal radiation of STWL(25)
!c
!c
!c
!c    I N P U T     V A R I A B L E S:
!c
!c       CC     :  C-sub-ij in Eq. SS(5), STWL(8b)
!c
!c       CMU    :  Abscissae for Gauss quadrature over angle cosine
!c
!c       OPRIM  :  Delta-M scaled single scattering albedo
!c
!c       XR0    :  Expansion coefficient b-sub-zero of thermal source
!c                   function, Eq. STWL(24c)
!c
!c       XR1    :  Expansion coefficient b-sub-one of thermal source
!c                   function Eq. STWL(24c)
!c
!c       (remainder are DISORT input variables)
!c
!c
!c    O U T P U T    V A R I A B L E S:
!c
!c       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16), STWL(26a)
!c
!c       Z1     :  Solution vectors Z-sub-one  of Eq. SS(16), STWL(26b)
!c
!c       ZPLK0, :  Permanent storage for Z0,Z1, but re-ordered
!c        ZPLK1
!c
!c
!c   I N T E R N A L    V A R I A B L E S:
!c
!c       ARRAY  :  Coefficient matrix in left-hand side of EQ. SS(16)
!c       IPVT   :  Integer vector of pivot indices required by LINPACK
!c       WK     :  Scratch array required by LINPACK
!c
!c   Called by- DISORT
!c   Calls- SGECO, ERRMSG, SGESL
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      INTEGER   MXCMU, NN, NSTR
      REAL      OPRIM, XR0, XR1
!c     ..
!c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ARRAY( MXCMU, MXCMU ), CC( MXCMU, MXCMU ), CMU( MXCMU ), &
               WK( MXCMU ), Z0( MXCMU ), Z1( MXCMU ), ZPLK0( MXCMU ), &
               ZPLK1( MXCMU )
!c     ..
!c     .. Local Scalars ..

      INTEGER   IQ, JQ
      REAL      RCOND
!c     ..
!c     .. External Subroutines ..

      EXTERNAL  SGECO, SGESL
!c     ..
	  
      DO 20 IQ = 1, NSTR

         DO 10 JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
   10    CONTINUE

         ARRAY( IQ, IQ ) = 1.0 + ARRAY( IQ, IQ )

         Z1( IQ ) = ( 1. - OPRIM ) * XR1

   20 CONTINUE
   
!c                       ** Solve linear equations: same as in UPBEAM,
!c                       ** except ZJ replaced by Z1 and Z0
      RCOND  = 0.0

      CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )

      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, Z1, 0 )

      DO 30 IQ = 1, NSTR
         Z0( IQ ) = ( 1. - OPRIM ) * XR0 + CMU( IQ ) * Z1( IQ )
   30 CONTINUE

      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, Z0, 0 )

      DO 40 IQ = 1, NN
		  
         ZPLK0( IQ + NN ) = Z0( IQ )
         ZPLK1( IQ + NN ) = Z1( IQ )
         ZPLK0( NN + 1 - IQ ) = Z0( IQ + NN )
         ZPLK1( NN + 1 - IQ ) = Z1( IQ + NN )
   40 CONTINUE
		
      END SUBROUTINE UPISOT



      SUBROUTINE USRINT( BPLANK, CMU, CWT, DELM0, DTAUCP, EMU, EXPBEA, &
                        FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL, &
                        LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT, NLYR, &
                        NN, NSTR, PLANK, NUMU, NTAU, PI, RMU, TAUCPR, &
                        TPLANK, UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U, &
                        ZZ, ZPLK0, ZPLK1, UUM )

!c       Computes intensity components at user output angles
!c       for azimuthal expansion terms in Eq. SD(2), STWL(6)
!c
!c
!c   I N P U T    V A R I A B L E S:
!c
!c       BPLANK :  Integrated Planck function for emission from
!c                 bottom boundary
!c
!c       CMU    :  Abscissae for Gauss quadrature over angle cosine
!c
!c       CWT    :  Weights for Gauss quadrature over angle cosine
!c
!c       DELM0  :  Kronecker delta, delta-sub-M0
!c
!c       EMU    :  Surface directional emissivity (user angles)
!c
!c       EXPBEA :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
!c
!c       GC     :  Eigenvectors at polar quadrature angles, SC(1)
!c
!c       GU     :  Eigenvectors interpolated to user polar angles
!c                    (i.e., G in Eq. SC(1) )
!c
!c       KK     :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
!c
!c       LAYRU  :  Layer number of user level UTAU
!c
!c       LL     :  Constants of integration in Eq. SC(1), obtained
!c                 by solving scaled version of Eq. SC(5);
!c                 exponential term of Eq. SC(12) not included
!c
!c       LYRCUT :  Logical flag for truncation of computational layer
!c
!c       MAZIM  :  Order of azimuthal component
!c
!c       NCUT   :  Total number of computational layers considered
!c
!c       NN     :  Order of double-Gauss quadrature (NSTR/2)
!c
!c       RMU    :  Surface bidirectional reflectivity (user angles)
!c
!c       TAUCPR :  Cumulative optical depth (delta-M-Scaled)
!c
!c       TPLANK :  Integrated Planck function for emission from
!c                 top boundary
!c
!c       UTAUPR :  Optical depths of user output levels in delta-M
!c                 coordinates;  equal to UTAU if no delta-M
!c
!c       Z0U    :  Z-sub-zero in Eq. SS(16) interpolated to user
!c                 angles from an equation derived from SS(16),
!c                 Y-sub-zero on STWL(26b)
!c
!c       Z1U    :  Z-sub-one in Eq. SS(16) interpolated to user
!c                 angles from an equation derived from SS(16),
!c                 Y-sub-one in STWL(26a)
!c
!c       ZZ     :  Beam source vectors in Eq. SS(19), STWL(24b)
!c
!c       ZPLK0  :  Thermal source vectors Z0, by solving Eq. SS(16),
!c                 Y-sub-zero in STWL(26)
!c
!c       ZPLK1  :  Thermal source vectors Z1, by solving Eq. SS(16),
!c                 Y-sub-one in STWL(26)
!c
!c       ZBEAM  :  Incident-beam source vectors
!c
!c       (Remainder are DISORT input variables)
!c
!c
!c    O U T P U T    V A R I A B L E S:
!c
!c       UUM    :  Azimuthal components of the intensity in EQ. STWJ(5),
!c                 STWL(6)
!c
!c
!c    I N T E R N A L    V A R I A B L E S:
!c
!c       BNDDIR :  Direct intensity down at the bottom boundary
!c       BNDDFU :  Diffuse intensity down at the bottom boundary
!c       BNDINT :  Intensity attenuated at both boundaries, STWJ(25-6)
!c       DTAU   :  Optical depth of a computational layer
!c       LYREND :  End layer of integration
!c       LYRSTR :  Start layer of integration
!c       PALINT :  Intensity component from parallel beam
!c       PLKINT :  Intensity component from planck source
!c       WK     :  Scratch vector for saving EXP evaluations
!c
!c       All the exponential factors ( EXP1, EXPN,... etc.)
!c       come from the substitution of constants of integration in
!c       Eq. SC(12) into Eqs. S1(8-9).  They all have negative
!c       arguments so there should never be overflow problems.
!c
!c   Called by- DISORT
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      LOGICAL, intent(in) ::   LAMBER, LYRCUT, PLANK
      INTEGER, intent(in) ::   MAZIM, MXCMU, MXULV, MXUMU, NCUT, NLYR, NN, NSTR, NTAU, NUMU
      REAL, intent(in) ::      BPLANK, DELM0, FBEAM, FISOT, PI, TPLANK, UMU0
!c     ..
!c     .. Array Arguments ..

      INTEGER, intent(inout) ::   LAYRU( * )
      REAL, intent(inout) ::  CMU( MXCMU ), CWT( MXCMU ), DTAUCP( * ), EMU( MXUMU ), &
               EXPBEA( 0:* ), GC( MXCMU, MXCMU, * ), &
               GU( MXUMU, MXCMU, * ), KK( MXCMU, * ), LL( MXCMU, * ), &
               RMU( MXUMU, 0:* ), TAUCPR( 0:* ), UMU( * ), &
               UTAUPR( MXULV ), UUM( MXUMU, MXULV ), WK( MXCMU ), &
               Z0U( MXUMU, * ), Z1U( MXUMU, * ), ZBEAM( MXUMU, * ), &
               ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ), ZZ( MXCMU, * )
!c     ..
!c     .. Local Scalars ..

      LOGICAL   NEGUMU
      INTEGER   IQ, IU, JQ, LC, LU, LYREND, LYRSTR, LYU
      REAL      BNDDFU, BNDDIR, BNDINT, DENOM, DFUINT, DTAU, DTAU1, &
               DTAU2, EXP0, EXP1, EXP2, EXPN, F0N, F1N, FACT, PALINT, &
               PLKINT, SGN

!c                          ** Incorporate constants of integration into
!c                          ** interpolated eigenvectors

	  EXP0 = 0.0
	  EXP1=0.0
	  EXP2=0.0
      DO 30 LC = 1, NCUT

         DO 20 IQ = 1, NSTR

			IU = 1
               GU( IU, IQ, LC ) = GU( IU, IQ, LC ) * LL( IQ, LC )

   20    CONTINUE

   30 CONTINUE
!c                           ** Loop over levels at which intensities
!c                           ** are desired ('user output levels')
      DO 160 LU = 1, NTAU

         IF( FBEAM.GT.0.0 ) EXP0  = EXP( -UTAUPR( LU ) / UMU0 )
         LYU  = LAYRU( LU )
!c                              ** Loop over polar angles at which
!c                              ** intensities are desired
			IU = 1

            IF( LYRCUT .AND. LYU.GT.NCUT ) GO TO  150


               LYRSTR = LYU + 1
               LYREND = NCUT
               SGN    = 1.0

!c                          ** For downward intensity, integrate from top
!c                          ** to LYU-1 in Eq. S1(8); for upward,
!c                          ** integrate from bottom to LYU+1 in S1(9)
            PALINT = 0.0
            PLKINT = 0.0

            DO 60 LC = LYRSTR, LYREND

               DTAU = DTAUCP( LC )
               EXP1 = EXP( ( UTAUPR(LU) - TAUCPR(LC-1) ) / UMU( IU ) )
               EXP2 = EXP( ( UTAUPR(LU) - TAUCPR(LC)   ) / UMU( IU ) )

               IF( PLANK .AND. MAZIM.EQ.0 ) THEN

!c                          ** Eqs. STWL(36b,c, 37b,c)
!c
                  F0N = SGN * ( EXP1 - EXP2 )

                  F1N = SGN * ( ( TAUCPR( LC-1 ) + UMU( IU ) ) * EXP1 - &
                               ( TAUCPR( LC )   + UMU( IU ) ) * EXP2 )

                  PLKINT = PLKINT + Z0U( IU,LC )*F0N + Z1U( IU,LC )*F1N

               END IF


               IF( FBEAM.GT.0.0 ) THEN

                  DENOM  = 1. + UMU( IU ) / UMU0

                  IF( ABS( DENOM ).LT.0.0001 ) THEN
!c                                                   ** L'Hospital limit
                     EXPN   = ( DTAU / UMU0 )*EXP0

                  ELSE

                     EXPN   = ( EXP1*EXPBEA( LC-1 ) - &
                               EXP2*EXPBEA( LC ) ) * SGN / DENOM

                  END IF

                  PALINT = PALINT + ZBEAM( IU, LC )*EXPN

               END IF

!c                                                   ** KK is negative
               DO 40 IQ = 1, NN

                  WK( IQ ) = EXP( KK( IQ,LC )*DTAU )
                  DENOM  = 1.0 + UMU( IU )*KK( IQ, LC )

                  IF( ABS( DENOM ).LT.0.0001 ) THEN
!c                                                   ** L'Hospital limit
                     EXPN   = DTAU / UMU( IU )*EXP2

                  ELSE

                     EXPN   = SGN*( EXP1*WK( IQ ) - EXP2 ) / DENOM

                  END IF

                  PALINT = PALINT + GU( IU, IQ, LC )*EXPN

   40          CONTINUE

!c                                                   ** KK is positive
               DO 50 IQ = NN + 1, NSTR

                  DENOM  = 1.0 + UMU( IU )*KK( IQ, LC )

                  IF( ABS( DENOM ).LT.0.0001 ) THEN
!c                                                   ** L'Hospital limit
                     EXPN  = -DTAU / UMU( IU )*EXP1

                  ELSE

                     EXPN  = SGN*( EXP1 - EXP2*WK( NSTR+1-IQ ) ) / DENOM

                  END IF

                  PALINT = PALINT + GU( IU, IQ, LC )*EXPN

   50          CONTINUE


   60       CONTINUE
!c                           ** Calculate contribution from user
!c                           ** output level to next computational level

            DTAU1  = UTAUPR( LU ) - TAUCPR( LYU - 1 )
            DTAU2  = UTAUPR( LU ) - TAUCPR( LYU )

            IF( ABS( DTAU2 ).LT.1.E-6  ) GO TO  90

			EXP2  = EXP( DTAU2/UMU( IU ) )

            IF( FBEAM.GT.0.0 ) THEN

               DENOM  = 1. + UMU( IU ) / UMU0

               IF( ABS( DENOM ).LT.0.0001 ) THEN

                  EXPN   = ( DTAU1 / UMU0 )*EXP0

               ELSE

                  EXPN  = ( EXP0 - EXPBEA( LYU )*EXP2 ) / DENOM

              END IF

               PALINT = PALINT + ZBEAM( IU, LYU )*EXPN

            END IF

!c                                                   ** KK is negative
            DTAU  = DTAUCP( LYU )

            DO 70 IQ = 1, NN

               DENOM  = 1. + UMU( IU )*KK( IQ, LYU )

               IF( ABS( DENOM ).LT.0.0001 ) THEN

                  EXPN = -DTAU2 / UMU( IU )*EXP2

               ELSE

                  EXPN = ( EXP( -KK( IQ,LYU ) * DTAU2 ) - EXP2 ) / DENOM

               END IF

               PALINT = PALINT + GU( IU, IQ, LYU )*EXPN

   70       CONTINUE

!c                                                   ** KK is positive
            DO 80 IQ = NN + 1, NSTR

               DENOM  = 1. + UMU( IU )*KK( IQ, LYU )

               IF( ABS( DENOM ).LT.0.0001 ) THEN

                  EXPN   = -DTAU1 / UMU( IU )*EXP1

               ELSE

                  EXPN = ( EXP( -KK( IQ,LYU ) * DTAU1 ) - &
	                     EXP( -KK( IQ,LYU ) * DTAU  ) * EXP2 ) / DENOM

               END IF

               PALINT = PALINT + GU( IU, IQ, LYU )*EXPN

   80       CONTINUE


            IF( PLANK .AND. MAZIM.EQ.0 ) THEN

!c                            ** Eqs. STWL (35-37) with tau-sub-n-1
!c                            ** replaced by tau for upward, and
!c                            ** tau-sub-n replaced by tau for downward
!c                            ** directions


                  EXPN  = EXP2
                  FACT  = TAUCPR( LYU ) + UMU( IU )


               F0N  = 1. - EXPN
               F1N  = UTAUPR( LU ) + UMU( IU ) - FACT * EXPN

               PLKINT = PLKINT + Z0U( IU, LYU )*F0N + Z1U( IU, LYU )*F1N

            END IF

!c                            ** Calculate intensity components
!c                            ** attenuated at both boundaries.
!c                            ** NOTE: no azimuthal intensity
!c                            ** component for isotropic surface
   90       CONTINUE
            BNDINT = 0.0


               IF( LYRCUT .OR. ( LAMBER.AND.MAZIM.GT.0 ) ) GO TO  140

               DO 100 JQ = NN + 1, NSTR
                  WK( JQ ) = EXP( -KK( JQ,NLYR )*DTAUCP( NLYR ) )
  100          CONTINUE

               BNDDFU = 0.0

               DO 130 IQ = NN, 1, -1

                  DFUINT = 0.0
                  DO 110 JQ = 1, NN
                     DFUINT = DFUINT + GC( IQ, JQ, NLYR )*LL( JQ, NLYR )
  110             CONTINUE

                  DO 120 JQ = NN + 1, NSTR
                     DFUINT = DFUINT + GC( IQ, JQ, NLYR )* &
                                      LL( JQ, NLYR )*WK( JQ )
  120             CONTINUE

                  IF( FBEAM.GT.0.0 ) DFUINT = DFUINT + &
                                          ZZ( IQ, NLYR )*EXPBEA( NLYR )

                  DFUINT = DFUINT + DELM0 * ( ZPLK0( IQ, NLYR ) + &
                                   ZPLK1( IQ,NLYR ) *TAUCPR( NLYR ) )
                  BNDDFU = BNDDFU + ( 1.+DELM0 ) * RMU(IU,NN+1-IQ) &
                                 * CMU(NN+1-IQ) * CWT(NN+1-IQ)* DFUINT
  130          CONTINUE

               BNDDIR = 0.0
               IF( FBEAM.GT.0.0 ) BNDDIR = UMU0*FBEAM / PI*RMU( IU, 0 )* &
                                          EXPBEA( NLYR )

               BNDINT = ( BNDDFU + BNDDIR + DELM0 * EMU(IU) * BPLANK ) &
                       * EXP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )


  140       CONTINUE

            UUM( IU, LU ) = PALINT + PLKINT + BNDINT

  150    CONTINUE

  160 CONTINUE


      END SUBROUTINE USRINT




!c ******************************************************************
!c ********** DISORT service routines ************************
!c ******************************************************************

     
      SUBROUTINE LEPOLY( NMU, M, MAXMU, TWONM1, MU, SQT, YLM )

!c       Computes the normalized associated Legendre polynomial,
!c       defined in terms of the associated Legendre polynomial
!c       Plm = P-sub-l-super-m as
!c
!c             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)
!c
!c       for fixed order m and all degrees from l = m to TWONM1.
!c       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
!c       from a prior call to the routine.
!c
!c       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!c                  High-Order Associated Legendre Polynomials,
!c                  J. Quant. Spectrosc. Radiat. Transfer 10,
!c                  557-562, 1970.  (hereafter D/A)
!c
!c       METHOD: Varying degree recurrence relationship.
!c
!c       NOTES:
!c       (1) The D/A formulas are transformed by setting M=n-1; L=k-1.
!c       (2) Assumes that routine is called first with  M = 0, then with
!c           M = 1, etc. up to  M = TWONM1.
!c
!c
!c  I N P U T     V A R I A B L E S:
!c
!c       NMU    :  Number of arguments of YLM
!c
!c       M      :  Order of YLM
!c
!c       MAXMU  :  First dimension of YLM
!c
!c       TWONM1 :  Max degree of YLM
!c
!c       MU(i)  :  Arguments of YLM (i = 1 to NMU)
!c
!c       SQT(k) :  Square root of k
!c
!c       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
!c       from a prior call.
!c
!c
!c  O U T P U T     V A R I A B L E:
!c
!c       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
!c                   polynomials evaluated at argument MU(i)
!c
!c   Called by- DISORT, ALBTRN
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      INTEGER, intent(in) :: M, MAXMU, NMU, TWONM1
!c     ..
!c     .. Array Arguments ..

      REAL, intent(inout) :: MU( * ), YLM( 0:MAXMU, * ), SQT( * )
!c     ..
!c     .. Local Scalars ..

      INTEGER   I, L
      REAL      TMP1, TMP2
!c     ..

      IF( M.EQ.0 ) THEN
!c                             ** Upward recurrence for ordinary
!c                             ** Legendre polynomials
         DO 20 I = 1, NMU
            YLM( 0, I ) = 1.0
            YLM( 1, I ) = MU( I )
   20    CONTINUE


         DO 40 L = 2, TWONM1

            DO 30 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) - &
                              ( L - 1 )*YLM( L-2, I ) ) / L
   30       CONTINUE

   40    CONTINUE


      ELSE

         DO 50 I = 1, NMU
!c                               ** Y-sub-m-super-m; derived from
!c                               ** D/A Eqs. (11,12), STWL(58c)

            YLM( M, I ) = - SQT( 2*M - 1 ) / SQT( 2*M )* &
                           SQRT( 1.- MU(I)**2 )*YLM( M-1, I )

!c                              ** Y-sub-(m+1)-super-m; derived from
!c                              ** D/A Eqs.(13,14) using Eqs.(11,12),
!c                              ** STWL(58f)

            YLM( M+1, I ) = SQT( 2*M + 1 )*MU( I )*YLM( M, I )

   50    CONTINUE

!c                                   ** Upward recurrence; D/A EQ.(10),
!c                                   ** STWL(58a)
         DO 70 L = M + 2, TWONM1

            TMP1  = SQT( L - M )*SQT( L + M )
            TMP2  = SQT( L - M - 1 )*SQT( L + M - 1 )

            DO 60 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) - &
                              TMP2*YLM( L-2, I ) ) / TMP1
   60       CONTINUE

   70    CONTINUE

      END IF

      END subroutine LEPOLY



      REAL FUNCTION PLKAVG( WNUMLO, WNUMHI, T )

!c        Computes Planck function integrated between two wavenumbers
!c
!c  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
!c
!c           WNUMHI : Upper wavenumber
!c
!c           T      : Temperature (K)
!c
!c  OUTPUT : PLKAVG : Integrated Planck function ( Watts/sq m )
!c                      = Integral (WNUMLO to WNUMHI) of
!c                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
!c                        (where h=Plancks constant, c=speed of
!c                         light, nu=wavenumber, T=temperature,
!c                         and k = Boltzmann constant)
!c
!c  Reference : Specifications of the Physical World: New Value
!c                 of the Fundamental Constants, Dimensions/N.B.S.,
!c                 Jan. 1974
!c

	  REAL, parameter :: C2 = 1.438786 , SIGMA = 5.67032E-8 , &
				VCP(7) =  (/ 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /)
	  
      REAL, intent(in) ::    T, WNUMHI, WNUMLO

      INTEGER   I, MMAX, M 
      REAL      EX, EXM, MV


      REAL      D( 2 ), V( 2 )

      V( 1 ) = C2*WNUMLO / T
      V( 2 ) = C2*WNUMHI / T

! the code for close wavenumbers has been deleted because in our particular case they are never that close

      DO 60 I = 1, 2

!c                      ** Use exponential series
            MMAX  = 0
!c                                ** Find upper limit of series
   40       CONTINUE
            MMAX  = MMAX + 1

            IF( V(I) .LT. VCP( MMAX ) ) GO TO  40

            EX     = EXP( - V(I) )
            EXM    = 1.0
            D( I ) = 0.0

            DO 50 M = 1, MMAX
               MV     = M*V( I )
               EXM    = EX*EXM
               D( I ) = D( I ) + EXM*( 6.+ MV*( 6.+ MV*( 3.+ MV ) ) ) / M**4
   50       CONTINUE

            D( I ) = PLKAVG_CONC*D( I )


   60 CONTINUE

      PLKAVG = PLKAVG_SIGDPI * T**4 * (D( 1 ) - D( 2 ))

      END function PLKAVG




      SUBROUTINE QGAUSN( M, GMU, GWT )

!c       Compute weights and abscissae for ordinary Gaussian quadrature
!c       on the interval (0,1);  that is, such that

!c           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

!c       is a good approximation to

!c           integral(0 to 1) ( f(x) dx )

!c   INPUT :    M       order of quadrature rule

!c   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
!c             GWT(I)   array of weights (I = 1 TO M)

!c   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!c                   Integration, Academic Press, New York, pp. 87, 1975

!c   METHOD:  Compute the abscissae as roots of the Legendre
!c            polynomial P-sub-M using a cubically convergent
!c            refinement of Newton's method.  Compute the
!c            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!c            that Newton's method can very easily diverge; only a
!c            very good initial guess can guarantee convergence.
!c            The initial guess used here has never led to divergence
!c            even for M up to 1000.

!c   ACCURACY:  relative error no better than TOL or computer
!c              precision (machine epsilon), whichever is larger

!c   INTERNAL VARIABLES:

!c    ITER      : number of Newton Method iterations
!c    MAXIT     : maximum allowed iterations of Newton Method
!c    PM2,PM1,P : 3 successive Legendre polynomials
!c    PPR       : derivative of Legendre polynomial
!c    P2PRI     : 2nd derivative of Legendre polynomial
!c    TOL       : convergence criterion for Legendre poly root iteration
!c    X,XI      : successive iterates in cubically-convergent version
!c                of Newtons Method (seeking roots of Legendre poly.)

!c   Called by- DREF, SETDIS, SURFAC
!c   Calls- D1MACH, ERRMSG
!c +-------------------------------------------------------------------+

!c     .. Scalar Arguments ..

      INTEGER, intent(in) ::  M
      REAL, intent(inout) ::  GMU(:), GWT(:)


      INTEGER   ITER, K, LIM, NN, NP1
      REAL      CONA, T
      DOUBLE PRECISION EN, NNP1, P, P2PRI, PM1, PM2, PPR, PROD, &
                      TMP, X, XI
		double precision :: tan_t, xang

	  P=0.0d0
	  

      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = FLOAT( M - 1 ) / ( 8*M**3 )

      LIM  = M / 2

      DO K = 1, LIM
!c                                        ** Initial guess for k-th root
!c                                        ** of Legendre polynomial, from
!c                                        ** Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
			if (t < approx_limit_tan) then 
				tan_t = t
			else
				tan_t = TAN(T)
			endif
			xang = T + CONA / tan_t
			if (xang < approx_limit) then 
				x = 1. - xang*xang/2.
			else
				x = cos(xang)
			endif
			
!         X  = COS( T + CONA / TAN( T ) )
         ITER = 0

!c                                        ** Upward recurrence for
!c                                        ** Legendre polynomials
   10    CONTINUE
         ITER   = ITER + 1
         PM2    = ONE
         PM1    = X

         DO 20 NN = 2, M
            P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
            PM2  = PM1
            PM1  = P
   20    CONTINUE
!c                                              ** Newton Method
         TMP    = ONE / ( ONE - X**2 )
         PPR    = EN*( PM2 - X*P )*TMP
         P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
         XI     = X - ( P / PPR )*( ONE + &
                 ( P / PPR )*P2PRI / ( TWO*PPR ) )

!c                                              ** Check for convergence
         IF( ABS( XI - X ).GT.TOL ) THEN

            X  = XI
            GO TO  10

         END IF
		 		 
!c                             ** Iteration finished--calculate weights,
!c                             ** abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
	  end do
!c                                        ** Convert from (-1,1) to (0,1)
      DO K = 1, M
         GMU( K ) = 0.5*GMU( K ) + 0.5
         GWT( K ) = 0.5*GWT( K )
	  end do

      END subroutine QGAUSN


      REAL FUNCTION RATIO( A, B )

!c        Calculate ratio  A/B  with over- and under-flow protection
!c        (thanks to Prof. Jeff Dozier for some suggestions here).
!c        Since this routine takes two logs, it is no speed demon,
!!c      but it is invaluable for comparing results from two runs
!c       of a program under development.
!c
!c       NOTE:  In Fortran90, built-in functions TINY and HUGE
!c               can replace the R1MACH calls.
!c
!c   Called by- DISORT
!c   Calls- R1MACH
!c +-------------------------------------------------------------------+

      REAL, intent(in) :: A, B
      REAL  ::    ABSA, ABSB


      IF( A == 0.0 ) THEN
         IF( B == 0.0 ) THEN
            RATIO  = 1.0
         ELSE
            RATIO  = 0.0
         END IF
      ELSE IF( B == 0.0 ) THEN
         RATIO  = SIGN( HUGE, A )
      ELSE
         ABSA   = ABS( A )
         ABSB   = ABS( B )
!c the chances of over/underflow like described here in our 
!c particular problem are pretty much zero. So...
         IF( ABSA < TINY .AND. ABSB < TINY ) THEN
            RATIO  = 1.0
         ELSE
            RATIO  = ABSA / ABSB
         END IF
!c                      ** DONT use old trick of determining sign
!c                      ** from A*B because A*B may (over/under)flow
         IF( ( A > 0.0 .AND. B < 0.0 ) .OR. &
            ( A < 0.0 .AND. B > 0.0 ) ) RATIO = -RATIO

      END IF
      END function ratio




end module DISORT_solver
