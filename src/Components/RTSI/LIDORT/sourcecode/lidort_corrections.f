C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LIDORT_SSCORR_NADIR (master)                     #
C #            LIDORT_DBCORRECTION (master)                     #
C #                                                             #
C #      Version 3.2 outgoing sphericity correction             #
C #              3.3.  partial-layer integration                #
C #                                                             #
C #             LIDORT_SSCORR_OUTGOING (master)                 #
C #                 outgoing_sphergeom_fine                     #
C #                 outgoing_integration_up                     #
C #                 outgoing_integration_dn                     #
C #                 multi_outgoing_adjustgeom                   #
C #                                                             #
C #      Version 3.3 outgoing Lambertian DB correction          #
C #                                                             #
C #            LIDORT_LAMBERTIAN_DBCORRECTION                   #
C #                                                             #
C ###############################################################
C #                                                             #
C #         Redundant routines................                  #
C #                 outgoing_sphergeom                          #
C #                 outgoing_integration_up_old                 #
C #                 outgoing_integration_dn_old                 #
C #                 outgoing_adjustgeom                         #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_SSCORR_NADIR (SSFLUX)

C  Single scatter exact calculation for the incoming solar beam
C   This is the regular pseudo-spherical calculation only

C   Programmed by R. Spurr, RT Solutions Inc.
C    Second Draft, April 14th 2005.
C    Third Draft,  May    6th 2005.
C       - additional code to deal with refraction.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file for storing geometrical indices

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  Argument
C  --------

      DOUBLE PRECISION SSFLUX

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, L, NM1
      INTEGER          UT, UTA, UM, IUM, IA, NC, NSAVE, IB, V

C  help variables (double precision)

      DOUBLE PRECISION FINAL_SOURCE, HELP, SS_LAYERSOURCE, DNM1, DNL1
      DOUBLE PRECISION MULTIPLIER, SSCORRECTION, VS, MUX, COSSCAT
      DOUBLE PRECISION TR_CUMSOURCE, SS_CUMSOURCE, LEGPOLY, F, FT
      DOUBLE PRECISION DF1(0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION DF2(0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION GK11(0:MAXMOMENTS_INPUT)

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)

      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  Set up operations
C  -----------------

C  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO
      NSAVE = 0

C  Create TMS factors, these get stored
C    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
          TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N)
        ENDIF
      ENDDO

C  Additional Delta-M scaling
C  New section. R. Spurr, 07 September 2007.
C   TMS gets modified by (1-F). Save the truncation factor.
C   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

C  save some geometrical quantities (cosines and sines)
C    Multiple layer quantities for the Refractive case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            MUX =  SUN_SZA_COSINES(N,IB)  
            CTHETA(N,IB) = MUX
            STHETA(N,IB) = DSQRT(ONE-MUX*MUX)
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            CTHETA(N,IB) = X0(IB)
            STHETA(N,IB) = SX0(IB)
          END DO
        END DO
      ENDIF

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

C  ====================================
C  Upwelling single scatter calculation
C  ====================================

      IF ( DO_UPWELLING ) THEN

C      VS = -1 for upwelling, +1 for downwelling

        VS = -1.0D0

C  Upwelling single scatter Phase matrices
C  ---------------------------------------

C  For each geometry (indexed by V), do the following:-
C     Develop the cosine scatter angle
C     Get Legendre polynomials for cos(THETA), save them
C     Form Exact Phase function
C     Scattering result multiplied by TMS factor

C  initialise and start layer loop

        NC = 0
        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N)) THEN
          NC = NC + 1

C  refractive geometry case - Always do the Legendre calculation

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO IA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + IA
              COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) +
     &                       STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
	      SS_PLEG_UP(V,N,0) = ONE
	      SS_PLEG_UP(V,N,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_UP(V,N,L) =
     &           DF1(L) * SS_PLEG_UP(V,N,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_UP(V,N,L-2)
	      ENDDO
             ENDDO
            ENDDO
           ENDDO

C  Non-refractive case, Just do for the first layer, then copy the rest

          ELSE IF ( .NOT. DO_REFRACTIVE_GEOMETRY ) THEN
           IF ( NC.EQ.1) THEN
            NSAVE = N
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              DO IA = 1, N_USER_RELAZMS
               V = UMOFF(IB,UM) + IA
               COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) +
     &                        STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
	       SS_PLEG_UP(V,N,0) = ONE
	       SS_PLEG_UP(V,N,1) = COSSCAT
	       DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_UP(V,N,L) =
     &            DF1(L) * SS_PLEG_UP(V,N,L-1) * COSSCAT  -
     &            DF2(L) * SS_PLEG_UP(V,N,L-2)
	       ENDDO
              ENDDO
             ENDDO
            ENDDO             
           ELSE
            DO V = 1, N_GEOMETRIES
	     DO L = 0, NMOMENTS_INPUT
	      SS_PLEG_UP(V,N,L) = SS_PLEG_UP(V,NSAVE,L)
             ENDDO
            ENDDO
           ENDIF
          ENDIF

C  Local moments

          GK11(0) = ONE
          IF ( DO_SSCORR_TRUNCATION ) THEN
            DO L = 1, NMOMENTS_INPUT
              DNL1 = DBLE(2*L + 1 )
              F    = SSFDEL(N) * DNL1
              FT   = ONE - SSFDEL(N)
              GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N)-F)/FT
            ENDDO
          ELSE
            DO L = 1, NMOMENTS_INPUT
              GK11(L) = PHASMOMS_TOTAL_INPUT(L,N)
            ENDDO
          ENDIF

C  Get the total phase function

          DO V = 1, N_GEOMETRIES   
            HELP = ZERO
            DO L = 0, NMOMENTS_INPUT
	     LEGPOLY = SS_PLEG_UP(V,N,L)
	     HELP = HELP + GK11(L)*LEGPOLY
	    ENDDO
	    EXACTSCAT_UP(V,N) = HELP * TMS(N)
          ENDDO

C  end layer loop

         ENDIF
        ENDDO

C  Upwelling single scatter recurrence
C  -----------------------------------

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_UP(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_UP(V,N)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE +
     &            T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,NC-1)
                ENDDO
              ENDDO
            ENDDO

C  end layer loop

          ENDDO

C  Offgrid output :
C  ----------------

C    Add additional partial layer source term = Exact Phase Func * Multiplier
C    Get final cumulative source and set the Single scatter results

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
                IUM = USEROUTPUT_INDEX(UM)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_UP(V,N)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
                  TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output :
C  ---------------

C     Set final cumulative source and single scatter intensity

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
                  SSCORRECTION = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  ======================================
C  Downwelling single scatter calculation
C  ======================================

      IF ( DO_DNWELLING ) THEN

C      VS = -1 for upwelling, +1 for downwelling

        VS = +1.0D0

C  Downwelling single scatter Phase matrices
C  -----------------------------------------

C  For each geometry (indexed by V), do the following:-
C     Develop the cosine scatter angle
C     Get Legendre polynomials for cos(THETA), save them
C     Form Exact Phase function
C     Scattering result multiplied by TMS factor

C  initialise and start layer loop

        NC = 0
        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N)) THEN
          NC = NC + 1

C  refractive geometry case - Always do the Legendre calculation

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO IA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + IA
              COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) +
     &                       STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
	      SS_PLEG_DN(V,N,0) = ONE
	      SS_PLEG_DN(V,N,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_DN(V,N,L) =
     &           DF1(L) * SS_PLEG_DN(V,N,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_DN(V,N,L-2)
	      ENDDO
             ENDDO
            ENDDO
           ENDDO

C  Non-refractive case, Just do for the first layer, then copy the rest

          ELSE IF ( .NOT. DO_REFRACTIVE_GEOMETRY ) THEN
           IF ( NC.EQ.1) THEN
            NSAVE = N
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              DO IA = 1, N_USER_RELAZMS
               V = UMOFF(IB,UM) + IA
               COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) +
     &                        STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
	       SS_PLEG_DN(V,N,0) = ONE
	       SS_PLEG_DN(V,N,1) = COSSCAT
	       DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_DN(V,N,L) =
     &            DF1(L) * SS_PLEG_DN(V,N,L-1) * COSSCAT  -
     &            DF2(L) * SS_PLEG_DN(V,N,L-2)
	       ENDDO
              ENDDO
             ENDDO
            ENDDO             
           ELSE
            DO V = 1, N_GEOMETRIES
	     DO L = 0, NMOMENTS_INPUT
	      SS_PLEG_DN(V,N,L) = SS_PLEG_DN(V,NSAVE,L)
             ENDDO
            ENDDO
           ENDIF
          ENDIF

C  Local moments

          IF ( DO_SSCORR_TRUNCATION ) THEN
            GK11(0) = ONE
            DO L = 1, NMOMENTS_INPUT
              DNL1 = DBLE(2*L + 1 )
              F    = SSFDEL(N) * DNL1
              FT   = ONE - SSFDEL(N)
              GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N)-F)/FT
            ENDDO
          ELSE
            DO L = 1, NMOMENTS_INPUT
              GK11(L) = PHASMOMS_TOTAL_INPUT(L,N)
            ENDDO
          ENDIF

C  Get the total phase function

          DO V = 1, N_GEOMETRIES   
           HELP = ZERO
           DO L = 0, NMOMENTS_INPUT
	    LEGPOLY = SS_PLEG_DN(V,N,L)
	    HELP = HELP + GK11(L)*LEGPOLY
	   ENDDO
	   EXACTSCAT_DN(V,N) = HELP * TMS(N)
          ENDDO

C  end layer loop

         ENDIF
        ENDDO

C  Downwelling single scatter recurrence
C  -------------------------------------

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT
           NC = N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_DN(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP =  EXACTSCAT_DN(V,N)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE +
     &             T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,NC-1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct the intensity

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
                IUM = USEROUTPUT_INDEX(UM)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_DN(V,N)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output :
C     Set final cumulative source and correct Stokes vector

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                IUM = USEROUTPUT_INDEX(UM)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
                  SSCORRECTION = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end optical depth loop and Downwelling clause

        ENDDO
      ENDIF

C  debug

c      nut = 97
c      if ( do_fdtest ) nut = 98
c       do v = 1, 8
c        write(nut,*)v,INTENSITY_SS(2,V,upIDX), INTENSITY_SS(2,V,dnIDX)   
c       enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_DBCORRECTION (FLUXMULT)

C  Prepares Exact Direct Beam reflection for the BRDF case

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Input arguments controlling type of surface
C  output arguments also go in here.

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include file for storing geometrical indices

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Output goes in the correction file

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, KL
      DOUBLE PRECISION FINAL_SOURCE, TR
      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN, REFL_ATTN, SUM

C  first stage
C  -----------

C  Initialize

      DO V = 1, N_GEOMETRIES
        EXACTDB_SOURCE(V) = ZERO
        DO UTA = 1, N_OUT_USERTAUS
          INTENSITY_DB(UTA,V)      = ZERO
        ENDDO
      ENDDO

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  Start geoemtry loops

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA

c  Beam attenuation

              IF ( DO_SSCORR_OUTGOING ) THEN
                X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                ATTN = BOA_ATTN(V)
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                ELSE
                  X0_BOA = X0(IB)
                ENDIF 
                ATTN  = SOLAR_BEAM_OPDEP(IB)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              ATTN    = X0_FLUX * ATTN
              ATTN_DB_SAVE(V) = ATTN

C  Start loop over albedo kernels

              DO KL = 1, N_BRDF_KERNELS
                REFL_ATTN = BRDF_FACTORS(KL) * ATTN
                SUM = EXACTDB_BRDFUNC(KL,UM,UA,IB)
                A_EXACTDB_SOURCE(V,KL) = REFL_ATTN * SUM
              ENDDO

C  Finish loops over geometries

            ENDDO
          ENDIF
        ENDDO
      ENDDO

C  Build together the kernels. Multiply result by solar Flux
C    Quantities are saved, because required later in linearization

      DO V = 1, N_GEOMETRIES
        DO KL = 1, N_BRDF_KERNELS
          EXACTDB_SOURCE(V) = EXACTDB_SOURCE(V) +
     &                   A_EXACTDB_SOURCE(V,KL)
        ENDDO
        EXACTDB_SOURCE(V) = EXACTDB_SOURCE(V) * FLUXMULT
      ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

      NC =  0
      DO V = 1, N_GEOMETRIES
        DB_CUMSOURCE(V,NC) = EXACTDB_SOURCE(V)
      ENDDO

C  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS(N,V)
                ELSE
                  TR = T_DELT_USERM(N,UM)
                ENDIF
                DB_CUMSOURCE(V,NC) = TR * DB_CUMSOURCE(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

C  end layer loop

        ENDDO

C  Offgrid output : partial layer transmittance, then set result

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS_UT(UT,V)
                ELSE
                  TR = T_UTUP_USERM(UT,UM)
                ENDIF
                INTENSITY_DB(UTA,V) = TR * DB_CUMSOURCE(V,NC)
              ENDDO
            ENDDO
          ENDDO

C  Ongrid output : Set final cumulative source directly

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                FINAL_SOURCE = DB_CUMSOURCE(V,NC)
                INTENSITY_DB(UTA,V) = FINAL_SOURCE
              ENDDO
            ENDDO
          ENDDO

        ENDIF

C  Check for updating the recursion 

        IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

C  end optical depth loop

      ENDDO

C  debug

c      um = 97
c      if ( do_fdtest ) um = 98
c      do v = 1, 8
c       write(um,*)v, intensity_db(1,v)
c      enddo
c      if ( do_fdtest) pause

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_SSCORR_OUTGOING
     &        ( SSFLUX, FAIL, MESSAGE)

C  Single scatter exact calculation for the outgoing LOS
C         - NEW for Versions 3.2 and 3.3

C   Programmed by R. Spurr, RT Solutions Inc.
C    First Draft, January 23rd 2007
C    Final Draft, October 3rd, 2007

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  Arguments
C  ---------

      DOUBLE PRECISION SSFLUX
      logical          fail
      character*(*)    message

C  Geometry routine inputs and outputs
C  -----------------------------------

c  control

      logical          do_fine
      LOGICAL          do_partials
      double precision alpha_boa, theta_boa, phi_boa

C  main outputs (geometry)

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level output (geometry)

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision radii_fine    (maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layer output (geometry)

      integer          ntraverse_ut(max_offgrid_usertaus)
      double precision sunpaths_ut (maxlayers,max_offgrid_usertaus)
      double precision radii_ut    (max_offgrid_usertaus)
      double precision alpha_ut    (max_offgrid_usertaus)

C  Other (incidental) geometrical output

      double precision lospaths(maxlayers)
      double precision lospaths_ut_up(max_offgrid_usertaus)
      double precision lospaths_ut_dn(max_offgrid_usertaus)
      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

C  Extinction

      double precision extinction ( maxlayers)

C  Partial layer heights

      double precision height_grid_ut(max_offgrid_usertaus)

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, L, NM1
      INTEGER          UT, UTA, UM, IA, NC, IB, V

C  help variables (double precision)

      DOUBLE PRECISION FINAL_SOURCE, HELP, SS_LAYERSOURCE, DNM1, DNL1
      DOUBLE PRECISION SSCORRECTION, COSSCAT, LEGPOLY, SS_CUMSOURCE
      DOUBLE PRECISION DF1(0:MAXMOMENTS_INPUT), XT
      DOUBLE PRECISION DF2(0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION GK11(0:MAXMOMENTS_INPUT), TRANS, F, FT

C  Set up operations
C  -----------------

C  Set up partials flag

      do_fine = .true.
      do_partials = ( n_offgrid_usertaus .gt. 0 )

C  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

C  Create TMS factors, these get stored
C    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
          TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N)
        ENDIF
      ENDDO

C  Create extinctions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

C  Create the heights of partial layers

      IF ( DO_PARTIALS ) THEN
        do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          xt = deltau_vert(n) - offgrid_utau_values(ut)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.
C   TMS gets modified by (1-F). Save the truncation factor.
C   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

C  Source term calculations
C  ========================

C  Start the main loop over all solar and viewing geometries
C   Use the adjusted values of the angles

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = BEAM_SZAS_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = UMOFF(IB,UM) + IA
    
C  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine
     i  ( maxlayers, maxfinelayers, max_offgrid_usertaus,
     i    do_fine, do_partials, nlayers, nfinelayers,
     i    n_offgrid_usertaus, offgrid_utau_layeridx,
     i    height_grid, height_grid_ut, earth_radius,
     i    alpha_boa, theta_boa, phi_boa,
     o    sunpaths,      radii,      ntraverse,      alpha_all, 
     o    sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,
     o    sunpaths_ut,   radii_ut,   ntraverse_ut,   alpha_ut,
     o    offgrid_utau_layerfineidx, lospaths,
     o    lospaths_ut_up, lospaths_ut_dn, 
     o    theta_all, phi_all, cosscat_up, cosscat_dn,
     o    fail, message )

C  debug geometry
c            do n = 0, nlayers
c              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
c            enddo
c            pause

C  Upwelling Source terms: Phase matrix, multipliers, transmittances
C  -----------------------------------------------------------------

            IF ( DO_UPWELLING ) THEN

C  Multipliers, transmittances

              call outgoing_integration_up
     i           ( nlayers, nfinelayers,
     i             do_fine, do_partials, extinction,
     i             n_offgrid_usertaus, offgrid_utau_layeridx,
     i             offgrid_utau_layerfineidx,
     i             sunpaths, radii, ntraverse, alpha_all,
     i             sunpaths_ut,   ntraverse_ut,   alpha_ut,
     i             sunpaths_fine, ntraverse_fine, alpha_fine,
     o             up_multipliers(1,v), up_lostrans(1,v), boa_attn(v),
     o             up_multipliers_ut(1,v), up_lostrans_ut(1,v) )

C  Debug
c              write(45,*)v
c              do n = 1, nlayers
c                write(45,'(1p2e18.10)')
c     &           up_lostrans(n,v),up_multipliers(n,v)
c              enddo
C  Debug: Multipliers (all layers), Upwelling OLD WAY
c              call outgoing_integration_old_up
c     i           ( nlayers, extinction, deltau_vert,
c     i             lospaths, sunpaths, radii, ntraverse, alpha_all,
c     o             up_multipliers(1,v), up_lostrans(1,v) )
c              write(46,*)v
c              do n = 1, nlayers
c                write(46,*)up_lostrans(n,v),up_multipliers(n,v)
c              enddo

C  legendre polynomials

              COSSCAT = COSSCAT_UP(NLAYERS)
	      SS_PLEG_UP(V,1,0) = ONE
	      SS_PLEG_UP(V,1,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_UP(V,1,L) =
     &           DF1(L) * SS_PLEG_UP(V,1,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_UP(V,1,L-2)
              ENDDO

C  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(N) ) THEN
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    GK11(0) = ONE
                    DO L = 1, NMOMENTS_INPUT
                      DNL1 = DBLE(2*L + 1 )
                      F    = SSFDEL(N) * DNL1
                      FT   = ONE - SSFDEL(N)
                      GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N)-F)/FT
                    ENDDO
                  ELSE
                    DO L = 1, NMOMENTS_INPUT
                      GK11(L) = PHASMOMS_TOTAL_INPUT(L,N)
                   ENDDO
                  ENDIF
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
	            LEGPOLY = SS_PLEG_UP(V,1,L)
	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*LEGPOLY
	          ENDDO
	          EXACTSCAT_UP(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

C  End upwelling clause

            ENDIF

C  Downwelling Source terms: Phase matrix, multipliers, transmittances
C  -------------------------------------------------------------------

            IF ( DO_DNWELLING ) THEN

C  Multipliers and transmittances

              call outgoing_integration_dn
     i           ( nlayers, nfinelayers,
     i             do_fine, do_partials, extinction,
     i             n_offgrid_usertaus, offgrid_utau_layeridx,
     i             offgrid_utau_layerfineidx,
     i             sunpaths, radii, ntraverse, alpha_all,
     i             sunpaths_ut,   ntraverse_ut,   alpha_ut,
     i             sunpaths_fine, ntraverse_fine, alpha_fine,
     o             dn_multipliers(1,v), dn_lostrans(1,v),
     o             dn_multipliers_ut(1,v), dn_lostrans_ut(1,v) )

C  Debug
c              write(65,*)v
c              do n = 1, nlayers
c                write(65,'(1p2e18.10)')
c     &           dn_lostrans(n,v),dn_multipliers(n,v)
c              enddo

C  Debug: Multipliers (all layers), Down welling OLD WAY
c              call outgoing_integration_old_dn
c     i           ( nlayers, extinction, deltau_vert,
c     i             lospaths, sunpaths, radii, ntraverse, alpha_all,
c     o             dn_multipliers(1,v), dn_lostrans(1,v) )
c              write(66,*)v
c              do n = 1, nlayers
c                write(66,*)dn_lostrans(n,v),dn_multipliers(n,v)
c              enddo

C  Legendre polynomials

              COSSCAT = COSSCAT_DN(NLAYERS)
	      SS_PLEG_DN(V,1,0) = ONE
	      SS_PLEG_DN(V,1,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_DN(V,1,L) =
     &           DF1(L) * SS_PLEG_DN(V,1,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_DN(V,1,L-2)
              ENDDO

C  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(N) ) THEN
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    GK11(0) = ONE
                    DO L = 1, NMOMENTS_INPUT
                      DNL1 = DBLE(2*L + 1 )
                      F    = SSFDEL(N) * DNL1
                      FT   = ONE - SSFDEL(N)
                      GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N)-F)/FT
                    ENDDO
                  ELSE
                    DO L = 1, NMOMENTS_INPUT
                      GK11(L) = PHASMOMS_TOTAL_INPUT(L,N)
                   ENDDO
                  ENDIF
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
	            LEGPOLY = SS_PLEG_DN(V,1,L)
	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*LEGPOLY
	          ENDDO
	          EXACTSCAT_DN(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

C  End Downwelling clause

            ENDIF

C   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

C  Recurrence relation for the UPWELLING intensity
C  ===============================================

      IF ( DO_UPWELLING ) THEN

C  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C  Multiplier using new integration scheme

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_UP(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS(N,V)
              SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE +
     &               UP_LOSTRANS(N,V) * SS_CUMSOURCE_UP(V,NC-1)
            ENDDO
          ENDDO

C  Offgrid output-------
C    Add additional partial layer source term = Exact Phase Func * Multiplier
C    Get final cumulative source and set the Single scatter results
C  Ongrid output--------
C     Set final cumulative source and single scatter intensity

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP           = EXACTSCAT_UP(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
              TRANS          = UP_LOSTRANS_UT(UT,V)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  Recurrence relation for the DOWNWELLING intensity
C  =================================================

      IF ( DO_DNWELLING ) THEN

C  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO
        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      Multiplier by new integration method

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              HELP =  EXACTSCAT_DN(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS(N,V)
              SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE +
     &                DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,NC-1)
            ENDDO
          ENDDO

C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct the intensity
C  Ongrid output :
C     Set final cumulative source and correct Stokes vector

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_DN(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
              TRANS          = DN_LOSTRANS_UT(UT,V)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end optical depth loop and Downwelling clause

        ENDDO
      ENDIF

C  Finish

      RETURN
      END

c

      subroutine outgoing_sphergeom_fine
     i  ( maxlayers, maxfine, maxpartials,
     i    do_fine, do_partials, nlayers, nfine,
     i    n_partials, partials_idx,
     i    heights, heights_p, eradius,
     i    alpha_boa, theta_boa, phi_boa,
     o    sunpaths,      radii,      ntraverse,      alpha_all,
     o    sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,
     o    sunpaths_p,    radii_p,    ntraverse_p,    alpha_p,
     o    partials_fineidx, lospaths,
     o    lospaths_p_up, lospaths_p_dn, 
     o    theta_all, phi_all, cosscat_up, cosscat_dn,
     o    fail, message )

C  Completely stand-alone geometry routine for the outgoing correction
C    starting inputs are the BOA values of SZA, VZA and PHI
C    need also the height grids, earth radius and control

C  This routine has the fine gridding treatment
C  Version 2.3. September 2007, Partial layer geometries added

C  inputs

      integer          maxlayers, maxfine, maxpartials
      integer          nlayers, nfine
      logical          do_fine, do_partials
      integer          n_partials
      integer          partials_idx    (maxpartials)
      double precision eradius, heights (0:maxlayers)
      double precision heights_p(maxpartials)
      double precision alpha_boa, theta_boa, phi_boa

C  main outputs (geometry)

      integer          ntraverse   (0:maxlayers)
      double precision sunpaths   (0:maxlayers,maxlayers)
      double precision radii      (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level output (geometry)

      integer          ntraverse_fine(maxlayers,maxfine)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfine)
      double precision radii_fine    (maxlayers,maxfine)
      double precision alpha_fine    (maxlayers,maxfine)

C  Partial layer output (geometry)

      integer          ntraverse_p(maxpartials)
      double precision sunpaths_p (maxpartials,maxlayers)
      double precision radii_p    (maxpartials)
      double precision alpha_p    (maxpartials)

C  Other (incidental) geometrical output

      double precision lospaths(maxlayers)
      double precision lospaths_p_up(maxpartials)
      double precision lospaths_p_dn(maxpartials)
      integer          partials_fineidx(maxpartials)

      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

C  Status output

      logical          fail
      character*(*)    message

C  Local

      logical          direct_sun
      integer          n, k, krad, n1, ut, np
      double precision deg_to_rad, ex, ey, ez, px, py, pz
      double precision salpha_boa, calpha_boa, sphi_boa
      double precision stheta_boa, ctheta_boa, cphi_boa
      double precision ksi, cksi, sksi, xicum, tangr, fac
      double precision ctheta, stheta, calpha, salpha, cphi
      double precision b, sth0, th0, ks1, sth1, th1, theta_p

C  Local arrays associated with fine grid output

      integer          maxlocalfine, j
      parameter        ( maxlocalfine = 20 )
      logical          direct_sunf(maxlocalfine)
      double precision difz, dfine1, saf, xicum0, path
      double precision thetaf(maxlocalfine), xicumf, difa
      double precision cthetaf(maxlocalfine)
      double precision sthetaf(maxlocalfine)
      double precision ksif(maxlocalfine)

C  Initialise output

      fail = .false.
      message = ' '

C  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif 

C  zero the sun paths
C  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo

C  Zero the partial paths, initialize the traverse number (nominal)

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          ntraverse_p(ut) = np
          do k = 1, np
            sunpaths_p(ut,k) = 0.0d0
          enddo
        enddo
      endif

C  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = 0.0d0
         enddo
        enddo
       enddo
      endif

C  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

C  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))
      cosscat_up (nlayers) = - calpha_boa * ctheta_boa +
     &                         salpha_boa * stheta_boa * cphi_boa 
      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa +
     &                         salpha_boa * stheta_boa * cphi_boa 

C  Radii
C  -----

C  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

c  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif
 
C  Partial levels
C   Find the fine level just below partial point = (0,1,2,...nfine)
 
      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          radii_p(ut) = eradius + heights_p(ut)
          j = 1
          do while (radii_p(ut).gt.radii_fine(np,j).and.j.le.nfine)
           j = j + 1
          enddo
          partials_fineidx(ut) = j - 1
        enddo
      endif

C  Special case. Direct nadir viewing
C  ==================================

C  Compute everything and Exit.
C    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

C  WHOLE LAYER and FINE divisions
C  ------------------------------

C  Start layer loop, working upwards

        do n = nlayers,1,-1

C  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = 0.0d0
            enddo
          endif

C  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

C  Non-overhead sun
C  Main output of solar paths
C  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = dsin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = dasin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = dsin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

C  End main layer loop

        enddo

C  PARTIAL LAYERS
C  --------------

        if ( do_partials ) then
          do ut = 1, n_partials
            np = partials_idx(ut)

C  Los angle and paths

            alpha_p(ut) = alpha_all(nlayers)
            lospaths_p_dn(ut) = radii(np-1)-radii_p(ut)
            lospaths_p_up(ut) = radii_p(ut)-radii(np)

C  Overhead sun

            if (stheta_boa.eq.0.0d0 ) then
              do k = 1, np-1
                sunpaths_p(ut,k) = radii(k-1)-radii(k)
              enddo
              sunpaths_p(ut,np) = radii(np-1)-radii_p(ut)
            endif

C  Non-overhead sun

            if (stheta_boa.gt.0.0d0 ) then
              sth0 = stheta_boa
              th0  = theta_all(np)
              sth1 = sth0*radii_p(ut)/radii(np-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
              sth0 = sth1
              th0  = th1
              do k = np-1, 1, -1
                sth1 = sth0*radii(k)/radii(k-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
                sth0 = sth1
                th0  = th1
              enddo
            endif

C  Finish partial layer stuff

          enddo
        endif

C  Return, as everything now done

        return

C  end regular pseudo-spherical clause, LOS is zero

      endif

C  Outgoing spehricity geometry
C  ============================

C  define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

C  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

C  Check single illumination
C      --Not required, now we have the tangent point treatment
c      if (stheta_boa.gt.0.0d0 ) then
c        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
c        xicum = alpha_all(nlayers)-xicum
c        px = - radii(0) * dsin(xicum)
c        py = 0.0d0
c        pz =   radii(0) * dcos(xicum)
c        b = ex*px + ey*py + ez*pz
c        ctheta = -b/radii(0)
c        if ( ctheta.le.0.0d0 ) then
c          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
c        endif
c      endif

C  initialise los cumulative angle

      xicum  = 0.0d0

C  set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

C  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

C  Next level up

        n1 = n + 1
  
C  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

C  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

C  Fine grid lospath output (angle and radius)
C    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = dsin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

C  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = dcos(thetaf(j))
             sthetaf(j) = dsqrt(1.0d0-ctheta*ctheta)
           enddo
         endif
        endif

C  Sun angles for the general case
C    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * dsin(xicumf)
             py = 0.0d0
             pz =   radii_fine(n1,j) * dcos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = dsqrt(1.0d0-cthetaf(j)*cthetaf(j))
             thetaf(j)  = dacos(cthetaf(j))
           enddo
         endif
        endif

C  Unit vector f2(i) perpendicular to OP but in plane of path
C  projection of f2(i) on solar path gives the relative azimuth at P
c        f2x = dsin(xicum)
c        f2y = 0.0d0
c        f2z = dcos(xicum)
c        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
c        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
c        if ( cphi.gt.1.0d0)  cphi = 1.0d0
c        if ( cphi.lt.-1.0d0) cphi = -1.0d0
C ********************************************* Apparently not correct

C  Fix phi by using constancy of scatter angle
C  Only for the scattering up directions..................

        cosscat_up(n) = cosscat_up(n+1)
        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

C  Sun paths, Direct sun at layer top
C  ==================================

C   Means that the SZA at layer top is < 90.
C    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

C  Work up from level n to TOA
C    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

C DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

C  Fine grid calculation is always required
C  ----------------------------------------

C   Start at the grid point on LOS path and work upwards,
C     first across partial layer to upper boundary, then continue
C     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
C  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,-1)
          enddo

C  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

C  Complete direct sun computations

        endif

C  Sun paths, Not direct sun , with tangent point
C  ==============================================

C  Although layer top has a tangent point, not all of the fine-grid
C   points will have a tangent point.

        if (.not.direct_sun ) then

C  First do the layer-top calculation.
C  -----------------------------------

C  TANGR = tangent point radius.

         tangr = stheta*radii(n)

C  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

C  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

C  Work downwards from TOA (sine rule) to level immediately above
C  the tangent layer. Don't forget to double the path length for
C  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

C  Tangent layer path length. Twice again. The following check is good.
c  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)   

C DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1,-1)

C  Fine layer development (slightly different)
C  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

C  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = dasin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

C DBG           write(*,*)'regular',n1,j,
C     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

C  Fine grid Calculation with tangent point
C  ----------------------------------------

           else

C  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

C  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = dasin(sth0)

C  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = dasin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

C  In layer below level n, Work down to level of fine grid radius
C     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = dasin(sth1)
            ks1 = th1-th0
            path = dsin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

C  In layer below level n, Complete the path down to the tangent point
C    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + 2.0d0*radii_fine(n1,j)*dcos(th1)
              sunpaths_fine(n1,krad+1,j) = path

C  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*dcos(th1)

C  In layer below level n, Need to go down further.
C    (from now on, we need double contributions)
C     --- Continue the path to the bottom of this layer
C     --- Continue down by whole-layer steps until reach level above Tangent
C     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = dasin(sth1)
              ks1 = th1-th0
              path = path + 2.0d0 * dsin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1              
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = dasin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=2.0d0*radii(krad)*dcos(th1)
C  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*dcos(th1)
            endif

c DBG           write(*,*)'tangent',n1,j,
c     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

C  Complete tangent point clause

           endif

C  Complete fine-grid loop

          enddo

c  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

C  Complete tangent point calculation

         endif

        endif

C  End layer loop

      enddo

C  PARTIALS Section
C  ----------------

C  Partial angles and lospaths

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)

C  establish the LOS angle, and lospaths

          th0 = alpha_all(np)
          sth0 = dsin(th0)
          salpha = radii(np)*sth0 / radii_p(ut)
          alpha_p(ut) = dasin(salpha)
          calpha = dsqrt(1.0d0-salpha*salpha)

          ksi = alpha_all(np) - alpha_p(ut)
          lospaths_p_up(ut) = dsin(ksi)*radii(np)/salpha
          ksi = alpha_p(ut) - alpha_all(np-1)
          lospaths_p_dn(ut) = dsin(ksi)*radii(np-1)/salpha

C  Check debug
c          lospaths_p_dn(ut) = lospaths(np)-lospaths_p_up(ut)
c          write(*,*)lospaths_p_up(ut),lospaths_p_dn(ut),
c     &              radii_p(ut)-radii(np)
c          pause

C  Establish the solar angle, both naidr =sun-at-BOA and generally

          xicum = alpha_all(nlayers) - alpha_p(ut)
          direct_sun = .true.
          if (stheta_boa.eq.0.0d0 ) then
            theta_p = xicum
            ctheta = dcos(theta_p)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
          else
            px = - radii_p(ut) * dsin(xicum)
            py = 0.0d0
            pz =   radii_p(ut) * dcos(xicum)
            b = ex*px + ey*py + ez*pz
            ctheta = -b/radii_p(ut)
            direct_sun = (direct_sun.and.ctheta.ge.0.d0)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            theta_p = dacos(ctheta)
          endif

c         write(*,*)direct_sun,theta_all(Np-1),theta_p,theta_all(np)

C  Sun paths, Direct sun at layer top
C   First calculate sunpath for layer containing partial point
C   Then work upwards to TOA using the sine rule.

          if ( direct_sun ) then
            sth0 = stheta
            th0  = theta_p
            sth1 = sth0*radii_p(ut)/radii(np-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
            sth0 = sth1
            th0  = th1
            do k = np-1, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
          endif

C  Debug
c 14       format(a11,1x,10f10.6)
c          j = partials_fineidx(ut)
c          write(*,*)ut,ntraverse_p(ut)
c          write(*,*)ut,j,np,lospaths_p_up(ut)/lospaths(np)
c           write(*,14)'regular feb',(sunpaths(np,k),k=np,1,-1)
c           do j = 1, nfine
c             write(*,14)'regular feb',(sunpaths_fine(np,k,j),k=np,1,-1)
c             if (j.eq.partials_fineidx(ut))
c     &        write(*,14)'regular ut ',(sunpaths_p(ut,k),k=np,1,-1)
c           enddo
c           write(*,14)'regular feb',0.0d0,(sunpaths(np-1,k),k=np-1,1,-1)
c          
C  Sun paths, Not direct sun , with tangent point. Code similar to above.
C  TANGR = tangent point radius for this ray
C   ntraverse_p(ut) = Number of layers traversed by ray.

C  Work downwards from TOA (sine rule) to level immediately above
C  the tangent layer. Don't forget to double the path length for
C  any layers which are traversed twice.

C  Tangent layer path length. Twice again. The following check is good.
c  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

          if (.not.direct_sun ) then
            tangr = stheta*radii(n)
            krad = nlayers
            do while (tangr.gt.radii(krad))
              krad = krad - 1
            enddo
            ntraverse_p(ut) = krad + 1
            sth0 = tangr/radii(0)
            th0 = dasin(sth0)
            do k = 1, krad
              sth1 = radii(k-1)*sth0/radii(k)
              th1 = dasin(sth1)
              ks1 = th1-th0
              fac = 1.0d0
              if ( k.gt.np) fac = 2.0d0
              sunpaths_p(ut,k) = fac*dsin(ks1)*radii(k-1)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            sunpaths_p(ut,krad+1)=2.0d0*radii_p(ut)*dcos(th1)   
          endif

C  Finish partials loop

        enddo
      endif

C  Finish

      return
      end

c

      subroutine outgoing_integration_up
     i  ( nlayers, nfinelayers, do_fine, do_partials, extinction,
     i    n_partials, partials_idx, partials_fineidx,
     i    sunpaths, radii, ntraverse, alpha_all,
     i    sunpaths_p,    ntraverse_p,    alpha_p,
     i    sunpaths_fine, ntraverse_fine, alpha_fine,
     o    multipliers, lostrans, boa_attn,
     o    multipliers_p, lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/LIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

C  control

      logical          do_fine, do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_offgrid_usertaus)
      integer          partials_fineidx(max_offgrid_usertaus)

C  Whole layers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layers

      integer          ntraverse_p(max_offgrid_usertaus)
      double precision sunpaths_p (max_offgrid_usertaus,maxlayers)
      double precision alpha_p    (max_offgrid_usertaus)

C  Extinction 

      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)
      double precision multipliers_p (max_offgrid_usertaus)
      double precision lostrans_p    (max_offgrid_usertaus)
      double precision boa_attn

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn      ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_offgrid_usertaus )

C  help variables
C  --------------

      integer          n, j, k, ut, np, nfine_p
      double precision tau, sum, salpha, calpha, dfine1, raycon, kn
      double precision csq_1, csq_2, cot_1, cot_2
      double precision func_1, func_2, tran_1, tran, step, func
      double precision term1, term2, dfine_p, step_f, step_p
      
C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in VLIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

C  Whole layers

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Partial layers
 
      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = 0.0d0
          lostrans_p(ut)    = 0.0d0
        enddo
      endif

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

C  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = 0.0d0
       attn(n) = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo
      boa_attn = attn(nlayers)
    
C  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = 0.0d0
          attn_p(ut) = 0.0d0
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif

C  Fine grid attenuations

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = 0.0d0
        attn_fine(n,j) = 0.0d0
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

C  Work up from the bottom of the atmosphere
C  =========================================

C  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

C  Start layer loop

      do n = nlayers, 1, -1

C  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration

C   Surely this is wrong.................???? 27 September 2007
c        func_1 = attn(n-1)   * csq_1 
c        tran_2 = dexp ( - kn * ( cot_2 - cot_1 ) )
c        func_2 = attn(n) * csq_2 * tran_2

        func_2 = attn(n-1)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
        func_1 = attn(n) * csq_1 * tran_1
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          tran = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
          func = attn_fine(n,j) * tran * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

c        if (n.eq.6)write(19,*)'whole',n,multipliers(n),lostrans(n)

C  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

C  Finish layer

      ENDDO

C  Finish if no partials

      if ( .not. do_partials ) return
      
C  Partial layer functions
C  =======================

C  start partials loop

      do ut = 1, n_partials
      
C  layer of occurrence and fine-layer cutoff

        np     = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .gt. 0 ) then
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_fine(np,nfine_p)-alpha_p(ut))
        else
          step_p = (alpha_all(np)-alpha_p(ut))
        endif
        
C  Initialize for layer of occurrence

        salpha = dsin(alpha_all(np))
        calpha = dcos(alpha_all(np))
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

C  Top of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

C  saved fine-grid quantities as far as cutoff

        do j = 1, nfine_p
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration
C  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
        func_1 = attn(np) * csq_1 * tran_1

        if ( nfine_p.gt.0) then
          sum = 0.5d0 * func_1
          do j = 1, nfine_p
            tran = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
            func = attn_fine(np,j) * tran * csq_fine(j)
            if ( j.lt.nfine_p ) sum = sum + func
            if ( j.eq.nfine_p ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

c        write(19,*)'partial up',ut,np,nfine_p,
c     &              multipliers_p(ut),lostrans_p(ut)

C  Finish partials loop

      ENDDO

C  finish

      RETURN
      END

C

      subroutine outgoing_integration_dn
     i  ( nlayers, nfinelayers, do_fine, do_partials, extinction,
     i    n_partials, partials_idx, partials_fineidx,
     i    sunpaths, radii, ntraverse, alpha_all,
     i    sunpaths_p,   ntraverse_p,   alpha_p,
     i    sunpaths_fine, ntraverse_fine, alpha_fine,
     o    multipliers, lostrans,
     o    multipliers_p, lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/LIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

C  control

      logical          do_fine, do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_offgrid_usertaus)
      integer          partials_fineidx(max_offgrid_usertaus)

C  Whole layers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layers

      integer          ntraverse_p(max_offgrid_usertaus)
      double precision sunpaths_p (max_offgrid_usertaus,maxlayers)
      double precision alpha_p    (max_offgrid_usertaus)

C  Extinction 

      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)
      double precision multipliers_p (max_offgrid_usertaus)
      double precision lostrans_p    (max_offgrid_usertaus)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn      ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_offgrid_usertaus )

C  help variables
C  --------------

      integer          n, j, k, ut, np, nfine_p
      double precision tau, sum, salpha, calpha, dfine1, raycon, kn
      double precision csq_1, csq_2, cot_1, cot_2
      double precision func_1, func_2, tran_1, tran, step, func
      double precision term1, term2, dfine_p, step_f, step_p
      
C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in VLIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

C  Whole layers

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Partial layers

      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = 0.0d0
          lostrans_p(ut)    = 0.0d0
        enddo
      endif

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

C  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = 0.0d0
       attn(n) = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo

C  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = 0.0d0
          attn_p(ut) = 0.0d0
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif

C  Attenuations for fine layer stuff

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = 0.0d0
        attn_fine(n,j) = 0.0d0
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

C  Work Down from the top of the atmosphere
C  ========================================

C  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

C  Start layer loop

      do n = 1, nlayers

C  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1
        do j = nfinelayers, 1, -1
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration

        tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
        func_1 = attn(n-1) * csq_1  * tran_1
        func_2 = attn(n)   * csq_2
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          tran = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
          func = attn_fine(n,j) * tran * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

c        if (n.eq.6)write(20,*)'whole',n,multipliers(n),lostrans(n)

C  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

C  Finish layer

      ENDDO

C  Finish if no partials

      if ( .not. do_partials ) return
      
C  Partial layer functions
C  =======================

C  start partials loop

      do ut = 1, n_partials
      
C  layer of occurrence and fine-layer cutoff

        np     = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .eq. nfinelayers ) then
          step_p = (alpha_p(ut)-alpha_all(np-1))
        else if ( nfine_p.eq.0 ) then
          step_f = (alpha_all(np) - alpha_all(np-1))/dfine1
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        else
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        endif

C  Top of layer of occurrence

        salpha = dsin(alpha_all(np-1))
        calpha = dcos(alpha_all(np-1))
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

C  End of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

C  saved fine-grid quantities as far as cutoff

        do j = nfinelayers, nfine_p+1, -1
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration
C  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
        func_1 = attn(np-1) * csq_1 * tran_1
        if ( nfine_p.lt.nfinelayers) then
          sum = 0.5d0 * func_1
          do j = nfinelayers, nfine_p+1, -1
            tran = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
            func = attn_fine(np,j) * tran * csq_fine(j)
            if ( j.gt.nfine_p+1 ) sum = sum + func
            if ( j.eq.nfine_p+1 ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

c        write(20,*)'partial dn',ut,np,nfine_p,
c     &           multipliers_p(ut),lostrans_p(ut)

C  Finish partials loop

      ENDDO

C  finish

      RETURN
      END

c
      
      subroutine multi_outgoing_adjustgeom
     i   ( max_vza, max_sza, max_azm, n_vza, n_sza, n_azm,
     i     hsurface, eradius, do_adjust_surface,
     i     alpha_boa, theta_boa, phi_boa,
     o     alpha_ssa, theta_ssa, phi_ssa, fail, message )

C stand-alone geometry routine for adjusting the outgoing correction
C    starting inputs are the BOA values of SZA, VZA and PHI
C    need also the height of new surface, earth radius.

C  Height grid here is artificial

C  inputs

      integer          max_vza, max_sza, max_azm
      integer          n_vza, n_sza, n_azm
      double precision eradius, hsurface
      logical          do_adjust_surface
      double precision alpha_boa (max_vza)
      double precision theta_boa (max_sza)
      double precision phi_boa   (max_azm)

C  outputs

      double precision alpha_ssa (max_vza)
      double precision theta_ssa (max_vza,max_sza,max_azm)
      double precision phi_ssa   (max_vza,max_sza,max_azm)
      logical          fail
      character*(*)    message

C  Local

      integer          j, i, k
      double precision deg_to_rad, ex, ey, ez, px, py, pz
      double precision salpha_boa, calpha_boa, sphi_boa
      double precision stheta_boa, ctheta_boa, cphi_boa
      double precision ksi, cksi, sksi, xicum, cosscat_up
      double precision phi_all, alpha_all, theta_all
      double precision ctheta, stheta, calpha, salpha, cphi
      double precision b,rssa

C  Initialise output

      fail = .false.
      message = ' '

C  check range of inputs

      do j = 1, n_vza
       if ( alpha_boa(j).ge.90.0d0.or.alpha_boa(j).lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
       endif
      enddo

      do k = 1, n_azm
        if ( phi_boa(k).lt.0.0d0   ) phi_boa(k) = - phi_boa(k)
        if ( phi_boa(k).gt.180.0d0 ) phi_boa(k) = 360.0d0 - phi_boa(k)
      enddo

      do i = 1, n_sza
       if ( theta_boa(i).ge.90.0d0.or.theta_boa(i).lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
       endif
      enddo

C  No adjustment, just copy and exit

      if ( .not. do_adjust_surface ) then
        do j = 1, n_vza
          alpha_ssa(j)   = alpha_boa(j)
          do i = 1, n_sza
            do k = 1, n_azm
              theta_ssa(j,i,k) = theta_boa(i)
              phi_ssa(j,i,k)   = phi_boa(k)
            enddo
          enddo
        enddo
        return
      endif

C  conversion

      deg_to_rad = dacos(-1.0d0) / 180.0d0

C  Radius of surface

      rssa = hsurface + eradius

C  Start VZA loop

      do j = 1, n_vza

        alpha_all = alpha_boa(j) * deg_to_rad
        salpha_boa = dsin(alpha_all)
        calpha_boa = dcos(alpha_all)

C  Special case. Direct nadir viewing. Compute everything and Exit.
C    (This is the same as the regular pseudo-spherical )

        if ( salpha_boa.eq.0.0d0 ) then
          alpha_ssa(j)   = alpha_boa(j)
          do i = 1, n_sza
            do k = 1, n_azm
              theta_ssa(j,i,k) = theta_boa(i)
              phi_ssa(j,i,k)   = phi_boa(k)
            enddo
          enddo
          go to 567
        endif

C  Los angle

        salpha       = eradius * salpha_boa / rssa
        alpha_ssa(j) = dasin(salpha)
        calpha       = dcos(alpha_ssa(j))

C  Lospaths

        ksi  = alpha_all - alpha_ssa(j)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        xicum = ksi

C  output angle in degrees

        alpha_ssa(j) = alpha_ssa(j) / deg_to_rad

C  Los vector

        px = - rssa * dsin(xicum)
        py = 0.0d0
        pz =   rssa * dcos(xicum)

C  Start SZA loop

        do i = 1, n_sza

          theta_all = theta_boa(i) * deg_to_rad
          stheta_boa = dsin(theta_all)
          ctheta_boa = dcos(theta_all)

C  Start azimuth loop

          do k = 1, n_azm

            phi_all   = phi_boa(k)   * deg_to_rad
            cphi_boa  = dcos(phi_all)
            sphi_boa  = dsin(phi_all)

C  define Unit solar vector

            ex = - stheta_boa * cphi_boa
            ey = - stheta_boa * sphi_boa
            ez = - ctheta_boa
  
C  Sun angle

            b = ex*px + ey*py + ez*pz
            ctheta = -b/rssa
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            theta_ssa(j,i,k) = dacos(ctheta)/deg_to_rad
            if ( ctheta.lt.0.0d0 ) stop'failed'

C  scattering angle

            cosscat_up  = - calpha_boa * ctheta_boa +
     &                      salpha_boa * stheta_boa * cphi_boa 

C  Fix phi by using constancy of scatter angle

            if ( phi_boa(k).eq.180.0d0 ) then
              phi_ssa(j,i,k) = phi_all
            else if ( phi_boa(k) .eq. 0.0d0 ) then
              phi_ssa(j,i,k) = 0.0d0
            else
              cphi = (cosscat_up+calpha*ctheta)/stheta/salpha
              if ( cphi.gt.1.0d0) cphi = 1.0d0
              if ( cphi.lt.-1.0d0) cphi = -1.0d0
              phi_ssa(j,i,k) = dacos(cphi)
            endif
            phi_ssa(j,i,k) = phi_ssa(j,i,k) / deg_to_rad

C  End Azimuth and solar loops

          enddo
        enddo

C  Continuation point

 567    continue

C  End los loop

      enddo

C  Finish

      return
      end

C

      SUBROUTINE LIDORT_LAMBERTIAN_DBCORRECTION (FLUXMULT)

C  Prepares Exact Direct Beam reflection for the Lambertian case

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Input arguments controlling type of surface

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Result goes in the Single scatter output file

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V
      DOUBLE PRECISION FINAL_SOURCE, TR, FACTOR
      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN

C  first stage
C  -----------

C  Initialize

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_OUT_USERTAUS
          INTENSITY_DB(UTA,V)  = ZERO
        ENDDO
      ENDDO

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Reflection of solar beam
C  ------------------------

C  Must use adjusted values here.

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  Start geometry loops

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA

c  Beam attenuation

              IF ( DO_SSCORR_OUTGOING ) THEN
                X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                ATTN = BOA_ATTN(V)               
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                ELSE
                  X0_BOA = X0(IB)
                ENDIF 
                ATTN = SOLAR_BEAM_OPDEP(IB)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              ATTN    = ATTN * X0_FLUX
              ATTN_DB_SAVE(V) = ATTN

C  Finish loops over geometries

            ENDDO
          ENDIF
        ENDDO
      ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term = F.A.mu_0.T/pi
C    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
C    Only require the Stokes total intensity component

      NC =  0
      FACTOR = FLUXMULT * LAMBERTIAN_ALBEDO
      DO V = 1, N_GEOMETRIES
        DB_CUMSOURCE(V,NC) = ATTN_DB_SAVE(V) * FACTOR
      ENDDO

C  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT.
C     I-Stokes component only

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS(N,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR = T_DELT_USERM(N,UM)
                ENDIF
                DB_CUMSOURCE(V,NC) = TR * DB_CUMSOURCE(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

C  end layer loop

        ENDDO

C  Offgrid output : partial layer transmittance, then set result

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS_UT(UT,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR = T_UTUP_USERM(UT,UM)
                ENDIF
                INTENSITY_DB(UTA,V) = TR * DB_CUMSOURCE(V,NC)
              ENDDO
            ENDDO
          ENDDO

C  Ongrid output : Set final cumulative source directly

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                FINAL_SOURCE = DB_CUMSOURCE(V,NC)
                INTENSITY_DB(UTA,V) = FINAL_SOURCE
              ENDDO
            ENDDO
          ENDDO

        ENDIF

C  Check for updating the recursion 

        IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

C  end optical depth loop

      ENDDO

C  Finish

      RETURN
      END

C
C  REDUNDANT ROUTINES
C

      subroutine outgoing_sphergeom
     &  ( maxlayers, nlayers, heights, eradius,
     &    alpha_boa, theta_boa, phi_boa,
     o    lospaths, sunpaths, radii, ntraverse,
     o    alpha_all, theta_all, phi_all, cosscat_up, cosscat_dn,
     o    fail, message )

C  Completely stand-alone geometry routine for the outgoing correction
C    starting inputs are the BOA values of SZA, VZA and PHI
C    need also the height grid, earth radius and control

C  inputs

      integer maxlayers, nlayers
      double precision eradius, heights (0:maxlayers)
      double precision alpha_boa, theta_boa, phi_boa

C  outputs

      integer          ntraverse(0:maxlayers)
      double precision lospaths(maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)

      double precision alpha_all  (0:maxlayers)
      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

      logical          fail
      character*(*)    message

C  Local

      logical          direct_sun
      integer          n, k, krad
      double precision deg_to_rad, ex, ey, ez, px, py, pz
      double precision salpha_boa, calpha_boa, sphi_boa
      double precision stheta_boa, ctheta_boa, cphi_boa
      double precision ksi, cksi, sksi, xicum, tangr, fac
      double precision ctheta, stheta, calpha, salpha, cphi
      double precision c, b, rtoasq, sth0, th0, ks1, sth1, th1

C  Initialise output

      fail = .false.
      message = ' '

C  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif

C  zero the sun paths
C  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
       ntraverse(n) = n
       do k = 1, nlayers
        sunpaths(n,k) = 0.0d0
       enddo
      enddo

C  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

C  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))
      cosscat_up (nlayers) = - calpha_boa * ctheta_boa +
     &                         salpha_boa * stheta_boa * cphi_boa 
      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa +
     &                         salpha_boa * stheta_boa * cphi_boa 

C  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo
      rtoasq =  radii(0)*radii(0)

C  Special case. Direct nadir viewing. Compute everything and Exit.
C    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then
        do n = nlayers,1,-1
          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up (n-1) = cosscat_up(n)
          cosscat_dn (n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)
          if (stheta_boa.eq.0.0d0 ) then
           do k = n, 1, -1
             sunpaths(n,k) = radii(k-1)-radii(k)
           enddo
          else
           sth0 = stheta_boa
           th0  = theta_all(n)
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
          endif
        enddo
        return
      endif

C  define Unit solar vector

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

C  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

C  Check single illumination
C      --Not required, now we have the tangent point treatment
c      if (stheta_boa.gt.0.0d0 ) then
c        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
c        xicum = alpha_all(nlayers)-xicum
c        px = - radii(0) * dsin(xicum)
c        py = 0.0d0
c        pz =   radii(0) * dcos(xicum)
c        b = ex*px + ey*py + ez*pz
c        ctheta = -b/radii(0)
c        if ( ctheta.le.0.0d0 ) then
c          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
c        endif
c      endif

C  initialise los cumulative angle

      xicum  = 0.0d0

C  set TOA direct illumination flag

      direct_sun = .true.

C  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1
  
C  Los angles

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

C  Lospaths

        ksi = alpha_all(n+1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n+1) = sksi * radii(n+1) / salpha
        xicum = xicum + ksi

C  Sun angle

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(xicum)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
        else
         c = rtoasq - radii(n)*radii(n)
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
        endif

C  Unit vector f2(i) perpendicular to OP but in plane of path
C  projection of f2(i) on solar path gives the relative azimuth at P
c        f2x = dsin(xicum)
c        f2y = 0.0d0
c        f2z = dcos(xicum)
c        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
c        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
c        if ( cphi.gt.1.0d0)  cphi = 1.0d0
c        if ( cphi.lt.-1.0d0) cphi = -1.0d0
C ********************************************* Apparently not correct

C  Fix phi by using constancy of scatter angle

        cosscat_up(n) = cosscat_up(n+1)
        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
         phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
c         cphi = (cosscat_dn(n)-calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

C  Sun paths, full illumination at TOA ---> stop

        if ( n.eq.0 .and. direct_sun )  go to 345

C  Sun paths, full illumination at other (non-boa, non_toa) levels

        if ( direct_sun ) then
          if ( stheta_boa.eq.0.0d0.and.n.gt.0 ) then
            do k = n, 1, -1
             sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
          else
            sth0 = stheta
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
          endif
          go to 345
        endif

C  Sun paths, non-boa level, with tangent point

        if (.not.direct_sun ) then
          tangr = stheta*radii(n)
          krad = nlayers
          do while (tangr.gt.radii(krad))
            krad = krad - 1
          enddo
          ntraverse(n) = krad + 1
          sth0 = tangr/radii(0)
          th0 = dasin(sth0)
c          write(*,'(a,2i3,2f12.5)')'tp',n,krad,tangr,th0/deg_to_rad
          do k = 1, krad
            sth1 = radii(k-1)*sth0/radii(k)
            th1 = dasin(sth1)
            ks1 = th1-th0
            fac = 1.0d0
            if ( k.gt.n) fac = 2.0d0
            sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
            sth0 = sth1
            th0  = th1
          enddo
c          write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)
          sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)   
        endif

c  Continuation point

 345    continue

C  End

      enddo

C  Finish

      return
      end


C

      subroutine outgoing_integration_old_up
     i  ( nlayers, extinction, deltaus,
     i    lospaths, sunpaths, radii, ntraverse, alpha_all,
     o    multipliers, lostrans )

C  Does the optical depth integration over layers.

C  include files

      include '../includes/LIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

      integer          nlayers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision lospaths(maxlayers)
      double precision radii   (0:maxlayers)

      double precision alpha_all  (0:maxlayers)

C  Extinction 

      double precision deltaus    (maxlayers)
      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)

C  local variables
C  ---------------

      INTEGER          N, K, G2_NILLUM
      DOUBLE PRECISION SP, S_T_0, S_T_1
      DOUBLE PRECISION FAC0, FAC1, TAU

C  local arrays

      DOUBLE PRECISION TSUN ( MAXLAYERS )
      DOUBLE PRECISION TMU  ( MAXLAYERS )
      DOUBLE PRECISION SUN_ATRAN ( MAXLAYERS )
      DOUBLE PRECISION SUN_AVSEC ( MAXLAYERS )
      DOUBLE PRECISION LOS_AVSEC ( MAXLAYERS )
      DOUBLE PRECISION TAU_SPHER ( 0:MAXLAYERS )

C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in VLIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  attenuation functions

      do n = 0, nlayers
       tau     = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       tau_spher(n) = tau
      enddo

C  average secants and initial transmittance factors, solar beam
C  -------------------------------------------------------------

      S_T_0 = 1.0d0
      G2_NILLUM = NLAYERS
      DO N = 1, NLAYERS
        LOS_AVSEC(N) = LOSPATHS(N) * EXTINCTION(N) / DELTAUS(N)
        IF  (N.LE.G2_NILLUM ) THEN
          IF ( TAU_SPHER(N) .GT. LOCAL_CUTOFF ) THEN
            G2_NILLUM = N
          ELSE
            S_T_1 = DEXP ( - TAU_SPHER(N) )
          ENDIF
          SUN_AVSEC(N) = (TAU_SPHER(N)-TAU_SPHER(N-1))/DELTAUS(N)
          SUN_ATRAN(N) = S_T_0
          S_T_0    = S_T_1
        ELSE
          SUN_AVSEC(N) = ZERO
          SUN_ATRAN(N) = ZERO
        ENDIF
      ENDDO

C  path transmittances TSUN and TMU (los)
C  --------------------------------------

       DO N = 1, NLAYERS
        IF ( N. GT. G2_NILLUM ) THEN
          TSUN(N) = ZERO
        ELSE
          TAU = DELTAUS(N) * SUN_AVSEC(N)
          IF ( TAU .GT. LOCAL_CUTOFF ) THEN
            TSUN(N) = ZERO
          ELSE
            TSUN(N) = DEXP ( - TAU )
          ENDIF
        ENDIF
        TMU(N) = DEXP ( - DELTAUS(N) * LOS_AVSEC(N) )
      ENDDO

C  Multiplier - Work up from the bottom of the atmosphere
C  ------------------------------------------------------

      DO N = NLAYERS, 1, -1
        SP   = SUN_AVSEC(N) + LOS_AVSEC(N)
        FAC0 = SUN_ATRAN(N) * LOS_AVSEC(N) / SP
        FAC1 = 1.0d0 - TMU(N) * TSUN(N)
        MULTIPLIERS(N) = FAC0 * FAC1
        LOSTRANS(N)   = TMU(N)
      ENDDO

C  finish

      RETURN
      END
      
C

      subroutine outgoing_integration_old_dn
     i  ( nlayers, extinction, deltaus,
     i    lospaths, sunpaths, radii, ntraverse, alpha_all,
     o    multipliers, lostrans )

C  Does the optical depth integration over layers.

C  include files

      include '../includes/LIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

      integer          nlayers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision lospaths(maxlayers)
      double precision radii   (0:maxlayers)

      double precision alpha_all  (0:maxlayers)

C  Extinction 

      double precision deltaus    (maxlayers)
      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)

C  local variables
C  ---------------

      INTEGER          N, K, G2_NILLUM
      DOUBLE PRECISION SP, S_T_0, S_T_1
      DOUBLE PRECISION FAC0, FAC1, TAU

C  local arrays

      DOUBLE PRECISION TSUN ( MAXLAYERS )
      DOUBLE PRECISION TMU  ( MAXLAYERS )
      DOUBLE PRECISION SUN_ATRAN ( MAXLAYERS )
      DOUBLE PRECISION SUN_AVSEC ( MAXLAYERS )
      DOUBLE PRECISION LOS_AVSEC ( MAXLAYERS )
      DOUBLE PRECISION TAU_SPHER ( 0:MAXLAYERS )

C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in VLIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  attenuation functions

      do n = 0, nlayers
       tau     = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       tau_spher(n) = tau
      enddo

C  average secants and initial transmittance factors, solar beam
C  -------------------------------------------------------------

      S_T_0 = 1.0d0
      G2_NILLUM = NLAYERS
      DO N = 1, NLAYERS
        LOS_AVSEC(N) = LOSPATHS(N) * EXTINCTION(N) / DELTAUS(N)
        IF  (N.LE.G2_NILLUM ) THEN
          IF ( TAU_SPHER(N) .GT. LOCAL_CUTOFF ) THEN
            G2_NILLUM = N
          ELSE
            S_T_1 = DEXP ( - TAU_SPHER(N) )
          ENDIF
          SUN_AVSEC(N) = (TAU_SPHER(N)-TAU_SPHER(N-1))/DELTAUS(N)
          SUN_ATRAN(N) = S_T_0
          S_T_0    = S_T_1
        ELSE
          SUN_AVSEC(N) = ZERO
          SUN_ATRAN(N) = ZERO
        ENDIF
      ENDDO

C  path transmittances TSUN and TMU (los)
C  --------------------------------------

       DO N = 1, NLAYERS
        IF ( N. GT. G2_NILLUM ) THEN
          TSUN(N) = ZERO
        ELSE
          TAU = DELTAUS(N) * SUN_AVSEC(N)
          IF ( TAU .GT. LOCAL_CUTOFF ) THEN
            TSUN(N) = ZERO
          ELSE
            TSUN(N) = DEXP ( - TAU )
          ENDIF
        ENDIF
        TMU(N) = DEXP ( - DELTAUS(N) * LOS_AVSEC(N) )
      ENDDO

C  Multiplier - Work down from the top of the atmosphere
C  -----------------------------------------------------

      DO N = 1, NLAYERS
        SP   = SUN_AVSEC(N) - LOS_AVSEC(N)
        FAC0 = SUN_ATRAN(N) * LOS_AVSEC(N) / SP
        FAC1 = TMU(N) - TSUN(N)
        MULTIPLIERS(N) = FAC0 * FAC1
        LOSTRANS(N)   = TMU(N)
      ENDDO

C  finish

      RETURN
      END

c

      subroutine outgoing_adjustgeom
     &  ( hsurface, eradius,
     &    alpha_boa, theta_boa, phi_boa,
     o    alpha_ssa, theta_ssa, phi_ssa, fail, message )

C stand-alone geometry routine for adjusting the outgoing correction
C    starting inputs are the BOA values of SZA, VZA and PHI
C    need also the height of new surface, earth radius.

C  Height grid here is artificial

C  inputs

      double precision eradius, hsurface
      double precision alpha_boa, theta_boa, phi_boa

C  outputs

      double precision alpha_ssa, theta_ssa, phi_ssa
      logical          fail
      character*(*)    message

C  Local

      double precision deg_to_rad, ex, ey, ez, px, py, pz
      double precision salpha_boa, calpha_boa, sphi_boa
      double precision stheta_boa, ctheta_boa, cphi_boa
      double precision ksi, cksi, sksi, xicum, cosscat_up
      double precision phi_all, alpha_all, theta_all
      double precision ctheta, stheta, calpha, salpha, cphi
      double precision b,rssa

C  Initialise output

      fail = .false.
      message = ' '

C  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif

C  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all = alpha_boa * deg_to_rad
      theta_all = theta_boa * deg_to_rad
      phi_all   = phi_boa   * deg_to_rad

C  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all)
      calpha_boa = dcos(alpha_all)
      stheta_boa = dsin(theta_all)
      ctheta_boa = dcos(theta_all)
      cphi_boa   = dcos(phi_all)
      sphi_boa   = dsin(phi_all)
      cosscat_up  = - calpha_boa * ctheta_boa +
     &                salpha_boa * stheta_boa * cphi_boa 

C  Radii

      rssa = hsurface + eradius

C  Special case. Direct nadir viewing. Compute everything and Exit.
C    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then
        alpha_ssa   = alpha_boa
        theta_ssa   = theta_boa
        phi_ssa     = phi_boa
        return
      endif

C  define Unit solar vector

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa
  
C  Los angle

      salpha     = eradius * salpha_boa / rssa
      alpha_ssa  = dasin(salpha)
      calpha     = dcos(alpha_ssa)

C  Lospaths

      ksi  = alpha_all - alpha_ssa
      sksi = dsin(ksi)
      cksi = dcos(ksi)
      xicum = ksi

C  Sun angle

      px = - rssa * dsin(xicum)
      py = 0.0d0
      pz =   rssa * dcos(xicum)
      b = ex*px + ey*py + ez*pz
      ctheta = -b/rssa
      stheta = dsqrt(1.0d0-ctheta*ctheta)
      theta_ssa = dacos(ctheta)

C  Fix phi by using constancy of scatter angle

      if ( phi_boa.eq.180.0d0 ) then
        phi_ssa = phi_all
      else if ( phi_boa .eq. 0.0d0 ) then
        phi_ssa = 0.0d0
      else
        cphi = (cosscat_up+calpha*ctheta)/stheta/salpha
        if ( cphi.gt.1.0d0) cphi = 1.0d0
        if ( cphi.lt.-1.0d0) cphi = -1.0d0
        phi_ssa     = dacos(cphi)
      endif
      
C  to degrees

      phi_ssa   = phi_ssa   / deg_to_rad
      alpha_ssa = alpha_ssa / deg_to_rad
      theta_ssa = theta_ssa / deg_to_rad
      
C  Finish

      return
      end


