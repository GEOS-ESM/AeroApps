C ###############################################################
C #                                                             #
C #                    THE VLIDORT  MODEL                       #
C #                                                             #
C #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
C #  -          --         -        -        -         -        #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, Inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #      Version 2.0. Nadir single scatter correction           #
C #                                                             #
C #            VLIDORT_SSCORR_NADIR (master)                    #
C #              VLIDORTSS_FMATRICES                            #
C #              VLIDORTSS_FMATRICES_MULTI                      #
C #                                                             #
C #      Version 2.2. outgoing sphericity correction            #
C #              2.3.  partial-layer integration                #
C #                                                             #
C #            VLIDORT_SSCORR_OUTGOING (master)                 #
C #                 SSCORR_OUTGOING_ZMATRIX                     #
C #                 outgoing_sphergeom_fine                     #
C #                 outgoing_integration_up                     #
C #                 outgoing_integration_dn                     #
C #                 multi_outgoing_adjustgeom                   #
C #                                                             #
C #      Version 2.3 outgoing Lambertian DB correction          #
C #                                                             #
C #            VLIDORT_LAMBERTIAN_DBCORR                        #
C #                                                             #
C #      Version 2.4R ------- Notes -------                     #
C #                                                             #
C #            Additional correction to ZMATRIX and FMATRIX     #
C #            routines for Azimuth > 180. HvdM, Eqs(94/95)     #
C #            Implemented by V. Natraj and R. Spurr, 5/1/09    #
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

C ###############################################################

      SUBROUTINE VLIDORT_SSCORR_NADIR (SSFLUX)

C  Single scatter exact calculation. Nadir viewing output
C   Programmed by R. Spurr, RT Solutions Inc.
C    Second Draft, April 14th 2005.
C    Third Draft,  May    6th 2005.
C       - only generates SS Stokes vector.
C       - additional code to deal with refraction.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  Argument
C  --------

      DOUBLE PRECISION SSFLUX

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1
      INTEGER          INDEX_11, INDEX_12, INDEX_34
      INTEGER          INDEX_22, INDEX_33, INDEX_44

C  help variables (double precision)

      DOUBLE PRECISION VIEWSIGN, MUX, DNM1
      DOUBLE PRECISION FINAL_SOURCE, HELP, SS_LAYERSOURCE
      DOUBLE PRECISION MULTIPLIER, SSCORRECTION
      DOUBLE PRECISION TR_CUMSOURCE, SS_CUMSOURCE
      DOUBLE PRECISION HELP2C1, HELP2S1, HELP3C1, HELP3S1

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)

      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  output Rotation cosines/sines from F-Matrix subroutine

      DOUBLE PRECISION C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S2 (MAX_GEOMETRIES,MAXLAYERS)

C  Sunlight parameters
C    - assumes the Flux vector is (F,0,0,0)
C   - drop this variables when further testing has been done

      LOGICAL           SUNLIGHT
      PARAMETER         ( SUNLIGHT = .TRUE. )
      INTEGER           NPOLAR

C  Set up operations
C  -----------------

C  Set npolar. For Natural light, there are only 3 Stokes parameters
C    ( no circular polarization for single scattering)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

C  indexing keys

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

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
        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DFLOAT(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
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
            CTHETA(N,IB) = COS_SZANGLES(IB)
            STHETA(N,IB) = SIN_SZANGLES(IB)
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

C  Upwelling single scatter Phase matrices
C  ---------------------------------------

      IF ( DO_UPWELLING ) THEN

C  Get Rotation cosine/sines and F matrices for all layers
C    Distinguish between refractive case and straightline case

        VIEWSIGN = -1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, 
     O    C1, S1, C2, S2, FMAT_UP )
        ELSE
          CALL VLIDORTSS_FMATRICES
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    C1, S1, C2, S2, FMAT_UP )
        ENDIF

C   Z matrices, Hovenier and van der Mee (1983), Eq 88.
C   ---------------------------------------------------

        IF ( NSTOKES .GT. 1 ) THEN

C  For sunlight, only need the first Column of the Z-matrix

         IF ( SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_UP(N)) THEN      
             ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11)
C  sign            ZMAT_UP(V,N,2,1) =   FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,2,1) =  -FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,3,1) =   FMAT_UP(V,N,INDEX_12) * S2(V,N)
             ZMAT_UP(V,N,4,1) =   ZERO
            ENDIF
           ENDDO
          ENDDO
         ENDIF

C  For general case. Not tested as of 15 April 2005.

         IF ( .NOT. SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_UP(N)) THEN      
             ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11)
C  sign!      ZMAT_UP(V,N,2,1) =   FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,2,1) =  -FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,3,1) =   FMAT_UP(V,N,INDEX_12) * S2(V,N)
             ZMAT_UP(V,N,4,1) =   ZERO
C  this next part remain untested. Need to check signs.
             HELP2C1 = FMAT_UP(V,N,INDEX_22) * C1(V,N)
             HELP2S1 = FMAT_UP(V,N,INDEX_22) * S1(V,N)
             HELP3C1 = FMAT_UP(V,N,INDEX_33) * C1(V,N)
             HELP3S1 = FMAT_UP(V,N,INDEX_33) * S1(V,N)
             ZMAT_UP(V,N,1,2) =   FMAT_UP(V,N,INDEX_12) * C1(V,N)
             ZMAT_UP(V,N,1,3) = - FMAT_UP(V,N,INDEX_12) * S1(V,N)
             ZMAT_UP(V,N,1,4) =   ZERO
             ZMAT_UP(V,N,2,2) =   C2(V,N) * HELP2C1 - S2(V,N) * HELP3S1
             ZMAT_UP(V,N,2,3) = - C2(V,N) * HELP2S1 - S2(V,N) * HELP3C1
             ZMAT_UP(V,N,2,4) = - FMAT_UP(V,N,INDEX_34) * S2(V,N)
             ZMAT_UP(V,N,3,2) =   S2(V,N) * HELP2C1 + C2(V,N) * HELP3S1
             ZMAT_UP(V,N,3,3) = - S2(V,N) * HELP2S1 + C2(V,N) * HELP3C1
             ZMAT_UP(V,N,3,4) =   FMAT_UP(V,N,INDEX_34) * C2(V,N)
             ZMAT_UP(V,N,4,2) = - FMAT_UP(V,N,INDEX_34) * S1(V,N)
             ZMAT_UP(V,N,4,3) = - FMAT_UP(V,N,INDEX_34) * C1(V,N)
             ZMAT_UP(V,N,4,4) =   FMAT_UP(V,N,INDEX_44)
            ENDIF
           ENDDO
          ENDDO
         ENDIF

C  end NSTOKES > 1

        ENDIF

C  debug
C         write(97,'(2i4,1pe15.7)')V,21,ZMAT_UP(V,21,1,1)
C         write(97,'(2i4,1pe15.7)')V,22,ZMAT_UP(V,22,1,1
C        pause

C  add TMS correction factor, scalar case
        
        IF ( NSTOKES .EQ. 1 ) THEN
         DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_UP(N)) THEN
           DO V = 1, N_GEOMETRIES
            ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11) * TMS(N)
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C  add TMS correction factor, vector case

        IF ( NSTOKES .GT. 1 ) THEN

C  sunlight only

         IF ( SUNLIGHT ) THEN
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_UP(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
c        if(n.eq.4)write(*,*)'b4tms',o1,ZMAT_UP(V,N,O1,1)
              ZMAT_UP(V,N,O1,1) = ZMAT_UP(V,N,O1,1) * TMS(N)
             ENDDO
            ENDDO
           ENDIF
          ENDDO

C  general case

         ELSE
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_UP(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              DO O2 = 1, NSTOKES
               ZMAT_UP(V,N,O1,O2) = ZMAT_UP(V,N,O1,O2) * TMS(N)
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF

C  finish this section

        ENDIF
      ENDIF

C  Upwelling single scatter recurrence
C  -----------------------------------

      IF ( DO_UPWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_UP(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N

C  sunlight case, no circular polarization

           IF ( SUNLIGHT ) THEN
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = EMULT_UP(UM,N,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NPOLAR
                HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* MULTIPLIER
                SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE +
     &             T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,O1,NC-1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO

C  general case

           ELSE
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = EMULT_UP(UM,N,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NSTOKES
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
                SS_LAYERSOURCE = HELP* MULTIPLIER
                SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE +
     &             T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,O1,NC-1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDIF

C  end layer loop

          ENDDO

C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct the intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C  sunlight case, no circular polarization

            IF ( SUNLIGHT ) THEN

             DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
               MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
               DO IA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + IA
                DO O1 = 1, NPOLAR
                 HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                 SS_LAYERSOURCE = HELP * MULTIPLIER
                 SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                 TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                 FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                 SSCORRECTION   = SSFLUX * FINAL_SOURCE
                 STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
             ENDDO

C  general case

            ELSE

             DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
               MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
               DO IA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + IA
                DO O1 = 1, NSTOKES
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                  HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                 ENDDO
                 SS_LAYERSOURCE = HELP * MULTIPLIER
                 SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                 TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                 FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                 SSCORRECTION   = SSFLUX * FINAL_SOURCE
                 STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
             ENDDO

            ENDIF

C  debug code to be inserted above
C                     write(98,'(a,3i4,1p2e15.7)')
C     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

C  Ongrid output :
C     Set final cumulative source and correct Stokes vector

          ELSE
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_UP(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDIF

C  insert this debug code in the above loops
C                     write(97,'(a,3i4,1p2e15.7)')
C     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  Downwelling single scatter Phase matrices
C  -----------------------------------------

      IF ( DO_DNWELLING ) THEN

C  Get Rotation cosine/sines and F matrices for all layers
C    Distinguish between refractive case and straightline case

        VIEWSIGN = +1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    C1, S1, C2, S2, FMAT_DN )
        ELSE
          CALL VLIDORTSS_FMATRICES
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, 
     O    C1, S1, C2, S2, FMAT_DN )
        ENDIF

C   Z matrices, Hovenier and van der Mee (1983), Eq 88.

        IF ( NSTOKES .GT. 1 ) THEN

C  For sunlight, only need the first column of Z matrix

         IF ( SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_DN(N)) THEN      
             ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11)
C  sign            ZMAT_DN(V,N,2,1) =   FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,2,1) =   -FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,3,1) =   FMAT_DN(V,N,INDEX_12) * S2(V,N)
             ZMAT_DN(V,N,4,1) =   ZERO
            ENDIF
           ENDDO
          ENDDO
         ENDIF

C  For general case. Not tested as of 15 April 2005.

         IF ( .NOT. SUNLIGHT ) THEN

          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_DN(N)) THEN
C  this code tested
             ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11)
C sign            ZMAT_DN(V,N,2,1) =   FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,2,1) =   -FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,3,1) =   FMAT_DN(V,N,INDEX_12) * S2(V,N)
             ZMAT_DN(V,N,4,1) =   ZERO
C  This code untested
             HELP2C1 = FMAT_DN(V,N,INDEX_22) * C1(V,N)
             HELP2S1 = FMAT_DN(V,N,INDEX_22) * S1(V,N)
             HELP3C1 = FMAT_DN(V,N,INDEX_33) * C1(V,N)
             HELP3S1 = FMAT_DN(V,N,INDEX_33) * S1(V,N)
             ZMAT_DN(V,N,1,2) =   FMAT_DN(V,N,INDEX_12) * C1(V,N)
             ZMAT_DN(V,N,1,3) = - FMAT_DN(V,N,INDEX_12) * S1(V,N)
             ZMAT_DN(V,N,1,4) =   ZERO
             ZMAT_DN(V,N,2,2) =   C2(V,N) * HELP2C1 - S2(V,N) * HELP3S1
             ZMAT_DN(V,N,2,3) = - C2(V,N) * HELP2S1 - S2(V,N) * HELP3C1
             ZMAT_DN(V,N,2,4) = - FMAT_DN(V,N,INDEX_34) * S2(V,N)
             ZMAT_DN(V,N,3,2) =   S2(V,N) * HELP2C1 + C2(V,N) * HELP3S1
             ZMAT_DN(V,N,3,3) = - S2(V,N) * HELP2S1 + C2(V,N) * HELP3C1
             ZMAT_DN(V,N,3,4) =   FMAT_DN(V,N,INDEX_34) * C2(V,N)
             ZMAT_DN(V,N,4,2) = - FMAT_DN(V,N,INDEX_34) * S1(V,N)
             ZMAT_DN(V,N,4,3) = - FMAT_DN(V,N,INDEX_34) * C1(V,N)
             ZMAT_DN(V,N,4,4) =   FMAT_DN(V,N,INDEX_44)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDIF

C  add TMS correction factor, scalar case
        
        IF ( NSTOKES .EQ. 1 ) THEN
         DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_DN(N)) THEN
           DO V = 1, N_GEOMETRIES
            ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11) * TMS(N)
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C  add TMS correction factor, vector case

        IF ( NSTOKES .GT. 1 ) THEN
         IF ( SUNLIGHT ) THEN

C  sunlight only

          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_DN(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              ZMAT_DN(V,N,O1,1) = ZMAT_DN(V,N,O1,1) * TMS(N)
             ENDDO
            ENDDO
           ENDIF
          ENDDO

C  general case

         ELSE
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_DN(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              DO O2 = 1, NSTOKES
               ZMAT_DN(V,N,O1,O2) = ZMAT_DN(V,N,O1,O2) * TMS(N)
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF

C  finish this section

        ENDIF
      ENDIF

C  Downwelling single scatter recurrence
C  -------------------------------------

      IF ( DO_DNWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_DN(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT
           NC = N

C  sunlight case

           IF ( SUNLIGHT ) THEN

            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = EMULT_DN(UM,N,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NPOLAR
                HELP =  ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* MULTIPLIER
                SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE +
     &             T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,O1,NC-1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO

C  general case

           ELSE

            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = EMULT_DN(UM,N,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NSTOKES
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
                SS_LAYERSOURCE = HELP* MULTIPLIER
                SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE +
     &             T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,O1,NC-1)
               ENDDO
              ENDDO
             ENDDO
            ENDDO

           ENDIF

C  end layer loop

          ENDDO

C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct the intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C  sunlight case

           IF ( SUNLIGHT ) THEN

            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NSTOKES
                HELP = ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP * MULTIPLIER
                SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
               ENDDO
              ENDDO
             ENDDO
            ENDDO

C  general case

           ELSE

            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
              DO IA = 1, N_USER_RELAZMS
               V = VZA_OFFSETS(IB,UM) + IA
               DO O1 = 1, NSTOKES
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
                SS_LAYERSOURCE = HELP * MULTIPLIER
                SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
               ENDDO
              ENDDO
             ENDDO
            ENDDO

           ENDIF

C  debug code to be inserted
C                     write(98,'(a,3i4,1p2e15.7)')
C     &                 'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

C  Ongrid output :
C     Set final cumulative source and correct Stokes vector

          ELSE

           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO IA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + IA
              DO O1 = 1, NSTOKES
               FINAL_SOURCE = SS_CUMSOURCE_DN(V,O1,NC)
               SSCORRECTION = SSFLUX * FINAL_SOURCE
               STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
              ENDDO
             ENDDO
            ENDDO
           ENDDO

C  debug code to be inserted above
C         write(97,'(a,3i4,1p2e15.7)')
C     &        'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

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

C

      SUBROUTINE VLIDORTSS_FMATRICES
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEKMOMS, 
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VSIGN, VZA_OFFSETS,
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI,
     O    C1, S1, C2, S2, FMAT )

C  include file of dimensions and numbers
C     Cannot use bookkeeping vars file

      INCLUDE '../includes/VLIDORT.PARS'

C  input
C  -----

C  Additional control for the truncation

      LOGICAL          DO_SSCORR_TRUNCATION

C  control integers

      INTEGER          NLAYERS, N_GEOMETRIES, NGREEKMOMS
      INTEGER          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS

C  +/- sign, offsets

      DOUBLE PRECISION VSIGN 
      INTEGER          VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

C  scattering input information

      INTEGER          LAYER_MAXMOMENTS(MAXLAYERS)
      DOUBLE PRECISION GREEKMAT_TOTAL_INPUT
     &      ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION SSFDEL( MAXLAYERS )

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  azimuth angles. Added V2.4R. for PHI > 180 case

      DOUBLE PRECISION PHI    (MAX_USER_RELAZMS)

C  output
C  ------

C  Rotation Angle cosines/sines

      DOUBLE PRECISION C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S2 (MAX_GEOMETRIES,MAXLAYERS)

C  F-Matrix

      DOUBLE PRECISION FMAT ( MAX_GEOMETRIES, MAXLAYERS, 6 )

C  Local
C  -----

      INTEGER          IB, UM, IA, V, N, GREEKMAT_INDEX(6)
      DOUBLE PRECISION COSSCAT, SINSCAT, HELP_SINSCAT, F, FT, DNL1
      DOUBLE PRECISION CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION UUU(MAX_GEOMETRIES)
      DOUBLE PRECISION P00(MAX_GEOMETRIES,2)
      DOUBLE PRECISION P02(MAX_GEOMETRIES,2)

      INTEGER          K, L, LNEW, LOLD, ITMP
      INTEGER          INDEX_11, INDEX_12, INDEX_34
      INTEGER          INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION DL, QROOT6, FAC1, FAC2, SQL4, SQL41, TMP1, TMP2
      DOUBLE PRECISION GK11(MAXLAYERS), GK12(MAXLAYERS)
      DOUBLE PRECISION GK44(MAXLAYERS), GK34(MAXLAYERS)

C  indexing key

C      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
C      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
C      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
C      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
C      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
C      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      SQL41 = ZERO

C  Geometrical quantities
C  ----------------------

      DO IB = 1, NBEAMS
        DO UM = 1, N_USER_STREAMS
          DO IA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + IA

C  cosine scatter angle (this is valid only for non-refracting atmosphere)
C  VSIGN = -1 for upwelling, +1 for downwelling

            COSSCAT = VSIGN * CTHETA(1,IB) * CALPHA(UM) +
     &                        STHETA(1,IB) * SALPHA(UM) * CPHI(IA)
            UUU(V)  = COSSCAT

C  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

C  a. safety
C    Watch for sin^2(scatter angle) less than zero (machine precision)
C    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

            HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
            IF ( HELP_SINSCAT.LE.ZERO ) THEN
              SINSCAT = 1.0D-12
            ELSE
              SINSCAT = DSQRT ( HELP_SINSCAT )
            ENDIF

C  b. necessary limit analyses - Hovenier limits.
C     R. Spurr and V. Natraj, 17 January 2006

            IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
              CSIG1 = ZERO
              CSIG2 = ZERO
            ELSE
             IF ( STHETA(1,IB) .EQ. ZERO ) THEN
              CSIG1 = -  CPHI(IA)
             ELSE
              CSIG1   = ( - VSIGN * CALPHA(UM) + CTHETA(1,IB)*COSSCAT )
     &                  / SINSCAT / STHETA(1,IB)
             ENDIF
             IF ( SALPHA(UM) .EQ. ZERO ) THEN
               CSIG2 = -  CPHI(IA)
             ELSE
               CSIG2   = ( - CTHETA(1,IB) + VSIGN*CALPHA(UM)*COSSCAT )
     &                    / SINSCAT / SALPHA(UM)
             ENDIF
            ENDIF
c    write(7,'(4i5,1p4e15.7)')V,IB,UM,IA,CSIG1,CSIG2,sinscat
            IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
            IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

C  output, H/VdM, Eqs. (89)-(94)
C    Rotation sines and cosines. Same for all layers

            CSIG1_2 = TWO * CSIG1
            CSIG2_2 = TWO * CSIG2
            IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
             SSIG1 = ZERO
            ELSE
             SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
            ENDIF
            IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
             SSIG2 = ZERO
            ELSE
             SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
            ENDIF

C  For relazm in [180,360), need sign reversal for S1 and S2
C  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

            C1(V,1) = CSIG1_2 * CSIG1 - ONE
            C2(V,1) = CSIG2_2 * CSIG2 - ONE

            IF (PHI(IA) .LE. 180.D0) THEN
              S1(V,1) = CSIG1_2 * SSIG1
              S2(V,1) = CSIG2_2 * SSIG2
            ELSE
              S1(V,1) = -CSIG1_2 * SSIG1
              S2(V,1) = -CSIG2_2 * SSIG2
            ENDIF

C  Copy to all other layers

            DO N = 2, NLAYERS
              C1(V,N) = C1(V,1)
              S1(V,N) = S1(V,1)
              C2(V,N) = C2(V,1)
              S2(V,N) = S2(V,1)
            ENDDO
c    write(8,'(i5,1p4e15.7)')V,C1(V,1),S1(V,1),C2(V,1),S2(V,1)

C  End geometry loops

          ENDDO
        ENDDO
      ENDDO

C  F-matrices
C  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

C initialise F-matrix

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO N = 1, NLAYERS
            FMAT(V,N,K) = ZERO
          END DO
        END DO
      END DO

C  Start loop over the coefficient index l
C  first update generalized spherical functions, then calculate coefs.
C  lold and lnew are pointer-like indices used in recurrence 

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)

C  Set the local Greek matrix elements that you need = 11 and 12.
c   44 and 34 are not required with natural sunlight (default here)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DO N = 1, NLAYERS
           DNL1 = DBLE(2*L + 1 )
           F    = SSFDEL(N) * DNL1
           FT   = ONE - SSFDEL(N)
           GK11(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))-F)/FT
           GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
           GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
           GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
           GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
           GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
           GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
          ENDDO
        ENDIF

C  First moment

        IF ( L .EQ. 0 ) THEN

C  Adding paper Eqs. (76) and (77) with m=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = ONE
            P00(V,LNEW) = ZERO
            P02(V,LOLD) = ZERO
            P02(V,LNEW) = ZERO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = (DL-ONE)/DL

C Adding paper Eq. (81) with m=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = FAC1*UUU(V)*P00(V,LNEW) - FAC2*P00(V,LOLD)
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = QROOT6*(ONE-UUU(V)*UUU(V))
            P02(V,LNEW) = ZERO
          END DO
          SQL41 = ZERO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = TMP1*UUU(V)*P02(V,LNEW) - TMP2*P02(V,LOLD)
          END DO

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

        DO N = 1, NLAYERS
         IF ( L.LE.LAYER_MAXMOMENTS(N) ) THEN
          DO V = 1, N_GEOMETRIES
           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11) + GK11(N)*P00(V,LNEW)
           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12) + GK12(N)*P02(V,LNEW)
           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44) + GK44(N)*P00(V,LNEW)
           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34) + GK34(N)*P02(V,LNEW)

C  Previous code...............................................
c           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11) +
c     &       GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))*P00(V,LNEW)
c           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12) +
c     &       GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))*P02(V,LNEW)
c           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44) +
c     &       GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))*P00(V,LNEW)
c           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34) +
c     &       GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))*P02(V,LNEW)
          ENDDO
         ENDIF
        END DO

C  end moment loop

      END DO

C  remaining symmetries for Mie particles

      DO V = 1, N_GEOMETRIES
        DO N = 1, NLAYERS
          FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_11)
          FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_44)
        END DO
      END DO

C  finish

      RETURN
      END

C

      SUBROUTINE VLIDORTSS_FMATRICES_MULTI
     I  ( DO_SSCORR_TRUNCATION, NLAYERS,
     I    N_GEOMETRIES, NGREEKMOMS, 
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL,
     I    VSIGN, VZA_OFFSETS,
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, 
     O    C1, S1, C2, S2, FMAT )

C  F-matrix routine with multiple-solar zenith angles (one for each layer)
C  Applicable with the Refractive Geometry case.

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  input
C  -----

C  Additional control for the truncation

      LOGICAL          DO_SSCORR_TRUNCATION

C  control integers

      INTEGER          NLAYERS, N_GEOMETRIES, NGREEKMOMS
      INTEGER          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS

C  +/- sign, offsets

      DOUBLE PRECISION VSIGN 
      INTEGER          VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

C  scattering input information

      INTEGER          LAYER_MAXMOMENTS(MAXLAYERS)
      DOUBLE PRECISION GREEKMAT_TOTAL_INPUT
     &      ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION SSFDEL( MAXLAYERS )

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  azimuth angles. Added V2.4R. for PHI > 180 case

      DOUBLE PRECISION PHI    (MAX_USER_RELAZMS)

C  output
C  ------

C  Rotation Angle cosines/sines

      DOUBLE PRECISION C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION S2 (MAX_GEOMETRIES,MAXLAYERS)

C  F-Matrix

      DOUBLE PRECISION  FMAT ( MAX_GEOMETRIES, MAXLAYERS, 6 )

C  Local
C  -----

      INTEGER          IB, UM, IA, V, N, GREEKMAT_INDEX(6)
      DOUBLE PRECISION COSSCAT, SINSCAT, HELP_SINSCAT, F, FT, DNL1
      DOUBLE PRECISION CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

C  now have layer dimensioning in the help arrays

      DOUBLE PRECISION UUU(MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION P00(MAX_GEOMETRIES,MAXLAYERS,2)
      DOUBLE PRECISION P02(MAX_GEOMETRIES,MAXLAYERS,2)

      INTEGER          K, L, LNEW, LOLD, ITMP
      INTEGER          INDEX_11, INDEX_12, INDEX_34
      INTEGER          INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION DL, QROOT6, FAC1, FAC2, SQL4, SQL41, TMP1, TMP2
      DOUBLE PRECISION GK11(MAXLAYERS), GK12(MAXLAYERS)
      DOUBLE PRECISION GK44(MAXLAYERS), GK34(MAXLAYERS)

C  indexing key

C      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
C      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
C      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
C      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
C      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
C      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      SQL41 = ZERO

C  Geometrical quantities
C  ----------------------

      DO IB = 1, NBEAMS
        DO UM = 1, N_USER_STREAMS
          DO IA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + IA
            DO N = 1, NLAYERS

C  cosine scatter angle (this is valid only for non-refracting atmosphere)
C  VSIGN = -1 for upwelling, +1 for downwelling

              COSSCAT = VSIGN * CTHETA(N,IB) * CALPHA(UM) +
     &                          STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
              UUU(V,N)  = COSSCAT

C  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

C  a. safety
C    Watch for sin^2(scatter angle) less than zero (machine precision)
C    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

              HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
              IF ( HELP_SINSCAT.LE.ZERO ) THEN
                SINSCAT = 1.0D-12
              ELSE
                SINSCAT = DSQRT ( HELP_SINSCAT )
              ENDIF

C  b. necessary limit analyses - Hovenier limits.
C     R. Spurr and V. Natraj, 17 January 2006

              IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
               CSIG1 = ZERO
               CSIG2 = ZERO
              ELSE
               IF ( STHETA(1,IB) .EQ. ZERO ) THEN
                CSIG1 = -  CPHI(IA)
               ELSE
                CSIG1 = ( - VSIGN*CALPHA(UM) + CTHETA(N,IB)*COSSCAT )
     &                       / SINSCAT / STHETA(N,IB)
               ENDIF
               IF ( SALPHA(UM) .EQ. ZERO ) THEN
                CSIG2 = -  CPHI(IA)
               ELSE
                CSIG2 = ( - CTHETA(N,IB) + VSIGN*CALPHA(UM)*COSSCAT )
     &                       / SINSCAT / SALPHA(UM)
               ENDIF
              ENDIF
              IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
              IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

C  output, H/VdM, Eqs. (89)-(94)

              CSIG1_2 = TWO * CSIG1
              CSIG2_2 = TWO * CSIG2
              IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
               SSIG1 = ZERO
              ELSE
               SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
              ENDIF
              IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
               SSIG2 = ZERO
              ELSE
               SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
              ENDIF

C  For relazm in [180,360), need sign reversal for S1 and S2
C  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

              C1(V,N) = CSIG1_2 * CSIG1 - ONE
              C2(V,N) = CSIG2_2 * CSIG2 - ONE

              IF (PHI(IA) .LE. 180.D0) THEN
                S1(V,N) = CSIG1_2 * SSIG1
                S2(V,N) = CSIG2_2 * SSIG2
              ELSE
                S1(V,N) = -CSIG1_2 * SSIG1
                S2(V,N) = -CSIG2_2 * SSIG2
              ENDIF

C  End layer loop

            ENDDO

C  End geometry loops

          ENDDO
        ENDDO
      ENDDO

C  F-matrices
C  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

C initialise F-matrix

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO N = 1, NLAYERS
            FMAT(V,N,K) = ZERO
          END DO
        END DO
      END DO

C  Start loop over the coefficient index l
C  first update generalized spherical functions, then calculate coefs.
C  lold and lnew are pointer-like indices used in recurrence 

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)

C  Set the local Greek matrix elements that you need = 11 and 12.
c   44 and 34 are not required with natural sunlight (default here)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DO N = 1, NLAYERS
           DNL1 = DBLE(2*L + 1 )
           F    = SSFDEL(N) * DNL1
           FT   = ONE - SSFDEL(N)
           GK11(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))-F)/FT
           GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
           GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
           GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
           GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
           GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
           GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
          ENDDO
        ENDIF

        IF ( L .EQ. 0 ) THEN

C  Adding paper Eqs. (76) and (77) with m=0

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P00(V,N,LOLD) = ONE
              P00(V,N,LNEW) = ZERO
              P02(V,N,LOLD) = ZERO
              P02(V,N,LNEW) = ZERO
            END DO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = (DL-ONE)/DL

C Adding paper Eq. (81) with m=0

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P00(V,N,LOLD) = FAC1*UUU(V,N)*P00(V,N,LNEW)
     &                        - FAC2*P00(V,N,LOLD)
            END DO
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P02(V,N,LOLD) = QROOT6*(ONE-UUU(V,N)*UUU(V,N))
              P02(V,N,LNEW) = ZERO
            END DO
          END DO
          SQL41 = ZERO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P02(V,N,LOLD) = TMP1*UUU(V,N)*P02(V,N,LNEW)
     &                        - TMP2*P02(V,N,LOLD)
            END DO
          END DO

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

        DO N = 1, NLAYERS
         IF ( L.LE.LAYER_MAXMOMENTS(N) ) THEN
          DO V = 1, N_GEOMETRIES
           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11)+GK11(N)*P00(V,N,LNEW)
           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12)+GK12(N)*P02(V,N,LNEW)
           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44)+GK44(N)*P00(V,N,LNEW)
           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34)+GK34(N)*P02(V,N,LNEW)
C  Previous code...............................................
c           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11) +
c     &    GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))*P00(V,N,LNEW)
c           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12) +
c     &    GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))*P02(V,N,LNEW)
c           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44) +
c     &    GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))*P00(V,N,LNEW)
c           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34) +
c     &    GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))*P02(V,N,LNEW)
          ENDDO
         ENDIF
        END DO

C  end moment loop

      END DO

C  remaining symmetries for Mie particles

      DO V = 1, N_GEOMETRIES
        DO N = 1, NLAYERS
          FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_11)
          FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_44)
        END DO
      END DO

C  finish

      RETURN
      END

c

      SUBROUTINE VLIDORT_SSCORR_OUTGOING
     &        ( SSFLUX, FAIL, MESSAGE)

C  Single scatter exact calculation for the outgoing LOS
C         - NEW for Version 2.2

C   Programmed by R. Spurr, RT Solutions Inc.
C    First Draft, January 30th 2007.
C    Validated against TOMRAD. 29 March 2007.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

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

      integer          ntraverse_ut(max_partlayers)
      double precision sunpaths_ut (max_partlayers,maxlayers)
      double precision radii_ut    (max_partlayers)
      double precision alpha_ut    (max_partlayers)

C  Other (incidental) geometrical output

      double precision lospaths(maxlayers)
      double precision lospaths_ut_up(max_partlayers)
      double precision lospaths_ut_dn(max_partlayers)
      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

C  Extinction

      double precision extinction ( maxlayers)

C  Partial layer heights

      double precision height_grid_ut(max_partlayers)

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, IA, NC, IB, V, O1, NM1

C  help variables (double precision)

      DOUBLE PRECISION FINAL_SOURCE, HELP, SS_CUMSOURCE
      DOUBLE PRECISION SS_LAYERSOURCE, SSCORRECTION, XT
      DOUBLE PRECISION CTHETA, STHETA, CALPHA, SALPHA, CPHI, CSA, TRANS
      DOUBLE PRECISION VSIGN, ZMAT_LOCAL(4), FMAT_LOCAL(6), DNM1

C  Sunlight parameters
C    - assumes the Flux vector is (F,0,0,0)
C   - drop this variables when further testing has been done

      LOGICAL           SUNLIGHT
      PARAMETER         ( SUNLIGHT = .TRUE. )
      INTEGER           NPOLAR

C  Set up operations
C  -----------------

C  Set up partials flag

      do_fine = .true.
      do_partials = ( n_partlayers .gt. 0 )

C  Set npolar. For Natural light, there are only 3 Stokes parameters
C    ( no circular polarization for single scattering)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

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
        do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          xt = deltau_vert(n) - partau_vert(ut)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.
C   TMS gets modified by (1-F). Save the truncation factor.
C   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DFLOAT(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

C  Source term calculations
C  ========================

C  Start the main loop over all solar and viewing geometries
C   Use the adjusted values of the angles

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_VZANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = SZANGLES_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = VZA_OFFSETS(IB,UM) + IA
    
C  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine
     i  ( maxlayers, maxfinelayers, max_partlayers,
     i    do_fine, do_partials, nlayers, nfinelayers,
     i    n_partlayers, PARTLAYERS_layeridx,
     i    height_grid, height_grid_ut, earth_radius,
     i    alpha_boa, theta_boa, phi_boa,
     o    sunpaths,      radii,      ntraverse,      alpha_all, 
     o    sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,
     o    sunpaths_ut,   radii_ut,   ntraverse_ut,   alpha_ut,
     o    PARTLAYERS_layerfineidx, lospaths,
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
     i             do_partials, extinction,
     i             n_partlayers, partlayers_layeridx,
     i             partlayers_layerfineidx,
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
c               pause

c              if (v.eq.1)write(35,*)v
c              do ut = 1, n_partlayers
c                if (v.eq.1)write(35,'(1p2e18.10)')
c     &           up_lostrans_ut(ut,v),up_multipliers_ut(ut,v)
c              enddo

C  Debug: Multipliers (all layers), Upwelling OLD WAY
c              call outgoing_integration_old_up
c     i           ( nlayers, extinction, deltau_vert,
c     i             lospaths, sunpaths, ntraverse,
c     o             up_multipliers(1,v), up_lostrans(1,v) )
c              write(46,*)v
c              do n = 1, nlayers
c                write(46,*)up_lostrans(n,v),up_multipliers(n,v)
c              enddo

C  Start layer loop

              DO N = NLAYERS, 1, -1

C  Trigonometry

                CTHETA = DCOS(THETA_ALL(N))
                STHETA = DSIN(THETA_ALL(N))
                CALPHA = DCOS(ALPHA_ALL(N))
                SALPHA = DSIN(ALPHA_ALL(N))
                CPHI   = DCOS(PHI_ALL(N))
                CSA    = COSSCAT_UP(N)
                VSIGN  = -1.0d0

C  Call to scattering law for phase matrix ZMAT

                CALL SSCORR_OUTGOING_ZMATRIX
     I          ( DO_SSCORR_TRUNCATION, 
     I            N, NSTOKES, NGREEK_MOMENTS_INPUT, SSFDEL(N),
     I            LAYER_MAXMOMENTS(N), GREEKMAT_TOTAL_INPUT,
     I            CTHETA, STHETA, CALPHA, SALPHA, CPHI,
     I            PHI_ALL(N), CSA, VSIGN,
     O            ZMAT_LOCAL, FMAT_LOCAL )

C  Phase matrix (multiplied by TMS factor). Save them.
C     Sunlight only, the first column of the matrix

                IF ( STERM_LAYERMASK_UP(N) ) THEN
                 DO O1 = 1, NSTOKES
                  ZMAT_UP(V,N,O1,1) = ZMAT_LOCAL(O1) * TMS(N)
                 ENDDO
                ENDIF

C  Finish the layer loop and end upwelling

              ENDDO
            ENDIF

C  Downwelling Source terms: Phase matrix, multipliers, transmittances
C  -------------------------------------------------------------------

            IF ( DO_DNWELLING ) THEN

C  Multipliers and transmittances

              call outgoing_integration_dn
     i           ( nlayers, nfinelayers,
     i             do_partials, extinction,
     i             n_partlayers, partlayers_layeridx,
     i             partlayers_layerfineidx,
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
c     i             lospaths, sunpaths, ntraverse,
c     o             dn_multipliers(1,v), dn_lostrans(1,v) )
c              write(66,*)v
c              do n = 1, nlayers
c                write(66,*)dn_lostrans(n,v),dn_multipliers(n,v)
c              enddo

C  start layer loop

              DO N = 1, NLAYERS

C  Trigonometry

                CTHETA = DCOS(THETA_ALL(N-1))
                STHETA = DSIN(THETA_ALL(N-1))
                CALPHA = DCOS(ALPHA_ALL(N-1))
                SALPHA = DSIN(ALPHA_ALL(N-1))
                CPHI   = DCOS(PHI_ALL(N-1))
                CSA    = COSSCAT_DN(N-1)
                VSIGN  = +1.0d0

C  Call to scattering law for phase matrix ZMAT

                CALL SSCORR_OUTGOING_ZMATRIX
     I          ( DO_SSCORR_TRUNCATION, 
     I            N, NSTOKES, NGREEK_MOMENTS_INPUT, SSFDEL(N),
     I            LAYER_MAXMOMENTS(N), GREEKMAT_TOTAL_INPUT,
     I            CTHETA, STHETA, CALPHA, SALPHA, CPHI, 
     I            PHI_ALL(N-1), CSA, VSIGN,
     O            ZMAT_LOCAL, FMAT_LOCAL )

C  Phase matrix (multiplied by TMS factor). Save them.
C     Sunlight only, the first column of the matrix

                IF ( STERM_LAYERMASK_DN(N) ) THEN
                 DO O1 = 1, NSTOKES
                  ZMAT_DN(V,N,O1,1) = ZMAT_LOCAL(O1) * TMS(N)
                 ENDDO
                ENDIF

C  End layer loop and downwelling

              ENDDO
            ENDIF

C   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

C  Recurrence relation for the UPWELLING Stokes vector
C  ===================================================

      IF ( DO_UPWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_UP(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      sunlight case, no circular polarization

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* UP_MULTIPLIERS(N,V)
                SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE +
     &             UP_LOSTRANS(N,V)*SS_CUMSOURCE_UP(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

C  sunlight case, no circular polarization
C  Offgrid output-------
C    Add additional partial layer source term = Exact Phase Func * Multiplier
C     Set final cumulative source and single scatter Stokes vector
C  Ongrid output--------
C     Set final cumulative source and single scatter Stokes vector

          IF ( partlayers_OUTFLAG(UTA) ) THEN
            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
                SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                TRANS = UP_LOSTRANS_UT(UT,V)
                FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_UP(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  Recurrence relation for the DOWNWELLING Stokes vector
C  =====================================================

      IF ( DO_DNWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_DN(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      sunlight case only

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP =  ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* DN_MULTIPLIERS(N,V)
                SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE +
     &             DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

C    sunlight case
C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C     Set final cumulative source and Correct Stokes vector
C  Ongrid output :
C     Set final cumulative source and Correct Stokes vector

          IF ( partlayers_OUTFLAG(UTA) ) THEN
            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
                SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                TRANS = DN_LOSTRANS_UT(UT,V)
                FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_DN(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
              ENDDO
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

C

      SUBROUTINE SSCORR_OUTGOING_ZMATRIX
     I  ( DO_SSCORR_TRUNCATION,
     I    LAYER, NSTOKES, NGREEKMOMS, TRUNCFACTOR,
     I    LAYER_MAXMOMENTS, GREEKMATRIX,
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, 
     I    PHI, COSSCAT, VSIGN,  
     O    ZMAT, FMAT )

C  Include file
C  ------------

      INCLUDE '../includes/VLIDORT.PARS'

C  input
C  -----

C  Control for using the local truncation

      LOGICAL          DO_SSCORR_TRUNCATION

C  control integers

      INTEGER          LAYER, NSTOKES, NGREEKMOMS

C  scattering input information

      INTEGER          LAYER_MAXMOMENTS
      DOUBLE PRECISION GREEKMATRIX (0:MAXMOMENTS_INPUT,MAXLAYERS,16)
      DOUBLE PRECISION TRUNCFACTOR

C  zenith angle, azimuth angle, scatter angle (sines/cosines)

      DOUBLE PRECISION CTHETA
      DOUBLE PRECISION STHETA
      DOUBLE PRECISION CALPHA
      DOUBLE PRECISION SALPHA
      DOUBLE PRECISION CPHI
      DOUBLE PRECISION COSSCAT
      DOUBLE PRECISION VSIGN

C  Azimuth angle. Added for Version 2.4R, required for PHI > 180.
C    R. Spurr and V. Natraj, 01 May 2009

      DOUBLE PRECISION PHI

C  output
C  ------

C  Z-Matrix, first column only. F-matrix
C    (this module only works with sunlight)

      DOUBLE PRECISION  ZMAT ( 4 )
      DOUBLE PRECISION  FMAT ( 6 )

C  Local
C  -----
       
      DOUBLE PRECISION  C1, S1, C2, S2, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION  CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2
      DOUBLE PRECISION  DL, QROOT6, UUU, P00(2), P02(2)
      DOUBLE PRECISION  FAC1, FAC2, SQL4, SQL41, TMP1, TMP2
      DOUBLE PRECISION  DNL1, FDNL1,  FACT, GK_11, GK_12

      INTEGER           GREEKMAT_INDEX(6), K, L, LNEW, LOLD, ITMP, N
      INTEGER           INDEX_11, INDEX_12, INDEX_34
      INTEGER           INDEX_22, INDEX_33, INDEX_44

C  indexing key

C      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
C      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
C      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
C      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
C      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
C      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      N = LAYER

C  Geometrical quantities
C  ----------------------

C  cosine scatter angle (this is valid only for non-refracting atmosphere)
C  VSIGN = -1 for upwelling, +1 for downwelling
C  Should be the same as the input value CSA

      UUU  = COSSCAT

C  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

C  a. safety
C    Watch for sin^2(scatter angle) less than zero (machine precision)
C    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

      HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
      IF ( HELP_SINSCAT.LE.ZERO ) THEN
        SINSCAT = 1.0D-12
      ELSE
        SINSCAT = DSQRT ( HELP_SINSCAT )
      ENDIF

C  b. necessary limit analyses - Hovenier limits.
C     R. Spurr and V. Natraj, 17 January 2006

      IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
        CSIG1 = ZERO
        CSIG2 = ZERO
      ELSE
        IF ( STHETA .EQ. ZERO ) THEN
          CSIG1 = -  CPHI
        ELSE
          CSIG1 = (-VSIGN*CALPHA+CTHETA*COSSCAT)/SINSCAT/STHETA
        ENDIF
        IF ( SALPHA .EQ. ZERO ) THEN
          CSIG2 = -  CPHI
        ELSE
          CSIG2 = (-CTHETA+VSIGN*CALPHA*COSSCAT)/SINSCAT/SALPHA
        ENDIF
      ENDIF
      IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
      IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

C  output, H/VdM, Eqs. (89)-(94)

      CSIG1_2 = TWO * CSIG1
      CSIG2_2 = TWO * CSIG2
      IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
        SSIG1 = ZERO
      ELSE
        SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
      ENDIF
      IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
        SSIG2 = ZERO
      ELSE
        SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
      ENDIF

C  For relazm in [180,360), need sign reversal for S1 and S2
C  See H/VdM, Eqs. 94-95

      C1 = CSIG1_2 * CSIG1 - ONE
      C2 = CSIG2_2 * CSIG2 - ONE
      IF ( PHI .LE. 180.0d0 ) THEN
        S1 = CSIG1_2 * SSIG1
        S2 = CSIG2_2 * SSIG2
      ELSE
        S1 = -CSIG1_2 * SSIG1
        S2 = -CSIG2_2 * SSIG2
      ENDIF

C  F-matrices
C  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

C initialise F-matrix

      DO K = 1, 6
        FMAT(K) = 0.0D0
      END DO

C  Start loop over the coefficient index l
C  first update generalized spherical functions, then calculate coefs.
C  lold and lnew are pointer-like indices used in recurrence 

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

C  Set the local Greek matrix elements that you need = 11 and 12.
c   44 and 34 are not required with natural sunlight (default here)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DNL1  = DBLE(2*L + 1 )
          FDNL1 = TRUNCFACTOR * DNL1
          FACT  = ONE - TRUNCFACTOR
          GK_11 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(1)) - FDNL1 ) / FACT
          GK_12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
c          GK_44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
c          GK_34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
        ELSE
          GK_11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
          GK_12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
c          GK_44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
c          GK_34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
        ENDIF

C  First moment

        IF ( L .EQ. 0 ) THEN

C  Adding paper Eqs. (76) and (77) with m=0

          P00(LOLD) = ONE
          P00(LNEW) = ZERO
          P02(LOLD) = ZERO
          P02(LNEW) = ZERO

        ELSE

          DL   = DBLE(L)
          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = (DL-ONE)/DL

C Adding paper Eq. (81) with m=0

          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          P02(LOLD) = QROOT6*(ONE-UUU*UUU)
          P02(LNEW) = ZERO
          SQL41 = ZERO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-4.0D0)
          TMP1  = (2.0D0*DL-1.0D0)/SQL41
          TMP2  = SQL4/SQL41
          P02(LOLD) = TMP1*UUU*P02(LNEW) - TMP2*P02(LOLD)

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

        IF ( L.LE.LAYER_MAXMOMENTS ) THEN
          FMAT(INDEX_11) = FMAT(INDEX_11) + GK_11 * P00(LNEW)
          FMAT(INDEX_12) = FMAT(INDEX_12) + GK_12 * P02(LNEW)
c          FMAT(INDEX_44) = FMAT(INDEX_44) + GK_44 * P00(LNEW)
c          FMAT(INDEX_34) = FMAT(INDEX_34) + GK_34 * P02(LNEW)

C  Previous code............................
C          FMAT(INDEX_11) = FMAT(INDEX_11) +
c     &         GREEKMATRIX(L,N,GREEKMAT_INDEX(1))*P00(LNEW)
c          FMAT(INDEX_12) = FMAT(INDEX_12) +
c     &         GREEKMATRIX(L,N,GREEKMAT_INDEX(3))*P02(LNEW)
C  do not need these with natural sunlight.............
c            FMAT(INDEX_44) = FMAT(INDEX_44) +
c     &        GREEKMATRIX(L,N,GREEKMAT_INDEX(6))*P00(LNEW)
c            FMAT(INDEX_34) = FMAT(INDEX_34) +
c     &        GREEKMATRIX(L,N,GREEKMAT_INDEX(5))*P02(LNEW)
        ENDIF

C  End moments

      END DO

C  remaining symmetries for Mie particles. Not needed

C      FMAT(INDEX_22) = FMAT(INDEX_11)
C      FMAT(INDEX_33) = FMAT(INDEX_44)

C  Z-matrix, first column only

      IF ( NSTOKES .EQ. 1 ) THEN
        ZMAT(1) =   FMAT(INDEX_11)
      ELSE
        ZMAT(1) =   FMAT(INDEX_11)
C  sign        ZMAT(2) =   FMAT(INDEX_12) * C2
C  sign        ZMAT(2) =  -FMAT(INDEX_12) * C2
        ZMAT(2) =  -FMAT(INDEX_12) * C2
        ZMAT(3) =   FMAT(INDEX_12) * S2
        ZMAT(4) =   0.0D0
      ENDIF

C  finish

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
          do while (radii_p(ut).gt.radii_fine(np,j))
             j = j + 1
             if (j.gt.nfine) go to 5778
          enddo
 5778     continue
          partials_fineidx(ut) = j - 1
        enddo
      endif
C          j = 1
C          do while (radii_p(ut).gt.radii_fine(np,j).and.j.le.nfine)
C           j = j + 1
C          enddo

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
     i  ( nlayers, nfinelayers, do_partials, extinction,
     i    n_partials, partials_idx, partials_fineidx,
     i    sunpaths, radii, ntraverse, alpha_all,
     i    sunpaths_p,    ntraverse_p,    alpha_p,
     i    sunpaths_fine, ntraverse_fine, alpha_fine,
     o    multipliers, lostrans, boa_attn,
     o    multipliers_p, lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/VLIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

C  control

      logical          do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_partlayers)
      integer          partials_fineidx(max_partlayers)

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

      integer          ntraverse_p(max_partlayers)
      double precision sunpaths_p (max_partlayers,maxlayers)
      double precision alpha_p    (max_partlayers)

C  Extinction 

      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)
      double precision multipliers_p   (max_partlayers)
      double precision lostrans_p      (max_partlayers)
      double precision boa_attn

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn      ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_partlayers )

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
     i  ( nlayers, nfinelayers, do_partials, extinction,
     i    n_partials, partials_idx, partials_fineidx,
     i    sunpaths, radii, ntraverse, alpha_all,
     i    sunpaths_p,   ntraverse_p,   alpha_p,
     i    sunpaths_fine, ntraverse_fine, alpha_fine,
     o    multipliers, lostrans,
     o    multipliers_p, lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/VLIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

C  control

      logical          do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_partlayers)
      integer          partials_fineidx(max_partlayers)

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

      integer          ntraverse_p(max_partlayers)
      double precision sunpaths_p (max_partlayers,maxlayers)
      double precision alpha_p    (max_partlayers)

C  Extinction 

      double precision extinction (maxlayers)

C  outputs
C  -------

      double precision multipliers (maxlayers)
      double precision lostrans    (maxlayers)
      double precision multipliers_p (max_partlayers)
      double precision lostrans_p    (max_partlayers)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn      ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_partlayers )

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

      SUBROUTINE VLIDORT_LAMBERTIAN_DBCORR (FLUXMULT)

C  Prepares Exact Direct Beam reflection for the Lambertian case

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Input arguments controlling type of surface

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Result goes in the Single scatter output file

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, O1
      DOUBLE PRECISION FINAL_SOURCE, TR, FACTOR
      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN

C  first stage
C  -----------

C  Initialize

      DO V = 1, N_GEOMETRIES
        DO O1 = 1, NSTOKES
          DO UTA = 1, N_USER_LEVELS
            STOKES_DB(UTA,V,O1)  = ZERO
          ENDDO
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
              V = VZA_OFFSETS(IB,UM) + UA

c  Beam attenuation

              IF ( DO_SSCORR_OUTGOING ) THEN
                X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                ATTN = BOA_ATTN(V)               
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                ELSE
                  X0_BOA = COS_SZANGLES(IB)
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
      O1 = 1
      FACTOR = FLUXMULT * LAMBERTIAN_ALBEDO
      DO V = 1, N_GEOMETRIES
        DB_CUMSOURCE(V,O1,NC) = ATTN_DB_SAVE(V) * FACTOR
      ENDDO

C  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

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
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS(N,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR = T_DELT_USERM(N,UM)
                ENDIF
                DB_CUMSOURCE(V,O1,NC) = TR * DB_CUMSOURCE(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

C  end layer loop

        ENDDO

C  Offgrid output : partial layer transmittance, then set result

        IF ( partlayers_OUTFLAG(UTA) ) THEN

          UT = partlayers_OUTINDEX(UTA)
          N  = partlayers_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS_UT(UT,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR = T_UTUP_USERM(UT,UM)
                ENDIF
                STOKES_DB(UTA,V,O1) = TR * DB_CUMSOURCE(V,O1,NC)
              ENDDO
            ENDDO
          ENDDO

C  Ongrid output : Set final cumulative source directly

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                FINAL_SOURCE = DB_CUMSOURCE(V,O1,NC)
                STOKES_DB(UTA,V,O1) = FINAL_SOURCE
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
     i    lospaths, sunpaths, ntraverse,
     o    multipliers, lostrans )

C  Does the optical depth integration over layers.

C  include files

      include '../includes/VLIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

      integer          nlayers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision lospaths(maxlayers)

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
     i    lospaths, sunpaths, ntraverse,
     o    multipliers, lostrans )

C  Does the optical depth integration over layers.

C  include files

      include '../includes/VLIDORT.PARS'

C  Geometry routine inputs
C  -----------------------

      integer          nlayers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision lospaths(maxlayers)

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

