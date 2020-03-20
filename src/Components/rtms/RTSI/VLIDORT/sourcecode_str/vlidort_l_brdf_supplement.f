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

C ##########################################################
C #                                                        #
C # Subroutines in this Module. The "L_BRDF SUPPLEMENT"    #
C #                                                        #
C #            VLIDORT_L_BRDF_CREATOR                      #
C #            VLIDORT_L_BRDF_FRTERMS                      #
C #                                                        #
C #    Basic BRDF setup, called by L_BRDF Creator:         #
C #                                                        #
C #            VLIDORT_BRDF_MAKER_PLUS                     #
C #                                                        #
C #  DB: Linearization w.r.t. atmospheric and surface      #
C #                                                        #
C #            VLIDORT_LAP_DBCORRECTION(master)            #
C #            VLIDORT_LAC_DBCORRECTION(master)            #
C #            VLIDORT_LS_DBCORRECTION(master)             #
C #                                                        #
C #    External calls in VLIDORT_L_FOURIER_BRDF:           #
C #    ----------------------------------------            #
C #                                                        #
C #  Linearization w.r.t. atmospheric parameter            #
C #                                                        #
C #            L_BVP_BRDF_SURFACE                          #
C #            L_BOA_BRDF_SOURCE                           #
C #                                                        #
C #  Linearization w.r.t. kernel factors                   #
C #                                                        #
C #            LS_BVP_BRDF_SURFACE                         #
C #            LS_BOA_BRDF_SOURCE                          #
C #                                                        #
C #  Linearization w.r.t. BRDF nonlinear parameters        #
C #                                                        #
C #            VLIDORT_KPARAMS_BRDFWF (master)             #
C #              KPARAMS_WFCOLSETUP                        #
C #              KPARAMS_BOA_SOURCE                        #
C #                                                        #
C ##########################################################


      SUBROUTINE VLIDORT_L_BRDF_CREATOR

C  Prepares (linearizations of) the bidirectional reflectance functions
C  necessary for LIDORT.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include files of VLIDORT linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Input arguments controlling type of surface
C     output arguments also go in here.

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  BRDF functions
C  --------------

C  ordinary BRDF without derivatives

      EXTERNAL       ROSSTHIN_VFUNCTION
      EXTERNAL       ROSSTHICK_VFUNCTION
      EXTERNAL       LISPARSE_VFUNCTION
      EXTERNAL       LIDENSE_VFUNCTION
      EXTERNAL       HAPKE_VFUNCTION
      EXTERNAL       ROUJEAN_VFUNCTION
      EXTERNAL       RAHMAN_VFUNCTION
      EXTERNAL       COXMUNK_VFUNCTION
      EXTERNAL       GISSCOXMUNK_VFUNCTION
      EXTERNAL       COXMUNK_VFUNCTION_DB
      EXTERNAL       GISSCOXMUNK_VFUNCTION_DB

C  Two new Land BRDFs, introduced 1 July 2008

      EXTERNAL       RHERMAN_VFUNCTION
      EXTERNAL       BREON_VFUNCTION

C  new for Version 2.4R, introduced 30 April 2009

      EXTERNAL       BREONVEG_VFUNCTION
      EXTERNAL       BREONSOIL_VFUNCTION

C  Hapke old uses exact DISORT code
C      EXTERNAL       HAPKE_VFUNCTION_OLD

C  BRDFs with derivatives

      EXTERNAL       LISPARSE_VFUNCTION_PLUS
      EXTERNAL       LIDENSE_VFUNCTION_PLUS
      EXTERNAL       HAPKE_VFUNCTION_PLUS
      EXTERNAL       RAHMAN_VFUNCTION_PLUS
      EXTERNAL       COXMUNK_VFUNCTION_PLUS
      EXTERNAL       GISSCOXMUNK_VFUNCTION_PLUS
      EXTERNAL       COXMUNK_VFUNCTION_PLUS_DB
      EXTERNAL       GISSCOXMUNK_VFUNCTION_PLUS_DB

C  Two new Land BRDFs, introduced 1 July 2008

      EXTERNAL       RHERMAN_VFUNCTION_PLUS
      EXTERNAL       BREON_VFUNCTION_PLUS

C  new for Version 2.4R, introduced 30 April 2009, 6 May 2009
C    2008 Veg/Soil functions based on Breon work 2008 as supplied to OCO
C    2009 function is final Kernel supplied by Breon, May 5 2009.

      EXTERNAL       BPDF2008VEG_VFUNCTION
      EXTERNAL       BPDF2008SOIL_VFUNCTION
      EXTERNAL       BPDF2009_VFUNCTION

C  local arguments
C  ---------------

      INTEGER          K, Q
      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL          LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

C  BRDF quadrature
C  ---------------

C  Save these quantities for efficient coding

      CALL BRDF_QUADRATURE

C  Fill BRDF arrays
C  ----------------

      DO K = 1, N_BRDF_KERNELS

C  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO Q = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(Q) = BRDF_PARAMETERS(K,Q)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO Q = 1, MAX_BRDF_PARAMETERS
            LOCAL_BRDF_DERIVS(Q) = DO_KERNEL_PARAMS_WFS(K,Q)
          ENDDO
        ENDIF

C  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROSSTHIN_VFUNCTION, ROSSTHIN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROSSTHICK_VFUNCTION, ROSSTHICK_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, LISPARSE_VFUNCTION_PLUS, LISPARSE_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, LISPARSE_VFUNCTION, LISPARSE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, LIDENSE_VFUNCTION_PLUS, LIDENSE_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, LIDENSE_VFUNCTION, LIDENSE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, HAPKE_VFUNCTION_PLUS, HAPKE_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, HAPKE_VFUNCTION, HAPKE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROUJEAN_VFUNCTION, ROUJEAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  Scalar-only original Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
         IF ( DO_COXMUNK_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I        ( K, COXMUNK_VFUNCTION_PLUS, COXMUNK_VFUNCTION_PLUS_DB,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I          LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I        ( K, COXMUNK_VFUNCTION, COXMUNK_VFUNCTION_DB,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I        ( K, COXMUNK_VFUNCTION_PLUS, COXMUNK_VFUNCTION_PLUS,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I          LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I        ( K, COXMUNK_VFUNCTION, COXMUNK_VFUNCTION,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ENDIF
        ENDIF

C  GISS Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX ) THEN
         IF ( DO_COXMUNK_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I   ( K, GISSCOXMUNK_VFUNCTION_PLUS, GISSCOXMUNK_VFUNCTION_PLUS_DB,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I          LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I        ( K, GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION_DB,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I    ( K, GISSCOXMUNK_VFUNCTION_PLUS, GISSCOXMUNK_VFUNCTION_PLUS,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I          LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I        ( K, GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION,
     I          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
         ENDIF
        ENDIF

C  Rahman kernel: (3 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, RAHMAN_VFUNCTION_PLUS, RAHMAN_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, RAHMAN_VFUNCTION, RAHMAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Rondeaux-Herman (vegetation) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RHERMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, RHERMAN_VFUNCTION_PLUS, RHERMAN_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, RHERMAN_VFUNCTION, RHERMAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  Breon et al (desert) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BREON_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS
     I          ( K, BREON_VFUNCTION_PLUS, BREON_VFUNCTION_PLUS,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I            LOCAL_BRDF_DERIVS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, BREON_VFUNCTION, BREON_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

C  BPDF 2008 Vegetation kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2008VEG_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, BPDF2008VEG_VFUNCTION, BPDF2008VEG_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  BPDF 2008 Soil kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2008SOIL_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, BPDF2008SOIL_VFUNCTION, BPDF2008SOIL_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

C  BPDF 2009 kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2009_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, BPDF2009_VFUNCTION, BPDF2009_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_BRDF_MAKER_PLUS
     I     ( BRDF_INDEX, BRDF_VFUNCTION_PLUS, BRDF_VFUNCTION_PLUS_DB,
     I       LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I       LOCAL_BRDF_DERIVS )

C  Prepares the bidirectional reflectance function derivatives

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of VLIDORT weighting function input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Include files of basic BRDF values

      INCLUDE '../includes/VLIDORT_BRDF.VARS'
      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Module arguments
C  ----------------

C  Index

      INTEGER          BRDF_INDEX

C  BRDF functions (external calls)

      EXTERNAL         BRDF_VFUNCTION_PLUS
      EXTERNAL         BRDF_VFUNCTION_PLUS_DB

C  Local number of parameters and local parameter array

      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL          LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )
            
C  local variables

      INTEGER          I, UI, J, K, Q, B, IB, O1
      DOUBLE PRECISION FUNC(MAXSTOKES_SQ)
      DOUBLE PRECISION DFUNC (MAXSTOKES_SQ, MAX_BRDF_PARAMETERS )
      DOUBLE PRECISION MUX, SZASURCOS(MAXBEAMS),SZASURSIN(MAXBEAMS)
      DOUBLE PRECISION PHIANG, COSPHI, SINPHI
      DOUBLE PRECISION SZAANG, COSSZA, SINSZA
      DOUBLE PRECISION VZAANG, COSVZA, SINVZA

C  shorthand

      B = BRDF_INDEX

C  Exact DB calculation
C  --------------------

      IF ( DO_DBCORRECTION ) THEN
        DO UI = 1, N_USER_STREAMS
          VZAANG = USER_VZANGLES_ADJUST(UI)
          COSVZA = DCOS(VZAANG*DEG_TO_RAD)
          SINVZA = DSIN(VZAANG*DEG_TO_RAD)        
          DO IB = 1, NBEAMS
            DO K = 1, N_USER_RELAZMS
              SZAANG = SZANGLES_ADJUST(UI,IB,K)
              COSSZA = DCOS(SZAANG*DEG_TO_RAD)
              SINSZA = DSIN(SZAANG*DEG_TO_RAD)
              PHIANG = USER_RELAZMS_ADJUST(UI,IB,K)
              COSPHI = DCOS(PHIANG*DEG_TO_RAD)
              SINPHI = DSIN(PHIANG*DEG_TO_RAD)
              CALL BRDF_VFUNCTION_PLUS_DB
     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          COSSZA, SINSZA, COSVZA, SINVZA,
     G          PHIANG, COSPHI, SINPHI,
     O          FUNC, DFUNC )
              DO O1 = 1, N_BRDF_STOKESSQ
               EXACTDB_BRDFUNC(O1,B,UI,K,IB) = FUNC(O1)
               DO Q = 1, LOCAL_BRDF_NPARS
                D_EXACTDB_BRDFUNC(O1,B,Q,UI,K,IB)  = DFUNC(O1,Q)
               ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF ( DO_SSFULL ) RETURN

C  Old code
c      IF ( DO_DBCORRECTION ) THEN
c        DO K = 1, N_USER_RELAZMS
c          PHIANG = USER_RELAZMS(K)
c          COSPHI = DCOS(PHIANG*DEG_TO_RAD)
c          SINPHI = DSIN(PHIANG*DEG_TO_RAD)
c          DO IB = 1, NBEAMS
c            DO UI = 1, N_USER_STREAMS
c              CALL BRDF_VFUNCTION_PLUS
c     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
c     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
c     C          LOCAL_BRDF_DERIVS,
c     G          SZASURCOS(IB), SZASURSIN(IB),
c     G          USER_STREAMS(UI), USER_SINES(UI),
c     G          PHIANG, COSPHI, SINPHI,
c     O          FUNC, DFUNC )
c              DO O1 = 1, N_BRDF_STOKESSQ
c               EXACTDB_BRDFUNC(O1,B,UI,K,IB) = FUNC(O1)
c               DO Q = 1, LOCAL_BRDF_NPARS
c                D_EXACTDB_BRDFUNC(O1,B,Q,UI,K,IB)  = DFUNC(O1,Q)
c               ENDDO
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDDO
c      ENDIF

C  Regular RT calculation
C  ======================

C  Usable solar beams
C  ------------------

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          MUX =  DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          SZASURCOS(IB) = COS_SZANGLES(IB)
          SZASURSIN(IB) = SIN_SZANGLES(IB)
        ENDDO
      ENDIF

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam

      DO IB = 1, NBEAMS
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF

C  Function call

          CALL BRDF_VFUNCTION_PLUS
     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          SZASURCOS(IB), SZASURSIN(IB), 
     G          QUAD_STREAMS(I), QUAD_SINES(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          FUNC, DFUNC )

C  copy to saved array

          DO O1 = 1, N_BRDF_STOKESSQ
           BRDFUNC_0(O1,B,I,IB,K) = FUNC(O1)
           DO Q = 1, LOCAL_BRDF_NPARS
            D_BRDFUNC_0(O1,B,Q,I,IB,K)  = DFUNC(O1,Q)
           ENDDO
          ENDDO

C  End loops

        ENDDO
       ENDDO
      ENDDO

C  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF

C  Function call

           CALL BRDF_VFUNCTION_PLUS
     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C          LOCAL_BRDF_DERIVS,
     G          QUAD_STREAMS(J), QUAD_SINES(J),
     G          QUAD_STREAMS(I), QUAD_SINES(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          FUNC, DFUNC )

C  copy to saved array

           DO O1 = 1, N_BRDF_STOKESSQ
            BRDFUNC(O1,B,I,J,K) = FUNC(O1)
            DO Q = 1, LOCAL_BRDF_NPARS
              D_BRDFUNC(O1,B,Q,I,J,K)  = DFUNC(O1,Q)
            ENDDO
           ENDDO

C  End loops

          ENDDO
        ENDDO
      ENDDO

C  Emissivity (optional) - BRDF quadrature input directions
C     PLACEHOLDER --> Linearized emissitivity not done here

C scalar code--------------------------------------------
c      IF ( DO_SURFACE_EMISSION ) THEN
c        NBRDF_HALF = NSTREAMS_BRDF / 2
c        DO I = 1, NSTREAMS
c          DO KE = 1, NBRDF_HALF
c            DO K = 1, NSTREAMS_BRDF
c              CALL BRDF_VFUNCTION_PLUS
c     C          ( MAX_BRDF_PARAMETERS,
c     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
c     C            LOCAL_BRDF_DERIVS,
c     G            CXE_BRDF(KE), SXE_BRDF(KE),
C     G            QUAD_STREAMS(I), QUAD_SINES(I),
c     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
c     O            FUNC, DFUNC )
c              EBRDFUNC(B,I,KE,K) = FUNC
c              DO Q = 1, LOCAL_BRDF_NPARS
c                D_EBRDFUNC(B,Q,I,KE,K)  = DFUNC(Q)
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDDO
c      ENDIF

C  User-streams outgoing directions
C  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

C  Incident Solar beam

        DO IB = 1, NBEAMS
         DO UI = 1, N_USER_STREAMS
          DO K = 1, NSTREAMS_BRDF

C  Function call

            CALL BRDF_VFUNCTION_PLUS
     C          ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C            LOCAL_BRDF_DERIVS,
     G            SZASURCOS(IB), SZASURSIN(IB),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            FUNC, DFUNC )

C  copy to saved array

            DO O1 = 1, N_BRDF_STOKESSQ
             USER_BRDFUNC_0(O1,B,UI,IB,K) = FUNC(O1)
             DO Q = 1, LOCAL_BRDF_NPARS
              D_USER_BRDFUNC_0(O1,B,Q,UI,IB,K)  = DFUNC(O1,Q)
             ENDDO
            ENDDO

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF

C  Function call

              CALL BRDF_VFUNCTION_PLUS
     C          ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     C            LOCAL_BRDF_DERIVS,
     G            QUAD_STREAMS(J), QUAD_SINES(J),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            FUNC, DFUNC )

C  copy to saved array

              DO O1 = 1, N_BRDF_STOKESSQ
               USER_BRDFUNC(O1,B,UI,J,K) = FUNC(O1)
               DO Q = 1, LOCAL_BRDF_NPARS
                D_USER_BRDFUNC(O1,B,Q,UI,J,K)  = DFUNC(O1,Q)
               ENDDO
              ENDDO

C  end loops

            ENDDO
          ENDDO
        ENDDO

C  Emissivity (optional) - BRDF quadrature input directions
C     PLACEHOLDER --> Linearized emissitivity not done here

C scalar code--------------------------------------------
c        IF ( DO_SURFACE_EMISSION ) THEN
c          DO UI = 1, N_USER_STREAMS
c            DO KE = 1, NBRDF_HALF
c              DO K = 1, NSTREAMS_BRDF
c                CALL BRDF_FUNCTION_PLUS
c     C            ( MAX_BRDF_PARAMETERS,
c     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
c     C              LOCAL_BRDF_DERIVS,
c     I              CXE_BRDF(KE), SXE_BRDF(KE),
c     I              USER_STREAMS(UI), USER_SINES(UI), 
c     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
c     O            FUNC, DFUNC )
c                USER_EBRDFUNC(B,UI,KE,K) = FUNC
c                DO Q = 1, LOCAL_BRDF_NPARS
c                  D_USER_EBRDFUNC(B,Q,UI,KE,K)  = DFUNC(Q)
c                ENDDO
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_BRDF_FRTERMS
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      SURFACE_FACTOR )

C  Prepares Fourier components of the bidirectional reflectance derivatives

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include file of VLIDORT weighting function input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Include files of basic BRDF values

      INCLUDE '../includes/VLIDORT_BRDF.VARS'
      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Module arguments
C  ----------------

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables

      INTEGER          I, UI, J, K, Q, B, IB, O1
      DOUBLE PRECISION SUM, HELP
      INTEGER          COSSIN_MASK(16)
      DATA COSSIN_MASK / 1,1,2,0,1,1,2,0,2,2,1,0,0,0,0,1 /

C  factor

      HELP = HALF * SURFACE_FACTOR

C  Quadrature outgoing directions
C  ------------------------------

C  Incident Solar beam (direct beam reflections)

      DO B = 1, N_BRDF_KERNELS
       IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
        DO IB = 1, NBEAMS
         IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          DO Q = 1, N_BRDF_PARAMETERS(B)
           IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
            DO I = 1, NSTREAMS
             DO O1 = 1, N_BRDF_STOKESSQ
              SUM = ZERO
              IF ( COSSIN_MASK(O1).EQ.1 ) THEN
               DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_BRDFUNC_0(O1,B,Q,I,IB,K)*BRDF_COSAZMFAC(K)
               ENDDO
              ELSE IF ( COSSIN_MASK(O1).EQ.2 ) THEN
               DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_BRDFUNC_0(O1,B,Q,I,IB,K)*BRDF_SINAZMFAC(K)
               ENDDO
              ENDIF
              D_BIREFLEC_0(O1,B,Q,I,IB) = SUM * HELP
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDIF
      ENDDO

C  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
       DO B = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
         DO Q = 1, N_BRDF_PARAMETERS(B)
          IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
           DO I = 1, NSTREAMS
            DO J = 1, NSTREAMS
             DO O1 = 1, N_BRDF_STOKESSQ
              SUM = ZERO
              IF ( COSSIN_MASK(O1).EQ.1 ) THEN
               DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_BRDFUNC(O1,B,Q,I,J,K)*BRDF_COSAZMFAC(K)
               ENDDO
              ELSE IF ( COSSIN_MASK(O1).EQ.2 ) THEN
               DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_BRDFUNC(O1,B,Q,I,J,K)*BRDF_SINAZMFAC(K)
               ENDDO
              ENDIF
              D_BIREFLEC(O1,B,Q,I,J) = SUM * HELP
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDIF

C  User-streams outgoing directions
C  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

C  Incident Solar beam (direct beam reflections)

       DO B = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
         DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
           DO Q = 1, N_BRDF_PARAMETERS(B)
            IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
             DO UI = 1, N_USER_STREAMS
              DO O1 = 1, N_BRDF_STOKESSQ
               SUM = ZERO
               IF ( COSSIN_MASK(O1).EQ.1 ) THEN
                DO K = 1, NSTREAMS_BRDF
                 SUM = SUM + 
     &             D_USER_BRDFUNC_0(O1,B,Q,UI,IB,K)*BRDF_COSAZMFAC(K)
                ENDDO
               ELSE IF ( COSSIN_MASK(O1).EQ.2 ) THEN
                DO K = 1, NSTREAMS_BRDF
                 SUM = SUM + 
     &             D_USER_BRDFUNC_0(O1,B,Q,UI,IB,K)*BRDF_SINAZMFAC(K)
                ENDDO
               ENDIF
               D_USER_BIREFLEC_0(O1,B,Q,UI,IB) = SUM * HELP
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO

C  incident quadrature directions (surface multiple reflections)

       IF ( DO_INCLUDE_SURFACE ) THEN
        DO B = 1, N_BRDF_KERNELS
         IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
          DO Q = 1, N_BRDF_PARAMETERS(B)
           IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
            DO UI = 1, N_USER_STREAMS
             DO J = 1, NSTREAMS
              DO O1 = 1, N_BRDF_STOKESSQ
               SUM = ZERO
               IF ( COSSIN_MASK(O1).EQ.1 ) THEN
                DO K = 1, NSTREAMS_BRDF
                 SUM = SUM + 
     &             D_USER_BRDFUNC(O1,B,Q,UI,J,K)*BRDF_COSAZMFAC(K)
                ENDDO
               ELSE IF ( COSSIN_MASK(O1).EQ.2 ) THEN
                DO K = 1, NSTREAMS_BRDF
                 SUM = SUM + 
     &             D_USER_BRDFUNC(O1,B,Q,UI,J,K)*BRDF_SINAZMFAC(K)
                ENDDO
               ENDIF
               D_USER_BIREFLEC(O1,B,Q,UI,J) = SUM * HELP
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDIF

      ENDIF

C  Emissivity
C  ----------

C  No linearization include here (PLACEHOLDER)
C    SCALAR CODE FOLLOWS (commented out)-------------------------------

c      IF ( DO_INCLUDE_SURFEMISS ) THEN
C  Start loop over kernels
c       DO B = 1, N_BRDF_KERNELS
c        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(B) ) THEN
c         DO Q = 1, N_BRDF_PARAMETERS(B)
c          IF ( DO_KERNEL_PARAMS_WFS(B,Q) ) THEN
C  Quadrature polar directions
c           DO I = 1, NSTREAMS
c            REFL = ZERO
c            DO KPHI= 1, NSTREAMS_BRDF
c             SUM = ZERO
c             DO K = 1, NBRDF_HALF
c              SUM = SUM + D_EBRDFUNC(B,Q,I,K,KPHI) * BAX_BRDF(K)
c             ENDDO
c             REFL = REFL + A_BRDF(KPHI) * SUM
c            ENDDO
c            D_EMISSIVITY(B,Q,I) = - REFL * BRDF_FACTORS(B)
c           ENDDO
C   user-defined polar directions
c           IF ( DO_USER_STREAMS ) THEN
c            DO UI = 1, N_USER_STREAMS
c             REFL = ZERO
c             DO KPHI= 1, NSTREAMS_BRDF
c              SUM = ZERO
c              DO K = 1, NBRDF_HALF
c               SUM = SUM + D_USER_EBRDFUNC(B,Q,UI,K,KPHI) * BAX_BRDF(K)
c              ENDDO
c              REFL = REFL + A_BRDF(KPHI) * SUM
c             ENDDO
c             D_USER_EMISSIVITY(B,Q,UI) = - REFL*BRDF_FACTORS(B)
c            ENDDO
c           ENDIF
C  parameter and kernel loop
c          ENDIF
c         ENDDO
c        ENDIF
c       ENDDO
C  end emissivity clause
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_LAP_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to PROFILE atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, KL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q, O1
      DOUBLE PRECISION L_A_EXACTDB_SOURCE, L_ATTN, SUM, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Start weighting function loop

      DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(K)

C  Initialize the output results

            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NSTOKES
                L_EXACTDB_SOURCE(Q,K,V,O1) = ZERO
                DO UTA = 1, N_USER_LEVELS
                  PROFILEWF_DB(Q,K,UTA,V,O1)      = ZERO
                ENDDO
              ENDDO
            ENDDO

C  Old code...................................................
C  For each solar beam
C    (a) Linearization factor L_ATTN for solar attenuation ATMOS_ATTN 
C    (b) for all geometries, inner um over BRDF kernels
C   R. Spurr, 24 May  2005. RT Solutions Inc. (Scalar)
C   R. Spurr, 19 July 2005. RT Solutions Inc. (Vector)
c            DO IB = 1, NBEAMS
c             IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
c              L_ATTN = L_DELTAU_SLANT(Q,NLAYERS,K,IB)   
c              DO UM = 1, N_USER_STREAMS
c               DO UA = 1, N_USER_RELAZMS
c                V = VZA_OFFSETS(IB,UM) + UA
c                DO O1 = 1, NSTOKES
c                 SUM = ZERO
c                 DO KL = 1, N_BRDF_KERNELS
c                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,O1,KL)* L_ATTN
c                  SUM = SUM + L_A_EXACTDB_SOURCE
c                 ENDDO
c                 L_EXACTDB_SOURCE(Q,K,V,O1) = SUM * FLUXMULT
c                ENDDO
c               ENDDO
c              ENDDO
c             ENDIF
c            ENDDO

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
                  L_ATTN = ZERO
                  IF (BOA_ATTN(V).GT.ZERO ) THEN
                    L_ATTN = L_BOA_ATTN(K,Q,V)/BOA_ATTN(V)
                  ENDIF
                ELSE
                  L_ATTN = - L_DELTAU_SLANT(Q,NLAYERS,K,IB)
                ENDIF

C  Loop over albedo kernels, assign reflection

                DO O1 = 1, NSTOKES
                 SUM = ZERO
                 DO KL = 1, N_BRDF_KERNELS
                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,O1,KL)* L_ATTN
                  SUM = SUM + L_A_EXACTDB_SOURCE
                 ENDDO
                 L_EXACTDB_SOURCE(Q,K,V,O1) = SUM * FLUXMULT
                ENDDO

C  Finish loops over geometries

               ENDDO
              ENDIF
             ENDDO
            ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

            NC =  0
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NSTOKES
                L_DB_CUMSOURCE(V,O1) = L_EXACTDB_SOURCE(Q,K,V,O1)
              ENDDO
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
C    loop over layers working upwards to level NUT

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N

C  If N  = K (varying layer), addition linearization of transmittance
C  If N ne K just transmitaance of the linearized source

                IF ( N.EQ.K) THEN
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = VZA_OFFSETS(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                      LTR = L_UP_LOSTRANS(N,Q,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                      LTR = L_T_DELT_USERM(N,UM,Q)
                    ENDIF
                    DO O1 = 1, NSTOKES
                     L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
     &                    + LTR * DB_CUMSOURCE(V,O1,NC-1) 
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                ELSE
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = VZA_OFFSETS(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                    ENDIF
                    DO O1 = 1, NSTOKES
                     L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF

C  end layer loop

              ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance
C  If N = K (varying layer), require additional linearization

              IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)

               IF ( N.EQ.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                     TR  = UP_LOSTRANS_UT(UT,V)
                     LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                   ELSE
                     TR  = T_UTUP_USERM(UT,UM)
                     LTR = L_T_UTUP_USERM(UT,UM,Q)
                   ENDIF
                   DO O1 = 1, NSTOKES
                    PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
     &                                     + LTR*DB_CUMSOURCE(V,O1,NC)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
               ELSE IF ( N.NE.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                     TR  = UP_LOSTRANS_UT(UT,V)
                   ELSE
                     TR  = T_UTUP_USERM(UT,UM)
                   ENDIF
                   DO O1 = 1, NSTOKES
                    PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
               ENDIF

C  Ongrid output : Set final cumulative source directly

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  DO O1 = 1, NSTOKES
                   PROFILEWF_DB(Q,K,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO

              ENDIF

C  Check for updating the recursion 

             IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
             NUT_PREV = NUT

C  end optical depth loop

            ENDDO

C  End weighting function loop

          ENDDO
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_LAC_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to COLUMN atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, KL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, K0, Q, O1
      DOUBLE PRECISION L_A_EXACTDB_SOURCE, L_ATTN, SUM, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Weighting function index

      K0 = 0

C  Start weighting function loop

      DO Q = 1, N_TOTALCOLUMN_WFS

C  Initialize the output results

        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            L_EXACTDBC_SOURCE(Q,V,O1) = ZERO
            DO UTA = 1, N_USER_LEVELS
              COLUMNWF_DB(Q,UTA,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

C  Old code...................................................
C  For each solar beam
C    (a) Linearization factor L_ATTN for solar attenuation ATMOS_ATTN 
C    (b) for all geometries, inner um over BRDF kernels
C   R. Spurr, 24 May  2005. RT Solutions Inc. (Scalar)
C   R. Spurr, 19 July 2005. RT Solutions Inc. (Vector)
c            DO IB = 1, NBEAMS
c             IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
c              L_ATTN = L_DELTAU_SLANT(Q,NLAYERS,K,IB)   
c              DO UM = 1, N_USER_STREAMS
c               DO UA = 1, N_USER_RELAZMS
c                V = VZA_OFFSETS(IB,UM) + UA
c                DO O1 = 1, NSTOKES
c                 SUM = ZERO
c                 DO KL = 1, N_BRDF_KERNELS
c                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,O1,KL)* L_ATTN
c                  SUM = SUM + L_A_EXACTDB_SOURCE
c                 ENDDO
c                 L_EXACTDB_SOURCE(Q,K,V,O1) = SUM * FLUXMULT
c                ENDDO
c               ENDDO
c              ENDDO
c             ENDIF
c            ENDDO

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
                  L_ATTN = ZERO
                  IF (BOA_ATTN(V).GT.ZERO ) THEN
                    L_ATTN = L_BOA_ATTN(K0,Q,V)/BOA_ATTN(V)
                  ENDIF
                ELSE
                  L_ATTN = ZERO
                  DO K = 1, NLAYERS
                    L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * 
     &                                   DELTAU_SLANT(NLAYERS,K,IB)
                  ENDDO   
                ENDIF

C  Loop over albedo kernels, assign reflection

                DO O1 = 1, NSTOKES
                  SUM = ZERO
                  DO KL = 1, N_BRDF_KERNELS
                    L_A_EXACTDB_SOURCE=A_EXACTDB_SOURCE(V,O1,KL)*L_ATTN
                    SUM = SUM + L_A_EXACTDB_SOURCE
                  ENDDO
                  L_EXACTDBC_SOURCE(Q,V,O1) = SUM * FLUXMULT
                ENDDO

C  Finish loops over geometries

              ENDDO
            ENDIF
          ENDDO
        ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            L_DB_CUMSOURCE(V,O1) = L_EXACTDBC_SOURCE(Q,V,O1)
          ENDDO
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
C    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

C  Linearization in very layer

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS(N,V)
                    LTR = L_UP_LOSTRANS(N,Q,V)
                  ELSE
                    TR  = T_DELT_USERM(N,UM)
                    LTR = L_T_DELT_USERM(N,UM,Q)
                  ENDIF
                  DO O1 = 1, NSTOKES
                   L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
     &                    + LTR * DB_CUMSOURCE(V,O1,NC-1) 
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

C  end layer loop

          ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                    LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                  ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                    LTR = L_T_UTUP_USERM(UT,UM,Q)
                  ENDIF
                  DO O1 = 1, NSTOKES
                    COLUMNWF_DB(Q,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
     &                                     + LTR*DB_CUMSOURCE(V,O1,NC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  DO O1 = 1, NSTOKES
                    COLUMNWF_DB(Q,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop

        ENDDO

C  End weighting function loop

      ENDDO

C  Finish

      RETURN
      END

c

      SUBROUTINE VLIDORT_LS_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to kernel amplitude and parameter variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, WF, O1, O2, OM
      INTEGER          PAR_INDEX, BRDF_INDEX
      DOUBLE PRECISION REFL_ATTN, SUM, L_KER, PAR, TR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  initialise output

      DO WF = 1, N_TOTALBRDF_WFS
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            LS_EXACTDB_SOURCE(WF,V,O1) = ZERO
            DO UTA = 1, N_USER_LEVELS
              SURFACEWF_DB(WF,UTA,V,O1)  = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  =========================================
C  1. Get the source function linearizations
C  =========================================

C  Initialise WF count

      WF = 0

C  Start loop over BRDF Kernels

      DO BRDF_INDEX = 1, N_BRDF_KERNELS

C  Amplitude weighting functions
C  -----------------------------

C  For each geometry, make Linearization w.r.t surface amplitude
C   R. Spurr, 24 May 2005. RT Solutions Inc.

        IF ( DO_KERNEL_FACTOR_WFS(BRDF_INDEX) ) THEN
          WF = WF + 1
          DO IB = 1, NBEAMS
           IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              DO O1 = 1, NSTOKES
               SUM= A_EXACTDB_SOURCE(V,O1,BRDF_INDEX)
               LS_EXACTDB_SOURCE(WF,V,O1) = SUM * FLUXMULT
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
        ENDIF

C  Parameter weighting functions
C  -----------------------------

C Old code..............................................................
C  For each geometry, make Linearization w.r.t nonlinear parameter
C   R. Spurr, 24 May 2005. RT Solutions Inc.
c        DO PAR_INDEX = 1, N_BRDF_PARAMETERS(BRDF_INDEX)
c         IF (DO_KERNEL_PARAMS_WFS(BRDF_INDEX,PAR_INDEX)) THEN
c          WF = WF + 1
c          PAR = BRDF_PARAMETERS(BRDF_INDEX,PAR_INDEX)
c          DO IB = 1, NBEAMS
c           IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
c            REFL_ATTN =  BRDF_FACTORS(BRDF_INDEX) * ATMOS_ATTN(IB)
c            DO UM = 1, N_USER_STREAMS
c             DO UA = 1, N_USER_RELAZMS
c              V = VZA_OFFSETS(IB,UM) + UA
c              DO O1 = 1, NSTOKES
c               SUM = ZERO
c               DO O2 = 1, NSTOKES
c                OM = MUELLER_INDEX(O1,O2)
c                L_KER = 
c     &           D_EXACTDB_BRDFUNC(O1,BRDF_INDEX,PAR_INDEX,UM,UA,IB)
c                SUM = SUM + FLUX_FACTOR * L_KER * PAR * FLUXVEC(O2)
c               ENDDO
c               LS_EXACTDB_SOURCE(WF,V,O1) = REFL_ATTN * SUM
c              ENDDO 
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
c         ENDIF
c        ENDDO

C  New code

        DO PAR_INDEX = 1, N_BRDF_PARAMETERS(BRDF_INDEX)
         IF (DO_KERNEL_PARAMS_WFS(BRDF_INDEX,PAR_INDEX)) THEN
          WF = WF + 1
          PAR = BRDF_PARAMETERS(BRDF_INDEX,PAR_INDEX)
          DO UM = 1, N_USER_STREAMS
           DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              REFL_ATTN =  BRDF_FACTORS(BRDF_INDEX) * ATTN_DB_SAVE(V)
              DO O1 = 1, NSTOKES
               SUM = ZERO
               DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                L_KER = 
     &           D_EXACTDB_BRDFUNC(O1,BRDF_INDEX,PAR_INDEX,UM,UA,IB)
                SUM = SUM + FLUX_FACTOR * L_KER * PAR * FLUXVEC(O2)
               ENDDO
               LS_EXACTDB_SOURCE(WF,V,O1) = REFL_ATTN * SUM
              ENDDO 
             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO

C  End kernel loop

      ENDDO

C  ====================================================
C  2. Transmittance of source term: upwelling recursion
C  ====================================================

      DO WF = 1, N_TOTALBRDF_WFS

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
         DO O1 = 1, NSTOKES
          LS_DB_CUMSOURCE(V,O1) = LS_EXACTDB_SOURCE(WF,V,O1)*FLUXMULT
         ENDDO
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
C    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS(N,V)
              ELSE
                TR  = T_DELT_USERM(N,UM)
              ENDIF
              DO O1 = 1, NSTOKES
                LS_DB_CUMSOURCE(V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO

C  Offgrid output
C  -------------- 

C  Require partial layer transmittance

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
           UT = PARTLAYERS_OUTINDEX(UTA)
           N  = PARTLAYERS_LAYERIDX(UT)
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,UM)
              ENDIF
              DO O1 = 1, NSTOKES
               SURFACEWF_DB(WF,UTA,V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              DO O1 = 1, NSTOKES
               SURFACEWF_DB(WF,UTA,V,O1) = LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop

        ENDDO

C  End Loop over all Surface  weighting functions

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BVP_BRDF_SURFACE
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BCL4,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Number of weighting functions

      INTEGER          N_LAYER_WFS

C  Fourier component and beam index

      INTEGER          FOURIER_COMPONENT

C  Flag for type of boundary condition

      LOGICAL          MODIFIED_BCL4

C  overall surface flag and surface factor

      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR

C  Output arguments
C  ----------------

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

C  Local variables
C  ---------------

      DOUBLE PRECISION PV_W
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION HV_P
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION HV_M
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      INTEGER          I, J, O1, O2, OM, KL, Q, N
      INTEGER          K, KO1, K0, K1, K2

      DOUBLE PRECISION H1R, H1I, REFL_P, REFL_B, REFL_M
      DOUBLE PRECISION FACTOR_KERNEL, H1, H2, HP, HM, BIREF

      DOUBLE PRECISION H1_CR,   H2_CR,   H1_CI,   H2_CI
      DOUBLE PRECISION H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI

C  Initial section
C  ---------------

C  Always zero the result to start

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES
        DO Q = 1, N_LAYER_WFS
          R2_L_BEAM(I,O1,Q)    = ZERO
          DO K = 1, NSTKS_NSTRMS
            R2_L_HOMP(I,O1,K,Q) = ZERO
            R2_L_HOMM(I,O1,K,Q) = ZERO
          ENDDO
        ENDDO
       ENDDO
      ENDDO

C  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Set up Auxiliary arrays
C  -----------------------

      N = NLAYERS

C  Need Beam arrays regardless
C   Always something from layers above the PBL and also from PBL

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QUAD_STRMWTS(J)
          ENDDO
        ENDDO
      ENDDO

C   Modified boundary condition
C     This applies when PBL is varying layer.
C       Require homogeneous solutions linearizations

      IF ( MODIFIED_BCL4 ) THEN

C  start loops

        DO J = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS

C  real homogeneous solution contributions

           DO K = 1, K_REAL(N)
            H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) +
     &             SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
            H2 = L_SOLB_XNEG(J,O1,K,N,Q)
            HV_P(J,O1,K,Q) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K,Q) = QUAD_STRMWTS(J)*H2
           ENDDO

C  Complex homogeneous solution contributions

           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N)
     &          - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N)
     &          +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q)
     &          -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
            H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N)
     &          + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N)
     &          +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q)
     &          +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
            HV_P(J,O1,K1,Q) = QUAD_STRMWTS(J)* H1R
            HV_P(J,O1,K2,Q) = QUAD_STRMWTS(J)* H1I
            HV_M(J,O1,K1,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K1,N,Q)
            HV_M(J,O1,K2,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K2,N,Q)
           ENDDO

C  End loops

          ENDDO
         ENDDO
        ENDDO

C  End modified BCL4 condition

      ENDIF

C  Integrated Downward reflection (Calculation)
C  --------------------------------------------

C  Start loop over BRDF kernels

      DO KL = 1, N_BRDF_KERNELS

C  Kernel amplitude

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(KL)

C  Lambertian Reflection
C  =====================

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
         IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, N_LAYER_WFS

C  Particular solution (only for the first Stokes component)

           O1 = 1
           REFL_B = ZERO
           DO J = 1, NSTREAMS
             REFL_B = REFL_B + PV_W(J,O1,Q)
           ENDDO
           REFL_B = REFL_B * FACTOR_KERNEL
           DO I = 1, NSTREAMS
             R2_L_BEAM(I,O1,Q) = R2_L_BEAM(I,O1,Q) + REFL_B
           ENDDO

C  Homogeneous solutions for the modified condition

           IF ( MODIFIED_BCL4 ) THEN

C  Homogeneous real solutions

             DO K = 1, K_REAL(NLAYERS)
              REFL_P = ZERO
              REFL_M = ZERO
              DO J = 1, NSTREAMS
               REFL_P = REFL_P + HV_P(J,O1,K,Q)
               REFL_M = REFL_M + HV_M(J,O1,K,Q)
              ENDDO
              REFL_P = REFL_P * FACTOR_KERNEL
              REFL_M = REFL_M * FACTOR_KERNEL
              DO I = 1, NSTREAMS
               R2_L_HOMP(I,O1,K,Q) = R2_L_HOMP(I,O1,K,Q) + REFL_P
               R2_L_HOMM(I,O1,K,Q) = R2_L_HOMM(I,O1,K,Q) + REFL_M
              ENDDO
             ENDDO

C  Homogeneous complex solutions

             KO1 = K_REAL(NLAYERS) + 1
             DO K = 1, K_COMPLEX(NLAYERS)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              H1_CR = ZERO
              H1_CI = ZERO
              H2_CR = ZERO
              H2_CI = ZERO
              DO J = 1, NSTREAMS
               H1_CR = H1_CR + HV_P(J,O1,K1,Q)
               H1_CI = H1_CI + HV_P(J,O1,K2,Q)
               H2_CR = H2_CR + HV_M(J,O1,K1,Q)
               H2_CI = H2_CI + HV_M(J,O1,K2,Q)
              ENDDO
              H1_CR = H1_CR * FACTOR_KERNEL
              H1_CI = H1_CI * FACTOR_KERNEL
              H2_CR = H2_CR * FACTOR_KERNEL
              H2_CI = H2_CI * FACTOR_KERNEL
              DO I = 1, NSTREAMS
               R2_L_HOMP(I,O1,K1,Q) = R2_L_HOMP(I,O1,K1,Q) + H1_CR
               R2_L_HOMP(I,O1,K2,Q) = R2_L_HOMP(I,O1,K2,Q) + H1_CI
               R2_L_HOMM(I,O1,K1,Q) = R2_L_HOMM(I,O1,K1,Q) + H2_CR
               R2_L_HOMM(I,O1,K2,Q) = R2_L_HOMM(I,O1,K2,Q) + H2_CI
              ENDDO
             ENDDO

C  end parameter loop

           ENDIF
          ENDDO
         ENDIF

C  Bidirectional reflecting kernel
C  ===============================

        ELSE

C  start loops

          DO I = 1, NSTREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, N_LAYER_WFS

C  Particular solution

             REFL_B = ZERO
             DO J = 1, NSTREAMS
              H1 = ZERO
              DO O2 = 1, NSTOKES
               OM = MUELLER_INDEX(O1,O2)
               H1 = H1 + PV_W(J,O2,Q)*BIREFLEC(OM,KL,J,I)
              ENDDO
              REFL_B = REFL_B + H1
             ENDDO
             REFL_B = REFL_B * FACTOR_KERNEL
             R2_L_BEAM(I,O1,Q) = R2_L_BEAM(I,O1,Q) + REFL_B

C  Homogeneous solutions for the modified condition

             IF ( MODIFIED_BCL4 ) THEN

C  Homogeneous real solutions

              DO K = 1, K_REAL(NLAYERS)
               REFL_P = ZERO
               REFL_M = ZERO
               DO J = 1, NSTREAMS
                HP = ZERO
                HM = ZERO
                DO O2 = 1, NSTOKES
                 OM = MUELLER_INDEX(O1,O2)
                 HP = HP + HV_P(J,O2,K,Q)*BIREFLEC(OM,KL,J,I)
                 HM = HP + HV_M(J,O2,K,Q)*BIREFLEC(OM,KL,J,I)
                ENDDO
                REFL_P = REFL_P + HP
                REFL_M = REFL_M + HM
               ENDDO
               REFL_P = REFL_P * FACTOR_KERNEL
               REFL_M = REFL_M * FACTOR_KERNEL
               R2_L_HOMP(I,O1,K,Q) = R2_L_HOMP(I,O1,K,Q) + REFL_P
               R2_L_HOMM(I,O1,K,Q) = R2_L_HOMM(I,O1,K,Q) + REFL_M
              ENDDO

C  homogeneous complex solutions

              KO1 = K_REAL(NLAYERS) + 1
              DO K = 1, K_COMPLEX(NLAYERS)
               K0 = 2*K - 2
               K1 = KO1 + K0
               K2 = K1  + 1
               H1_CR = ZERO
               H2_CR = ZERO
               H1_CI = ZERO
               H2_CI = ZERO
               DO J = 1, NSTREAMS
                H1_S_CR = ZERO
                H2_S_CR = ZERO
                H1_S_CI = ZERO
                H2_S_CI = ZERO
                DO O2 = 1, NSTOKES
                 OM    = MUELLER_INDEX(O1,O2)
                 BIREF = BIREFLEC(OM,KL,J,I)
                 H1_S_CR = H1_S_CR + HV_P(J,O2,K1,Q)*BIREF
                 H2_S_CR = H2_S_CR + HV_M(J,O2,K1,Q)*BIREF
                 H1_S_CI = H1_S_CI + HV_P(J,O2,K2,Q)*BIREF
                 H2_S_CI = H2_S_CI + HV_M(J,O2,K2,Q)*BIREF
                ENDDO
                H1_CR = H1_CR + H1_S_CR
                H2_CR = H2_CR + H2_S_CR
                H1_CI = H1_CI + H1_S_CI
                H2_CI = H2_CI + H2_S_CI
               ENDDO
               H1_CR = H1_CR * FACTOR_KERNEL
               H2_CR = H2_CR * FACTOR_KERNEL
               H1_CI = H1_CI * FACTOR_KERNEL
               H2_CI = H2_CI * FACTOR_KERNEL
               R2_L_HOMP(I,O1,K1,Q) = R2_L_HOMP(I,O1,K1,Q) + H1_CR
               R2_L_HOMM(I,O1,K1,Q) = R2_L_HOMM(I,O1,K1,Q) + H2_CR
               R2_L_HOMP(I,O1,K2,Q) = R2_L_HOMP(I,O1,K2,Q) + H1_CI
               R2_L_HOMM(I,O1,K2,Q) = R2_L_HOMM(I,O1,K2,Q) + H2_CI
              ENDDO

C  End modified condition

             ENDIF

C  end parameter and stream and Stokes loops

            ENDDO
           ENDDO
          ENDDO

C  End Kernel clause

        ENDIF

C  End kernels loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BOA_BRDF_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_DIRECTBEAM,
     I      SURFACE_FACTOR,
     I      FOURIER_COMPONENT, IBEAM,
     I      NV, NV_PARAMETERS,
     I      L_BOA_MSSOURCE, L_BOA_DBSOURCE)

C  Linearized Bottom of the atmosphere source term

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  albedo and Fourier/beam indices

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM

C  linearization control

      INTEGER          NV, NV_PARAMETERS

C  output arguments
C  ----------------

      DOUBLE PRECISION L_BOA_MSSOURCE
     &          (MAX_USER_STREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION L_BOA_DBSOURCE
     &          (MAX_USER_STREAMS,MAXSTOKES,MAX_ATMOSWFS)

C  local variables
C  ---------------

      INTEGER          M, N, J, I, UM
      INTEGER          Q, KL, IB, O1, O2, OM
      INTEGER          K, KO1, K0, K1, K2

      DOUBLE PRECISION L_DOWN (MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

      DOUBLE PRECISION REFLEC, S_REFLEC, L_BEAM, FAC, KMULT
      DOUBLE PRECISION SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION LXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION LXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION LXR2, NXR2, LLXR2

C  Starting section
C  ----------------

C  Fourier number, layer number

      M   = FOURIER_COMPONENT
      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  initialise linearized BOA source functions

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, NV_PARAMETERS
          DO O1 = 1, NSTOKES
            L_BOA_MSSOURCE(UM,O1,Q) = ZERO
            L_BOA_DBSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  reflectance from surface
C  ------------------------

C  reflectance integrand  a(j).x(j).L(downwelling)(-j)

      IF ( DO_INCLUDE_SURFACE ) THEN

C  Two cases:
C  (a) If  NV = N, this is also the layer that is varying --> Extras!
C  (b) If  N > NV with variations in layer NV above N

        IF ( NV .EQ. N ) THEN

C  stream and stokes loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N)
     &                       +  LXR   * L_T_DELT_EIGEN(K,N,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N)
     &                  - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N)
     &                  + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
     &                   - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
     &                   + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q)
     &                   - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N)
     &                   - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q)
     &                             - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR 

               ENDDO

C  real part and add particular solution

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I,O1,N,Q)
               L_DOWN(I,O1,Q) = QUAD_STRMWTS(I) * ( SPAR + SHOM )

C  End loops over Q, O1 and I

              ENDDO
            ENDDO
          ENDDO

C  otherwise the varying layer is above the boundary layer

        ELSE IF (NV.LT.N) THEN

C  stream and stokes loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_DELT_EIGEN(K,N)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N)
     &                   - NXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

C  real part homogeneous solution

               SHOM = SHOM_R + SHOM_CR

C  Add particular solution

               SPAR = L_WLOWER(I,O1,N,Q)

C  Final result

               L_DOWN(I,O1,Q) = QUAD_STRMWTS(I) * ( SPAR + SHOM )

C  End loops over Q, O1 and I

              ENDDO
            ENDDO
          ENDDO

C  End cases

        ENDIF

C  reflected multiple scatter intensity at user defined-angles
C  -----------------------------------------------------------

        DO KL = 1, N_BRDF_KERNELS

          KMULT = SURFACE_FACTOR * BRDF_FACTORS(KL)

C  ###### Lambertian reflectance (same for all user-streams)

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              O1 = 1
              DO Q = 1, NV_PARAMETERS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + L_DOWN(J,O1,Q)
                ENDDO
                REFLEC = KMULT * REFLEC
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q)
     &                             + REFLEC
                ENDDO
              ENDDO
            ENDIF

C  .. integrate with BRDF reflectance function at user angles 

          ELSE

            DO Q = 1, NV_PARAMETERS
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  REFLEC = ZERO
                  DO J = 1, NSTREAMS
                    S_REFLEC = ZERO
                    DO O2 = 1, NSTOKES
                      OM = MUELLER_INDEX(O1,O2)
                      S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * 
     &                     USER_BIREFLEC(OM,KL,UM,J)
                    ENDDO
                    REFLEC = REFLEC + S_REFLEC
                  ENDDO
                  REFLEC = KMULT * REFLEC
                  L_BOA_MSSOURCE(UM,O1,Q) = 
     &                   L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  end loop over kernels

        ENDDO
        
C  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              FAC = - USER_DIRECT_BEAM(UM,IB,O1)*DELTAU_SLANT(N,NV,IB)
              DO Q = 1, NV_PARAMETERS
                L_BEAM = L_DELTAU_VERT(Q,NV) * FAC
                L_BOA_DBSOURCE(UM,O1,Q) = L_BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  End inclusion of surface terms

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_BVP_BRDF_SURFACE
     I   ( DO_REGULAR_BVP,DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I     IBEAM_INDEX, BRDF_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  inputs
C  ------

C  .. regular BVP flag. Shoudl be true, not used here.

      LOGICAL          DO_REGULAR_BVP

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission inputs

      LOGICAL          DO_INCLUDE_SURFEMISS

C  Beam index

      INTEGER          IBEAM_INDEX
  
C  BRDF index (indicates which albedo kernel)

      INTEGER          BRDF_INDEX

C  local variables
C  ---------------

      INTEGER          N, I, C0, CM, N1, IB, B, IROW, IR, O1
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I

C  initialise (Should not be necessary!)

      DO I = 1, NTOTAL
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for first layer upper (TOA)

      N = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COL2_WFALB(IROW,1) = ZERO
        ENDDO
      ENDDO

C  boundary conditions not changed for all intermediate levels

      DO N = 2, NLAYERS - 1
        N1 = N - 1
        C0 = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,1) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      N  = NLAYERS
      C0 = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      B  = BRDF_INDEX
      IB = IBEAM_INDEX

C  Diffuse scatter contributions

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM   = C0 + IROW

C  Beam contribution for this Kernel

          L_BEAM = A_DIFFUSE_BEAM(I,O1,B)

C  Real homogeneous solutions for this Kernel

          L_HOM_R = ZERO
          DO K = 1, K_REAL(N)
            L_HOM_R = L_HOM_R
     &       + LCON(K,N) * A_DIFFUSE_HOMP(I,O1,K,B) * T_DELT_EIGEN(K,N) 
     &       + MCON(K,N) * A_DIFFUSE_HOMM(I,O1,K,B)
          ENDDO

C  Complex homogeneous solutions for this Kernel

          L_HOM_CR = ZERO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1 + 1
            H1R =   A_DIFFUSE_HOMP(I,O1,K1,B) * T_DELT_EIGEN(K1,N)
     &            - A_DIFFUSE_HOMP(I,O1,K2,B) * T_DELT_EIGEN(K2,N)
            H1I =   A_DIFFUSE_HOMP(I,O1,K1,B) * T_DELT_EIGEN(K2,N)
     &            + A_DIFFUSE_HOMP(I,O1,K2,B) * T_DELT_EIGEN(K1,N)
            L_HOM_CR = L_HOM_CR
     &                  + LCON(K1,N) * H1R  - LCON(K2,N) * H1I 
     &                  + MCON(K1,N) * A_DIFFUSE_HOMM(I,O1,K1,B)
     &                  - MCON(K2,N) * A_DIFFUSE_HOMM(I,O1,K2,B)
          ENDDO

C  Final contribution

          COL2_WFALB(CM,1) = L_BEAM + L_HOM_R + L_HOM_CR

C  End Streams and Stokes loops

        ENDDO
      ENDDO

C  Add direct beam variation of albedo

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,1) = 
     &         COL2_WFALB(CM,1) + A_DIRECT_BEAM(I,IB,O1,B)
          ENDDO
        ENDDO
      ENDIF

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX,N, COL2_WFALB(N,1)
c        ENDDO
c      ENDIF
           
C  Finish

      RETURN
      END

C

      SUBROUTINE LS_BOA_BRDF_SOURCE
     I        ( DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, 
     O          L_BOA_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM
C  kernel index

      INTEGER          BRDF_INDEX

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  output

      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)

C  Local variables
C  ---------------

      INTEGER          UM, I, J, KL, N, IB, O1, O2, OM
      INTEGER          K, KO1, K0, K1, K2      
      DOUBLE PRECISION INTEGRAND(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION SUM_R, SUM_CR, REFLEC, S_REFLEC
      DOUBLE PRECISION FACTOR_BRDF, H1, H2
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1

C  Initialise
C  ----------

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  initialise Derivative of BOA source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          L_BOA_SOURCE(UM,O1) = ZERO
        ENDDO
      ENDDO

C  Return if only quadrature

      IF ( .NOT. DO_USER_STREAMS ) RETURN

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

C  start loops

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES

C  Real homogeneous solutions

        SUM_R = ZERO
        DO K = 1, K_REAL(N)
         NXR = NCON_ALB(K,N) * SOLA_XPOS(I,O1,K,N)
         PXR = PCON_ALB(K,N) * SOLB_XNEG(I,O1,K,N)
         SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
        ENDDO

C  Complex solutions

        SUM_CR = ZERO
        DO K = 1, K_COMPLEX(N)
         K0 = 2 * K - 2
         K1 = KO1 + K0
         K2 = K1  + 1
         NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K1,N)
     &          - NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K2,N)
         NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K2,N)
     &          + NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K1,N)
         PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K1,N)
     &          - PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K2,N)
         H1 =  NXR1 *  * T_DELT_EIGEN(K1,N)
     &        -NXR2 *  * T_DELT_EIGEN(K2,N)
         H2 =  PXR1
         SUM_CR = SUM_CR + H1 + H2
        ENDDO

C  Final result

        INTEGRAND(I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

C  end loops

       ENDDO
      ENDDO

C  integrated BRDF reflectance term (loop over all kernels)
C  --------------------------------------------------------

      DO KL = 1, N_BRDF_KERNELS

        FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)

C  .. Either, integrate Lambertian case, same for all user-streams

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) .AND.
     &                 FOURIER_COMPONENT.EQ.0 ) THEN
          O1 = 1
          REFLEC = ZERO
          DO J = 1, NSTREAMS
             REFLEC = REFLEC + INTEGRAND(J,O1)
          ENDDO
          REFLEC = FACTOR_BRDF * REFLEC
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
          ENDDO

C  .. Or, integrate with BRDF kernel function at user angles

        ELSE
  
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + INTEGRAND(J,O2) * 
     &                           USER_BIREFLEC(OM,KL,UM,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              REFLEC = FACTOR_BRDF * REFLEC
              L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
            ENDDO
          ENDDO

        ENDIF

C  End loop over all kernels

      ENDDO

C  Contributions due to direct variation of kernel factor
C  ------------------------------------------------------

C  Add linearization term due to kernel factor variation 

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1)
     &              + A_USER_DIFFUSE_SOURCE  ( UM, O1, BRDF_INDEX )
     &              + A_USER_DIRECT_SOURCE ( UM, IB, O1, BRDF_INDEX )
        ENDDO
      ENDDO

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)
C      IF ( DO_INCLUDE_SURFEMISS ) THEN
C        DO UM = LOCAL_UM_START, N_USER_STREAMS
C          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) - 
C     &           FP_SURFBB * ALBEDO_USER_EMISSIVITY(UM,BRDF_INDEX)
C        ENDDO
C      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_KPARAMS_BRDFWF
     I       (  DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          FOURIER_COMPONENT, IBEAM,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER,
     I          BRDF_INDEX, PAR_INDEX, WF_INDEX,
     O          STATUS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of solution variables

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  local control flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Fourier surface factor, Fourier number, Beam index, Flux multiplier

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Weighting function index number

      INTEGER          WF_INDEX

C  BRDF index (indicates which BRDF kernel)

      INTEGER          BRDF_INDEX

C  Parameter index

      INTEGER          PAR_INDEX

C  Output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Other local variables

      INTEGER          N, C0, IROW, IROW1, IROW_S, IROW1_S
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)

C  Initialise status

       STATUS = VLIDORT_SUCCESS

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

      CALL KPARAMS_WFCOLSETUP
     I       (  DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR,
     I          IBEAM, BRDF_INDEX, PAR_INDEX )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (multilayer) in VLIDORT_KPARAMS_BRDFWF'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO N = 1, NLAYERS
          C0 = (N-1)*NSTKS_NSTRMS_2
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_ALB(K,N) = COL2_WFALB(C0+IROW,1)
            PCON_ALB(K,N) = COL2_WFALB(C0+IROW1,1)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_ALB(K1,N) = COL2_WFALB(C0+IROW,   1)
            NCON_ALB(K2,N) = COL2_WFALB(C0+IROW_S, 1)
            PCON_ALB(K1,N) = COL2_WFALB(C0+IROW1,  1)
            PCON_ALB(K2,N) = COL2_WFALB(C0+IROW1_S,1)
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2_WFALB, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (one layer) in VLIDORT_KPARAMS_BRDFWF'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        N = 1
        DO K = 1, K_REAL(N)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          NCON_ALB(K,N) = SCOL2_WFALB(IROW,1)
          PCON_ALB(K,N) = SCOL2_WFALB(IROW1,1)
        ENDDO
        KO1 = K_REAL(N) + 1
        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          NCON_ALB(K1,N) = SCOL2_WFALB(IROW,   1)
          NCON_ALB(K2,N) = SCOL2_WFALB(IROW_S, 1)
          PCON_ALB(K1,N) = SCOL2_WFALB(IROW1,  1)
          PCON_ALB(K2,N) = SCOL2_WFALB(IROW1_S,1)
        ENDDO

C  end clause
 
      ENDIF

C  debug------------------------------------------
        if ( do_debug_write) then
         if (fourier_component.eq.0 ) then
         DO N = 1, NLAYERS
          DO K = 1, K_REAL(N)
           write(86,'(4I3,1p4e20.10)')FOURIER_COMPONENT,IBEAM,K,N,
     &                LCON(K,N), MCON(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N)
          ENDDO
         ENDDO
        ENDIF
        ENDIF

C  Get the weighting functions
C  ---------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE) for these weighting functions

        CALL KPARAMS_BOA_SOURCE
     I        ( DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, PAR_INDEX,
     O          L_BOA_SOURCE )

C  Upwelling Surface weighting functions (w.r.t Kernel amplitudes)

        CALL UPUSER_SURFACEWF
     I      ( DO_INCLUDE_SURFEMISS,
     I        FLUX_MULTIPLIER,
     I        IBEAM, WF_INDEX,
     I        L_BOA_SOURCE )

      ENDIF

C  Downwelling Surface weighting functions (w.r.t Kernel amplitudes)

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF
     I       ( FLUX_MULTIPLIER, IBEAM,
     I         WF_INDEX )
      ENDIF

C  mean value output

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL VLIDORT_LS_INTEGRATED_OUTPUT
     I  ( DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER,
     I    IBEAM, WF_INDEX )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE KPARAMS_WFCOLSETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I      SURFACE_FACTOR, IBEAM_INDEX, BRDF_INDEX, PAR_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of linearized BRDF variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  inputs
C  ------

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission inputs

      LOGICAL          DO_INCLUDE_SURFEMISS

C  Beam index

      INTEGER          IBEAM_INDEX
 
C  .. BRDF index (indicates which albedo kernel)

      INTEGER          BRDF_INDEX

C  -- parameter index

      INTEGER          PAR_INDEX

C  factor

      DOUBLE PRECISION SURFACE_FACTOR

C  local variables
C  ---------------

      INTEGER          N,  N1, IB, B, Q
      INTEGER          K, KO1, K0, K1, K2
      INTEGER          I, J, O1, O2, OM, IR, IROW, C0, CM
      DOUBLE PRECISION REFL_B, L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I
      DOUBLE PRECISION PARVALUE, FACTOR_BRDF, REFL_ATTN, AWF_DIRECT
      DOUBLE PRECISION H_1,   H_2,   H_1_S,   H_2_S, H1, H2, DBR
      DOUBLE PRECISION H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

C  Help arrays

      DOUBLE PRECISION PV_W  ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION HV_P  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION HV_M  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

C  initialise. Vitally necessary
C    We found that when commented out, the whole thing didn't work.
C    R. Spurr and V. Natraj, 20 january 2006

      OM = 0
      DO I = 1, NTOTAL
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for first layer upper (TOA)

      N = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COL2_WFALB(IROW,1) = ZERO
        ENDDO
      ENDDO

C  boundary conditions not changed for all intermediate levels

      DO N = 2, NLAYERS - 1
        N1 = N - 1
        C0 = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,1) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      N  = NLAYERS
      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      Q  = PAR_INDEX
      B  = BRDF_INDEX
      IB = IBEAM_INDEX

      PARVALUE    = BRDF_PARAMETERS(B,Q)
      FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(B) * PARVALUE

C  Save some quantities
C  --------------------

C  ( This is a repetition of earlier code and could be stored)

C  start loops

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES

C  Beam

          PV_W(J,O1) = WLOWER(J,O1,N) * QUAD_STRMWTS(J)

C  real homogeneous solution contributions

          DO K = 1, K_REAL(N)
            H1 = SOLA_XPOS(J,O1,K,N)
            H2 = SOLB_XNEG(J,O1,K,N)
            HV_P(J,O1,K) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K) = QUAD_STRMWTS(J)*H2
          ENDDO

C  Complex homogeneous solution contributions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            HV_P(J,O1,K1) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K1,N)
            HV_P(J,O1,K2) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K2,N)
            HV_M(J,O1,K1) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K1,N)
            HV_M(J,O1,K2) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K2,N)
          ENDDO

C  End loops

        ENDDO
      ENDDO

C  Diffuse scatter contributions
C  -----------------------------

C  start loops

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM   = C0 + IROW

C  Beam contribution for this Kernel

          REFL_B = ZERO
          DO J = 1, NSTREAMS
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              REFL_B = REFL_B + PV_W(J,O2) * D_BIREFLEC(OM,B,Q,J,I)
            ENDDO
          ENDDO
          L_BEAM = REFL_B * FACTOR_BRDF

C  Real homogeneous solutions contribution to this Kernel

          L_HOM_R = ZERO
          DO K = 1, K_REAL(N)
            H_1 = ZERO
            H_2 = ZERO
            DO J = 1, NSTREAMS
              H_1_S = ZERO
              H_2_S = ZERO
              DO O2 = 1, NSTOKES
                OM  = MUELLER_INDEX(O1,O2)
                DBR = D_BIREFLEC(OM,B,Q,J,I)
                H_1_S = H_1_S + HV_P(J,O2,K) * DBR
                H_2_S = H_2_S + HV_M(J,O2,K) * DBR
              ENDDO
              H_1 = H_1 + H_1_S
              H_2 = H_2 + H_2_S
            ENDDO
            H_1 = FACTOR_BRDF * H_1
            H_2 = FACTOR_BRDF * H_2
            L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) 
     &                        + MCON(K,N) * H_2
          ENDDO

C  homogeneous complex solutions

          L_HOM_CR = ZERO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
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
                DBR = D_BIREFLEC(OM,B,Q,J,I)
                H_1_S_CR = H_1_S_CR + HV_P(J,O2,K1) * DBR
                H_2_S_CR = H_2_S_CR + HV_M(J,O2,K1) * DBR
                H_1_S_CI = H_1_S_CI + HV_P(J,O2,K2) * DBR
                H_2_S_CI = H_2_S_CI + HV_M(J,O2,K2) * DBR
              ENDDO
              H_1_CR = H_1_CR + H_1_S_CR
              H_2_CR = H_2_CR + H_2_S_CR
              H_1_CI = H_1_CI + H_1_S_CI
              H_2_CI = H_2_CI + H_2_S_CI
            ENDDO
            H_1_CR = FACTOR_BRDF  * H_1_CR
            H_1_CI = FACTOR_BRDF  * H_1_CI
            H_2_CR = FACTOR_BRDF  * H_2_CR
            H_2_CI = FACTOR_BRDF  * H_2_CI
            H1R =   H_1_CR * T_DELT_EIGEN(K1,N)
     &            - H_1_CI * T_DELT_EIGEN(K2,N)
            H1I =   H_1_CR * T_DELT_EIGEN(K2,N)
     &            + H_1_CI * T_DELT_EIGEN(K1,N)
            L_HOM_CR = L_HOM_CR
     &        + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I 
     &        + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI

c 4567  format   (3i4,1p2e20.10)  
c         write(21,4567) i,o1,k, H1I, H_2_CI

          ENDDO

C  Final contribution

          COL2_WFALB(CM,1) = L_BEAM + L_HOM_R + L_HOM_CR

C  End loops

        ENDDO
      ENDDO

C  Direct beam reflection

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        REFL_ATTN = BRDF_FACTORS(B) * ATMOS_ATTN(IB)
        REFL_B    = PARVALUE        * REFL_ATTN
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW
            OM = MUELLER_INDEX(O1,1)
            AWF_DIRECT = REFL_B * D_BIREFLEC_0(OM,B,Q,I,IB)
            COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + AWF_DIRECT
          ENDDO
        ENDDO
      ENDIF

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX,N, COL2_WFALB(N,1)
c        ENDDO
c      ENDIF
 
C  Finish

      RETURN
      END

C

      SUBROUTINE KPARAMS_BOA_SOURCE
     I        ( DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, BRDF_INDEX, PAR_INDEX,
     O          L_BOA_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_BRDF.VARS'

C  include files of Linearized solution and reflectance variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'

C  Module arguments
C  ----------------

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM
C  kernel index

      INTEGER          BRDF_INDEX

C  parameter index

      INTEGER          PAR_INDEX

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  output

      DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)

C  Local variables
C  ---------------

      INTEGER          UM, I, J, KL, N, IB, O1, O2, OM, B, Q
      INTEGER          K, KO1, K0, K1, K2      
      DOUBLE PRECISION INTEGRAND(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      DOUBLE PRECISION FACTOR_BRDF, H1, H2, PAR
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1

C  Initialise
C  ----------

C  Local indices

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      B   = BRDF_INDEX
      Q   = PAR_INDEX

C  initialise Derivative of BOA source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          L_BOA_SOURCE(UM,O1) = ZERO
        ENDDO
      ENDDO

C  Return if only quadrature

      IF ( .NOT. DO_USER_STREAMS ) RETURN

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

C  start loops

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES

C  Real homogeneous solutions

        SUM_R = ZERO
        DO K = 1, K_REAL(N)
         NXR = NCON_ALB(K,N) * SOLA_XPOS(I,O1,K,N)
         PXR = PCON_ALB(K,N) * SOLB_XNEG(I,O1,K,N)
         SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
        ENDDO

C  Complex solutions

        SUM_CR = ZERO
        DO K = 1, K_COMPLEX(N)
         K0 = 2 * K - 2
         K1 = KO1 + K0
         K2 = K1  + 1
         NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K1,N)
     &          - NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K2,N)
         NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K2,N)
     &          + NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K1,N)
         PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K1,N)
     &          - PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K2,N)
         H1 =  NXR1 *  * T_DELT_EIGEN(K1,N)
     &        -NXR2 *  * T_DELT_EIGEN(K2,N)
         H2 =  PXR1
         SUM_CR = SUM_CR + H1 + H2
        ENDDO

C  Final result

        INTEGRAND(I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

C  end loops

       ENDDO
      ENDDO

C  integrated BRDF reflectance term (loop over all kernels)
C  --------------------------------------------------------

      DO KL = 1, N_BRDF_KERNELS

        FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(KL)

C  .. Either, integrate Lambertian case, same for all user-streams

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) .AND.
     &                 FOURIER_COMPONENT.EQ.0 ) THEN
          O1 = 1
          REFLEC = ZERO
          DO J = 1, NSTREAMS
             REFLEC = REFLEC + INTEGRAND(J,O1)
          ENDDO
          REFLEC = FACTOR_BRDF * REFLEC
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
          ENDDO

C  .. Or, integrate with BRDF kernel function at user angles

        ELSE
  
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + INTEGRAND(J,O2) * 
     &                           USER_BIREFLEC(OM,KL,UM,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              REFLEC = FACTOR_BRDF * REFLEC
              L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
            ENDDO
          ENDDO

        ENDIF

C  End loop over all kernels

      ENDDO

C  Contributions due to direct variation of kernel parameter
C  ---------------------------------------------------------

C  Add linearization term for variation of the diffuse reflectance function

      PAR = BRDF_PARAMETERS(B,Q)
      FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(B) * PAR
      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            S_REFLEC = ZERO
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) *
     &                              D_USER_BIREFLEC(OM,B,Q,UM,J)
            ENDDO
            REFLEC = REFLEC + S_REFLEC
          ENDDO
          REFLEC = FACTOR_BRDF * REFLEC
          L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
        ENDDO
      ENDDO

C  Add linearization term for variation of direct beam reflectance 

      REFL_ATTN = BRDF_FACTORS(B) * ATMOS_ATTN(IB)
      FACTOR_BRDF = PAR * REFL_ATTN
      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          REFLEC = ZERO
          DO O2 = 1, NSTOKES
            OM = MUELLER_INDEX(O1,O2)
            REFLEC = REFLEC + 
     *              D_USER_BIREFLEC_0(OM,B,Q,UM,IB) * FLUXVEC(O2)
          ENDDO
          REFLEC = FACTOR_BRDF * REFLEC
          L_BOA_SOURCE(UM,O1) = L_BOA_SOURCE(UM,O1) + REFLEC
        ENDDO
      ENDDO

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)
C      IF ( DO_INCLUDE_SURFEMISS ) THEN
c        DO UM = LOCAL_UM_START, N_USER_STREAMS
c          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) - 
c     &                  FP_SURFBB * D_EMISSIVITY(B,Q,UM)
c        ENDDO
c      ENDIF

C  Finish

      RETURN
      END
