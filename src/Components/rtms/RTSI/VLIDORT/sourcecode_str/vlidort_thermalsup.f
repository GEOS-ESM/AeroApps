C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
C #                 9 Channing Street                           #
C #                 Cambridge, MA 02138, USA                    #
C #                 Tel: (617) 492 1183                         #
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
C #   setups                                                    #
C #          SUBROUTINE thermal_setup                           #
C #                                                             #
C #   discrete ordinate particular integral                     #
C #          SUBROUTINE thermal_clsolution                      #
C #                                                             #
C #   postprocessing source terms                               #
C #          SUBROUTINE thermal_sterms_up                       #
C #          SUBROUTINE thermal_sterms_dn                       #
C #                                                             #
C ###############################################################

      SUBROUTINE THERMAL_SETUP

C  set-up of thermal expansion coefficients, always done after delta-m.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of thermal setup variables (output to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  Local variables
C  ---------------

      INTEGER          n, s, ut, um, nt
      DOUBLE PRECISION HELP(MAXLAYERS)
      DOUBLE PRECISION xtau, sum, cosmum

C  Multipliers (whole layer)

      DOUBLE PRECISION T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)

C  powers of optical thickness
C  ---------------------------

C  Whole layer

      DO n = 1, nlayers
       deltau_power(n,1) = one
       DO s = 2, n_thermal_coeffs
        deltau_power(n,s) = deltau_vert(n) * deltau_power(n,s-1)
       END DO
      END DO

C  Partial layer

      DO ut = 1, n_PARTLAYERS
       xtau = partau_vert(ut)
       xtau_power(ut,1) = one
       DO s = 2, n_thermal_coeffs
        xtau_power(ut,s) = xtau * xtau_power(ut,s-1)
       END DO
      END DO

!  initial set of coefficients

      DO n = 1, nlayers
       thermcoeffs(n,1) = thermal_bb_input(n-1)
       help(n) = 
     &    (thermal_bb_input(n)-thermal_bb_input(n-1)) / deltau_vert(n)
      END DO

!  piecewise continuous for linear regime

      IF ( n_thermal_coeffs == 2 ) THEN
       DO n = 1, nlayers
        thermcoeffs(n,2) = help(n)
       END DO
      END IF

!  Derivative continuity for quadratic regime
!    ( first layer is linear; better than using a Free derivative at TOA)
!  IF ( n_thermal_coeffs == 3 ) THEN
!     thermcoeffs(1,3) = zero
!     thermcoeffs(1,2) = help(1)
!     DO n = 1, n_comp_layers - 1
!        n1 = n + 1
!        thermcoeffs(n1,2) = thermcoeffs(n,2) + two * deltaus(n) * thermcoeffs(n,3)
!        thermcoeffs(n1,3) = ( help(n1) - thermcoeffs(n1,2) ) / deltaus(n1)
!     END DO
!  END IF

!  alternative scheme: backward looking quadratics
!    ( first layer is linear; better than using a Free derivative at TOA)

!      IF ( n_thermal_coeffs == 3 ) THEN
!       thermcoeffs(1,3) = zero
!       thermcoeffs(1,2) = help(1)
!       DO n = 1, nlayers - 1
!        n1 = n + 1
!        sum = ( (thermcoeffs(n,1)-thermcoeffs(n1,1)) / deltau_vert(n) )
!     &         + help(n1)
!       thermcoeffs(n1,3) = sum / ( deltau_vert(n1) + deltau_vert(n) )
!        thermcoeffs(n1,2) = help(n1) - deltau_vert(n1)*thermcoeffs(n1,3)
!       END DO
!      END IF

!  debug check

c      xtau = zero
c      DO n = 1, nlayers
c       sum = deltau_vert(n) / 10.0
c       DO s = 1, 10
c        xtau = xtau + sum
c        WRITE(34,'(1p2e15.6)')
c     &      xtau,thermcoeffs(n,1)+thermcoeffs(n,2)*s*sum+
c     &         thermcoeffs(n,3)*s*sum*s*sum
c       END DO
c      END DO

!  auxiliary quantities

      if ( do_thermal_transonly ) then
        DO n = 1, nlayers
          DO s = 1, n_thermal_coeffs
           tcom1(n,s) = thermcoeffs(n,s)
          END DO
        END DO
      else
        DO n = 1, nlayers
          omegas1(n) = one - omega_total(n)
          DO s = 1, n_thermal_coeffs
           tcom1(n,s) = thermcoeffs(n,s) * omegas1(n)
          END DO
        END DO
      endif

! return if post-processing not flagged

      IF ( .NOT. do_user_streams ) RETURN

C  Short hand

      nt = n_thermal_coeffs

!  Upwelling Direct solution source terms
! ---------------------------------------

      IF ( do_upwelling ) THEN

C  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

        do n = 1, nlayers
         if ( sterm_layermask_up(n) ) then
          t_mult_up(n,nt) = tcom1(n,nt)
          do s = nt - 1, 1, -1
           t_mult_up(n,s) = tcom1(n,s) + s * cosmum * t_mult_up(n,s+1)
          enddo    
          sum = t_mult_up(n,1)
          do s = 2, nt
           sum = sum + t_mult_up(n,s) * deltau_power(n,s)
          enddo
          t_mult_up(n,0) = - sum
          t_direct_up(um,n) = t_mult_up(n,0) * t_delt_userm(n,um)
     &                      + t_mult_up(n,1)
         endif
        enddo

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_PARTLAYERS
          n  = partlayers_layeridx(ut)
          t_ut_direct_up(um,ut) = t_mult_up(n,0) * t_utup_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_up(um,ut) = t_ut_direct_up(um,ut) + sum
         END DO
        endif

C  End user stream loop, and upwelling

       enddo
      endif

!  Downwelling Direct solution source terms
! -----------------------------------------

      IF ( do_dnwelling ) THEN

C  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

        DO n = 1, nlayers
         IF ( sterm_layermask_dn(n) ) THEN
          t_mult_dn(n,nt) = tcom1(n,nt)
          DO s = nt - 1, 1, -1
           t_mult_dn(n,s) = tcom1(n,s) - s * cosmum * t_mult_dn(n,s+1)
          END DO    
          t_mult_dn(n,0) = - t_mult_dn(n,1)
          t_direct_dn(um,n) = t_mult_dn(n,0) * t_delt_userm(n,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * deltau_power(n,s)
          END DO    
          t_direct_dn(um,n) = t_direct_dn(um,n) + sum
         END IF
        END DO

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_PARTLAYERS
          n  = partlayers_layeridx(ut)
          t_ut_direct_dn(um,ut) = t_mult_dn(n,0) * t_utdn_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_dn(um,ut) = t_ut_direct_dn(um,ut) + sum
         END DO
        END IF

C  End user stream loop, and downwelling

       enddo
      endif

!  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_clsolution (status, message)

!  Classical function thermal particular integral, all layers.
!  Uses coefficient expansion of attenuation.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of Setup stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of Solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of thermal setup variables (output to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  Status argument
!  ---------------

      integer          status
      character*(*)    message

! Local variables
! ---------------

      INTEGER          aa, aa1, i, i1, j, j1, l, m, o1
      INTEGER          s, n, nt, info, um, ut
      double precision sum_m, sum_p, tk, k1, sd, su, a5
      double precision term1, term2, sum1, sum2
      double precision pos1, pos2, neg1, neg2

C  Local matrices and vectors

      DOUBLE PRECISION
     &   T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS),
     &   T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

      DOUBLE PRECISION HVEC1(MAXSTREAMS)
      DOUBLE PRECISION HVEC2(MAXSTREAMS)
      DOUBLE PRECISION JVEC1(MAXSTREAMS)
      DOUBLE PRECISION TVEC1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION TVEC2(MAXSTREAMS_2,MAXLAYERS)

      DOUBLE PRECISION T_HELP1(0:MAXMOMENTS)
      DOUBLE PRECISION T_HELP2(0:MAXMOMENTS)

      DOUBLE PRECISION TMAT(MAXSTREAMS,MAXSTREAMS)
      INTEGER          TPIVOT(MAXSTREAMS)

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

      do i = 1, nstreams_2
       do n = 1, nlayers
        t_wupper(i,n) = zero
        t_wlower(i,n) = zero
       enddo
      enddo

C  (1,1) component only

      o1 = 1

C  initialize exception handling

      status  = vlidort_success
      message = ' '

C  shorthand

      nt = n_thermal_coeffs

C  Fourier component = 0 (always)

      m = 0

C  Thermal Transmittance only, quadrature solutions
C  ================================================

      if ( do_thermal_transonly ) then

C  Whole layer solutions

       DO n = 1, nlayers
        DO aa = 1, nstreams
         aa1 = aa + nstreams
         k1 = quad_streams(aa)
         tk = t_delt_disords(aa,n)
         t_c_minus(aa,n,nt)  = k1 * thermcoeffs(n,nt)
         t_c_plus(aa,n,nt)   = k1 * thermcoeffs(n,nt)
         DO s = nt - 1, 1, -1
          t_c_minus(aa,n,s)= k1*(thermcoeffs(n,s)-s*t_c_minus(aa,n,s+1))
          t_c_plus(aa,n,s) = k1*(thermcoeffs(n,s)+s*t_c_plus (aa,n,s+1))
         END DO    
         sum_p = t_c_plus (aa,n,1)
         sum_m = t_c_minus(aa,n,1)
         DO s = 2, nt
          sum_m = sum_m + t_c_minus(aa,n,s) * deltau_power(n,s)
          sum_p = sum_p + t_c_plus(aa,n,s)  * deltau_power(n,s)
         END DO
         t_c_minus(aa,n,0) = - t_c_minus(aa,n,1)
         t_c_plus(aa,n,0)  = - sum_p
         t_wlower(aa,n)  = tk * t_c_minus(aa,n,0) + sum_m
         t_wupper(aa1,n) = tk * t_c_plus(aa,n,0)  + t_c_plus(aa,n,1)
        END DO
       END DO

!  Offgrid: only for quadrature or mean-value output

       IF ( do_quad_output .OR. do_additional_mvout ) THEN
        IF ( n_PARTLAYERS .gt. 0 ) THEN
         DO ut = 1, n_PARTLAYERS

C  regular off-grid solution

          n  = partlayers_layeridx(ut)
          DO aa = 1, nstreams
           aa1 = aa + nstreams
           sd = t_c_minus(aa,n,0) * t_disords_utdn(aa,ut)
           su = t_c_plus(aa,n,0)  * t_disords_utup(aa,ut)
           DO s = 1, nt
            sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
            su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
           END DO
           ut_t_partic(aa,ut)  = sd
           ut_t_partic(aa1,ut) = su
          END DO

C  end offgrid output 

         ENDDO
        ENDIF
       END IF

C  return thermal transmittance-only

       RETURN

C  End thermal transmittance-only clause

      endif

! ---------------------------------------
! Classical solutions for all layers
! ---------------------------------------

      DO n = 1, nlayers

C  Source constants

        TERM1 = two * tcom1(n,1)
        TERM2 = two * tcom1(n,2)

C  H solutions
C  -----------

C  solution matrix for the reduced problem
C  ( matrix should be saved in the LU decomposition form)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = SAB(I,J,O1,O1,N) * QUAD_STREAMS(I)
          ENDDO
        ENDDO

C  L-U decomposition of the solution matrix

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
        IF ( INFO .NE. 0 ) THEN
          MESSAGE = 
     & 'Classical Thermal H solution LU decomposition (DGETRF)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  H Vector_1 and solution by back-substitution

        DO I = 1, NSTREAMS
          HVEC1(I) = TERM1
        ENDDO
        CALL DGETRS
     &    ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,
     &     HVEC1,MAXSTREAMS,INFO)
        IF ( INFO .NE. 0 ) THEN
          MESSAGE = 
     & 'Classical Thermal solution H1 back-substitution (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  H Vector_2 and solution by back-substitution

        DO I = 1, NSTREAMS
          HVEC2(I) = TERM2
        ENDDO
        CALL DGETRS
     &    ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,
     &     HVEC2,MAXSTREAMS,INFO)
        IF ( INFO .NE. 0 ) THEN
          MESSAGE = 
     & 'Classical Thermal solution H2 back-substitution (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  J solution
C  ----------

C  solution matrix and vector for the reduced problem
C  ( matrix should be saved in the LU decomposition form)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = - DAB(I,J,O1,O1,N)
          ENDDO
          JVEC1(I) = HVEC2(I)
        ENDDO

C  L-U decomposition of the solution matrix

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
        IF ( INFO .NE. 0 ) THEN
          MESSAGE = 
     & 'Classical Thermal J solution LU decomposition (DGETRF)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  J Vector_1 solution by back-substitution

        CALL DGETRS
     &    ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,
     &     JVEC1,MAXSTREAMS,INFO)
        IF ( INFO .NE. 0 ) THEN
          MESSAGE = 
     & 'Classical Thermal solution J1 back-substitution (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  Set solution
C  ============

C  Expansion

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TVEC1(I,N)  = HALF * (HVEC1(I) + JVEC1(I))
          TVEC1(I1,N) = HALF * (HVEC1(I) - JVEC1(I))
          TVEC2(I,N)  = HALF * HVEC2(I)
          TVEC2(I1,N) = TVEC2(I,N)
        ENDDO
 
C  Values at the layer boundaries

        DO I = 1, NSTREAMS_2
          T_WUPPER(I,N) = TVEC1(I,N)
          T_WLOWER(I,N) = TVEC1(I,N) + TVEC2(I,N) * DELTAU_POWER(N,2)
c          write(16,'(2i4,1p2e10.10)')I,N,TVEC1(I,N),TVEC2(I,N)
        ENDDO

C  User solutions

        DO L = M, NMOMENTS
          SUM1 = ZERO
          SUM2 = ZERO
          DO  J = 1, NSTREAMS
            A5 = QUAD_HALFWTS(J)
            J1 = J + NSTREAMS
            POS1 = TVEC1(J1,N) * A5 * PI_XQP    (L,J,O1,O1)
            NEG1 = TVEC1(J,N)  * A5 * PI_XQM_PRE(L,J,O1,O1)
            SUM1 = SUM1 + POS1 + NEG1
            POS2 = TVEC2(J,N)  * A5 * PI_XQP    (L,J,O1,O1)
            NEG2 = TVEC2(J1,N) * A5 * PI_XQM_PRE(L,J,O1,O1)
            SUM2 = SUM2 + POS2 + NEG2
          ENDDO
          T_HELP1(L) = SUM1 * OMEGA_GREEK(L,N,O1,O1)
          T_HELP2(L) = SUM2 * OMEGA_GREEK(L,N,O1,O1)
        ENDDO

C  UPWELLING: sum over all harmonic contributions, each user stream

        IF ( DO_UPWELLING.and.STERM_LAYERMASK_UP(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            POS1 = ZERO
            POS2 = ZERO
            DO L = M, NMOMENTS
              POS1 = POS1 + T_HELP1(L) * PI_XUP(L,UM,O1,O1)
              POS2 = POS2 + T_HELP2(L) * PI_XUP(L,UM,O1,O1)
            ENDDO
            U_TPOS1(UM,N) = POS1
            U_TPOS2(UM,N) = POS2
          ENDDO
        ENDIF

C  DNWELLING: sum over all harmonic contributions, each user stream

        IF ( DO_DNWELLING.and.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            NEG1 = ZERO
            NEG2 = ZERO
            DO L = M, NMOMENTS
              NEG1 = NEG1 + T_HELP1(L)*PI_XUM(L,UM,O1,O1)
              NEG2 = NEG2 + T_HELP2(L)*PI_XUM(L,UM,O1,O1)
            ENDDO
            U_TNEG1(UM,N) = NEG1
            U_TNEG2(UM,N) = NEG2
          ENDDO
        ENDIF

C  End layer loop

      ENDDO

!  Offgrid: only for quadrature or mean-value output

      IF ( do_quad_output .OR. do_additional_mvout ) THEN
        IF ( n_PARTLAYERS .gt. 0 ) THEN
          DO ut = 1, n_PARTLAYERS
            n  = partlayers_layeridx(ut)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = 
     &            TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO 
          ENDDO
        ENDIF
      END IF

C  Finish

      RETURN
      END

c

      SUBROUTINE thermal_sterms_up

C  thermal contributions to layer source terms (upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables (input/output to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  Local variables

      INTEGER          um, n, ut, s, nt
      DOUBLE PRECISION spar, cosmum, fac, sum

C  Multipliers (whole layer)

      DOUBLE PRECISION T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)

C  Particular solution layer source terms ( Green's function solution )
C  --------------------------------------------------------------------

C  Initial modulus = 4.pi if solar sources are included

      fac = one
      if ( do_solar_sources ) fac = pi4

C  shorthand

      nt = n_thermal_coeffs

C  Start user angle loop

      DO um = local_um_start, n_user_streams

C  local cosine

       cosmum = user_streams(um)

C  Direct terms to start

       do n = n_alllayers_up, nlayers
         layer_tsup_up(um,n) = fac * t_direct_up(um,n)
       enddo
       if ( do_partlayers ) then
        do ut = 1, n_PARTLAYERS
         layer_tsup_utup(um,ut) = fac*t_ut_direct_up(um,ut)
        enddo
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  classical section particular integral
C  -------------------------------------

!  Whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

       do n = 1, nlayers
         if ( sterm_layermask_up(n) ) then
          t_mult_up(n,2) = u_tpos2(um,n)
          t_mult_up(n,1) = u_tpos1(um,n) + cosmum * t_mult_up(n,2)
          sum = t_mult_up(n,1)
          do s = 2, nt
           sum = sum + t_mult_up(n,s) * deltau_power(n,s)
          enddo
          t_mult_up(n,0)   = - sum
          spar = t_mult_up(n,0) * t_delt_userm(n,um)
     &         + t_mult_up(n,1)
          layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar * fac
         endif
       enddo

C  Partial layer sources - add thermal diffuse term

       if ( do_partlayers ) then
         do ut = 1, n_PARTLAYERS
           n = partlayers_layeridx(ut)
           sum = t_mult_up(n,0) * t_utup_userm(ut,um)
           DO s = 1, nt
            sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
           END DO    
           layer_tsup_utup(um,ut) = 
     &           layer_tsup_utup(um,ut) + sum*fac
         enddo
       endif

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      enddo

C  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_sterms_dn

!  thermal contributions to layer source terms (downwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables (input/output to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  Local variables

      INTEGER          um, n, ut, s, nt
      DOUBLE PRECISION spar, cosmum, fac, sum

C  Multipliers (whole layer)

      DOUBLE PRECISION T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)

C  Particular solution layer source terms ( Green's function solution )
C  --------------------------------------------------------------------

C  Initial modulus = 4.pi if solar sources are included

      fac = one
      if ( do_solar_sources ) fac = pi4

C  shorthand

      nt = n_thermal_coeffs

C  Start user angle loop

      DO um = local_um_start, n_user_streams

C  Direct terms to start

       DO n = 1, n_alllayers_dn
         layer_tsup_dn(um,n) = fac * t_direct_dn(um,n)
       enddo
       if ( do_partlayers ) then
        do ut = 1, n_PARTLAYERS
         layer_tsup_utdn(um,ut) = fac*t_ut_direct_dn(um,ut)
        enddo
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  local cosine

       cosmum = user_streams(um)

C  classical section particular integral

!  whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

       do n = 1, nlayers
         if ( sterm_layermask_dn(n) ) then
          t_mult_dn(n,2) = u_tneg2(um,n)
          t_mult_dn(n,1) = u_tneg1(um,n) - cosmum * t_mult_dn(n,2)
          t_mult_dn(n,0) = - t_mult_dn(n,1)
          spar = t_mult_dn(n,0) * t_delt_userm(n,um)
          spar = spar + t_mult_dn(n,1)
          spar = spar + t_mult_dn(n,2) * deltau_power(n,2)
          layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar * fac
         endif
       enddo

C  Partial layer source terms, add thermal diffuse part

       if ( do_partlayers ) then
         DO ut = 1, n_PARTLAYERS
           n  = partlayers_layeridx(ut)
           sum = t_mult_dn(n,0) * t_utdn_userm(ut,um)
           DO s = 1, nt
             sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
           END DO    
           layer_tsup_utdn(um,ut) = 
     &           layer_tsup_utdn(um,ut) + sum*fac
         END DO
       END IF

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      enddo

C  Finish

      RETURN
      END

