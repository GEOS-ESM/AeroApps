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
C #          SUBROUTINE thermal_setup_plus                      #
C #                                                             #
C #   discrete ordinate particular integrals                    #
C #          SUBROUTINE thermal_clsolution_plus                 #
C #                                                             #
C #   postprocessing source terms                               #
C #          SUBROUTINE thermal_sterms_up_plus                  #
C #          SUBROUTINE thermal_sterms_dn_plus                  #
C #                                                             #
C #   LTE Linearization   (new 9/11/09)                         #
C #          SUBROUTINE thermal_lte_linearization               #
C #                                                             #
C ###############################################################

      SUBROUTINE THERMAL_SETUP_PLUS

C  set-up of thermal expansion coefficients, always done after delta-m.
!  set-up of Linearization of thermal expansion coefficients.

!  Linearization w.r.t profile variables, not BB input.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  regular files
C  -------------

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of regular set up stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of thermal setup variables (input to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  linearized variables
C  --------------------

C  Include file of linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include file of linearizedd  set up stuff (input to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Local variables
C  ---------------

      INTEGER          n, s, ut, um, q, npars, nt
      DOUBLE PRECISION xtau, sum, cosmum, l_d, l_x
      DOUBLE PRECISION HELP(MAXLAYERS)
      DOUBLE PRECISION T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION L_T_MULT_UP
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_T_MULT_DN
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

C  powers of optical thickness
C  ---------------------------

C  Whole layer

      do n = 1, nlayers
       deltau_power(n,1) = one
       DO s = 2, n_thermal_coeffs
        deltau_power(n,s) = deltau_vert(n) * deltau_power(n,s-1)
       END DO
      enddo

      if ( do_atmos_linearization ) then
       do n = 1, nlayers
        do q = 1, layer_vary_number(n)
         l_deltau_power(n,1,q) = zero
         l_d = l_deltau_vert(q,n) * deltau_vert(n)
         do s = 2, n_thermal_coeffs
          l_deltau_power(n,s,q) = dble(s-1) * l_d * deltau_power(n,s-1)
         enddo
        enddo
       enddo
      endif

C  Partial layer

      if ( do_partlayers ) then
       DO ut = 1, n_partlayers

        xtau = partau_vert(ut)
        xtau_power(ut,1) = one
        DO s = 2, n_thermal_coeffs
         xtau_power(ut,s) = xtau * xtau_power(ut,s-1)
        END DO

        if ( do_atmos_linearization ) then
         n = partlayers_layeridx(ut)
         do q = 1, layer_vary_number(n)
          l_x = l_deltau_vert(q,n) *  partau_vert(ut)
          l_xtau_power(ut,1,q) = zero
          do s = 2, n_thermal_coeffs
           l_xtau_power(ut,s,q) = dble(s-1) * l_x * xtau_power(ut,s-1)
          enddo
         enddo
        endif

       enddo
      endif

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

      if ( n_thermal_coeffs .ne. 2 ) stop' number of coefficients not 2'

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

! linearized coefficients and auxiliary quantities

      if ( do_thermal_transonly ) then
       DO n = 1, nlayers
        if ( layer_vary_flag(n) ) then
         do q = 1, layer_vary_number(n) 
          l_thermcoeffs(n,1,q) = zero
          l_thermcoeffs(n,2,q) = - l_deltau_vert(q,n)*thermcoeffs(n,2)
          do s = 1, n_thermal_coeffs
           l_tcom1(n,s,q) =  l_thermcoeffs(n,s,q)
          enddo
         enddo
        else
         do s = 1, n_thermal_coeffs
           l_tcom1(n,s,1) = zero
         enddo
        endif
       enddo
      else
       DO n = 1, nlayers
        if ( layer_vary_flag(n) ) then
         do q = 1, layer_vary_number(n) 
          l_thermcoeffs(n,1,q) = zero
          l_thermcoeffs(n,2,q) = - l_deltau_vert(q,n)*thermcoeffs(n,2)
          do s = 1, n_thermal_coeffs
           l_tcom1(n,s,q) =  omegas1(n) * l_thermcoeffs(n,s,q)
     &       - thermcoeffs(n,s) * l_omega_total(q,n) * omega_total(n)
          enddo
         enddo
        else
         do s = 1, n_thermal_coeffs
           l_tcom1(n,s,1) = zero
         enddo
        endif
       enddo
      endif

! return if post-processing not flagged

      IF ( .NOT. do_user_streams ) RETURN

C  Short hand

      nt = n_thermal_coeffs

!  Upwelling Direct solution source terms
!  --------------------------------------

      IF ( do_upwelling ) THEN

C  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

C  start layer loop

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

C  Linearization of direct solution: whole layer source terms

         IF ( do_atmos_linearization ) then
          IF ( sterm_layermask_up(n).and.layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_mult_up(n,nt,q) = l_tcom1(n,nt,q)
            DO s = n_thermal_coeffs - 1, 1, -1
             l_t_mult_up(n,s,q) = l_tcom1(n,s,q) 
     &                  + s * cosmum * l_t_mult_up(n,s+1,q)
            END DO    
            sum = l_t_mult_up(n,1,q)
            DO s = 2, n_thermal_coeffs
             sum = sum + l_t_mult_up(n,s,q) *   deltau_power(n,s)
     &                 +   t_mult_up(n,s)   * l_deltau_power(n,s,q)
            END DO
            l_t_mult_up(n,0,q) = - sum
            l_t_direct_up(um,n,q) = l_t_mult_up(n,1,q)
     &            + l_t_mult_up(n,0,q) *   t_delt_userm(n,um)
     &            +   t_mult_up(n,0)   * l_t_delt_userm(n,um,q)
           END DO
          endif
         endif

C  End whole layer loop

        enddo

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          t_ut_direct_up(um,ut) = t_mult_up(n,0) * t_utup_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_up(um,ut) = t_ut_direct_up(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          IF ( layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_ut_direct_up(um,ut,q) =
     &            l_t_mult_up(n,0,q) *   t_utup_userm(ut,um)
     &           +  t_mult_up(n,0)   * l_t_utup_userm(ut,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_up(n,s,q) *   xtau_power(ut,s)
     &                 +   t_mult_up(n,s)   * l_xtau_power(ut,s,q)
            END DO    
            l_t_ut_direct_up(um,ut,q) = l_t_ut_direct_up(um,ut,q) + sum
           enddo
          END IF
         END DO
        endif

C  End user stream loop, and upwelling

       enddo
      endif

!  Downwelling Direct solution source terms
!  ----------------------------------------

      IF ( do_dnwelling ) THEN

C  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

C  start whole layer loop

        do n = 1, nlayers

         if ( sterm_layermask_dn(n) ) then
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
         endif

C  Linearization of direct solution: whole layer source terms

         if ( do_atmos_linearization ) then
          IF ( sterm_layermask_dn(n).and.layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_mult_dn(n,nt,q) = l_tcom1(n,nt,q)
            do s = nt - 1, 1, -1
             l_t_mult_dn(n,s,q) = l_tcom1(n,s,q) -
     &                  s * cosmum * l_t_mult_dn(n,s+1,q)
            enddo  
            l_t_mult_dn(n,0,q) = - l_t_mult_dn(n,1,q)
            l_t_direct_dn(um,n,q) =
     &           l_t_mult_dn(n,0,q) *   t_delt_userm(n,um)
     &           + t_mult_dn(n,0)   * l_t_delt_userm(n,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_dn(n,s,q) *   deltau_power(n,s)
     &                 +   t_mult_dn(n,s)   * l_deltau_power(n,s,q)
            END DO    
            l_t_direct_dn(um,n,q) = l_t_direct_dn(um,n,q) + sum
           enddo
          endif
         endif

C  End whole layer loop

        enddo

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          t_ut_direct_dn(um,ut) = t_mult_dn(n,0) * t_utdn_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_dn(um,ut) = t_ut_direct_dn(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          IF ( layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_ut_direct_dn(um,ut,q) =
     &             l_t_mult_dn(n,0,q) *   t_utdn_userm(ut,um)
     &             + t_mult_dn(n,0)   * l_t_utdn_userm(ut,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_dn(n,s,q) *   xtau_power(ut,s)
     &                 +   t_mult_dn(n,s)   * l_xtau_power(ut,s,q)
            END DO    
            l_t_ut_direct_dn(um,ut,q) = l_t_ut_direct_dn(um,ut,q) + sum
           ENDDO
          ENDIF
         ENDDO
        ENDIF

C  End user stream loop, and downwelling

       enddo
      endif

!  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_clsolution_plus (status,message)

!  Layer linearizations of Classical thermal particular integrals

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of Setup stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include files of Solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of thermal solution variables (input to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include file of linearized thermal solution variables (OUTPUT)

      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

!  Status arguments
!  ----------------

      integer          status
      character*(*)    message

! Local variables
! ---------------

      logical          do_vary
      INTEGER          aa, aa1, i, i1, j, j1, l, m, o1
      INTEGER          s, n, ut, nt, info, um, q, nvary
      double precision sum_m, sum_p, tk, k1, sd, su, sum
      double precision pos1, pos2, neg1, neg2, sum1, sum2, a5
      double precision term1, l_term1(max_atmoswfs), L_sum1, L_sum2
      double precision term2, l_term2(max_atmoswfs)
      double precision lsu, lsd, l_m, l_p, l_tk

C  Local matrices and vectors

      DOUBLE PRECISION HVEC1(MAXSTREAMS)
      DOUBLE PRECISION HVEC2(MAXSTREAMS)
      DOUBLE PRECISION JVEC1(MAXSTREAMS)
      DOUBLE PRECISION TVEC1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION TVEC2(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION TMAT(MAXSTREAMS,MAXSTREAMS)
      INTEGER          TPIVOT(MAXSTREAMS)

      DOUBLE PRECISION T_HELP1(0:MAXMOMENTS)
      DOUBLE PRECISION T_HELP2(0:MAXMOMENTS)
      DOUBLE PRECISION L_T_HELP1(0:MAXMOMENTS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_T_HELP2(0:MAXMOMENTS,MAX_ATMOSWFS)

      DOUBLE PRECISION LHVEC1(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION LHVEC2(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION LJVEC1(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_TVEC1(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_TVEC2(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION
     &   T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS),
     &   T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

      DOUBLE PRECISION L_T_C_MINUS
     &  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION L_T_C_PLUS
     &  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

      do i = 1, nstreams_2
       do n = 1, nlayers
        t_wupper(i,n) = zero
        t_wlower(i,n) = zero
       enddo
      enddo

      if ( do_atmos_linearization ) then
       do i = 1, nstreams_2
        do n = 1, nlayers
         do q = 1, layer_vary_number(n)
           l_t_wupper(i,n,q) = zero
           l_t_wlower(i,n,q) = zero
         enddo
        enddo
       enddo
      endif

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

C  Linearized transmittance-only thermal solutions

        if ( do_atmos_linearization ) then
         do q = 1, layer_vary_number(n)
          DO aa = 1, nstreams
           aa1 = aa + nstreams
           k1 = quad_streams(aa)
           tk = t_delt_disords(aa,n)
           l_tk = l_t_delt_disords(aa,n,q)
           l_t_c_minus(aa,n,nt,q)  =  k1 * l_thermcoeffs(n,nt,q)
           l_t_c_plus(aa,n,nt,q)   = l_t_c_minus(aa,n,nt,q)
           DO s = nt - 1, 1, -1
            l_t_c_minus(aa,n,s,q) =
     &       + k1 * (l_thermcoeffs(n,s,q) - s * l_t_c_minus(aa,n,s+1,q))
            l_t_c_plus(aa,n,s,q)  = 
     &       + k1 * (l_thermcoeffs(n,s,q) + s * l_t_c_plus(aa,n,s+1,q))
           END DO    
           l_p   = l_t_c_plus (aa,n,1,q)
           l_m   = l_t_c_minus(aa,n,1,q)
           DO s = 2, nt
            l_m = l_m + l_t_c_minus(aa,n,s,q) *   deltau_power(n,s)
     &                +   t_c_minus(aa,n,s)   * l_deltau_power(n,s,q) 
            l_p = l_p + l_t_c_plus(aa,n,s,q)  *   deltau_power(n,s)
     &                +   t_c_plus(aa,n,s)    * l_deltau_power(n,s,q)
           END DO
           l_t_c_minus(aa,n,0,q) = - l_t_c_minus(aa,n,1,q)
           l_t_c_plus(aa,n,0,q)  = - l_p
           l_t_wlower(aa,n,q) =
     &        (l_tk*t_c_minus(aa,n,0)+tk*l_t_c_minus(aa,n,0,q)+l_m)
           l_t_wupper(aa1,n,q) =
     &          ( l_tk*t_c_plus(aa,n,0) + tk*l_t_c_plus(aa,n,0,q)
     &                          + l_t_c_plus(aa,n,1,q) )
          END DO
         enddo
        endif

C  end layer loop

       END DO

!  Offgrid: only for quadrature or mean-value output

       IF ( do_quad_output .OR. do_additional_mvout ) THEN
        IF ( n_partlayers .gt. 0 ) THEN
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)

C  regular

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

C  Linearized

          if ( do_atmos_linearization ) then
           if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
             DO aa = 1, nstreams
              aa1 = aa + nstreams
              lsd = l_t_c_minus(aa,n,0,q) *   t_disords_utdn(aa,ut)
     &              + t_c_minus(aa,n,0)   * l_t_disords_utdn(aa,ut,q)
              lsu = l_t_c_plus(aa,n,0,q)  *   t_disords_utup(aa,ut)
     &              + t_c_plus(aa,n,0)    * l_t_disords_utup(aa,ut,q)
              DO s = 1, nt
               lsd = lsd + l_t_c_minus(aa,n,s,q) *   xtau_power(ut,s)
     &                   +   t_c_minus(aa,n,s)   * l_xtau_power(ut,s,q) 
               lsu = lsu + l_t_c_plus(aa,n,s,q)  * xtau_power(ut,s)
     &                   +   t_c_plus(aa,n,s)    * l_xtau_power(ut,s,q)
              END DO
              l_ut_t_partic(aa,ut,q)  = lsd 
              l_ut_t_partic(aa1,ut,q) = lsu 
             END DO
            enddo
           endif
          endif

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

C  Variation conditions

        do_vary = do_atmos_linearization.and.layer_vary_flag(n)
        nvary   = layer_vary_number(n)

C  Source constants

        TERM1 = two * tcom1(n,1)
        TERM2 = two * tcom1(n,2)
        if ( do_vary ) then
          do q = 1, nvary
            L_TERM1(q) = two * l_tcom1(n,1,q)
            L_TERM2(q) = two * l_tcom1(n,2,q)
          enddo
        endif

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

C  Linearization H Vector_1

        if ( do_vary ) then
          do q = 1, nvary
            do i = 1, nstreams
              sum = zero
              do j = 1, nstreams
                sum = sum -
     &               l_sab(i,j,o1,o1,n,q) * quad_streams(i) * hvec1(j)
              enddo
              lhvec1(i,q) = l_term1(q) + sum
            enddo
          enddo
          CALL DGETRS
     &    ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT,
     &     LHVEC1,MAXSTREAMS,INFO)            
          IF ( INFO .NE. 0 ) THEN
            MESSAGE = 
     & 'Classical Thermal solution LH1 back-substitution (DGETRS)'
            CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
            IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF
        endif

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

C  Linearization H Vector_2

        if ( do_vary ) then
          do q = 1, nvary
            do i = 1, nstreams
              sum = zero
              do j = 1, nstreams
                sum = sum -
     &               l_sab(i,j,o1,o1,n,q) * quad_streams(i) * hvec2(j)
              enddo
              lhvec2(i,q) = l_term2(q) + sum
            enddo
          enddo
          CALL DGETRS
     &    ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT,
     &     LHVEC2,MAXSTREAMS,INFO)            
          IF ( INFO .NE. 0 ) THEN
            MESSAGE = 
     & 'Classical Thermal solution LH2 back-substitution (DGETRS)'
            CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
            IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF
        endif

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

C  Linearization J Vector_1

        if ( do_vary ) then
          do q = 1, nvary
            do i = 1, nstreams
              sum = zero
              do j = 1, nstreams
                sum = sum + l_dab(i,j,o1,o1,n,q) * jvec1(j)
              enddo
              ljvec1(i,q) = lhvec2(i,q) + sum
            enddo
          enddo
          CALL DGETRS
     &    ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT,
     &     LJVEC1,MAXSTREAMS,INFO)            
          IF ( INFO .NE. 0 ) THEN
            MESSAGE = 
     & 'Classical Thermal solution LJ1 back-substitution (DGETRS)'
            CALL VLIDORT_LAPACK_ERROR ( INFO, N, MESSAGE, STATUS )
            IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF
        endif

C  Set solution
!  ------------

C  Expansion

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TVEC1(I,N)  = HALF * (HVEC1(I) + JVEC1(I))
          TVEC1(I1,N) = HALF * (HVEC1(I) - JVEC1(I))
          TVEC2(I,N)  = HALF * HVEC2(I)
          TVEC2(I1,N) = TVEC2(I,N)
        ENDDO

C  Linearized expansion

        if ( do_vary ) then
          do q = 1, nvary
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              L_TVEC1(I,N,Q)  = HALF * (LHVEC1(I,Q) + LJVEC1(I,Q))
              L_TVEC1(I1,N,Q) = HALF * (LHVEC1(I,Q) - LJVEC1(I,Q))
              L_TVEC2(I,N,Q)  = HALF * LHVEC2(I,Q)
              L_TVEC2(I1,N,Q) = L_TVEC2(I,N,Q)
            ENDDO
          enddo
        endif

C  Values at the layer boundaries

        DO I = 1, NSTREAMS_2
          T_WUPPER(I,N) = TVEC1(I,N)
          T_WLOWER(I,N) = TVEC1(I,N) + TVEC2(I,N) * DELTAU_POWER(N,2)
        ENDDO

C  Linearized Values at the layer boundaries

        if ( do_vary ) then
          do q = 1, nvary
            DO I = 1, NSTREAMS_2
              L_T_WUPPER(I,N,Q) = L_TVEC1(I,N,Q)
              L_T_WLOWER(I,N,Q) = L_TVEC1(I,N,Q)
     &             + L_TVEC2(I,N,Q) *   DELTAU_POWER(N,2)
     &             +   TVEC2(I,N)   * L_DELTAU_POWER(N,2,Q)
            ENDDO
          enddo
        endif

C  User solutions - help vectors

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

C  Linearized User solutions, help vectors

          if ( do_vary ) then
            do q = 1, nvary
              L_SUM1 = ZERO
              L_SUM2 = ZERO
              DO  J = 1, NSTREAMS
                J1 = J + NSTREAMS
                A5 = QUAD_HALFWTS(J)
                POS1 = L_TVEC1(J1,N,Q) * A5 * PI_XQP    (L,J,O1,O1)
                NEG1 = L_TVEC1(J,N,Q)  * A5 * PI_XQM_PRE(L,J,O1,O1)
                L_SUM1 = L_SUM1 + POS1 + NEG1
                POS2 = L_TVEC2(J,N,Q)  * A5 * PI_XQP    (L,J,O1,O1)
                NEG2 = L_TVEC2(J1,N,Q) * A5 * PI_XQM_PRE(L,J,O1,O1)
                L_SUM2 = L_SUM2 + POS2 + NEG2
              ENDDO
              L_T_HELP1(L,Q) = L_SUM1 *   OMEGA_GREEK(L,N,O1,O1) 
     &                         + SUM1 * L_OMEGA_GREEK(L,N,O1,O1,Q)
              L_T_HELP2(L,Q) = L_SUM2 *   OMEGA_GREEK(L,N,O1,O1)
     &                         + SUM2 * L_OMEGA_GREEK(L,N,O1,O1,Q)
            ENDDO
          ENDIF
        enddo

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
            if ( do_vary ) then
              do q = 1, nvary
                POS1 = ZERO
                POS2 = ZERO
                DO L = M, NMOMENTS
                  POS1 = POS1 + L_T_HELP1(L,Q) * PI_XUP(L,UM,O1,O1)
                  POS2 = POS2 + L_T_HELP2(L,Q) * PI_XUP(L,UM,O1,O1)
                ENDDO
                L_U_TPOS1(UM,N,Q) = POS1
                L_U_TPOS2(UM,N,Q) = POS2
              enddo
            endif
          ENDDO
        ENDIF

C  DNWELLING: sum over all harmonic contributions, each user stream

        IF ( DO_DNWELLING.and.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            NEG1 = ZERO
            NEG2 = ZERO
            DO L = M, NMOMENTS
              NEG1 = NEG1 + T_HELP1(L) * PI_XUM(L,UM,O1,O1)
              NEG2 = NEG2 + T_HELP2(L) * PI_XUM(L,UM,O1,O1)
            ENDDO
            U_TNEG1(UM,N) = NEG1
            U_TNEG2(UM,N) = NEG2
            if ( do_vary ) then
              do q = 1, nvary
                NEG1 = ZERO
                NEG2 = ZERO
                DO L = M, NMOMENTS
                  NEG1 = NEG1 + L_T_HELP1(L,Q) * PI_XUM(L,UM,O1,O1)
                  NEG2 = NEG2 + L_T_HELP2(L,Q) * PI_XUM(L,UM,O1,O1)
                ENDDO
                L_U_TNEG1(UM,N,Q) = NEG1
                L_U_TNEG2(UM,N,Q) = NEG2
              enddo
            endif
          ENDDO
        ENDIF

C  End layer loop

      ENDDO

!  Offgrid: only for quadrature or mean-value output
!  =================================================

      IF ( do_quad_output .OR. do_additional_mvout ) THEN
        IF ( n_PARTLAYERS .gt. 0 ) THEN
          DO ut = 1, n_PARTLAYERS
            n  = partlayers_layeridx(ut)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = 
     &            TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO 
            if ( layer_vary_flag(n) ) then
              do q = 1, layer_vary_number(n)
                DO I = 1, NSTREAMS_2
                  L_UT_T_PARTIC(I,UT,Q) = L_TVEC1(I,N,Q) 
     &              + L_TVEC2(I,N,Q) *   XTAU_POWER(UT,2)
     &              +   TVEC2(I,N)   * L_XTAU_POWER(UT,2,Q)
                ENDDO
              enddo
            endif
          ENDDO
        ENDIF
      END IF



C  Finish

      RETURN
      END

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_sterms_up_plus

C  thermal contributions to layer source terms (upwelling)
C  Linearized thermal contributions to layer source terms (upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of Linearized setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of thermal variables (input to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Local variables

      INTEGER          um, n, ut, s, nt, q
      DOUBLE PRECISION spar, cosmum, fac, sum

C  Multipliers (whole layer)

      DOUBLE PRECISION T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION L_T_MULT_UP
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

C  Particular solution layer source terms ( Green's function solution )
C  --------------------------------------------------------------------

C  Initial modulus = 4.pi if solar sources are included

      fac   = one
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
        do ut = 1, n_partlayers
         layer_tsup_utup(um,ut) = fac*t_ut_direct_up(um,ut)
        enddo
       endif

C  Linearized direct terms to start

       if ( do_atmos_linearization ) then
        do n = n_alllayers_up, nlayers
         do q = 1, layer_vary_number(n)
          l_layer_tsup_up(um,n,q) = fac * l_t_direct_up(um,n,q)
         enddo
        enddo
        if ( do_partlayers ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utup(um,ut,q) = fac*l_t_ut_direct_up(um,ut,q)
          enddo
         enddo
        endif
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  classical section particular integral
C  -------------------------------------

C  Layer loop

       do n = 1, nlayers
         if ( sterm_layermask_up(n) ) then

!  Whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

          t_mult_up(n,2) = u_tpos2(um,n)
          t_mult_up(n,1) = u_tpos1(um,n) + cosmum * t_mult_up(n,2)
          sum = t_mult_up(n,1)
          sum = sum + t_mult_up(n,2) * deltau_power(n,2)
          t_mult_up(n,0)   = - sum
          spar = t_mult_up(n,0) * t_delt_userm(n,um)
     &         + t_mult_up(n,1)
          layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar * fac

C  LInearized whole layer source terms

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            L_t_mult_up(n,2,q) = L_u_tpos2(um,n,q)
            L_t_mult_up(n,1,q) = L_u_tpos1(um,n,q) 
     &             + cosmum * L_t_mult_up(n,2,q)
            sum = L_t_mult_up(n,1,q)
            sum = sum + L_t_mult_up(n,2,q) *   deltau_power(n,2) +
     &                +   t_mult_up(n,2)   * L_deltau_power(n,2,q)
            L_t_mult_up(n,0,q)   = - sum
            spar = t_mult_up(n,0)   * L_t_delt_userm(n,um,q)
     &         + L_t_mult_up(n,0,q) *   t_delt_userm(n,um)
     &         + L_t_mult_up(n,1,q)
            L_layer_tsup_up(um,n,q) =
     &          L_layer_tsup_up(um,n,q) + spar * fac
           enddo
          endif
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
           if ( layer_vary_flag(n) ) then
             do q = 1, layer_vary_number(n)
              sum = L_t_mult_up(n,0,q) *   t_utup_userm(ut,um)
     &              + t_mult_up(n,0)   * L_t_utup_userm(ut,um,q)
              DO s = 1, nt
                sum = sum + L_t_mult_up(n,s,q) *   xtau_power(ut,s)
     &                    +   t_mult_up(n,s)   * L_xtau_power(ut,s,q)
              END DO    
              L_layer_tsup_utup(um,ut,q) = 
     &           L_layer_tsup_utup(um,ut,q) + sum*fac
             enddo
           endif
         enddo
       endif

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      END DO

C  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_sterms_dn_plus

!  thermal contributions to layer source terms (downwelling)
C  Linearized thermal contributions to layer source terms (downwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of Linearized setup/solution stuff (input to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of thermal variables (input to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Local variables

      INTEGER          um, n, ut, s, nt, q
      DOUBLE PRECISION spar, cosmum, fac, sum

C  Multipliers (whole layer)

      DOUBLE PRECISION T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION L_T_MULT_DN
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

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

       do n = 1, n_alllayers_dn
         layer_tsup_dn(um,n) = fac * t_direct_dn(um,n)
       enddo
       if ( do_partlayers ) then
        do ut = 1, n_partlayers
         layer_tsup_utdn(um,ut) = fac*t_ut_direct_dn(um,ut)
        enddo
       endif

C  Linearized direct terms to start

       if ( do_atmos_linearization ) then
        do n = 1, n_alllayers_dn
         do q = 1, layer_vary_number(n)
          l_layer_tsup_dn(um,n,q) = fac * l_t_direct_dn(um,n,q)
         enddo
        enddo
        if ( do_partlayers ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utdn(um,ut,q) = fac*l_t_ut_direct_dn(um,ut,q)
          enddo
         enddo
        endif
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  classical section particular integral
C  =====================================

       do n = 1, nlayers
         if ( sterm_layermask_dn(n) ) then

!  whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

          t_mult_dn(n,2) = u_tneg2(um,n)
          t_mult_dn(n,1) = u_tneg1(um,n) - cosmum * t_mult_dn(n,2)
          t_mult_dn(n,0) = - t_mult_dn(n,1)
          spar = t_mult_dn(n,0) * t_delt_userm(n,um)
          spar = spar + t_mult_dn(n,1)
          spar = spar + t_mult_dn(n,2) * deltau_power(n,2)
          layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar * fac

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            L_t_mult_dn(n,2,q) = L_u_tneg2(um,n,q)
            L_t_mult_dn(n,1,q) = L_u_tneg1(um,n,q) 
     &             - cosmum * L_t_mult_dn(n,2,q)
            L_t_mult_dn(n,0,q) = - L_t_mult_dn(n,1,q)
            spar =     t_mult_dn(n,0)   * L_t_delt_userm(n,um,q)
     &             + L_t_mult_dn(n,0,q) *   t_delt_userm(n,um)
            spar = spar + L_t_mult_dn(n,1,q)
            spar = spar + L_t_mult_dn(n,2,q) *   deltau_power(n,2)
     &                  +   t_mult_dn(n,2)   * L_deltau_power(n,2,q)
            L_layer_tsup_dn(um,n,q) = 
     &            L_layer_tsup_dn(um,n,q) + spar * fac
           enddo
          endif
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
           if ( layer_vary_flag(n) ) then
             do q = 1, layer_vary_number(n)
              sum = L_t_mult_dn(n,0,q) *   t_utdn_userm(ut,um)
     &              + t_mult_dn(n,0)   * L_t_utdn_userm(ut,um,q)
              DO s = 1, nt
                sum = sum + L_t_mult_dn(n,s,q) *   xtau_power(ut,s)
     &                    +   t_mult_dn(n,s)   * L_xtau_power(ut,s,q)
              END DO
              L_layer_tsup_utdn(um,ut,q) = 
     &             L_layer_tsup_utdn(um,ut,q) + sum*fac
             enddo
           endif
         END DO
       END IF

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      enddo

C  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE THERMAL_LTE_LINEARIZATION
     &   ( DO_INCLUDE_SURFACE, SURFACE_FACTOR, FLUX_MULTIPLIER )

!  set-up of LTE Linearization of thermal expansion coefficients.

!  Linearization w.r.t profile optical depth AND the BB input.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  regular files
C  -------------

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of regular setup and solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of thermal setup variables (input to this module)

      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  linearized variables
C  --------------------

C  Include file of linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Include file of linearized output variables

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Input arguments

      logical          do_include_surface
      double precision flux_multiplier
      double precision surface_factor

C  Local variables
C  ---------------

      INTEGER          k, i, ui, q, uta
      INTEGER          N, NC, NUT, NSTART, NUT_PREV, NLEVEL
      DOUBLE PRECISION mu, xi, term1, term2, kt, kt1, amb, bd, kd
      DOUBLE PRECISION LTE_kt, LTE_amb, LTE_bd, final_source
      DOUBLE PRECISION kmult, reflec

!  Local cumulative source

      DOUBLE PRECISION LTE_CUMSOURCE ( MAX_USER_STREAMS )

!  Linearized coefficients

      DOUBLE PRECISION LTE_TCOEFFS(2,MAXLAYERS,MAX_THERMAL_COEFFS)

!  Linearized transmittances

      DOUBLE PRECISION LTE_TRANS_DISORDS(2,MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION LTE_TRANS_USERM  (2,MAXLAYERS,MAX_USER_STREAMS)

!  Linearized source terms (user streams)

      DOUBLE PRECISION LTE_T_DIRECT_UP(2,MAXLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION LTE_T_DIRECT_DN(2,MAXLAYERS,MAX_USER_STREAMS)

!  Source terms (and linearized) for discrete ordinates

      DOUBLE PRECISION T_DISORDS_DN    (MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION LTE_T_DISORDS_DN(2,MAXLAYERS,MAXSTREAMS)

!  Surface reflection stuff

      DOUBLE PRECISION DOWNCUMS     (0:MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION LTE_DOWNSURF (MAXSTREAMS)
      DOUBLE PRECISION LTE_BOA_SOURCES (0:MAXLAYERS,MAXSTREAMS)

!  LTE-linearized coefficients
!  ---------------------------

      DO n = 1, nlayers
         lte_tcoeffs(1,n,1) = lte_thermal_bb_input(n-1)
         lte_tcoeffs(2,n,1) = zero
         lte_tcoeffs(1,n,2) = - lte_thermal_bb_input(n-1)
     &              - thermcoeffs(n,2) * lte_deltau_vert_input(1,n)
         lte_tcoeffs(2,n,2) = + lte_thermal_bb_input(n)
     &              - thermcoeffs(n,2) * lte_deltau_vert_input(2,n)
         lte_tcoeffs(1,n,2) = lte_tcoeffs(1,n,2)/deltau_vert_input(n)
         lte_tcoeffs(2,n,2) = lte_tcoeffs(2,n,2)/deltau_vert_input(n)
      ENDDO

!  LTE-linearization of transmittances
!  -----------------------------------

!  User streams, always need this

      DO ui = 1, n_user_streams
        mu = user_streams(ui)
        do n = 1, nlayers
          kt  = t_delt_userm(n,ui)
          kd = - kt / mu
          do q = 1, 2   
            LTE_trans_userm(q,n,ui) = kd * LTE_deltau_vert_input(q,n)
          enddo
        enddo
      enddo

!  Discrete ordinate stsreams, only for surface contribution for Upwelling

      if ( do_upwelling .and. do_include_surface ) then
        DO i = 1, nstreams
          xi = quad_streams(i)
          do n = 1, nlayers
            kt  = t_delt_disords(i,n)
            kd = - kt / xi
            do q = 1, 2   
              LTE_trans_disords(q,n,i) = kd * LTE_deltau_vert_input(q,n)
            enddo
          enddo
        enddo
      endif

!  Upwelling Direct solution source terms
!  --------------------------------------

      IF ( do_upwelling ) THEN
        DO ui = 1, n_user_streams
          mu = user_streams(ui)
          do n = 1, nlayers
            if ( sterm_layermask_up(n) ) then
              amb = thermcoeffs(n,1) + mu * thermcoeffs(n,2)
              bd  = thermcoeffs(n,2) * deltau_vert_input(n)
              kt  = t_delt_userm(n,ui)
              kt1 = one - kt
              do q = 1, 2
                LTE_kt  = LTE_trans_userm(q,n,ui)
                LTE_amb = LTE_tcoeffs(q,n,1) + mu * LTE_tcoeffs(q,n,2)
                LTE_bd  = thermcoeffs(n,2) * LTE_deltau_vert_input(q,n)
     &                + LTE_tcoeffs(q,n,2) * deltau_vert_input(n)
                term1 = - kt  * LTE_bd  - LTE_kt * bd
                term2 = + kt1 * LTE_amb - LTE_kt * amb
                LTE_t_direct_up(q,n,ui) = term1 + term2
              enddo
            endif
          enddo
        enddo
      endif

!  Downwelling Discrete Ordinate source terms
!  ------------------------------------------

!  These are required if there is surface reflectance

      if ( do_upwelling .and. do_include_surface ) then
        DO i = 1, nstreams
          xi = quad_streams(i)
          do n = 1, nlayers
            amb = thermcoeffs(n,1) - xi * thermcoeffs(n,2)
            bd  = thermcoeffs(n,2) * deltau_vert_input(n)
            kt  = t_delt_disords(i,n)
            kt1 = one - kt
            do q = 1, 2
              LTE_kt  = LTE_trans_disords(q,n,i)
              LTE_amb = LTE_tcoeffs(q,n,1) - xi * LTE_tcoeffs(q,n,2)
              LTE_bd  = thermcoeffs(n,2)   * LTE_deltau_vert_input(q,n)
     &                + LTE_tcoeffs(q,n,2) * deltau_vert_input(n)
              term1 = + LTE_bd 
              term2 = + kt1 * LTE_amb - LTE_kt * amb
              LTE_t_disords_dn(q,n,i) = term1 + term2
            enddo
            t_disords_dn(n,i) = amb * kt1 + bd
          enddo
        enddo
      endif

!  Downwelling Direct solution source terms
!  ----------------------------------------

      IF ( do_dnwelling ) THEN
        DO ui = 1, n_user_streams
          mu = user_streams(ui)
          do n = 1, nlayers
            if ( sterm_layermask_dn(n) ) then
              amb = thermcoeffs(n,1) - mu * thermcoeffs(n,2)
              bd  = thermcoeffs(n,2) * deltau_vert_input(n)
              kt  = t_delt_userm(n,ui)
              kt1 = one - kt
              do q = 1, 2
                LTE_kt  = LTE_trans_userm(q,n,ui)
                LTE_amb = LTE_tcoeffs(q,n,1) - mu * LTE_tcoeffs(q,n,2)
                LTE_bd  = thermcoeffs(n,2)  * LTE_deltau_vert_input(q,n)
     &                 + LTE_tcoeffs(q,n,2) * deltau_vert_input(n)
                term1 = + LTE_bd 
                term2 = + kt1 * LTE_amb - LTE_kt * amb
                LTE_t_direct_dn(q,n,ui) = term1 + term2
              enddo
            endif
          enddo
        enddo
      endif

!  ==============================
!  Now do the weighting functions
!  ==============================

!  For upwelling, have to do Linearization of STOKES_DOWNSURF first

      if ( do_upwelling ) then
        if ( do_include_surface ) then

! develop the cumulative sources

          do i = 1, nstreams
            downcums(0,i) = zero
          enddo
          do n = 1, nlayers
            do i = 1, nstreams
              downcums(n,i) = t_disords_dn(n,i) +
     &                   t_delt_disords(i,n) * downcums(n-1,i)
            enddo
          enddo

!  For each LTE weighting function

          do k = 0, nlayers

!  Linearize the downwelling solution

            do i = 1, nstreams
              LTE_downsurf(i) = zero
            enddo
            do n = 1, nlayers
              if ( k .eq. n ) then
                do i = 1, nstreams
                  LTE_downsurf(i) = LTE_t_disords_dn(2,n,i) +
     &                  t_delt_disords(i,n)  * LTE_downsurf(i) +
     &              LTE_trans_disords(2,n,i) * downcums(n-1,i)
                enddo
              else if ( k .eq. n-1 ) then
                do i = 1, nstreams
                  LTE_downsurf(i) = LTE_t_disords_dn(1,n,i) +
     &                  t_delt_disords(i,n)  * LTE_downsurf(i) +
     &              LTE_trans_disords(1,n,i) * downcums(n-1,i)
                enddo
              else
                do i = 1, nstreams
                  LTE_downsurf(i) = 
     &                  t_delt_disords(i,n)   * LTE_downsurf(i)
                enddo
              endif
            enddo

!  reflect off the surface (Lambertian only)

            kmult = surface_factor * lambertian_albedo
            reflec = 0.0d0
            do i = 1, nstreams
              reflec = reflec + quad_strmwts(I) * LTE_downsurf(i) 
            enddo
            reflec = reflec * kmult

!  Set the BOA source terms (Lambertian only)

            do ui = 1, n_user_streams
              LTE_boa_sources(k,ui) = reflec
            enddo

!  End LTE weighting function loop

          enddo

!  End surface reflection clause

        endif
      endif

!  Zero this contribution if no surface

      if ( .not. do_include_surface ) then
        do k = 0, nlayers
          do ui = 1, n_user_streams
            LTE_boa_sources(k,ui) = 0.0d0
          enddo
        enddo
      endif

!  Start Upwelling LTE weighting functions 

      if ( do_upwelling ) then

        do k = 0, nlayers

C  initialise cumulative source term loop

          NC  = 0
          NUT = 0
          NSTART   = NLAYERS
          NUT_PREV = NSTART + 1

C  initialise source

          DO UI = LOCAL_UM_START, N_USER_STREAMS
            LTE_CUMSOURCE(UI) = LTE_BOA_SOURCES(K,UI)
          ENDDO

C  loop over all output optical depths
C  -----------------------------------

          DO UTA = N_USER_LEVELS, 1, -1
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            NUT = NLEVEL + 1

C  Cumulative layer processing

            DO N = NSTART, NUT, -1
              NC = NLAYERS + 1 - N
              IF ( K .EQ. N ) THEN           
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_UP(2,N,UI) + 
     &           +     T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI)
     &           + LTE_TRANS_USERM(2,N,UI) * CUMSOURCE_UP(UI,1,NC-1)
                ENDDO
              ELSE IF ( K .EQ. N -1 ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_UP(1,N,UI) + 
     &           +     T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI)
     &           + LTE_TRANS_USERM(1,N,UI) * CUMSOURCE_UP(UI,1,NC-1)
                ENDDO
              ELSE
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = 
     &               + T_DELT_USERM(N,UI) * LTE_CUMSOURCE(UI)
                ENDDO
              ENDIF
            ENDDO

C  User-defined stream output, just set to the cumulative source term

            DO UI = LOCAL_UM_START, N_USER_STREAMS
              FINAL_SOURCE = FLUX_MULTIPLIER * LTE_CUMSOURCE(UI)
              LTE_ATMOSWF(K,UTA,UI,UPIDX) = FINAL_SOURCE
            ENDDO

C  Check for updating the recursion

            IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT

C  end loop over optical depth

          ENDDO

!  End upwelling LTE weighting functions

        ENDDO
      ENDIF

!  Start Downwelling LTE weighting functions
!  -----------------------------------------

      if ( do_dnwelling ) then

        do k = 0, nlayers

C  initialise cumulative source term loop

          NC  = 0
          NUT = 0
          NSTART   = 1
          NUT_PREV = NSTART - 1

C  initialise source

          DO UI = LOCAL_UM_START, N_USER_STREAMS
            LTE_CUMSOURCE(UI) = zero
          ENDDO

C  loop over all output optical depths

          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

C  Cumulative layer processing

            DO N = NSTART, NUT
              NC = N
              IF ( K .EQ. N ) THEN           
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_DN(2,N,UI) + 
     &           +     T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI)
     &           + LTE_TRANS_USERM(2,N,UI) * CUMSOURCE_DN(UI,1,NC-1)
                ENDDO
              ELSE IF ( K .EQ. N -1 ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_DN(1,N,UI) + 
     &           +     T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI)
     &           + LTE_TRANS_USERM(1,N,UI) * CUMSOURCE_DN(UI,1,NC-1)
                ENDDO
              ELSE
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = 
     &               + T_DELT_USERM(N,UI) * LTE_CUMSOURCE(UI)
                ENDDO
              ENDIF
            ENDDO

C  User-defined stream output, just set to the cumulative source term

            DO UI = LOCAL_UM_START, N_USER_STREAMS
              FINAL_SOURCE = FLUX_MULTIPLIER * LTE_CUMSOURCE(UI)
              LTE_ATMOSWF(K,UTA,UI,DNIDX) = FINAL_SOURCE
            ENDDO

C  Check for updating the recursion

            IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT

C  end loop over optical depth

          ENDDO

!  End Downwelling LTE weighting functions

        ENDDO
      ENDIF

!  Finish

      RETURN
      END

