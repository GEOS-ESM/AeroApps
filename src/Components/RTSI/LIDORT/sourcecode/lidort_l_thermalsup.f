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
C #   setups                                                    #
C #          SUBROUTINE thermal_setup_plus                      #
C #                                                             #
C #   discrete ordinate solution                                #
C #          SUBROUTINE thermal_gfsolution_plus                 #
C #                                                             #
C #   postprocessing source terms                               #
C #          SUBROUTINE thermal_sterms_up_plus                  #
C #          SUBROUTINE thermal_sterms_dn_plus                  #
C #                                                             #
C ###############################################################

      SUBROUTINE THERMAL_SETUP_PLUS

C  set-up of thermal expansion coefficients, always done after delta-m.
!  set-up of Linearization of thermal expansion coefficients.

!  Linearization w.r.t profile variables, not BB input.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  regular files
C  -------------

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of regular set up stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of thermal setup variables (input to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  linearized variables
C  --------------------

C  Include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include file of linearizedd  set up stuff (input to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'


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

      if ( do_user_taus ) then
       DO ut = 1, n_offgrid_usertaus

        xtau = offgrid_utau_values(ut)
        xtau_power(ut,1) = one
        DO s = 2, n_thermal_coeffs
         xtau_power(ut,s) = xtau * xtau_power(ut,s-1)
        END DO

        if ( do_atmos_linearization ) then
         n = offgrid_utau_layeridx(ut)
         do q = 1, layer_vary_number(n)
          l_x = l_deltau_vert(q,n) *  offgrid_utau_values(ut)
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

        if ( do_user_taus ) then
         DO ut = 1, n_offgrid_usertaus
          n  = offgrid_utau_layeridx(ut)
          t_ut_direct_up(um,ut) = t_mult_up(n,0) * t_utup_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_up(um,ut) = t_ut_direct_up(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_user_taus ) then
         DO ut = 1, n_offgrid_usertaus
          n  = offgrid_utau_layeridx(ut)
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

        if ( do_user_taus ) then
         DO ut = 1, n_offgrid_usertaus
          n  = offgrid_utau_layeridx(ut)
          t_ut_direct_dn(um,ut) = t_mult_dn(n,0) * t_utdn_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_dn(um,ut) = t_ut_direct_dn(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_user_taus ) then
         DO ut = 1, n_offgrid_usertaus
          n  = offgrid_utau_layeridx(ut)
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

      SUBROUTINE thermal_gfsolution_plus

!  Layer linearizations of Green's function thermal particular integrals

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of Setup stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of Solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of thermal solution variables (input to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include file of linearized thermal solution variables (OUTPUT)

      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

! Local variables
! ---------------

      INTEGER          aa, aa1, i, i1, s, n, ut, nt, q
      DOUBLE PRECISION sum, help, s_p_u, s_p_l, s_m_u, s_m_l
      double precision sum_m, sum_p, tk, k1, tt
      double precision par1, par2, spar, sovern, tterm, norm
      double precision l_k1, l_tk, l_tt, l_sum, l_norm
      double precision su, sd, lsu, lsd, l_m, l_p

C  Multipliers (discrete ordinate solutions)

      DOUBLE PRECISION T_GMULT_UP(MAXSTREAMS)
      DOUBLE PRECISION T_GMULT_DN(MAXSTREAMS)
      DOUBLE PRECISION UT_T_GMULT_UP(MAXSTREAMS)
      DOUBLE PRECISION UT_T_GMULT_DN(MAXSTREAMS)
      DOUBLE PRECISION L_T_GMULT_UP(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_T_GMULT_DN(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_UT_T_GMULT_UP(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_UT_T_GMULT_DN(MAXSTREAMS,MAX_ATMOSWFS)

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

C  shorthand

      nt = n_thermal_coeffs

C  Thermal Transmittance only, quadrature solutions
C  ================================================

      if ( do_thermal_transonly ) then

C  Whole layer solutions

       DO n = 1, nlayers
        DO aa = 1, nstreams
         aa1 = aa + nstreams
         k1 = x(aa)
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
           k1 = x(aa)
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
        IF ( n_offgrid_usertaus .gt. 0 ) THEN
         DO ut = 1, n_offgrid_usertaus
          n  = offgrid_utau_layeridx(ut)

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

! --------------------------------------------------
! Linearized Green function solutions for all layers
! --------------------------------------------------

C  start layer loop

      DO n = 1, nlayers

!  For each eigenstream, get the constant term

       DO aa = 1, nstreams
        sum = zero
        DO i = 1, nstreams
          i1 = i + nstreams
          help = a(i)*(xpos(i,aa,n)+xpos(i1,aa,n))
          sum  = sum + help
        END DO
        tterm_save(aa,n) = omegas1(n) * sum / norm_saved(n,aa)
       END DO

!  For each eigenstream, get the linearization of tterm_save

       IF ( layer_vary_flag(n) ) THEN
        DO aa = 1, nstreams
         tterm  = tterm_save(aa,n) 
         norm   = norm_saved(n,aa)
         sovern = omega_total(n) * tterm / omegas1(n)
         do q = 1, layer_vary_number(n)
          l_norm = l_norm_saved(aa,n,q)
          l_sum  = zero
          DO i = 1, nstreams
           i1 = i + nstreams
           help  = a(i)*(l_xpos(i,aa,n,q)+l_xpos(i1,aa,n,q))
           l_sum = l_sum + help
          END DO
          l_tterm_save(aa,n,q) = - l_omega_total(q,n) * sovern +
     &       ( omegas1(n)*l_sum - tterm*l_norm) / norm
         END DO
        enddo
       endif

C  Multiplier section
C  ------------------

C  start index loop

       DO aa = 1, nstreams

C  Regular Green function multipliers

        k1 = one / keigen(aa,n)
        tk = t_delt_eigen(aa,n)
        t_c_minus(aa,n,nt)  = k1 * thermcoeffs(n,nt)
        t_c_plus(aa,n,nt)   = k1 * thermcoeffs(n,nt)
        DO s = nt - 1, 1, -1
         t_c_minus(aa,n,s) = k1*(thermcoeffs(n,s)-s*t_c_minus(aa,n,s+1))
         t_c_plus(aa,n,s)  = k1*(thermcoeffs(n,s)+s*t_c_plus (aa,n,s+1))
        END DO    
        sum_p = t_c_plus (aa,n,1)
        sum_m = t_c_minus(aa,n,1)
        DO s = 2, nt
         sum_m = sum_m + t_c_minus(aa,n,s) * deltau_power(n,s)
         sum_p = sum_p + t_c_plus(aa,n,s)  * deltau_power(n,s)
        END DO
        tt = tterm_save(aa,n)
        t_c_minus(aa,n,0) = - t_c_minus(aa,n,1)
        t_c_plus(aa,n,0)  = - sum_p
        t_gmult_dn(aa) = tt*(tk*t_c_minus(aa,n,0) + sum_m)
        t_gmult_up(aa) = tt*(tk*t_c_plus(aa,n,0) + t_c_plus(aa,n,1))

C  Linearized Green function multipliers

        do q = 1, layer_vary_number(n)
          l_k1 = - k1 * k1 * l_keigen(aa,n,q)
          l_tk = l_t_delt_eigen(aa,n,q)
          l_t_c_minus(aa,n,nt,q)  = l_k1 *   thermcoeffs(n,nt)
     &                              + k1 * l_thermcoeffs(n,nt,q)
          l_t_c_plus(aa,n,nt,q)   = l_t_c_minus(aa,n,nt,q)
          DO s = nt - 1, 1, -1
           l_t_c_minus(aa,n,s,q) =
     &       l_k1 * (  thermcoeffs(n,s)   - s *   t_c_minus(aa,n,s+1))
     &       + k1 * (l_thermcoeffs(n,s,q) - s * l_t_c_minus(aa,n,s+1,q))
           l_t_c_plus(aa,n,s,q)  = 
     &       l_k1 * (  thermcoeffs(n,s)   + s *   t_c_plus(aa,n,s+1))
     &       + k1 * (l_thermcoeffs(n,s,q) + s * l_t_c_plus(aa,n,s+1,q))
          END DO    
          l_p   = l_t_c_plus (aa,n,1,q)
          l_m   = l_t_c_minus(aa,n,1,q)
          DO s = 2, nt
           l_m = l_m + l_t_c_minus(aa,n,s,q) *   deltau_power(n,s)
     &               +   t_c_minus(aa,n,s)   * l_deltau_power(n,s,q) 
           l_p = l_p + l_t_c_plus(aa,n,s,q)  *   deltau_power(n,s)
     &               +   t_c_plus(aa,n,s)    * l_deltau_power(n,s,q)
          END DO
          l_tt = l_tterm_save(aa,n,q)
          l_t_c_minus(aa,n,0,q) = - l_t_c_minus(aa,n,1,q)
          l_t_c_plus(aa,n,0,q)  = - l_p
          l_t_gmult_dn(aa,q) =
     &       l_tt * ( tk*t_c_minus(aa,n,0) + sum_m ) +
     &       tt * (l_tk*t_c_minus(aa,n,0)+tk*l_t_c_minus(aa,n,0,q)+l_m)
          l_t_gmult_up(aa,q) =
     &        l_tt * ( tk*t_c_plus(aa,n,0) + t_c_plus(aa,n,1) ) +
     &          tt * ( l_tk*t_c_plus(aa,n,0) + tk*l_t_c_plus(aa,n,0,q)
     &                          + l_t_c_plus(aa,n,1,q) )
        END DO

C  end index loop

       enddo

! -----------------------------------------------------
! Set particular integral from Green function expansion
! -----------------------------------------------------
 
       DO i = 1, nstreams
        i1 = i + nstreams
        s_p_u = zero
        s_p_l = zero
        s_m_u = zero
        s_m_l = zero
        DO aa = 1, nstreams
         s_p_u = s_p_u + t_gmult_up(aa)*xpos(i1,aa,n)
         s_m_u = s_m_u + t_gmult_up(aa)*xpos(i,aa,n)
         s_p_l = s_p_l + t_gmult_dn(aa)*xpos(i,aa,n)
         s_m_l = s_m_l + t_gmult_dn(aa)*xpos(i1,aa,n)
        END DO
        t_wupper(i,n)  = s_p_u
        t_wupper(i1,n) = s_m_u
        t_wlower(i1,n) = s_m_l
        t_wlower(i,n)  = s_p_l
       END DO

C  Set linearized form of particular integral at boundaries

       DO Q = 1, LAYER_VARY_NUMBER(N)
        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         S_P_U = ZERO
         S_P_L = ZERO
         S_M_U = ZERO
         S_M_L = ZERO
         DO AA = 1, NSTREAMS
          S_P_U = S_P_U + L_t_gmult_UP(AA,Q) *   XPOS(I1,AA,N) +
     &                      t_gmult_UP(AA)   * L_XPOS(I1,AA,N,Q)
          S_M_U = S_M_U + L_t_gmult_UP(AA,Q) *   XPOS(I,AA,N) +
     &                      t_gmult_UP(AA)    * L_XPOS(I,AA,N,Q)
          S_P_L = S_P_L + L_t_gmult_DN(AA,Q) *   XPOS(I,AA,N) +
     &                      t_gmult_DN(AA)   * L_XPOS(I,AA,N,Q)
          S_M_L = S_M_L + L_t_gmult_DN(AA,Q) *   XPOS(I1,AA,N) +
     &                      t_gmult_DN(AA)   * L_XPOS(I1,AA,N,Q)
         ENDDO
         L_T_WUPPER(I,N,Q)  = S_P_U
         L_T_WUPPER(I1,N,Q) = S_M_U
         L_T_WLOWER(I1,N,Q) = S_M_L
         L_T_WLOWER(I,N,Q)  = S_P_L
        ENDDO
       ENDDO

C  End layer loop

      END DO

! ------------------------------------------------------------
! Offgrid Green's function Linearized multipliers and solution
! ------------------------------------------------------------

!  only for quadrature or mean-value output

      IF ( do_quad_output .OR. do_additional_mvout ) THEN
       IF ( n_offgrid_usertaus .gt. 0 ) THEN

!  start loop over offgrid optical depths and parameters

        DO ut = 1, n_offgrid_usertaus
         n  = offgrid_utau_layeridx(ut)

!  Solutions
!  ---------

!  multipliers

         DO aa = 1, nstreams
          sd = t_c_minus(aa,n,0) * t_utdn_eigen(aa,ut)
          su = t_c_plus(aa,n,0)  * t_utup_eigen(aa,ut)
          DO s = 1, nt
           sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
           su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
          END DO
          ut_t_gmult_dn(aa) = sd * tterm_save(aa,n)
          ut_t_gmult_up(aa) = su * tterm_save(aa,n)
         END DO

!  upwelling solution

         IF ( do_upwelling ) THEN
          DO i = 1, nstreams
           i1 = i + nstreams
           spar = zero
           DO aa = 1, nstreams
            par1 = xpos(i,aa,n)  * ut_t_gmult_up(aa)
            par2 = xpos(i1,aa,n) * ut_t_gmult_dn(aa)
            spar = spar + par1 + par2
           END DO
           ut_t_partic(i1,ut) = spar
          END DO
         END IF

!  Downwelling solution

         IF ( do_dnwelling ) THEN
          DO i = 1, nstreams
           i1 = i + nstreams
           spar = zero
           DO aa = 1, nstreams
            par1 = xpos(i1,aa,n) * ut_t_gmult_up(aa)
            par2 = xpos(i,aa,n)  * ut_t_gmult_dn(aa)
            spar = spar + par1 + par2
           END DO
           ut_t_partic(i,ut) = spar
          END DO
         END IF

!  Linearized solutions
!  --------------------

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)

!  multipliers

           DO aa = 1, nstreams
            sd =  t_c_minus(aa,n,0) * t_utdn_eigen(aa,ut)
            su =  t_c_plus(aa,n,0) * t_utup_eigen(aa,ut)
            lsd = l_t_c_minus(aa,n,0,q) *   t_utdn_eigen(aa,ut)
     &            + t_c_minus(aa,n,0)   * l_t_utdn_eigen(aa,ut,q)
            lsu = l_t_c_plus(aa,n,0,q)  *   t_utup_eigen(aa,ut)
     &            + t_c_plus(aa,n,0)    * l_t_utup_eigen(aa,ut,q)
            DO s = 1, nt
             sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
             su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
             lsd = lsd + l_t_c_minus(aa,n,s,q) *   xtau_power(ut,s)
     &                 +   t_c_minus(aa,n,s)   * l_xtau_power(ut,s,q) 
             lsu = lsu + l_t_c_plus(aa,n,s,q)  * xtau_power(ut,s)
     &                 +   t_c_plus(aa,n,s)    * l_xtau_power(ut,s,q)
            END DO
            l_ut_t_gmult_dn(aa,q) = lsd *   tterm_save(aa,n) +
     &                               sd * l_tterm_save(aa,n,q)
            l_ut_t_gmult_up(aa,q) = lsu *   tterm_save(aa,n) +
     &                               su * l_tterm_save(aa,n,q)
           END DO

!  upwelling solution

           IF ( do_upwelling ) THEN
            DO i = 1, nstreams
             i1 = i + nstreams
             spar = zero
             DO aa = 1, nstreams
              par1 = l_xpos(i,aa,n,q)  *   ut_t_gmult_up(aa)
     &               + xpos(i,aa,n)    * l_ut_t_gmult_up(aa,q)
              par2 = l_xpos(i1,aa,n,q) *   ut_t_gmult_dn(aa)
     &               + xpos(i1,aa,n)   * l_ut_t_gmult_dn(aa,q)
              spar = spar + par1 + par2
             END DO
             l_ut_t_partic(i1,ut,q) = spar
            END DO
           END IF

!  Downwelling solution

           IF ( do_dnwelling ) THEN
            DO i = 1, nstreams
             i1 = i + nstreams
             spar = zero
             DO aa = 1, nstreams
              par1 = l_xpos(i1,aa,n,q) *   ut_t_gmult_up(aa)
     &               + xpos(i1,aa,n)   * l_ut_t_gmult_up(aa,q)
              par2 = l_xpos(i,aa,n,q)  *   ut_t_gmult_dn(aa)
     &               + xpos(i,aa,n)    * l_ut_t_gmult_dn(aa,q)
              spar = spar + par1 + par2
             END DO
             l_ut_t_partic(i,ut,q) = spar
            END DO
           END IF

C  Finish linearization

          enddo
         endif

! finish off-grid solutions

        END DO
       END IF
      END IF

C  Finish

      RETURN
      END

C

      SUBROUTINE thermal_sterms_up_plus

C  thermal contributions to layer source terms (upwelling)
C  Linearized thermal contributions to layer source terms (upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of Linearized setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  include file of thermal variables (input to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Local variables

      INTEGER          aa, um, n, ut, s, nt, q
      DOUBLE PRECISION sum_m, sum_p, su, sd, spar, cosmum, fac
      double precision sp, sm, lspar, lsu, lsd, su0, sd0

C  Local multipliers

      DOUBLE PRECISION l_tsgm_uu
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION l_tsgm_ud
     &         (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION tsgm_uu(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION tsgm_ud(MAXLAYERS,0:MAX_THERMAL_COEFFS)

C  Particular solution layer source terms ( Green's function solution )
C  --------------------------------------------------------------------

C  Initial modulus = 4.pi if solar sources are included

      fac   = one
      lspar = zero
      if ( do_solar_sources ) fac = pi4

C  shorthand

      nt = n_thermal_coeffs

C  Start user angle loop

      DO um = local_um_start, n_user_streams

C  Direct terms to start

       do n = n_alllayers_up, nlayers
         layer_tsup_up(um,n) = fac * t_direct_up(um,n)
       enddo
       if ( do_user_taus ) then
        do ut = 1, n_offgrid_usertaus
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
        if ( do_user_taus ) then
         do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utup(um,ut,q) = fac*l_t_ut_direct_up(um,ut,q)
          enddo
         enddo
        endif
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  local cosine

       cosmum = user_streams(um)

C  Start index loop

       do aa = 1, nstreams

C  Start layer loop

        DO n = n_alllayers_up, nlayers

C  Multipliers

         tsgm_uu(n,nt) = t_c_plus(aa,n,nt)
         tsgm_ud(n,nt) = t_c_minus(aa,n,nt)
         DO s = nt - 1, 1, -1
          tsgm_uu(n,s) = t_c_plus(aa,n,s)  + s * cosmum * tsgm_uu(n,s+1)
          tsgm_ud(n,s) = t_c_minus(aa,n,s) + s * cosmum * tsgm_ud(n,s+1)
         END DO
         sum_p = zero  
         sum_m = zero  
         DO s = 1, nt
          sum_p = sum_p + tsgm_uu(n,s) * deltau_power(n,s)
          sum_m = sum_m + tsgm_ud(n,s) * deltau_power(n,s)
         END DO
         tsgm_uu(n,0) = - sum_p
         tsgm_ud(n,0) = - sum_m

c  Linearized multipliers

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)
           l_tsgm_uu(n,nt,q) = l_t_c_plus(aa,n,nt,q)
           l_tsgm_ud(n,nt,q) = l_t_c_minus(aa,n,nt,q)
           DO s = nt - 1, 1, -1
            l_tsgm_uu(n,s,q) = l_t_c_plus(aa,n,s,q)
     &            + s * cosmum * l_tsgm_uu(n,s+1,q)
            l_tsgm_ud(n,s,q) = l_t_c_minus(aa,n,s,q)
     &            + s * cosmum * l_tsgm_ud(n,s+1,q)
           END DO
           sum_p = zero  
           sum_m = zero  
           DO s = 1, nt
            sum_p = sum_p + l_tsgm_uu(n,s,q) *   deltau_power(n,s)
     &                    +   tsgm_uu(n,s)   * l_deltau_power(n,s,q)
            sum_m = sum_m + l_tsgm_ud(n,s,q) *   deltau_power(n,s)
     &                    +   tsgm_ud(n,s)   * l_deltau_power(n,s,q)
           END DO
           l_tsgm_uu(n,0,q) = - sum_p
           l_tsgm_ud(n,0,q) = - sum_m
          END DO
         endif

C  Compute thermal diffuse term, add to WHOLE layer source

         su0 = t_c_plus(aa,n,0)  * hmult_1(aa,um,n)
         su0 = su0 + tsgm_uu(n,0) * t_delt_userm(n,um) + tsgm_uu(n,1)
         su  = su0 * tterm_save(aa,n)
         sd0 = t_c_minus(aa,n,0) * hmult_2(aa,um,n)
         sd0 = sd0 + tsgm_ud(n,0) * t_delt_userm(n,um) + tsgm_ud(n,1)
         sd  = sd0 * tterm_save(aa,n)
         spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
         layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar

C  Compute Linearized thermal diffuse term, add to WHOLE layer source

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)
           lsu = l_t_c_plus(aa,n,0,q) *   hmult_1(aa,um,n) +
     &             t_c_plus(aa,n,0)   * l_hmult_1(aa,um,n,q)
           lsu = lsu + l_tsgm_uu(n,0,q) *   t_delt_userm(n,um)
     &               +   tsgm_uu(n,0)   * l_t_delt_userm(n,um,q)
     &               + l_tsgm_uu(n,1,q)
           lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
           lsd = l_t_c_minus(aa,n,0,q) *   hmult_2(aa,um,n) +
     &             t_c_minus(aa,n,0)   * l_hmult_2(aa,um,n,q)
           lsd = lsd + l_tsgm_ud(n,0,q) *   t_delt_userm(n,um)
     &               +   tsgm_ud(n,0)   * l_t_delt_userm(n,um,q)
     &               + l_tsgm_ud(n,1,q)
           lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
           lspar =   l_u_xpos(um,aa,n,q)*sd + u_xpos(um,aa,n)*lsd
     &             + l_u_xneg(um,aa,n,q)*su + u_xneg(um,aa,n)*lsu
           lspar = lspar * fac
           l_layer_tsup_up(um,n,q) = l_layer_tsup_up(um,n,q) + lspar
          enddo
         endif

C  End whole layer loop

        enddo

C  Partial terms
C  =============

        if ( do_user_taus ) then
         do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)

C  Compute thermal diffuse term, add to PARTIAL layer source term

          sum_p = zero  
          sum_m = zero  
          do s = 1, nt
           sum_p = sum_p + tsgm_uu(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_ud(n,s) * xtau_power(ut,s)
          enddo
          su0 = t_c_plus(aa,n,0)  * ut_hmult_uu(aa,um,ut)
          su0 = su0 + tsgm_uu(n,0) * t_utup_userm(ut,um) + sum_p
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * ut_hmult_ud(aa,um,ut)
          sd0 = sd0 + tsgm_ud(n,0) * t_utup_userm(ut,um) + sum_m
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
          layer_tsup_utup(um,ut) = layer_tsup_utup(um,ut) + spar

C  Linearized terms

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            sp = zero  
            sm = zero  
            do s = 1, nt
             sp = sp + l_tsgm_uu(n,s,q) *   xtau_power(ut,s)
     &               +   tsgm_uu(n,s)   * l_xtau_power(ut,s,q)
             sm = sm + l_tsgm_ud(n,s,q) *   xtau_power(ut,s)
     &               +   tsgm_ud(n,s)   * l_xtau_power(ut,s,q)
            enddo
            lsu = l_t_c_plus(aa,n,0,q)  *   ut_hmult_uu(aa,um,ut)
     &            + t_c_plus(aa,n,0)    * l_ut_hmult_uu(aa,um,ut,q)
            lsu = lsu + sp + l_tsgm_uu(n,0,q) *   t_utup_userm(ut,um)
     &                     +   tsgm_uu(n,0)   * l_t_utup_userm(ut,um,q)
            lsu = lsu * tterm_save(aa,n)+ su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q)  *   ut_hmult_ud(aa,um,ut)
     &            + t_c_minus(aa,n,0)    * l_ut_hmult_ud(aa,um,ut,q)
            lsd = lsd + sm + l_tsgm_ud(n,0,q) *   t_utup_userm(ut,um)
     &                     +   tsgm_ud(n,0)   * l_t_utup_userm(ut,um,q)
            lsd = lsd * tterm_save(aa,n)+ sd0 * l_tterm_save(aa,n,q)
            lspar =   l_u_xpos(um,aa,n,q)*sd + u_xpos(um,aa,n)*lsd
     &              + l_u_xneg(um,aa,n,q)*su + u_xneg(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_utup(um,ut,q) =
     &                l_layer_tsup_utup(um,ut,q) + lspar
           enddo
          endif

C  Finish partial terms

         enddo
        endif

C  Finish index loop

       enddo

C  This checks out............. 1 March 2007
c       do n = 1, nlayers
c        if(um.eq.1)write(*,*)ut,
c     &        l_layer_tsup_utup(um,1,1),l_t_ut_direct_up(um,1,1),
c     &                layer_tsup_utup(um,1),t_ut_direct_up(um,1)
c       enddo

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

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of Linearized setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  include file of thermal variables (input to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include file of linearized thermal variables (output to this module)

      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  Local variables

      INTEGER          aa, um, n, ut, s, nt, q
      DOUBLE PRECISION sum_m, sum_p, su, sd, spar, cosmum, fac
      double precision sp, sm, lspar, lsu, lsd, su0, sd0

C  Local multipliers

      DOUBLE PRECISION L_TSGM_DU
     &       (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_TSGM_DD
     &       (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION TSGM_DU(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION TSGM_DD(MAXLAYERS,0:MAX_THERMAL_COEFFS)

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

       do n = 1, n_alllayers_dn
         layer_tsup_dn(um,n) = fac * t_direct_dn(um,n)
       enddo
       if ( do_user_taus ) then
        do ut = 1, n_offgrid_usertaus
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
        if ( do_user_taus ) then
         do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utdn(um,ut,q) = fac*l_t_ut_direct_dn(um,ut,q)
          enddo
         enddo
        endif
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  local cosine

       cosmum = user_streams(um)

C  Start index loop

       DO aa = 1, nstreams

C  Start layer loop

        DO n = 1, n_alllayers_dn

C  Multipliers

         tsgm_du(n,nt) = t_c_plus(aa,n,nt)
         tsgm_dd(n,nt) = t_c_minus(aa,n,nt)
         DO s = nt - 1, 1, -1
          tsgm_du(n,s) = t_c_plus(aa,n,s) 
     &           - s * cosmum * tsgm_du(n,s+1)
          tsgm_dd(n,s) = t_c_minus(aa,n,s)
     &           - s * cosmum * tsgm_dd(n,s+1)
         END DO
         tsgm_du(n,0) = - tsgm_du(n,1)
         tsgm_dd(n,0) = - tsgm_dd(n,1)

c  Linearized multipliers

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)
           l_tsgm_du(n,nt,q) = l_t_c_plus(aa,n,nt,q)
           l_tsgm_dd(n,nt,q) = l_t_c_minus(aa,n,nt,q)
           DO s = nt - 1, 1, -1
            l_tsgm_du(n,s,q) = l_t_c_plus(aa,n,s,q)
     &            - s * cosmum * l_tsgm_du(n,s+1,q)
            l_tsgm_dd(n,s,q) = l_t_c_minus(aa,n,s,q)
     &            - s * cosmum * l_tsgm_dd(n,s+1,q)
           END DO
           l_tsgm_du(n,0,q) = - l_tsgm_du(n,1,q)
           l_tsgm_dd(n,0,q) = - l_tsgm_dd(n,1,q)
          END DO
         endif

C  Compute thermal diffuse term, add to WHOLE layer source

         sum_p = zero  
         sum_m = zero  
         DO s = 1, nt
          sum_p = sum_p + tsgm_du(n,s) * deltau_power(n,s)
          sum_m = sum_m + tsgm_dd(n,s) * deltau_power(n,s)
         END DO
         su0 = t_c_plus(aa,n,0)  * hmult_2(aa,um,n)
         su0 = su0 + tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
         su  = su0 * tterm_save(aa,n)
         sd0 = t_c_minus(aa,n,0) * hmult_1(aa,um,n)
         sd0 = sd0 + tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
         sd  = sd0 * tterm_save(aa,n)
         spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
         layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar

C  Compute Linearized thermal diffuse term, add to WHOLE layer source

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)
           sp = zero  
           sm = zero  
           DO s = 1, nt
            sp = sp + l_tsgm_du(n,s,q) *   deltau_power(n,s)
     &              +   tsgm_du(n,s)   * l_deltau_power(n,s,q)
            sm = sm + l_tsgm_dd(n,s,q) *   deltau_power(n,s)
     &              +   tsgm_dd(n,s)   * l_deltau_power(n,s,q)
           END DO
           lsu = l_t_c_plus(aa,n,0,q) *   hmult_2(aa,um,n) +
     &             t_c_plus(aa,n,0)   * l_hmult_2(aa,um,n,q)
           lsu = lsu + sp + l_tsgm_du(n,0,q) *   t_delt_userm(n,um)
     &                    +   tsgm_du(n,0)   * l_t_delt_userm(n,um,q)
           lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
           lsd = l_t_c_minus(aa,n,0,q) *   hmult_1(aa,um,n) +
     &             t_c_minus(aa,n,0)   * l_hmult_1(aa,um,n,q)
           lsd = lsd + sm + l_tsgm_dd(n,0,q) *   t_delt_userm(n,um)
     &                    +   tsgm_dd(n,0)   * l_t_delt_userm(n,um,q)
           lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
           lspar =   l_u_xneg(um,aa,n,q)*sd + u_xneg(um,aa,n)*lsd
     &             + l_u_xpos(um,aa,n,q)*su + u_xpos(um,aa,n)*lsu
           lspar = lspar * fac
           l_layer_tsup_dn(um,n,q) = l_layer_tsup_dn(um,n,q) + lspar
          enddo
         endif

C  End whole layer loop

        enddo

C  Partial terms
C  =============

        if ( do_user_taus ) then
         do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)

C  Compute thermal diffuse term, add to PARTIAL layer source term

          sum_p = zero  
          sum_m = zero  
          DO s = 1, nt
           sum_p = sum_p + tsgm_du(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_dd(n,s) * xtau_power(ut,s)
          END DO
          su0 = t_c_plus(aa,n,0)  * ut_hmult_du(aa,um,ut)
          su0 = su0 + tsgm_du(n,0) * t_utdn_userm(ut,um) + sum_p
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * ut_hmult_dd(aa,um,ut)
          sd0 = sd0 + tsgm_dd(n,0) * t_utdn_userm(ut,um) + sum_m
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
          layer_tsup_utdn(um,ut) = layer_tsup_utdn(um,ut) + spar

C  Linearized terms

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            sp = zero  
            sm = zero  
            do s = 1, nt
             sp = sp + l_tsgm_du(n,s,q) *   xtau_power(ut,s)
     &               +   tsgm_du(n,s)   * l_xtau_power(ut,s,q)
             sm = sm + l_tsgm_dd(n,s,q) * xtau_power(ut,s)
     &               +   tsgm_dd(n,s)   * l_xtau_power(ut,s,q)
            enddo
            lsu = l_t_c_plus(aa,n,0,q)  *   ut_hmult_du(aa,um,ut)
     &            + t_c_plus(aa,n,0)    * l_ut_hmult_du(aa,um,ut,q)
            lsu = lsu + sp + l_tsgm_du(n,0,q) *   t_utdn_userm(ut,um)
     &                     +   tsgm_du(n,0)   * l_t_utdn_userm(ut,um,q) 
            lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q)  *   ut_hmult_dd(aa,um,ut)
     &            + t_c_minus(aa,n,0)    * l_ut_hmult_dd(aa,um,ut,q)
            lsd = lsd + sm + l_tsgm_dd(n,0,q) *   t_utdn_userm(ut,um)
     &                     +   tsgm_dd(n,0)   * l_t_utdn_userm(ut,um,q) 
            lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
            lspar =  l_u_xneg(um,aa,n,q)*sd + u_xneg(um,aa,n)*lsd
     &             + l_u_xpos(um,aa,n,q)*su + u_xpos(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_utdn(um,ut,q) =
     &                l_layer_tsup_utdn(um,ut,q) + lspar
           enddo
          endif

C  Finish partial terms

         enddo
        endif

C  Finish index loop

       enddo

c       do n = 1, nlayers
c        if(um.eq.1.and.n.eq.4)write(*,*)
c     &     n,l_layer_tsup_dn(um,n,1),l_t_direct_dn(um,n,1),
c     &                layer_tsup_dn(um,n),t_direct_dn(um,n)
c       enddo

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      END DO

C  Finish

      RETURN
      END
