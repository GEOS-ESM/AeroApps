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
C #          SUBROUTINE thermal_setup                           #
C #                                                             #
C #   discrete ordinate solution                                #
C #          SUBROUTINE thermal_gfsolution                      #
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

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of thermal setup variables (output to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  Local variables
C  ---------------

      INTEGER          n, n1, s, ut, um, nt
      DOUBLE PRECISION HELP(MAXLAYERS)
      DOUBLE PRECISION T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION xtau, sum, cosmum

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

      DO ut = 1, n_offgrid_usertaus
       xtau = offgrid_utau_values(ut)
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

      IF ( n_thermal_coeffs == 3 ) THEN
       thermcoeffs(1,3) = zero
       thermcoeffs(1,2) = help(1)
       DO n = 1, nlayers - 1
        n1 = n + 1
        sum = ( (thermcoeffs(n,1)-thermcoeffs(n1,1)) / deltau_vert(n) )
     &         + help(n1)
        thermcoeffs(n1,3) = sum / ( deltau_vert(n1) + deltau_vert(n) )
        thermcoeffs(n1,2) = help(n1) - deltau_vert(n1)*thermcoeffs(n1,3)
       END DO
      END IF

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
        END IF

C  End user stream loop, and downwelling

       enddo
      endif

!  Finish

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE thermal_gfsolution

!  Green's function thermal particular integral, all layers.
!  Uses coefficient expansion of attenuation.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of Setup stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of Solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of thermal setup variables (output to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

! Local variables
! ---------------

      INTEGER          aa, aa1, i, i1, s, n, ut, nt
      DOUBLE PRECISION sum, help, s_p_u, s_p_l, s_m_u, s_m_l
      double precision sum_m, sum_p, tk, k1, tt
      double precision par1, par2, sd, su, spar

C  Multipliers (discrete ordinate solutions)

      DOUBLE PRECISION T_GMULT_UP(MAXSTREAMS)
      DOUBLE PRECISION T_GMULT_DN(MAXSTREAMS)
      DOUBLE PRECISION UT_T_GMULT_UP(MAXSTREAMS)
      DOUBLE PRECISION UT_T_GMULT_DN(MAXSTREAMS)

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

      do i = 1, nstreams_2
       do n = 1, nlayers
        t_wupper(i,n) = zero
        t_wlower(i,n) = zero
       enddo
      enddo

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
       END DO

!  Offgrid: only for quadrature or mean-value output

       IF ( do_quad_output .OR. do_additional_mvout ) THEN
        IF ( n_offgrid_usertaus .gt. 0 ) THEN
         DO ut = 1, n_offgrid_usertaus

C  regular off-grid solution

          n  = offgrid_utau_layeridx(ut)
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
! Green function solutions for all layers
! ---------------------------------------

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

! --------------------------
! Green function multipliers
! --------------------------

       DO aa = 1, nstreams
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
       END DO
 
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

C  End layer loop

      END DO

! -------------------------------------------------
! Offgrid Green's function multipliers and solution
! -------------------------------------------------

!  only for quadrature or mean-value output

      IF ( do_quad_output .OR. do_additional_mvout ) THEN
       IF ( n_offgrid_usertaus .gt. 0 ) THEN

!  start loop over offgrid optical depths

        DO ut = 1, n_offgrid_usertaus
         n  = offgrid_utau_layeridx(ut)

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

! finish off-grid solutions

        END DO
       END IF
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

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables (input/output to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  Local variables

      INTEGER          aa, um, n, ut, s, nt
      DOUBLE PRECISION sum_m, sum_p, su, sd, spar, cosmum, fac

C  Local multipliers

      DOUBLE PRECISION tsgm_uu(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION tsgm_ud(MAXLAYERS,0:MAX_THERMAL_COEFFS)

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

       do n = n_alllayers_up, nlayers
         layer_tsup_up(um,n) = fac * t_direct_up(um,n)
       enddo
       if ( do_user_taus ) then
        do ut = 1, n_offgrid_usertaus
         layer_tsup_utup(um,ut) = fac*t_ut_direct_up(um,ut)
        enddo
       endif

C  Finish if transmittance only

       if ( do_thermal_transonly ) go to 678

C  local cosine

       cosmum = user_streams(um)

C  Start index loop

       DO aa = 1, nstreams

C  Start layer loop

        DO n = n_alllayers_up, nlayers

C  Multipliers

         tsgm_uu(n,nt) = t_c_plus(aa,n,nt)
         tsgm_ud(n,nt) = t_c_minus(aa,n,nt)
         DO s = nt - 1, 1, -1
          tsgm_uu(n,s) = t_c_plus(aa,n,s)
     &          + s * cosmum * tsgm_uu(n,s+1)
          tsgm_ud(n,s) = t_c_minus(aa,n,s)
     &          + s * cosmum * tsgm_ud(n,s+1)
         END DO
         sum_p = zero  
         sum_m = zero  
         DO s = 1, nt
          sum_p = sum_p + tsgm_uu(n,s) * deltau_power(n,s)
          sum_m = sum_m + tsgm_ud(n,s) * deltau_power(n,s)
         END DO
         tsgm_uu(n,0) = - sum_p
         tsgm_ud(n,0) = - sum_m

C  Compute thermal diffuse term, add to WHOLE layer source

         su = t_c_plus(aa,n,0)  * hmult_1(aa,um,n)
         su = su + tsgm_uu(n,0) * t_delt_userm(n,um) 
     &           + tsgm_uu(n,1)
         su = su * tterm_save(aa,n)
         sd = t_c_minus(aa,n,0) * hmult_2(aa,um,n)
         sd = sd + tsgm_ud(n,0) * t_delt_userm(n,um) 
     &           + tsgm_ud(n,1)
         sd = sd * tterm_save(aa,n)
         spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
         layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar

C  End whole layer loop

        enddo

C  Compute thermal diffuse term, add to PARTIAL layer source term

        if ( do_user_taus ) then
         do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          sum_p = zero  
          sum_m = zero  
          do s = 1, nt
           sum_p = sum_p + tsgm_uu(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_ud(n,s) * xtau_power(ut,s)
          enddo
          su = t_c_plus(aa,n,0)  * ut_hmult_uu(aa,um,ut)
          su = su + tsgm_uu(n,0) * t_utup_userm(ut,um) + sum_p
          su = su * tterm_save(aa,n)
          sd = t_c_minus(aa,n,0) * ut_hmult_ud(aa,um,ut)
          sd = sd + tsgm_ud(n,0) * t_utup_userm(ut,um) + sum_m
          sd = sd * tterm_save(aa,n)
          spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
          layer_tsup_utup(um,ut) = layer_tsup_utup(um,ut) + spar
         enddo
        endif

C  End index loop 

       enddo

c       if ( do_user_taus ) then
c        do ut = 1, n_offgrid_usertaus
C   This checks out
c         do n = 1, nlayers
c          if ( um.eq.1)write(*,*)ut,
c     &        layer_tsup_utup(um,1),t_ut_direct_up(um,1)
c         enddo
c        enddo
c       endif

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

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of thermal variables (input/output to this module)

      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  Local variables

      INTEGER          aa, um, n, ut, s, nt
      DOUBLE PRECISION sum_m, sum_p, su, sd, spar, cosmum, fac

C  Local multipliers

      DOUBLE PRECISION tsgm_du(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION tsgm_dd(MAXLAYERS,0:MAX_THERMAL_COEFFS)

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
       if ( do_user_taus ) then
        do ut = 1, n_offgrid_usertaus
         layer_tsup_utdn(um,ut) = fac*t_ut_direct_dn(um,ut)
        enddo
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

C  Compute thermal diffuse term, add to WHOLE layer source

         sum_p = zero  
         sum_m = zero  
         DO s = 1, nt
          sum_p = sum_p + tsgm_du(n,s) * deltau_power(n,s)
          sum_m = sum_m + tsgm_dd(n,s) * deltau_power(n,s)
         END DO
         su = t_c_plus(aa,n,0)  * hmult_2(aa,um,n)
         su = su + tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
         su = su * tterm_save(aa,n)
         sd = t_c_minus(aa,n,0) * hmult_1(aa,um,n)
         sd = sd + tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
         sd = sd * tterm_save(aa,n)
         spar = fac * (  u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
         layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar

C  End whole layer loop

        enddo

C  Compute thermal diffuse term, add to PARTIAL layer source term

        if ( do_user_taus ) then
         DO ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          sum_p = zero  
          sum_m = zero  
          DO s = 1, nt
           sum_p = sum_p + tsgm_du(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_dd(n,s) * xtau_power(ut,s)
          END DO
          su = t_c_plus(aa,n,0)  * ut_hmult_du(aa,um,ut)
          su = su + tsgm_du(n,0) * t_utdn_userm(ut,um) + sum_p
          su = su * tterm_save(aa,n)
          sd = t_c_minus(aa,n,0) * ut_hmult_dd(aa,um,ut)
          sd = sd + tsgm_dd(n,0) * t_utdn_userm(ut,um) + sum_m
          sd = sd * tterm_save(aa,n)
          spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
          layer_tsup_utdn(um,ut) = layer_tsup_utdn(um,ut) + spar
         END DO
        endif

C  End index loop

       enddo

C   This checks out
c       do n = 1, nlayers
c        if (um.eq.1)write(*,*)n,layer_tsup_dn(um,n),t_direct_dn(um,n)
c       enddo

C  Continuation point for avoiding scattering calculations

 678   continue

C  End user-stream loop

      enddo

C  Finish

      RETURN
      END

