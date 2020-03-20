
      PROGRAM test_2sJL

C  THIS IS A PURE LIDORT CALCULATION with the same 3-layer
C  inputs as the test_2sj program. Validation of 2-stream code.

      implicit none

C  LIDORT Include files
C  --------------------

C  Include file of dimensions and numbers

	INCLUDE '../includes/LIDORT.PARS'

C  Input variables

	INCLUDE '../includes/LIDORT_INPUTS.VARS'
	INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
	INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Result variables

	INCLUDE '../includes/LIDORT_RESULTS.VARS'
	INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  output status

      INTEGER        STATUS_INPUTCHECK
      INTEGER        STATUS_CALCULATION
      INTEGER        STATUS_INPUTREAD

C  Other variables

      integer n
      double precision abs_save(3),sca_save(3),asy_save(3),stream_value
      double precision abs(3), sca(3), eps, epsfac, column_save, column
      data abs_save / 0.05, 0.02, 0.01 /
      data sca_save / 0.15, 0.28, 0.34 /
      data asy_save / 0.70, 0.76, 0.77 /
      data column_save / 1.0d0 /

C  Call lidort input master

      call lidort_l_input_master
     +     ('test_2sjl.inp','logfile',status_inputread)
      if (status_inputread .ne. 0) STOP 'badread test 2sjl'

C  Epsfac

      eps = 1.0d-04
      epsfac = 1.0d0 + eps

      height_grid(0) = 36.0d0
      height_grid(1) = 20.0d0
      height_grid(2) = 10.0d0
      height_grid(3) = 0.0d0

          stream_value = 0.5d0
c        stream_value = dsqrt(1.0d0 / 3.0d0 )
        do_full_quadrature =  ( stream_value.ne.0.5d0) 

C  Baseline calculation
C  ====================

C  Optical property inputs

      column = column_save
      do n = 1, nlayers
        abs(n) = column * abs_save(n)
        sca(n) = sca_save(n)
        deltau_vert_input(n) = abs(n) + sca(n)
        omega_total_input(n) = sca(n) /  deltau_vert_input(n)
        phasmoms_total_input(0,n) = 1.0d0
        phasmoms_total_input(1,n) = 3.0d0 * asy_save(n)
       enddo

C  Linearization inputs, column WFs

       if ( do_column_linearization ) then
        do n = 1, nlayers
         l_deltau_vert_input(1,n) = abs_save(n)
         l_omega_total_input(1,n) =
     &        -abs_save(n)*omega_total_input(n)/deltau_vert_input(n)
         l_phasmoms_total_input(1,0,n) = 0.0d0
         l_phasmoms_total_input(1,1,n) = 0.0d0
         layer_vary_flag(n)   = .true.
         layer_vary_number(n) = 1
         n_totalcolumn_wfs = 1
        enddo
       endif

C  Linearization inputs, profile WFs

       if ( do_profile_linearization ) then
        do n = 1, nlayers
         l_deltau_vert_input(1,n) = column
         l_deltau_vert_input(2,n) = 1.0d0
         l_omega_total_input(1,n) = - column * 
     &                    omega_total_input(n)/deltau_vert_input(n)
         l_omega_total_input(2,n) =
     &       (1.0d0-omega_total_input(n))/deltau_vert_input(n)

         l_deltau_vert_input(1,n) = l_deltau_vert_input(1,n) 
     &                      * abs_save(n) / deltau_vert_input(n)
         l_deltau_vert_input(2,n) = l_deltau_vert_input(2,n) 
     &                      * sca_save(n) / deltau_vert_input(n)

         l_omega_total_input(1,n) = l_omega_total_input(1,n)
     &                      * abs_save(n) / omega_total_input(n)
         l_omega_total_input(2,n) = l_omega_total_input(2,n)
     &                      * sca_save(n) / omega_total_input(n)

         l_phasmoms_total_input(1,0,n) = 0.0d0
         l_phasmoms_total_input(1,1,n) = 0.0d0
         l_phasmoms_total_input(2,0,n) = 0.0d0
         l_phasmoms_total_input(2,1,n) = 0.0d0
         layer_vary_flag(n)   = .true.
         layer_vary_number(n) = 2
         n_totalcolumn_wfs  = 0
         n_totalprofile_wfs = 2
        enddo
       endif

C  Lidort call

       call lidort_l_master(status_inputcheck,status_calculation)
       call lidort_status(status_inputcheck,status_calculation)
       if (status_inputcheck .ne. 0) STOP 
     +		'internal input check failed'
       if (status_calculation .ne. 0) STOP 'calculation failed'

c      write(*,*)INTENSITY(1,1,1), INTENSITY(2,1,2)
c      write(*,*)COLUMNWF(1,1,1,1), COLUMNWF(1,2,1,2)
      write(*,*)PROFILEWF(1,1,1,1,1), PROFILEWF(1,1,2,1,2)

C  Finish

      STOP
      END PROGRAM test_2sJL
