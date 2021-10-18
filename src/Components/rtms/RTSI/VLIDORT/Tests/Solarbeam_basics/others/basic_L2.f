      PROGRAM basic_L2

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  input variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  result variables

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'
      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Arguments to the IOP preparation routine (all inputs)
C  -----------------------------------------------------

      double precision depol, aerssalb0, aerasymm0
      double precision deltau_mol(maxlayers)
      double precision deltau_ray(maxlayers)
      double precision deltau_aer(maxlayers)
      double precision temp(0:maxlayers), surftemp
      double precision greekmat_moms(0:500,maxlayers,16)

C  output arrays

      DOUBLE PRECISION RAD_BS
     &    (max_user_levels,max_geometries,4)
      DOUBLE PRECISION MRAD_BS (max_user_levels,maxbeams,4)
      DOUBLE PRECISION FLUX_BS (max_user_levels,maxbeams,4)

C  Exception handling variables (from VLIDORT routines)

      INTEGER          STATUS_INPUTREAD
      INTEGER          STATUS_INPUTCHECK
      INTEGER          STATUS_CALCULATION

C  help variables when needed

      INTEGER          ndum, n, k, j, ut, o1
      integer          ib, um, az, v, idr, lu
      double precision ddum
      character*80     filename, filename1
      logical          do_scalar_testing

C  Additional variables for timing tests

c      real             e1,e2,eb1,eb2,efd2,ebas

C  LAMBERTIAN ONLY

C  Read the input configuration file
C   number of layers will be set from the atmosphere data read

      CALL VLIDORT_L_INPUT_MASTER
     &    ('vlidort_v2p4_testing.inp_L2',
     &     'vlidort_v2p4_testing.log_L2',status_inputread)
      if ( status_inputread.ne.vlidort_success)
     &     stop'BAD READ of file: vlidort_v2p4_testing.inp_L2'

C  Set scalar testing flag

      do_scalar_testing =  ( nstokes .eq. 1 )

C  Upwelling or Downwelling

      if ( DO_UPWELLING ) IDR = 1
      if ( DO_DNWELLING ) IDR = 2

C  using Tomrad atmospheric input data:
C     26 layers, Rayleigh + O3 absorption, 1 UV wavelength
C     9 layers with additional aerosol

      open(2,file = 'aermoms_test.dat_rev',status='old')
      read(2,'(2I5)')    nlayers, ngreek_moments_input
      read(2,'(1pe14.5)')  depol
      read(2,'(1pe14.5)')  aerssalb0
      read(2,'(1pe14.5)')  aerasymm0
      read(2,'(i3,1pe14.5,44x,1pe10.2)')  ndum,height_grid(0),temp(0)
      do n = 1, nlayers
        read(2,'(i3,1pe14.5,1p2e15.7,1pe14.5,1pe10.2)')
     &      ndum, height_grid(n),
     &      deltau_mol(n), deltau_ray(n), deltau_aer(n),temp(n)
      enddo
      read(2,*)surftemp
      close(2)

C  Moments (Vector testing only)

      if ( .not.do_scalar_testing ) then
       open(2,file = 'aermoms_test.dat_moms',status='old')
       do n = 1, nlayers
        read(2,*) ndum, ddum, ddum
        do k = 1, 16
          read (2,*) (greekmat_moms(j,n,k),j = 0, ngreek_moments_input)
        enddo
       enddo
       close(2)
      endif

C  Open file, write headers

      if ( .not. do_mvout_only ) then
        filename = 'vlidort_v2p4_testing.out_L2'
        open(98,file=filename, status = 'unknown')
      endif
      if ( do_mvout_only .or. do_additional_mvout ) then
        filename1 = 'vlidort_v2p4_testing.mfout_L2'
        open(78,file=filename1, status = 'unknown')
      endif

C  BS CALCULATION 1
C  ================

C  Get the optical properties

      call prepare_2p4_testing_iops
     I   ( do_scalar_testing,  depol, aerssalb0, aerasymm0, 
     I     greekmat_moms, deltau_mol, deltau_ray, deltau_aer )

C  Call to VLIDORT_LMASTER with full linearization

c      call get_elapsed_time(e1)
      CALL VLIDORT_L_MASTER_BRDF
     &          ( STATUS_INPUTCHECK, STATUS_CALCULATION )
      CALL VLIDORT_STATUS ( STATUS_INPUTCHECK, STATUS_CALCULATION )
c      call get_elapsed_time(e2)
c      eb1 = e2-e1

C  Progress

      write(lu,'(a,50i3)')
     &    '   --- Done baseline calculation 2; # Fourier terms = ', 
     &        (FOURIER_SAVED(IB),IB = 1, N_SZANGLES)

C  Exception handling

      if (status_inputcheck.eq.vlidort_serious ) then
        write(*,*)'baseline calculation 2, input check failed'
        write(*,*)'ACTION: Check Log file for messages'
        write(*,*)'Stop'
        stop
      endif

C  store the Baseline output, standard radiances + Jacobians

      if ( .not. do_mvout_only ) then
       do o1 = 1, nstokes
        DO V = 1, N_GEOMETRIES
         DO UT = 1, N_USER_LEVELS
          RAD_BS(UT,V,O1) = STOKES(UT,V,O1,IDR)
         ENDDO
        ENDDO
       enddo
      endif

C  Store the baseline output (mean-I/Flux values + associated WFs)

      if ( do_mvout_only .or. do_additional_mvout ) then
       DO O1 = 1, NSTOKES
        DO V = 1, N_SZANGLES
         DO UT = 1, N_USER_LEVELS
          MRAD_BS(UT,V,O1) = MEAN_STOKES(UT,V,O1,IDR)
          FLUX_BS(UT,V,O1) = FLUX_STOKES(UT,V,O1,IDR)
         ENDDO
        ENDDO
       ENDDO
      endif

C  Write radiances

      if ( .not.do_mvout_only ) then
        DO IB = 1, N_SZANGLES
          DO UM = 1, N_USER_VZANGLES
           DO AZ = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + AZ
            write(98,'(a)')' '
            DO UT = 1, N_USER_LEVELS
             DO O1 = 1, NSTOKES
                write(98,677) O1,UT, USER_LEVELS(UT),'|',
     &    V, SZANGLES(IB),USER_VZANGLES(UM), USER_RELAZMS(AZ),
     &     '  |',RAD_BS(UT,V,O1)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDDO
      endif

C  Write mean radiances/Fluxes 

      if ( do_mvout_only . or. do_additional_mvout ) then
        DO V = 1, N_SZANGLES
          write(78,'(a)')' '
          do UT = 1, N_USER_LEVELS
           DO O1 = 1, NSTOKES
              write(78,678) O1,UT, USER_LEVELS(UT),'|', V, SZANGLES(V),
     &       '|',MRAD_BS(UT,V,O1), 
     &       '|',FLUX_BS(UT,V,O1)
           enddo
          ENDDO
        ENDDO
      endif


C  Output format statements

 677  format(4x,i1,3x,i2,1x,f7.3,1x,a1,i5,3f9.3,a3,1pe15.6)
 678  format
     &(4x,i1,3x,i2,1x,f7.3,1x,a1,i3,f9.3,1x,a1,1pe15.6,1x,a1,1pe15.6)

C  Close files

      if ( .not. do_mvout_only  )                   close(98)
      if ( do_mvout_only .or. do_additional_mvout ) close(78)

C  Stop

      stop
      end


      SUBROUTINE GET_ELAPSED_TIME (ELAPSED_TIME)

C       ----------------------------------------------------------------
C       output arguments
C       ----------------------------------------------------------------
        REAL            ELAPSED_TIME

C       ----------------------------------------------------------------
C       external functions
C       ----------------------------------------------------------------
        REAL            ETIME
        EXTERNAL        ETIME

C       ----------------------------------------------------------------
C       local variables
C       ----------------------------------------------------------------
        REAL            USER_SYS_TIME

C       ----------------------------------------------------------------
C       get elapsed time
C       ----------------------------------------------------------------
        ELAPSED_TIME = ETIME(USER_SYS_TIME)

C       ----------------------------------------------------------------
C       return
C       ----------------------------------------------------------------
        RETURN
        END
C

      subroutine prepare_2p4_testing_iops
     I   ( do_scalar_testing, depol, aerssalb0, aerasymm0,
     I     greekmat_moms, deltau_mol, deltau_ray, deltau_aer )

C  VLIDORT include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  VLIDORT Bookkeeping file

      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  VLIDORT include files with IOP variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Arguments (all inputs)
C  ----------------------

      logical          do_scalar_testing
      double precision depol, aerssalb0, aerasymm0
      double precision greekmat_moms(0:500,maxlayers,16)

      double precision deltau_mol(maxlayers)
      double precision deltau_ray(maxlayers)
      double precision deltau_aer(maxlayers)

C  Local variables
C  ---------------

      integer          n, l, k, ncoeffs, km(16), k1

      double precision deltau_gas(maxlayers)
      double precision colgas0, colray0, colaer0
      double precision colgas,  colray,  colaer
      double precision delgas, delray, delaer, delsca, deltau
      double precision aerswt, rayswt, raysmom2, omega, factor
      double precision aermom(0:maxmoments_input)
      double precision rayvmoms(0:2,16)
      double precision aerssalb, aerasymm, scaling

C  start of code
C  -------------

C  Ncoeffs and mask

      ncoeffs = 16
      do k = 1, ncoeffs
       km(k) = k
      enddo

C  Find total Column optical depths for GAS, RAY and AER

      colgas0 = 0.0d0
      colray0 = 0.0d0
      colaer0 = 0.0d0
      do n = 1, nlayers
        deltau_gas(n) = deltau_mol(n)-deltau_ray(n)
        colgas0 = colgas0 + deltau_gas(n)
        colray0 = colray0 + deltau_ray(n)
        colaer0 = colaer0 + deltau_aer(n)
      enddo

C  Column scaling
C    Set to Unity for consistency with column AND profile FD results
c      scaling = 1.1d0

      scaling = 1.0d0
      colgas = colgas0 * scaling
      colray = colray0 * scaling
      colaer = colaer0 * scaling

C  Perturbations for the Column Finite-Difference testing

      aerssalb = aerssalb0
      if ( do_scalar_testing ) then
        aerasymm = aerasymm0
      endif

C  Start layer loop

      do n = 1, nlayers

C  Set up basic IOP inputs
C  =======================

C  set up the layer optical densities. 

        delgas = deltau_gas(n)
        delray = deltau_ray(n)
        delaer = deltau_aer(n)

C  Perturbations for the profile finite difference testing

        aerssalb = aerssalb0 
        if ( do_scalar_testing ) then
          aerasymm = aerasymm0
        endif

C  total optical depths for extinction and scattering

        deltau = delgas + delray + delaer
        delsca = delray + aerssalb * delaer

C  single scattering albedo (total) - use toggle

        omega = delsca / deltau
        if ( omega.gt.0.99999d0)omega = 0.99999d0

C  Assign VLIDORT optical properties

        deltau_vert_input(n) = deltau
        omega_total_input(n) = omega

C  Scalar testing Rayleigh second moment

        if ( do_scalar_testing ) then
          raysmom2 = ( 1.0d0 - depol ) / (2.0d0 + depol )
        endif

C  Vector Rayleigh coefficients

        if ( .not.do_scalar_testing ) then
          do k = 1, 16
            do l = 0, 2
              rayvmoms(l,k) = greekmat_moms(l,1,k)
            enddo
          enddo
        endif

C  Phase function moments (Rayleigh only)

        if ( do_scalar_testing ) then
         if ( delaer .eq. 0.0d0 ) then
          greekmat_total_input(0,n,1) = 1.0d0
          greekmat_total_input(1,n,1) = 0.0d0
          greekmat_total_input(2,n,1) = raysmom2
          do l = 3, ngreek_moments_input
            greekmat_total_input(l,n,1) = 0.0d0
          enddo
         endif
        endif

C  SCALAR TESTING PHase function moments
C     (Mixture of Rayleigh + Henyey-Greenstein aerosol)

        if ( do_scalar_testing ) then
         if ( delaer .ne. 0.0d0 ) then
          aerswt  = aerssalb * delaer / delsca
          rayswt  = delray / delsca
          aermom(0) = 1.0d0
          aermom(1) = 3.0d0 * aerasymm
          aermom(2) = 5.0d0 * aerasymm * aermom(1) / 3.0d0
          greekmat_total_input(0,n,1) = 1.0d0
          greekmat_total_input(1,n,1) = aermom(1) * aerswt
          greekmat_total_input(2,n,1) = raysmom2*rayswt+aermom(2)*aerswt
          do l = 3, ngreek_moments_input
            factor = dble(2*l+1) / dble(2*l-1)
            aermom(l) = factor * aerasymm * aermom(l-1) 
            greekmat_total_input(l,n,1) = aermom(l) * aerswt
          enddo
         endif
        endif

C  VECTOR TESTING - use the coefficients from file

        if ( .not.do_scalar_testing ) then
         if ( delaer .eq. 0 ) then
          do k = 1, ncoeffs
            do l = 0, 2
              greekmat_total_input(l,n,k) = rayvmoms(l,km(k))
            enddo
           enddo
           greekmat_total_input(0,n,1) = 1.0d0
         else
          aerswt  = aerssalb * delaer / delsca
          rayswt  = delray / delsca
          do k = 1, ncoeffs
            k1 = km(k)
            do l = 0, 2
              greekmat_total_input(l,n,k) =
     &            rayvmoms(l,k1)*rayswt + greekmat_moms(l,n,k1)*aerswt
            enddo
            do l = 3, ngreek_moments_input
              greekmat_total_input(l,n,k) = greekmat_moms(l,n,k1)*aerswt
            enddo
          enddo
          greekmat_total_input(0,n,1) = 1.0d0
         endif
        endif

C  End layer loop

      enddo

      
C  FInish

      return
      end

