	PROGRAM test_lidort_2stream

        implicit none

C  LIDORT Include files
C  --------------------

C  Include file of dimensions and numbers

	INCLUDE '../includes/LIDORT.PARS'

C  Input variables

	INCLUDE '../includes/LIDORT_INPUTS.VARS'
	INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Result variables

	INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Parameters

        integer ndat,ntypes
        parameter(ndat = 26261,ntypes = 2)
        integer nspars
        parameter(nspars = 1)
        integer surftype
        parameter(surftype = 1)

C  Variables for read_prop

        integer ncoefsa(ntypes,2)
        character*15 scatfile(ntypes)
        double precision depol
        double precision wavnum(ndat)
        double precision taug(maxlayers,ndat)
        double precision taudp(maxlayers,ndat)
        double precision omega(maxlayers,ndat)
        double precision coefsr(0:2,6)
        double precision coefsa(ntypes,0:maxmoments_input,6,2)
        double precision spars(nspars)
        double precision heights(0:maxlayers)
        double precision theta,theta0,phi
        double precision fr(maxlayers,ndat)
        double precision fa(ntypes,maxlayers,ndat)

C  Variables for lidort_input_master

	integer status_inputread

C  Variables for lidort_master

	integer status_inputcheck,status_calculation

C  Variables for twostream_master

        double precision stream_value

C  Optical properties

        double precision DELTAU_VERT(MAXLAYERS,1)
        double precision OMEGA_TOTAL(MAXLAYERS,1)
        double precision ASYM(MAXLAYERS,1)
        double precision ALB(1)

C  Output
C  ------

C  Results

        double precision INTENSITY_TOA(1,1)
        double precision INTENSITY_BOA(1,1)

C  Other variables

        integer i,j,k1,l,n
        double precision coefsai(ntypes,0:maxmoments_input,6)
        double precision coefs(0:maxmoments_input,maxlayers,6)
c        double precision asymi(ntypes)

        depol = 0.0279d0

        scatfile(1) = 'cont.mom'
        scatfile(2) = 'oceanic.mom'

C  Call lidort input master

      	call lidort_input_master('test_lidort_2stream.inp','logfile',
     +				 status_inputread)
      	if (status_inputread .ne. 0) STOP 'badread test lidort 2stream'

C  Two choices of stream value.

C   LIDORT setting for single Gauss scheme is do_full_quadrature
C    If Stream value = 0.5,       then using double Gauss scheme over [0,1]
C    If Stream value = 1/SQRT(3), then using single Gauss scheme over [-1,1]

c          stream_value = 0.5d0
        stream_value = dsqrt(1.0d0 / 3.0d0 )
        do_full_quadrature =  ( stream_value.ne.0.5d0) 

C  Read in optical properties

c        call read_prop
c     +        (nlayers,ndat,ntypes,nspars,nmoments_input+1,depol,
c     +         scatfile,wavnum,taug,taudp,omega,coefsr,coefsa,ncoefsa,
c     +         spars,heights,theta,theta0,phi,fr,fa)

C  Set up lidort inputs

        beam_szas(1)         = 50.0
        user_angles_input(1) = 30.0
        user_relazms(1)      = 10.0
	lambertian_albedo    = 0.2
	brdf_factors(1)      = lambertian_albedo

C  Set up for two stream calculation

        ALB(1) = 0.2
        height_grid(0) = 60.0d0

C  wavelength loop
C  ===============

c        do i = 1, ndat
        do i = 1, 1

          do n = 1, nlayers
            read(87,'(1p5e24.12)')heights(n),taudp(n,1),omega(n,1),
     7        coefs(0,n,1),coefs(1,n,1)
          enddo

          do n = 1, nlayers
            height_grid(n) = heights(n)
            deltau_vert_input(nlayers-n+1) = taudp(n,i)
	    omega_total_input(nlayers-n+1) = omega(n,i)
	    if (dabs(omega_total_input(nlayers-n+1)-1.d0) .le. 1.d-6)
     +		omega_total_input(nlayers-n+1) = 0.999999d0
            DELTAU_VERT(NLAYERS-N+1,1) = taudp(n,i)
	    OMEGA_TOTAL(NLAYERS-N+1,1) = omega(n,i)
	    if (dabs(OMEGA_TOTAL(NLAYERS-N+1,1)-1.d0) .le. 1.d-6)
     +		OMEGA_TOTAL(NLAYERS-N+1,1) = 0.999999d0
            phasmoms_total_input(0:nmoments_input,nlayers-n+1) =
     +               coefs(0:nmoments_input,n,1)
          enddo
          phasmoms_total_input(0,:) = 1.d0

	  call lidort_master(status_inputcheck,status_calculation)
	  call lidort_status(status_inputcheck,status_calculation)
	  if (status_inputcheck .ne. 0) STOP 
     +		'internal input check failed'
	  if (status_calculation .ne. 0) STOP 'calculation failed'

          do n = 1, nlayers
            asym(nlayers-n+1,1) = coefs(1,n,1)/3.d0
          enddo

C  2stream call

	  CALL TWOSTREAM_MASTER
     +        ( DO_UPWELLING, DO_DNWELLING,
     +         DO_PLANE_PARALLEL,
     +         NLAYERS, 2*NLAYERS, NFINELAYERS, 1, NMOMENTS_INPUT,
     +         1, 1, 1, 1,
     +         STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     +         BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,
     +         SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, 1,
     +         ALB, DELTAU_VERT, OMEGA_TOTAL, ASYM,
     +         INTENSITY_TOA, INTENSITY_BOA,
     +         STATUS_INPUTCHECK, STATUS_CALCULATION)
	  if (status_inputcheck .ne. 0) STOP 
     +		'internal input check failed'
	  if (status_calculation .ne. 0) STOP 'calculation failed'

          write(*,*)'Stream value = ', stream_value
	  write(0,*) ' --lidort intensity   ',intensity(1,1,upidx)
	  write(0,*) ' --twostream intensity',intensity_toa(1,1)

	enddo

C  Finish

	STOP
	END PROGRAM test_lidort_2stream
