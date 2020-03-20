
      PROGRAM test_2sJ

      implicit none

C  integers

      integer nthreads
      integer nbeams, n_user_angles, n_user_relazms, ngeoms
      integer nlayers, nfinelayers,nmoments_input
      integer npars, nspars, surftype

C  Parameters

      parameter ( nthreads = 1 )

c      parameter ( nbeams = 5, n_user_angles = 2, n_user_relazms = 3)
      parameter ( nbeams = 1, n_user_angles = 1, n_user_relazms = 1)
      parameter ( ngeoms = nbeams * n_user_angles*n_user_relazms)

      parameter(nlayers = 3,nfinelayers = 2)
      parameter(nmoments_input = 1)

      parameter(npars = 2)
      parameter(nspars = 1)
      parameter(surftype = 1)

C  Flags

      LOGICAL          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL

C  Linearization flags

      LOGICAL          DO_PROFILE_WFS,  DO_COLUMN_WFS,  DO_SURFACE_WFS
      LOGICAL          DO_PROFILE_WFFD, DO_COLUMN_WFFD, DO_SURFACE_WFFD
      LOGICAL          DO_SIMULATION_ONLY

C  Linearization numbers

      LOGICAL          LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER          LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER          N_COLUMN_WFS, N_SURFACE_WFS

C  Geometry

      DOUBLE PRECISION BEAM_SZAS    ( nbeams )
      DOUBLE PRECISION USER_ANGLES  ( n_user_angles )
      DOUBLE PRECISION USER_RELAZMS ( n_user_relazms )

C  height and earth radius

      DOUBLE PRECISION EARTH_RADIUS
      DOUBLE PRECISION HEIGHT_GRID ( 0:NLAYERS )

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Geometry specification height

      DOUBLE PRECISION GEOMETRY_SPECHEIGHT

C  Optical properties

      double precision DELTAU_VERT(NLAYERS,NTHREADS)
      double precision OMEGA_TOTAL(NLAYERS,NTHREADS)
      double precision ASYM       (NLAYERS,NTHREADS)
      double precision ALB(NTHREADS)

C  Linearized optical properties

      DOUBLE PRECISION L_DELTAU_VERT(NLAYERS, NPARS, NTHREADS)
      DOUBLE PRECISION L_OMEGA_TOTAL(NLAYERS, NPARS, NTHREADS)
      DOUBLE PRECISION L_ASYM       (NLAYERS, NPARS, NTHREADS)

C  Output
C  ------

C  Results

      double precision INTENSITY_TOA(NGEOMS, NTHREADS)
      double precision INTENSITY_BOA(NGEOMS, NTHREADS)
      double precision INTENSITY_TOA_FD(NGEOMS, NTHREADS)
      double precision INTENSITY_BOA_FD(NGEOMS, NTHREADS)

      DOUBLE PRECISION PROFILEWF_TOA(NGEOMS,NLAYERS,NPARS,NTHREADS)
      DOUBLE PRECISION PROFILEWF_BOA(NGEOMS,NLAYERS,NPARS,NTHREADS)

      DOUBLE PRECISION COLUMNWF_TOA(NGEOMS,NPARS,NTHREADS)
      DOUBLE PRECISION COLUMNWF_BOA(NGEOMS,NPARS,NTHREADS)

      DOUBLE PRECISION SURFACEWF_TOA(NGEOMS,NSPARS,NTHREADS)
      DOUBLE PRECISION SURFACEWF_BOA(NGEOMS,NSPARS,NTHREADS)

C  output status

      INTEGER        STATUS_INPUTCHECK
      INTEGER        STATUS_CALCULATION

C  Output strings

      character*3  c3
      character*24 cstream

C  Other variables

      integer n, q, k, g
      double precision abs_save(3), sca_save(3), asy_save(3)
      double precision abs(3), sca(3), eps, epsfac, column
      double precision column_save, albedo_save
      data abs_save / 0.05, 0.02, 0.01 /
      data sca_save / 0.15, 0.28, 0.34 /
      data asy_save / 0.70, 0.76, 0.77 /
      data column_save / 1.0d0 /
      data albedo_save / 0.2d0 /

C  Set up inputs

      DO_UPWELLING      = .true.
      DO_DNWELLING      = .true.

c      DO_PLANE_PARALLEL = .true.
      DO_PLANE_PARALLEL = .false.

      flux_factor = 1.0d0
      geometry_specheight = 0.0d0

      beam_szas(1)   = 57.0d0
c      beam_szas(2)   = 57.0d0
c      beam_szas(3)   = 57.0d0
c      beam_szas(4)   = 57.0d0
c      beam_szas(5)   = 57.0d0
      user_angles(1) = 35.0d0
c      user_angles(2) = 35.0d0
c      user_angles(3) = 35.0d0
      user_relazms(1)= 10.0d0
c      user_relazms(2)= 10.0d0
c      user_relazms(3)= 10.0d0

      earth_radius = 6371.0d0
      height_grid(0) = 36.0d0
      height_grid(1) = 20.0d0
      height_grid(2) = 10.0d0
      height_grid(3) = 0.0d0

C  Two choices of stream value................ CHOOSE One !!!!!

      stream_value = dsqrt(1.0d0 / 3.0d0 )
c      stream_value = 0.5d0

C  Headers

      if ( stream_value .eq. 0.5d0 ) then
        c3 = 'SV1'
        cstream = '(Stream value = 0.50000)'
      else
        c3 = 'SV2'
        cstream = '(Stream value = 0.57735)'
      endif

C  Epsfac for finite differencing

      eps = 1.0d-04
      epsfac = 1.0d0 + eps

C  Baseline calculation
C  ====================

C  Linearization control

C  Choose this.......................................
      DO_PROFILE_wfs     = .true.
      DO_COLUMN_wfs      = .false.

C  Or choose this....................................
c      DO_PROFILE_wfs     = .false.
c      DO_COLUMN_wfs      = .true.

      DO_SURFACE_wfs     = .true.
c      DO_SURFACE_wfs     = .false.

      DO_SIMULATION_ONLY = .false.
      n_column_wfs  = 0
      n_surface_wfs = 0

C  FD values

      do_profile_wffd = do_profile_wfs
      do_column_wffd  = do_column_wfs
      do_surface_wffd = do_surface_wfs

C  Abort if both profile and column WFs are flagged

      if ( do_profile_wfs .and. do_column_wfs ) stop
     & ' * Cannot have both Profile and Column WFs in same run !'

C  Profile WFS
C  -----------

      if ( do_profile_wfs ) then

C  Optical property inputs

       alb(1) = albedo_save
       column = column_save
       do n = 1, nlayers
         abs(n) = column * abs_save(n)
         sca(n) = sca_save(n)
         deltau_vert(n,1) = abs(n) + sca(n)
         omega_total(n,1) = sca(n) /  deltau_vert(n,1)
         asym(n,1)        = asy_save(n)
       enddo

C  Linearization inputs, profile WFs

       do n = 1, nlayers
         l_deltau_vert(n,1,1) = column
         l_deltau_vert(n,2,1) = 1.0d0
         l_omega_total(n,1,1) = - column * 
     &                    omega_total(n,1)/deltau_vert(n,1)
         l_omega_total(n,2,1) =(1.0d0-omega_total(n,1))/deltau_vert(n,1)
         l_asym(n,1,1) = 0.0d0
         l_asym(n,2,1) = 0.0d0
         layer_vary_flag(n)   = .true.
         layer_vary_number(n) = 2
       enddo

C  Normalize, profile WFs

       do n = 1, nlayers
         l_deltau_vert(n,1,1) = l_deltau_vert(n,1,1) * abs_save(n)
         l_deltau_vert(n,2,1) = l_deltau_vert(n,2,1) * sca_save(n)
         l_omega_total(n,1,1) = l_omega_total(n,1,1) * abs_save(n)
         l_omega_total(n,2,1) = l_omega_total(n,2,1) * sca_save(n)
       enddo

C  end clause

      endif

C  Column WFs
C  -----------

      if ( do_column_wfs ) then

C  Optical property inputs

       alb(1) = albedo_save
       column = column_save
       do n = 1, nlayers
         abs(n) = column * abs_save(n)
         sca(n) = sca_save(n)
         deltau_vert(n,1) = abs(n) + sca(n)
         omega_total(n,1) = sca(n) /  deltau_vert(n,1)
         asym(n,1)        = asy_save(n)
       enddo

C  Linearization inputs, column WFs

       do n = 1, nlayers
         l_deltau_vert(n,1,1) = abs_save(n)
         l_omega_total(n,1,1) =
     &        -abs_save(n)*omega_total(n,1)/deltau_vert(n,1)
         l_asym(n,1,1) = 0.0d0
         layer_vary_flag(n)   = .true.
         layer_vary_number(n) = 1
         n_column_wfs = 1
       enddo

C  Normalize, column WFs

       do n = 1, nlayers
         l_deltau_vert(n,1,1) = l_deltau_vert(n,1,1) * column
         l_omega_total(n,1,1) = l_omega_total(n,1,1) * column
       enddo

C  end clause

      endif

C  basic calculation with linearizations
C  -------------------------------------

      call TWOSTREAM_L_MASTER
     I   ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     I     DO_PROFILE_WFS, DO_COLUMN_WFS,
     I     DO_SURFACE_WFS, DO_SIMULATION_ONLY,
     I     NLAYERS, 2*NLAYERS, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     I     NBEAMS, N_USER_ANGLES, N_USER_RELAZMS,
     I     Ngeoms, npars, nspars,
     I     LAYER_VARY_FLAG, LAYER_VARY_NUMBER,
     I     n_column_wfs, n_surface_wfs,
     I     STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     I     BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     I     SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, 1,
     I     ALB, DELTAU_VERT, OMEGA_TOTAL, ASYM,
     I     L_DELTAU_VERT, L_OMEGA_TOTAL, L_ASYM,
     O     INTENSITY_TOA, INTENSITY_BOA,
     O     PROFILEWF_TOA, PROFILEWF_BOA,
     O     COLUMNWF_TOA,  COLUMNWF_BOA,
     O     SURFACEWF_TOA, SURFACEWF_BOA,
     O     STATUS_INPUTCHECK, STATUS_CALCULATION )

C  Debug
c      write(*,*)COLUMNWF_TOA(1,1,1), COLUMNWF_BOA(1,1,1)
c      write(*,*)PROFILEWF_TOA(1,1,1,1), PROFILEWF_BOA(1,1,1,1)

C  FD testing initialization
C  =========================

C  Linearization control

      DO_PROFILE_wfs     = .false.
      DO_COLUMN_wfs      = .false.
      DO_SURFACE_wfs     = .false.
      DO_SIMULATION_ONLY = .true.
      do n = 1, nlayers
        do q = 1, npars
         l_deltau_vert(n,q,1) = 0.0d0
         l_omega_total(n,q,1) = 0.0d0
         l_asym(n,q,1)        = 0.0d0
        enddo
        layer_vary_flag(n)   = .false.
        layer_vary_number(n) = 0
        n_column_wfs  = 0
        n_surface_wfs = 0
      enddo

C  Open file

      if ( do_profile_wffd .and. do_surface_wffd ) then
       open(55,file='profile_surface_wfs_'//c3//'.out',status='unknown')
      else if ( do_column_wffd .and. do_surface_wffd ) then
       open(55,file='column_surface_wfs_'//c3//'.out',status='unknown')
      endif

C  FD testing for PROFILE WFS
C  ==========================

C  Skip to column section if not flagged

      if ( .not. do_profile_wffd ) go to 3456

C  write header

      write(55,'(a//a,a/)')  'Profile WFS '//cstream//':',
     &  '      geom#  WF# Layer#  ',
     &  '  2s-analytic       2s-Findiff.       %diff'

C  Start FD loops

      do k = 1, nlayers
       do q = 1, 2
c      do k = 1, 1
c       do q = 1, 1
    
C  Optical property inputs

        alb(1) = albedo_save
        column = column_save
        do n = 1, nlayers
         abs(n) = column * abs_save(n)
         sca(n) = sca_save(n)
         if ( q.eq.1 .and. n.eq.k) abs(n) = abs(n) * epsfac
         if ( q.eq.2 .and. n.eq.k) sca(n) = sca(n) * epsfac
         deltau_vert(n,1) = abs(n) + sca(n)
         omega_total(n,1) = sca(n) /  deltau_vert(n,1)
         asym(n,1)        = asy_save(n)
        enddo

C  2stream call

        CALL TWOSTREAM_MASTER
     +      (DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     +       NLAYERS, 2*NLAYERS, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     +       NBEAMS, N_USER_ANGLES, N_USER_RELAZMS, Ngeoms,
     +       STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     +       BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     +       SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, 1,
     +       ALB, DELTAU_VERT, OMEGA_TOTAL, ASYM,
     +       INTENSITY_TOA_FD, INTENSITY_BOA_FD,
     +       STATUS_INPUTCHECK, STATUS_CALCULATION)

C  write results

        do g = 1, ngeoms
         write(55,34)'uptoa',g,q,k,profilewf_toa(g,k,q,1),
     &   (intensity_toa_fd(g,1)-intensity_toa(g,1))/eps,
     &  100.0d0*((profilewf_toa(g,k,q,1)*eps/
     &     (intensity_toa_fd(g,1)-intensity_toa(g,1)))-1.0d0)
         write(55,34)'dnboa',g,q,k,profilewf_boa(g,k,q,1),
     &   (intensity_boa_fd(g,1)-intensity_boa(g,1))/eps,
     &  100.0d0*((profilewf_boa(g,k,q,1)*eps/
     &     (intensity_boa_fd(g,1)-intensity_boa(g,1)))-1.0d0)
        enddo

C  End FD loops

       enddo
      enddo

C  Continuation point for avoiding profile FD testing

 3456 continue

C  FD testing for COLUMN WFS
C  =========================

C  Skip next section if not doing column WFS

      if ( .not. do_column_wffd ) go to 4567

C  write header

      write(55,'(a//a,a/)')  'Column WFS '//cstream//':',
     &  '      geom#  WF#    ',
     &  '  2s-analytic       2s-Findiff.       %diff'

C  Start FD loop

       do q = 1, 1
    
C  Optical property inputs

        alb(1) = albedo_save
        column = column_save * epsfac
        do n = 1, nlayers
         abs(n) = column * abs_save(n)
         sca(n) = sca_save(n)
         deltau_vert(n,1) = abs(n) + sca(n)
         omega_total(n,1) = sca(n) /  deltau_vert(n,1)
         asym(n,1)        = asy_save(n)
        enddo

C  2stream call

        CALL TWOSTREAM_MASTER
     +     ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     +       NLAYERS, 2*NLAYERS, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     +       NBEAMS, N_USER_ANGLES, N_USER_RELAZMS, Ngeoms,
     +       STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     +       BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     +       SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, 1,
     +       ALB, DELTAU_VERT, OMEGA_TOTAL, ASYM,
     +       INTENSITY_TOA_FD, INTENSITY_BOA_FD,
     +       STATUS_INPUTCHECK, STATUS_CALCULATION)

C  write results

        do g = 1, ngeoms
         write(55,35)'uptoa',g,q,columnwf_toa(g,q,1),
     &   (intensity_toa_fd(g,1)-intensity_toa(g,1))/eps,
     &  100.0d0*((columnwf_toa(g,q,1)*eps/
     &     (intensity_toa_fd(g,1)-intensity_toa(g,1)))-1.0d0)
         write(55,35)'dnboa',g,q,columnwf_boa(g,q,1),
     &   (intensity_boa_fd(g,1)-intensity_boa(g,1))/eps,
     &  100.0d0*((columnwf_boa(g,q,1)*eps/
     &     (intensity_boa_fd(g,1)-intensity_boa(g,1)))-1.0d0)
        enddo

C  End FD loop

      enddo

C  Continuation point for avoiding column FD testing

 4567 continue

C  FD testing for SURFACE WFS
C  ==========================

C  Skip next section if not doing surface WFS

      if ( .not. do_surface_wffd ) go to 5678

C  write header

      write(55,'(/a//a,a/)')  'Surface WFS '//cstream//':',
     &  '      geom#  WF#    ',
     &  '  2s-analytic       2s-Findiff.       %diff'

C  Optical property inputs

      alb(1) = albedo_save * epsfac
      column = column_save
      do n = 1, nlayers
        abs(n) = column * abs_save(n)
        sca(n) = sca_save(n)
        deltau_vert(n,1) = abs(n) + sca(n)
        omega_total(n,1) = sca(n) /  deltau_vert(n,1)
        asym(n,1)        = asy_save(n)
      enddo

C  2stream call

      CALL TWOSTREAM_MASTER
     +     ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     +       NLAYERS, 2*NLAYERS, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     +       NBEAMS, N_USER_ANGLES, N_USER_RELAZMS, Ngeoms,
     +       STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     +       BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     +       SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, 1,
     +       ALB, DELTAU_VERT, OMEGA_TOTAL, ASYM,
     +       INTENSITY_TOA_FD, INTENSITY_BOA_FD,
     +       STATUS_INPUTCHECK, STATUS_CALCULATION)

C  write results

      do g = 1, ngeoms
        write(55,35)'uptoa',g,1,surfacewf_toa(g,1,1),
     &   (intensity_toa_fd(g,1)-intensity_toa(g,1))/eps,
     &  100.0d0*((surfacewf_toa(g,1,1)*eps/
     &     (intensity_toa_fd(g,1)-intensity_toa(g,1)))-1.0d0)
        write(55,35)'dnboa',g,1,surfacewf_boa(g,1,1),
     &   (intensity_boa_fd(g,1)-intensity_boa(g,1))/eps,
     &  100.0d0*((surfacewf_boa(g,1,1)*eps/
     &     (intensity_boa_fd(g,1)-intensity_boa(g,1)))-1.0d0)
      enddo

C  Continuation point for avoiding  SURFACE WFs

 5678 continue

C  Formats

 34   format(a,3i5,2x,1p3e18.9)
 35   format(a,2i5,2x,1p3e18.9)

C  CLose file

      close(55)

C  Finish

      STOP
      END PROGRAM test_2sJ
