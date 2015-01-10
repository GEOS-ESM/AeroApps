      PROGRAM Siewert2000_validation

C  include file of dimensions and numbers

      INCLUDE '../../includes/VLIDORT.PARS'

C  input variables

      INCLUDE '../../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../../includes/VLIDORT_BOOKKEEP.VARS'

C  result variables

      INCLUDE '../../includes/VLIDORT_RESULTS.VARS'

C  local dimensions

      integer          maxmom,maxlay
      parameter        ( maxmom = 100, maxlay = 20 )

C  help variables

      INTEGER          STATUS_INPUTREAD
      INTEGER          STATUS_INPUTCHECK
      INTEGER          STATUS_CALCULATION

      integer          um,ua,l,n,v,ut,o1
      double precision diff, tautotal, taudiff, topheight

C  when needed

c      double precision s_i,s_q,s_u,s_v,pzn,dpzn, aerosol_tau
c      double precision ds_i,ds_q,ds_u,ds_v,dpzn1,dpzn2,depol
c      INTEGER          ib,um,ua,v,l,n,m,o1,drn,lout,mdum
c      double precision air(22),hdum,pdum(18),agas,sigray,difgas

C  aerosol data problem 2a

      DOUBLE PRECISION PROBLEM_IIA(6,0:11)
      DATA (PROBLEM_IIA(1,L),L = 0,11)/
     &   0.000000D0,  0.000000D0,  3.726079D0,  2.202868D0,
     &   1.190694D0,  0.391203D0,  0.105556D0,  0.020484D0,
     &   0.003097D0,  0.000366D0,  0.000035D0,  0.000003D0 /
      DATA (PROBLEM_IIA(2,L),L = 0,11)/
     &   1.000000D0,  2.104031D0,  2.095158D0,  1.414939D0,
     &   0.703593D0,  0.235001D0,  0.064039D0,  0.012837D0,
     &   0.002010D0,  0.000246D0,  0.000024D0,  0.000002D0 /
      DATA (PROBLEM_IIA(3,L),L = 0,11)/
     &   0.000000D0,  0.000000D0, -0.116688D0, -0.209370D0,
     &  -0.227137D0, -0.144524D0, -0.052640D0, -0.012400D0,
     &  -0.002093D0, -0.000267D0, -0.000027D0, -0.000002D0 /
      DATA (PROBLEM_IIA(4,L),L = 0,11)/
     &   0.915207D0,  2.095727D0,  2.008624D0,  1.436545D0,
     &   0.706244D0,  0.238475D0,  0.056448D0,  0.009703D0,
     &   0.001267D0,  0.000130D0,  0.000011D0,  0.000001D0 /
      DATA (PROBLEM_IIA(5,L),L = 0,11)/
     &   0.000000D0,  0.000000D0,  0.065456D0,  0.221658D0,
     &   0.097752D0,  0.052458D0,  0.009239D0,  0.001411D0,
     &   0.000133D0,  0.000011D0,  0.000001D0,  0.000000D0 /
      DATA (PROBLEM_IIA(6,L),L = 0,11)/
     &   0.000000D0,  0.000000D0,  3.615946D0,  2.240516D0,
     &   1.139473D0,  0.365605D0,  0.082779D0,  0.013649D0,
     &   0.001721D0,  0.000172D0,  0.000014D0,  0.000001D0 /

C  initialize Greek matrix

      DO N = 1, MAXLAYERS
        DO L = 0, MAXMOMENTS
          DO O1 = 1, NSTOKES_SQ
            GREEKMAT_TOTAL_INPUT(L,N,O1) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Test 2
C  ======

C  description
C     Problem IIA Slab, Siewert (2000)
C     Single layer of optical thickness 1, Lambertian albedo 0.0
C     Aerosol scattering problem IIA, SSA = 0.973257
C     plane parallel, classical solution
C     1 solar angles, 11 viewing angles, 3 azimuths (33 geometries)
C     Upwelling/downwelling radiance at 7 Tau values
C     WITH POLARIZATION

      CALL VLIDORT_INPUT_MASTER
     &        ( 'Siewert2000_validation.inp',
     &          'Siewert2000_validation.log',status_inputread)
      if ( status_inputread.ne.vlidort_success)stop'badread input data'

C  Siewert value

      SZANGLES(1) = DACOS(0.6d0)/DEG_TO_RAD

C  Set up of optical properties

      topheight = 10.0d0
      tautotal  = 1.0d0
      diff = dble(nlayers)
      height_grid(0) = topheight
      taudiff = tautotal / diff
      DO n = 1, nlayers
       height_grid(n) = topheight - dble(n)*topheight/diff
       deltau_vert_input(n) = taudiff
       omega_total_input(n) = 0.973527d0
       DO L = 0, NGREEK_MOMENTS_INPUT
        GREEKMAT_TOTAL_INPUT(L,n,1)  = PROBLEM_IIA(2,L)
        IF ( L.EQ.0)GREEKMAT_TOTAL_INPUT(L,n,1) = 1.0D0
        IF ( NSTOKES .GT. 1 ) THEN
          GREEKMAT_TOTAL_INPUT(L,n,2)  =  PROBLEM_IIA(3,L)
          GREEKMAT_TOTAL_INPUT(L,n,5)  =  PROBLEM_IIA(3,L)
          GREEKMAT_TOTAL_INPUT(L,n,6)  =  PROBLEM_IIA(1,L)
          GREEKMAT_TOTAL_INPUT(L,n,11) =  PROBLEM_IIA(6,L)
          GREEKMAT_TOTAL_INPUT(L,n,12) = - PROBLEM_IIA(5,L)
          GREEKMAT_TOTAL_INPUT(L,n,15) =   PROBLEM_IIA(5,L)
          GREEKMAT_TOTAL_INPUT(L,n,16) =   PROBLEM_IIA(4,L)
        ENDIF
       ENDDO
      ENDDO

C  Call to VLIDORT_MASTER

       write(*,'(/a)')'Siewert2000 validation status:'
       CALL VLIDORT_MASTER_LAMBERTIAN
     &        ( STATUS_INPUTCHECK, STATUS_CALCULATION )
       CALL VLIDORT_STATUS ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  output

       IF ( STATUS_CALCULATION .NE. VLIDORT_SUCCESS .OR. 
     &     STATUS_INPUTCHECK   .EQ. VLIDORT_SERIOUS ) THEN
        WRITE(*,*)'Failure: No output written to file'
       ELSE
        WRITE(*,*)
     &      'Success: output written to file Siewert2000_validation.out'
        OPEN(1,FILE='Siewert2000_validation.out',STATUS='unknown')

        UA = 1
        write(1,*)'Table 2'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        write(1,*)'Table 3'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
           V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO

        UA = 2
        write(1,*)'Table 4'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        write(1,*)'Table 5'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
           V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        write(1,*)'Table 6'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,3,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,3,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        write(1,*)'Table 7'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,4,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
           V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,4,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO

        UA = 3
        write(1,*)'Table 8'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
          V = VZA_OFFSETS(1,UM) + UA
          write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,1,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        write(1,*)'Table 9'
        DO UM = 1, N_USER_VZANGLES
          V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')-USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,UPIDX),UT = 1, N_USER_LEVELS)
        ENDDO
        DO UM = N_USER_VZANGLES, 1, -1
           V = VZA_OFFSETS(1,UM) + UA
           write(1,'(f8.2,1p7e15.5)')USER_STREAMS(UM),
     &        (PIE*STOKES(UT,V,2,DNIDX),UT = 1, N_USER_LEVELS)
        ENDDO

        CLOSE(1)
       ENDIF

C  Stop

       stop
       end
