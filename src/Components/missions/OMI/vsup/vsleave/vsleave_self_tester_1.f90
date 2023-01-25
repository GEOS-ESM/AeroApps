

      program VSLEAVE_Self_Tester_2p8p5

!  New test of VSLEAVE supplement, Version 2.8
!  --------------------------------------------

!   This is designed to produce multigeometrical output for
!      1. VSLEAVE Exact (direct-bounce) reflectance
!      2. diffuse-field (mutliply-scattered) VSLEAVE Fourier components
      
!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VSLEAVE_SUP_AUX_m, Only : VSLEAVE_READ_ERROR
      USE VSLEAVE_SUP_MOD_m

!  Implicit none

      IMPLICIT NONE

!  VSLEAVE supplement file inputs status structure

      TYPE(VSLEAVE_Input_Exception_Handling)   :: VSLEAVE_Sup_InputStatus

!  VSLEAVE supplement input structures

      TYPE(VSLEAVE_Sup_Inputs)             :: VSLEAVE_Sup_In

!  VSLEAVE supplement output structures

      TYPE(VSLEAVE_Sup_Outputs)            :: VSLEAVE_Sup_Out
      TYPE(VSLEAVE_Output_Exception_Handling)  :: VSLEAVE_Sup_OutputStatus

!  Saved Results
!  =============

!  22 June 2022. Add ESLEAVE_S0

      REAL(fpk) :: ISLEAVE_EX    ( MAXSTOKES, MAXBEAMS )
      REAL(fpk) :: ESLEAVE_EX    ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk) :: ESLEAVE_S0    ( MAXSTOKES, MAXBEAMS )
      REAL(fpk) :: MSLEAVE_SD    ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: MSLEAVE_SU    ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Local Variables
!  ===============

      INTEGER           :: i, ia, ib, j, m, o1, um
      INTEGER           :: nstokes, nbeams, nusers, nazims, ndisos, nmoms
      CHARACTER(LEN=80) :: infile, outfile
      CHARACTER(LEN=99) :: trace
      REAL(fpk)         :: Lw_basic

!  7/7/21. AOS_MODEl control

      INTEGER :: AOS_MODEL

!  7/7/21. ABTOT_CONTROL
!    -- If True,  use R2 = btot/(atot+btot)
!    -- If False, use R2 = btot/atot

      LOGICAL :: ABTOT_CONTROL

!  New inputs
      NAMELIST /VSLEAVE/ VSLEAVE_Sup_In
      AOS_MODEL = 99 
      ABTOT_CONTROL = .false.

!  Define input file

      infile = 'vsleave_self_tester_water.cfg'
!      write(*,*) 'trim(infile) = |',trim(infile),'|'

!  Get the VSLEAVE inputs

      CALL VSLEAVE_INPUTMASTER (   &
        infile,                    & ! Input
        VSLEAVE_Sup_In,            & ! Outputs
        VSLEAVE_Sup_InputStatus )    ! Outputs

      IF ( VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VSLEAVE_READ_ERROR ( 'vsleave_self_tester_water.log', VSLEAVE_Sup_InputStatus )

!  Define some inputs not handled by VSLEAVE_LIN_INPUTMASTER

      VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = 'data'

!  Setup Proxies

      nstokes = VSLEAVE_Sup_In%SL_nstokes
      nbeams  = VSLEAVE_Sup_In%SL_nbeams
      nusers  = VSLEAVE_Sup_In%SL_n_user_streams
      nazims  = VSLEAVE_Sup_In%SL_n_user_relazms
      ndisos  = VSLEAVE_Sup_In%SL_nstreams
      nmoms   = 2*ndisos - 1

!  VSLEAVE baseline call


      CALL VSLEAVE_MAINMASTER ( AOS_MODEL, ABTOT_CONTROL, &
        VSLEAVE_Sup_In,          & ! Inputs
        VSLEAVE_Sup_Out,         & ! Outputs
        VSLEAVE_Sup_OutputStatus ) ! Output Status

      if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
         TRACE = 'VSLEAVE_MAINMASTER failed' ; go to 678
      endif

!  Make Output. 22 June 2022. Add ESLEAVE_S0

      o1 = 1
      Lw_basic =  VSLEAVE_Sup_Out%SL_LW_BASIC(o1,1)
      print *, Lw_basic
      open(87,file="VSLEAVE_input.nml")
      write(87,nml=VSLEAVE)

      stop
      do ib = 1, nbeams
         ISLEAVE_EX(o1,ib) = Lw_basic * VSLEAVE_Sup_Out%SL_KSCALE_ISOTROPIC(o1,ib)
         ESLEAVE_S0(o1,ib) = Lw_basic * VSLEAVE_Sup_Out%SL_KSCALE_SOLARVERT(o1,ib)
         do um = 1, nusers
            do ia = 1, nazims
               ESLEAVE_EX(o1,um,ia,ib) = Lw_basic * VSLEAVE_Sup_Out%SL_KSCALE_USERANGLES(o1,um,ia,ib)
            enddo
         enddo
         do m = 0, nmoms
            do j = 1, ndisos
               MSLEAVE_SD(m,o1,j,ib)  = Lw_basic *VSLEAVE_Sup_Out%SL_KSCALE_F_0(m,o1,j,ib)
            enddo
            do um = 1, nusers
               MSLEAVE_SU(m,o1,um,ib) = Lw_basic *VSLEAVE_Sup_Out%SL_KSCALE_USER_F_0(m,o1,um,ib)
            enddo
         enddo
      enddo

!  Define & open output file

      outfile = 'results_vsleave_self_test_water.all'
      open(36,file=trim(outfile),status='unknown')

!  Write results

      write(36,*)
      write(36,'(A)') '****************************************'
      write(36,'(T9,A)')'Water surface-leaving'
      write(36,'(A)') '****************************************'

      do o1 = 1, nstokes

      write(36,'(/T10,A,i1)')'For Stokes: ',o1

      write(36,'(/T10,A,es15.6)')'LW_BASIC = ', Lw_basic

!  Basic output. 22 June 2022. Add ESLEAVE_S0

      write(36,'(/T10,A)')'STANDARD QUANTITIES'

      write(36,'(/A//T15,A/)')'Isotropic SLEAVE, all geometries','sza     Iso SLEAVE value'
      do ib = 1, nbeams
        write(36,'(A,T15,i2,5x,es15.6)')' ',ib,ISLEAVE_EX(o1,ib)
      enddo

      write(36,'(/A//T15,A/)')'Exact SLEAVE, Direct Line-of-sight','sza vza azm    Exact SLEAVE value'
      do ib = 1, nbeams
         do um = 1, nusers
            do ia = 1, nazims
             write(36,'(A,T15,3(i2,2x),es16.6)')' ',ib,um,ia,ESLEAVE_EX(o1,um,ia,ib)
            enddo
         enddo
      enddo

      write(36,'(/A//T15,A/)')'Exact SLEAVE, Solar-Vertical ','sza    Exact SLEAVE value'
      do ib = 1, nbeams
         write(36,'(A,T15,i2,2x,es16.6)')' ',ib,ESLEAVE_S0(o1,ib)
      enddo

      write(36,'(/A)')'Fourier Components:'
      DO M = 0, NMOMS
         write(36,'(/T15,A,i2//T15,A/)')'Fourier component # ',M,'sza str    SLEAVE value'
         do ib = 1, nbeams
            do j = 1, ndisos
             write(36,'(A,T15,2(i2,2x),es14.6)')'Solars-Stream',ib,j,MSLEAVE_SD(M,o1,j,ib)
            enddo
         enddo

         write(36,'(/T16,A/)')'sza vza    SLEAVE value'
         do ib = 1, nbeams
            do um = 1, nusers
             write(36,'(A,T15,2(i2,2x),es14.6)')'Solars-Users ',ib,um,MSLEAVE_SU(M,o1,um,ib)
            enddo
         enddo
      ENDDO

      enddo

!  Close output file

      close(36)

! 98   format(3I4,1p10e16.7)
! 99   format(I4,1p10e16.7)

      write(*,*)
      write(*,*) 'Main program finished successfully'

      stop

!  Error finish

678   continue

      write(*,*)
      write(*,'(1x,a)') trim(TRACE)
      write(*,*)'Number of error messages from VSLEAVE calculation = ',&
      VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
      do i = 1, VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
       write(*,'(a,i2,a,a)')' Message # ', i, ': ',&
         adjustl(trim(VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)))
      enddo

!  Finish

end program VSLEAVE_Self_Tester_2p8p5

