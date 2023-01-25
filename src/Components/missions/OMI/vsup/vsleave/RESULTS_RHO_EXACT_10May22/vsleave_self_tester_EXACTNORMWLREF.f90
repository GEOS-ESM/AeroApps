

      program VSLEAVE_Self_Tester_ExactNormWLRef

!  6 May 2022. Computation of Exact Normalized Water-leaving reflectances
!    -- New formalism.

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

      integer, parameter :: maxwavs = 501
      REAL(fpk) :: WLExactNorm ( 3, maxwavs ), WLUnNorm ( 3, maxwavs ), wavs ( maxwavs ), savconc(3), cos50, res, start
      data savconc / 0.05d0, 0.5d0, 5.0d0 /

!  Local Variables
!  ===============

      INTEGER           :: w, ic, nwavs, i
      CHARACTER(LEN=80) :: infile, outfile
      CHARACTER(LEN=99) :: trace

!  7/7/21. AOS_MODEL control. ABTOT_CONTROL
!    -- If True,  use R2 = btot/(atot+btot)
!    -- If False, use R2 = btot/atot

      INTEGER :: AOS_MODEL     = 99
!      LOGICAL :: ABTOT_CONTROL = .true.
      LOGICAL :: ABTOT_CONTROL = .false.

!  Define input file

      infile = 'vsleave_self_tester_ExactNormWLRef_503090.cfg'

!  Get the VSLEAVE inputs

      CALL VSLEAVE_INPUTMASTER (   &
        infile,                    & ! Input
        VSLEAVE_Sup_In,            & ! Outputs
        VSLEAVE_Sup_InputStatus )    ! Outputs

      IF ( VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VSLEAVE_READ_ERROR ( 'vsleave_self_tester_ExactNormWLRef_503090.log', VSLEAVE_Sup_InputStatus )

!  Define some inputs not handled by VSLEAVE_LIN_INPUTMASTER

      VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = 'data'

!  3 chlorophyll loops and wavelengths every nm from 300-800 nm

      nwavs = 501 ; res = 1.0d0 ; start = 300.0d0
      do ic = 1, 3
         VSLEAVE_Sup_In%SL_CHLORCONC = savconc(ic)
write(*,*)'Doing ic = ',ic
         do w = 1, nwavs
            wavs(w) = start + res * dble(w-1)
            VSLEAVE_Sup_In%SL_WAVELENGTH = wavs(w)  * 0.001d0 
if ( mod(w,50).eq.0) write(*,*)'Doing wav # = ',w


!  VSLEAVE baseline call

            CALL VSLEAVE_MAINMASTER ( AOS_MODEL, ABTOT_CONTROL, &
              VSLEAVE_Sup_In,          & ! Inputs
              VSLEAVE_Sup_Out,         & ! Outputs
              VSLEAVE_Sup_OutputStatus ) ! Output Status

            if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
               TRACE = 'VSLEAVE_MAINMASTER failed' ; go to 678
            endif

!  save results

            cos50 = 0.6427876d0
            WLExactNorm(ic,w) = pie * VSLEAVE_Sup_Out%SL_LW_BASIC(1,1)
            WLUnNorm(ic,w)    = pie * VSLEAVE_Sup_Out%SL_LW_BASIC(1,1) * VSLEAVE_Sup_Out%SL_KSCALE_USERANGLES(1,1,1,1) / cos50

!  end doloops

         enddo
      enddo

!  Define & open output file

      outfile = 'EXACTNORM_WLREF.out'
      open(36,file=trim(outfile),status='unknown')
      do w = 1, nwavs
         write(36,98)w,wavs(w),WLExactNorm(1:3,w),WLUnNorm(1:3,w)
      enddo
      close(36)
 98   format(I5,f10.4,2(2x,1p3e16.7))

!  Normal finish

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

end program VSLEAVE_Self_Tester_ExactNormWLRef

