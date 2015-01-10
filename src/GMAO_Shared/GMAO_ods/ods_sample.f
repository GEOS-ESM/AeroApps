
*====================================================================
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  ods_sample -- Main Driver
!
! !DESCRIPTION: This routine reads a list of GEOS/DAS innovation files
!               in ODS format.
!
! !INTERFACE:
      program ods_sample
!
! !REVISION HISTORY: 
!
!     12Dec2000  C. Redder  Original code.
!
! EOP
!-------------------------------------------------------------------------

      use m_ODS_Sample
      use m_ods,         only : ods_vect,
     .                          ODS_Get,
     .                          ODS_Put
      use m_die,         only : die
      use m_stdio,       only : stdout
      use m_StrTemplate, only : StrTemplate

      implicit NONE
      integer, external     :: ODS_CalDat
      integer, external     :: ODS_Julian

      type ( ods_vect )     :: Obs
      type ( ods_vect ), pointer, dimension (:) :: Binned_Obs

      integer               :: iret
      integer               :: NSyn1,  NSyn2,   iSyn,    NAllSyn
      integer               :: JHr1,   JHrBeg1, JHrEnd1, DJHr1
      integer               :: JHr2,   JHrBeg2, JHrEnd2, DJHr2
      integer               :: nymd1,  nhms1
      integer               :: nymd2,  nhms2
      real,   dimension (2) :: f
      integer               :: NFiles, iFile
      integer, parameter    :: NFiles_Max = 4096
      character (len = 255), dimension ( NFiles_Max )
     .                      :: InFiles
      character (len = 255) :: InFile, OutFile
      character (len = 255) :: TPlate
      character (len = 255) :: ExpID
      character (len = 255) :: Type
      character (len = *), parameter :: MyName = 'ODS_Sample::Main'
      integer, parameter    :: NLines_Max   = 4
      character (len = 100), dimension ( NLines_Max )
     .                      :: ErrorMessage

!     Initialize package
!     ------------------
      nullify ( Binned_Obs )
      call Init ( ExpID, NSyn2, f, TPlate, NFiles, InFiles ) 

!     For each input file ...
!     -----------------------
      do iFile = 1, NFiles

!        ... read relevant input file parameters
!        ---------------------------------------
         InFile = InFiles ( iFile )
         call Inquire_ODS ( InFile, NSyn1, JHrBeg1, JHrEnd1 )

!        ... and for each synoptic time in the input file ...
!        ----------------------------------------------------
         DJHr1  = 24 / NSyn1
         do JHr1 = JHrBeg1, JHrEnd1, DJHr1
            nymd1 = ODS_CalDat ( JHr1 / 24 + 1 )
            nhms1 =        mod ( JHr1,  24 ) * 10000

!           ... read the observation data
!           -----------------------------
            call ODS_Get ( InFile, nymd1, nhms1, Type, Obs, iret )
            if ( iret .ne. 0 ) then
               write ( ErrorMessage, 901 ) trim(InFile), iret
               call die ( MyName, ErrorMessage (:1) )

            end if

            if ( Obs%data%nobs<1 ) cycle

!           ... determine the synoptic times for the output file
!           ----------------------------------------------------
            call Set_SynTimes ( JHr1,  Obs % Data % Time, f,
     .                          NSyn2, JHrBeg2, JHrEnd2 )

!           Allocate space for the binned observations ...
!           ----------------------------------------------
            NAllSyn = ( JHrEnd2 - JHrBeg2 ) * NSyn2 / 24 + 1
            if ( associated ( Binned_Obs ) ) deallocate ( Binned_Obs )
            allocate ( Binned_Obs ( NAllSyn ), stat = iret )
            if ( iret .ne. 0 ) then
               write ( ErrorMessage, 902 ) 'Binned_Obs', 
     .                                      trim ( ItoA ( NAllSyn )),
     .                                      trim ( ItoA ( iret    ))
               call die ( MyName, ErrorMessage ( :1 ) )

            end if

!           ... bin and sort the observation vectors
!               according to synoptic index time.
!           ----------------------------------------
            call Bin_Obs      ( JHr1, f, NSyn2, JHrBeg2, Type,
     .                          Obs,  Binned_Obs )

!           ... and for each synoptic time in the output file
!           -------------------------------------------------
            do iSyn = 1, NAllSyn
               JHr2  = JHrBeg2 + ( iSyn - 1 ) * 24 / NSyn2
               nymd2 = ODS_CalDat ( JHr2 / 24 + 1 )
               nhms2 =        mod ( JHr2,  24 ) * 10000

!              ... reset the time attribute to be consistent
!                  with the synoptic times for the output file
!              -----------------------------------------------
               call NewTime ( JHr1, JHr2,
     .                        Binned_Obs ( iSyn ) % Data % Time )

!              ... generate the name of the output file
!              ----------------------------------------
               call StrTemplate ( OutFile, TPlate,
     .                            xid  = ExpID,
     .                            nymd = nymd2,
     .                            nhms = nhms2 )

!              ... write the results to the output file
!              ----------------------------------------
               call ODS_Put     ( OutFile, Type,
     .                            nymd2,   nhms2,
     .                            Binned_Obs ( iSyn ),
     .                            iret,
     .                            append = .true. )
               if ( iret .ne. 0 ) then
                  write ( ErrorMessage, 903 ) trim (OutFile), iret
                  call die ( MyName, ErrorMessage (:1) )

               end if

            end do
         end do
         write ( stdout, 801 ) MyName, InFile

      end do

      call exit (0)
 801  format ( /, a, ': Processed input file, ', a )
 901  format ( ' Error in reading the ODS input file, ',  a,
     .          '; iret = ', i2 )
 902  format ( ' Error in allocating the array ', a, '(', a,
     .          '); iret = ', a )
 903  format ( ' Error in writing to the ODS output file, ', a,
     .          '; iret = ', i2 )

      end program ods_sample
*====================================================================
