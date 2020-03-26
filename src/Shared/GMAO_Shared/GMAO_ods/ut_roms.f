
      use m_ods
      use m_roms
      type(ods_vect) :: ods 

       integer only_Ts, rc

      character*255 fname(512), ftype(512)
      character(len=3) :: ch

      print *, 'Enter iopt: '
      read  *, iopt

      nymd = 20020101
      do nhms = 0, 0*60000,  60000

         if ( nhms .eq. 000000 ) ch = '00z'
         if ( nhms .eq. 060000 ) ch = '06z'
         if ( nhms .eq. 120000 ) ch = '12z'
         if ( nhms .eq. 180000 ) ch = '18z'

         call  filenames_ ( ch, iopt )

         print *
         print *, 'nhms = ', nhms
         print *
         call roms_get ( nfiles, fname, nymd, nhms, ftype, ods, rc, only_Ts )
!!!         call ods_tally ( 6, ods, ods%data%nobs, rc )
         print *, 'NOBS = ', ods%data%nobs
         if ( iopt .eq. 1 ) then
            open(10,file='/tmp/obs01.txt',form='formatted')
         else
            open(10,file='/tmp/obs16.txt',form='formatted')
         endif
            do i = 1, ods%data%nobs 
!!!               if ( ods%data%qcexcl(i) .eq. 0 ) then
                   write(10,'(i5,f10.1,i2)') ods%data%kx(i), 
     &                     ods%data%obs(i), ods%data%qcexcl(i) 
!!!               end if
             end do
             close(10)
         call ods_clean ( ods, ier )
      end do

      stop

      CONTAINS

      subroutine filenames_ ( ch, iopt )
      integer iopt
      character(len=3) :: ch
      i = 1

      if ( iopt .eq. 1 ) then
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms'; i=i+1
      nfiles = i - 1
      return
      end if

!
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.1'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.2'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.3'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.4'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.5'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.6'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.7'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.8'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.9'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.10'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.11'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.12'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.13'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.14'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nj.20020101_'//ch//'.oms.15'; i=i+1

!      nfiles = i - 1
!      return

!
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.1'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.2'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.3'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.4'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.5'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.6'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.7'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.8'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.9'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.10'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.11'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.12'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.13'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.14'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nk.20020101_'//ch//'.oms.15'; i=i+1
!
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.1'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.2'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.3'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.4'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.5'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.6'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.7'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.8'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.9'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.10'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.11'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.12'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.13'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.14'; i=i+1
      fname(i) = 'T13r2_16cpuc.level2.tovs_nl.20020101_'//ch//'.oms.15'; i=i+1

      nfiles = i - 1

!      do i = 1, nfiles
!         print *, 'File: ', trim(fname(i))
!      end do

      end subroutine filenames_

      end

