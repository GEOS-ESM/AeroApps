!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ReadDel --- Reads innovation data from del-file.
! 
! !INTERFACE:
!
      subroutine Read_Del ( fname, nsto, iopt,
     &                      lat, lon, lev, del, kx, kt, qc, nobs,
     &                      nymd, nhms, ier )


! !INPUT PARAMETERS: 
!
      use m_stdio, only : stderr
      use m_stdio, only : stdout

      use m_ioutil, only : luavail

      implicit NONE
      character*(*)  fname         ! del file name
      integer        nsto          ! space for data vectors as
                                   !  declared in calling program
      integer        iopt          ! Option:
                                   !   0 - all observations on file
                                   !   1 - only obs which passed QC
!
! !OUTPUT PARAMETERS:
!
      real           lat(nsto)     ! latitude (degrees)
      real           lon(nsto)     ! longitude (degrees)
      real           lev(nsto)     ! level (hPa)
      real           del(nsto)     ! inovation
      integer        kx(nsto)      ! data source index
      integer        kt(nsto)      ! data type   index
      real           qc(nsto)      ! GEOS-1 Quality Control flag
      integer        nobs          ! actual number of observations

      integer        nymd          ! year-month-day, e.g., 970301
      integer        nhms          ! hour-min-sec, e.g., 120000
      integer        ier           ! Error return code
                                   !   0  successfull
                                   !   1  end-of-file reached
                                   !   2  error reading file
                                   !   3  not enough work space
        real zz
!
! !DESCRIPTION: 
!
!            Reads next synoptic time on file. Returns surface, moisture
!  and upper air innovations. User can select only those observations which 
!  passed quality control. All data that passed QC are stored in an 1-D array
!  with other proper information are stored in several other arrays. 
!
!
! !REVISION HISTORY: 
!
!  06Oct96   da Silva/Lamich   Initial code. Derived from Read_Hobs().
!  17Oct97   G. P. Lou         Modifications.
!  15Mar98   da Silva          Minor touch ups (comments only)
!  10Jun04   Todling           Using m_ioutil
!
!EOP
!-------------------------------------------------------------------------
!BOC

      integer     OK
      parameter ( OK = 0 )
      real        GOOD,      VERY_GOOD
      parameter ( GOOD = 1., VERY_GOOD = 3. )

      integer EOF, ERR, STORAGE
      parameter ( EOF=1, ERR=2, STORAGE=3 )

      real AMISS 
      parameter ( AMISS = 1.E15 )

      logical let_u, let_p, let_z, let_q
      
      character      quant*3, Equant*3

      integer        nob, ios
      integer        i, j, nrobs, nsobs, lf
      real           rqc, rKx, u, v, z, d, rmix, slu, slv, slp
      real           ylat, xlon, p

      save            previous_fname
      character*(255) previous_fname

      save     lu
      integer  lu

      data lu / -1 /
      data previous_fname / '!@#$%^&**()_+|' /   ! just garbage


!     Initialize counters
!     -------------------
      nrobs = 0
      nsobs = 0
      nobs = 0

!     Open innovation file
!     --------------------
      lf = len(fname)
      if ( fname .ne. previous_fname(1:lf) ) then
         if ( lu .ne. -1 ) close(lu) 
         lu = luavail()
         open(lu,
     &        file=fname,
     &        form='unformatted',
     &        status='old',
     &        iostat=ios)
         if ( ios .ne. OK ) then
            write(stdout,*) 'ReadDel: Cannot read file ', fname
            ier = ERR
            return
         else
            write(stdout,*) 
            write(stdout,*) 'ReadDel: Just opened file ', fname(1:lf)
         end if
         previous_fname = fname
      else
            write(stdout,*) 'ReadDel: Re-using already opened file ', 
     &                      fname(1:lf)
      end if

      nobs = 0


!     Read header for surface observations
!     ------------------------------------
      Equant = 'SLP'
      read(lu,end=999,err=998) nob, nymd, nhms, quant         
      if ( quant .eq. 'MIX' ) go to 2000
      if ( quant .eq. 'ZUV' ) go to 3000
      if ( quant .ne. Equant ) then
           write(stderr,*) 'ReadDel: quantity ', Equant, 
     &                     ' expected but ', 
     &                     quant, ' found on file ', fname
           ier = ERR
           return
      end if

!     Read all sfc obs for this synoptic time
!     ---------------------------------------
 1000 do 10 i = 1, nob
         read(lu,iostat=ios) ylat, xlon, rqc, p, rKx, 
     &        slu, slv, slp
         if ( ios .ne. OK ) then
            print *, 'ReadDel: problems with input file, slp '
            ier = ERR
            return
         end if

!        Save the report
!        ---------------
         if ( nobs+3 .gt. nsto ) then
            write(stderr,*) 'ReadDel: not enough space to read obs'
            write(stderr,*) 'ReadDel: nobs = ', nobs+3
            write(stderr,*) 'ReadDel: nsto = ', nsto
            ier = STORAGE
            return
         end if

!        Filter obs not used by analysis if iopt=1
!        -----------------------------------------
         let_u = .true.
         let_p = .true.
         if ( iopt .eq. 1 ) then
            let_u = .false.
            let_p = .false.
            if ( nint(rqc) .eq. 1 ) then
               let_p = .true.
            else if ( nint(rqc) .eq. 2 ) then 
               let_u = .true.
            else if ( nint(rqc) .eq. 3 ) then 
               let_u = .true.
               let_p = .true.
            end if
         end if               

         do j = nobs+1, nobs+3
            lat(j) = ylat
            lon(j) = xlon
            lev(j) = 2000.
            kx(j)  = nint(rKx)
            qc(j)  = rqc
         end do

         j = nobs
         if ( let_u ) then
            if ( abs(slu-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 1
               del(j) = slu
            end if
            
            if ( abs(slv-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 2
               del(j) = slv
            end if
         end if

         if ( let_p ) then
            if ( abs(slp-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 3
               del(j) = slp
            end if
         end if


         nobs = j           ! update obs counter

 10   continue



!     Read header of mixing ratio reports
!     -----------------------------------
      Equant = 'MIX'
      read(lu,iostat=ios) nob, nymd, nhms, quant         
      if ( ios .ne. OK ) then
         print *, 'ReadDel: problems with input file, mix 1'
         call exit(7)
      end if
      if ( quant .eq. 'ZUV' ) go to 3000
      if ( quant .ne. Equant ) then
         write(stderr,*) 'ReadDel: quantity ', Equant, 
     &        ' expected but ', 
     &        quant, ' found.'
         ier = ERR
         return
      end if

!     Read all mixing ratio reports
!     -----------------------------      
 2000 do 20 i = 1, nob

         read(lu,iostat=ios) ylat, xlon, rqc, p, rKx, 
     &        d, d, rmix
         if ( ios .ne. OK ) then
            print *, 'ReadDel: problems with input file, mix 2'
            ier = ERR
            return
         end if


!        Filter obs not used by analysis if iopt=1
!        -----------------------------------------
         let_q = .true.
         if ( iopt .eq. 1 .and. nint(rqc) .ne. 1 ) let_q = .false.

!        Save the report
!        ---------------
         if ( nobs+1 .gt. nsto ) then
            write(stderr,*) 'ReadDel: not enough space to read obs'
            write(stderr,*) 'ReadDel: nobs = ', nobs+1
            write(stderr,*) 'ReadDel: nsto = ', nsto
            ier = STORAGE
            return
         end if
        
         if ( let_q ) then
            if ( abs(rmix-amiss)/amiss .gt. 0.001 ) then
               lat(nobs+1) = ylat
               lon(nobs+1) = xlon
               lev(nobs+1) = p
               kx(nobs+1)  = nint(rKx)
               qc(nobs+1)  = rqc
               kt(nobs+1)  = 7
               del(nobs+1) = rmix
               nobs = nobs + 1  ! update obs counter
            end if
         end if

 20   continue


!     Read (u,v,z)
!     ------------
      Equant = 'ZUV'
      read(lu,iostat=ios) nob, nymd, nhms, quant         
      if ( ios .ne. OK ) then
         print *, 'ReadDel: problems with input file, zuv 1'
         ier = ERR
         return
      end if
      if ( quant .ne. Equant ) then
         write(stderr,*) 'ReadDel: quantity ', Equant, 
     &        ' expected but ', 
     &        quant, ' found.'
         ier = ERR
         return
      end if
      
      
 3000 do 30 i = 1, nob
         
!        Read the upper level report
!        ---------------------------
         read(lu,iostat=ios) ylat, xlon, rqc, p, rKx, u, v, z
         if ( ios .ne. OK ) then
            print *, 'ReadDel: problems with input file, zuv 2'
            ier = ERR
            return
         end if
         
!        Filter obs not used by analysis if iopt=1
!        -----------------------------------------
         let_u = .true.
         let_z = .true.
         if ( iopt .eq. 1 ) then
            let_u = .false.
            let_z = .false.
            if ( nint(rqc) .eq. 1 ) then
               let_z = .true.
            else if ( nint(rqc) .eq. 2 ) then 
               let_u = .true.
            else if ( nint(rqc) .eq. 3 ) then
               let_z = .true.
               let_u = .true.
            end if
         end if               


!        Save the report
!        ---------------
         if ( nobs+3 .gt. nsto ) then
            write(stderr,*) 'ReadDel: not enough space to read obs'
            write(stderr,*) 'ReadDel: nobs = ', nobs+3
            write(stderr,*) 'ReadDel: nsto = ', nsto
            ier = STORAGE
            return
         end if
         
         do j = 1, 3
            lat(nobs+j) = ylat
            lon(nobs+j) = xlon
            lev(nobs+j) = p
            kx(nobs+j)  = nint(rKx)
            qc(nobs+j)  = rqc
         end do

!    missing data are evaluated 
         j = nobs
         if ( let_u ) then
            if ( abs(u-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 4
               del(j) = u
            end if
            if ( abs(v-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 5
               del(j) = v
            end if
         end if

         if ( let_z ) then
            if ( abs(z-amiss)/amiss .gt. 0.001 ) then
               j = j + 1
               kt(j) = 6
               del(j) = z
               zz = abs(z - nint(z) )
            end if
         end if
         
         nobs = j          ! update obs counter
         
 30   continue
      

!     All is well
!     -----------
      ier = OK
      return


!     Error
!     -----
 998  continue
      ier = ERR
      return

!     End of file reached
!     -------------------
 999  continue
      ier = EOF
      return

      end
                                  
!EOC
