!-------------------------------------------------------------------------
!          NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Output --- Write del grids to IEEE GrADS file
! 
! !INTERFACE:
!
       subroutine GrADS_Output  ( grid_ofname, 
     &                            var_list, var_descr, listsz,
     &                            l2d, list_2d, bias_2d, stdv_2d, nobs_2d,
     &                            l3d, list_3d, bias_3d, stdv_3d, nobs_3d,
     &                            nstat )
!
! !USES:
!
      use m_gritas_grids, only : im, jm, km
      use m_gritas_grids, only : AMISS
      use m_gritas_grids, only : LonMin, dLon
      use m_gritas_grids, only : LatMin, dLat
      use m_gritas_grids, only : nlevs
      use m_gritas_grids, only : plevs

      use m_ioutil, only : luavail

      implicit NONE
!
! !INPUT PARAMETERS: 
!

      character(len=*)  grid_ofname      ! gridded O-F (del) file

      integer       listsz               ! size of the list (< lmax)
      character*11  var_list(listsz)     ! variable name for output 
      character*256 var_descr(listsz)    ! ... and its description
      character*11  var_num     ! variable name of number of obs for output 

      integer         l2d          ! actual number of 2d grids
      integer         l3d          ! actual number of 3d grids

      integer         list_2d(l2d) ! list item for 2d grids
      integer         list_3d(l3d) ! list item for 3d grids

                                   ! Surface (2D) grids:
      real  bias_2d(im,jm,l2d)     !   time mean
      real  stdv_2d(im,jm,l2d)     !   stdv
      real  nobs_2d(im,jm,l2d)     !   number of obs per grid box

                                   ! Upper-air (3D) grids:
      real  bias_3d(im,jm,km,l3d)  !   time mean
      real  stdv_3d(im,jm,km,l3d)  !   stdv
      real  nobs_3d(im,jm,km,l3d)  !   number of obs per grid box
      integer nstat                ! No of stats to write to file:
                                   ! nstat = 1   ...  bias
                                   ! nstat = 2   ...  bias, rms
                                   ! nstat = 3   ...  bias, rms, nobs

!
! !DESCRIPTION: Write del-grids to file.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  17OCt97  G. P. Lou  Modifications.
!  26Mar98  da Silva   Added nstat parameter for GrISAS's sake.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      integer  i, k, l, lu, ldescr


      if ( nstat .lt. 1 .or. nstat .gt. 3 ) then
         call gritas_perror ( 'GrITAS_Output: invalid nstat' )
      endif

!     Open descriptor file
!     --------------------
      lu = luavail()
      open(lu,file=trim(grid_ofname)//'.ctl',
     &status='unknown',form='formatted' )

      print *
      print *, ' [] Writing descriptor file ', trim(grid_ofname)//'.ctl'

      write(lu,61) trim(grid_ofname)
      write(lu,62) 
      write(lu,63)  AMISS
      write(lu,64) im, LonMin, dLon
      write(lu,65) jm, LatMin, dLat
      write(lu,66) nlevs, (plevs(k), k=1, nlevs)
      write(lu,67) nstat
      write(lu,68) listsz
      do i = 1, l2d
         l      = list_2d(i)
         ldescr = max ( len_trim (var_descr(l)),1)
         write(lu,69) var_list(l), var_descr(l)(:ldescr)
      end do
      do i = 1, l3d
         l      = list_3d(i)
         ldescr = max ( len_trim (var_descr(l)),1)
         write(lu,70) var_list(l), nlevs,
     &                var_descr(l)(:ldescr)
      end do
      write(lu,71) 
61    format('DSET ^' ,100a)
62    format('OPTIONS sequential big_endian')
63    format('UNDEF ', g15.10)
64    format('XDEF ',i4,' LINEAR ', 2f9.2)
65    format('YDEF ',i4,' LINEAR ', 2f9.2)
66    format('ZDEF ',i4,' LEVELS ', 4f10.2,/,6(10x,6f10.2/))
67    format('TDEF ',i4,' LINEAR 01DEC1900 1mon')
68    format('VARS ',i5)
69    format(a, '  0 99 ',  a)
70    format(a, i3, ' 99 ', a)
71    format('ENDVARS')

      close(lu)


!     Open output data file
!     ---------------------
      lu = luavail()
      open(lu,file=trim(grid_ofname),
     &status='unknown',form='unformatted')

      print *, ' [] Writing data file ', trim(grid_ofname)
      print *

!     Means (uses GrADS t=1 slot)
!     ---------------------------
      do l = 1, l2d
         call write_2d ( lu, bias_2d(1,1,l), im * jm )
      end do
      do l = 1, l3d
         do k = 1, nlevs
            call write_2d ( lu, bias_3d(1,1,k,l), im * jm )
         end do
      end do
      
!     Stdv  (uses GrADS t=2 slot)
!     ---------------------------
      if ( nstat .ge. 2 ) then
      do l = 1, l2d
         call write_2d ( lu, stdv_2d(1,1,l), im * jm )
      end do
      do l = 1, l3d
         do k = 1, nlevs
            call write_2d ( lu, stdv_3d(1,1,k,l), im * jm )
         end do
      end do
      end if
      
!     Nobs  (uses GrADS t=3 slot)
!     ---------------------------
      if ( nstat .ge. 3 ) then
      do l = 1, l2d
         call write_2d ( lu, nobs_2d(1,1,l), im * jm )
      end do
      do l = 1, l3d
         do k = 1, nlevs
            call write_2d ( lu, nobs_3d(1,1,k,l), im * jm )
         end do
      end do
      end if

      close(lu)

      end

!...............................................................


      subroutine write_2d ( lu, a, n )
      implicit NONE
      integer lu, n
      real a(n)
      real*4 a4(n)    ! _RT: added since all compiles r8 now
      a4 = a
      write(lu) a4
      end
!EOC


