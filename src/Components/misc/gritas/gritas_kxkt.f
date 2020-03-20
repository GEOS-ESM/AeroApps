!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_KxKt --- Set kx-kt table from user defined list
! 
! !INTERFACE:
!
      subroutine GrITAS_KxKt ( KTMAX, KXMAX, L2D_MAX, L3D_MAX,
     &                         listsz, kt_list, kx_list, 
     &                         kxkt_table, list_2d, list_3d, l2d, l3d )
!
! !INPUT PARAMETERS: 
!
      implicit NONE
      integer         KTMAX                    ! max number of kt's
      integer         KXMAX                    ! max number of kx's
      integer         L2D_MAX                  ! max number of 2d grids
      integer         L3D_MAX                  ! max number of 3d grids
      integer         listsz                   ! actual list size 
      integer         kt_list(listsz)          ! data type list
      integer         kx_list(KXMAX,listsz)    ! data source list
!
! !OUTPUT PARAMETERS:
!
      integer         kxkt_table(KXMAX,KTMAX)  ! kx-kt table
      integer         list_2d(l2d_max)         ! corresponding list item for
                                               !  2d grids
      integer         list_3d(l3d_max)         ! corresponding list item for
                                               !  3d grids
      integer         l2d                      ! actual number of 2d grids
      integer         l3d                      ! actual number of 3d grids

! !DESCRIPTION: From the list of kt/kx's for each grid, this routine
!               constructs a kx-kt table which contains the grid number
!  for a pair of (kx,kt). These grid numbers correspond to a 2d or 3d
!  grid depending on the kt value:
! \bv
!       kt = 1, 2, 3, ...  -  surface   (2d) grids; see is_surface()
!       kt = 4, 5, 6, 7    -  upper-air (3d) grids
! \ev
!  
!  kt refers to variable type. For example, when kt is 1, 2, and 3,
!  they represent surface measurements of u, v, and p, respectively.
!  When kt is 4, 5, 6, and 7, they are upper air u, v, h, and q (mixing
!  ratio.)
!
!  kx refers to instrument type. There are many different instruments
!  that will measure same variable type. For example, when kx are
!  3 through 6, it represent ship measurements; when kx are 7 through
!  13, it represent rawinsondes, so on and so forth.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!           GPLou      Modifications
!  14Jul99  da Silva   Added land-surface kt's
!  28May09  Redder     Changed code to define negative value in kt_list
!                      and kxkt_table as indices for surface variables.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      integer i, j, kx, kt, l
      logical is_surface

!     Initialize table, etc.
!     ----------------------
      l2d = 0
      l3d = 0
      do kt = 1, KTMAX
         do kx = 1, KXMAX
            kxkt_table(kx,kt) = 0
         end do
      end do

!     Loop over list and set kx-kt table
!     ----------------------------------
      do j = 1, listsz
         kt = abs ( kt_list(j))
         if ( kt .gt. KTMAX ) 
     &        call gritas_perror ('KxKt: invalid kt - too large' ) 
         if ( kt .lt. 1 ) 
     &        call gritas_perror ('KxKt: invalid kt - negative' )
         if ( is_surface(kt_list(j)) ) then  ! surface variables 
              l2d = l2d + 1
              if ( l2d .gt. L2D_MAX )
     &           call gritas_perror ('KxKt: invalid l2d - too large' )
              l = -1 * l2d          ! set to negative to denote surface var
              list_2d(l2d) = j      ! for easily getting to variable names
         else                       ! upper aur
              l3d = l3d + 1 
              if ( l3d .gt. L3D_MAX )
     &           call gritas_perror ('KxKt: invalid l3d - too large' )
              l = l3d 
              list_3d(l3d) = j ! for easily getting to variable names
         end if
         do i = 1, KXMAX
            kx = kx_list(i,j) 
            if ( kx .lt. 0 ) go to 11
            if ( kx .gt. KXMAX ) 
     &           call gritas_perror ('KxKt: invalid kx - too large' ) 
            kxkt_table(kx,kt) = l                
         end do
 11      continue
      end do

      print *
      print *, 'Number of 2D SURFACE   grids: ', l2d, '    ==> ',
     &          ( list_2d(i), i = 1, l2d )
      print *, 'Number of 3D UPPER-AIR grids: ', l3d, '    ==> ',
     &          ( list_3d(i), i = 1, l3d )
      print *, '                                ---'
      print *, '                              ', listsz
      print *

      end
!EOC



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  is_surface --- determines whether this is a surface kt
! 
! !INTERFACE:

      logical function is_surface ( kt )

      use m_odsmeta, only : ktSurfAll
!
! !INPUT PARAMETERS: 
!
      integer kt   ! data type index

!
! !DESCRIPTION: Examines data type index {\tt kt} and determines
!               whether it is a surface data type. Indicies less than
!               zero are defined as surface variables.
!
! !REVISION HISTORY: 
!
!  14jul99   da Silva   First version
!  24jun04   Todling    Generalized.
!  28may09   Redder     Added case when kt < 0 which is defined as a 
!                       surface variable.
!
!EOP
!-------------------------------------------------------------------------

      integer ic(1)

      is_surface = .false.

      ic(1) = count(ktSurfAll.eq.kt)
      if ( ic(1) > 0 .or. kt < 0 ) is_surface = .true.

      return
      end







