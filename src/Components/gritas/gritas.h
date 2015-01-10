!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !INCLUDE: gritas_dim.h --- GRITAS internal parameters
! 
! !DESCRIPTION: 
!
!  Defines internal parameters for GrITAS.
!
! !DEFINED PARAMETERS: 
!

!     Maximum number of obs per synoptic time
!     Note: only needed in support of "del files"
!     -------------------------------------------
      integer     nobsmax
      parameter ( nobsmax = 700 * 1000 )

!     Number of grids
!     ---------------
      integer     l2d_max                    ! surface grids       
      parameter ( l2d_max = 20 )
      integer     l3d_max                    ! upper-air grids       
      parameter ( l3d_max = 20 )
      integer     lmax                       ! sfc + upper-air
      parameter ( lmax = l2d_max + l3d_max )


!     Max number of observation files
!     -------------------------------
      integer          NDEL_MAX
      parameter       (NDEL_MAX = 1024)

      real, parameter :: eps = 1.e-10

!
! !SEE ALSO: 
!
!  gritas.f        -   GrITAS driver.
!  GrITAS_Init()   -   GrITAS initializer; the vertical pressure levels
!                      are set by this routine.
!  GrITAS_Grids()  -   Initialize grid parameters.
!
! !REVISION HISTORY: 
!
!  05sep97   da Silva   Original code.
!  25mar99   da Silva   Increased ktmax/kxmax
!  09Jun04   Todling    Cleaned up in view of direct use of libods.a
!  16Jul09   Redder     Increased l2d_max and l3d_max to 20
!
!EOP
!-------------------------------------------------------------------------





