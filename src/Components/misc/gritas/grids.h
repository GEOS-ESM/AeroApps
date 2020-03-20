!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !INCLUDE: grids.h --- GRITAS dimension and grid parameters
! 
! !DESCRIPTION: 
!
!  Defines dimensions and other relevant parameters.
!
! !DEFINED PARAMETERS: 
!

!     Dimensions of O-F Grids (2x2.5, 19 levels)
!     -----------------------------------------      
      integer     im                         ! no. zonal      gridpoints
      integer     jm                         ! no. meridional gridpoints
      integer     km                         ! no. vertical   levels
!!!      parameter ( im=360, jm=181, km=26 ) 
!!!      parameter ( im=144, jm=91, km=26 ) 
      parameter ( im=72, jm=46, km=26 )


!     Longitudinal grid
!     -----------------
      real        LonMin                     ! = lambda(1)
      real        dLon                       ! mesh size (deg)
      real        glon(im)                   ! zonal gridpoints

      parameter ( LonMin = -180. )
      parameter ( dLon = 360. / im )
   

!     Latitudinal grid
!     ----------------
      real        LatMin                     ! = phi(1)
      real        dLat                       ! mesh size (deg)
      real        glat(jm)                   ! neridional gridpoints

      parameter ( LatMin = -90. )
      parameter ( dLat = 180. / ( jm - 1 ) )



!     Missing value
!     -------------
      real        AMISS
      parameter  (AMISS = 1.0E15 )

      common / grid_parms / glat, glon


! !REVISION HISTORY: 
!
!  05sep97   da Silva   Original code.
!
!EOP
!-------------------------------------------------------------------------

