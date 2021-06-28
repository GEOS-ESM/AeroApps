      subroutine binObs2D(lon,lat,obs,nobs,gObs,im,jm,missing)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm)
!
! Bins obss of no-fire pixels. It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,10)
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, n
      real :: dLon, dLat, xLon
      real :: nsample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0
      nSample = 0.0
      do n = 1, nobs

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = 1 + nint((xlon + 180.  ) / dLon)
         j = 1 + nint((lat(n) + 90. ) / dLat)

         if ( i>im ) i = i - im
         if ( i<1  ) i = i + im
         
         if (abs(obs(n)-missing) > 0.1) then
         gObs(i,j) = gObs(i,j) + obs(n)
         nSample(i,j) = nSample(i,j) + 1.0
         end if
      end do

      where ( nsample > 0 ) 
              gObs = gObs / nSample
      elsewhere
              gObs = missing
      end where

      end

!-----------

      subroutine binObs3D(lon,lat,obs,nobs,gObs,im,jm,km,missing)
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
! Bins obss of no-fire pixels. It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, k, n
      real :: dLon, dLat, xLon
      real :: nSample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0

      do k = 1, km
         nSample = 0.0
         
         do n = 1,nobs
                 
            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
       
            if (abs(obs(n,k)-missing) > 0.1) then
               gObs(i,j,k) = gObs(i,j,k) + obs(n,k)
               nSample(i,j) = nSample(i,j) + 1.0
            end if

         end do

         where ( nsample > 0 ) 
            gObs(:,:,k) = gObs(:,:,k) / nSample
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end
