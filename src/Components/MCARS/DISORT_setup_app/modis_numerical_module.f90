   module modis_numerical_module

implicit none

private

public :: bisectionsearch, linearinterpolation, linearinterpolation_3d, bilinear_interpolation, BisectionSimple

contains


SUBROUTINE BisectionSimple(theta,thetaArray,ntheta,iAngLow, iAngHi)
	
		INTEGER,INTENT(IN)::ntheta
		REAL,INTENT(IN)::theta,thetaArray(1:ntheta)
		INTEGER,INTENT(OUT)::iAngLow, iAngHi
	
		INTEGER:: iAng
	

		iAngLow = 1
		iAngHi = ntheta

		DO
			IF((iAngHi - iAngLow) == 1)EXIT
			iAng = (iAngHi + iAngLow)/2
			IF(thetaArray(iAng) < theta)THEN
				iAngLow = iAng
			ELSE
				iAngHi = iAng
			ENDIF
		ENDDO

		if (iAngLow < 1) iAngLow = 1
		if (ianghi < 1) ianghi = 1
		if (iAngHi > ntheta) iAngHi = ntheta
		if (ianglow > ntheta) ianglow = ntheta

END SUBROUTINE BisectionSimple



! currently bisectionsearch works for ascending arrays only
subroutine bisectionsearch(xx,x,jlo,jhi)
 
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: jlo, jhi
REAL, INTENT(IN) :: x
REAL, DIMENSION(:), INTENT(IN) :: xx
INTEGER :: n,inc,jm
LOGICAL :: ascnd
n=size(xx)
ascnd = (xx(n) >= xx(1))
if (jlo <= 0 .or. jlo > n) then
   jlo=0
   jhi=n+1
else
   inc=1
   if (x >= xx(jlo) .eqv. ascnd) then
      do
         jhi=jlo+inc
         if (jhi > n) then
            jhi=n+1
            exit
         else
            if (x < xx(jhi) .eqv. ascnd) exit
            jlo=jhi
            inc=inc+inc
         end if
      end do
   else
      jhi=jlo
      do
         jlo=jhi-inc
         if (jlo < 1) then
            jlo=0
            exit
         else
            if (x >= xx(jlo) .eqv. ascnd) exit
            jhi=jlo
            inc=inc+inc
         end if
      end do
   end if
end if
do
   if (jhi-jlo <= 1) then
!     NOTE: the following two lines assume that the xx array MUST have been sorted in ascending order for this to work as designed
      if (x >= xx(n)) jlo=n-1
      if (x <= xx(1)) jlo=1
      exit
   else
      jm=(jhi+jlo)/2
      if (x >= xx(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      end if
   end if
end do

! G.Wind 12.16.05 inserted logic to make sure that the values can not go outside
! the valid array bounds
if (jlo < 1) jlo = 1
if (jlo > n) jlo = n
if (jhi < 1) jhi = 1
if (jhi > n) jhi = n

end subroutine bisectionsearch


logical function real_s_equal(x,y)
   real :: x, y
   real_s_equal = (abs(x-y) <= epsilon(x)) 
end function real_s_equal

logical function realsingle_s_equal(x,y)
   real :: x, y
   realsingle_s_equal = (abs(x-y) <= epsilon(x)) 
end function realsingle_s_equal

subroutine realsingle_s_where_equal(x,y)
   real ,intent(inout) :: x(:)
   real :: y

   where( abs(x - y) <= epsilon(y) )
      x = 1.
   elsewhere
      x = 0.
   endwhere

end subroutine realsingle_s_where_equal


real function linearinterpolation(X,Y,XX)


!       XX    R     Interpolation POINT
!       X   R(NN)   INDEPENDENT VARIABLE
!       Y   R(NN)   DEPENDENT   VARIABLE
!
!       Written by      
!        Mark A Gray
!        SM&A SSG .
!        Code 913, NASA/GSFC
!        Greenbelt, MD 20771

   implicit none

   real, intent(in)  ::   X(2), xx 
   real, intent(in)  ::   y(2)

   if (realsingle_s_equal(x(1),x(2))) then
     LinearInterpolation = y(1)
     return
   elseif(realsingle_s_equal(x(1),xx)) then 
     LinearInterpolation = y(1)
     return
   elseif (realsingle_s_equal(x(2),xx)) then 
     LinearInterpolation = y(2)
     return
   else 
     LinearInterpolation = y(1) + ( ( xx - x(1) ) / ( x(2) - x(1)) ) * ( y(2) -y(1) )
   endif

end function linearinterpolation


real function bilinear_interpolation( X, Y, XX, YY, source, method ) 

	implicit none

	real, dimension(2), intent(in) :: X, Y
	real, intent(in) :: XX, YY
	integer, intent(in) :: method
	real, dimension(:,:), intent(in) :: source
	
	real :: mid1, mid2


	real :: area1, area2, area3, area4
	real :: deltaX1, deltaX2, deltaY1, deltaY2
	
	real :: num1, num2
	
	if (method == 1) then 	
	
		deltaX1 = X(1) - XX
		deltaX2 = X(2) - XX
	
		deltaY1 = Y(1) - YY
		deltaY2 = Y(2) - YY
	
		area1 = abs(deltaX1 * deltaY1)
		area2 = abs(deltaX2 * deltaY1)
		area3 = abs(deltaX1 * deltaY2)
		area4 = abs(deltaX2 * deltaY2)
	
		num2 = area1 + area2 + area3 + area4
	
		num1 = source(2, 2) * area1 + &
			source(1, 1) * area4 + &
			source(1, 2) * area2 + &
			source(2, 1) * area3
		   
		bilinear_interpolation = num1 / num2
	else
	
		mid1 = linearinterpolation(X, (/source(1,1), source(1,2)/) , XX)
		mid2 = linearinterpolation(X, (/source(2,1), source(2,2)/) , XX)
	
		bilinear_interpolation = linearinterpolation(Y, (/mid1, mid2/), YY)
	endif


end function bilinear_interpolation


real function LinearInterpolation_3d(index_theta01,index_theta02,  &
                                       index_theta1, index_theta2,  &
                                       index_phi1,   index_phi2,&
                                       theta,theta0,phi,&
                                       theta_grid,  &
                                       theta0_grid,&
                                       phi_grid,&
                                       debug, &
                                       reflectance)
       
       
!Description: 
!       Performa a linear interpolation within a 3d grid for a defined
!       point bounded by defined points.
! 
!Input parameters:
!
!Output parameters:
!
!Revision history: 
!
!Team-unique header: 
!
!References and credits:
!    written by; Mark Gray
!                  
   
   implicit none

   integer, intent(in)  :: index_theta1, index_theta2, index_phi1,index_phi2, &
                           index_theta01, index_theta02
       
   real, intent(in)     :: reflectance(:,:,:)
   real, intent(in)     :: theta_grid(:), &
                             theta0_grid(:),phi_grid(:),theta,theta0,phi
   
   logical, intent(in)    :: debug
      
   real                 :: corner1, corner2, corner3, corner4,  &
                                   corner5, corner6, corner7, corner8
   real                 :: volume1, volume2, volume3, volume4,  &
                                   volume5, volume6, volume7, volume8
   real                 :: num_1, num_2    
   real                 :: theta0_val_1, theta0_val_2, theta_val_1, theta_val_2,phi_val_1, phi_val_2


   theta0_val_1 = theta0_grid(index_theta01)-theta0
   theta0_val_2 = theta0_grid(index_theta02)-theta0

   theta_val_1 = theta_grid(index_theta1)-theta
   theta_val_2 = theta_grid(index_theta2)-theta
 
   phi_val_1 = phi_grid(index_phi1)-phi
   phi_val_2 = phi_grid(index_phi2)-phi

   volume2= abs(theta0_val_1 * theta_val_1 * phi_val_2)
   volume3= abs(theta0_val_1 * theta_val_2 * phi_val_2)
   volume4= abs(theta0_val_1 * theta_val_2 * phi_val_1)
   volume5= abs(theta0_val_2 * theta_val_1 * phi_val_1)
   volume6= abs(theta0_val_2 * theta_val_1 * phi_val_2)
   volume7= abs(theta0_val_2 * theta_val_2 * phi_val_2)
   volume8= abs(theta0_val_2 * theta_val_2 * phi_val_1)
   volume1= abs(theta0_val_1 * theta_val_1 * phi_val_1)

   num_1 =     reflectance(index_theta01,index_theta1,index_phi1) * volume7 + &
               reflectance(index_theta01,index_theta1,index_phi2) * volume8 +  &
               reflectance(index_theta02,index_theta1,index_phi2) * volume4 +  &
               reflectance(index_theta02,index_theta1,index_phi1) * volume3 +  &
               reflectance(index_theta01,index_theta2,index_phi1) * volume6 +  &
               reflectance(index_theta01,index_theta2,index_phi2) * volume5 +  &
               reflectance(index_theta02,index_theta2,index_phi2) * volume1 +  &
               reflectance(index_theta02,index_theta2,index_phi1) * volume2 

   num_2 =     volume1 + volume2 +  &
               volume3 + volume4 +  &
               volume5 + volume6 +  &
               volume7 + volume8   

   LinearInterpolation_3d = num_1/num_2

end function LinearInterpolation_3d

end module modis_numerical_module
