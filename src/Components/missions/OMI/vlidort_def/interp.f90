MODULE interp

IMPLICIT NONE

contains

subroutine interp1( xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: xData(:), yData(:)

  real, intent(in) :: xVal
  real(KIND=8),intent(out) :: yVal
  integer :: inputIndex, dataIndex
  real :: minXdata, minYdata, xRange, weight, maxXData

  ! Possible checks on inputs could go here
  ! Things you may want to check:
  !   monotonically increasing xData
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)

  minXData = xData(1)
  maxXData = xData(size(xData))
  xRange = xData(2)-xData(1)


  ! possible checks for out of range xVal could go here
  
  ! this will work if x is uniformly spaced, otherwise increment
  ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
  dataIndex = floor((xVal-minXData)/xRange) + 1
  weight = (xVal - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
  yVal = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)


end subroutine

END MODULE interp
