MODULE rd_fp

! Description: read in omglerpre datasets for ocean

! Revision History:
!     initial version: 2018-Feb-19
!     2nd version: 2018-Feb-23
!     3rd version: 2019-Jun-10
!     complete history: svn log file://omi/cm/svn/app/OMPREler

! Author: Wenhan Qin (SSAI)
!*******************************************************************************

use HDF5


IMPLICIT NONE

! key dimensions of datasets:
!
  integer :: hdferr
  integer(HID_T) :: fileId, dataId
  integer(HSIZE_T), dimension(1) :: dims_wave
  integer(HSIZE_T), dimension(3) :: dims

  integer(HSIZE_T), dimension(5) :: dim3

  type aer_data ! for one pixel
    integer ::       nlayer = 72
    integer ::       nPol = 6     !numbers of moments
    integer ::       nMom = 300   !numbers of legendre coef
    real    ::       aerflv = 9.9999999E14 ! fill value for aerosol data

    real (kind = 4), allocatable, dimension (:)   :: wvs, p0
    real (kind = 4), allocatable, dimension (:,:)   :: h_mid,delp, T
    real (kind = 4), allocatable, dimension (:,:,:)   :: dtau, ssa
    real (kind = 4), allocatable, dimension (:,:,:)   :: pmom
  endtype aer_data

contains

subroutine rd_aerpfl (aerpfname,flval,aer_data_rec,nwv,wvs)
! read in MERRA-2 aerosol profiles

  implicit none

  type (aer_data)        :: aer_data_rec
  real                   :: flval, wave
  character(len=*)       :: aerpfname

! local
  integer ::       nlay, i,j,ilay,iwv, nsample=63, ipx, iw
  integer, intent(in) :: nwv
  
  integer(HSIZE_T), dimension(1) :: dim_px, dim_wv
  integer(HSIZE_T), dimension(2) :: dim2
  real (kind = 4), intent(in) :: wvs(:)
  real (kind = 4), allocatable, dimension (:,:,:)   :: tau, ssa
  
  nlay=aer_data_rec%nlayer

  dims(3)=nlay; dims(2)=15; dims(1) = nsample

  if (allocated(aer_data_rec%p0)) deallocate(aer_data_rec%p0)
  allocate(aer_data_rec%p0(nsample))

  if (allocated(aer_data_rec%h_mid)) deallocate(aer_data_rec%h_mid)
  allocate(aer_data_rec%h_mid(nsample,nlay))


  if (allocated(aer_data_rec%T)) deallocate(aer_data_rec%T)
  allocate(aer_data_rec%T(nsample,nlay))

  if (allocated(aer_data_rec%delp)) deallocate(aer_data_rec%delp)
  allocate(aer_data_rec%delp(nsample,nlay))

  if (allocated(aer_data_rec%wvs)) deallocate(aer_data_rec%wvs)
  allocate(aer_data_rec%wvs(dims(2)))

  if (allocated(tau)) deallocate(tau)
  if (allocated(ssa)) deallocate(ssa)
  allocate(tau(nsample,dims(2),nlay),ssa(nsample,dims(2),nlay))


  if (allocated(aer_data_rec%dtau)) deallocate(aer_data_rec%dtau)
  if (allocated(aer_data_rec%ssa)) deallocate(aer_data_rec%ssa)
  allocate(aer_data_rec%dtau(nsample,nwv,nlay),aer_data_rec%ssa(nsample,nwv,nlay))

! initializing the output parameters
  aer_data_rec%p0=flval
  aer_data_rec%h_mid=flval

  aer_data_rec%delp=flval
  aer_data_rec%T=flval
  aer_data_rec%dtau=flval
  aer_data_rec%ssa=flval
  aer_data_rec%wvs=flval




  CALL h5open_f(hdferr)
  call h5fopen_f(aerpfname, H5F_ACC_RDONLY_F, fileId, hdferr)
  if(hdferr /= 0) stop "h5fopen aerpfname failed"


  dim2(1)=nsample
  dim2(2)=nlay

  dim_wv(1) = nwv
  dim_px(1) = nsample


  call h5dopen_f(fileId,'PS',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen PS failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%p0,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread PS failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'H',dataId,hdferr) ! for h_mid
  if(hdferr /= 0) stop "h5dopen H failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%h_mid,dim2,hdferr)
  if(hdferr /= 0) stop "h5dread H failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

! need to use aerflv to screen the data first before assigning it


  call h5dopen_f(fileId,'T',dataId,hdferr) ! for T
  if(hdferr /= 0) stop "h5dopen T failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%T,dim2,hdferr)
  if(hdferr /= 0) stop "h5dread PS failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  call h5dopen_f(fileId,'DELP',dataId,hdferr) ! for delp
  if(hdferr /= 0) stop "h5dopen delp failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%delp,dim2,hdferr)
  if(hdferr /= 0) stop "h5dread delp failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  call h5dopen_f(fileId,'wavelength',dataId,hdferr) ! for wvs
  if(hdferr /= 0) stop "h5dopen Wavelength failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%wvs,dim_wv,hdferr)
  if(hdferr /= 0) stop "h5dread Wavelength failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset
! write(*,*) aer_data_rec%wvs


  call h5dopen_f(fileId,'ssa',dataId,hdferr) ! for ssa
  if(hdferr /= 0) stop "h5dopen ssa failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,ssa,dims,hdferr)
  if(hdferr /= 0) stop "h5dread ssa failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'tau',dataId,hdferr) ! for dtau
  if(hdferr /= 0) stop "h5dopen tau failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,tau,dims,hdferr)
  call h5dclose_f(dataId, hdferr) ! close the dataset

  do ipx = 1, nsample
     do i=1,nlay
        do iw=1,nwv

           call interp1(aer_data_rec%wvs,ssa(ipx,:,i),wvs(iw), aer_data_rec%ssa(ipx,iw,i))
           call interp1(aer_data_rec%wvs,tau(ipx,:,i),wvs(iw), aer_data_rec%dtau(ipx,iw,i))
        enddo
     enddo
  enddo


  if(hdferr /= 0) stop "h5dread tau failed"

  call h5fclose_f(fileId, hdferr) ! close the hdf file
  CALL h5close_f(hdferr)
! stop 'so far so good rd_aer'

end subroutine rd_aerpfl

subroutine interp1( xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: xData(:), yData(:)

  real, intent(in) :: xVal
  real(KIND=4),intent(out) :: yVal
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


subroutine rd_phsfunc (pcoefname,flval,aer_data_rec,wave,ipx)
! read in aerosol phase function (legendre coefs)
  implicit none

  type (aer_data)        :: aer_data_rec
  real                   :: flval , wave
  character(len=*)       :: pcoefname
  integer ::       nPol,nMom, nlay, nwav=15, i,j, k, nsample=63, ipx
  real (kind = 4), allocatable, dimension (:)   :: wavelength
  real (kind = 4), allocatable, dimension (:,:,:,:,:)   :: pmom


! read in legendre coefficients


  CALL h5open_f(hdferr)
  call h5fopen_f(pcoefname, H5F_ACC_RDONLY_F, fileId, hdferr)
  if(hdferr /= 0) stop "h5fopen pcoefname failed"

  dim3(1)=npol; dim3(2)=nmom; dim3(3)=nsample; dim3(4) = nwav; dim3(5) = nlay
  dims_wave(1) = nwav
 if (allocated(wavelength)) deallocate(wavelength)
  allocate(wavelength(nwav))

  call h5dopen_f(fileId,'wavelength',dataId,hdferr) ! for ssa
  if(hdferr /= 0) stop "h5dopen waves failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,wavelength,dims_wave,hdferr)
  if(hdferr /= 0) stop "h5dread waves failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  !do i=1,nwav
  !   PRINT *, wavelength(i),wave
  !   if (wavelength(i) .EQ. wave ) THEN
  !      wave_ind = i
  !   endif
  !enddo
  !print *, wave, wavelength(wave_ind)

  nlay=aer_data_rec%nlayer
  npol=aer_data_rec%nPol
  nmom=aer_data_rec%nMom

  if (allocated(pmom)) deallocate(pmom)
  allocate(pmom(npol, nmom,nsample, nwav,nlay))

  if (allocated(aer_data_rec%pmom)) deallocate(aer_data_rec%pmom)
  allocate(aer_data_rec%pmom(npol, nmom,nlay))

  aer_data_rec%pmom=flval

  call h5dopen_f(fileId,'pmom',dataId,hdferr) ! for ssa
  if(hdferr /= 0) stop "h5dopen pmom failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,pmom,dim3,hdferr)
  if(hdferr /= 0) stop "h5dread pmom failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset
  call h5fclose_f(fileId, hdferr) ! close the hdf file
  CALL h5close_f(hdferr)
  do i=1,npol
     do j = 1,nmom
        do k = 1,nlay
           call interp1(wavelength,pmom(i,j,ipx,:,k),wave, aer_data_rec%pmom(i,j,k))
        enddo
     enddo
  enddo
  !aer_data_rec%pmom = pmom


end subroutine rd_phsfunc

!  End modul



end module rd_fp
