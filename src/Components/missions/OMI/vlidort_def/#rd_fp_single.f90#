MODULE rd_fp_single

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
  integer(HSIZE_T), dimension(2) :: dims

  integer(HSIZE_T), dimension(4) :: dim3

  type aer_data ! for one pixel
    integer ::       nlayer = 72, nwv
    integer ::       nPol = 6     !numbers of moments
    integer ::       nMom = 300   !numbers of legendre coef
    real    ::       aerflv = 9.9999999E14, p0, sza, raa, vza, chl, wsp, lat, lon
    integer    ::       year, month, day, orbit, line, row 

    real (kind = 4), allocatable, dimension (:)   :: wvs
    real (kind = 4), dimension (72)   :: h_mid,delp, T, oz_layer
    real (kind = 4), allocatable, dimension (:,:)   :: dtau, ssa
    real (kind = 4), allocatable, dimension (:,:)   :: o4_xsec, oz_xsec
    real (kind = 4), dimension (6,300,72)   :: pmom
  endtype aer_data

contains

subroutine rd_aerpfl (aerpfname,flval,aer_data_rec)
! read in MERRA-2 aerosol profiles

  implicit none

  type (aer_data)        :: aer_data_rec
  real                   :: flval, wave
  character(len=*)       :: aerpfname

! local
  integer ::       nlay, i,j,ilay,iwv,  iw, m2_nwv, ozxsec_nwv, o4xsec_nwv, vl_nwv 

  real (kind = 4), dimension (:,:), allocatable   :: dtau, ssa
  integer(HSIZE_T), dimension(1) :: dim_px, dim_wv
  integer(HSIZE_T), dimension(1) :: dim2
  integer(HSIZE_T), dimension(2) :: dim_xsec
  real (kind = 4), allocatable, dimension (:)   ::  m2_wvs

  dim_px(1) = 1

  CALL h5open_f(hdferr)
  call h5fopen_f(aerpfname, H5F_ACC_RDONLY_F, fileId, hdferr)
  if(hdferr /= 0) stop "h5fopen aerpfname failed"


  call h5dopen_f(fileId,'nwv_merra2',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen Num Waves failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,m2_nwv,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread Num Waves failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'nwv_vl',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen Num Waves failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,vl_nwv,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread Num Waves failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  !call h5dopen_f(fileId,'nwv_ozxsec',dataId,hdferr) ! for p0
  !if(hdferr /= 0) stop "h5dopen Num Waves failed"
  !call h5dread_f(dataId,H5T_NATIVE_INTEGER,ozxsec_nwv,dim_px,hdferr)
  !if(hdferr /= 0) stop "h5dread Num Waves failed"
  !call h5dclose_f(dataId, hdferr) ! close the dataset

  !call h5dopen_f(fileId,'nwv_o4xsec',dataId,hdferr) ! for p0
  !if(hdferr /= 0) stop "h5dopen Num Waves failed"
  !call h5dread_f(dataId,H5T_NATIVE_INTEGER,o4xsec_nwv,dim_px,hdferr)
  !if(hdferr /= 0) stop "h5dread Num Waves failed"
  !call h5dclose_f(dataId, hdferr) ! close the dataset


  
  nlay=aer_data_rec%nlayer

  dims(2)=nlay; dims(1)=m2_nwv

  dim_xsec(2)=nlay; dim_xsec(1)=31


  if (allocated(aer_data_rec%ssa)) deallocate(aer_data_rec%ssa)
  allocate(aer_data_rec%ssa(vl_nwv,nlay))


  if (allocated(aer_data_rec%dtau)) deallocate(aer_data_rec%dtau)
  allocate(aer_data_rec%dtau(vl_nwv,nlay))

  !if (allocated(aer_data_rec%o4_xsec)) deallocate(aer_data_rec%o4_xsec)
  !allocate(aer_data_rec%o4_xsec(o4xsec_nwv,nlay))

  !if (allocated(aer_data_rec%oz_xsec)) deallocate(aer_data_rec%oz_xsec)
  !allocate(aer_data_rec%oz_xsec(ozxsec_nwv,nlay))


  if (allocated(aer_data_rec%wvs)) deallocate(aer_data_rec%wvs)
  allocate(aer_data_rec%wvs(vl_nwv))


  if (allocated(m2_wvs)) deallocate(m2_wvs)
  allocate(m2_wvs(m2_nwv))


  if (allocated(dtau)) deallocate(dtau)
  if (allocated(ssa)) deallocate(ssa)
  allocate(dtau(vl_nwv,nlay),ssa(vl_nwv,nlay))


! initializing the output parameters
  aer_data_rec%p0=flval
  aer_data_rec%h_mid=flval

  aer_data_rec%delp=flval
  aer_data_rec%T=flval
  aer_data_rec%dtau=flval
  aer_data_rec%ssa=flval
  aer_data_rec%wvs=flval


  dim2(1)=nlay
  dim_wv(1) = vl_nwv



  call h5dopen_f(fileId,'PS',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen PS failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%p0,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread PS failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Lat',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen lat failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%lat,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread lat failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Lon',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen lon failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%lon,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread lon failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Year',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen year failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%year,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread year failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Month',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen month failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%month,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread month failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Day',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen day failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%day,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread day failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Orbit',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen orb failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%orbit,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread orb failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Line',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen line failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%line,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread line failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'Row',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen row failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,aer_data_rec%row,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread row failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'VZA',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen VZA failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%vza,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread VZA failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'RAA',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen RAA failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%raa,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread RAA failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'SZA',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen SZA failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%sza,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread SZA failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'WSP',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen WSP failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%wsp,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread WSP failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'CHL',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen CHL failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%chl,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread CHL failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  call h5dopen_f(fileId,'H',dataId,hdferr) ! for h_mid
  if(hdferr /= 0) stop "h5dopen H failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%h_mid,dim2,hdferr)
  if(hdferr /= 0) stop "h5dread H failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'LayerOzone',dataId,hdferr) ! for h_mid
  if(hdferr /= 0) stop "h5dopen H failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%oz_layer,dim2,hdferr)
  if(hdferr /= 0) stop "h5dread H failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

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

  call h5dopen_f(fileId,'m2_waves',dataId,hdferr) ! for wvs
  if(hdferr /= 0) stop "h5dopen Wavelength failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,m2_wvs,dim_wv,hdferr)
  if(hdferr /= 0) stop "h5dread Wavelength failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset

  call h5dopen_f(fileId,'vl_waves',dataId,hdferr) ! for wvs
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
  call h5dread_f(dataId,H5T_NATIVE_REAL,dtau,dims,hdferr)
  call h5dclose_f(dataId, hdferr) ! close the dataset

  !call h5dopen_f(fileId,'o4_xsec',dataId,hdferr) ! for dtau
  !if(hdferr /= 0) stop "h5dopen tau failed"
  !call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%o4_xsec,dim_xsec,hdferr)
  !call h5dclose_f(dataId, hdferr) ! close the dataset

  !call h5dopen_f(fileId,'ozone_xsec',dataId,hdferr) ! for dtau
  !if(hdferr /= 0) stop "h5dopen tau failed"
  !call h5dread_f(dataId,H5T_NATIVE_REAL,aer_data_rec%oz_xsec,dim_xsec,hdferr)
  !call h5dclose_f(dataId, hdferr) ! close the dataset

  aer_data_rec%ssa = ssa
  aer_data_rec%dtau = dtau

  aer_data_rec%nwv = vl_nwv

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


subroutine rd_phsfunc (pcoefname,flval,aer_data_rec,wave)
! read in aerosol phase function (legendre coefs)
  implicit none

  type (aer_data)        :: aer_data_rec
  real                   :: flval , wave
  character(len=*)       :: pcoefname
  integer ::       nPol,nMom, nlay, nwav, i,j, k
  real (kind = 4), allocatable, dimension (:)   :: wavelength
  real (kind = 4), allocatable, dimension (:,:,:,:)   :: pmom
  integer(HSIZE_T), dimension(1) :: dim_px 
  integer :: wave_ind

  dim_px(1) = 1

! read in legendre coefficients


  CALL h5open_f(hdferr)
  call h5fopen_f(pcoefname, H5F_ACC_RDONLY_F, fileId, hdferr)
  if(hdferr /= 0) stop "h5fopen pcoefname failed"

  call h5dopen_f(fileId,'nwv_vl',dataId,hdferr) ! for p0
  if(hdferr /= 0) stop "h5dopen Num Waves failed"
  call h5dread_f(dataId,H5T_NATIVE_INTEGER,nwav,dim_px,hdferr)
  if(hdferr /= 0) stop "h5dread Num Waves failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  dim3(1)=npol; dim3(2)=nmom;  dim3(3) = nwav; dim3(4) = nlay
  dims_wave(1) = nwav

 if (allocated(wavelength)) deallocate(wavelength)
  allocate(wavelength(nwav))

  call h5dopen_f(fileId,'vl_waves',dataId,hdferr) ! for ssa
  if(hdferr /= 0) stop "h5dopen waves failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,wavelength,dims_wave,hdferr)
  if(hdferr /= 0) stop "h5dread waves failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset


  do i=1,nwav
     PRINT *, wavelength(i),wave
     if (wavelength(i) .EQ. wave ) THEN
        wave_ind = i
     endif
  enddo


  nlay=aer_data_rec%nlayer
  npol=aer_data_rec%nPol
  nmom=aer_data_rec%nMom

  if (allocated(pmom)) deallocate(pmom)
  allocate(pmom(npol, nmom, nwav,nlay))

  !if (allocated(aer_data_rec%pmom)) deallocate(aer_data_rec%pmom)
  !allocate(aer_data_rec%pmom(npol, nmom,nlay))

  aer_data_rec%pmom=flval

  call h5dopen_f(fileId,'pmom',dataId,hdferr) ! for ssa
  if(hdferr /= 0) stop "h5dopen pmom failed"
  call h5dread_f(dataId,H5T_NATIVE_REAL,pmom,dim3,hdferr)
  if(hdferr /= 0) stop "h5dread pmom failed"
  call h5dclose_f(dataId, hdferr) ! close the dataset
  call h5fclose_f(fileId, hdferr) ! close the hdf file
  CALL h5close_f(hdferr)
  aer_data_rec%pmom = pmom(:,:,wave_ind,:)
  




end subroutine rd_phsfunc

!  End modul



end module rd_fp_single
