module SOLAR
!
! Module to read in high resolution solar reference spectrum
! Apply instrument operator and interpolate to desired wavelength
!
! July 2016 P. Castellanos
!.............................................................................
  implicit None

  PUBLIC read_SOLAR_reference
  PUBLIC get_SOLAR_nwavenum

  contains

  function get_SOLAR_nwavenum(solar_path,nheader)

    character(len=*),    intent(in)  :: solar_path
    integer,             intent(in)  :: nheader
    integer                          :: get_SOLAR_nwavenum

    integer                          :: i

    open(unit=500, file=solar_path)
    ! Read in header
    do i=1,nheader
      read(500,*,end=10)
    end do

    ! Get number of wavenumber in file
    get_SOLAR_nwavenum = 0
    do 
      read(500,*,end=10)
      get_SOLAR_nwavenum = get_SOLAR_nwavenum + 1
    end do
10  close(500)

  end function get_SOLAR_nwavenum
!.............................................................................

  subroutine read_SOLAR_reference(solar_path,nheader,nwavenum,wavenum,irradiance)

    character(len=*),    intent(in)  :: solar_path
    integer,             intent(in)  :: nheader
    integer,             intent(in)  :: nwavenum
    real*8,              intent(out) :: wavenum(nwavenum)
    real*8,              intent(out) :: irradiance(nwavenum)

    integer                          :: i

    open(unit=500,file=solar_path)
    ! Read in header
    do i=1,nheader
      read(500,*,end=20)
    end do

    do i=1,nwavenum
      read(500,*,end=20) wavenum(i), irradiance(i)
    end do

20  close(500)

  end subroutine read_Solar_reference
!.............................................................................

  subroutine convolve_SOLAR_instrument(wavenum,irradiance,fwhm,cwavenum,cirradiance)

    real*8,             intent(in) :: wavenum(:)
    real*8,             intent(in) :: irradiance(:)
    real*8,             intent(in) :: fwhm
    real*8,             intent(in) :: cwavenum(:)
    real*8,             intent(out):: cirradiance(:)

    integer                        :: i, midx, sidx, eidx
    real*8                         :: dwave
    real*8                         :: hw1e
    real*8                         :: hw1esq
    integer                        :: nhalf
    integer                        :: ncwavenum, nwavenum  
    real*8, allocatable            :: slit(:)  
    real*8                         :: ssum                 


    if (fwhm == 0.0) then
      return
    endif

    dwave  = wavenum(2) - wavenum(1) 
    hw1e   = fwhm / 1.66551; 
    hw1esq = hw1e ** 2
    nhalf  = hw1e / ABS(dwave) * 2.65    
    ncwavenum = SIZE(cwavenum)
    nwavenum  = SIZE(wavenum)

    allocate( slit(nwavenum) )

    do i = 1,ncwavenum
      ! Find the closest wavenumber
      midx = MINVAL(MAXLOC(wavenum, MASK=(wavenum <= cwavenum(i)))) + 1

      sidx = MAX(midx - nhalf, 1)
      eidx = MIN(nwavenum, midx + nhalf)
      slit(sidx:eidx) = EXP(-(cwavenum(i) - wavenum(sidx:eidx))**2 / hw1esq )

      ssum           = SUM(slit(sidx:eidx))
      cirradiance(i) = SUM(irradiance(sidx:eidx) * slit(sidx:eidx)) / ssum
    end do

  end subroutine convolve_SOLAR_instrument
end module SOLAR


