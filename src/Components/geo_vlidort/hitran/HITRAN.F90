module HITRAN
! Module to calculate spectra convolved with instrument functions from HITRAN line parameters
!
!
! Comments from Kelly's original code - with some clarifications from P. Castellanos
!
! (1) molnum is the hitran molecule number.
!
! (2) start, step, and npoints define the calculation grid for the high resolution spectrum. 
!     in order to avoid undersampling the lines the grid should be at least as fine as 1/3 
!     of the smallest gaussian fwhm of a spectral line, OR 0.5550 times the smallest gaussian hw1e.
!     this insures that the maximum sampling error is <1.24e-4 of the full-scale
!     line-shape. later: add automatic check for undersampling wrt voigt widths
!
! (3) Press is the pressure in atmospheres (millibars / 1013.25). 
!
! (4) Temp is the temperature in degrees kelvin.
!
! (5) nvoigt is the number of grid points to each side of a spectral line for
!     performing the voigt calculation. assuming calculation to <= 1e-4 is desired,
!     nvoigt should be the greater of (A) 100 * hwhm_max / step, where hwhm_max is
!     the largest lorentzian hwhm of a line and step is the grid spacing;
!     (B) 3.035 * hw1e_max / step, where hw1e_max is the largest gaussian hw1e of a
!     line.
!
! (6) hw1e is the gaussian slit with at 1/e intensity. hw1e = the full width at
!     half-maximum, fwhm, * 0.60056 = fwhm / (2 * sqrt (ln 2)).
!


! Notes: 
! 1. Qpower needs to be more accurate
! 2. pressure-induced shift is not considered
! 3. Intensity for different isotopes is weighted by their fraction in the atmosphere
! 4. Is nvoigt enough?
! 5. Need to speed up Voigt calculations?
! 6. Sinc function?

! Updates
! Feb. 2012: add HUMLIK function, which improves the speed of the calculation by a factor of 5 compared to voigt
! Feb. 25, 2012: add crsdt
! Aug. 27, 2012: Update Kelly's new partition function

! RT Solutions 12 March 2013
!   * Make into a module
!   * Introduce HITRAN path as variable

! RT Solutions 07 July 2013
!   * Linearized Version of code
!   * Modularize code more, add exception handling

! P. Castellanos July 2016
! (1) Simplicf code.  Add comments and formatting. Change subroutine and variable naming.
! (2) real*8

  use Linespec_aux_m  , only : gauss_f2c, ibin, StrLowCase, StrUpCase, Reverse
  use Linespec_basic_m, only : q_load, q_lookup, hitran_setup, voigt, Humlik 
  use Linespec_plus_m , only : q_lookup_lin, Lin_Humlik 

  implicit None

  PUBLIC  get_HITRAN_cross_section

  contains

  SUBROUTINE get_HITRAN_cross_section &
         ( Hitran_path, the_molecule, nlambda, lambda, &
           is_wavenum, nz, ps, ts, fwhm, do_dP, do_dT, &
           crs, errstat, message, crsdt, crsdp )



  ! Input arguments
  CHARACTER(LEN=*),                  INTENT(IN)      :: Hitran_path   ! directory for hitran line parameter files
  CHARACTER (LEN=6),                 INTENT(IN)      :: the_molecule  !  
  INTEGER,                           INTENT(IN)      :: nlambda       ! number of instrument wavelengths
  INTEGER,                           INTENT(IN)      :: nz            ! number of vertical layers
  LOGICAL,                           INTENT(IN)      :: is_wavenum    ! for instrument *** T: wavenumber, F: nm ***
  REAL*8,                            INTENT(IN)      :: fwhm          ! instrument resolution *** in cm^-1 or nm ***
  REAL*8, DIMENSION(nlambda),        INTENT(IN)      :: lambda        ! instrument wavelength/wavenumber grid
  REAL*8, DIMENSION(nz),             INTENT(IN)      :: ps, ts        ! Pressure and temperature vertical profiles
  LOGICAL,                           INTENT(IN)      :: do_dP, do_dT  ! Flags for P and T derivatives

  !  Output arguments
  INTEGER,                           INTENT(OUT)     :: errstat       ! Error status
  CHARACTER*(*),                     INTENT(OUT)     :: message       ! output message
  REAL*8, DIMENSION(nlambda, nz),    INTENT(OUT)     :: crs           ! output coarse cross sections
  REAL*8, DIMENSION(nlambda, nz),    INTENT(OUT),optional     :: crsdt         ! output coarse cross sections temperature derivative
  REAL*8, DIMENSION(nlambda, nz),    INTENT(OUT),optional     :: crsdp         ! output corase cross sections pressure derivative

  !  Local dimensions and names
  INTEGER, PARAMETER :: maxlines  = 110000       ! spectral lines
  INTEGER, PARAMETER :: maxmols   = 42           ! number of molecules in hitran
  INTEGER, PARAMETER :: maxiso    = 8            ! maximum # of isotopes
  INTEGER, PARAMETER :: maxpoints = 210001       ! number of points in spectrum
  CHARACTER (LEN=6), DIMENSION(maxmols), PARAMETER :: molnames = (/ &
       'H2O   ', 'CO2   ', 'O3    ', 'N2O   ', 'CO    ', 'CH4   ', 'O2    ', &
       'NO    ', 'SO2   ', 'NO2   ', 'NH3   ', 'HNO3  ', 'OH    ', 'HF    ', &
       'HCL   ', 'HBR   ', 'HI    ', 'CLO   ', 'OCS   ', 'H2CO  ', 'HOCL  ', &
       'N2    ', 'HCN   ', 'CH3CL ', 'H2O2  ', 'C2H2  ', 'C2H6  ', 'PH3   ', &
       'COF2  ', 'SF6   ', 'H2S   ', 'HCOOH ', 'HO2   ', 'O     ', 'CLONO2', &
       'NO+   ', 'HOBR  ', 'C2H4  ', 'CH3OH ', 'CH3BR ', 'CH3CN ', 'CF4   '/)

  ! constants
  REAL*8, PARAMETER :: pi       = 3.14159265358979d0
  REAL*8, PARAMETER :: one_rtpi = 0.56418958354776d0
  REAL*8, PARAMETER :: c        = 2.99792458d10
  REAL*8, PARAMETER :: p0       = 1013.25d0
  REAL*8, PARAMETER :: t0       = 296.d0          ! hitran standard

  ! codata 2002 constants
  REAL*8, PARAMETER :: h  = 6.6260693d-27
  REAL*8, PARAMETER :: an = 6.0221415d23
  REAL*8, PARAMETER :: r  = 82.057463d0      ! derived
  REAL*8, PARAMETER :: rk = 1.3806505d-16
  REAL*8, PARAMETER :: du = 2.6867773d16    ! derived
  REAL*8, PARAMETER :: c2 = 1.4387752d0

  ! Local arrays
  ! ------------
  INTEGER, DIMENSION(maxlines)             :: mol, iso
  REAL*8, DIMENSION(maxlines)              :: sigma0, strnth, einstein, alpha, &
                                              elow, coeff, selfbrdn, pshift
  REAL*8, DIMENSION(maxmols, maxiso)       :: amu, q296, q
  LOGICAL, DIMENSION(maxmols, maxiso)      :: if_q
  REAL*8, DIMENSION(maxpoints)             :: pos, spec, voigtx, v, posnm


  !  --Linearized------->
  REAL*8, DIMENSION(maxmols, maxiso)       :: dqdT
  REAL*8, DIMENSION(maxpoints)             :: dT_voigtx, dP_voigtx, dT_spec, dP_spec
  REAL*8, DIMENSION(maxpoints,2)           :: dX, dv

  !  Local variables
  !  ---------------
  CHARACTER(LEN=256)                       :: hitran_filename, q_file

  INTEGER          :: molnum, npoints, nvoigt, ntemp
  REAL*8           :: wstart, wend, step, press, temp, minline_vg, &
                      minslit_fwhm, min_fwhm, maxline_vg, maxline_hwhm, voigt_extra

  INTEGER          :: i, mol_temp, iso_temp, nlines, nvlo, nvhi, idx, MEND, iz
  REAL*8           :: sigma0_temp, strnth_temp, einstein_temp, alpha_temp, &
                      selfbrdn_temp, elow_temp, coeff_temp, pshift_temp, sigma_temp
  REAL*8           :: vg, voigta, ratio1, ratio2, ratio, vnorm, rt0t, rc2t, rc2t0
  CHARACTER(LEN=2) :: molc
  CHARACTER(LEN=6) :: the_moleculeU, c6

  !  Linearized

  INTEGER          :: m, mt, mp, nv, moli, isoi
  REAL*8           :: vg_00, dP_vg, dT_vg, voigta_00, term1, term2, term3, term4, term5, term6, dY(2)
  REAL*8           :: dP_term3, dT_term1, dT_term2, dT_term3, dP_term5
  REAL*8           :: dP_ratio1, dP_ratio2, dP_ratio, dP_vnorm, dP_voigta, dT_rt0t, dT_rc2t
  REAL*8           :: dT_ratio1, dT_ratio, dT_vnorm, dT_voigta, amwt

  !  Local Control
  !  =============

  LOGICAL          :: write_diagnostics = .true.
  LOGICAL          :: use_humlik        = .true.

  !  Saved variables for Q-input control
  REAL*8, DIMENSION(maxmols, maxiso, 148:342), SAVE :: q_input 
  LOGICAL,                                     SAVE :: first_qload = .TRUE.                                       


  ! Initialize error status
  errstat = 0 ; message = ' '

!   !  Initialize output

  crs   = 0.0d0
  if (do_dP)   crsdp = 0.0d0
  if (do_dT)   crsdt = 0.0d0

  !  Cannot use non-Humlicek routines with linearization

  if ( (do_dP.or.do_dT) .and..not.use_humlik ) then
     message = 'Failure: Cannot do linearizations with NON-Humlicek Functions!!!'
     errstat = 1 ; RETURN
  endif


  m = 0
  if ( do_dT ) m = m + 1 ; mt = m
  if ( do_dP ) m = m + 1 ; mp = m

  ! Determine which molecule and database file
  the_moleculeU = StrUpCase(the_molecule)
  molnum = 0
  DO i = 1, maxmols
     IF (TRIM(the_moleculeU) == TRIM(molnames(i)) ) THEN 
        molnum = i; EXIT
     ENDIF
  ENDDO
  IF (molnum == 0) THEN
     message = 'Failure: This molecule '//the_molecule//' is not found in HITRAN!!!'
     errstat = 1; RETURN
  ENDIF

  ! Determin start and end wavenumber
  IF (is_wavenum) THEN
     wstart = lambda(1); wend = lambda(nlambda)
  ELSE
     wstart = 1.0D7/lambda(nlambda); wend = 1.0D7/lambda(1)
  ENDIF
  IF (write_diagnostics) WRITE(*, *) 'wstart = ', wstart, ' wend = ', wend

  ! Get HITRAN line parameters file name
  WRITE(molc, '(I2.2)') molnum
  hitran_filename = adjustl(trim(HITRAN_path)) // molc // '_hit16.par'

  IF (write_diagnostics) WRITE(*, *) TRIM(ADJUSTL(hitran_filename))
  OPEN(unit = 22, file = TRIM(ADJUSTL(hitran_filename)), err = 455, status = 'old')

  go to 456       ! successfully opened

  !  Open error status --> Return

  455 continue
  message = 'Failure: Hitran File not found - Look at hitran_filename input !!!'
  errstat = 1 ; return

  !  File found successfully

  456 continue
  ! setup hitrans
  CALL hitran_setup (maxmols, maxiso, amu)

  ! read lines (15 cm^-1 extra on both sides)
  i = 1
  DO 
     READ (22, '(i2, i1, f12.6, 2e10.3, f5.4, f5.3, f10.4, f4.2, f8.6)', IOSTAT = MEND) mol_temp, &
          iso_temp, sigma0_temp, strnth_temp, einstein_temp, alpha_temp, selfbrdn_temp, &
          elow_temp, coeff_temp, pshift_temp
     IF (MEND < 0 .OR. sigma0_temp > wend + 15.0) EXIT
     IF ( (mol_temp .EQ. 2 .AND. iso_temp .EQ. 9) .OR. (mol_temp .EQ. 2 .AND. iso_temp .EQ. 10) &
          .OR. (mol_temp .EQ. 6 .AND. iso_temp .EQ. 4) .OR. (mol_temp .EQ. 27 .AND. iso_temp .EQ. 2) &
          .OR. (mol_temp .EQ. 40) ) CYCLE

     ! only count lines for the specified molecule
     write(*,*),'+++',mol_temp,molnum
     write(*,*),'===',sigma0_temp,wstart-15
     IF (mol_temp == molnum .AND. sigma0_temp > wstart - 15.0 ) THEN
        if_q(mol_temp, iso_temp) = .TRUE.
        mol(i) = mol_temp
        iso(i) = iso_temp
        sigma0(i) = sigma0_temp
        strnth(i) = strnth_temp
        einstein(i) = einstein_temp
        alpha(i) = alpha_temp
        selfbrdn(i) = selfbrdn_temp
        elow(i)  = elow_temp
        coeff(i) = coeff_temp
        pshift(i) = pshift_temp
        i = i + 1
     ENDIF
  ENDDO
  CLOSE(unit = 22)
  nlines = i - 1
  IF (write_diagnostics) WRITE (*, *) 'nlines = ', nlines

  IF (nlines > maxlines) THEN
     message = 'Failure: Nlines > maxlines, need to increase maxlines.!!!'
     errstat = 1; RETURN
  ELSE IF (nlines == 0) THEN
     message =  'Warning: No absorption lines are found in this spectral range!!!'
     errstat = 2; RETURN
  ENDIF

  ! Determine step size for calculating spectrum to avoid undersampling
  minline_vg = MINVAL(4.30140d-7 * sigma0(1:nlines) * dsqrt (t0 / amu(molnum, 1)))   ! it is hw1e
  IF (is_wavenum) THEN
     minslit_fwhm = fwhm 
  ELSE
     minslit_fwhm = (1.0D7/lambda(nlambda) - 1.0D7/(lambda(nlambda) + fwhm))
  ENDIF
  !min_fwhm = SQRT(minslit_fwhm ** 2.0 + minline_vg ** 2.0)
  min_fwhm = minline_vg 
  step = min_fwhm * 0.555

  ! Step size should at least be as fine as the input grid
  ! IF fwhm == 0, then set step size to be the same as input spectral grid
  IF (fwhm == 0.0) THEN
      IF ((wend - wstart) / (nlambda - 1) > step) THEN
        WRITE(*, *) 'Input spectral interval might be too large to resolve the spectral lines!!!'
      ENDIF
      step = (wend - wstart) / (nlambda - 1)    
  ELSE
      IF (nlambda > 1) THEN
        step = MIN(step, (wend - wstart) / (nlambda - 1)) 
      ENDIF
  ENDIF

  IF (write_diagnostics) WRITE(*, *) 'minline_vg = ', minline_vg, ' minslit_fwhm = ', minslit_fwhm
  IF (write_diagnostics) WRITE(*, *) 'Wavenumber step = ', step

  ! Determine nvoigt: number of grid points to each side of a spectral line for
  !     performing the voigt calculation
  maxline_vg = MAXVAL(4.30140d-7 * sigma0(1:nlines) * dsqrt (t0 / amu(molnum, 1)))    ! it is hwle
  maxline_hwhm = MAXVAL(ps) * MAXVAL(alpha(1:nlines))
  IF (write_diagnostics) WRITE(*, *) 'maxline_vg = ', maxline_vg, ' maxline_hwhm = ', maxline_hwhm
  voigt_extra = MAX(maxline_hwhm * 100.0, maxline_vg * 3.035)   ! for hw1e
  nvoigt = INT(voigt_extra / step)
  voigt_extra = nvoigt * step
  IF (write_diagnostics) WRITE(*, *) 'voigt_extra = ', voigt_extra, ' nvoigt = ', nvoigt

  ! Get wavenumbers for high resolution spectra
  ! IF fwhm == 0.0, include the exact original input wavenumber grid with extra edges
  IF (fwhm > 0.0d0) THEN
     npoints = (wend - wstart) / step + 2 * nvoigt
     IF (npoints > maxpoints) THEN
        write(c6,'(I6)')npoints
        message = 'Failure: Npoints = '//c6//', > maxpoints, fwhm>0.0,  need to increase maxpoints.!!!'
        errstat = 1; RETURN
     ENDIF
     
     DO i = 1, npoints
        pos(i) = wstart - voigt_extra + (i - 1) * step
     ENDDO
  ELSE
     npoints = nlambda + 2 * nvoigt
     IF (npoints > maxpoints) THEN
        write(c6,'(I7)')npoints
        message = 'Failure: Npoints = '//c6//', > maxpoints, fwhm=0.0, need to increase maxpoints.!!!'
        errstat = 1; RETURN
     ENDIF
    
     IF (is_wavenum) THEN
        pos(nvoigt + 1:nvoigt + nlambda) = lambda(1:nlambda)
     ELSE
        pos(nvoigt + 1:nvoigt + nlambda) = 1.0D7/lambda(1:nlambda)
        CALL REVERSE(pos(nvoigt + 1:nvoigt + nlambda), nlambda)
     ENDIF

     DO i = nvoigt, 1, -1
        pos(i) = pos(i + 1) - step
     ENDDO

     DO i = nvoigt + nlambda + 1, npoints
        pos(i) = pos(i - 1) + step
     ENDDO
  ENDIF

  posnm(1:npoints) = 1.0D7 / pos(1:npoints)
  CALL REVERSE(posnm(1:npoints), npoints)

  IF (write_diagnostics) THEN
     WRITE(*, *) 'npoints = ', npoints, ' nlambda = ', nlambda
     WRITE(*, *) 'pos(1) = ', pos(1), ' pos(npoints) = ', pos(npoints)
     WRITE(*, *) 'posnm(1) = ', posnm(1), ' posnm(npoints) = ', posnm(npoints)
  ENDIF    

  ! Read in Parition Function Database
  IF (first_qload) THEN
     q_file = adjustl(trim(HITRAN_path)) //'hitran08-parsum.resorted'
     first_qload = .FALSE.
     CALL q_load (maxmols, maxiso, q_input, q_file)   ! @@@ RTS
  ENDIF
  temp = 296.0
  CALL q_lookup (maxmols, maxiso, q_input, molnum, temp, if_q, q296)

  ! Loop over altitude
  DO iz = 1, nz
    ! initialize cross sections and calculate the spectrum grid.
    spec(1:npoints) = 0.d0
    dP_spec(1:npoints) = 0.d0
    dT_spec(1:npoints) = 0.d0
    temp  = ts(iz)
    press = ps(iz)

    if ( .not. do_dT ) then
       CALL q_lookup     (maxmols, maxiso, q_input, molnum, temp, if_q, q)
    else
       CALL q_lookup_lin (maxmols, maxiso, q_input, molnum, temp, if_q, q, dqdt)
    endif

    rt0t = t0 / temp; rc2t = c2 / temp; rc2t0 = c2 / t0
    if ( do_dT ) then
       dT_rt0t = - rt0t / temp ; dT_rc2t = -rc2t / temp
    endif

    ! loop over lines to fill out cross section array
    DO i = 1, nlines
      voigtx = 0.0d0 ; dX = 0.0d0 ; v = 0.0d0 ; dv = 0.0d0

      moli = mol(i) ; isoi = iso(i)
      sigma_temp = sigma0(i) + pshift(i) * press  ! Add pressure induced shift
      amwt  = amu (moli,isoi)

      vg_00 = 4.30140d-7 * sqrt (temp/amwt)
      vg    = vg_00 * sigma_temp
      voigta_00 = alpha(i) * (rt0t)**coeff(i) / vg
      voigta    =  voigta_00 * press

      term1     = dexp(-elow(i) * rc2t)
      term2     = dexp(- sigma_temp * rc2t)
      term3     = term1 * term2
      ratio1    = term1 - term3

      term4     = dexp ( - elow(i) * rc2t0 )
      term5     = dexp ( - sigma_temp * rc2t0 )
      term6     = term4 * term5
      ratio2    = term4 - term6

      ratio = ratio1 / ratio2 * q296(moli, isoi) / q(moli,isoi)
      vnorm    = ratio * strnth(i) / vg

      if ( do_dT ) then
        dT_vg     = 0.5d0 * vg / temp
        dT_voigta = - voigta * ( (coeff(i)/temp) + (dT_vg/vg) )
        dY(mt)    = dT_voigta 
        dT_term1  = - term1 * elow(i) * dT_rc2t
        dT_term2  = - term2 * sigma_temp * dT_rc2t
        dT_term3  = term1 * dT_term2 + term2 * dT_term1
        dT_ratio1 = dT_term1 - dT_term3
        dT_ratio  = ratio * ( (dT_ratio1/ratio1) - (dqdT(moli,isoi)/q(moli,isoi)) )
        dT_vnorm  = vnorm * ( (dT_ratio/ratio) - (dT_vg/vg) ) 
      endif

      if ( do_dP ) then
        dP_vg     = vg_00 * pshift(i)
        dP_voigta =  voigta_00 - voigta * (dP_vg/vg)
        dY(mp)    = dP_voigta
        dP_term3  = - term3 * pshift(i) * rc2t
        dP_ratio1 =  - dP_term3
        dP_term5  = - term5 * pshift(i) * rc2t0
        dP_ratio2 = - term4 * dP_term5
        dP_ratio  = ratio * ( (dP_ratio1/ratio1) - (dP_ratio2/ratio2) )
        dP_vnorm  = vnorm * ( (dP_ratio/ratio) - (dP_vg/vg) ) 
      endif


      idx = ibin (sigma_temp, pos, npoints)
      IF (idx == 0) THEN
        IF (sigma_temp < pos(1) - voigt_extra .OR. sigma_temp > pos(npoints) + voigt_extra) THEN
           nvlo = 0; nvhi = 0
        ELSE IF (sigma_temp < pos(1)) THEN
           nvlo = 1; nvhi = nvoigt
        ELSE
           nvlo = npoints - nvoigt + 1
           nvhi = npoints
        ENDIF
      ELSE
        nvlo = MAX(1, idx - nvoigt)
        nvhi = MIN(npoints, idx + nvoigt)
      ENDIF
      ntemp = nvhi - nvlo + 1

      IF (ntemp > 0.and.nvlo.ne.0) THEN  

        term2 = - pshift(i) / vg
        do nv = nvlo, nvhi 
           voigtx(nv)    = ( pos(nv) - sigma_temp ) / vg
        enddo

        if ( do_dP ) then
           term1 = - dP_vg / vg
           do nv = nvlo, nvhi 
              dP_voigtx(nv) = term1 * voigtx(nv) + term2  ; dX(nv,mp) = dP_voigtx(nv)
           enddo
        endif

        if ( do_dT ) then
           term3 = - dT_vg / vg
           do nv = nvlo, nvhi 
              dT_voigtx(nv) = term3 * voigtx(nv)          ; dX(nv,mt) = dT_voigtx(nv)
           enddo
        endif

        IF (.NOT. use_humlik) THEN
           CALL voigt (voigtx(1:npoints), voigta, v(1:npoints), npoints, nvlo, nvhi)   
           spec(nvlo:nvhi)    = spec(nvlo:nvhi) + vnorm * v(nvlo:nvhi)
        ELSE
           vnorm    = vnorm     * one_rtpi
           if ( do_dP .or. do_dT ) then
              CALL Lin_HUMLIK ( ntemp, m, voigtx(nvlo:nvhi), voigta, dX(nvlo:nvhi,1:2), dY, &
                                v(nvlo:nvhi), dv(nvlo:nvhi,1:2) )
              spec(nvlo:nvhi)    = spec(nvlo:nvhi) + vnorm * v(nvlo:nvhi)
              if ( do_dP ) then
                 dP_vnorm = dP_vnorm  * one_rtpi
                 dP_spec(nvlo:nvhi) = dP_spec(nvlo:nvhi) + dP_vnorm * v(nvlo:nvhi)+ vnorm * dv(nvlo:nvhi,mp)
              endif
              if ( do_dT ) then
                 dT_vnorm = dT_vnorm  * one_rtpi
                 dT_spec(nvlo:nvhi) = dT_spec(nvlo:nvhi) + dT_vnorm * v(nvlo:nvhi)+ vnorm * dv(nvlo:nvhi,mt)
              endif
           else
              CALL HUMLIK   ( ntemp, voigtx(nvlo:nvhi), voigta, v(nvlo:nvhi) ) 
              spec(nvlo:nvhi)    = spec(nvlo:nvhi) + vnorm * v(nvlo:nvhi)
           endif

        ENDIF
      ENDIF

    ENDDO

    ! convolve with instrument function (original)
    IF (fwhm > 0.0d0) THEN
       IF (is_wavenum) THEN
          CALL gauss_f2c (pos(1:npoints), spec(1:npoints), npoints, 1,  &
               fwhm, lambda(1:nlambda), crs(1:nlambda, iz), nlambda)   
       ELSE
          CALL REVERSE(spec(1:npoints), npoints)
          CALL gauss_f2c (posnm(1:npoints), spec(1:npoints), npoints, 1, &
               fwhm, lambda(1:nlambda), crs(1:nlambda, iz), nlambda) 
       ENDIF
    ELSE
       crs(1:nlambda, iz) = spec(nvoigt + 1:nvoigt + nlambda)
       IF (.NOT. is_wavenum) CALL REVERSE(crs(1:nlambda, iz), nlambda)
    ENDIF

!  Convolve linearized arrays with instrument function
    if ( do_dT ) then
      IF (fwhm > 0.0d0 ) THEN
        IF (is_wavenum) THEN
          CALL gauss_f2c (pos(1:npoints), dT_spec(1:npoints), npoints, 1,  &
               fwhm, lambda(1:nlambda), crsdT(1:nlambda, iz), nlambda)   
        ELSE
          CALL REVERSE(dT_spec(1:npoints), npoints)
          CALL gauss_f2c (posnm(1:npoints), dT_spec(1:npoints), npoints, 1, &
               fwhm, lambda(1:nlambda), crsdT(1:nlambda, iz), nlambda) 
        ENDIF
      ELSE
        crsdT(1:nlambda, iz) = dT_spec(nvoigt + 1:nvoigt + nlambda)
        IF (.NOT. is_wavenum) CALL REVERSE(crsdT(1:nlambda, iz), nlambda)
      ENDIF
    endif

    if ( do_dP ) then
      IF (fwhm > 0.0d0 ) THEN
        IF (is_wavenum) THEN
          CALL gauss_f2c (pos(1:npoints), dP_spec(1:npoints), npoints, 1,  &
               fwhm, lambda(1:nlambda), crsdP(1:nlambda, iz), nlambda)   
        ELSE
          CALL REVERSE(dP_spec(1:npoints), npoints)
          CALL gauss_f2c (posnm(1:npoints), dP_spec(1:npoints), npoints, 1, &
               fwhm, lambda(1:nlambda), crsdP(1:nlambda, iz), nlambda) 
        ENDIF
      ELSE
        crsdP(1:nlambda, iz) = dP_spec(nvoigt + 1:nvoigt + nlambda)
        IF (.NOT. is_wavenum) CALL REVERSE(crsdP(1:nlambda, iz), nlambda)
      ENDIF
    endif

  ENDDO  !  End level loop

  RETURN

END SUBROUTINE  get_HITRAN_cross_section

!  End module

END MODULE  HITRAN
