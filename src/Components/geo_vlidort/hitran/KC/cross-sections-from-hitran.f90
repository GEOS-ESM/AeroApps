program cross_sections_from_hitran

implicit real*8 (a - h, o - z)
character*132 input, output
! spectral lines
parameter (maxlines = 60000)
! number of molecules in hitran
parameter (maxmols = 39)
! maximum number of isotopologues per molecule
parameter (maxiso = 8)
! number of points in spectrum
parameter (maxpoints = 600000)

dimension mol (maxlines), iso (maxlines), sigma0 (maxlines), &
  strnth (maxlines), einstein (maxlines), alpha (maxlines), &
  elow (maxlines), coeff (maxlines)
! maximum # species and isotopomers
dimension qpower (maxmols), amu (maxmols, maxiso)
dimension pos (maxpoints), spec (maxpoints), specmod (maxpoints)
dimension voigtx (maxpoints), v (maxpoints)
! constants
pi = 3.14159265358979d0
c = 2.99792458d10
p0 = 1013.25d0
t0 = 296.d0 ! hitran standard
! codata 2002 constants
h = 6.6260693d-27
an = 6.0221415d23
r = 82.057463d0 ! derived
rk = 1.3806505d-16
du = 2.6867773d16 ! derived
c2 = 1.4387752d0

write (*, '(5x, a)') 'enter hitran input file'
read (*, '(a)') input
write (*, '(5x, a)') 'enter cross section output file'
read (*, '(a)') output
open (unit = 1, file = input, status = 'old')
open (unit = 2, file = output, status = 'unknown')

! read logicals and parameters controlling calculations
!
! molnum is the hitran molecule number.
!
! start, step, and npoints define the calculation grid. in order to avoid
! undersampling the grid should be at least as fine as 1/3 of the smallest
! gaussian fwhm of a spectral line, or 0.5550 times the smallest gaussian hw1e.
! this insures that the maximum sampling error is <1.24e-4 of the full-scale
! line-shape. later: add automatic check for undersampling wrt voigt widths
!
! press is the pressure in atmospheres (millibars / 1013.25). temp is the
! temperature in degrees kelvin.
!
! nvoigt is the number of grid points to each side of a spectral line for
! performing the voigt calculation. assuming calculation to <= 1e-4 is desired,
! nvoigt should be the greater of (1) 100 * hwhm_max / step, where hwhm_max is
! the largest lorentzian hwhm of a line and step is the grid spacing;
! (2) 3.035 * hw1e_max / step, where hw1e_max is the largest gaussian hw1e of a
! line.
!
! hw1e is the gaussian slit with at 1/e intensity. hw1e = the full width at
! half-maximum, fwhm, * 0.60056 = fwhm / (2 * sqrt (ln 2)).
!
! nmod gives the option not to write out humongous spectral files by printing
! every nmod^th spectral value
write (*, '(a)') &
  'read molnum, start, step, npoints, press, temp, nvoigt, hw1e, nmod'
read (*, *) molnum, start, step, npoints, press, temp, nvoigt, hw1e, nmod
if (npoints .gt. maxpoints) pause 'npoints .gt. maxpoints'

! initialize cross sections and calculate the spectrum grid.
spec = 0.d0
specmod = 0.d0
do i = 1, npoints
  pos (i) = start + (i - 1) * step
end do

! setup hitran
call hitran_setup (qpower, amu)

! read lines
i = 1
40 read (1, '(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2)', end = 50) mol_temp, &
iso_temp, sigma0_temp, strnth_temp, einstein_temp, alpha_temp, dum, elow_temp, &
coeff_temp
write (*, *) mol_temp, molnum
! only count lines for the specified molecule
  if (mol_temp .eq. molnum) then
    mol (i) = mol_temp
    iso (i) = iso_temp
    sigma0 (i) = sigma0_temp
    strnth (i) = strnth_temp
    einstein (i) = einstein_temp
    alpha (i) = alpha_temp
    elow (i) = elow_temp
    coeff (i) = coeff_temp
    i = i + 1
  end if
  go to 40
50 nlines = i - 1
write (*, *) 'nlines = ', nlines
if (nlines .gt. maxlines) pause 'nlines .gt. maxlines'

! loop over lines to fill out cross section array
do i = 1, nlines
  vg = 4.30140d-7 * sigma0 (i) * dsqrt (temp / amu (mol (i), iso (i)))
  voigta = press * alpha (i) * (t0 / temp)**coeff (i) / vg
  ratio1 = dexp (-elow (i) * c2 / temp) - dexp (-(sigma0 (i) + &
  elow (i)) * c2 / temp)
  ratio2 = dexp (-elow (i) * c2 / t0) - dexp (-(sigma0 (i) + &
  elow (i)) * c2 / t0)
  ratio = ratio1 / ratio2 * ((t0 / temp)**qpower (mol (i)))
  vnorm = ratio * strnth (i) / vg
  nvlo = max (1, ibin (sigma0 (i), pos, npoints) - nvoigt)
  nvhi = min (npoints, ibin (sigma0 (i), pos, npoints) + nvoigt)
  do j = nvlo, nvhi
    voigtx (j) = (pos (j) - sigma0 (i)) / vg
  end do
  call voigt (voigtx, voigta, v, nvlo, nvhi)
  do j = nvlo, nvhi
    spec (j) = spec (j) + vnorm * v (j)
  end do
end do

! convolve with instrument function
call gauss (pos, spec, specmod, npoints, hw1e)

! output results
do i = 1, npoints
  if (nmod .eq. 1 .or. mod (i, nmod) .eq. 1) write (2, '(f10.3, 1p2e13.5)') &
  pos (i), spec (i), specmod (i)
end do

close (unit = 1)
close (unit = 2)
stop
end
!
subroutine hitran_setup (qpower, amu)

implicit real*8 (a - h, o - z)
parameter (maxmols = 39)
parameter (maxiso = 8)

dimension qpower (maxmols), amu (maxmols, maxiso)
qpower = 1.5d0
amu = 0.d0

! hitran numbers:
!   h2o (1)
!   co2 (2)
!    o3 (3)
!   n2o (4)
!    co (5)
!   ch4 (6)
!    o2 (7)
!    no (8)
!   so2 (9)
!   no2 (10)
!   nh3 (11)
!  hno3 (12)
!    oh (13)
!    hf (14)
!   hcl (15)
!   hbr (16)
!    hi (17)
!   clo (18)
!   ocs (19)
!  h2co (20)
!  hocl (21)
!    n2 (22)
!   hcn (23)
! ch3cl (24)
!  h2o2 (25)
!  c2h2 (26)
!  c2h6 (27)
!   ph3 (28)
!  cof2 (29)
!   sf6 (30)
!   h2s (31)
! hcooh (32)
!   ho2 (33)
!     o (34)
!clono2 (35)
!   no+ (36)
!  hobr (37)
!  c2h4 (38)
! ch3oh (39)
! ch3br (40)
! ch3cn (41)
!   cf4 (42)

do i = 1, maxmols
  if (i .eq. 34) qpower (i) = 0.d0
  if (i .eq. 2 .or. i .eq. 4 .or. i .eq. 5 .or. i .eq. 7 .or. i .eq. 8 &
  .or. i .eq. 13 .or. i .eq. 14 .or. i .eq. 15 .or. i .eq. 16 .or. i &
  .eq. 17 .or. i .eq. 18 .or. i .eq. 19 .or. i .eq. 22 .or. i .eq. 23 &
  .or. i .eq. 26 .or. i .eq. 36) qpower (i) = 1.d0
end do
amu (1, 1) = 18.010565  
amu (1, 2) = 20.014811  
amu (1, 3) = 19.014780
amu (1, 4) = 19.016740  
amu (1, 5) = 21.020985  
amu (1, 6) = 20.020956  
amu (2, 1) = 43.989830  
amu (2, 2) = 44.993185  
amu (2, 3) = 45.994076  
amu (2, 4) = 44.994045  
amu (2, 5) = 46.997431  
amu (2, 6) = 45.997400  
amu (2, 7) = 47.998322  
amu (2, 8) = 46.998291  
amu (3, 1) = 47.984745  
amu (3, 2) = 49.988991  
amu (3, 3) = 49.988991  
amu (3, 4) = 48.988960  
amu (3, 5) = 48.988960  
amu (4, 1) = 44.001062  
amu (4, 2) = 44.998096  
amu (4, 3) = 44.998096  
amu (4, 4) = 46.005308  
amu (4, 5) = 45.005278  
amu (5, 1) = 27.994915  
amu (5, 2) = 28.998270  
amu (5, 3) = 29.999161  
amu (5, 4) = 28.999130  
amu (5, 5) = 31.002516  
amu (5, 6) = 30.002485  
amu (6, 1) = 16.031300  
amu (6, 2) = 17.034655  
amu (6, 3) = 17.037475  
amu (7, 1) = 31.989830  
amu (7, 2) = 33.994076  
amu (7, 3) = 32.994045  
amu (8, 1) = 29.997989  
amu (8, 2) = 30.995023  
amu (8, 3) = 32.002234  
amu (9, 1) = 63.961901  
amu (9, 2) = 65.957695  
amu (10, 1) = 45.992904  
amu (11, 1) = 17.026549  
amu (11, 2) = 18.023583  
amu (12, 1) = 62.995644  
amu (13, 1) = 17.002740  
amu (13, 2) = 19.006986  
amu (13, 3) = 18.008915  
amu (14, 1) = 20.006229  
amu (15, 1) = 35.976678  
amu (15, 2) = 37.973729  
amu (16, 1) = 79.926160  
amu (16, 2) = 81.924115  
amu (17, 1) = 127.912297  
amu (18, 1) = 50.963768  
amu (18, 2) = 52.960819  
amu (19, 1) = 59.966986  
amu (19, 2) = 61.962780  
amu (19, 3) = 60.970341  
amu (19, 4) = 60.966371  
amu (19, 5) = 61.971231  
amu (20, 1) = 30.010565  
amu (20, 2) = 31.013920  
amu (20, 3) = 32.014811  
amu (21, 1) = 51.971593  
amu (21, 2) = 53.968644  
amu (22, 1) = 28.006147  
amu (23, 1) = 27.010899  
amu (23, 2) = 28.014254  
amu (23, 3) = 28.007933  
amu (24, 1) = 49.992328  
amu (24, 2) = 51.989379  
amu (25, 1) = 34.005480  
amu (26, 1) = 26.015650  
amu (26, 2) = 27.019005  
amu (27, 1) = 30.046950  
amu (28, 1) = 33.997238  
amu (29, 1) = 65.991722  
amu (30, 1) = 145.962492
amu (31, 1) = 33.987721
amu (31, 2) = 35.983515
amu (31, 3) = 34.987105
amu (32, 1) = 46.005480
amu (33, 1) = 32.997655
amu (34, 1) = 15.994915
amu (35, 1) = 96.956672
amu (35, 2) = 98.953723
amu (26, 1) = 29.997989
amu (37, 1) = 95.921076
amu (37, 2) = 97.919027
amu (38, 1) = 28.031300
amu (38, 2) = 29.034655
amu (39, 1) = 32.026215

return
end
!
subroutine voigt (x, a, v, nvlo, nvhi)

! the following calculated voigt values at all grid values for each point
! subroutine voigt (x, a, v, nx)
! voigt first and second derivatives commented out
! subroutine voigt (x, a, v, dv, d2v, nx)

implicit real*8 (a - h, o - z)
integer capn
real*8 lamda
logical b
parameter (ndim = 600000)
dimension x (ndim), v (ndim)
! dimension x (ndim), v (ndim), dv (ndim), d2v (ndim)
save
data sqrtpi /1.77245385090551d0/, twooverpi /0.63661977236758d0/, &
  fouroverpi /1.27323954473516d0/

! a = 0.
if (a .lt. 1.0d-8) then
  do i = nvlo, nvhi
! do i = 1, nx
    v (i) = dexp (-x (i)**2) / sqrtpi
!   dv (i) = -2.0d0 * x (i) * v (i)
!   d2v (i) = (4.0d0 * x (i) * x (i) - 2.0d0) * v (i)
  end do
! add lorentzian check here, for speed
else
! coefficient for second derivative
! c2v = 4.0d0 * a * a + 2.d0
  sfac = 1.0d0 - a / 4.29d0
  do i = nvlo, nvhi
! do i = 1, nx
    absx = dabs (x (i))
    if ((a .lt. 4.29d0) .and. (absx .lt. 5.33d0)) then
      s = sfac * dsqrt (1.d0 - (x (i) / 5.33d0)**2)
      h = 1.6d0 * s
      h2 = 2.0d0 * h
      capn = 6.d0 + 23.0d0 * s
      lamda = h2**capn
      nu = 9.0d0 + 21.0d0 * s
          else
      h = 0.0d0
      capn = 0
      nu = 8
    end if
    b = (h .eq. 0.0d0) .or. (lamda .eq. 0.0d0)
    r1 = 0.0d0
    r2 = 0.0d0
    s1 = 0.0d0
    s2 = 0.0d0
    n = nu
    do in = 1, nu + 1
      np1 = n + 1
      t1 = a + h + dfloat (np1) * r1
      t2 = absx - dfloat (np1) * r2
      c = .5d0 / (t1 * t1 + t2 * t2)
      r1 = c * t1
      r2 = c * t2
      if ((h .gt. 0.0d0) .and. (n .le. capn)) then
        t1 = lamda + s1
        s1 = r1 * t1 - r2 * s2
        s2 = r2 * t1 + r1 * s2
        lamda = lamda / h2
      end if
      n = n - 1
    end do
    if (b) then
      v (i) = twooverpi * r1
!     dv (i) = fouroverpi * (a * r2 - absx * r1)
    else
      v (i) = twooverpi * s1
!     dv (i) = fouroverpi * (a * s2 - absx * s1)
    end if
!   dv (i) = -dsign (dv (i), x (i))
!   d2v (i) = fouroverpi * a - (c2v + 4.d0 * x (i) * x (i)) * &
!   v (i) - 4.d0 * x (i) * dv (i)
  end do
end if

return
end
!
integer function ibin (target, array, nentries)

! binary search in an array of real numbers in increasing order.
! returned is the number of the last entry which is less than target, or
! 0 if not within array. (this was written to find values enclosing
! target for a linear interpolation scheme.) 4/9/84 john lavagnino;
! adapted from jon bentley, cacm february 1984, vol. 27, no. 2, p. 94.

implicit real*8 (a - h, o - z)
parameter (maxpoints = 600000)
dimension array (maxpoints)
integer upper

lower = 0
upper = nentries + 1
do while (lower + 1 .ne. upper)
  middle = (lower + upper) / 2
  if (array (middle) .lt. target) then
    lower = middle
  else
    upper = middle
  end if
end do

! at this point, either array (lower) <= target <= array (upper), or
! lower = 0, or upper = nentries + 1 (initial values).
if (lower .gt. 0 .and. upper .ne. nentries + 1) then
  ibin = lower
else
  ibin = 0
end if

end function ibin
!
subroutine gauss (pos, spec, specmod, npoints, hw1e)

! convolves input spectrum with gaussian slit function of specified hw1e
! (half-width at 1/e intensity). Assumes input spectrum has constant
! spacing (will need to modify or replace with nm version when and if
! needed).

implicit real*8 (a - h, o - z)
parameter (maxpts = 600000)
parameter (maxslit = 10000)
dimension pos (maxpts), spec (maxpts), specmod (maxpts), slit (maxpts)

write (*, *) 'hw1e = ', hw1e
if (hw1e .eq. 0.d0) then
  write (*, *) ' no gaussian convolution applied.'
  do i = 1, npoints
    specmod (i) = spec (i)
  end do
  return
end if

emult = - 1.d0 / (hw1e**2)
delpos = pos (2) - pos (1)

! apply slit function convolution
write (*, *) 'applying slit function convolution'
! calculate slit function values out to 0.001 times x0 value, normalize
! so that sum = 1.

slitsum = 1.d0
slit0 = 1.d0
do i = 1, maxslit
  slit (i) = exp (emult * (delpos * i)**2)
  slitsum = slitsum + 2.d0 * slit (i)
  if (slit (i) / slit0 .le. 0.001) then
    nslit = i
    write (*, *) ' nslit = ', nslit
    go to 40
  end if
end do
40 if (i .eq. maxslit) write (*, *) 'maxslit is too small'

slit0 = slit0 / slitsum
do i = 1, nslit
  slit (i) = slit (i) / slitsum
end do

! convolve spectrum. don't reflect at endpoints (for now).
do i = 1, npoints
  specmod (i) = slit0 * spec (i)
  do j = 1, nslit
    nlo = i - j
    nhi = i + j
    if (nlo .ge. 1) specmod (i) = specmod (i) + slit (j) * spec (nlo)
    if (nhi .le. npoints) specmod (i) = specmod (i) + slit (j) * &
    spec (nhi)
  end do
end do

return
end
