
module Linespec_aux_m

!  Contains the following auxiliary routines:

!      FUNCTION ibin
!      FUNCTION StrUpCase
!      FUNCTION StrLowCase

!      SUBROUTINE reverse
!      SUBROUTINE reverse_idxs
!      SUBROUTINE gauss
!      SUBROUTINE gaussio
!      SUBROUTINE gauss_f2c
!      SUBROUTINE gauss_f2cio


!   - All Subroutines and functions are stand-alone, no dependencies
!   - All Subroutines and functions are Implicit none, and public.

INTEGER, PARAMETER :: dp = KIND(1.0D0)

public
private :: dp

contains

FUNCTION ibin (vtarget, array, nentries) RESULT(idx)

! binary search in an array of real numbers in increasing order.
! returned is the number of the last entry which is less than target, or
! 0 if not within array. (this was written to find values enclosing
! target for a linear interpolation scheme.) 4/9/84 john lavagnino;
! adapted from jon bentley, cacm february 1984, vol. 27, no. 2, p. 94.

IMPLICIT NONE
INTEGER, INTENT(IN)      :: nentries
REAL(KIND=dp), INTENT(IN) :: vtarget
REAL(KIND=dp), DIMENSION(nentries), INTENT(IN) :: array

INTEGER :: upper, lower, middle
INTEGER :: idx

lower = 0
upper = nentries + 1

DO WHILE (lower + 1 /= upper)
  middle = (lower + upper) / 2
  IF (array(middle) < vtarget) THEN
     lower = middle
  ELSE
     upper = middle
  ENDIF
ENDDO

! at this point, either array (lower) <= target <= array (upper), or
! lower = 0, or upper = nentries + 1 (initial values).
IF (lower > 0 .and. upper /= nentries + 1) THEN
   idx = lower
ELSE
   idx = 0
ENDIF

END FUNCTION ibin

!!

!-----------------------------------------------------------------------------
!S+
! NAME:
!       StrUpCase
!
! PURPOSE:
!       Function to convert an input string to upper case.
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       Result = StrUpCase( String )
!
! INPUT ARGUMENTS:
!       String:  Character string to be converted to upper case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( * )
!                DIMENSION:  Scalar
!                ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result:  The input character string converted to upper case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( LEN(String) )
!                DIMENSION:  Scalar
!
! CALLS:
!       None.
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! EXAMPLE:
!       string = 'this is a string'
!       WRITE( *, '( a )' ) StrUpCase( string )
!   THIS IS A STRING
!
! PROCEDURE:
!       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!       1995 Springer-Verlag, New York.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )
  
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  
  ! -- Local variables
  INTEGER :: i, n
  
  ! -- Copy input string
  Output_String = Input_String
  
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
     
     ! -- Find location of letter in lower case constant string
     n = INDEX( LOWER_CASE, Output_String( i:i ) )
     
     ! -- If current substring is a lower case letter, make it upper case
     IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
     
  END DO
  
END FUNCTION StrUpCase


!------------------------------------------------------------------------------
!S+
! NAME:
!       StrLowCase
!
! PURPOSE:
!       Function to convert an input string to lower case.
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       Result = StrLowCase( String )
!
! INPUT ARGUMENTS:
!       String: Character string to be converted to lower case.
!               UNITS:      N/A
!               TYPE:       CHARACTER( * )
!               DIMENSION:  Scalar
!               ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result:  The input character string converted to lower case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( LEN(String) )
!                DIMENSION:  Scalar
!
! CALLS:
!       None.
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! EXAMPLE:
!       string = 'THIS IS A STRING'
!       WRITE( *, '( a )' ) StrLowCase( string )
!   this is a string
!
! PROCEDURE:
!       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!       1995 Springer-Verlag, New York.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String )
  
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  
  ! -- Local variables
  INTEGER :: i, n
  
  
  ! -- Copy input string
  Output_String = Input_String
  
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
     
     ! -- Find location of letter in upper case constant string
     n = INDEX( UPPER_CASE, Output_String( i:i ) )
     
     ! -- If current substring is an upper case letter, make it lower case
     IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
     
  END DO
  
END FUNCTION StrLowCase

SUBROUTINE reverse ( inarr, num )
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  INTEGER, INTENT(IN) :: num
  INTEGER             :: i
  REAL (KIND=dp), DIMENSION(1: num), INTENT(INOUT) :: inarr
  REAL (KIND=dp), DIMENSION(1: num)                :: temp

  DO i = 1, num
     temp(i) = inarr(num - i + 1)
  ENDDO
  inarr = temp

  RETURN
END SUBROUTINE reverse

SUBROUTINE reverse_idxs (num, idxs )
  IMPLICIT NONE

  INTEGER,                  INTENT(IN) :: num
  INTEGER, DIMENSION(num), INTENT(OUT) :: idxs
  INTEGER                              :: i

  DO i = 1, num
     idxs(i) = num - i + 1
  ENDDO

  RETURN
END SUBROUTINE reverse_idxs

! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity)
!
! The :symetric Gaussian g(x) is defined as
!                   _                _
!                  |         x^2      |
!      g(x) =  EXP | - -------------- |
!                  |_     hw1e ^2    _|
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================

! Assume wavelength grid is evenly spaced and output/input has the same spectral grid

SUBROUTINE gauss (wvlarr, specarr, specmod, npoints, fwhm)
  
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                            INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: fwhm
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                            :: nhi, nlo, i, j, num_slit
  REAL (KIND=dp)                      :: emult, delwvl, slitsum, slit0, hw1e
  REAL (KIND=dp), DIMENSION (npoints) :: slit

  hw1e = fwhm / 1.66551

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0 ) RETURN

  emult  = -1.0 / ( hw1e * hw1e )
  delwvl = wvlarr(2) - wvlarr(1)
 
  !  Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.

  slitsum = 1.0 ;  slit0 = 1.0
  i = 1  ;  num_slit = 0
  DO WHILE ( num_slit <= npoints )
     slit (i) = EXP (emult * (delwvl * i)**2)
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit (i) <= 0.001 ) EXIT 
     i = i + 1
  ENDDO
  num_slit = i
  
  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! Convolve spectrum.  reflect at endpoints.
  ! Doesn't look right
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = i - j 
        IF (nlo < 1) nlo = -nlo + 2 
        nhi = i + j
        IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE gauss

SUBROUTINE gaussio (wvlarr, specarr, i0, specmod, npoints, fwhm, scalex)

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                            INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: fwhm, scalex
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr, i0
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=dp), DIMENSION (npoints) :: i0mod, abspec, abspecmod

  specmod = specarr
  if (fwhm == 0.0) RETURN

  CALL gauss(wvlarr, i0, i0mod, npoints, fwhm)
  abspec = i0 * exp(-specarr * scalex)
  CALL gauss(wvlarr, abspec, abspecmod, npoints, fwhm)

  specmod = -LOG(abspecmod / i0mod) / scalex

  RETURN

END SUBROUTINE gaussio

! convolve high-resolution spectra to low resolution spectra
! fwave: high-resolution wavelength grid
! fspec: high-resolution reference spectra
! nf:    number of wavelengths at high-resolution
! nspec: number of spectra
! fwhm:  slit function in terms of FWHM (nm)
! cwave: low-resolution wavelength grid (could be the same as high-resolution grid)
! cspec: spectra at fwhm
! nc:    number of wavelengths for the low-resolution grid
SUBROUTINE gauss_f2c (fwave, fspec, nf, nspec, fwhm, cwave, cspec, nc)
  
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                               INTENT (IN) :: nc, nf, nspec
  REAL (KIND=dp),                         INTENT (IN) :: fwhm  
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                       :: i, j, midx, sidx, eidx, nhalf
  REAL (KIND=dp)                 :: hw1esq, dfw, ssum, hw1e
  REAL (KIND=dp), DIMENSION (nf) :: slit

  if (fwhm == 0.0) then
     return
  endif

  dfw  = fwave(2) - fwave(1) 
  hw1e = fwhm / 1.66551; hw1esq = hw1e ** 2
  nhalf  = hw1e / ABS(dfw) * 2.65
  
  DO i = 1, nc
     ! Find the closest pixel
     if (dfw > 0) then
        midx = MINVAL(MAXLOC(fwave, MASK=(fwave <= cwave(i)))) + 1
     else
        midx = MINVAL(MINLOC(fwave, MASK=(fwave <= cwave(i)))) + 1
     endif
     
     sidx = MAX(midx - nhalf, 1)
     eidx = MIN(nf, midx + nhalf)
     slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / hw1esq )
     
     ssum = SUM(slit(sidx:eidx))
     DO j = 1, nspec
        cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
     ENDDO
  ENDDO
  
  RETURN
  
END SUBROUTINE gauss_f2c

! Same as above, except for correcting solar I0 effect using high resolution
! solar reference spectrum following the Beer's Law
! i0: high resolution solar reference
! scalex: scaling, normally number of molecules
SUBROUTINE gauss_f2ci0(fwave, fspec, i0, nf, nspec, scalex, fwhm, cwave, cspec, nc)

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                               INTENT (IN) :: nc, nf, nspec
  REAL (KIND=dp),                         INTENT (IN) :: fwhm, scalex  
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave, i0
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, ntemp
  REAL (KIND=dp), DIMENSION(nc)        :: ci0
  REAL (KIND=dp), DIMENSION(nc, nspec) :: cabspec
  REAL (KIND=dp), DIMENSION(nf, nspec) :: abspec

  if (fwhm == 0.0) then
     return
  endif

  ntemp = 1
  CALL gauss_f2c (fwave, i0, nf, ntemp, fwhm, cwave, ci0, nc)
  
  ! Follow Beer's Law
  DO i = 1, nspec
     abspec(:, i) = i0 * EXP(-fspec(:, i) * scalex)
  ENDDO
  
  CALL gauss_f2c (fwave, abspec, nf, nspec, fwhm, cwave, cabspec, nc)

  DO i = 1, nspec
     cspec(:, i) = - LOG(cabspec(:, i)/ ci0) / scalex
  ENDDO

  RETURN

END SUBROUTINE gauss_f2ci0

!  End module

end module Linespec_aux_m
