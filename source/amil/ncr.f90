SUBROUTINE ncr(n, r, ncomb, ier)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-20  Time: 18:08:52

!     Calculate the number of different combinations of r objects out of n.

!     ier = 0 if no error is detected
!         = 1 if n < 1
!         = 2 if r < 0
!         = 3 if r > n
!         = 4 if nCr > 1.e+308, i.e. if it overflows.  In this case, the
!                natural log of nCr is returned.

!     Programmer: Alan.Miller @ cmis.csiro.au
!     Latest revision - 28 July 1988 (Fortran 77 version)

IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: r
REAL (dp), INTENT(OUT)  :: ncomb
INTEGER, INTENT(OUT)    :: ier

INTEGER :: rr, i, nn

INTERFACE
  FUNCTION lngamma(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION lngamma
END INTERFACE

IF (n < 1) THEN
  ier = 1
ELSE IF (r < 0) THEN
  ier = 2
ELSE IF (r > n) THEN
  ier = 3
ELSE
  ier = 0
END IF
IF (ier /= 0) RETURN

IF (r <= n-r) THEN
  rr = r
ELSE
  rr = n - r
END IF

IF (rr == 0) THEN
  ncomb = 1.0_dp
  RETURN
END IF

IF (rr > 25) THEN
  ncomb = lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1))
  IF (ncomb > 709._dp) THEN
    ier = 4
  ELSE
    ncomb = EXP(ncomb)
  END IF
  RETURN
END IF

ncomb = n
i = 1
nn = n
DO
  IF (i == rr) RETURN
  nn = nn - 1
  i = i + 1
  ncomb = (ncomb * nn) / REAL(i)
END DO

RETURN
END SUBROUTINE nCr



PROGRAM test_nCr
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTERFACE
  SUBROUTINE ncr(n, r, ncomb, ier)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: r
    REAL (dp), INTENT(OUT)  :: ncomb
    INTEGER, INTENT(OUT)    :: ier
  END SUBROUTINE ncr
END INTERFACE

INTEGER    :: n, r, ier
REAL (dp)  :: result

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter n, r : '
  READ(*, *) n, r
  CALL nCr(n, r, result, ier)
  IF (ier /= 0) THEN
    WRITE(*, *) ' Error, IER = ', ier
    IF (ier == 4) WRITE(*, '(a, f12.5)') ' Ln(nCr) = ', result
  ELSE
    WRITE(*, '(a, g16.8)') ' nCr = ', result
  END IF
END DO

STOP
END PROGRAM test_nCr



FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL(dp), INTENT(IN)  :: z
REAL(dp)              :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                              -1259.139216722289D0, 771.3234287757674D0, &
                              -176.6150291498386D0, 12.50734324009056D0, &
                              -0.1385710331296526D0, 0.9934937113930748D-05, &
                               0.1659470187408462D-06 /), zero = 0.D0,   &
                               one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
                               half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER          :: j

IF (z <= zero) THEN
  WRITE(*, *) 'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma


