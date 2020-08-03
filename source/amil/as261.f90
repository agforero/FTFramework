MODULE rsquared
IMPLICIT NONE

INTEGER, SAVE :: ifault

PUBLIC :: ifault, qr2

CONTAINS

!------------------------------------------------------------

! The following code for calculating quantiles of R^2 (qr2.for)
! is by the author of AS 261 and uses an improved algorithm.

FUNCTION qr2(m, size, rho2, p) RESULT(quantile)

!  Computes the quantile of the distribution of the square of the
!  sample multiple correlation coefficient for given number of
!  random variables M, sample size SIZE, square of the population
!  multiple correlation coefficient RHO2, and lower tail area P

!  Reference:
!  Ding, C.G. (1996) `On the computation of the distribution of
!  the square of the sample multiple correlation coefficient',
!  Comput. Statist. & Data Analysis, vol. 22, 345-350.

!  IFAULT is a fault indicator:
!  = 1 if there is no convergence after 10 Newton's iterations
!  = 2 if any of the input values is illegal
!  = 0 otherwise

!  No auxiliary algorithm is required

!  ELF90 compatible version by Alan Miller
!  N.B. The error indicator (ifault) has been put into the module contents
!       as functions can only return one result in ELF90.
!  Latest revision - 20 May 1997

INTEGER, INTENT(IN) :: m, size
REAL, INTENT(IN)    :: rho2, p
REAL                :: quantile

! Local variables
REAL               :: n, a, b, ab, ga, q0, yp, y, coeff, q, s, t, cdf, pdf,  &
                      gb, gab, v, bndcdf, bndpdf, ynew, diff
INTEGER            :: na, i, nab, nb, iter
REAL, PARAMETER    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      eps = 1.E-06, delta = 1.E-04, rp = 1.772453850905516028
INTEGER, PARAMETER :: itrmax = 10

quantile = p

!        Test for admissibility of arguments

ifault = 2
IF (m <= 1 .OR. size <= m .OR. rho2 < zero .OR.  &
    rho2 > one .OR. p < zero .OR. p > one) RETURN
ifault = 0
IF (p == zero .OR. p == one) RETURN

!        Calculate the constants needed for each Newton's iteration

a = (m - 1) / two
b = (size - m) / two
ab = (size - 1) / two
IF (MOD(m + 1, 2) == 0) THEN
  na = a + half
  ga = one
  DO i = 1, na
    ga = ga * i
  END DO
ELSE
  na = a + one
  ga = rp
  DO i = 1, na
    ga = ga * (i - half)
  END DO
END IF
IF (MOD(size - m, 2) == 0) THEN
  nb = b - half
  gb = one
  DO i = 1, nb
    gb = gb * i
  END DO
ELSE
  nb = b
  gb = rp
  DO i = 1, nb
    gb = gb * (i - half)
  END DO
END IF
IF (MOD(size - 1, 2) == 0) THEN
  nab = ab - half
  gab = one
  DO i = 1, nab
    gab = gab * i
  END DO
ELSE
  nab = ab
  gab = rp
  DO i = 1, nab
    gab = gab * (i - half)
  END DO
END IF
q0 = (one - rho2) ** ab
coeff = gab / ga / gb

!        Use 0.5 as a starting value for Newton's iterations

y = half

!        Perform Newton's iterations

DO iter = 1, itrmax
  
!        Evaluate the first terms of the series for CDF (distribution
!        function) and PDF (density)
  
  n = one
  yp = one - y
  t = coeff * y ** a * yp ** b
  s = a * t / y / yp
  q = q0
  v = q
  cdf = v * t
  pdf = q * s
  
!        Check if a + n > (a + b + n)y
  
  80 IF (a + n > (a + b + n) * y) GO TO 90
  
!        Evaluate the next terms of two series and then the
!        partial sums
  
  q = q * (a + b + n - one) * rho2 / n
  v = v + q
  s = t * (a + b + n - one) / yp
  t = t * y * (a + b + n - one) / (a + n)
  cdf = cdf + v * t
  pdf = pdf + q * s
  n = n + one
  GO TO 80
  
!        Find the error bounds and check for convergence for both
!        series
  
  90 bndcdf = t * y * (a + b + n - one) / (a + n - (a + b + n) * y)
  bndpdf = t * (a + b + n - one) * (one - v) / yp
  100 IF (bndcdf <= eps .AND. bndpdf <= eps) GO TO 110
  
!        Continue to update the terms and then accumulate
  
  q = q * (a + b + n - one) * rho2 / n
  v = v + q
  IF (bndcdf <= eps) THEN
    s = s * y * (a + b + n - one) / (a + n - one)
    pdf = pdf + q * s
    n = n + one
    bndpdf = s * y * (a + b + n - one) * (one - v) / (a + n - one)
    GO TO 100
  ELSE IF (bndpdf <= eps) THEN
    t = t * y * (a + b + n - one) / (a + n)
    cdf = cdf + v * t
    n = n + one
    bndcdf = t * y * (a+b+n-one) / (a+n - (a+b+n)*y)
    GO TO 100
  ELSE
    s = t * (a + b + n - one) / yp
    t = t * y * (a + b + n - one) / (a + n)
    cdf = cdf + v * t
    pdf = pdf + q * s
    n = n + one
    GO TO 90
  END IF
  
!        Obtain a new Y and make changes if it is illegal
  
  110 diff = (cdf - p) / pdf
  ynew = y - diff
  IF (ynew <= zero) THEN
    y = y / two
  ELSE IF (ynew >= one) THEN
    y = (y + one) / two
  ELSE
    y = ynew
  END IF
  
!        Check for convergence of Newton's iterations
  
  IF (ABS(diff) <= delta * y) THEN
    quantile = y
    RETURN
  END IF
END DO
ifault = 1
RETURN
END FUNCTION qr2

END MODULE rsquared



PROGRAM main

!        This is a driver program that calls QR2 and produces output

USE rsquared
IMPLICIT NONE

INTEGER :: m, size, icode, k
REAL    :: rho2, p, y

10 WRITE (*, 11)
11 FORMAT (/' ENTER M (>1), N (>M), RHO2 (BETWEEN 0 AND 1), and',  &
           /' P (BETWEEN 0 AND 1) FOR QR2 ==> ')
READ (*,*) m, size, rho2, p
y = qr2(m, size, rho2, p)
icode = ifault + 1
SELECT CASE (icode)
  CASE (1)
    WRITE (*, 21) m, size, rho2, p, y
    21 FORMAT (/' qr2(', i3, ',', i3, ',', g13.4, ',', g13.4, ') = ', g13.5)
  CASE (2)
    WRITE (*, 31)
    31 FORMAT (/' NO convergence after 10 Newton ITERATIONS')
  CASE (3)
    WRITE (*, 41)
    41 FORMAT (/' THE INPUT VALUE IS ILLEGAL!')
END SELECT
WRITE (*, 61)
61 FORMAT (/ ' ENTER 1 TO CONTINUE OR 0 TO QUIT ==> ')
READ (*,*) k
IF (k == 1) GO TO 10

STOP
END PROGRAM main
