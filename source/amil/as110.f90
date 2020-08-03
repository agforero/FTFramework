SUBROUTINE lpest(n, p, x, y, maxit, a, b, sd, r, rate, it, npo, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-17  Time: 15:07:09

! FORMAL PARAMETERS
!    N     INTEGER        input : the number of points
!    P     REAL           input : p in the Lp norm
!    X     REAL ARRAY (N) input : the observed values x(i)
!    Y     REAL ARRAY (N) input : the observed values y(i)
!    MAXIT INTEGER        input : maximum allowable number of iterations
!    A     REAL           output : the estimate of alfa
!    B     REAL           output : the estimate of beta
!    SD    REAL           output : the Lp norm
!    R     REAL ARRAY (N) output : the signed residuals
!    RATE  REAL           output : abs(S(k+1)-S(k))/S(k+1), where S(k) is
!                                  the Lp norm of the fit on the kth iteration
!    IT    INTEGER        output : the number of iterations
!    NPO   INTEGER        output : the number of the points on the line
!   IFAULT INTEGER        output :
!                 IFAULT = 0 if the routine converged.
!                        = 1 if return was due to an increase in the norm.
!                        = 2 if the maximum iterations
!                            specified was less than 2.
!                        = 3 if the weighted sample variance of x is 0.
!                        = 4 if N given less than 2.
!                        = 5 if convergence has not be achieved within
!                            the maximum number of iterations specified.

!  ALGORITHM AS 110 APPL. STATIST. (1977), VOL.26, NO.1

!  Lp-NORM FIT OF STRAIGHT LINE BY EXTENSION OF SCHLOSSMACHER,
!  particularly for 1 <= p <= 2.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: p
REAL, INTENT(IN)      :: x(n)
REAL, INTENT(IN)      :: y(n)
INTEGER, INTENT(IN)   :: maxit
REAL, INTENT(OUT)     :: a
REAL, INTENT(OUT)     :: b
REAL, INTENT(OUT)     :: sd
REAL, INTENT(OUT)     :: r(n)
REAL, INTENT(OUT)     :: rate
INTEGER, INTENT(OUT)  :: it
INTEGER, INTENT(OUT)  :: npo
INTEGER, INTENT(OUT)  :: ifault

INTEGER  :: i, isw
REAL     :: absri, a2, b2, div, dx, dxy, eps2,   &
    res, sd2, sw, ssx, spxy, w, xiw, xmean, wp, xi, yi, ymean
REAL, PARAMETER  :: eps = 1.0E-6

IF(maxit < 2) GO TO 9
IF(n < 2) GO TO 10
ifault = 0
wp = p - 2.0
eps2 = 2.0 * eps
sd = 0.0
r(1:n) = 1.0

DO  it = 1, maxit
  npo = 0
  
!          CALCULATE A AND B BY LEAST SQUARES ON WEIGHTED DATA,
!          USING THE HERRAMAN ALGORITHM.
!          OMIT OBSERVATIONS WITH SMALL RESIDUALS.
  
  sw = 0.0
  xmean = 0.0
  ymean = 0.0
  ssx = 0.0
  spxy = 0.0
  DO  i = 1, n
    absri = ABS(r(i))
    IF(absri <= eps) GO TO 2
    w = absri ** wp
    sw = sw + w
    div = w / sw
    xi = x(i) - xmean
    yi = y(i) - ymean
    xiw = xi * w
    dx = xi * xiw
    dxy = yi * xiw
    ssx = ssx + dx - dx * div
    spxy = spxy + dxy - dxy * div
    xmean = xmean + xi * div
    ymean = ymean + yi * div
    CYCLE
    2 npo = npo + 1
  END DO
  IF(ssx < eps) GO TO 11
  b = spxy / ssx
  a = ymean - b * xmean
  
!          FORM RESIDUALS AND TEST CONVERGENCE
  
  sd2 = 0.0
  isw = 0
  DO  i = 1, n
    res = y(i) - a - b * x(i)
    absri = ABS(res)
    IF(ABS(absri - ABS(r(i))) > eps2) isw = 1
    sd2 = sd2 + absri ** p
    r(i) = res
  END DO
  rate = ABS(sd2 - sd) / sd2
  IF(isw == 0) RETURN
  IF(it == 1) GO TO 5
  
!        TEST FOR INCREASE IN NORM
  
  IF(sd2 > sd) GO TO 7
  5 sd = sd2
  a2 = a
  b2 = b
END DO
!     FAILED TO CONVERGE IN MAXIT ITERATIONS

ifault = 5
RETURN

!      NORM INCREASED, RESTORE A, B, AND R, THEN RETURN

7 ifault = 1
a = a2
b = b2
DO  i = 1, n
  r(i) = y(i) - a - b * x(i)
END DO
RETURN

!      MAXIT SPECIFIED LESS THAN 2

9 ifault = 2
RETURN

!      N LESS THAN 2

10 ifault = 4
RETURN

!      VARIANCE OF WEIGHTED X IS ZERO

11 ifault = 3
RETURN
END SUBROUTINE lpest
