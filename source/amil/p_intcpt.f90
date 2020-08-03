SUBROUTINE Poly_Intercept (a, b, x, y, n, u, v, m, num, ierr)
!-----------------------------------------------------------------------

!                 INTERSECTION OF A STRAIGHT LINE
!                       AND POLYGONAL PATH

!-----------------------------------------------------------------------

! The polygon is defined by the set of points (xi, yi), i = 1, 2, ..., n.
! The straight line is from (a1,a2) to (b1,b2).
! On exit, the arrays U and V contain the num points at which the line
! crosses the polygon in order, provided that num <= m.
! Error indicator:
! ierr = 0 no error detected
!      = 1 if a = b
!      = 2 U and V require more storage, i.e. num > m.
!      = -i if the ith segment of the polygon is coincident with part of the
!           line.

! Based upon routine PFIND from the NSWC Mathematics Library.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:01

IMPLICIT NONE

REAL, INTENT(IN)      :: a(2)
REAL, INTENT(IN)      :: b(2)
REAL, INTENT(IN)      :: x(:)
REAL, INTENT(IN)      :: y(:)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: u(:)
REAL, INTENT(OUT)     :: v(:)
INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(OUT)  :: num
INTEGER, INTENT(OUT)  :: ierr

! Local variables

INTEGER  :: i, ind, nm1
REAL     :: d, diff, diff1, eps, h, hi, k, ki, onem, onep, p, q, s,  &
            t, tmax, tmin, tol, tol0
!----------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0 .

eps = EPSILON(1.0)

!----------------------
num = 0
IF (n < 2) GO TO 200
h = b(1) - a(1)
k = b(2) - a(2)
IF (h == 0.0 .AND. k == 0.0) GO TO 200

ierr = 0
nm1 = n - 1
tol = 4.0*eps
tol0 = 2.0*eps
onep = 1.0 + tol
onem = 0.5 + (0.5 - tol0)

ind = 0
DO  i = 1, nm1
  hi = x(i + 1) - x(i)
  ki = y(i + 1) - y(i)
  IF (hi == 0.0 .AND. ki == 0.0) CYCLE
  ind = 1
  
!           CHECK IF THE LINE FROM A TO B AND THE I-TH
!                 LINE IN THE PATH ARE PARALLEL
  
  s = hi*k
  t = h*ki
  d = s - t
  IF (ABS(d) <= tol*MAX(ABS(s), ABS(t))) GO TO 40
!-----------------------------------------------------------------------
!                   THE LINES ARE NOT PARALLEL
!-----------------------------------------------------------------------
  p = x(i) - a(1)
  q = y(i) - a(2)
  s = hi*q
  t = ki*p
  diff = s - t
  IF (ABS(diff) <= tol*MAX(ABS(s),ABS(t))) diff = 0.0
  s = h*q
  t = k*p
  diff1 = s - t
  IF (ABS(diff1) <= tol*MAX(ABS(s),ABS(t))) diff1 = 0.0
  
  s = diff/d
  t = diff1/d
  IF (s < 0.0 .OR. s > onep) CYCLE
  IF (t < 0.0 .OR. t > onep) CYCLE
  IF (num > 0 .AND. t == 0.0) CYCLE
  IF (s > 0.0) GO TO 20
  
!                   POINT A IS ON THE I-TH LINE
  
  10 num = num + 1
  IF (num > m) GO TO 210
  u(num) = a(1)
  v(num) = a(2)
  CYCLE
  
!                   POINT B IS ON THE I-TH LINE
  
  20 IF (s < onem) GO TO 30
  21 num = num + 1
  IF (num > m) GO TO 210
  u(num) = b(1)
  v(num) = b(2)
  CYCLE
  
!              THE INTERIOR OF THE LINE FROM A TO B
!                 INTERSECTS WITH THE I-TH LINE
  
  30 num = num + 1
  IF (num > m) GO TO 210
  u(num) = a(1) + s*h
  v(num) = a(2) + s*k
  CYCLE
!-----------------------------------------------------------------------
!                     THE LINES ARE PARALLEL
!-----------------------------------------------------------------------
  40 IF (ABS(hi) > ABS(ki)) GO TO 50
  
  d = a(2) - y(i)
  IF (ABS(d) <= tol0*MAX(ABS(a(2)),ABS(y(i)))) d = 0.0
  s = d/ki
  
  p = x(i) + s*hi
  IF (ABS(a(1) - p) > tol*MAX(ABS(a(1)),ABS(p))) CYCLE
  
  d = b(2) - y(i)
  IF (ABS(d) <= tol0*MAX(ABS(b(2)),ABS(y(i)))) d = 0.0
  t = d/ki
  GO TO 60
  
  50 d = a(1) - x(i)
  IF (ABS(d) <= tol0*MAX(ABS(a(1)),ABS(x(i)))) d = 0.0
  s = d/hi
  
  p = y(i) + s*ki
  IF (ABS(p - a(2)) > tol*MAX(ABS(p),ABS(a(2)))) CYCLE
  
  d = b(1) - x(i)
  IF (ABS(d) <= tol0*MAX(ABS(b(1)),ABS(x(i)))) d = 0.0
  t = d/hi
  
!              THE 2 LINES ARE PORTIONS OF THE SAME
!                     STRAIGHT INFINITE LINE
  
  60 IF (s > 0.0 .AND. s < onem) GO TO 220
  IF (t > 0.0 .AND. t < onem) GO TO 220
  tmin = MIN(s,t)
  tmax = MAX(s,t)
  IF (tmax <= 0.0) GO TO 70
  IF (tmin >= onem) GO TO 80
  GO TO 220
  
  70 IF (tmax < 0.0) CYCLE
  IF (num > 0) CYCLE
  IF (tmax == s) GO TO 10
  GO TO 21
  
  80 IF (tmin > 1.0) CYCLE
  IF (tmin == s) GO TO 10
  GO TO 21
  
END DO
IF (ind == 0) GO TO 200

IF (num < 2) RETURN
IF (u(num) == x(1) .AND. v(num) == y(1)) num = num - 1
RETURN

!                          ERROR RETURN

200 ierr = 1
RETURN

210 ierr = 2
num = num - 1
RETURN

220 ierr = -i
RETURN
END SUBROUTINE Poly_Intercept
