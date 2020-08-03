SUBROUTINE locpt (x0, y0, x, y, n, l, m)
!-----------------------------------------------------------------------
! GIVEN A POLYGONAL LINE CONNECTING THE VERTICES (X(I),Y(I)) (I = 1,...,N)
! TAKEN IN THIS ORDER.  IT IS ASSUMED THAT THE POLYGONAL PATH IS A LOOP,
! WHERE (X(N),Y(N)) = (X(1),Y(1)) OR THERE IS AN ARC FROM (X(N),Y(N)) TO
! (X(1),Y(1)).  N.B. The polygon may cross itself any number of times.

! (X0,Y0) IS AN ARBITRARY POINT AND L AND M ARE VARIABLES.
! On output, L AND M ARE ASSIGNED THE FOLLOWING VALUES ...

!    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
!    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
!    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH

! M = 0 IF (X0,Y0) IS ON OR OUTSIDE THE PATH.  IF (X0,Y0) IS INSIDE THE
! PATH THEN M IS THE WINDING NUMBER OF THE PATH AROUND THE POINT (X0,Y0).

! Fortran 66 version by A.H. Morris
! Converted to ELF90 compatibility by Alan Miller, 15 February 1997

!-----------------------

IMPLICIT NONE
REAL, INTENT(IN)     :: x0, y0, x(:), y(:)
INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: l, m

!     Local variables
INTEGER :: i, n0
REAL    :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
n0 = n
IF (x(1) == x(n) .AND. y(1) == y(n)) n0 = n - 1
pi = ATAN2(0.0, -1.0)
pi2 = 2.0*pi
tol = 4.0*eps*pi
l = -1
m = 0

u = x(1) - x0
v = y(1) - y0
IF (u == 0.0 .AND. v == 0.0) GO TO 20
IF (n0 < 2) RETURN
theta1 = ATAN2(v, u)

sum = 0.0
theta = theta1
DO i = 2, n0
  u = x(i) - x0
  v = y(i) - y0
  IF (u == 0.0 .AND. v == 0.0) GO TO 20
  thetai = ATAN2(v, u)
  
  angle = ABS(thetai - theta)
  IF (ABS(angle - pi) < tol) GO TO 20
  IF (angle > pi) angle = angle - pi2
  IF (theta > thetai) angle = -angle
  sum = sum + angle
  theta = thetai
END DO

angle = ABS(theta1 - theta)
IF (ABS(angle - pi) < tol) GO TO 20
IF (angle > pi) angle = angle - pi2
IF (theta > theta1) angle = -angle
sum = sum + angle

!     SUM = 2*PI*M WHERE M IS THE WINDING NUMBER

m = ABS(sum)/pi2 + 0.2
IF (m == 0) RETURN
l = 1
IF (sum < 0.0) m = -m
RETURN

!     (X0, Y0) IS ON THE BOUNDARY OF THE PATH

20 l = 0
RETURN
END SUBROUTINE locpt
