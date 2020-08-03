FUNCTION student(t, ndf, normal) RESULT(fn_val)

!  Calculates the two-tailed probability of Student's t.
!  The user must supply a function `normal' to calculate the area under the
!  normal distribution curve, such as Applied Statistics algorithm AS 66.

!  Translation from Algol by Alan Miller, CSIRO Division of Mathematics
!  & Statistics, Clayton, Victoria 3169, Australia of:

!  Algorithm 395: Student's t-distribution by G.W. Hill
!  Comm. A.C.M., vol.13(10), 617-619, October 1970

!  N.B. The number of degrees of freedom is not necessarily an integer

!  Arguments:
!  t       double precision    Value of Student's t
!  ndf     double precision    Number of degrees of freedom
!  normal  double precision    External function to return the
!                              area under the standard normal curve
!  ier     integer             Error indicator = 0 for normal exit
!                                              = 1 if ndf < 1

!  ELF90-compatible version by Alan Miller
!  Latest revision - 7 August 1997

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN) :: t, ndf
REAL (dp)             :: fn_val

INTERFACE
  FUNCTION normal(x) RESULT(fx)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION normal
END INTERFACE

!     Local variables

REAL (dp)            :: a, b, n, t2, y, z
REAL (dp), PARAMETER :: eps = 1.d-12, one = 1.d0, twoonpi = 0.63661977236758D0
INTEGER              :: j

IF (ndf < one) THEN
  WRITE(*, *) 'Error in student: Number of degrees of freedom < 1'
  RETURN
END IF

t2 = t*t
y = t2/ndf
b = one + y
z = one

IF ( (ndf - INT(ndf)) > eps .OR. (ndf > 20.d0 - eps .AND. ndf > t2)  &
     .OR.  ndf > 200.d0 ) THEN
  
!     Asymptotic series for large or non-integer ndf
  
  IF (y > 1.d-06) y = LOG(b)
  a = ndf - 0.5D0
  b = 48.d0 * a**2
  y = a * y
  y = (((((-0.4D0*y - 3.3D0)*y - 24.d0)*y - 85.5D0) /  &
      (0.8D0*y**2 + 100.d0 + b) + y +3.d0) / b + one) * SQRT(y)
  fn_val = 2.d0 * normal(-y)
  RETURN
  
ELSE IF ( ndf < 20.d0 .AND. t2 < 4.d0) THEN
  
!     Nested summation of cosine series
  
  y = SQRT(y)
  a = y
  IF (ndf == one) a = 0.d0
  n = ndf
  
!     Tail series expansion for large t-values
  
ELSE
  a = SQRT(b)
  y = a * ndf
  j = 0
  10   j = j + 2
  IF (ABS(a-z) > eps) THEN
    z = a
    y = y * (j - 1) / (b * j)
    a = a + y / (ndf + j)
    GO TO 10
  END IF
  n = ndf + 2.d0
  z = 0.d0
  y = 0.d0
  a = - a
END IF

!     This is the `loop' from the Algol code.
!     It required jumping between different parts of an IF () THEN
!     ELSE IF () .. block.

20 n = n - 2.d0
IF (n > 1.d0) THEN
  a = (n - 1.d0) * a / (b * n) + y
  GO TO 20
END IF

IF ( ABS(n) < eps ) THEN
  a = a / SQRT(b)
ELSE
  a = (ATAN(y) + a/b) * twoonpi
END IF
fn_val = z - a

RETURN
END FUNCTION student
