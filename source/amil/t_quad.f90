PROGRAM test_quad
!     Three tests of the basic quadruple precision module.
!     1.       (A + B)^2          = (A - B)^2 + 4.A.B
!     2.  (A^2 - B^2) / (A - B)   = A + B
!     3.  SQRT(A^2 + 2.A.B + B^2) = A + B

USE quadruple_precision
IMPLICIT NONE

TYPE (quad)          :: a, b, lhs, rhs, diff
REAL (dp)            :: half = 0.5_dp, small = EPSILON(half)
INTEGER, ALLOCATABLE :: seed(:)
INTEGER              :: k, i

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

DO i = 1, 10
  CALL RANDOM_NUMBER(a%hi)
  a%hi = (a%hi - half) / a%hi
  a%lo = a%hi * small
  CALL RANDOM_NUMBER(b%hi)
  b%hi = (b%hi - half) / b%hi
  b%lo = b%hi * small
                             ! (A + B)^2 = (A - B)^2 + 4.A.B
  lhs = (a + b) * (a + b)
  rhs = a * b
  rhs = 4._dp * rhs
  rhs = (a - b) * (a - b) + rhs
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! (A^2 - B^2) / (A - B) = (A + B)
  lhs = (a*a - b*b) / (a - b)
  rhs = a + b
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! SQRT(A^2 + 2.A.B + B^2) = A + B
  lhs = a * b
  lhs = 2._dp * lhs
  lhs = SQRT(a*a + lhs + b*b)
  IF (rhs%hi < 0._dp) THEN             ! Force +ve square root
    rhs = -rhs
  END IF
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi

END DO

STOP
END PROGRAM test_quad
