PROGRAM test_longlog_longexp
! Tests:
! 1. log(x^2) = 2.*log(x)
! 2. e^(2x)   = (e^x)**2

USE quadruple_precision
IMPLICIT NONE

TYPE (quad)          :: x, lhs, rhs, diff
INTEGER, ALLOCATABLE :: seed(:)
INTEGER              :: k, i
REAL (dp)            :: half = 0.5_dp, small = 1.e-16_dp

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

DO i = 1, 10
  CALL RANDOM_NUMBER(x%hi)
  x%hi = (x%hi - half) / x%hi
  x%lo = x%hi * small
                                       ! log(x^2) = 2.*log(x)
  lhs = LOG(x * x)
  IF (x%hi > 0._dp) THEN
    rhs = LOG(x)
  ELSE
    rhs = LOG(-x)
  END IF
  rhs = 2._dp * rhs
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                                       ! e^(2x) = (e^x)**2
  lhs = 2._dp * x
  lhs = EXP(lhs)
  rhs = EXP(x)
  rhs = rhs * rhs
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
END DO

STOP
END PROGRAM test_longlog_longexp
