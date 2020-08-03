PROGRAM test_cst
!     Test of sine, cosine & tangent
! 1. sin(2.x)  = 2.sin(x).cos(x)
! 2. cos(2.x)  = (cos(x) + sin(x)).(cos(x) - sin(x))
! 3. sec(x)**2 = 1 + tan(x)**2

USE quadruple_precision
IMPLICIT NONE

TYPE (quad)          :: x, twox, sine, cosine, tangent, lhs, rhs, diff
REAL (dp)            :: half = 0.5_dp, small = EPSILON(half), r
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

DO i = 1, 100
  CALL RANDOM_NUMBER(x%hi)
  CALL RANDOM_NUMBER(r)
  x%hi = (x%hi - half) / r
  x%lo = x%hi * small
                             ! sin(2.x)  = 2.sin(x).cos(x)
  twox = 2._dp * x
  lhs = SIN(twox)
  sine = SIN(x)
  cosine = COS(x)
  rhs = sine * cosine
  rhs = 2._dp * rhs
  diff = lhs - rhs
  WRITE(*, '(" sin(2x)   lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! cos(2.x)  = (cos(x) + sin(x)).(cos(x) - sin(x))
  lhs = COS(twox)
  rhs = (cosine + sine) * (cosine - sine)
  diff = lhs - rhs
  WRITE(*, '(" cos(2x)   lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! sec(x)**2 = 1 + tan(x)**2
  diff%hi = 1._dp
  diff%lo = 0._dp
  lhs = diff / cosine
  lhs = lhs * lhs
  tangent = TAN(x)
  rhs = diff + tangent * tangent
  diff = lhs - rhs
  WRITE(*, '(" sec(x)^2  lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
END DO

STOP
END PROGRAM test_cst
