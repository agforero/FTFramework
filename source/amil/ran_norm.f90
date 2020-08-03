PROGRAM time_random_normal
IMPLICIT NONE

INTERFACE
  FUNCTION random_normal() RESULT (ran_norm)
    IMPLICIT NONE
    REAL :: ran_norm
  END FUNCTION random_normal

  SUBROUTINE gasdev_s(harvest)
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
  END SUBROUTINE gasdev_s
END INTERFACE

INTEGER              :: i, k, start(8), finish(8)
REAL                 :: x, time_taken
INTEGER, ALLOCATABLE :: seed(:)

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

!     Generate 10 million random normals using TOMS 712

CALL DATE_AND_TIME( VALUES=start )
DO i = 1, 10000000
  x = random_normal()
END DO
CALL DATE_AND_TIME( VALUES=finish )
time_taken = 3600.*(finish(5) - start(5)) + 60.*(finish(6) - start(6)) +  &
             (finish(7) - start(7)) + 0.001*(finish(8) - start(8))
WRITE(*, '(a, f8.2, a)') ' Time taken using TOMS 712 = ', time_taken, ' sec.'

!     Generate 10 million random normals using NR

CALL DATE_AND_TIME( VALUES=start )
DO i = 1, 10000000
  CALL gasdev_s(x)
END DO
CALL DATE_AND_TIME( VALUES=finish )
time_taken = 3600.*(finish(5) - start(5)) + 60.*(finish(6) - start(6)) +  &
             (finish(7) - start(7)) + 0.001*(finish(8) - start(8))
WRITE(*, '(a, f8.2, a)') ' Time taken using Num. Rec.= ', time_taken, ' sec.'

STOP
END PROGRAM time_random_normal



FUNCTION random_normal() RESULT (ran_norm)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.   This version uses the default
!  uniform random number generator which is in your fortran library.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

!  Fortran 90 version by Alan Miller (alan @ mel.dms.csiro.au)

IMPLICIT NONE
REAL :: ran_norm

!     Local variables
REAL, PARAMETER :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,   &
                   half = 0.5, r1 = 0.27597, r2 = 0.27846
REAL            :: u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
ran_norm = v/u
RETURN

END FUNCTION random_normal



SUBROUTINE gasdev_s(harvest)
! Numerical Recipes routine for generating a single normal random deviate,
! adapted to use the compiler's random number generator.

IMPLICIT NONE
REAL, INTENT(OUT) :: harvest

! Local variables
REAL          :: rsq, v1, v2
REAL, SAVE    :: g
LOGICAL, SAVE :: gaus_stored = .false.

IF (gaus_stored) THEN
   harvest = g
   gaus_stored = .false.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1 = 2.0*v1 - 1.0
      v2 = 2.0*v2 - 1.0
      rsq = v1**2 + v2**2
      if (rsq > 0.0 .and. rsq < 1.0) EXIT
   END DO
   rsq = SQRT(-2.0*LOG(rsq)/rsq)
   harvest = v1*rsq
   g = v2*rsq
   gaus_stored = .true.
END IF

RETURN
END SUBROUTINE gasdev_s
