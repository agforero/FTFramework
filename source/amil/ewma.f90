MODULE Update_EWMA

! Update exponentially weighted moving averages (EWMA) & variances.
! New average = lambda.x + (1 - lambda).(old average)
! The variance of the estimated mean = [lambda/(2 - lambda)].variance(x)

! Code by Alan Miller
! amiller@bigpond.net.au
! Web sites:
! http://www.ozemail.com.au/~milleraj/
! http://users.bigpond.net.au/amiller

! Latest revision - 24 October 2001

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), SAVE, PRIVATE  :: lambda, one_lambda, deg_freedom
LOGICAL, SAVE, PRIVATE    :: initialized = .FALSE.

CONTAINS


SUBROUTINE initialize(l)

REAL (dp), INTENT(IN)  :: l

lambda = l
initialized = .TRUE.
one_lambda = 1.0_dp - lambda
deg_freedom = 2.0_dp * one_lambda / (lambda*(2.0_dp - lambda))

RETURN
END SUBROUTINE initialize



SUBROUTINE update(x, average, sumsq)

! New average = lambda.x + (1 - lambda).(old average)

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN OUT)  :: average, sumsq

! Local variable
REAL (dp)  :: dev

! For speed, the following lines are not used
! IF (.NOT. initialized) THEN
!   WRITE(*, *) ' *** EWMA not initialized ***'
!   STOP
! END IF

dev = x - average
average = average + lambda * dev
sumsq = one_lambda * sumsq + dev * (x - average)

RETURN
END SUBROUTINE update



FUNCTION variance(sumsq) RESULT(fn_val)

! Calculate the variance from the sum of squares.

REAL (dp), INTENT(IN)  :: sumsq
REAL (dp)              :: fn_val

fn_val = sumsq / deg_freedom

RETURN
END FUNCTION variance

END MODULE Update_EWMA



PROGRAM test_EWMA

! Generate random normal values with mean = 7.0 and variance = 1.0
! Calculate and output the weighted average and variance after every 10 values.

USE Update_EWMA
IMPLICIT NONE

REAL (dp)  :: lambda, x, average, sumsq
INTEGER    :: i

WRITE(*, '(a)', ADVANCE='NO') ' Enter lambda between 0 and 1: '
READ(*, *) lambda
CALL initialize(lambda)

! Set arbitrary starting values for average & sumsq
average = 5.0_dp
sumsq = 1.0_dp

WRITE(*, '(a, f6.2)') ' Lambda = ', lambda
WRITE(*, *) '  I     Avge.  Variance'
DO i = 1, 1000
  x = rnorm() + 7.0
  CALL update(x, average, sumsq)
  IF (i == 10*(i/10)) THEN
    WRITE(*, '(i5, f7.1, f7.1)') i, average, variance(sumsq)
  END IF
END DO

STOP


CONTAINS


FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
REAL  :: fn_val

! Local variables

REAL            :: u, sum
REAL, SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
  fn_val = u*sln
END IF

RETURN
END FUNCTION rnorm

END PROGRAM test_EWMA
