PROGRAM se_logistic
! Generate artificial data to see how much bias there is in
! the estimated standard errors.

USE Logistic_Regression
IMPLICIT NONE

REAL (dp)  :: chisq, devnce, beta(0:1), se_beta(0:1), p, z(10), dev,  &
              aver(0:1) = 0.0_dp, sumbb(0:1) = 0.0_dp, aver_se(0:1) = 0.0_dp
INTEGER    :: k, i, case, n(5), s(5), ngroups = 5, ndf, ier
INTEGER, ALLOCATABLE  :: seed(:)
REAL (dp), PARAMETER  :: x(5,1) = RESHAPE(  &
                                  (/ -2.1972_dp, -0.8473_dp, 0.0_dp,  &
                                      0.8473_dp,  2.1972_dp /), (/ 5, 1 /) )

CALL RANDOM_SEED(SIZE=k)
ALLOCATE( seed(k) )

WRITE(*, '(a, i2, a)')' Enter ', k, ' integers for random no. seeds: '
READ(*, *) seed
WRITE(*, '(a, (7i10))') ' Random no. seeds: ', seed
CALL RANDOM_SEED(PUT=seed)

n = 10
DO case = 1, 10000
!     Generate binomial samples from samples of size 10,
!     with p = 0.1, 0.3, 0.5, 0.7 & 0.9
  DO i = 1, 5
    p = 0.1_dp * (2*i - 1)
    CALL RANDOM_NUMBER( z )
    s(i) = COUNT( z < p )
  END DO

  CALL logistic(ngroups, x, 1, s, n, chisq, devnce, ndf, beta, se_beta, ier)

!     Update the average
  dev = beta(0) - aver(0)
  aver(0) = aver(0) + dev/case
  sumbb(0) = sumbb(0) + dev*(beta(0) - aver(0))
  dev = beta(1) - aver(1)
  aver(1) = aver(1) + dev/case
  sumbb(1) = sumbb(1) + dev*(beta(1) - aver(1))
  aver_se = aver_se + (se_beta - aver_se) / case
END DO

sumbb = SQRT( sumbb / 9999. )
WRITE(*, '(a, 2f9.5)') ' Average estimates of beta = ', aver
WRITE(*, '(a, 2f9.5)') ' Std. devn. of sample beta = ', sumbb
WRITE(*, '(a, 2f9.5)') ' Average se_beta           = ', aver_se

STOP
END PROGRAM se_logistic
