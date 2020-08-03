FUNCTION ran_gamma(s, first) RESULT(fn_val)
USE random
IMPLICIT NONE

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS).

! Generates a random gamma deviate for shape parameter s >= 1.

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (s < 1.0) THEN
   WRITE(*, *) 'Shape parameter must be >= 1'
   STOP
END IF

IF (first) THEN
  d = s - 1./3.
  c = 1.0/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (1.0 + c*x)**3
    IF (v > 0.0) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1.0 - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < 0.5*x**2 + d*(1.0 - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION ran_gamma



PROGRAM test_rgamma
USE random
IMPLICIT NONE

INTERFACE
  FUNCTION ran_gamma(s, first) RESULT(fn_val)
    USE random
    IMPLICIT NONE
    REAL, INTENT(IN)    :: s
    LOGICAL, INTENT(IN) :: first
    REAL                :: fn_val
  END FUNCTION ran_gamma
END INTERFACE

REAL, PARAMETER  :: shape(8) = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0 /)
INTEGER          :: i, j
REAL             :: aver, dev, sumsq, start, finish, x

DO i = 1, 8
  WRITE(*, *)
  WRITE(*, '(a, f7.1)') ' Shape parameter = ', shape(i)

! First time the new algorithm

  CALL CPU_TIME(start)
  x = ran_gamma(shape(i), .TRUE.)
  aver = x
  sumsq = 0.0
  DO j = 2, 1000000
    x = ran_gamma(shape(i), .FALSE.)
    dev = x - aver
    aver = aver + dev/j
    sumsq = sumsq + dev*(x-aver)
  END DO
  CALL CPU_TIME(finish)
  sumsq = SQRT(sumsq/999999.)
  WRITE(*, *) 'Using Marsaglia & Tsang algorithm'
  WRITE(*, '(a, f7.2, a, 2f7.3)')  &
        ' Time = ', finish - start, 'sec.   Mean & st.devn. = ', aver, sumsq

! Repeat using algorithm in module random

  CALL CPU_TIME(start)
  x = random_gamma(shape(i), .TRUE.)
  aver = x
  sumsq = 0.0
  DO j = 2, 1000000
    x = random_gamma(shape(i), .FALSE.)
    dev = x - aver
    aver = aver + dev/j
    sumsq = sumsq + dev*(x-aver)
  END DO
  CALL CPU_TIME(finish)
  sumsq = SQRT(sumsq/999999.)
  WRITE(*, *) 'Using algorithm in module random'
  WRITE(*, '(a, f7.2, a, 2f7.3)')  &
        ' Time = ', finish - start, 'sec.   Mean & st.devn. = ', aver, sumsq
END DO

STOP
END PROGRAM test_rgamma
