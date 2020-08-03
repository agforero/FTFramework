FUNCTION rngpi(t, n) RESULT(fn_val)

! N.B. Argument IFAULT has been removed.

! Code converted using TO_F90 by Alan Miller
! Date: 2003-02-09  Time: 14:24:28

!        Algorithm AS 126  Appl. Statist. (1978) Vol. 27, No. 2

!        Computes the probability of the normal range given t, the
!        upper limit of integration, and n, the sample size.

!     Auxiliary function required: ALNORM = algorithm AS66

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN)  :: t
INTEGER, INTENT(IN)    :: n
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: a, b, c, xl, y
INTEGER    :: i

REAL (dp), PARAMETER  :: g(8) = (/  &
    0.4947004675D0, 0.4722875115D0, 0.4328156012D0, 0.3777022042D0,  &
    0.3089381222D0, 0.2290083888D0, 0.1408017754D0, 0.04750625492D0 /)

REAL (dp), PARAMETER  :: h(8) = (/  &
    0.01357622971D0, 0.03112676197D0, 0.04757925584D0, 0.06231448563D0,  &
    0.07479799441D0, 0.08457825969D0, 0.09130170752D0, 0.09472530523D0 /)

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, two = 2.0_dp,  &
                         eight = 8.0_dp

INTERFACE
  FUNCTION alnorm( x, upper ) RESULT( fn_val )
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 100)
    REAL(DP), INTENT(IN)   ::  x
    LOGICAL,   INTENT(IN)  ::  upper
    REAL(DP)               ::  fn_val
  END FUNCTION alnorm
END INTERFACE

fn_val = zero
IF (t <= zero .OR. n <= 1) RETURN

xl = half * t
a = half * (eight+xl)
b = eight - xl
y = zero
DO  i = 1, 8
  c = b * g(i)
  y = y + h(i) * (risf(a+c) + risf(a-c))
END DO
fn_val = (two*(alnorm(xl,.false.)-half)) ** n + two * b * y * n
RETURN


CONTAINS

FUNCTION risf(x) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = 0.3989422804D0 * EXP(-half*x*x) * (alnorm(x, .false.) -  &
         alnorm(x-t, .false.)) ** (n-1)
RETURN
END FUNCTION risf

END FUNCTION rngpi
