FUNCTION chi2nc(x, f, theta) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-16  Time: 00:34:40

! ALGORITHM AS 275 APPL.STATIST. (1992), VOL.41, NO.2

! Computes the noncentral chi-square distribution function with positive
! real degrees of freedom f and nonnegative noncentrality parameter theta

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp), INTENT(IN)  :: f
REAL (dp), INTENT(IN)  :: theta
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: ifault
LOGICAL    :: flag
REAL (dp)  :: lam, n, u, v, x2, f2, t, term, bound

! EXTERNAL alogam

INTEGER, PARAMETER    :: itrmax = 50
REAL (dp), PARAMETER  :: errmax = 1.0E-12, zero = 0.0_dp, one = 1.0_dp,  &
                         two = 2.0_dp

fn_val = x
ifault = 2
IF (f <= zero .OR. theta < zero) GO TO 999
ifault = 3
IF (x < zero) GO TO 999
ifault = 0
IF (x == zero) RETURN
lam = theta / two

!       Evaluate the first term

n = one
u = EXP(-lam)
v = u
x2 = x / two
f2 = f / two
t = x2 ** f2 * EXP(-x2) / EXP(lngamma(f2 + one))
term = v * t
fn_val = term

!       Check if (f+2n) is greater than x

flag = .false.
10 IF (f + two * n - x <= zero) GO TO 30

!       Find the error bound and check for convergence

flag = .true.
20 bound = t * x / (f + two * n - x)
IF (bound > errmax .AND. INT(n) <= itrmax) GO TO 30
IF (bound > errmax) ifault = 1
GO TO 999

!       Evaluate the next term of the expansion and then the partial sum

30 u = u * lam / n
v = v + u
t = t * x / (f + two * n)
term = v * t
fn_val = fn_val + term
n = n + one
IF (flag) GO TO 20
GO TO 10

999 WRITE(*, *) 'Error in CHI2NC, IFAULT = ', ifault
RETURN


CONTAINS


FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

REAL (dp), INTENT(IN)  :: z
REAL (dp)              :: lanczos

! Local variables

REAL (dp), PARAMETER   :: a(9) = (/   &
                          0.9999999999995183_dp, 676.5203681218835_dp, &
                          -1259.139216722289_dp, 771.3234287757674_dp, &
                          -176.6150291498386_dp, 12.50734324009056_dp, &
                          -0.1385710331296526_dp, 0.9934937113930748D-05, &
                          0.1659470187408462D-06 /),  &
                          lnsqrt2pi = 0.9189385332046727_dp, half = 0.5_dp, &
                          sixpt5 = 6.5_dp, seven = 7.0_dp
INTEGER    :: j
REAL (dp)  :: tmp

IF (z <= zero) THEN
  WRITE(*, *) 'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma


END FUNCTION chi2nc
