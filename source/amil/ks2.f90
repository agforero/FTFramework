MODULE Kolmogorov_Smirnov
! Calculate 2-tail probabilities for the max abs deviation of the sample
! distribution function from the hypothetical.
! Includes KS1 which calculates exact single-tail probabilities.
! N = sample size
! D = max. abs. deviation

! N.B. The hypothesized distribution function MUST NOT use ANY parameters
!      estimated from the data used to form the empirical distribution.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

PRIVATE
PUBLIC  :: ks1, ks2

CONTAINS


FUNCTION KS2(n, d) RESULT( prob )

! Uses a combination of 3 different algorithms for different N & D.

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: d
REAL (dp)              :: prob

IF (n <= 50) THEN
  prob = 1.0_dp - pks2(n, d)
ELSE IF (n > 500) THEN
  prob = KS2_asympt(n, d)
ELSE
  IF (d <= 1.2_dp / SQRT( DBLE(n) )) THEN
    IF (n <= 100) THEN
      prob = 1.0_dp - pks2(n, d)
    ELSE
      prob = KS2_asympt(n, d)
    END IF
  ELSE
    prob = 2. * KS1(n, d)
  END IF
END IF

RETURN
END FUNCTION KS2



! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-28  Time: 18:34:23

!     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12, P. 703.

FUNCTION pks2(n, d) RESULT( prob )

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: d
REAL (dp)              :: prob

! N IS THE SAMPLE SIZE USED.

! D IS THE MAXIMUM MAGNITUDE (OF THE DISCREPANCY BETWEEN THE EMPIRICAL AND
! PROPOSED DISTRIBUTIONS) IN EITHER THE POSITIVE OR NEGATIVE DIRECTION.
! PKS2 IS THE EXACT PROBABILITY OF OBTAINING A DEVIATION NO LARGER THAN D.
! THESE FORMULAS APPEAR AS (23) AND (24) IN J. DURBIN.  THE PROBABILITY THAT
! THE SAMPLE DISTRIBUTION FUNCTION LIES BETWEEN TWO PARALLEL STRAIGHT LINES.
! ANNALS OF MATHEMATICAL STATISTICS 39, 2(APRIL 1968), 398-411.

REAL (dp)               :: sum, ci, fn, fnd, ft, fu, fv, sign
INTEGER                 :: i, j, jmax, k, nd, ndd, nddp, ndp, ndt
REAL (dp), ALLOCATABLE  :: q(:), fact(:)
REAL (dp), PARAMETER    :: zero = 0.0_dp, one = 1.0_dp

IF (n == 1) GO TO 90
fn = n
fnd = fn*d
ndt = 2.*fnd
IF (ndt < 1) GO TO 100

ALLOCATE( q(n+1), fact(n+1) )
nd = fnd
ndd = MIN(2*nd,n)
ndp = nd + 1
nddp = ndd + 1
fact(1) = one
ci = one
DO  i=1,n
  fact(i+1) = fact(i)*ci
  ci = ci + one
END DO
q(1) = one
IF (ndd == 0) GO TO 50

ci = one
DO  i=1,ndd
  q(i+1) = ci**i/fact(i+1)
  ci = ci + one
END DO
IF (ndp > n) GO TO 80

fv = ndp - fnd
jmax = fv + 1
DO  i=ndp,ndd
  sum = zero
  ft = fnd
  k = i
  fu = fv
  DO  j=1,jmax
    sum = sum + ft**(j-2)/fact(j)*fu**k/ fact(k+1)
    ft = ft + one
    fu = fu - one
    k = k - 1
  END DO
  q(i+1) = q(i+1) - 2.*fnd*sum
  jmax = jmax + 1
  fv = fv + one
END DO
IF (ndd == n) GO TO 80

50 DO  i=nddp,n
  sum = zero
  SIGN = one
  ft = 2.*fnd
  DO  j=1,ndt
    ft = ft - one
    k = i - j + 1
    sum = sum + SIGN*ft**j/fact(j+1)*q(k)
    SIGN = -SIGN
  END DO
  q(i+1) = sum
END DO

80 prob = q(n+1)*fact(n+1) / fn**n
DEALLOCATE( q, fact )
RETURN

90 prob = 2.*d - one
RETURN

100 prob = zero
RETURN
END FUNCTION pks2



FUNCTION KS1(n, d) RESULT( prob )
! Calculates the single-tailed probability for the on sample KS test.

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: d
REAL (dp)              :: prob

! Local variables
REAL (dp)  :: eps, lncj, lncj0, term
INTEGER    :: j, j0, jupper
REAL (dp), PARAMETER  :: one = 1.0_dp

eps = 10.0 * EPSILON(one)
jupper = n * (one - d)

! Evaluate one of the largest terms in the sum,
! and then sum up & down from term j0 until terms are negligible

j0 = MIN(jupper, n/2)
lncj0 = log_nCr(n, j0)
prob = EXP( lncj0 + (n-j0)*LOG(one - d - j0/DBLE(n)) +   &
                     (j0-1)*LOG(d + j0/DBLE(n)) )

lncj = lncj0
DO j = j0+1, jupper
  lncj = lncj + LOG(DBLE(n+1-j) / DBLE(j))
  term = EXP( lncj + (n-j)*LOG(one - d - j/DBLE(n)) +   &
                     (j-1)*LOG(d + j/DBLE(n)) )
  prob = prob + term
  IF (term < eps*prob) EXIT
END DO

lncj = lncj0
DO j = j0-1, 0, -1
  lncj = lncj + LOG(DBLE(j+1) / DBLE(n-j))
  term = EXP( lncj + (n-j)*LOG(one - d - j/DBLE(n)) +   &
                     (j-1)*LOG(d + j/DBLE(n)) )
  prob = prob + term
  IF (term < eps*prob) EXIT
END DO

prob = d * prob

RETURN
END FUNCTION KS1



FUNCTION log_nCr(n, r) RESULT(fn_val)
! Evaluate log of the number of combinations of r out of n.

INTEGER, INTENT(IN)  :: n, r
REAL (dp)            :: fn_val

! Local variables

INTEGER  :: i, j

! If MIN(r, n-r) < 10, use multiplications & divisions, otherwise
! use the gamma function.

j = MIN(r, n-r)
IF (j < 10) THEN
  fn_val = 1.0_dp
  DO i = 1, j
    fn_val = fn_val * (n+1-i) / DBLE(i)
  END DO
  fn_val = LOG(fn_val)
ELSE
  fn_val = lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n+1-r))
END IF

RETURN
END FUNCTION log_nCr



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

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
                      -1259.139216722289_dp, 771.3234287757674_dp, &
                      -176.6150291498386_dp, 12.50734324009056_dp, &
                      -0.1385710331296526_dp, 0.9934937113930748D-05, &
                       0.1659470187408462D-06 /), zero = 0._dp,   &
                       one = 1._dp, lnsqrt2pi = 0.9189385332046727_dp, &
                       half = 0.5_dp, sixpt5 = 6.5_dp, seven = 7._dp, tmp
INTEGER  :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
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



FUNCTION KS2_asympt(n, d) RESULT( prob )
! Using Smirnov's asymptotic formaula for large N

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: d
REAL (dp)              :: prob

! Local variables

REAL (dp)  :: sgn, sum, term, z
INTEGER    :: i
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

z = d * SQRT( DBLE(n) )
sum = zero
i = 1
sgn = one
DO
  term = EXP( -2.0_dp * (i*z)**2 )
  sum = sum + sgn*term
  IF (term < EPSILON(one)) EXIT
  i = i + 1
  sgn = -sgn
END DO
prob = 2.0_dp*sum

RETURN
END FUNCTION KS2_asympt

END MODULE Kolmogorov_Smirnov
