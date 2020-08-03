MODULE Logistic_Regression
! A package for fitting linear logistic models by iteratively re-weighted
! least squares.

! The model fitted is that the probability of a `success' is:

!      p = 1 - 1/[1 + exp(b0 + b1.X1 + b2.X2 + ... + bk.Xk)]

! where X1, X2, ... , Xk are the predictor variables, and the coefficients
! b0, b1, b2, ..., bk are to be determined.

! N.B. The residual variance used here in estimating standard errors is the
!      larger of 1.0 and that from the weighted least squares calculations;
!      it is not the theoretical residual variance (1.0) assuming a binomial
!      distribution about the logistic curve.   If the fit of the logistic
!      is poor, the standard errors of the coefficients in the logistic will
!      be much larger.

! The calculation of chi-squared was corrected on 9 August 2003.
! My thanks to Zhang Guan

! By Alan Miller
! amiller @ bigpond.net.au
! users.bigpond.net.au/amiller/

! Latest revision - 9 August 2003

USE lsq
IMPLICIT NONE

PUBLIC  :: logistic

CONTAINS


SUBROUTINE logistic(ngroups, x, k, s, n, chisq, devnce, ndf, beta, se_beta,  &
                    ier, cov_beta, fit, stdres)

! Input arguments:

! ngroups     The number of groups of observations.

! x(:,:)      x(i,j) contains the value of the j-th predictor for group i.
!             The X-values are the same for all the n(i) trials in group i.

! k           The number of predictors.   N.B. A constant will be fitted in
!             the model; do not include it in the k predictor variables.

! s(:), n(:)  s(i) is the number of `successes' out of the n(i) `trials' in
!             the i-th group.   In many applications, each n(i) = 1, that is
!             each case will be treated individually, and then s(i) = 0 or 1
!             for each of the two possible outcomes.

INTEGER, INTENT(IN)    :: ngroups, k, s(:), n(:)
REAL (dp), INTENT(IN)  :: x(:,:)

! Output arguments:

! chisq       The value of chi-squared on exit, when a model has been fitted.

! devnce      The deviance on exit, when a model has been fitted.

! ndf         Number of degrees of freedom.

! beta(0:)    The fitted coefficients in the logistic model.

! se_beta(0:) Approximate standard errors of the beta coefficients.

! ier         Error indicator
!             = 0 successful termination
!             = 1 if ngroups < 2 or ndf < 0
!             = 2 if any n(i) < 0
!             = 3 if any r(i) < 0
!             = 4 if any r(i) > n(i)
!             = 5 if any X-variable is constant
!             = 6 if a singularity is detected
!             = 7 if any beta(i) is tending to +/- infinity
!             = 8 failure to converge in 20 iterations

REAL (dp), INTENT(OUT)  :: chisq, devnce, beta(0:), se_beta(0:)
INTEGER, INTENT(OUT)    :: ndf, ier

! Optional output arguments:

! cov_beta(0:,0:)     Approximate covariance matrix of the fitted coefficients.

! fit(:)      The fitted probabilities of a success for each group.

! stdres(:)   Vector of standardized residuals.

REAL (dp), INTENT(OUT), OPTIONAL  :: cov_beta(0:,0:), fit(:), stdres(:)


! Local variables

INTEGER    :: i, iter, j, ncov, pos
REAL (dp)  :: propn(ngroups), p(ngroups), wt(ngroups), xrow(0:k), db(0:k), &
              bnew(0:k), dev_new, xb, pnew(ngroups), wnew(ngroups), a, b,  &
              range(k), var, e(ngroups), hii
LOGICAL    :: lindep(0:k)
REAL (dp), ALLOCATABLE  :: covmat(:)

! Initial checks

ier = 0
ndf = ngroups - k - 1
IF (ngroups < 2 .OR. ndf < 0) THEN
  ier = 1
  RETURN
END IF
IF (ANY(n(1:ngroups) < 0)) THEN
  ier = 2
  RETURN
END IF
IF (ANY(s(1:ngroups) < 0)) THEN
  ier = 3
  RETURN
END IF
IF (ANY(s(1:ngroups) > n(1:ngroups))) THEN
  ier = 4
  RETURN
END IF

! Calculate ranges of the X-variables to use in testing whether any beta
! is tending to +/- infinity.  Also test that no variable is constant.

DO i = 1, k
  a = MAXVAL(x(1:ngroups,i))
  b = MINVAL(x(1:ngroups,i))
  range(i) = a - b
  IF (range(i) < EPSILON(0.0_dp) * (ABS(a) + ABS(b))) THEN
    ier = 5
    RETURN
  END IF
END DO

! Start with all beta's = 0 and weights = 1.

beta(0:k) = 0.0_dp
wt(1:ngroups) = 1.0_dp
p(1:ngroups) = 0.5_dp

! propn stores the sample proportions, i.e. s(i) / n(i)
propn(1:ngroups) = REAL(s(1:ngroups), KIND=dp) / n(1:ngroups)
iter = 1

! Start of iterative cycle

DO
  CALL startup(k, .TRUE.)
  DO i = 1, ngroups
    IF (iter == 1) THEN
      xrow(0) = 0.25_dp
      xrow(1:k) = 0.25_dp*x(i, 1:k)
    ELSE
      xrow(0) = p(i)*(1.0_dp - p(i))
      xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
    END IF
    CALL includ(wt(i), xrow, propn(i)-p(i))
  END DO

! Test for a singularity

  CALL sing(lindep, ier)
  IF (ier /= 0) THEN
    DO i = 1, k
      IF (lindep(i)) THEN
        WRITE(*, '(a, i6, a)') ' Variable number ', i,  &
                               ' is linearly dependent upon earlier variables'
      END IF
    END DO
    ier = 6
    RETURN
  END IF

  CALL regcf(db, k+1, ier)
  10 bnew = beta(0:k) + db

! Calculate new p(i)'s, weights & deviance

  dev_new = 0.0_dp
  DO i = 1, ngroups
    xb = DOT_PRODUCT( x(i,1:k), bnew(1:k) ) + bnew(0)
    xb = EXP(xb)
    pnew(i) = xb / (1.0_dp + xb)
    wnew(i) = REAL(n(i), KIND=dp) / (pnew(i)*(1.0_dp - pnew(i)))
    IF (iter == 1) wnew(i) = SQRT(wnew(i))
    IF (s(i) > 0) dev_new = dev_new + s(i)*LOG(propn(i)/pnew(i))
    IF (s(i) < n(i)) dev_new = dev_new +   &
                           (n(i)-s(i))*LOG((1.0_dp-propn(i))/(1.0_dp-pnew(i)))
  END DO
  dev_new = 2 * dev_new

! If deviance has increased, reduce the step size.

  IF (iter > 2) THEN
    IF (dev_new > devnce*1.0001_dp) THEN
      db = 0.5_dp * db
      GO TO 10
    END IF
  END IF

! Replace betas, weights & p's with new values

  beta(0:k) = bnew(0:k)
  wt = wnew
  p(1:ngroups) = pnew

! Test for convergence

  IF (iter > 2 .AND. devnce - dev_new < 0.0001_dp) EXIT
  devnce = dev_new
  iter = iter + 1
  IF (iter > 20) THEN
    ier = 8
    RETURN
  END IF

! Test for a very large beta

  DO i = 1, k
    IF (ABS(beta(i))*range(i) > 30.0_dp) THEN
      WRITE(*, '(a, i4, a)') ' Coefficient for variable no.', i,  &
                             ' tending to infinity'
      ier = 7
      RETURN
    END IF
  END DO

END DO

e = n(1:ngroups)*p(1:ngroups)
chisq = SUM( (s(1:ngroups) - e)**2 / (e * (1.0_dp - p(1:ngroups))) )
devnce = dev_new

! Calculate the approximate covariance matrix for the beta's, if ndf > 0.

IF (ndf > 0) THEN
  ncov = (k+1)*(k+2)/2
  ALLOCATE( covmat(ncov) )
  CALL cov(k+1, var, covmat, ncov, se_beta, ier)
  IF (var < 1.0_dp) THEN
    covmat = covmat / var
    se_beta = se_beta / SQRT(var)
  END IF

  IF(PRESENT(cov_beta)) THEN
    pos = 1
    DO i = 0, k
      cov_beta(i,i) = covmat(pos)
      pos = pos + 1
      DO j = i+1, k
        cov_beta(i,j) = covmat(pos)
        cov_beta(j,i) = covmat(pos)
        pos = pos + 1
      END DO
    END DO
  END IF
END IF

IF(PRESENT(fit)) fit(1:ngroups) = p

IF (PRESENT(stdres)) THEN
  DO i = 1, ngroups
    xrow(0) = p(i)*(1.0_dp - p(i))
    xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
    CALL hdiag(xrow, k+1, hii, ier)
    stdres(i) = (s(i) - n(i)*p(i)) /   &
                SQRT(n(i)*p(i)*(1.0_dp - p(i))*(1.0_dp - hii))
  END DO
END IF

IF (ALLOCATED(covmat)) DEALLOCATE( covmat )

RETURN
END SUBROUTINE logistic



FUNCTION chi_squared(ndf, chi2) RESULT(prob)
! Calculate the chi-squared distribution function
! ndf  = number of degrees of freedom
! chi2 = chi-squared value
! prob = probability of a chi-squared value <= chi2 (i.e. the left-hand
!        tail area)

INTEGER, INTENT(IN)    :: ndf
REAL (dp), INTENT(IN)  :: chi2
REAL (dp)              :: prob

! Local variables
REAL (dp) :: half = 0.5_dp, x, p

x = half * chi2
p = half * REAL(ndf)
prob = gammad(x, p)
RETURN

END FUNCTION chi_squared



FUNCTION gammad(x, p) RESULT(gamma_prob)

!  ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!  Computation of the Incomplete Gamma Integral

!  Auxiliary functions required: ALNORM = algorithm AS66 (included) & LNGAMMA

!  Converted to be compatible with ELF90 by Alan Miller
!  N.B. The return parameter IFAULT has been removed as ELF90 allows only
!  one output parameter from functions.   An error message is issued instead.

! This revision - 15 October 1996

REAL (dp), INTENT(IN) :: x, p
REAL (dp)             :: gamma_prob

!     Local variables

REAL (dp) :: pn1, pn2, pn3, pn4, pn5, pn6, tol = 1.d-14, oflo = 1.d+37,  &
             xbig = 1.d+8, arg, c, rn, a, b, one = 1._dp, zero = 0._dp, an, &
             two = 2._dp, elimit = -88._dp, plimit = 1000._dp, three = 3._dp, &
             nine = 9._dp

gamma_prob = zero

!      Check that we have valid values for X and P

IF (p <= zero .OR. x < zero) THEN
  WRITE(*, *)'Error: Function gammad.  1st argument < 0 or 2nd argument <= 0'
  RETURN
END IF
IF (x == zero) RETURN

!      Use a normal approximation if P > PLIMIT

IF (p > plimit) THEN
  pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one / (nine * p) - one)
  gamma_prob = alnorm(pn1, .false.)
  RETURN
END IF

!      If X is extremely large compared to P then set gamma_prob = 1

IF (x > xbig) THEN
  gamma_prob = one
  RETURN
END IF

IF (x <= one .OR. x < p) THEN

!      Use Pearson's series expansion.
!      (Note that P is not large enough to force overflow in LNGAMMA)

  arg = p * LOG(x) - x - lngamma(p + one)
  c = one
  gamma_prob = one
  a = p
  DO
    a = a + one
    c = c * x / a
    gamma_prob = gamma_prob + c
    IF (c < tol) EXIT
  END DO
  arg = arg + LOG(gamma_prob)
  gamma_prob = zero
  IF (arg >= elimit) gamma_prob = EXP(arg)

ELSE

!      Use a continued fraction expansion

  arg = p * LOG(x) - x - lngamma(p)
  a = one - p
  b = a + x + one
  c = zero
  pn1 = one
  pn2 = x
  pn3 = x + one
  pn4 = x * b
  gamma_prob = pn3 / pn4
  DO
    a = a + one
    b = b + two
    c = c + one
    an = a * c
    pn5 = b * pn3 - an * pn1
    pn6 = b * pn4 - an * pn2
    IF (ABS(pn6) > zero) THEN
      rn = pn5 / pn6
      IF (ABS(gamma_prob - rn) <= MIN(tol, tol * rn)) EXIT
      gamma_prob = rn
    END IF

    pn1 = pn3
    pn2 = pn4
    pn3 = pn5
    pn4 = pn6
    IF (ABS(pn5) >= oflo) THEN

  !      Re-scale terms in continued fraction if terms are large

      pn1 = pn1 / oflo
      pn2 = pn2 / oflo
      pn3 = pn3 / oflo
      pn4 = pn4 / oflo
    END IF
  END DO
  arg = arg + LOG(gamma_prob)
  gamma_prob = one
  IF (arg >= elimit) gamma_prob = one - EXP(arg)
END IF

RETURN
END FUNCTION gammad



FUNCTION alnorm(x, upper) RESULT(norm_prob)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

REAL (dp), INTENT(IN)  :: x
LOGICAL, INTENT(IN)    :: upper
REAL (dp)              :: norm_prob


! Local variables
REAL (dp) :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
REAL (dp) :: con = 1.28_dp, z, y, ltone = 7.0_dp, utzero = 18.66_dp
REAL (dp) :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
             r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
             a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
             b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
             c1 = -3.8052D-8, c2 = 3.98064794D-4,         &
             c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
             c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
             d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
             d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
             d5 = 30.789933034_dp
LOGICAL   :: up

up = upper
z = x
IF(z >=  zero) GO TO 10
up = .NOT. up
z = -z
10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
norm_prob = zero
GO TO 40
20 y = half*z*z
IF(z > con) GO TO 30

norm_prob = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
GO TO 40
30 norm_prob = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40 IF(.NOT. up) norm_prob = one - norm_prob
RETURN

END FUNCTION alnorm



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
INTEGER   :: j

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


END MODULE Logistic_Regression
