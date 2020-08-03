PROGRAM test_chi_squared
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER             :: ndf
REAL (dp)           :: chi2

INTERFACE
  FUNCTION chi_squared(ndf, chi2) RESULT(prob)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    INTEGER, INTENT(IN)    :: ndf
    REAL (dp), INTENT(IN)  :: chi2
    REAL (dp)              :: prob
  END FUNCTION chi_squared
END INTERFACE

DO
  WRITE(*, *) 'Enter no. of degrees of freedom: '
  READ(*, *) ndf
  WRITE(*, *) 'Enter chi-squared value: '
  READ(*, *) chi2
  WRITE(*, *) 'Chi-squared probability = ', chi_squared(ndf, chi2)
END DO
STOP

END PROGRAM test_chi_squared



FUNCTION chi_squared(ndf, chi2) RESULT(prob)
! Calculate the chi-squared distribution function
! ndf  = number of degrees of freedom
! chi2 = chi-squared value
! prob = probability of a chi-squared value <= chi2 (i.e. the left-hand
!        tail area)

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER, INTENT(IN)    :: ndf
REAL (dp), INTENT(IN)  :: chi2
REAL (dp)              :: prob

INTERFACE
  FUNCTION gammad(x, p) RESULT(gamma_prob)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)  :: x, p
    REAL (dp)              :: gamma_prob
  END FUNCTION gammad
END INTERFACE

! Local variables
REAL (dp)  :: half = 0.5d0, x, p

x = half * chi2
p = half * REAL(ndf)
prob = gammad(x, p)
RETURN

END FUNCTION chi_squared



FUNCTION gammad(x, p) RESULT(fn_val)

!   ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!   Computation of the Incomplete Gamma Integral

!   Auxiliary functions required: lngamma = logarithm of the gamma
!   function, and ALNORM = algorithm AS66

! ELF90-compatible version by Alan Miller
! amiller@bigpond.net.au
! Latest revision - 28 October 2000

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)  :: x, p
REAL (dp)              :: fn_val

INTERFACE
  FUNCTION alnorm( x, upper ) RESULT( fn_val )
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
    REAL(dp), INTENT(IN)  :: x
    LOGICAL, INTENT(IN)   :: upper
    REAL(DP)              :: fn_val
  END FUNCTION alnorm

  FUNCTION lngamma(z) RESULT(lanczos)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 60)
    REAL(dp), INTENT(IN)  :: z
    REAL(dp)              :: lanczos
  END FUNCTION lngamma
END INTERFACE

! Local variables
REAL (dp)             :: pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, &
                         oflo = 1.d+37, three = 3.0_dp, nine = 9.0_dp, &
                         tol = 1.d-14, xbig = 1.d+8, plimit = 1000.0_dp, &
                         elimit = -88.0_dp
! EXTERNAL lngamma, alnorm

fn_val = zero

!      Check that we have valid values for X and P

IF (p <= zero .OR. x < zero) THEN
  WRITE(*, *) 'AS239: Either p <= 0 or x < 0'
  RETURN
END IF
IF (x == zero) RETURN

!      Use a normal approximation if P > PLIMIT

IF (p > plimit) THEN
  pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one /(nine * p) - one)
  fn_val = alnorm(pn1, .false.)
  RETURN
END IF

!      If X is extremely large compared to P then set fn_val = 1

IF (x > xbig) THEN
  fn_val = one
  RETURN
END IF

IF (x <= one .OR. x < p) THEN
  
!      Use Pearson's series expansion.
!      (Note that P is not large enough to force overflow in lngamma).
!      No need to test IFAULT on exit since P > 0.
  
  arg = p * LOG(x) - x - lngamma(p + one)
  c = one
  fn_val = one
  a = p
  DO
    a = a + one
    c = c * x / a
    fn_val = fn_val + c
    IF (c <= tol) EXIT
  END DO
  arg = arg + LOG(fn_val)
  fn_val = zero
  IF (arg >= elimit) fn_val = EXP(arg)
  
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
  fn_val = pn3 / pn4
  DO
    a = a + one
    b = b + two
    c = c + one
    an = a * c
    pn5 = b * pn3 - an * pn1
    pn6 = b * pn4 - an * pn2
    IF (ABS(pn6) > zero) THEN
      rn = pn5 / pn6
      IF (ABS(fn_val - rn) <= MIN(tol, tol * rn)) EXIT
      fn_val = rn
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

  arg = arg + LOG(fn_val)
  fn_val = one
  IF (arg >= elimit) fn_val = one - EXP(arg)
END IF


RETURN
END FUNCTION gammad



!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 2001

FUNCTION alnorm( x, upper ) RESULT( fn_val )
   IMPLICIT NONE
   INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
   REAL(dp), INTENT(IN)  :: x
   LOGICAL, INTENT(IN)   :: upper
   REAL(dp)              :: fn_val

   !  Local variables
   REAL(dp), PARAMETER   ::  zero=0.0_dp, one=1.0_dp, half=0.5_dp, con=1.28_dp
   REAL(dp)              ::  z, y
   LOGICAL               ::  up

   !  Machine dependent constants
   REAL(dp), PARAMETER  :: ltone = 7.0_dp, utzero = 18.66_dp
   REAL(dp), PARAMETER  :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
                           r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
                           a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
                           b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
                           c1 = -3.8052D-8, c2 = 3.98064794D-4,           &
                           c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
                           c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
                           d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
                           d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
                           d5 = 30.789933034_dp

   up = upper
   z = x
   IF( z < zero ) THEN
      up = .NOT. up
      z = -z
   END IF
   IF( z <= ltone  .OR.  (up  .AND.  z <= utzero) ) THEN
      y = half*z*z
      IF( z > con ) THEN
         fn_val = r*EXP( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
      ELSE
         fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
      END IF
   ELSE
      fn_val = zero
   END IF

   IF( .NOT. up ) fn_val = one - fn_val
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

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 60)
REAL(dp), INTENT(IN)  :: z
REAL(dp)              :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                       -1259.139216722289D0, 771.3234287757674D0, &
                       -176.6150291498386D0, 12.50734324009056D0, &
                       -0.1385710331296526D0, 0.9934937113930748D-05, &
                        0.1659470187408462D-06 /), zero = 0.D0,   &
             one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
             half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER   :: j

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


