FUNCTION ncbeta (a, b, lambda, x, errmax) RESULT(fn_val)

!       ALGORITHM AS 310 APPL. STATIST. (1997), VOL. 46, NO. 1

!       Computes the cumulative distribution function of a
!       non-central beta random variable

! ELF90-compatible version by Alan Miller
! This version includes all of the other functions called.
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: a, b, lambda, x, errmax
REAL (dp)             :: fn_val

! Local variables
INTEGER               :: xj, m, i, j, iter1, iter2, iterlo, iterhi

!       Local variable XJ gives the number of iterations taken

REAL (dp)             :: fx, gx, temp, ftemp, ebd, errbd, q, r, s, t, s0, &
                         s1, t0, t1, sum, psum, c, beta
REAL (dp), PARAMETER  :: zero = 0.0, half = 0.5, one = 1.0, five = 5.0

! EXTERNAL alngam, betain, betanc, gammad

fn_val = x

!       Check for admissibility of parameters

IF (lambda <= zero .OR. a <= zero .OR. b <= zero) THEN
  WRITE(*, *) 'Illegal value for A, B or LAMBDA in calling FUNCTION NCBETA'
  STOP
END IF
IF (x < zero .OR. x > one) THEN
  WRITE(*, *) 'Illegal value for X in calling FUNCTION NCBETA'
  STOP
END IF
IF (x == zero .OR. x == one) RETURN

c = lambda * half
xj = zero

IF (lambda < 54.0_dp) THEN
  
!       AS 226 as it stands is sufficient in this situation
  
  fn_val = betanc(x, a, b, lambda)
  RETURN
ELSE
  m = INT(c + half)
  iterlo = m - five * SQRT(REAL(m))
  iterhi = m + five * SQRT(REAL(m))
  t = - c + m * LOG(c) - alngam(m + one)
  q = EXP(t)
  r = q
  psum = q
  
  beta = alngam(a+ m) + alngam(b) - alngam(a + m + b)
  s1 = (a + m) * LOG(x) + b * LOG(one - x) - LOG(a + m) - beta
  gx = EXP(s1)
  fx = gx
  temp = betain(x, a + m, b, beta)
  ftemp= temp
  xj = xj + one
  sum = q - temp
  iter1= m
  
!      The first set of iterations starts from M and goes downwards
  
  20 IF (iter1 < iterlo) GO TO 30
  IF (q < errmax) GO TO 30
  q = q - iter1 / c
  xj = xj + one
  gx = (a + iter1) / (x * (a + b + iter1 - one)) * gx
  iter1= iter1- one
  temp =temp + gx
  psum = psum + q
  sum = sum + q * temp
  GO TO 20
  30 t0 = alngam(a + b) - alngam(a + one) - alngam(b)
  s0 = a * LOG(x) + b * LOG(one - x)
  
  DO i=1, iter1
    j = i - one
    s = s + EXP(t0 + s0 + j * LOG(x))
    t1 = LOG(a + b + j) - LOG(a + one + j) + t0
    t0 = t1
  END DO
  
!       Compute the first part of error bound
  
  errbd = (one - gammad(c, DBLE(iter1)))*(temp + s)
  q = r
  temp = ftemp
  gx = fx
  iter2 = m
  DO
    ebd = errbd + (one - psum) * temp
    IF (ebd < errmax .OR. iter2 >= iterhi) EXIT
    iter2 = iter2 + one
    xj = xj + one
    q = q * c / iter2
    psum = psum + q
    temp = temp - gx
    gx = x * (a + b + iter2 - one) / (a + iter2) * gx
  sum = sum + q * temp
  END DO
END IF
fn_val = sum

RETURN

CONTAINS


FUNCTION betanc(x, a, b, lambda) RESULT(fn_val)

!     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
!     Incorporates modification AS R84 from AS vol. 39, pp 311-2, 1990

!     Returns the cumulative probability of X for the non-central beta
!     distribution with parameters A, B and non-centrality LAMBDA

!     Auxiliary routines required: ALNGAM - log-gamma function (ACM
!     291 or AS 245), BETAIN - incomplete-beta function (AS 63), and
!     the incomplete gamma function (gammad, AS 239).

! N.B. Argument IFAULT has been removed

REAL (dp), INTENT(IN) :: x, a, b, lambda
REAL (dp)             :: fn_val

! Local variables
REAL (dp)             :: ax, beta, c, errbd, gx, q, sumq, temp, xj, a0, x0

!     Change ERRMAX and ITRMAX if desired ...

REAL (dp), PARAMETER  :: errmax = 1.0E-6_dp, ualpha = 5.0_dp, zero = 0.0_dp, &
                         half = 0.5_dp, one = 1.0_dp
INTEGER, PARAMETER    :: itrmax = 100

fn_val = x

IF (lambda < zero .OR. a <= zero .OR. b <= zero) THEN
  WRITE(*, *) 'Illegal value for A, B or LAMBDA in calling FUNCTION BETANC'
  STOP
END IF
IF (x < zero .OR. x > one) THEN
  WRITE(*, *) 'Illegal value for X in calling FUNCTION BETANC'
  STOP
END IF
IF (x == zero .OR. x == one) RETURN

c = lambda * half

!     Initialize the series ...

x0 = INT( MAX(c - ualpha*SQRT(c), zero) )
a0 = a + x0
beta = alngam(a0) + alngam(b) - alngam(a0+b)
temp = betain(x, a0, b, beta)
gx = EXP(a0 * LOG(x) + b * LOG(one - x) - beta - LOG(a0))
IF (a0 > a) THEN
  q = EXP(-c + x0*LOG(c)) - alngam(x0 + one)
ELSE
  q = EXP(-c)
END IF
xj = zero
ax = q * temp
sumq = one - q
fn_val = ax

!     Recur over subsequent terms until convergence is achieved...

10 xj = xj + one
temp = temp - gx
gx = x * (a + b + xj - one) * gx / (a + xj)
q = q * c / xj
sumq = sumq - q
ax = temp * q
fn_val = fn_val + ax

!     Check for convergence and act accordingly...

errbd = (temp - gx) * sumq
IF ((INT(xj) < itrmax) .AND. (errbd > errmax)) GO TO 10
IF (errbd > errmax) WRITE(*, *) 'Error bound too large in BETANC'

RETURN
END FUNCTION betanc


FUNCTION betain(x, p, q, beta) RESULT(fn_val)

!     algorithm as 63  appl. statist. (1973), vol.22, no.3

!     computes incomplete beta function ratio for arguments
!     x between zero and one, p and q positive.
!     log of complete beta function, beta, is assumed to be known

! ELF90-compatible version by Alan Miller
! N.B. Argument IFAULT has been removed

! Latest revision - 29 November 1997

REAL (dp), INTENT(IN) :: x, p, q, beta
REAL (dp)             :: fn_val

! Local variables
LOGICAL               :: indx
INTEGER               :: ns
REAL (dp)             :: psq, cx, xx, pp, qq, term, ai, rx, temp

!     define accuracy and initialise

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, acu = 1.0E-14_dp

fn_val = x

!     test for admissibility of arguments

IF(p <= zero .OR. q <= zero) THEN
  WRITE(*, *) 'AS63: Either p or q <= 0'
  RETURN
END IF
IF(x < zero .OR. x > one) THEN
  WRITE(*, *) 'AS63: Argument x outside range (0, 1)'
  RETURN
END IF
IF(x == zero .OR. x == one) RETURN

!     change tail if necessary and determine s

psq = p + q
cx = one - x
IF(p < psq*x) THEN
  xx = cx
  cx = x
  pp = q
  qq = p
  indx = .true.
ELSE
  xx = x
  pp = p
  qq = q
  indx = .false.
END IF
term = one
ai = one
fn_val = one
ns = qq + cx*psq

!     user soper's reduction formulae.

rx = xx/cx
3 temp = qq - ai
IF(ns == 0) rx = xx
4 term = term*temp*rx/(pp+ai)
fn_val = fn_val + term
temp = ABS(term)
IF(temp <= acu .AND. temp <= acu*fn_val) GO TO 5
ai = ai + one
ns = ns - 1
IF(ns >= 0) GO TO 3
temp = psq
psq = psq+one
GO TO 4

!     calculate result

5 fn_val = fn_val*EXP(pp*LOG(xx) + (qq-one)*LOG(cx) - beta)/pp
IF(indx) fn_val = one - fn_val

RETURN
END FUNCTION betain



FUNCTION gammad(x, p) RESULT(fn_val)

!      ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!      Computation of the Incomplete Gamma Integral

!      Auxiliary functions required: ALOGAM = logarithm of the gamma
!      function, and ALNORM = algorithm AS66

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: x, p
REAL (dp)             :: fn_val

! Local variables
REAL (dp)             :: pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an
REAL (dp), PARAMETER  :: zero = 0.d0, one = 1.d0, two = 2.d0, &
                         oflo = 1.d+37, three = 3.d0, nine = 9.d0, &
                         tol = 1.d-14, xbig = 1.d+8, plimit = 1000.d0, &
                         elimit = -88.d0
! EXTERNAL alogam, alnorm

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
!      (Note that P is not large enough to force overflow in ALNGAM).
  
  arg = p * LOG(x) - x - alngam(p + one)
  c = one
  fn_val = one
  a = p
  40   a = a + one
  c = c * x / a
  fn_val = fn_val + c
  IF (c > tol) GO TO 40
  arg = arg + LOG(fn_val)
  fn_val = zero
  IF (arg >= elimit) fn_val = EXP(arg)
  
ELSE
  
!      Use a continued fraction expansion
  
  arg = p * LOG(x) - x - alngam(p)
  a = one - p
  b = a + x + one
  c = zero
  pn1 = one
  pn2 = x
  pn3 = x + one
  pn4 = x * b
  fn_val = pn3 / pn4
  60   a = a + one
  b = b + two
  c = c + one
  an = a * c
  pn5 = b * pn3 - an * pn1
  pn6 = b * pn4 - an * pn2
  IF (ABS(pn6) > zero) THEN
    rn = pn5 / pn6
    IF (ABS(fn_val - rn) <= MIN(tol, tol * rn)) GO TO 80
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
  GO TO 60
  80   arg = arg + LOG(fn_val)
  fn_val = one
  IF (arg >= elimit) fn_val = one - EXP(arg)
END IF

RETURN
END FUNCTION gammad


FUNCTION alngam(xvalue) RESULT(fn_val)

!     ALGORITHM AS245  APPL. STATIST. (1989) VOL. 38, NO. 2

!     Calculation of the logarithm of the gamma function

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: xvalue
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x, x1, x2, y

!     Coefficients of rational functions

REAL (dp), PARAMETER :: r1(9) = (/ -2.66685511495_dp, -24.4387534237_dp,  &
                                   -21.9698958928_dp,  11.1667541262_dp,  &
                                    3.13060547623_dp,  0.607771387771_dp, &
                                    11.9400905721_dp,  31.4690115749_dp,  &
                                    15.2346874070_dp /)
REAL (dp), PARAMETER :: r2(9) = (/ -78.3359299449_dp, -142.046296688_dp,  &
                                    137.519416416_dp,  78.6994924154_dp,  &
                                    4.16438922228_dp,  47.0668766060_dp,  &
                                    313.399215894_dp,  263.505074721_dp,  &
                                    43.3400022514_dp /)
REAL (dp), PARAMETER :: r3(9) = (/ -2.12159572323E5_dp,  2.30661510616E5_dp,  &
                                    2.74647644705E4_dp, -4.02621119975E4_dp,  &
                                   -2.29660729780E3_dp, -1.16328495004E5_dp,  &
                                   -1.46025937511E5_dp, -2.42357409629E4_dp,  &
                                   -5.70691009324E2_dp /)
REAL (dp), PARAMETER :: r4(5) = (/ 0.279195317918525_dp, 0.4917317610505968_dp, &
                                   0.0692910599291889_dp, 3.350343815022304_dp, &
                                   6.012459259764103_dp /)

!     Fixed constants

REAL (dp), PARAMETER :: alr2pi = 0.918938533204673_dp, four = 4._dp,  &
                        half = 0.5_dp, one = 1._dp, onep5 = 1.5_dp,   &
                        twelve = 12._dp, zero = 0._dp

!     Machine-dependant constants.
!     A table of values is given at the top of page 399 of the paper.
!     These values are for the IEEE double-precision format for which
!     B = 2, t = 53 and U = 1023 in the notation of the paper.

REAL (dp), PARAMETER :: xlge = 5.10E6_dp, xlgst = HUGE(1.0_dp)

x = xvalue
fn_val = zero

!     Test for valid function argument

IF (x >= xlgst) THEN
  WRITE(*, *) 'AS 245: Argument x too large'
  RETURN
END IF
IF (x <= zero) THEN
  WRITE(*, *) 'AS 245: Argument x <= 0'
  RETURN
END IF

!     Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined

IF (x < onep5) THEN
  IF (x < half) THEN
    fn_val = -LOG(x)
    y = x + one
    
!     Test whether X < machine epsilon
    
    IF (y == one) RETURN
  ELSE
    fn_val = zero
    y = x
    x = (x - half) - half
  END IF
  fn_val = fn_val + x * ((((r1(5)*y + r1(4))*y + r1(3))*y + r1(2))*y + r1(1)) / &
                    ((((y + r1(9))*y + r1(8))*y+ r1(7))*y + r1(6))
  RETURN
END IF

!     Calculation for 1.5 <= X < 4.0

IF (x < four) THEN
  y = (x - one) - one
  fn_val = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x + r2(1)) /  &
               ((((x + r2(9))*x + r2(8))*x + r2(7))*x+ r2(6))
  RETURN
END IF

!     Calculation for 4.0 <= X < 12.0

IF (x < twelve) THEN
  fn_val = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) /  &
           ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
  RETURN
END IF

!     Calculation for X >= 12.0

y = LOG(x)
fn_val = x * (y - one) - half * y + alr2pi
IF (x > xlge) RETURN
x1 = one / x
x2 = x1 * x1
fn_val = fn_val + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) /  &
         ((x2 + r4(5))*x2 + r4(4))
RETURN
END FUNCTION alngam



FUNCTION alnorm(x, upper) RESULT(fn_val)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

REAL (dp), INTENT(IN) :: x
LOGICAL, INTENT(IN)   :: upper
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp, &
                        con = 1.28_dp
REAL (dp) :: z, y
LOGICAL   :: up

!*** machine dependent constants
REAL (dp), PARAMETER :: ltone = 7.0_dp, utzero = 18.66_dp

REAL (dp), PARAMETER :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
                        r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
                        a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
                        b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
                        c1 = -3.8052E-8_dp, c2 = 3.98064794E-4_dp,     &
                        c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
                        c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
                        d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
                        d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
                        d5 = 30.789933034_dp

up = upper
z = x
IF(z >=  zero) GO TO 10
up = .NOT. up
z = -z
10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
fn_val = zero
GO TO 40
20 y = half*z*z
IF(z > con) GO TO 30

fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
GO TO 40
30 fn_val = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40 IF(.NOT. up) fn_val = one - fn_val

RETURN
END FUNCTION alnorm

END FUNCTION ncbeta
