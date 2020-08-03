MODULE Owens_T
! Calculates Owen's T-function which is useful in calculating bivariate
! normal probabilities.   Also included is function biv_norm for calculating
! bivariate normal probabilities.

! From the paper:
! Patefield, Mike `Fast and accurate calculation of Owen's T-function',
! J. of Statistical Software, vol. 5, 2000
! http://www.stat.ucla.edu/journals/jss/v05

! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-21  Time: 21:16:36

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PUBLIC   :: t, alnorm, znorm1, znorm2
PRIVATE  :: tf

CONTAINS


FUNCTION t(h, a) RESULT(fn_val)
 
!        COMPUTES OWEN'S T-FUNCTION OF H AND A
!        WITH H,A ANY REAL (dp) NUMBERS

REAL (dp), INTENT(IN)  :: h
REAL (dp), INTENT(IN)  :: a
REAL (dp)              :: fn_val

REAL (dp)  :: absa, absh, ah, normh, normah
REAL (dp), PARAMETER  :: zero = 0.0_dp, quart = 0.25_dp, half = 0.5_dp,  &
                         cut = 0.67_dp, one = 1.0_dp

absh=ABS(h)
absa=ABS(a)
ah = absa * absh
IF (absa <= one) THEN
  fn_val = tf(absh, absa, ah)
ELSE
  IF (absh <= cut) THEN
    fn_val = quart - znorm1(absh) * znorm1(ah) - tf(ah, one / absa, absh)
  ELSE
    normh = znorm2(absh)
    normah = znorm2(ah)
    fn_val = half * ( normh + normah ) - normh * normah -  &
             tf(ah, one / absa, absh)
  END IF
END IF
IF (a < zero) fn_val = - fn_val

RETURN
END FUNCTION t



FUNCTION tf(h, a, ah) RESULT(fn_val)

!        COMPUTES OWEN'S T-FUNCTION OF H AND A
!        H >= 0 AND 0 <= A <= 1 ; INPUT AH MUST EQUAL A * H

REAL (dp), INTENT(IN)  :: h
REAL (dp), INTENT(IN)  :: a
REAL (dp), INTENT(IN)  :: ah
REAL (dp)              :: fn_val

REAL (dp)  :: z, zi, x, vi, hs, dhs, as, aj, yi, dj, gj, r, ai, y, normh
INTEGER    :: ihint, iaint, m, j, jj, i, maxii, ii, icode
REAL (dp), PARAMETER  :: rtwopi = 0.15915494309189533577_dp,   &
                         rrtpi = 0.39894228040143267794_dp,  &
                         zero = 0.0_dp, mhalf = -0.5_dp, half = 0.5_dp,  &
                         one = 1.0_dp
REAL (dp), PARAMETER  :: c2(21) = (/ 0.99999999999999987510D+00,  &
       -0.99999999999988796462D+00,  0.99999999998290743652D+00,  &
       -0.99999999896282500134D+00,  0.99999996660459362918D+00,  &
       -0.99999933986272476760D+00,  0.99999125611136965852D+00,  &
       -0.99991777624463387686D+00,  0.99942835555870132569D+00,  &
       -0.99697311720723000295D+00,  0.98751448037275303682D+00,  &
       -0.95915857980572882813D+00,  0.89246305511006708555D+00,  &
       -0.76893425990463999675D+00,  0.58893528468484693250D+00,  &
       -0.38380345160440256652D+00,  0.20317601701045299653D+00,  &
       -0.82813631607004984866D-01,  0.24167984735759576523D-01,  &
       -0.44676566663971825242D-02,  0.39141169402373836468D-03 /)
REAL (dp), PARAMETER  :: pts(13) = (/ 0.35082039676451715489D-02,  &
         0.31279042338030753740D-01,  0.85266826283219451090D-01,  &
         0.16245071730812277011D+00,  0.25851196049125434828D+00,  &
         0.36807553840697533536D+00,  0.48501092905604697475D+00,  &
         0.60277514152618576821D+00,  0.71477884217753226516D+00,  &
         0.81475510988760098605D+00,  0.89711029755948965867D+00,  &
         0.95723808085944261843D+00,  0.99178832974629703586D+00 /)
REAL (dp), PARAMETER  :: wts(13) = (/ 0.18831438115323502887D-01,  &
         0.18567086243977649478D-01,  0.18042093461223385584D-01,  &
         0.17263829606398753364D-01,  0.16243219975989856730D-01,  &
         0.14994592034116704829D-01,  0.13535474469662088392D-01,  &
         0.11886351605820165233D-01,  0.10070377242777431897D-01,  &
         0.81130545742299586629D-02,  0.60419009528470238773D-02,  &
         0.38862217010742057883D-02,  0.16793031084546090448D-02 /)
INTEGER, PARAMETER    :: meth(18) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 2,  &
                                       2, 2, 3, 4, 4, 4, 4, 5, 6 /)
INTEGER, PARAMETER    :: ord(18) = (/  2, 3, 4, 5, 7,10,12,18,10,  &
                                      20,30,20, 4, 7, 8,20,13, 0 /)
REAL (dp), PARAMETER  :: hrange(14) = (/ 0.02_dp, 0.06_dp, 0.09_dp, 0.125_dp, &
                         0.26_dp, 0.4_dp,  0.6_dp, 1.6_dp, 1.7_dp, 2.33_dp,  &
                         2.4_dp,  3.36_dp, 3.4_dp, 4.8_dp /)
REAL (dp), PARAMETER  :: arange(7) = (/ 0.025_dp, 0.09_dp, 0.15_dp, 0.36_dp,  &
                                        0.5_dp,   0.9_dp,  0.99999_dp /)
INTEGER, PARAMETER  :: select(15,8) = RESHAPE(  &
 (/ 1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9,  &
    1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9,  &
    2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10,  &
    2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10,  &
    2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11,  &
    2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12,  &
    2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12,  &
    2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12 /), (/ 15, 8 /) )

!  DETERMINE APPROPRIATE METHOD FROM T1...T6

DO  ihint=1,14
  IF (h <= hrange(ihint)) GO TO 20
END DO
ihint=15
20 DO  iaint=1,7
  IF (a <= arange(iaint)) GO TO 40
END DO
iaint=8
40 icode = select(ihint,iaint)
m = ord(icode)
SELECT CASE ( meth(icode) )
  CASE (    1)
    GO TO 100
  CASE (    2)
    GO TO 200
  CASE (    3)
    GO TO 300
  CASE (    4)
    GO TO 400
  CASE (    5)
    GO TO 500
  CASE (    6)
    GO TO 600
END SELECT

!  T1(H, A, M) ; M = 2, 3, 4, 5, 7, 10, 12 OR 18
!  JJ = 2J - 1 ; GJ = EXP(-H*H/2) * (-H*H/2)**J / J!
!  AJ = A**(2J-1) / (2*PI)

100 hs = mhalf * h * h
dhs = EXP(hs)
as = a * a
j = 1
jj = 1
aj = rtwopi * a
fn_val = rtwopi * ATAN(a)
dj = dhs - one
gj = hs * dhs
110 fn_val = fn_val + dj * aj / jj
IF (j >= m) RETURN
j = j + 1
jj = jj + 2
aj = aj * as
dj = gj - dj
gj = gj * hs / j
GO TO 110

!  T2(H, A, M) ; M = 10, 20 OR 30
!  Z = (-1)**(I-1) * ZI ; II = 2I - 1
!  VI = (-1)**(I-1) * A**(2I-1) * EXP[-(A*H)**2/2] / SQRT(2*PI)

200 maxii = m + m + 1
ii = 1
fn_val = zero
hs = h * h
as = -a * a
vi = rrtpi * a * EXP(mhalf * ah * ah)
z = znorm1(ah) / h
y = one / hs
210 fn_val = fn_val + z
IF (ii >= maxii) GO TO 220
z = y * (vi - ii * z)
vi = as * vi
ii = ii + 2
GO TO 210
220 fn_val = fn_val * rrtpi * EXP (mhalf * hs)
RETURN

!  T3(H, A, M) ; M = 20
!  II = 2I - 1
!  VI = A**(2I-1) * EXP[-(A*H)**2/2] / SQRT(2*PI)

300 i = 1
ii = 1
fn_val = zero
hs = h * h
as = a * a
vi = rrtpi * a * EXP(mhalf * ah * ah)
zi = znorm1(ah) / h
y = one / hs
310 fn_val = fn_val + zi * c2(i)
IF (i > m) GO TO 320
zi = y  * (ii * zi - vi)
vi = as * vi
i = i + 1
ii = ii + 2
GO TO 310
320 fn_val = fn_val * rrtpi * EXP(mhalf * hs)
RETURN

!  T4(H, A, M) ; M = 4, 7, 8 OR 20;  II = 2I + 1
!  AI = A * EXP[-H*H*(1+A*A)/2] * (-A*A)**I / (2*PI)

400 maxii = m + m + 1
ii = 1
hs = h * h
as = -a * a
fn_val = zero
ai = rtwopi * a * EXP(mhalf * hs * (one - as))
yi = one
410 fn_val = fn_val + ai * yi
IF (ii >= maxii) RETURN
ii = ii + 2
yi = (one - hs * yi) / ii
ai = ai * as
GO TO 410

!  T5(H, A, M) ; M = 13
!  2M - POINT GAUSSIAN QUADRATURE

500 fn_val = zero
as = a * a
hs = mhalf * h * h
DO  i = 1 , m
  r = one + as * pts(i)
  fn_val = fn_val + wts(i) * EXP(hs * r) / r
END DO
fn_val = a * fn_val
RETURN

!  T6(H, A);  APPROXIMATION FOR A NEAR 1, (A <= 1)

600 normh = znorm2(h)
fn_val = half * normh * (one - normh)
y = one - a
x = ATAN(y / (one + a))
fn_val = fn_val - rtwopi * x * EXP(mhalf * y * h * h / x)

RETURN
END FUNCTION tf



FUNCTION alnorm(x, upper) RESULT(fn_val)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

REAL (dp), INTENT(IN)  :: x
LOGICAL, INTENT(IN)    :: upper
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp, &
                         con = 1.28_dp
REAL (dp)  :: z, y
LOGICAL    :: up

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



FUNCTION znorm1(x) RESULT(fn_val)
! Area under the standard normal curve between 0 and x

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = alnorm(x, .FALSE.) - 0.5_dp

RETURN
END FUNCTION znorm1



FUNCTION znorm2(x) RESULT(fn_val)
! Area under the standard normal curve to the right of x

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = alnorm(x, .TRUE.)

RETURN
END FUNCTION znorm2



FUNCTION biv_norm(h, xmean, xsd, k, ymean, ysd, rho) RESULT(fn_val)
! Calculate bivariate normal probabilities,
! i.e. Prob(X > h, Y > k) where X, Y have a bivariate normal distribution.

! Reference:
! Baughman, A.L. (1988) `A Fortran function for the bivariate normal integral',
! Computer Methods and Programs in Biomedicine, Vol.27, 169-174.

! Evaluates L(h,k,r) = Prob(X > h, Y > k)
! where X and Y are standardized normal variables and R is their correlation.
! Uses Owen's relationship:
! L(h,k,r) = 1.0 - (PH + PK)/2 - T(h, ah) - T(k, ak)
! where PH = tail area of the standardized normal curve from minus infinity
!            up to h,
!       PK = tail area of the standardized normal curve from minus infinity
!            up to k,
!   T(h,a) = integral from 0 to a of
!            (1/(2.pi)).EXP(-0.5.h^2.(1 + x^2)) / (1 + x^2)
!       ah = (h - r.k) / (h.sqrt(1 - r^2))
!       ak = (k - r.h) / (k.sqrt(1 - r^2))

! Uses the TFUNC of Mike Patefield

! Arguments:

! h      lower limit for X
! xmean  mean of X
! xsd    standard deviation of X

! k      lower limit for Y
! ymean  mean of Y
! ysd    standard deviation of Y

! rho    correlation coefficient

REAL (dp), INTENT(IN)  :: h, xmean, xsd, k, ymean, ysd, rho
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: rr, horig, korig, ah, ak, hh, kk, thah, tkak
REAL (dp), PARAMETER  :: zero = 0.0_dp, quart = 0.25_dp, half = 0.5_dp,  &
                         one = 1.0_dp, eight = 8.0_dp

! Check input arguments

IF (ABS(rho) > one) THEN
  WRITE(*, *) '** Error in FUNCTION biv_norm **'
  WRITE(*, '(a, g13.5, a)') ' Argument RHO = ', rho, ' It must be between -1 and 1'
END IF
IF (xsd <= zero) THEN
  WRITE(*, *) '** Error in FUNCTION biv_norm **'
  WRITE(*, '(a, g13.5, a)') ' Argument XSD = ', xsd, ' It must be between > 0'
END IF
IF (ysd <= zero) THEN
  WRITE(*, *) '** Error in FUNCTION biv_norm **'
  WRITE(*, '(a, g13.5, a)') ' Argument YSD = ', ysd, ' It must be between > 0'
END IF

! Standardize h and k

hh = (h - xmean) / xsd
kk = (k - ymean) / ysd

! Calculate probability only for positive hh or kk;
! adjust to the right quadrant later.

horig = hh
korig = kk
hh = ABS(hh)
kk = ABS(kk)

! The relations between L(h,k,r) and PX, the tail area up to X for the
! standard normal curve are equations (22.4) to (22.6) in Johnson, N.L.
! and Kotz, S. `Distributions in Statistics: Continuous Multivariate
! Distributions, Ch. 36, Wiley: New York (1972).

IF (rho == zero) THEN
  fn_val = alnorm(horig, .TRUE.) * alnorm(korig, .TRUE.)
  GO TO 100

ELSE IF (rho == -one) THEN
  fn_val = zero
  IF (horig + korig < zero) fn_val = one - alnorm(horig, .FALSE.) -  &
                                           alnorm(korig, .FALSE.)
  GO TO 100

ELSE IF (rho == one) THEN
  fn_val = alnorm(MAX(horig, korig), .TRUE.)
  GO TO 100

ELSE

! For negative h and positive k, or vice-vera, reverse the sign of rho.
! Reference equation (22.2) of J & K.

  IF ((horig < zero .AND. korig >= zero) .OR.  &
      (horig >= zero .AND. korig < zero) ) THEN
    rr = -rho
    IF(hh == zero .AND. kk == zero) THEN
      fn_val = ASIN(rr) / (eight * ATAN(one)) + quart
      GO TO 100
    ELSE IF (hh == zero .AND. kk /= zero) THEN
      thah = quart
      ak = (hh - kk*rr) / (kk*SQRT(one - rr**2))
      tkak = t(k,ak)
    ELSE IF (hh /= zero .AND. kk == zero) THEN
      ak = (kk - hh*rr) / (hh*SQRT(one - rr**2))
      thah = t(h,ah)
      tkak = quart
    ELSE
      ah = (kk - hh*rr) / (hh*SQRT(one - rr**2))
      thah = t(h,ah)
      ak = (hh - kk*rr) / (kk*SQRT(one - rr**2))
      tkak = t(k,ak)
    END IF

  ELSE
    IF(hh == zero .AND. kk == zero) THEN
      fn_val = ASIN(rho) / (eight * ATAN(one)) + quart
      GO TO 100
    ELSE IF (hh == zero .AND. kk /= zero) THEN
      thah = quart
      ak = (hh - kk*rho) / (kk*SQRT(one - rho**2))
      tkak = t(k,ak)
    ELSE IF (hh /= zero .AND. kk == zero) THEN
      ak = (kk - hh*rho) / (hh*SQRT(one - rho**2))
      thah = t(h,ah)
      tkak = quart
    ELSE
      ah = (kk - hh*rho) / (hh*SQRT(one - rho**2))
      thah = t(h,ah)
      ak = (hh - kk*rho) / (kk*SQRT(one - rho**2))
      tkak = t(k,ak)
    END IF
  END IF
END IF

fn_val = one - half * (alnorm(hh, .FALSE.) + alnorm(kk, .FALSE.))  &
             - thah - tkak

IF (horig < zero .AND. korig >= zero) THEN
  fn_val = - fn_val + alnorm(kk, .TRUE.)
END IF

IF (horig >= zero .AND. korig < zero) THEN
  fn_val = - fn_val + alnorm(hh, .TRUE.)
END IF

IF (horig < zero .AND. korig < zero) THEN
  fn_val = fn_val + alnorm(hh, .FALSE.) + alnorm(kk, .FALSE.) - one
END IF

100 fn_val = MAX(zero, fn_val)
fn_val = MIN(one, fn_val)

RETURN
END FUNCTION biv_norm

END MODULE Owens_T
