MODULE Complex_Bessel
! SUBROUTINE CBSSLJ evaluates Bessel functions of complex ORDER
! SUBROUTINE CBSSLJ(x, nu, w)
! x = argument, nu = order, w = function value

! Latest revision - 20 July 2001
! amiller @ bigpond.net.au

! WARNING (from the NSWC Math. Library Manual)
! Precision:
! The real and imaginary parts are normally accurate to 11-12 significant
! digits when the function value is not near zero.   The exception is when
! Real(z) is small, Imag(nu) is small, and 14 < |z| < 17.5 + 0.5|nu|^2
! In this region, the real part of the result is still accurate, but all
! accuracy may be lost in the imaginary part.

USE constants_NSWC
IMPLICIT NONE

CONTAINS


FUNCTION cpabs(x, y) RESULT(fn_val)
!     --------------------------------------
!     EVALUATION OF SQRT(X*X + Y*Y)
!     --------------------------------------
REAL (dp), INTENT(IN)  :: x, y
REAL (dp)              :: fn_val

! Local variable
REAL (dp)  :: a

IF (ABS(x) > ABS(y)) THEN
  a = y / x
  fn_val = ABS(x) * SQRT(1.0_dp + a*a)
  RETURN
END IF
IF (y /= 0.0_dp) THEN
  a = x / y
  fn_val = ABS(y) * SQRT(1.0_dp + a*a)
  RETURN
END IF
fn_val = 0.0_dp
RETURN
END FUNCTION cpabs



SUBROUTINE dcrec(x, y, u, v)
!-----------------------------------------------------------------------
!             COMPLEX RECIPROCAL U + I*V = 1/(X + I*Y)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)   :: x, y
REAL (dp), INTENT(OUT)  :: u, v

! Local variables
REAL (dp)  :: t, d

IF (ABS(x) <= ABS(y)) THEN
  t = x / y
  d = y + t * x
  u = t / d
  v = -1.0_dp / d
  RETURN
END IF
t = y / x
d = x + t * y
u = 1.0_dp / d
v = -t / d
RETURN
END SUBROUTINE dcrec



FUNCTION cdiv(a,b) RESULT(fn_val)
!-----------------------------------------------------------------------
!              COMPLEX DIVISION A/B WHERE B IS NONZERO
!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)  :: a, b
COMPLEX (dp)              :: fn_val


REAL (dp)  :: ai, ar, bi, br, d, t
REAL (dp)  :: u, v

ar = REAL(a, KIND=dp)
ai = AIMAG(a)
br = REAL(b, KIND=dp)
bi = AIMAG(b)

IF (ABS(br) >= ABS(bi)) THEN
  t = bi / br
  d = br + t * bi
  u = (ar+ai*t) / d
  v = (ai-ar*t) / d
  fn_val = CMPLX(u,v, KIND=dp)
  RETURN
END IF
t = br / bi
d = bi + t * br
u = (ar*t+ai) / d
v = (ai*t-ar) / d
fn_val = CMPLX(u,v, KIND=dp)
RETURN
END FUNCTION cdiv



FUNCTION dsin1(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!                REAL (dp) EVALUATION OF SIN(PI*X)

!                             --------------

!     THE EXPANSION FOR SIN(PI*A) (ABS(A) <= PI/4) USING A1,...,A13
!     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!     THE EXPANSION FOR COS(PI*A) (ABS(A) <= PI/4) USING B1,...,B13
!     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)            :: a, t, w
REAL (dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197D+00
REAL (dp), PARAMETER :: a1 = -.1028083791780141522795259479153765743002D+00,   &
      a2  = .3170868848763100170457042079710451905600D-02,   &
      a3  = -.4657026956105571623449026167864697920000D-04,  &
      a4  = .3989844942879455643410226655783424000000D-06,   &
      a5  = -.2237397227721999776371894030796800000000D-08,  &
      a6  = .8847045483056962709715066675200000000000D-11,   &
      a7  = -.2598715447506450292885585920000000000000D-13,  &
      a8  = .5893449774331011070033920000000000000000D-16 ,  &
      a9  = -.1062975472045522550784000000000000000000D-18,   &
      a10 = .1561182648301780992000000000000000000000D-21,    &
      a11 = -.1903193516670976000000000000000000000000D-24,   &
      a12 = .1956617650176000000000000000000000000000D-27,    &
      a13 = -.1711276032000000000000000000000000000000D-30
REAL (dp), PARAMETER :: b1 = -.3084251375340424568385778437461297229882D+00, &
      b2  = .1585434424381550085228521039855226435920D-01,   &
      b3  = -.3259918869273900136414318317506279360000D-03,  &
      b4  = .3590860448591510079069203991239232000000D-05,   &
      b5  = -.2461136950494199754009084061808640000000D-07,  &
      b6  = .1150115912797405152263195572224000000000D-09,   &
      b7  = -.3898073171259675439899172864000000000000D-12,  &
      b8  = .1001886461636271969091584000000000000000D-14,   &
      b9  = -.2019653396886572027084800000000000000000D-17,  &
      b10 = .3278483561466560512000000000000000000000D-20,   &
      b11 = -.4377345082051788800000000000000000000000D-23,  &
      b12 = .4891532381388800000000000000000000000000D-26,   &
      b13 = -.4617089843200000000000000000000000000000D-29
INTEGER  :: max, n
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
t = MAX
IF (a >= t) THEN
  fn_val = 0.0_dp
  RETURN
END IF

n = a
t = n
a = a - t
IF (a <= 0.75_dp) THEN
  IF (a < 0.25_dp) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25_dp + (0.25_dp-a)
  t = 16._dp * a * a
  fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t  &
           + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t +  &
           0.5_dp) + 0.5_dp
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25_dp + (0.75_dp-a)
10 t = 16._dp * a * a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t + a8)*t + a7)*t  &
    + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5_dp) + 0.5_dp
fn_val = pi * a * w

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dsin1



FUNCTION dcos1 (x) RESULT(fn_val)
 
!-----------------------------------------------------------------------

!                REAL (dp) EVALUATION OF COS(PI*X)

!                             --------------

!     THE EXPANSION FOR SIN(PI*A) (ABS(A) .LE. PI/4) USING A1,...,A13
!     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!     THE EXPANSION FOR COS(PI*A) (ABS(A) .LE. PI/4) USING B1,...,B13
!     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp)  :: a, t, w
INTEGER    :: MAX, n
!------------------------
REAL (dp), PARAMETER  :: pi = 3.141592653589793238462643383279502884197_dp
!------------------------
REAL (dp), PARAMETER  :: &
    a1  = -.1028083791780141522795259479153765743002D+00,  &
    a2  =  .3170868848763100170457042079710451905600D-02,  &
    a3  = -.4657026956105571623449026167864697920000D-04,  &
    a4  =  .3989844942879455643410226655783424000000D-06,  &
    a5  = -.2237397227721999776371894030796800000000D-08,  &
    a6  =  .8847045483056962709715066675200000000000D-11,  &
    a7  = -.2598715447506450292885585920000000000000D-13,  &
    a8  =  .5893449774331011070033920000000000000000D-16,  &
    a9  = -.1062975472045522550784000000000000000000D-18,  &
    a10 =  .1561182648301780992000000000000000000000D-21,  &
    a11 = -.1903193516670976000000000000000000000000D-24,  &
    a12 =  .1956617650176000000000000000000000000000D-27,  &
    a13 = -.1711276032000000000000000000000000000000D-30
!------------------------
REAL (dp), PARAMETER  :: &
    b1  = -.3084251375340424568385778437461297229882D+00,  &
    b2  =  .1585434424381550085228521039855226435920D-01,  &
    b3  = -.3259918869273900136414318317506279360000D-03,  &
    b4  =  .3590860448591510079069203991239232000000D-05,  &
    b5  = -.2461136950494199754009084061808640000000D-07,  &
    b6  =  .1150115912797405152263195572224000000000D-09,  &
    b7  = -.3898073171259675439899172864000000000000D-12,  &
    b8  =  .1001886461636271969091584000000000000000D-14,  &
    b9  = -.2019653396886572027084800000000000000000D-17,  &
    b10 =  .3278483561466560512000000000000000000000D-20,  &
    b11 = -.4377345082051788800000000000000000000000D-23,  &
    b12 =  .4891532381388800000000000000000000000000D-26,  &
    b13 = -.4617089843200000000000000000000000000000D-29
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

MAX = HUGE(0)

!------------------------
a = ABS(x)
t = MAX
IF (a < t) GO TO 10
fn_val = 1.d0
RETURN

10 n = a
t = n
a = a - t
IF (a > 0.75D0) GO TO 20
IF (a < 0.25D0) GO TO 21

!                    0.25 .LE. A .LE. 0.75

a = 0.25D0 + (0.25D0 - a)
t = 16.d0*a*a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t +  &
    a8)*t + a7)*t + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5D0) + 0.5D0
fn_val = pi*a*w
GO TO 30

!                 A .LT. 0.25  OR  A .GT. 0.75

20 a = 0.25D0 + (0.75D0 - a)
n = n - 1
21 t = 16.d0*a*a
fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t + &
         b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t + 0.5D0) + 0.5D0

!                        TERMINATION

30 IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dcos1




FUNCTION drexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp) :: e, w, z
REAL (dp) :: a0 = .248015873015873015873016D-04,   &
    a1 = -.344452080605731005808147D-05, a2 = .206664230430046597475413D-06,  &
    a3 = -.447300111094328162971036D-08, a4 = .114734027080634968083920D-11,  &
    b1 = -.249994190011341852652396D+00, b2 = .249987228833107957725728D-01,  &
    b3 = -.119037506846942249362528D-02, b4 = .228908693387350391768682D-04
REAL (dp) :: c1 = .1666666666666666666666666666666667D+00,   &
             c2 = .4166666666666666666666666666666667D-01,   &
             c3 = .8333333333333333333333333333333333D-02,   &
             c4 = .1388888888888888888888888888888889D-02,   &
             c5 = .1984126984126984126984126984126984D-03
!---------------------------
IF (ABS(x) <= 0.15_dp) THEN

!        Z IS A MINIMAX APPROXIMATION OF THE SERIES

!                C6 + C7*X + C8*X**2 + ....

!        THIS APPROXIMATION IS ACCURATE TO WITHIN
!        1 UNIT OF THE 23-RD SIGNIFICANT DIGIT.
!        THE RESULTING VALUE FOR W IS ACCURATE TO
!        WITHIN 1 UNIT OF THE 33-RD SIGNIFICANT DIGIT.

  z = ((((a4*x + a3)*x + a2)*x + a1)*x + a0) /  &
      ((((b4*x + b3)*x + b2)*x + b1)*x + 1._dp)
  w = ((((((z*x + c5)*x + c4)*x + c3)*x + c2)*x + c1)*x + 0.5_dp)*x + 1._dp
  fn_val = x * w
  RETURN
END IF

IF (x >= 0.0_dp) THEN
  e = EXP(x)
  fn_val = e * (0.5_dp + (0.5_dp - 1.0_dp/e))
  RETURN
END IF
IF (x >= -77._dp) THEN
  fn_val = (EXP(x)-0.5_dp) - 0.5_dp
  RETURN
END IF
fn_val = -1._dp
RETURN
END FUNCTION drexp



SUBROUTINE dcgama(mo, z, w)
!-----------------------------------------------------------------------

!        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

!                        ---------------

!     MO IS AN INTEGER.  Z AND W ARE INTERPRETED AS REAL (dp)
!     COMPLEX NUMBERS.  IT IS ASSUMED THAT Z(1) AND Z(2) ARE THE REAL
!     AND IMAGINARY PARTS OF THE COMPLEX NUMBER Z, AND THAT W(1) AND
!     W(2) ARE THE REAL AND IMAGINARY PARTS OF W.

!                 W = GAMMA(Z)       IF MO = 0
!                 W = LN(GAMMA(Z))   OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)        :: mo
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: w

! Local variables
REAL (dp), PARAMETER :: c0(30)  &
        = (/ .8333333333333333333333333333333333333333D-01,  &
        -.2777777777777777777777777777777777777778D-02,  &
         .7936507936507936507936507936507936507937D-03,  &
        -.5952380952380952380952380952380952380952D-03,  &
         .8417508417508417508417508417508417508418D-03,  &
        -.1917526917526917526917526917526917526918D-02,  &
         .6410256410256410256410256410256410256410D-02,  &
        -.2955065359477124183006535947712418300654D-01,  &
         .1796443723688305731649384900158893966944D+00,  &
        -.1392432216905901116427432216905901116427D+01,  &
         .1340286404416839199447895100069013112491D+02,  &
        -.1568482846260020173063651324520889738281D+03,  &
         .2193103333333333333333333333333333333333D+04,  &
        -.3610877125372498935717326521924223073648D+05,  &
         .6914722688513130671083952507756734675533D+06,  &
        -.1523822153940741619228336495888678051866D+08,  &
         .3829007513914141414141414141414141414141D+09,  &
        -.1088226603578439108901514916552510537473D+11,  &
         .3473202837650022522522522522522522522523D+12,  &
        -.1236960214226927445425171034927132488108D+14,  &
         .4887880647930793350758151625180229021085D+15,  &
        -.2132033396091937389697505898213683855747D+17,  &
         .1021775296525700077565287628053585500394D+19,  &
        -.5357547217330020361082770919196920448485D+20,  &
         .3061578263704883415043151051329622758194D+22,  &
        -.1899991742639920405029371429306942902947D+24,  &
         .1276337403382883414923495137769782597654D+26,  &
        -.9252847176120416307230242348347622779519D+27,  &
         .7218822595185610297836050187301637922490D+29,  &
        -.6045183405995856967743148238754547286066D+31 /),  &
        dlpi = 1.144729885849400174143427351353058711647_dp,  &
        hl2p =  .9189385332046727417803297364056176398614_dp,  &
        pi = 3.141592653589793238462643383279502884197_dp,  &
        pi2 = 6.283185307179586476925286766559005768394_dp
REAL (dp) :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, q1, q2, s, sn,  &
             s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
INTEGER   :: j, k, max, n, nm1
!---------------------------
!     DLPI = LOG(PI)
!     HL2P = 0.5 * LOG(2*PI)
!---------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST NUMBER SUCH
!            THAT  1._dp + EPS > 1._dp.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0_dp)

!---------------------------
x = REAL(z, KIND=dp)
y = AIMAG(z)
IF (x < 0._dp) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi * y
  et = EXP(t)
  e2t = et * et

!     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

  a1 = 0.5_dp * (1._dp+e2t)
  t2 = t + t
  IF (t2 >= -0.15_dp) THEN
    a2 = -0.5_dp * drexp(t2)
  ELSE
    a2 = 0.5_dp * (0.5_dp+(0.5_dp-e2t))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X)

  u = MAX
  IF (ABS(x) >= MIN(u,1._dp/eps)) GO TO 80
  k = ABS(x)
  u = x + k
  k = MOD(k,2)
  IF (u <= -0.5_dp) THEN
    u = 0.5_dp + (0.5_dp+u)
    k = k + 1
  END IF
  u = pi * u
  sn = SIN(u)
  cn = COS(u)
  IF (k == 1) THEN
    sn = -sn
    cn = -cn
  END IF

!     SET  H1 + H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))

  a1 = sn * a1
  a2 = cn * a2
  a = a1 * a1 + a2 * a2
  IF (a == 0._dp) GO TO 80
  IF (mo == 0) THEN

    h1 = a1 / a
    h2 = -a2 / a
    c = pi * et
    h1 = c * h1
    h2 = c * h2
  ELSE

    h1 = (dlpi+t) - 0.5_dp * LOG(a)
    h2 = -ATAN2(a2,a1)
  END IF
  IF (AIMAG(z) >= 0._dp) THEN
    x = 1.0 - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0 - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
w1 = 0._dp
w2 = 0._dp
n = 0
t = x
y2 = y * y
a = t * t + y2
cut = 225._dp
IF (eps > 1.d-30) cut = 144._dp
IF (eps > 1.d-20) cut = 64._dp
IF (a < cut) THEN
  IF (a == 0._dp) GO TO 80
  10 n = n + 1
  t = t + 1._dp
  a = t * t + y2
  IF (a < cut) GO TO 10

!     LET S1 + S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)

  u1 = (x*t+y2) / a
  u2 = y / a
  s1 = u1
  s2 = n * u2
  IF (n >= 2) THEN
    u = t / a
    nm1 = n - 1
    DO j = 1, nm1
      v1 = u1 + j * u
      v2 = (n-j) * u2
      c = s1 * v1 - s2 * v2
      d = s1 * v2 + s2 * v1
      s1 = c
      s2 = d
    END DO
  END IF

!     SET  W1 + W2*I = LOG(S1 + S2*I)  WHEN MO IS NONZERO

  s = s1 * s1 + s2 * s2
  IF (mo /= 0) THEN
    w1 = 0.5_dp * LOG(s)
    w2 = ATAN2(s2,s1)
  END IF
END IF

!     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

t1 = 0.5_dp * LOG(a) - 1._dp
t2 = ATAN2(y,t)
u = x - 0.5_dp
v1 = (u*t1-0.5_dp) - y * t2
v2 = u * t2 + y * t1

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

u1 = t / a
u2 = -y / a
q1 = u1 * u1 - u2 * u2
q2 = 2._dp * u1 * u2
a1 = 0._dp
a2 = 0._dp
DO j = 1, 30
  t1 = a1
  t2 = a2
  a1 = a1 + c0(j) * u1
  a2 = a2 + c0(j) * u2
  IF (a1 == t1) THEN
    IF (a2 == t2) GO TO 40
  END IF
  t1 = u1 * q1 - u2 * q2
  t2 = u1 * q2 + u2 * q1
  u1 = t1
  u2 = t2
END DO
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
40 w1 = (((a1+hl2p)-w1)+v1) - n
w2 = (a2-w2) + v2
IF (REAL(z, KIND=dp) < 0.0_dp) GO TO 60
IF (mo == 0) THEN

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

  a = EXP(w1)
  w1 = a * COS(w2)
  w2 = a * SIN(w2)
  IF (n == 0) GO TO 70
  c = (s1*w1+s2*w2) / s
  d = (s1*w2-s2*w1) / s
  w1 = c
  w2 = d
  GO TO 70
END IF

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
!     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

50 IF (w2 <= pi) THEN
  k = 0.5_dp - w2 / pi2
  w2 = w2 + pi2 * k
  GO TO 70
END IF
k = w2 / pi2 - 0.5_dp
u = k + 1
w2 = w2 - pi2 * u
IF (w2 <= -pi) w2 = pi
GO TO 70

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

60 IF (mo /= 0) THEN
  w1 = h1 - w1
  w2 = h2 - w2
  GO TO 50
END IF

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO = 0

a = EXP(-w1)
t1 = a * COS(-w2)
t2 = a * SIN(-w2)
w1 = h1 * t1 - h2 * t2
w2 = h1 * t2 + h2 * t1
IF (n /= 0) THEN
  c = w1 * s1 - w2 * s2
  d = w1 * s2 + w2 * s1
  w1 = c
  w2 = d
END IF

!     TERMINATION

70 w = CMPLX(w1, w2, KIND=dp)
RETURN
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
80 w = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
RETURN
END SUBROUTINE dcgama




FUNCTION cgam0(z) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF 1/GAMMA(1 + Z)  FOR ABS(Z) < 1.0
!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)  :: z
COMPLEX (dp)              :: fn_val

COMPLEX (dp)  :: w
INTEGER       :: i, k, n
!-----------------------
REAL (dp)  :: x, y
REAL, PARAMETER :: a(25) = (/ .577215664901533_dp, -.655878071520254_dp,  &
        -.420026350340952D-01, .166538611382291_dp, -.421977345555443D-01,  &
        -.962197152787697D-02, .721894324666310D-02, -.116516759185907D-02,  &
        -.215241674114951D-03, .128050282388116D-03, -.201348547807882D-04,  &
        -.125049348214267D-0, .113302723198170D-05, -.205633841697761D-0,  &
        .611609510448142D-08, .500200764446922D-08, -.118127457048702D-08,  &
        .104342671169110D-09, .778226343990507D-11, -.369680561864221D-11,  &
        .510037028745448D-12, -.205832605356651D-13, -.534812253942302D-14,  &
        .122677862823826D-14, -.118125930169746D-15 /)
!-----------------------
n = 25
x = REAL(z, KIND=dp)
y = AIMAG(z)
IF (x*x+y*y <= 0.04D0) n = 14

k = n
w = a(n)
DO  i = 2, n
  k = k - 1
  w = a(k) + z * w
END DO
fn_val = 1.0D0 + z * w
RETURN
END FUNCTION cgam0



FUNCTION dgamma(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!                EVALUATION OF THE GAMMA FUNCTION FOR
!                     REAL (dp) ARGUMENTS

!                           -----------

!     DGAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.41893853320467274178032973640562_dp,  &
                        pi = 3.14159265358979323846264338327950_dp
REAL (dp) :: s, t, x, w
INTEGER   :: j, n
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
fn_val = 0.0_dp
x = a
IF (ABS(a) <= 20._dp) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DGAMMA(A) FOR ABS(A) <= 20
!-----------------------------------------------------------------------
  t = 1.0_dp
  n = x
  n = n - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (n < 0) THEN
    GO TO 40
  ELSE IF (n == 0) THEN
    GO TO 30
  END IF

  DO j = 1, n
    x = x - 1._dp
    t = x * t
  END DO
  30 x = x - 1._dp
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0._dp) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1._dp
        t = x * t
      END DO
    END IF
    x = (x+0.5_dp) + 0.5_dp
    t = x * t
    IF (t == 0._dp) RETURN
  END IF

!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.d-33) THEN
    IF (ABS(t)*HUGE(1.0_dp) <= 1.000000001_dp) RETURN
    fn_val = 1._dp / t
    RETURN
  END IF

!     COMPUTE DGAMMA(1 + X) FOR 0 <= X < 1

  60 fn_val = 1._dp / (1._dp + dgam1(x))

!     TERMINATION

  IF (a >= 1._dp) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DGAMMA(A) FOR ABS(A) > 20
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.d3) RETURN
IF (a <= 0._dp) THEN
  s = dsin1(a) / pi
  IF (s == 0._dp) RETURN
  x = -a
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

w = dpdel(x)

!     FINAL ASSEMBLY

w = (d+w) + (x-0.5_dp) * (LOG(x)-1._dp)
IF (w > dxparg(0)) RETURN
fn_val = EXP(w)
IF (a < 0._dp) fn_val = (1._dp/(fn_val*s)) / x

RETURN
END FUNCTION dgamma



FUNCTION glog(x) RESULT(fn_val)
!     -------------------
!     EVALUATION OF LN(X) FOR X >= 15
!     -------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c1 = .286228750476730_dp, c2 = .399999628131494+dp,  &
                        c3 = .666666666752663_dp
REAL (dp), PARAMETER :: w(163) = (/ .270805020110221007D+01, .277258872223978124D+01,  &
  .283321334405621608D+01, .289037175789616469D+01, .294443897916644046D+01, .299573227355399099D+01, .304452243772342300D+01, &
  .309104245335831585D+01, .313549421592914969D+01, .317805383034794562D+01, .321887582486820075D+01, .325809653802148205D+01, &
  .329583686600432907D+01, .333220451017520392D+01, .336729582998647403D+01, .340119738166215538D+01, .343398720448514625D+01, &
  .346573590279972655D+01, .349650756146648024D+01, .352636052461616139D+01, .355534806148941368D+01, .358351893845611000D+01, &
  .361091791264422444D+01, .363758615972638577D+01, .366356164612964643D+01, .368887945411393630D+01, .371357206670430780D+01, &
  .373766961828336831D+01, .376120011569356242D+01, .378418963391826116D+01, .380666248977031976D+01, .382864139648909500D+01, &
  .385014760171005859D+01, .387120101090789093D+01, .389182029811062661D+01, .391202300542814606D+01, .393182563272432577D+01, &
  .395124371858142735D+01, .397029191355212183D+01, .398898404656427438D+01, .400733318523247092D+01, .402535169073514923D+01, &
  .404305126783455015D+01, .406044301054641934D+01, .407753744390571945D+01, .409434456222210068D+01, .411087386417331125D+01, &
  .412713438504509156D+01, .414313472639153269D+01, .415888308335967186D+01, .417438726989563711D+01, .418965474202642554D+01, &
  .420469261939096606D+01, .421950770517610670D+01, .423410650459725938D+01, .424849524204935899D+01, .426267987704131542D+01, &
  .427666611901605531D+01, .429045944114839113D+01, .430406509320416975D+01, .431748811353631044D+01, .433073334028633108D+01, &
  .434380542185368385D+01, .435670882668959174D+01, .436944785246702149D+01, .438202663467388161D+01, .439444915467243877D+01, &
  .440671924726425311D+01, .441884060779659792D+01, .443081679884331362D+01, .444265125649031645D+01, .445434729625350773D+01, &
  .446590811865458372D+01, .447733681447820647D+01, .448863636973213984D+01, .449980967033026507D+01, .451085950651685004D+01, &
  .452178857704904031D+01, .453259949315325594D+01, .454329478227000390D+01, .455387689160054083D+01, .456434819146783624D+01, &
  .457471097850338282D+01, .458496747867057192D+01, .459511985013458993D+01, .460517018598809137D+01, .461512051684125945D+01, &
  .462497281328427108D+01, .463472898822963577D+01, .464439089914137266D+01, .465396035015752337D+01, .466343909411206714D+01, &
  .467282883446190617D+01, .468213122712421969D+01, .469134788222914370D+01, .470048036579241623D+01, .470953020131233414D+01, &
  .471849887129509454D+01, .472738781871234057D+01, .473619844839449546D+01, .474493212836325007D+01, .475359019110636465D+01, &
  .476217393479775612D+01, .477068462446566476D+01, .477912349311152939D+01, .478749174278204599D+01, .479579054559674109D+01, &
  .480402104473325656D+01, .481218435537241750D+01, .482028156560503686D+01, .482831373730230112D+01, .483628190695147800D+01, &
  .484418708645859127D+01, .485203026391961717D+01, .485981240436167211D+01, .486753445045558242D+01, .487519732320115154D+01, &
  .488280192258637085D+01, .489034912822175377D+01, .489783979995091137D+01, .490527477843842945D+01, .491265488573605201D+01, &
  .491998092582812492D+01, .492725368515720469D+01, .493447393313069176D+01, .494164242260930430D+01, .494875989037816828D+01, &
  .495582705760126073D+01, .496284463025990728D+01, .496981329957600062D+01, .497673374242057440D+01, .498360662170833644D+01, &
  .499043258677873630D+01, .499721227376411506D+01, .500394630594545914D+01, .501063529409625575D+01, .501727983681492433D+01, &
  .502388052084627639D+01, .503043792139243546D+01, .503695260241362916D+01, .504342511691924662D+01, .504985600724953705D+01, &
  .505624580534830806D+01, .506259503302696680D+01, .506890420222023153D+01, .507517381523382692D+01, .508140436498446300D+01, &
  .508759633523238407D+01, .509375020080676233D+01, .509986642782419842D+01, .510594547390058061D+01, .511198778835654323D+01, &
  .511799381241675511D+01, .512396397940325892D+01, .512989871492307347D+01, .513579843705026176D+01, .514166355650265984D+01, &
  .514749447681345304D+01, .515329159449777895D+01, .515905529921452903D+01, .516478597392351405D+01, .517048399503815178D+01, &
  .517614973257382914D+01 /)
REAL (dp) :: t, t2, z
INTEGER   :: n
!     -------------------

IF (x < 178.0_dp) THEN
  n = x
  t = (x-n) / (x+n)
  t2 = t * t
  z = (((c1*t2 + c2)*t2 + c3)*t2 + 2.0) * t
  fn_val = w(n-14) + z
  RETURN
END IF

fn_val = LOG(x)
RETURN
END FUNCTION glog



SUBROUTINE cbsslj(z,cnu,w)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEX BESSEL FUNCTION J   (Z)
!                                                    CNU
!-----------------------------------------------------------------------

!     WRITTEN BY
!         ANDREW H. VAN TUYL AND ALFRED H. MORRIS, JR.
!         NAVAL SURFACE WARFARE CENTER
!         OCTOBER, 1991

!     A MODIFICATION OF THE PROCEDURE DEVELOPED BY ALLEN V. HERSHEY
!     (NAVAL SURFACE WARFARE CENTER) IN 1978 FOR HANDLING THE DEBYE
!     APPROXIMATION IS EMPLOYED.

!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(IN)   :: cnu
COMPLEX (dp), INTENT(OUT)  :: w

COMPLEX (dp)  :: c, nu, s, sm1, sm2, t, tsc, w0, w1, zn, zz
!-----------------------
REAL (dp) :: a, cn1, cn2, e, fn
REAL (dp) :: pn, qm, qn, qnp1
REAL (dp) :: r, rn2, r2, sn, t1, t2
REAL (dp) :: u, v, x, y
INTEGER   :: i, k, m, n
REAL (dp), PARAMETER  :: pi = 3.141592653589793238462643383279502884197_dp
!-----------------------
x = REAL(z, KIND=dp)
y = AIMAG(z)
r = cpabs(x,y)
cn1 = REAL(cnu, KIND=dp)
cn2 = AIMAG(cnu)
rn2 = cn1 * cn1 + cn2 * cn2
pn = INT(cn1)
fn = cn1 - pn
sn = 1.0_dp

!          CALCULATION WHEN ORDER IS AN INTEGER

IF (fn == 0.0_dp .AND. cn2 == 0.0_dp) THEN
  n = pn
  pn = ABS(pn)
  cn1 = pn
  IF (n < 0 .AND. n /= (n/2)*2) sn = -1.0_dp
END IF

!          SELECTION OF METHOD

IF (r > 17.5D0) THEN
  IF (r > 17.5D0 + 0.5D0*rn2) GO TO 10
  GO TO 20
END IF

!          USE MACLAURIN EXPANSION AND RECURSION

IF (cn1 < 0.0D0) THEN
  qn = -1.25D0 * (r + 0.5D0*ABS(cn2) - ABS(y-0.5D0*cn2))
  IF (cn1 < qn) THEN
    qn = 1.25D0 * (r - MAX(1.2D0*r,ABS(y-cn2)))
    IF (cn1 < qn) THEN
      qn = MIN(pn, REAL(-INT(1.25D0*(r-ABS(cn2))), KIND=dp))
      GO TO 60
    END IF
  END IF
END IF

r2 = r * r
qm = 0.0625D0 * r2 * r2 - cn2 * cn2
qn = MAX(pn, REAL(INT(SQRT(MAX(0.0D0,qm))), KIND=dp))
GO TO 60

!          USE ASYMPTOTIC EXPANSION

10 CALL cbja(z,cnu,w)
RETURN

!          CALCULATION FOR 17.5 < ABS(Z) <= 17.5 + 0.5*ABS(CNU)**2

20 n = 0
IF (ABS(cn2) < 0.8D0*ABS(y)) THEN
  qm = -1.25D0 * (r + 0.5D0*ABS(cn2) - ABS(y-0.5D0*cn2))
  IF (cn1 < qm) THEN
    qm = 1.25D0 * (r - MAX(1.2D0*r, ABS(y-cn2)))
    IF (cn1 < qm) n = 1
  END IF
END IF

qn = pn
a = 4.d-3 * r * r
zz = z
IF (x < 0.0D0) zz = -z

!          CALCULATION OF ZONE OF EXCLUSION OF DEBYE APPROXIMATION

30 nu = CMPLX(qn+fn,cn2, KIND=dp)
zn = nu / z
t2 = AIMAG(zn) * AIMAG(zn)
u = 1.0D0 - REAL(zn, KIND=dp)
t1 = u * u + t2
u = 1.0D0 + DBLE(zn)
t2 = u * u + t2
u = t1 * t2
v = a * u / (t1*t1 + t2*t2)
IF (u*v*v <= 1.0D0) THEN
  
!          THE ARGUMENT LIES INSIDE THE ZONE OF EXCLUSION
  
  qn = qn + 1.0D0
  IF (n == 0) GO TO 30
  
!          USE MACLAURIN EXPANSION WITH FORWARD RECURRENCE
  
  qn = MIN(pn, REAL(-INT(1.25D0*(r-ABS(cn2))), KIND=dp))
ELSE
  
!          USE BACKWARD RECURRENCE STARTING FROM THE ASYMPTOTIC EXPANSION
  
  qnp1 = qn + 1.0D0
  IF (ABS(qn) < ABS(pn)) THEN
    IF (r >= 17.5D0 + 0.5D0*(qnp1*qnp1 + cn2*cn2)) THEN
      
      nu = CMPLX(qn+fn,cn2, KIND=dp)
      CALL cbja(zz,nu,sm1)
      nu = CMPLX(qnp1+fn,cn2, KIND=dp)
      CALL cbja(zz,nu,sm2)
      GO TO 40
    END IF
  END IF
  
!          USE BACKWARD RECURRENCE STARTING FROM THE DEBYE APPROXIMATION
  
  nu = CMPLX(qn+fn,cn2, KIND=dp)
  CALL cbdb(zz,nu,fn,sm1)
  IF (qn == pn) GO TO 50
  nu = CMPLX(qnp1+fn,cn2, KIND=dp)
  CALL cbdb(zz,nu,fn,sm2)
  
  40 nu = CMPLX(qn+fn,cn2, KIND=dp)
  tsc = 2.0D0 * nu * sm1 / zz - sm2
  sm2 = sm1
  sm1 = tsc
  qn = qn - 1.0D0
  IF (qn /= pn) GO TO 40
  
  50 w = sm1
  IF (sn < 0.0D0) w = -w
  IF (x >= 0.0D0) RETURN
  
  nu = pi * CMPLX(-cn2,cn1, KIND=dp)
  IF (y < 0.0D0) nu = -nu
  w = EXP(nu) * w
  RETURN
END IF

!          USE MACLAURIN EXPANSION WITH FORWARD OR BACKWARD RECURRENCE.

60 m = qn - pn
IF (ABS(m) <= 1) THEN
  nu = CMPLX(cn1,cn2, KIND=dp)
  CALL cbjm(z,nu,w)
ELSE
  nu = CMPLX(qn+fn,cn2, KIND=dp)
  CALL cbjm(z,nu,w1)
  w0 = 0.25D0 * z * z
  IF (m <= 0) THEN
    
!          FORWARD RECURRENCE
    
    m = ABS(m)
    nu = nu + 1.0D0
    CALL cbjm(z,nu,w)
    DO  i = 2, m
      c = nu * (nu+1.0D0)
      t = (c/w0) * (w-w1)
      w1 = w
      w = t
      nu = nu + 1.0D0
    END DO
  ELSE
    
!          BACKWARD RECURRENCE
    
    nu = nu - 1.0D0
    CALL cbjm(z,nu,w)
    DO  i = 2, m
      c = nu * (nu+1.0D0)
      t = (w0/c) * w1
      w1 = w
      w = w - t
      nu = nu - 1.0D0
    END DO
  END IF
END IF

!          FINAL ASSEMBLY

IF (fn == 0.0D0 .AND. cn2 == 0.0D0) THEN
  k = pn
  IF (k == 0) RETURN
  e = sn / dgamma(pn+1.0D0)
  w = e * w * (0.5D0*z) ** k
  RETURN
END IF

s = cnu * LOG(0.5D0*z)
w = EXP(s) * w
IF (rn2 <= 0.81D0) THEN
  w = w * cgam0(cnu)
  RETURN
END IF
CALL dcgama(0,cnu,t)
w = cdiv(w,cnu*t)
RETURN
END SUBROUTINE cbsslj



SUBROUTINE cbjm(z,cnu,w)
!-----------------------------------------------------------------------

!       COMPUTATION OF  (Z/2)**(-CNU) * GAMMA(CNU + 1) * J(CNU,Z)

!                           -----------------

!     THE MACLAURIN EXPANSION IS USED.  IT IS ASSUMED THAT CNU IS NOT
!     A NEGATIVE INTEGER.

!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(IN)   :: cnu
COMPLEX (dp), INTENT(OUT)  :: w

COMPLEX (dp)  :: nu, nup1, p, s, sn, t, ti
!--------------------------
REAL (dp)  :: a, a0, eps, inu, m, rnu
INTEGER    :: i, imin, k, km1, km2

!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0 .

eps = EPSILON(0.0_dp)

!--------------------------
s = -0.25D0 * (z*z)
nu = cnu
rnu = REAL(nu, KIND=dp)
inu = AIMAG(nu)
a = 0.5D0 + (0.5D0+rnu)
nup1 = CMPLX(a,inu, KIND=dp)

IF (a > 0.0D0) THEN
  m = 1.0D0
  t = s / nup1
  w = 1.0D0 + t
ELSE
  
!     ADD 1.0 AND THE FIRST K-1 TERMS
  
  k = INT(-a) + 2
  km1 = k - 1
  w = (1.0D0,0.0D0)
  t = w
  DO  i = 1, km1
    m = i
    t = t * (s/(m*(nu+m)))
    w = w + t
    IF (anorm(t) <= eps*anorm(w)) GO TO 20
  END DO
  GO TO 50
  
!     CHECK IF THE (K-1)-ST AND K-TH TERMS CAN BE IGNORED.
!     IF SO THEN THE SUMMATION IS COMPLETE.
  
  20 IF (i /= km1) THEN
    imin = i + 1
    IF (imin < k-5) THEN
      ti = t
      
      m = km1
      t = s / (nu+m)
      a0 = anorm(t) / m
      t = t * (s/(nu+(m+1.0D0)))
      a = anorm(t) / (m*(m+1.0D0))
      a = MAX(a,a0)
      
      t = (1.0D0,0.0D0)
      km2 = k - 2
      DO  i = imin, km2
        m = i
        t = t * (s/(m*(nu+m)))
        IF (a*anorm(t) < 0.5D0) RETURN
      END DO
      t = t * ti
      imin = km2
    END IF
    
!     ADD THE (K-1)-ST TERM
    
    a = 1.0D0
    p = (1.0D0,0.0D0)
    sn = p
    DO  i = imin, km1
      m = i
      a = a * m
      p = p * (nu+m)
      sn = s * sn
    END DO
    t = t * (cdiv(sn,p)/a)
    w = w + t
  END IF
END IF

!     ADD THE REMAINING TERMS

50 m = m + 1.0D0
t = t * (s/(m*(nu+m)))
w = w + t
IF (anorm(t) > eps*anorm(w)) GO TO 50

RETURN
END SUBROUTINE cbjm



SUBROUTINE cbdb(cz,cnu,fn,w)
!-----------------------------------------------------------------------

!         CALCULATION OF J   (CZ) BY THE DEBYE APPROXIMATION
!                         CNU
!                         ------------------

!     IT IS ASSUMED THAT REAL(CZ) .GE. 0 AND THAT REAL(CNU) = FN + K
!     WHERE K IS AN INTEGER.

!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)   :: cz, cnu
REAL (dp), INTENT(IN)      :: fn
COMPLEX (dp), INTENT(OUT)  :: w

! Local variables
REAL (dp)     :: is, inu, izn
COMPLEX (dp)  :: c1, c2, eta, nu, p, p1, q, r, s, s1, s2, sm, t, z, zn
!----------------------
!     C = 1/SQRT(2*PI)
!     BND = PI/3
!----------------------
REAL (dp), PARAMETER  :: c = .398942280401433_dp, pi = 3.14159265358979_dp,  &
                         pi2 = 6.28318530717959_dp, bnd = 1.04719755119660_dp
COMPLEX (dp), PARAMETER :: j = (0.0, 1.0)
REAL (dp)  :: alpha, am, aq, ar
REAL (dp)  :: phi, sgn, theta
REAL (dp)  :: u, v, x, y
INTEGER    :: ind, k, l, m

!----------------------
!             COEFFICIENTS OF THE FIRST 16 POLYNOMIALS
!                   IN THE DEBYE APPROXIMATION
!----------------------

REAL (dp)  :: a(136) = (/ 1.0_dp, -.208333333333333_dp, .125000000000000_dp, .334201388888889_dp, &
  -.401041666666667_dp, .703125000000000D-01,-.102581259645062D+01, .184646267361111D+01, &
  -.891210937500000_dp, .732421875000000D-01, .466958442342625D+01,-.112070026162230D+02, &
   .878912353515625D+01,-.236408691406250D+01, .112152099609375_dp,-.282120725582002D+02, &
   .846362176746007D+02,-.918182415432400D+02, .425349987453885D+02,-.736879435947963D+01, &
   .227108001708984_dp, .212570130039217D+03,-.765252468141182D+03, .105999045252800D+04, &
  -.699579627376133D+03, .218190511744212D+03,-.264914304869516D+02, .572501420974731_dp, &
  -.191945766231841D+04, .806172218173731D+04,-.135865500064341D+05, .116553933368645D+05, &
  -.530564697861340D+04, .120090291321635D+04,-.108090919788395D+03, .172772750258446D+01, &
   .202042913309661D+05,-.969805983886375D+05, .192547001232532D+06,-.203400177280416D+06, &
   .122200464983017D+06,-.411926549688976D+05, .710951430248936D+04,-.493915304773088D+03, &
   .607404200127348D+01,-.242919187900551D+06, .131176361466298D+07,-.299801591853811D+07, &
   .376327129765640D+07,-.281356322658653D+07, .126836527332162D+07,-.331645172484564D+06, &
   .452187689813627D+05,-.249983048181121D+04, .243805296995561D+02, .328446985307204D+07, &
  -.197068191184322D+08, .509526024926646D+08,-.741051482115327D+08, .663445122747290D+08, &
  -.375671766607634D+08, .132887671664218D+08,-.278561812808645D+07, .308186404612662D+06, &
  -.138860897537170D+05, .110017140269247D+03,-.493292536645100D+08, .325573074185766D+09, &
  -.939462359681578D+09, .155359689957058D+10,-.162108055210834D+10, .110684281682301D+10, &
  -.495889784275030D+09, .142062907797533D+09,-.244740627257387D+08, .224376817792245D+07, &
  -.840054336030241D+05, .551335896122021D+03, .814789096118312D+09,-.586648149205185D+10, &
   .186882075092958D+11,-.346320433881588D+11, .412801855797540D+11,-.330265997498007D+11, &
   .179542137311556D+11,-.656329379261928D+10, .155927986487926D+10,-.225105661889415D+09, &
   .173951075539782D+08,-.549842327572289D+06, .303809051092238D+04,-.146792612476956D+11, &
   .114498237732026D+12,-.399096175224466D+12, .819218669548577D+12,-.109837515608122D+13, &
   .100815810686538D+13,-.645364869245377D+12, .287900649906151D+12,-.878670721780233D+11, &
   .176347306068350D+11,-.216716498322380D+10, .143157876718889D+09,-.387183344257261D+07, &
   .182577554742932D+05, .286464035717679D+12,-.240629790002850D+13, .910934118523990D+13, &
  -.205168994109344D+14, .305651255199353D+14,-.316670885847852D+14, .233483640445818D+14, &
  -.123204913055983D+14, .461272578084913D+13,-.119655288019618D+13, .205914503232410D+12, &
  -.218229277575292D+11, .124700929351271D+10,-.291883881222208D+08, .118838426256783D+06, &
  -.601972341723401D+13, .541775107551060D+14,-.221349638702525D+15, .542739664987660D+15, &
  -.889496939881026D+15, .102695519608276D+16,-.857461032982895D+15, .523054882578445D+15, &
  -.232604831188940D+15, .743731229086791D+14,-.166348247248925D+14, .248500092803409D+13, &
  -.229619372968246D+12, .114657548994482D+11,-.234557963522252D+09, .832859304016289D+06 /)

z = cz
nu = cnu
inu = AIMAG(cnu)
IF (inu < 0.0D0) THEN
  z = CONJG(z)
  nu = CONJG(nu)
END IF
x = REAL(z, KIND=dp)
y = AIMAG(z)

!          TANH(GAMMA) = SQRT(1 - (Z/NU)**2) = W/NU
!          T = EXP(NU*(TANH(GAMMA) - GAMMA))

zn = z / nu
izn = AIMAG(zn)
IF (ABS(izn) <= 0.1D0*ABS(REAL(zn, KIND=dp))) THEN
  
  s = (1.0D0-zn) * (1.0D0+zn)
  eta = 1.0D0 / s
  q = SQRT(s)
  s = 1.0D0 / (nu*q)
  t = zn / (1.0D0 + q)
  t = EXP(nu*(q + LOG(t)))
ELSE
  
  s = (nu-z) * (nu+z)
  eta = (nu*nu) / s
  w = SQRT(s)
  q = w / nu
  IF (REAL(q, KIND=dp) < 0.0D0) w = -w
  s = 1.0D0 / w
  t = z / (nu+w)
  t = EXP(w + nu*LOG(t))
END IF

is = AIMAG(s)
r = SQRT(s)
c1 = r * t
ar = REAL(r, KIND=dp) * REAL(r, KIND=dp) + AIMAG(r) * AIMAG(r)
aq = -1.0D0 / (REAL(q, KIND=dp)*REAL(q, KIND=dp) + AIMAG(q)*AIMAG(q))

phi = ATAN2(y,x) / 3.0D0
q = nu - z
theta = ATAN2(AIMAG(q),REAL(q, KIND=dp)) - phi
ind = 0
IF (ABS(theta) >= 2.0D0*bnd) THEN
  
  ind = 1
  CALL dcrec(REAL(t, KIND=dp), AIMAG(t),u,v)
  c2 = -j * r * CMPLX(u,v, KIND=dp)
  IF (is >= 0.0D0) THEN
    IF (is > 0.0D0) GO TO 10
    IF (REAL(s, KIND=dp) <= 0.0D0) GO TO 10
  END IF
  c2 = -c2
END IF

!          SUMMATION OF THE SERIES S1 AND S2

10 sm = s * s
p = (a(2)*eta + a(3)) * s
p1 = ((a(4)*eta + a(5))*eta + a(6)) * sm
s1 = (1.0D0 + p) + p1
IF (ind /= 0) s2 = (1.0D0-p) + p1
sgn = 1.0D0
am = ar * ar
m = 4
l = 6

!          P = VALUE OF THE M-TH POLYNOMIAL

20 l = l + 1
alpha = a(l)
p = CMPLX(a(l),0.0D0, KIND=dp)
DO  k = 2, m
  l = l + 1
  alpha = a(l) + aq * alpha
  p = a(l) + eta * p
END DO

!          ONLY THE S1 SUM IS FORMED WHEN IND = 0

sm = s * sm
p = p * sm
s1 = s1 + p
IF (ind /= 0) THEN
  sgn = -sgn
  s2 = s2 + sgn * p
END IF
am = ar * am
IF (1.0D0 + alpha*am /= 1.0D0) THEN
  m = m + 1
  IF (m <= 16) GO TO 20
END IF

!          FINAL ASSEMBLY

s1 = c * c1 * s1
IF (ind == 0) THEN
  w = s1
ELSE
  
  s2 = c * c2 * s2
  q = nu + z
  theta = ATAN2(AIMAG(q),REAL(q, KIND=dp)) - phi
  IF (ABS(theta) <= bnd) THEN
    w = s1 + s2
  ELSE
    
    alpha = pi2
    IF (izn < 0.0D0) alpha = -alpha
    t = alpha * CMPLX(ABS(inu),-fn, KIND=dp)
    alpha = EXP(REAL(t))
    u = AIMAG(t)
    r = CMPLX(COS(u),SIN(u), KIND=dp)
    t = s1 - (alpha*r) * s1
    IF (x == 0.0D0 .AND. inu == 0.0D0) t = -t
    
    IF (y < 0.0D0) THEN
      IF (izn >= 0.0D0 .AND. theta <= SIGN(pi,theta)) s2 =  &
                                                      s2 * ( CONJG(r)/alpha)
      IF (x == 0.0D0) GO TO 40
      IF (izn >= 0.0D0) THEN
        IF (is < 0.0D0) GO TO 40
      END IF
    END IF
    
    w = s2 + t
    GO TO 50
    40 w = s2 - t
  END IF
END IF

50 IF (inu < 0.0D0) w = CONJG(w)
RETURN
END SUBROUTINE cbdb



SUBROUTINE cbja(cz,cnu,w)
!-----------------------------------------------------------------------
!        COMPUTATION OF J(NU,Z) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------

COMPLEX (dp), INTENT(IN)   :: cz
COMPLEX (dp), INTENT(IN)   :: cnu
COMPLEX (dp), INTENT(OUT)  :: w

! Local variables
REAL (dp)     :: eps, inu, m
COMPLEX (dp)  :: a, a1, arg, e, eta, nu, p, q, t, z, zr, zz
!--------------------------
REAL (dp) :: r, rnu, tol, u, v
REAL (dp) :: x, y
INTEGER   :: i, ind

!--------------------------
!     PIHALF = PI/2
!     C = 2*PI**(-1/2)
!--------------------------
REAL (dp), PARAMETER    :: pihalf = 1.5707963267949_dp, c = 1.12837916709551_dp
COMPLEX (dp), PARAMETER :: j = (0.0_dp, 1.0_dp)
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0 .

eps = EPSILON(0.0_dp)

!--------------------------
z = cz
x = REAL(z, KIND=dp)
y = AIMAG(z)
nu = cnu
ind = 0
IF (ABS(x) <= 1.d-2*ABS(y)) THEN
  IF (AIMAG(nu) < 0.0D0 .AND. ABS(REAL(nu)) < 1.d-2*ABS(AIMAG(nu))) THEN
    ind = 1
    nu = CONJG(nu)
    z = CONJG(z)
    y = -y
  END IF
END IF

IF (x < -1.d-2*y) z = -z
zz = z + z
CALL dcrec(REAL(zz, KIND=dp),AIMAG(zz),u,v)
zr = CMPLX(u,v, KIND=dp)
eta = -zr * zr

p = (0.0D0,0.0D0)
q = (0.0D0,0.0D0)
a1 = nu * nu - 0.25D0
a = a1
t = a1
m = 1.0D0
tol = eps * anorm(a1)
DO  i = 1, 16
  a = a - 2.0D0 * m
  m = m + 1.0D0
  t = t * a * eta / m
  p = p + t
  a = a - 2.0D0 * m
  m = m + 1.0D0
  t = t * a / m
  q = q + t
  IF (anorm(t) <= tol) GO TO 20
END DO

20 p = p + 1.0D0
q = (q+a1) * zr
w = z - pihalf * nu
IF (ABS(AIMAG(w)) <= 1.0D0) THEN
  arg = w - 0.5D0 * pihalf
  w = c * SQRT(zr) * (p*COS(arg) - q*SIN(arg))
ELSE
  e = EXP(-j*w)
  t = q - j * p
  IF (AIMAG(z) > 0.0D0 .AND. REAL(z, KIND=dp) <= 1.d-2*AIMAG(z).AND.  &
      ABS(REAL(nu, KIND=dp)) < 1.d-2*AIMAG(nu)) t = 0.5D0 * t
  CALL dcrec(REAL(e, KIND=dp),AIMAG(e),u,v)
  w = 0.5D0 * c * SQRT(j*zr) * ((p-j*q)*e + t*CMPLX(u,v, KIND=dp))
END IF

IF (x < -1.d-2*y) THEN
  IF (y < 0.0D0) nu = -nu
  
!     COMPUTATION OF EXP(I*PI*NU)
  
  rnu = REAL(nu, KIND=dp)
  inu = AIMAG(nu)
  r = EXP(-2.0D0*pihalf*inu)
  u = r * dcos1(rnu)
  v = r * dsin1(rnu)
  w = w * CMPLX(u,v, KIND=dp)
END IF

IF (ind /= 0) w = CONJG(w)
RETURN
END SUBROUTINE cbja



FUNCTION anorm(z) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

COMPLEX (dp), INTENT(IN)  :: z
REAL (dp)                 :: fn_val

fn_val = MAX( ABS( REAL(z, KIND=dp)), ABS(AIMAG(z) ) )
RETURN
END FUNCTION anorm



FUNCTION dgam1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF 1/GAMMA(1 + X) - 1  FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------

!     THE FOLLOWING ARE THE FIRST 49 COEFFICIENTS OF THE MACLAURIN
!     EXPANSION FOR 1/GAMMA(1 + X) - 1. THE COEFFICIENTS ARE
!     CORRECT TO 40 DIGITS. THE COEFFICIENTS WERE OBTAINED BY
!     ALFRED H. MORRIS JR. (NAVAL SURFACE WARFARE CENTER) AND ARE
!     GIVEN HERE FOR REFERENCE. ONLY THE FIRST 14 COEFFICIENTS ARE
!     USED IN THIS CODE.

!                           -----------

!     DATA A(1)  / .5772156649015328606065120900824024310422D+00/,
!    *     A(2)  /-.6558780715202538810770195151453904812798D+00/,
!    *     A(3)  /-.4200263503409523552900393487542981871139D-01/,
!    *     A(4)  / .1665386113822914895017007951021052357178D+00/,
!    *     A(5)  /-.4219773455554433674820830128918739130165D-01/,
!    *     A(6)  /-.9621971527876973562114921672348198975363D-02/,
!    *     A(7)  / .7218943246663099542395010340446572709905D-02/,
!    *     A(8)  /-.1165167591859065112113971084018388666809D-02/,
!    *     A(9)  /-.2152416741149509728157299630536478064782D-03/,
!    *     A(10) / .1280502823881161861531986263281643233949D-03/
!     DATA A(11) /-.2013485478078823865568939142102181838229D-04/,
!    *     A(12) /-.1250493482142670657345359473833092242323D-05/,
!    *     A(13) / .1133027231981695882374129620330744943324D-05/,
!    *     A(14) /-.2056338416977607103450154130020572836513D-06/,
!    *     A(15) / .6116095104481415817862498682855342867276D-08/,
!    *     A(16) / .5002007644469222930055665048059991303045D-08/,
!    *     A(17) /-.1181274570487020144588126565436505577739D-08/,
!    *     A(18) / .1043426711691100510491540332312250191401D-09/,
!    *     A(19) / .7782263439905071254049937311360777226068D-11/,
!    *     A(20) /-.3696805618642205708187815878085766236571D-11/
!     DATA A(21) / .5100370287454475979015481322863231802727D-12/,
!    *     A(22) /-.2058326053566506783222429544855237419746D-13/,
!    *     A(23) /-.5348122539423017982370017318727939948990D-14/,
!    *     A(24) / .1226778628238260790158893846622422428165D-14/,
!    *     A(25) /-.1181259301697458769513764586842297831212D-15/,
!    *     A(26) / .1186692254751600332579777242928674071088D-17/,
!    *     A(27) / .1412380655318031781555803947566709037086D-17/,
!    *     A(28) /-.2298745684435370206592478580633699260285D-18/,
!    *     A(29) / .1714406321927337433383963370267257066813D-19/,
!    *     A(30) / .1337351730493693114864781395122268022875D-21/
!     DATA A(31) /-.2054233551766672789325025351355733796682D-21/,
!    *     A(32) / .2736030048607999844831509904330982014865D-22/,
!    *     A(33) /-.1732356445910516639057428451564779799070D-23/,
!    *     A(34) /-.2360619024499287287343450735427531007926D-25/,
!    *     A(35) / .1864982941717294430718413161878666898946D-25/,
!    *     A(36) /-.2218095624207197204399716913626860379732D-26/,
!    *     A(37) / .1297781974947993668824414486330594165619D-27/,
!    *     A(38) / .1180697474966528406222745415509971518560D-29/,
!    *     A(39) /-.1124584349277088090293654674261439512119D-29/,
!    *     A(40) / .1277085175140866203990206677751124647749D-30/
!     DATA A(41) /-.7391451169615140823461289330108552823711D-32/,
!    *     A(42) / .1134750257554215760954165259469306393009D-34/,
!    *     A(43) / .4639134641058722029944804907952228463058D-34/,
!    *     A(44) /-.5347336818439198875077418196709893320905D-35/,
!    *     A(45) / .3207995923613352622861237279082794391090D-36/,
!    *     A(46) /-.4445829736550756882101590352124643637401D-38/,
!    *     A(47) /-.1311174518881988712901058494389922190237D-38/,
!    *     A(48) / .1647033352543813886818259327906394145400D-39/,
!    *     A(49) /-.1056233178503581218600561071538285049997D-40/

!                           -----------

!     C = A(1) - 1 IS ALSO FREQUENTLY NEEDED. C HAS THE VALUE ...

!     DATA C /-.4227843350984671393934879099175975689578D+00/

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: d, t, w, z
REAL (dp), PARAMETER :: a0 = .611609510448141581788D-08, a1  &
        = .624730830116465516210D-08, b1 = .203610414066806987300D+00, b2  &
        = .266205348428949217746D-01, b3 = .493944979382446875238D-03, b4  &
        = -.851419432440314906588D-05, b5 = -.643045481779353022248D-05, b6  &
        = .992641840672773722196D-06, b7 = -.607761895722825260739D-07, b8  &
        = .195755836614639731882D-09
REAL (dp), PARAMETER :: p0 = .6116095104481415817861D-08, p1  &
        = .6871674113067198736152D-08, p2 = .6820161668496170657, p3  &
        = .4686843322948848031080D-10, p4 = .1572833027710446286995D-11, p5  &
        = -.1249441572276366213222D-12, p6 = .4343529937408594255178D-14, q1  &
        = .3056961078365221025009D+00, q2 = .5464213086042296536016D-01, q3  &
        = .4956830093825887312, q4 = .2692369466186361192876D-03
REAL (dp), PARAMETER :: c = -.422784335098467139393487909917598D+00, c0  &
        = .577215664901532860606512090082402D+00, c1  &
        = -.655878071520253881077019515145390D+00, c2  &
        = -.420026350340952355290039348754298D-01, c3  &
        = .166538611382291489501700795102105D+00, c4  &
        = -.421977345555443367482083012891874D-01, c5  &
        = -.962197152787697356211492167234820D-02, c6  &
        = .721894324666309954239501034044657D-02, c7  &
        = -.116516759185906511211397108401839D-02, c8  &
        = -.215241674114950972815729963053648D-03, c9  &
        = .128050282388116186153198626328164D-03, c10  &
        = -.201348547807882386556893914210218D-04, c11  &
        = -.125049348214267065734535947383309D-05, c12  &
        = .113302723198169588237412962033074D-05, c13  &
        = -.205633841697760710345015413002057D-06
!----------------------------
t = x
d = x - 0.5_dp
IF (d > 0._dp) t = d - 0.5_dp
IF (t < 0.0_dp) THEN
  GO TO 30
ELSE IF (t > 0.0_dp) THEN
  GO TO 20
END IF

fn_val = 0._dp
RETURN
!------------

!                CASE WHEN 0 < T <= 0.5

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
20 w = ((((((p6*t + p5)*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) /   &
       ((((q4*t+q3)*t + q2)*t + q1)*t + 1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c0

IF (d <= 0._dp) THEN
  fn_val = x * z
  RETURN
END IF
fn_val = (t/x) * ((z-0.5_dp)-0.5_dp)
RETURN
!------------

!                CASE WHEN -0.5 <= T < 0

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
30 w = (a1*t + a0) / ((((((((b8*t + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t+1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c

IF (d <= 0._dp) THEN
  fn_val = x * ((z+0.5_dp)+0.5_dp)
  RETURN
END IF
fn_val = t * z / x
RETURN
END FUNCTION dgam1



FUNCTION dpdel(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF THE FUNCTION DEL(X) FOR  X >= 10  WHERE
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X)

!                         --------

!     THE SERIES FOR DPDEL ON THE INTERVAL 0.0 TO 1.0 DERIVED BY
!     A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY
!     OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(15) = (/ .833333333333333333333333333333D-01,  &
        -.277777777777777777777777752282D-04,  &
        .793650793650793650791732130419D-07,  &
        -.595238095238095232389839236182D-09,  &
        .841750841750832853294451671990D-11,  &
        -.191752691751854612334149171243D-12,  &
        .641025640510325475730918472625D-14,  &
        -.295506514125338232839867823991D-15,  &
        .179643716359402238723287696452D-16,  &
        -.139228964661627791231203060395D-17,  &
        .133802855014020915603275339093D-18,  &
        -.154246009867966094273710216533D-19,  &
        .197701992980957427278370133333D-20,  &
        -.234065664793997056856992426667D-21,  &
        .171348014966398575409015466667D-22 /)
REAL (dp) :: t, w
INTEGER   :: i, k
!-----------------------------------------------------------------------
t = (10._dp/x) ** 2
w = a(15)
DO i = 1, 14
  k = 15 - i
  w = t * w + a(k)
END DO
fn_val = w / x
RETURN
END FUNCTION dpdel

END MODULE Complex_Bessel
