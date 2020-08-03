MODULE LogBeta
! Log of the beta function
! Includes log of the gamma function

IMPLICIT NONE


CONTAINS


FUNCTION betaln(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: e = .918938533204673_dp
REAL (dp)             :: a, b, c, h, u, v, w, z
INTEGER               :: i, n
REAL (dp), PARAMETER  :: half = 0.5_dp, one = 1.0_dp, two = 2.0_dp,  &
                         eight = 8.0_dp
!--------------------------
a = MIN(a0, b0)
b = MAX(a0, b0)
IF (a < eight) THEN
  IF (a < one) THEN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A < 1
!-----------------------------------------------------------------------
    IF (b < eight) THEN
      fn_val = gamln(a) + (gamln(b) - gamln(a+b))
      RETURN
    END IF
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 <= A < 8
!-----------------------------------------------------------------------
  IF (a <= two) THEN
    IF (b <= two) THEN
      fn_val = gamln(a) + gamln(b) - gsumln(a,b)
      RETURN
    END IF
    w = 0.0
    IF (b < eight) GO TO 20
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                REDUCTION OF A WHEN B <= 1000

  IF (b > 1000.0_dp) GO TO 40
  n = a - one
  w = one
  DO i = 1, n
    a = a - one
    h = a / b
    w = w * (h/(one + h))
  END DO
  w = LOG(w)
  IF (b >= eight) THEN
    fn_val = w + gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                 REDUCTION OF B WHEN B < 8

  20 n = b - one
  z = one
  DO i = 1, n
    b = b - one
    z = z * (b/(a+b))
  END DO
  fn_val = w + LOG(z) + (gamln(a) + (gamln(b) - gsumln(a,b)))
  RETURN

!                REDUCTION OF A WHEN B > 1000

  40 n = a - one
  w = one
  DO i = 1, n
    a = a - one
    w = w * (a/(one+a/b))
  END DO
  fn_val = (LOG(w) - n*LOG(b)) + (gamln(a) + algdiv(a,b))
  RETURN
END IF
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A >= 8
!-----------------------------------------------------------------------
w = bcorr(a,b)
h = a / b
c = h / (one + h)
u = -(a-half) * LOG(c)
v = b * alnrel(h)
IF (u > v) THEN
  fn_val = (((-half*LOG(b) + e) + w) - v) - u
  RETURN
END IF
fn_val = (((-half*LOG(b) + e) + w) - u) - v
RETURN
END FUNCTION betaln



FUNCTION gamln(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.418938533204673_dp, c0 = 0.833333333333333D-01,   &
              c1 = -0.277777777760991D-02, c2 = 0.793650666825390D-03,   &
              c3 = -0.595202931351870D-03, c4 = 0.837308034031215D-03,   &
              c5 = -0.165322962780713D-02
REAL (dp)  :: t, w
INTEGER    :: i, n
!-----------------------------------------------------------------------
IF (a <= 0.8_dp) THEN
  fn_val = gamln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.25_dp) THEN
  t = (a-0.5_dp) - 0.5_dp
  fn_val = gamln1(t)
  RETURN
END IF

IF (a < 10.0_dp) THEN
  n = a - 1.25_dp
  t = a
  w = 1.0_dp
  DO i = 1, n
    t = t - 1.0_dp
    w = t * w
  END DO
  fn_val = gamln1(t-1.0) + LOG(w)
  RETURN
END IF

t = (1.0/a) ** 2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a
fn_val = (d+w) + (a-0.5) * (LOG(a)-1.0)
RETURN
END FUNCTION gamln



FUNCTION algdiv(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
                    c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
                    c4 = .837308034031215D-03, c5 = -.165322962780713D-02
REAL (dp) :: c, d, h, s3, s5, s7, s9, s11, t, u, v, w, x, x2
!------------------------
IF (a > b) THEN
  h = b / a
  c = 1.0 / (1.0_dp + h)
  x = h / (1.0_dp + h)
  d = a + (b - 0.5_dp)
ELSE
  h = a / b
  c = h / (1.0_dp + h)
  x = 1.0 / (1.0_dp + h)
  d = b + (a - 0.5_dp)
END IF

!                SET SN = (1 - X**N)/(1 - X)

x2 = x * x
s3 = 1.0 + (x + x2)
s5 = 1.0 + (x + x2*s3)
s7 = 1.0 + (x + x2*s5)
s9 = 1.0 + (x + x2*s7)
s11 = 1.0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0_dp/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                    COMBINE THE RESULTS

u = d * alnrel(a/b)
v = a * (LOG(b) - 1.0_dp)
IF (u > v) THEN
  fn_val = (w-v) - u
  RETURN
END IF
fn_val = (w-u) - v
RETURN
END FUNCTION algdiv



FUNCTION gsumln(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 <= A <= 2  AND  1 <= B <= 2
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x

x = a + b - 2.0_dp
IF (x <= 0.25_dp) THEN
  fn_val = gamln1(1.0_dp + x)
  RETURN
END IF
IF (x <= 1.25_dp) THEN
  fn_val = gamln1(x) + alnrel(x)
  RETURN
END IF
fn_val = gamln1(x - 1.0_dp) + LOG(x*(1.0_dp + x))
RETURN
END FUNCTION gsumln



FUNCTION bcorr(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.

!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
                   c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
                   c4 = .837308034031215D-03, c5 = -.165322962780713D-02
REAL (dp) :: a, b, c, h, s3, s5, s7, s9, s11, t, w, x, x2
!------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)

h = a / b
c = h / (1.0_dp + h)
x = 1.0 / (1.0_dp + h)
x2 = x * x

!                SET SN = (1 - X**N)/(1 - X)

s3 = 1.0 + (x + x2)
s5 = 1.0 + (x + x2*s3)
s7 = 1.0 + (x + x2*s5)
s9 = 1.0 + (x + x2*s7)
s11 = 1.0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0_dp/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                   COMPUTE  DEL(A) + W

t = (1.0_dp/a) ** 2
fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a + w
RETURN
END FUNCTION bcorr



FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: p1 = -.129418923021993D+01, p2 = .405303492862024D+00,  &
        p3 = -.178874546012214D-01, q1 = -.162752256355323D+01,   &
        q2 =  .747811014037616D+00, q3 = -.845104217945565D-01,  &
        zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp, two = 2.0_dp
REAL (dp) :: t, t2, w, x
!--------------------------
IF (ABS(a) <= 0.375_dp) THEN
  t = a / (a + two)
  t2 = t * t
  w = (((p3*t2 + p2)*t2 + p1)*t2 + one) / (((q3*t2 + q2)*t2 + q1)*t2 + one)
  fn_val = two * t * w
  RETURN
END IF

x = one + a
IF (a < zero) x = (a + half) + half
fn_val = LOG(x)
RETURN
END FUNCTION alnrel



FUNCTION gamln1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: p0 = .577215664901533D+00, p1 = .844203922187225D+00,  &
        p2 = -.168860593646662D+00, p3 = -.780427615533591D+00,   &
        p4 = -.402055799310489D+00, p5 = -.673562214325671D-01,   &
        p6 = -.271935708322958D-02, q1 = .288743195473681D+01,   &
        q2  = .312755088914843D+01, q3 = .156875193295039D+01,   &
        q4  = .361951990101499D+00, q5 = .325038868253937D-01,   &
        q6  = .667465618796164D-03, r0 = .422784335098467D+00,   &
        r1  = .848044614534529D+00, r2 = .565221050691933D+00,   &
        r3  = .156513060486551D+00, r4 = .170502484022650D-01,   &
        r5  = .497958207639485D-03, s1 = .124313399877507D+01,   &
        s2  = .548042109832463D+00, s3 = .101552187439830D+00,   &
        s4  = .713309612391000D-02, s5 = .116165475989616D-03
REAL (dp)  :: w, x
!----------------------
IF (a < 0.6_dp) THEN
  w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0) /   &
      ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0)
  fn_val = -a * w
  RETURN
END IF

x = (a - 0.5_dp) - 0.5_dp
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0) / (((((s5*x + s4)*x +  &
    s3)*x + s2)*x + s1)*x + 1.0_dp)
fn_val = x * w
RETURN
END FUNCTION gamln1

END MODULE LogBeta



FUNCTION betain(x, p, q, beta) RESULT(fn_val)

!   Algorithm AS 63  Appl. Statist. (1973), vol.22, no.3

!   Computes incomplete beta function ratio for arguments
!   x between zero and one, p and q positive.
!   Log of complete beta function, beta, is assumed to be known

! ELF90-compatible version by Alan Miller
! N.B. Argument IFAULT has been removed

! Latest revision - 5 July 2003

USE LogBeta
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: x, p, q, beta
REAL (dp)             :: fn_val

! Local variables
LOGICAL               :: indx
INTEGER               :: ns
REAL (dp)             :: psq, cx, xx, pp, qq, term, ai, rx, temp

!   Define accuracy and initialise

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, acu = 1.0E-14_dp

fn_val = x

!   Test for admissibility of arguments

IF(p <= zero .OR. q <= zero) THEN
  WRITE(*, *) 'AS63: Either p or q <= 0'
  RETURN
END IF
IF(x < zero .OR. x > one) THEN
  WRITE(*, *) 'AS63: Argument x outside range (0, 1)'
  RETURN
END IF
IF(x == zero .OR. x == one) RETURN

!   Change tail if necessary and determine s

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

!     Use Soper's reduction formulae.

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

!     Calculate result

5 fn_val = fn_val*EXP(pp*LOG(xx) + (qq-one)*LOG(cx) - beta)/pp
IF(indx) fn_val = one - fn_val

RETURN
END FUNCTION betain
