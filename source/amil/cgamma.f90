SUBROUTINE cgamma(mo, z, w)
!-----------------------------------------------------------------------

!        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

!                        ---------------

!     MO IS AN INTEGER, Z A COMPLEX ARGUMENT, AND W A COMPLEX VARIABLE.

!                 W = GAMMA(Z)       IF MO = 0
!                 W = LN(GAMMA(Z))   OTHERWISE

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA

!     This version, in a subset of Fortran 90, prepared by
!     Alan.Miller @ vic.cmis.csiro.au
!     http://www.ozemail.com.au/~milleraj

!     This version is accurate to within 5 in the 14th significant
!     decimal digit.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

INTEGER, INTENT(IN)       :: mo
COMPLEX (dp), INTENT(IN)  :: z
COMPLEX (dp), INTENT(OUT) :: w

! Local variables
COMPLEX (dp) :: eta, eta2, sum
REAL (dp), PARAMETER :: c0(12) = (/ .833333333333333E-01_dp,  &
          -.277777777777778E-02_dp, .793650793650794E-03_dp,  &
          -.595238095238095E-03_dp, .841750841750842E-03_dp,  &
          -.191752691752692E-02_dp, .641025641025641E-02_dp,  &
          -.295506535947712E-01_dp, .179644372368831_dp,      &
          -1.39243221690590_dp,     13.4028640441684_dp,      &
          -156.848284626002_dp /), pi = 3.14159265358979_dp,  &
    pi2  = 6.28318530717959_dp, alpi = 1.14472988584940_dp,  &
    hl2p = .918938533204673_dp, half = 0.5_dp
REAL (dp)  :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
              s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
INTEGER    :: j, k, l, m, max, n, nm1
!---------------------------
!     ALPI = LOG(PI)
!     HL2P = 0.5 * LOG(2*PI)
!---------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            SUCH THAT 1.0 + EPS > 1.0.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0_dp)

!---------------------------
x = REAL(z, KIND=dp)
y = AIMAG(z)
IF (x < 0.0_dp) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi * y
  et = EXP(t)
  e2t = et * et

!     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

  a1 = half * (1.0_dp + e2t)
  t2 = t + t
  IF (t2 >= -0.15_dp) THEN
    a2 = -half * rexp(t2)
  ELSE
    a2 = half * (half + (half - e2t))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X)

  IF (ABS(x) >= MIN(REAL(MAX), 1.0_dp/eps)) GO TO 70
  k = ABS(x)
  u = x + k
  k = MOD(k,2)
  IF (u <= -half) THEN
    u = half + (half + u)
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
  IF (a == 0.0_dp) GO TO 70
  IF (mo == 0) THEN

    h1 = a1 / a
    h2 = -a2 / a
    c = pi * et
    h1 = c * h1
    h2 = c * h2
  ELSE

    h1 = (alpi+t) - half * LOG(a)
    h2 = -ATAN2(a2,a1)
  END IF
  IF (AIMAG(z) >= 0.0_dp) THEN
    x = 1.0_dp - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0_dp - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
w1 = 0.0_dp
w2 = 0.0_dp
n = 0
t = x
y2 = y * y
a = t * t + y2
cut = 36.0_dp
IF (eps > 1.e-8_dp) cut = 16.0_dp
IF (a < cut) THEN
  IF (a == 0.0_dp) GO TO 70
  10 n = n + 1
  t = t + 1.0_dp
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
    w1 = half * LOG(s)
    w2 = ATAN2(s2,s1)
  END IF
END IF

!     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

t1 = half * LOG(a) - 1.0_dp
t2 = ATAN2(y,t)
u = x - half
v1 = (u*t1-half) - y * t2
v2 = u * t2 + y * t1

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

eta = CMPLX(t/a, -y/a, KIND=dp)
eta2 = eta * eta
m = 12
IF (a >= 289.0_dp) m = 6
IF (eps > 1.e-8) m = m / 2
sum = CMPLX(c0(m), 0.0_dp, KIND=dp)
l = m
DO j = 2, m
  l = l - 1
  sum = CMPLX(c0(l), 0.0_dp, KIND=dp) + sum * eta2
END DO
sum = sum * eta
a1 = REAL(sum, KIND=dp)
a2 = AIMAG(sum)
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
w1 = (((a1 + hl2p) - w1) + v1) - n
w2 = (a2 - w2) + v2
IF (REAL(z, KIND=dp) < 0.0_dp) GO TO 50
IF (mo == 0) THEN

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

  a = EXP(w1)
  w1 = a * COS(w2)
  w2 = a * SIN(w2)
  IF (n == 0) GO TO 60
  c = (s1*w1 + s2*w2) / s
  d = (s1*w2 - s2*w1) / s
  w1 = c
  w2 = d
  GO TO 60
END IF

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
!     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

40 IF (w2 <= pi) THEN
  k = half - w2 / pi2
  w2 = w2 + pi2 * k
  GO TO 60
END IF
k = w2 / pi2 - half
w2 = w2 - pi2 * REAL(k+1)
IF (w2 <= -pi) w2 = pi
GO TO 60

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

50 IF (mo /= 0) THEN
  w1 = h1 - w1
  w2 = h2 - w2
  GO TO 40
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

60 w = CMPLX(w1, w2, KIND=dp)
RETURN
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
70 w = (0.0_dp, 0.0_dp)
RETURN

CONTAINS


FUNCTION rexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER           :: p1 =  .914041914819518E-09_dp,  &
    p2 = .238082361044469E-01_dp, q1 = -.499999999085958_dp,      &
    q2 = .107141568980644_dp,     q3 = -.119041179760821E-01_dp,  &
    q4 = .595130811860248E-03_dp
REAL (dp) :: e
!-----------------------
IF (ABS(x) <= 0.15_dp) THEN
  fn_val = x * (((p2*x + p1)*x + 1.0_dp) /  &
           ((((q4*x + q3)*x + q2)*x + q1)*x + 1.0_dp))
  RETURN
END IF

IF (x >= 0.0_dp) THEN
  e = EXP(x)
  fn_val = e * (half + (half - 1.0_dp/e))
  RETURN
END IF
IF (x >= -37.0_dp) THEN
  fn_val = (EXP(x) - half) - half
  RETURN
END IF
fn_val = -1.0_dp
RETURN
END FUNCTION rexp

END SUBROUTINE cgamma
