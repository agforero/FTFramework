SUBROUTINE cexpli(mo, z, w)

!-----------------------------------------------------------------------
!           EVALUATION OF THE COMPLEX EXPONENTIAL INTEGRAL
!-----------------------------------------------------------------------

! From the NSWC Mathematics Library
! Original code by Hershey, A.V. & Alfred Morris.

! Arguments:
! If mo = 0, then w = Ei(z),
!       otherwise w = exp(-z).Ei(z)

! Accuracy:
! If mo = 0, Re(w) and Im(w) are accurate to within 2 units in the 12th
! significant digit, except when z is near a zero of either the real or
! imaginary part of Ei(z).

IMPLICIT NONE
INTEGER, PARAMETER   :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)       :: mo
COMPLEX (dp), INTENT(IN)  :: z
COMPLEX (dp), INTENT(OUT) :: w

! Local variables
REAL (dp), PARAMETER :: pi = 3.14159265358979_dp, euler = .577215664901533_dp
REAL (dp), PARAMETER :: cd(18) = (/ 0.0_dp, .311105957086528E-01_dp,  &
     .103661260539112E+00_dp, .216532335244554E+00_dp, .369931427960192E+00_dp,  &
     .566766259990589E+00_dp, .814042066324748E+00_dp, .112384247540813E+01_dp,  &
     .151400478148512E+01_dp, .200886795032284E+01_dp, .264052411823592E+01_dp,  &
     .345098449933392E+01_dp, .449583360763202E+01_dp, .585058263409822E+01_dp,  &
     .762273501463380E+01_dp, .997814501584578E+01_dp, .132122064896408E+02_dp,  &
     .180322948376021E+02_dp /),   &
                        ce(18) = (/ .850156516121093E-02_dp,  &
     .505037465849058E-01_dp, .836817368956407E-01_dp, .107047582417607E+00_dp,  &
     .120424719029462E+00_dp, .125096631582229E+00_dp, .122314435224685E+00_dp,  &
     .112621417553907E+00_dp, .963419407392582E-01_dp, .747398422757511E-01_dp,  &
     .508596135953441E-01_dp, .290822706773628E-01_dp, .132201640530101E-01_dp,  &
     .443802939829067E-02_dp, .992612478987576E-03_dp, .126579795112011E-03_dp,  &
     .702150908253350E-05_dp, .910281532564632E-07_dp /)
REAL (dp)    :: qf(2), sm(2), tm(2), ts(2), g0(2), gn(2), h0(2), hn(2), wn(2),  &
                c, cy, d, e, eps, n, np1, qm, r, ss, sy, tol, u, x, y
INTEGER      :: i
LOGICAL      :: ind

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.D0 + EPS > 1.D0.

eps = EPSILON(1.0_dp)

!-------------------------

x = REAL(z, KIND=dp)
y = AIMAG(z)
r = dcpabs(x,y)
eps = MAX(eps, 1.e-15_dp)

IF (r > 1.0) THEN
  IF (r >= 40.0) GO TO 40
  IF (r >= 4.0) THEN
    IF (x <= 0.0 .OR. ABS(y) > 8.0) GO TO 40
    IF (r < 10.0 .AND. ABS(y) > 1.8*x) GO TO 40
  ELSE
    IF (x < 0.09*y*y) GO TO 20
    IF (r > 3.6 .AND. ABS(y) > 1.8*x) GO TO 40
  END IF
END IF

!                        TAYLOR SERIES

sm(1) = 0.0
sm(2) = 0.0
tm(1) = x
tm(2) = y
n = 1.0
10 n = n + 1.0
ts(1) = tm(1) * x - tm(2) * y
ts(2) = tm(1) * y + tm(2) * x
tm(1) = ts(1) / n
tm(2) = ts(2) / n
ts(1) = tm(1) / n
ts(2) = tm(2) / n
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
IF (anorm(ts(1), ts(2)) > eps*anorm(sm(1), sm(2))) GO TO 10
sm(1) = x + sm(1)
sm(2) = y + sm(2)

sm(1) = (euler + LOG(r)) + sm(1)
sm(2) = ATAN2(-y, -x) + sm(2)
GO TO 70

!                      RATIONAL EXPANSION

20 sm(1) = 0.0
sm(2) = 0.0
DO i = 1, 18
  ts(1) = x - cd(i)
  ts(2) = y
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  sm(1) = sm(1) + ce(i) * ts(1) / ss
  sm(2) = sm(2) - ce(i) * ts(2) / ss
END DO
GO TO 60

!         PADE APPROXIMATION FOR THE ASYMPTOTIC EXPANSION
!                       FOR EXP(-Z)*EI(Z)

40 x = -x
y = -y
d = 4.0 * r
IF (r < 10.0) d = 32.0
g0(1) = 1.0
g0(2) = 0.0
gn(1) = (1.0+x) / d
gn(2) = y / d
h0(1) = 1.0
h0(2) = 0.0
u = x + 2.0
hn(1) = u / d
hn(2) = gn(2)
w = CMPLX(1.0+x, y, KIND=dp) / CMPLX(u, y, KIND=dp)
wn(1) = REAL(w, KIND=dp)
wn(2) = AIMAG(w)
np1 = 1.0
tol = 4.0 * eps

50 n = np1
np1 = n + 1.0
e = (n*np1) / d
u = u + 2.0
tm(1) = ((u*gn(1)-y*gn(2))-e*g0(1)) / d
tm(2) = ((u*gn(2)+y*gn(1))-e*g0(2)) / d
g0(1) = gn(1)
g0(2) = gn(2)
gn(1) = tm(1)
gn(2) = tm(2)
tm(1) = ((u*hn(1)-y*hn(2))-e*h0(1)) / d
tm(2) = ((u*hn(2)+y*hn(1))-e*h0(2)) / d
h0(1) = hn(1)
h0(2) = hn(2)
hn(1) = tm(1)
hn(2) = tm(2)

tm(1) = wn(1)
tm(2) = wn(2)
w = CMPLX(gn(1), gn(2), KIND=dp) / CMPLX(hn(1), hn(2), KIND=dp)
wn(1) = REAL(w, KIND=dp)
wn(2) = AIMAG(w)
IF (anorm(tm(1)-wn(1),tm(2)-wn(2)) > tol*anorm(wn(1),wn(2)))  &
GO TO 50

x = REAL(z, KIND=dp)
y = AIMAG(z)
w = w / z
sm(1) = REAL(w, KIND=dp)
sm(2) = AIMAG(w)

!                         TERMINATION

60 ind = x <= 0.0 .OR. ABS(y) > 1.e-2
IF (ind .AND. mo /= 0) GO TO 90
c = pi
IF (y > 0.0) c = -pi
qm = EXP(x)
cy = COS(y)
sy = SIN(y)
qf(1) = qm * cy
qf(2) = qm * sy
IF (mo == 0) GO TO 80

r = c / qm
sm(1) = sm(1) + r * sy
sm(2) = sm(2) + r * cy
GO TO 90

70 IF (mo == 0) GO TO 90
ind = .true.
qm = EXP(-x)
qf(1) = qm * COS(-y)
qf(2) = qm * SIN(-y)

80 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
IF (.NOT.ind) sm(2) = sm(2) + c

90 w = CMPLX(sm(1), sm(2), KIND=dp)
RETURN


CONTAINS


FUNCTION dcpabs(x, y) RESULT(fn_val)
!     --------------------------------------
!     EVALUATION OF SQRT(X*X + Y*Y)
!     --------------------------------------
REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

! Local variable
REAL (dp) :: a

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
END FUNCTION dcpabs



FUNCTION anorm(x, y) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

fn_val = MAX(ABS(x), ABS(y))
RETURN
END FUNCTION anorm


END SUBROUTINE cexpli
