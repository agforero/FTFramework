MODULE double_complex_erf

! Double precision complex error function
! Adapted from the Naval Surface Warfare Center Mathematics Library
! by Alan.Miller @ vic.cmis.csiro.au
! http://www.ozemail.com.au/~milleraj

! Latest revision - 18 September 2000

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION dnorm(x, y) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

fn_val = MAX(ABS(x), ABS(y))
RETURN
END FUNCTION dnorm



SUBROUTINE dcrec (x, y, u, v)
!-----------------------------------------------------------------------
!             COMPLEX RECIPROCAL U + I*V = 1/(X + I*Y)
!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: x, y
REAL (dp), INTENT(OUT)  :: u, v

REAL (dp) :: d, t

IF (ABS(x) > ABS(y)) GO TO 10
t = x/y
d = y + t*x
u = t/d
v = -1._dp/d
RETURN

10 t = y/x
d = x + t*y
u = 1._dp/d
v = -t/d

RETURN
END SUBROUTINE dcrec



SUBROUTINE cdivid (ar, ai, br, bi, cr, ci)
!-----------------------------------------------------------------------
!     REAL (dp) COMPLEX DIVISION C = A/B AVOIDING OVERFLOW
!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: ar, ai, br, bi
REAL (dp), INTENT(OUT)  :: cr, ci

REAL (dp) :: d, t, u, v

IF (ABS(br) <= ABS(bi)) GO TO 10
t = bi/br
d = br + t*bi
u = (ar + ai*t)/d
v = (ai - ar*t)/d
cr = u
ci = v
RETURN

10 IF (bi == 0._dp) GO TO 20
t = br/bi
d = bi + t*br
u = (ar*t + ai)/d
v = (ai*t - ar)/d
cr = u
ci = v
RETURN

!     DIVISION BY ZERO. C = INFINITY

20 cr = HUGE(1.0_dp)
ci = cr

RETURN
END SUBROUTINE cdivid



SUBROUTINE dcerf (mo, z, w)
!-----------------------------------------------------------------------

!               COMPUTATION OF THE COMPLEX ERROR FUNCTION

!                          -----------------

!                       W = ERF(Z)    IF MO = 0
!                       W = ERFC(Z)   OTHERWISE

!-----------------------------------------------------------------------

INTEGER, INTENT(IN)     :: mo
REAL (dp), INTENT(IN)   :: z(2)
REAL (dp), INTENT(OUT)  :: w(2)

REAL (dp) :: m, n, n2, n4, np1
REAL (dp) :: c2, d, d2, e, eps, r, sn, tol, x, y
REAL (dp) :: a0(2), an(2), b0(2), bn(2)
REAL (dp) :: g0(2), gn(2), h0(2), hn(2)
REAL (dp) :: qf(2), sm(2), sz(2), tm(2), ts(2), w0(2), wn(2)
!------------------------
!     C = 1/SQRT(PI)
!------------------------
REAL (dp), PARAMETER :: c = .56418958354775628694807945156077_dp
!------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp.

eps = EPSILON(1.0_dp)

!------------------------
x = z(1)
y = z(2)
sn = 1._dp
IF (x >= 0._dp) GO TO 10
x = -x
y = -y
sn = -1._dp

10 r = x*x + y*y
sz(1) = x*x - y*y
sz(2) = 2._dp*x*y

IF (r <= 1._dp) GO TO 20
IF (r >= 144._dp) GO TO 100
IF (ABS(y) > 31.8_dp*x) GO TO 50
IF (ABS(y) > 7.0_dp*x .AND. r < 64._dp) GO TO 50
IF (ABS(y) > 3.2_dp*x .AND. r < 49._dp) GO TO 50
IF (ABS(y) > 2.0_dp*x .AND. r < 36._dp) GO TO 50
IF (ABS(y) > 1.2_dp*x .AND. r < 25._dp) GO TO 50
IF (ABS(y) > 0.9_dp*x .AND. r < 16._dp) GO TO 50
IF (r >= 6.25_dp) GO TO 80
IF (ABS(y) > 0.6_dp*x) GO TO 50
IF (r >= 4.0_dp) GO TO 40

d = x - 2._dp
IF (d*d + y*y < 1._dp) GO TO 40
GO TO 50

!                          TAYLOR SERIES

20 c2 = c + c
tm(1) = c2*x
tm(2) = c2*y
sm(1) = tm(1)
sm(2) = tm(2)
tol = 2._dp*eps
m = 0._dp

21 m = m + 1._dp
d = m + m + 1._dp
ts(1) = tm(1)*sz(1) - tm(2)*sz(2)
ts(2) = tm(1)*sz(2) + tm(2)*sz(1)
tm(1) = -ts(1)/m
tm(2) = -ts(2)/m
ts(1) = tm(1)/d
ts(2) = tm(2)/d
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
IF (dnorm(ts(1),ts(2)) > tol*dnorm(sm(1),sm(2))) GO TO 21

IF (mo /= 0) GO TO 30
w(1) = sn*sm(1)
w(2) = sn*sm(2)
RETURN

30 IF (sn == 1._dp) GO TO 31
w(1) = 1._dp + sm(1)
w(2) = sm(2)
RETURN

31 w(1) = 0.5_dp + (0.5_dp - sm(1))
w(2) = -sm(2)
RETURN

!                  TAYLOR SERIES AROUND Z0 = 2

40 tm(1) = x
tm(2) = y
CALL erfcm2 (0, tm, w)
IF (mo /= 0) GO TO 41
w(1) = sn*(0.5_dp + (0.5_dp - w(1)))
w(2) = - sn*w(2)
RETURN
41 IF (sn > 0._dp) RETURN
w(1) = 2._dp - w(1)
w(2) = - w(2)
RETURN

!            PADE APPROXIMATION FOR THE TAYLOR SERIES
!                    FOR  (EXP(Z*Z)/Z)*ERF(Z)

50 d = 4._dp
IF (r > 16._dp) d = 16._dp
IF (r > 64._dp) d = 64._dp
d2 = d*d
CALL dcrec (sz(1), sz(2), w(1), w(2))
a0(1) = 1._dp
a0(2) = 0._dp
an(1) = (w(1) + 4._dp/15._dp)*d
an(2) = w(2)*d
b0(1) = 1._dp
b0(2) = 0._dp
bn(1) = (w(1) - 0.4_dp)*d
bn(2) = w(2)*d
CALL cdivid (an(1), an(2), bn(1), bn(2), wn(1), wn(2))
tol = 10._dp*eps
n4 = 0._dp

60 n4 = n4 + 4._dp
e = (n4 + 1._dp)*(n4 + 5._dp)
tm(1) = d*(w(1) - 2._dp/e)
tm(2) = d*w(2)
e = d2*(n4*(n4 + 2.0))/((n4 - 1.0)*(n4 + 3.0)*(n4 + 1.0)**2)

qf(1) = (tm(1)*an(1) - tm(2)*an(2)) + e*a0(1)
qf(2) = (tm(1)*an(2) + tm(2)*an(1)) + e*a0(2)
a0(1) = an(1)
a0(2) = an(2)
an(1) = qf(1)
an(2) = qf(2)
qf(1) = (tm(1)*bn(1) - tm(2)*bn(2)) + e*b0(1)
qf(2) = (tm(1)*bn(2) + tm(2)*bn(1)) + e*b0(2)
b0(1) = bn(1)
b0(2) = bn(2)
bn(1) = qf(1)
bn(2) = qf(2)

w0(1) = wn(1)
w0(2) = wn(2)
CALL cdivid (an(1), an(2), bn(1), bn(2), wn(1), wn(2))
IF (dnorm(wn(1) - w0(1), wn(2) - w0(2)) > tol*dnorm(wn(1), wn(2))) GO TO 60

c2 = c + c
sm(1) = c2*(x*wn(1) - y*wn(2))
sm(2) = c2*(x*wn(2) + y*wn(1))
e = EXP(-sz(1))
qf(1) = e*COS(-sz(2))
qf(2) = e*SIN(-sz(2))
tm(1) = qf(1)*sm(1) - qf(2)*sm(2)
tm(2) = qf(1)*sm(2) + qf(2)*sm(1)

w(1) = sn*tm(1)
w(2) = sn*tm(2)
IF (mo == 0) RETURN
w(1) = 1._dp - w(1)
w(2) = - w(2)
RETURN

!         PADE APPROXIMATION FOR THE ASYMPTOTIC EXPANSION
!                    FOR  Z*EXP(Z*Z)*ERFC(Z)

80 d = 4._dp*r
IF (r < 16._dp) d = 16._dp*r
d2 = d*d
tm(1) = sz(1) + sz(1)
tm(2) = sz(2) + sz(2)
g0(1) = 1._dp
g0(2) = 0._dp
gn(1) = (2._dp + tm(1))/d
gn(2) = tm(2)/d
h0(1) = 1._dp
h0(2) = 0._dp
tm(1) = 3._dp + tm(1)
hn(1) = tm(1)/d
hn(2) = tm(2)/d
CALL cdivid (gn(1), gn(2), hn(1), hn(2), wn(1), wn(2))
np1 = 1._dp
tol = 10._dp*eps

90 n = np1
np1 = n + 1._dp
n2 = n + n
e = (n2*(n2 + 1._dp))/d2
tm(1) = tm(1) + 4._dp
qf(1) = (tm(1)*gn(1) - tm(2)*gn(2))/d - e*g0(1)
qf(2) = (tm(1)*gn(2) + tm(2)*gn(1))/d - e*g0(2)
g0(1) = gn(1)
g0(2) = gn(2)
gn(1) = qf(1)
gn(2) = qf(2)
qf(1) = (tm(1)*hn(1) - tm(2)*hn(2))/d - e*h0(1)
qf(2) = (tm(1)*hn(2) + tm(2)*hn(1))/d - e*h0(2)
h0(1) = hn(1)
h0(2) = hn(2)
hn(1) = qf(1)
hn(2) = qf(2)

w0(1) = wn(1)
w0(2) = wn(2)
CALL cdivid (gn(1), gn(2), hn(1), hn(2), wn(1), wn(2))
IF (dnorm(wn(1) - w0(1), wn(2) - w0(2)) > tol*dnorm(wn(1), wn(2))) GO TO 90

tm(1) = x*hn(1) - y*hn(2)
tm(2) = x*hn(2) + y*hn(1)
CALL cdivid (c*gn(1), c*gn(2), tm(1), tm(2), sm(1), sm(2))
GO TO 130

!                      ASYMPTOTIC EXPANSION

100 CALL dcrec (x, y, tm(1), tm(2))
sm(1) = tm(1)
sm(2) = tm(2)
qf(1) = tm(1)*tm(1) - tm(2)*tm(2)
qf(2) = 2._dp*tm(1)*tm(2)
tol = 2._dp*eps
d = -0.5_dp
110 d = d + 1._dp
ts(1) = tm(1)*qf(1) - tm(2)*qf(2)
ts(2) = tm(1)*qf(2) + tm(2)*qf(1)
tm(1) = -d*ts(1)
tm(2) = -d*ts(2)
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (dnorm(tm(1), tm(2)) > tol*dnorm(sm(1), sm(2))) GO TO 110
sm(1) = c*sm(1)
sm(2) = c*sm(2)
IF (x < 1.d-2) GO TO 200

!                       TERMINATION

130 e = EXP(-sz(1))
qf(1) = e*COS(-sz(2))
qf(2) = e*SIN(-sz(2))
ts(1) = qf(1)*sm(1) - qf(2)*sm(2)
ts(2) = qf(1)*sm(2) + qf(2)*sm(1)
sm(1) = ts(1)
sm(2) = ts(2)

IF (mo /= 0) GO TO 140
w(1) = sn*(0.5_dp + (0.5_dp - sm(1)))
w(2) = - sn*sm(2)
RETURN

140 IF (sn == 1._dp) GO TO 141
w(1) = 2._dp - sm(1)
w(2) = -sm(2)
RETURN

141 w(1) = sm(1)
w(2) = sm(2)
RETURN

!               MODIFIED ASYMPTOTIC EXPANSION

200 e = EXP(-sz(1))
qf(1) = e*COS(-sz(2))
qf(2) = e*SIN(-sz(2))
w(1) = qf(1)*sm(1) - qf(2)*sm(2)
w(2) = qf(1)*sm(2) + qf(2)*sm(1)
IF (mo == 0) GO TO 210
w(1) = 1._dp + sn*w(1)
w(2) = sn*w(2)
RETURN

210 IF (sn < 0.0) RETURN
w(1) = - w(1)
w(2) = - w(2)

RETURN
END SUBROUTINE dcerf



SUBROUTINE erfcm2 (mo, z, w)
!-----------------------------------------------------------------------
!           CALCULATION OF ERFC(Z) USING THE TAYLOR SERIES
!                          AROUND Z0 = 2
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)     :: mo
REAL (dp), INTENT(IN)   :: z(2)
REAL (dp), INTENT(OUT)  :: w(2)

REAL (dp) :: eps, h(2), t(2), tol, x, y
INTEGER   :: j, n
!------------------------------
!     C = (2/SQRT(PI))*EXP(-4)
!     E = ERFC(2)
!------------------------------
REAL (dp), PARAMETER :: c = .020666985354092053857068941306585476_dp,  &
                        e = .0046777349810472658379307436327470714_dp
!------------------------------
REAL (dp) :: a(63) = (/ .20000000000000000000000000000000000D+01,  &
   .23333333333333333333333333333333333D+01,  .16666666666666666666666666666666667D+01,  &
   .63333333333333333333333333333333333D+00, -.22222222222222222222222222222222222D-01,  &
  -.16349206349206349206349206349206349D+00, -.76984126984126984126984126984126984D-01,  &
  -.24250440917107583774250440917107584D-02,  .12716049382716049382716049382716049D-01,  &
   .50208433541766875100208433541766875D-02, -.25305969750414194858639303083747528D-03,  &
  -.78593217482106370995259884148773038D-03, -.19118154038788959423880058800693721D-03,  &
   .46324144207742091339974937858535742D-04,  .33885549097189308829520469732109944D-04,  &
   .28637897646612243562134629672756034D-05, -.29071891082127275370004560446169188D-05,  &
  -.89674405786490646425523560263096103D-06,  .96069103941908684338469767911200105D-07,  &
   .99432863129093191401848891268744113D-07,  .97610310501460621303387795457283579D-08,  &
  -.65557500375673133822289344530892436D-08, -.18706782059105426900361744016236561D-08,  &
   .20329898993447386223176373714372370D-09,  .16941915827254374668448114614201210D-09,  &
   .10619149520827430973786114446699534D-10, -.10136148256511788733365237088810952D-10,  &
  -.21042890133669970575386166675721692D-11,  .37186985840699828780916522245407087D-12,  &
   .17921843632701679986488128324051002D-12, -.89823991804248069863542565948598397D-16,  &
  -.10533182313660970970232171410372199D-13, -.12340742690978398320850088252659714D-14,  &
   .44315624546581333350568244777175883D-15,  .11584041639989442481950487524296214D-15,  &
  -.10765703619385988116658460442219647D-16, -.70653158723054941879586082239984222D-17,  &
  -.18708903154917138727191793341667090D-18,  .32549879318817103966053527398133297D-18,  &
   .40654116689599228385911733319215613D-19, -.11250074516817311101947327325293424D-19,  &
  -.28923865378584966737386008432031980D-20,  .23653053641701517160704870522922706D-21,  &
   .14665384680061888088099002254334292D-21,  .26971039707314316218154193225264469D-23,  &
  -.58753834789274356433279671015522650D-23, -.59960357240498652932299485494869633D-24,  &
   .18586826578121663981412155416486531D-24,  .38364131854692721721867481914852428D-25,  &
  -.41342210492630142578080062451711039D-26, -.17646283105274988992381528904600860D-26,  &
   .19828685934364181151988692232131608D-28,  .65592252170840353572672782446212733D-28,  &
   .40626551379996340638338449938639730D-29, -.20097984104191034713653294173834095D-29,  &
  -.28104226475997460044096389060743131D-30,  .48705319298749358709127987806547949D-31,  &
   .12664655832830787747161769929972617D-31, -.75168312488894341862391776330113688D-33,  &
  -.45760473722605993842481669806804415D-33, -.56725491529575395930156379514718000D-35,  &
   .13932664042920082608489441616061541D-34,  .10452448992516358449586503951463322D-35 /)
!------------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS .GT. 1._dp

eps = EPSILON(1.0_dp)

!------------------------------
tol = eps*1.d+12
h(1) = 1._dp + (1._dp - z(1))
h(2) = - z(2)

x = 1._dp
y = 0._dp
w(1) = a(30)
w(2) = 0._dp
DO n = 31,63
  t(1) = x*h(1) - y*h(2)
  t(2) = x*h(2) + y*h(1)
  x = t(1)
  y = t(2)
  t(1) = a(n)*x
  t(2) = a(n)*y
  w(1) = w(1) + t(1)
  w(2) = w(2) + t(2)
  IF (dnorm(t(1), t(2)) <= tol*dnorm(w(1), w(2))) EXIT
END DO

DO j = 1,29
  n = 30 - j
  x = h(1)*w(1) - h(2)*w(2)
  w(2) = h(1)*w(2) + h(2)*w(1)
  w(1) = a(n) + x
END DO
x = h(1)*w(1) - h(2)*w(2)
w(2) = h(1)*w(2) + h(2)*w(1)
w(1) = 1._dp + x

x = c*(h(1)*w(1) - h(2)*w(2))
w(2) = c*(h(1)*w(2) + h(2)*w(1))
w(1) = e + x
IF (mo == 0) RETURN

!                     COMPUTE EXP(Z*Z)*ERFC(Z)

x = z(1)*z(1) - z(2)*z(2)
y = 2._dp*z(1)*z(2)
x = EXP(x)
t(1) = x*COS(y)
t(2) = x*SIN(y)
x = t(1)*w(1) - t(2)*w(2)
y = t(1)*w(2) + t(2)*w(1)
w(1) = x
w(2) = y

RETURN
END SUBROUTINE erfcm2

END MODULE double_complex_erf
