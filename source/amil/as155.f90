MODULE Nonc_chisq
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL (dp), SAVE  :: aintl, ersm
REAL, SAVE       :: pi, aln28, sigsq, almax, almin, amean, c
INTEGER, SAVE    :: icount, ir, lim
LOGICAL, SAVE    :: ndtsrt, fail


CONTAINS


FUNCTION iround(x) RESULT(ival)
IMPLICIT NONE

REAL, INTENT(IN)  :: x
INTEGER           :: ival

REAL, PARAMETER  :: half = 0.5

ival = INT(x + SIGN(half,x))
RETURN
END FUNCTION iround



FUNCTION qf(alb, anc, n, irr, sigma, cc, lim1, acc, ith, trace, ifault)  &
   RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-01-06  Time: 08:37:41

!     ALGORITHM AS 155  APPL. STATIST. (1980) VOL.29, NO.3

!     Distribution of a linear combination of non-central chi-squared
!     random variables.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: irr
REAL, INTENT(IN)      :: alb(irr)
REAL, INTENT(IN)      :: anc(irr)
INTEGER, INTENT(IN)   :: n(irr)
REAL, INTENT(IN)      :: sigma
REAL, INTENT(IN)      :: cc
INTEGER, INTENT(IN)   :: lim1
REAL, INTENT(IN)      :: acc
INTEGER, INTENT(OUT)  :: ith(irr)
REAL, INTENT(OUT)     :: trace(7)
INTEGER, INTENT(OUT)  :: ifault
REAL                  :: fn_val

INTEGER  :: j, nj, nt, ntm
REAL     :: acc1, almx, xlim, xnt, xntm
REAL     :: utx, tausq, sd, aintv, aintv1, x, up, un, d1, d2, alj, ancj
REAL (dp) :: aintl, ersm
REAL     :: pi, aln28, sigsq, almax, almin, amean, c
INTEGER  :: icount, ir, lim
LOGICAL  :: ndtsrt, fail
REAL, PARAMETER  :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, four = 4.0, &
    sixtn = 16.0, fourp5 = 4.5, pt07 = 0.07, pt2 = 0.2, quart = 0.25,  &
    ten = 10.0, pt33 = 0.33, pt67 = 0.67, pt75 = 0.75, onept5 = 1.5,  &
    three = 3.0, onept1 = 1.1

!     Setting constants in COMMON.   ALN28 = ln(2) / 8.

pi = 3.14159265358979
aln28 = 0.0866

c = cc
ir = irr
lim = lim1
trace(1:7) = zero
ifault = 0
icount = 0
aintl = zero
ersm = zero
fn_val = -one
acc1 = acc
ndtsrt = .true.
fail = .false.

!     Find mean, sd, max & min of ALB.
!     Check that parameter values are valid.

xlim = lim
sigsq = sigma**2
sd = sigsq
almax = zero
almin = zero
amean = zero
j = 1

20 IF (.NOT.(j <= ir)) GO TO 60
nj = n(j)
alj = alb(j)
ancj = anc(j)
IF (.NOT.(nj < 0 .OR. ancj < zero)) GO TO 30
ifault = 3
GO TO 260

30 sd = sd + alj**2*(2*nj + four*ancj)
amean = amean + alj*(nj + ancj)
IF (.NOT.(almax < alj)) GO TO 40
almax = alj
GO TO 50

40 IF (.NOT.(almin > alj)) GO TO 50
almin = alj
50 j = j + 1
GO TO 20

60 IF (.NOT.(sd == zero)) GO TO 80
IF (.NOT.(c > zero)) GO TO 70
fn_val = one
GO TO 260

70 fn_val = zero
GO TO 260
80 IF (.NOT.(almin == zero .AND. almax == zero .AND. sigma == zero)) GO TO 90
ifault = 3
GO TO 260

90 sd = SQRT(sd)
IF (.NOT.(almax < -almin)) GO TO 100
almx = -almin
GO TO 110
100 almx = almax

!     Starting values for FINDU * CTFF.

110 utx = sixtn/sd
up = fourp5/sd
un = -up

!     Truncation point with no convergence factor.

CALL findu (n, alb, anc, utx, half*acc1)

!     Does convergence factor help ?

IF (.NOT.(c /= zero .AND. almx > pt07*sd)) GO TO 130
tausq = quart*acc1/cfe(n, alb, anc, ith, c)
IF (.NOT.(fail)) GO TO 120
fail = .false.
GO TO 130

120 IF (.NOT.(truncn(n, alb, anc, utx, tausq) < pt2*acc1)) GO TO 130
sigsq = sigsq + tausq
CALL findu (n, alb, anc, utx, quart*acc1)
trace(6) = SQRT(tausq)
130 trace(5) = utx
acc1 = half*acc1

!     Find 'range' of distribution, quit if outside of this.

140 d1 = ctff(n, alb, anc, acc1, up) - c
IF (.NOT.(d1 < zero)) GO TO 150
fn_val = one
GO TO 260
150 d2 = c - ctff(n, alb, anc, acc1, un)
IF (.NOT.(d2 < zero)) GO TO 160
fn_val = zero
GO TO 260

!     Find integration interval.

160 IF (.NOT.(d1 > d2)) GO TO 170
aintv = d1
GO TO 180
170 aintv = d2
180 aintv = two*pi/aintv

!     Calculate number of terms required for main & auxiliary
!     integrations.

xnt = utx/aintv
xntm = three/SQRT(acc1)
IF (.NOT.(xnt > xntm*onept5)) GO TO 220
IF (.NOT.(xntm > xlim)) GO TO 190
ifault = 1
GO TO 260

!     Parameters for auxiliary integration.

190 ntm = iround(xntm)
aintv1 = utx/xntm
x = two*pi/aintv1
IF (.NOT.(x <= ABS(c))) GO TO 200
GO TO 220

!     Calculate convergence factor.

200 tausq = cfe(n, alb, anc, ith, c - x) + cfe(n, alb, anc, ith, c + x)
tausq = pt33*acc1/(onept1*tausq)
IF (.NOT.(fail)) GO TO 210
GO TO 220
210 acc1 = pt67*acc1

!     Auxiliary integration.

CALL integr (n, alb, anc, ntm, aintv1, tausq, .false.)
xlim = xlim - xntm
sigsq = sigsq + tausq
trace(3) = trace(3) + 1
trace(2) = trace(2) + ntm + 1

!     Find truncation point with new convergence factor.

CALL findu (n, alb, anc, utx, quart*acc1)
acc1 = pt75*acc1
GO TO 140

!     Main integration.

220 trace(4) = aintv
IF (.NOT.(xnt > xlim)) GO TO 230
ifault = 1
GO TO 260

230 nt = iround(xnt)
CALL integr (n, alb, anc, nt, aintv, zero, .true.)
trace(3) = trace(3) + 1
trace(2) = trace(2) + nt + 1
fn_val = half - aintl
trace(1) = ersm
up = ersm

!     Test whether round-off error could be significant.
!     Allow for radix 8 or 16 machines.

x = up + acc/ten
j = 1
240 IF (.NOT.(j <= 8)) GO TO 260
IF (.NOT.(j*x == j*up)) GO TO 250
ifault = 2
250 j = j*2
GO TO 240

260 trace(7) = icount
RETURN
END FUNCTION qf



SUBROUTINE countr()

!     Count number of calls to ERRBD, TRUNCN & CFE.

IMPLICIT NONE
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim

icount = icount + 1
IF (.NOT.(icount > lim)) GO TO 20
WRITE (6, 10)
10 FORMAT (' qf: cannot locate integration parameters'/)
STOP
20 RETURN
END SUBROUTINE countr



FUNCTION alog1(x, first) RESULT(fn_val)

!     If FIRST then return ln(1 + x) else ln(1 + x) - x.

IMPLICIT NONE
REAL, INTENT(IN)     :: x
LOGICAL, INTENT(IN)  :: first
REAL                 :: fn_val

REAL    :: s, s1, term, y, ak
REAL, PARAMETER  :: pt1 = 0.1, one = 1.0, two = 2.0, three = 3.0

IF (.NOT.(ABS(x) > pt1)) GO TO 20
IF (.NOT.(first)) GO TO 10
fn_val = LOG(one + x)
GO TO 70

10 fn_val = LOG(one + x) - x
GO TO 70

20 y = x/(two + x)
term = two*y**3
ak = three
IF (.NOT.(first)) GO TO 30
s = two
GO TO 40

30 s = -x
40 s = s*y
y = y**2
s1 = s + term/ak

50 IF (.NOT.(s1 /= s)) GO TO 60
ak = ak + two
term = term*y
s = s1
s1 = s + term/ak
GO TO 50

60 fn_val = s
70 RETURN
END FUNCTION alog1



FUNCTION exp1(x) RESULT(fn_val)

IMPLICIT NONE
REAL, INTENT(IN)  :: x
REAL              :: fn_val

REAL, PARAMETER  :: zero = 0.0, neg50 = -50.0

IF (.NOT.(x < neg50)) GO TO 10
fn_val = zero
GO TO 20
10 fn_val = EXP(x)
20 RETURN
END FUNCTION exp1



SUBROUTINE order (alb, ith)

!     Find order of absolute values of ALB.

IMPLICIT NONE
REAL, INTENT(IN)         :: alb(*)
INTEGER, INTENT(IN OUT)  :: ith(*)

INTEGER  :: j, k, k1, ithk
REAL     :: alj
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim

j = 1
10 IF (.NOT.(j <= ir)) GO TO 60
alj = ABS(alb(j))
k = j - 1
20 IF (.NOT.(k > 0)) GO TO 40
ithk = ith(k)
k1 = k + 1
IF (.NOT.(alj > ABS(alb(ithk)))) GO TO 50
ith(k1) = ithk
k = k - 1
GO TO 20

40 k = 0
k1 = 1
50 ith(k1) = j
j = j + 1
GO TO 10

60 ndtsrt = .false.
RETURN
END SUBROUTINE order



FUNCTION errbd(n, alb, anc, uu, cx) RESULT(fn_val)

!     Find bound on tail probability using mgf.
!     Cut-off point returned to CX.

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n(*)
REAL, INTENT(IN)     :: alb(*)
REAL, INTENT(IN)     :: anc(*)
REAL, INTENT(IN)     :: uu
REAL, INTENT(OUT)    :: cx
REAL                 :: fn_val

REAL     :: u
REAL     :: sum1, alj, ancj, x, y, const
INTEGER  :: j, nj
! COMMON /qfcom/ aintl, ersm, pi, aln28, sigsq, almax, almin, amean, c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: half = 0.5, one = 1.0, two = 2.0

CALL countr()
u = uu
const = u*sigsq
sum1 = u*const
u = two*u
j = ir
10 IF (.NOT.(j > 0)) GO TO 20
nj = n(j)
alj = alb(j)
ancj = anc(j)
x = u*alj
y = one - x
const = const + alj*(ancj/y + nj)/y
sum1 = sum1 + ancj*(x/y)**2
sum1 = sum1 + nj*(x**2/y + alog1(-x,.false.))
j = j - 1
GO TO 10
20 fn_val = exp1(-half*sum1)
cx = const
RETURN
END FUNCTION errbd



FUNCTION ctff(n, alb, anc, accx, upn) RESULT(fn_val)

!     Find CTFF so that P(QF > CTFF) < ACCX if UPN > 0;
!     P(QF < CTFF) < ACCX otherwise.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n(*)
REAL, INTENT(IN)      :: alb(*)
REAL, INTENT(IN)      :: anc(*)
REAL, INTENT(IN)      :: accx
REAL, INTENT(IN OUT)  :: upn
REAL                  :: fn_val

REAL     :: u1, u2, u, rb, const, c1, c2
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: zero = 0.0, one = 1.0, two = 2.0

u2 = upn
u1 = zero
c1 = amean
IF (.NOT.(u2 > zero)) GO TO 10
rb = almax
GO TO 20
10 rb = almin
20 rb = two*rb
u = u2/(one + u2*rb)
30 IF (.NOT.(errbd(n, alb, anc, u, c2) > accx)) GO TO 40
u1 = u2
c1 = c2
u2 = two*u2
u = u2/(one + u2*rb)
GO TO 30

40 u = (c1 - amean)/(c2 - amean)
50 IF (.NOT.(u < 0.9)) GO TO 80
u = (u1 + u2)/two
IF (.NOT.(errbd(n, alb, anc, u/(one + u*rb), const) > accx)) GO TO 60
u1 = u
c1 = const
GO TO 70

60 u2 = u
c2 = const
70 u = (c1 - amean)/(c2 - amean)
GO TO 50

80 fn_val = c2
upn = u2
RETURN
END FUNCTION ctff



FUNCTION truncn(n, alb, anc, uu, tausq) RESULT(fn_val)

!     Bound integration error due to truncation at U.

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n(*)
REAL, INTENT(IN)     :: alb(*)
REAL, INTENT(IN)     :: anc(*)
REAL, INTENT(IN)     :: uu
REAL, INTENT(IN)     :: tausq
REAL                 :: fn_val

REAL     :: u
REAL     :: sum1, sum2, prod1, prod2, prod3, alj, ancj, x, y, err1, err2
INTEGER  :: j, nj, ns
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: zero = 0.0, quart = 0.25, half = 0.5, one = 1.0, two = 2.0

CALL countr()
u = uu
sum1 = zero
prod2 = zero
prod3 = zero
ns = 0
sum2 = (sigsq + tausq)*u**2
prod1 = two*sum2
u = two*u
j = 1

10 IF (.NOT.(j <= ir)) GO TO 40
alj = alb(j)
ancj = anc(j)
nj = n(j)
x = (u*alj)**2
sum1 = sum1 + ancj*x/(one + x)
IF (.NOT.(x > one)) GO TO 20
prod2 = prod2 + nj*LOG(x)
prod3 = prod3 + nj*alog1(x,.true.)
ns = ns + nj
GO TO 30

20 prod1 = prod1 + nj*alog1(x,.true.)
30 j = j + 1
GO TO 10

40 sum1 = half*sum1
prod2 = prod1 + prod2
prod3 = prod1 + prod3
x = exp1(-sum1 - quart*prod2)/pi
y = exp1(-sum1 - quart*prod3)/pi
IF (.NOT.(ns == 0)) GO TO 50
err1 = one
GO TO 60

50 err1 = x*two/ns
60 IF (.NOT.(prod3 > one)) GO TO 70
err2 = 2.5*y
GO TO 80

70 err2 = one
80 IF (.NOT.(err2 < err1)) GO TO 90
err1 = err2
90 x = half*sum2
IF (.NOT.(x <= y)) GO TO 100
err2 = one
GO TO 110

100 err2 = y/x
110 IF (.NOT.(err1 < err2)) GO TO 120
fn_val = err1
GO TO 130

120 fn_val = err2
130 RETURN
END FUNCTION truncn



SUBROUTINE findu (n, alb, anc, utx, accx)

!     Find U such that TRUNCN(U) < ACCX & TRUNCN(U / 1.2) > ACCX.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n(*)
REAL, INTENT(IN)      :: alb(*)
REAL, INTENT(IN)      :: anc(*)
REAL, INTENT(IN OUT)  :: utx
REAL, INTENT(IN)      :: accx

REAL     :: u, ut
INTEGER  :: i
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: divis(4) = (/ 2.0, 1.4, 1.2, 1.1 /)
REAL, PARAMETER  :: four = 4.0, zero = 0.0

ut = utx
u = ut/four
IF (.NOT.(truncn(n, alb, anc, u, zero) > accx)) GO TO 20
u = ut
10 IF (.NOT.(truncn(n, alb, anc, u, zero) > accx)) GO TO 40
ut = ut*four
u = ut
GO TO 10

20 ut = u
u = u/four
30 IF (.NOT.(truncn(n, alb, anc, u, zero) <= accx)) GO TO 40
ut = u
u = u/four
GO TO 30

40 DO  i = 1,4
  u = ut/divis(i)
  IF (.NOT.(truncn(n, alb, anc, u, zero) <= accx)) CYCLE
  ut = u
END DO
utx = ut
RETURN
END SUBROUTINE findu



SUBROUTINE integr (n, alb, anc, nterm, aintrv, tausq, main)

!     Carry out integration with NTERM terms, at interval AINTRV.
!     If not MAIN then multiply integrand by  1 - exp(-0.5 * TAUSQ * U**2).

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n(*)
REAL, INTENT(IN)     :: alb(*)
REAL, INTENT(IN)     :: anc(*)
INTEGER, INTENT(IN)  :: nterm
REAL, INTENT(IN)     :: aintrv
REAL, INTENT(IN)     :: tausq
LOGICAL, INTENT(IN)  :: main

REAL      :: ainpi, u, sum1, sum2, sum3, x, y, z
INTEGER   :: k, j, nj
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: quart = 0.25, half = 0.5, one = 1.0, two = 2.0

ainpi = aintrv/pi
k = nterm
10 IF (.NOT.(k >= 0)) GO TO 50
u = (k + half)*aintrv
sum1 = -two*u*c
sum2 = ABS(sum1)
sum3 = -half*sigsq*u**2
j = ir

20 IF (.NOT.(j > 0)) GO TO 30
nj = n(j)
x = two*alb(j)*u
y = x**2
sum3 = sum3 - quart*nj*alog1(y,.true.)
y = anc(j)*x/(one + y)
z = nj*ATAN(x) + y
sum1 = sum1 + z
sum2 = sum2 + ABS(z)
sum3 = sum3 - half*x*y
j = j - 1
GO TO 20

30 x = ainpi*exp1(sum3)/u
IF (.NOT.(.NOT.main)) GO TO 40
x = x*(one - exp1(-half*tausq*u**2))
40 sum1 = SIN(half*sum1)*x
sum2 = half*sum2*x
aintl = aintl + sum1
ersm = ersm + sum2
k = k - 1
GO TO 10

50 RETURN
END SUBROUTINE integr



FUNCTION cfe(n, alb, anc, ith, x) RESULT(fn_val)

!     Coefficient of TAUSQ in error when convergence factor of
!     exp(-0.5 * TAUSQ * U**2) is used when df is evaluated at X.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n(*)
REAL, INTENT(IN)      :: alb(*)
REAL, INTENT(IN)      :: anc(*)
INTEGER, INTENT(OUT)  :: ith(*)
REAL, INTENT(IN)      :: x
REAL                  :: fn_val

REAL      :: axl, axl1, axl2, sxl, sum1, alj
INTEGER   :: j, k, it, itk
! COMMON /qfcom/ aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,  &
!    icount,ir,ndtsrt,fail,lim
REAL, PARAMETER  :: zero = 0.0, one = 1.0, four = 4.0, hundrd = 100.0

CALL countr()
IF (.NOT.(ndtsrt)) GO TO 10
CALL order (alb, ith)
10 axl = ABS(x)
sxl = SIGN(one,x)
sum1 = zero
j = ir

20 IF (.NOT.(j > 0)) GO TO 70
it = ith(j)
IF (.NOT.(alb(it)*sxl > zero)) GO TO 60
alj = ABS(alb(it))
axl1 = axl - alj*(n(it) + anc(it))
axl2 = alj/aln28
IF (.NOT.(axl1 > axl2)) GO TO 30
axl = axl1
GO TO 60

30 IF (.NOT.(axl > axl2)) GO TO 40
axl = axl2
40 sum1 = (axl - axl1)/alj
k = j - 1
50 IF (.NOT.(k > 0)) GO TO 70
itk = ith(k)
sum1 = sum1 + (n(itk) + anc(itk))
k = k - 1
GO TO 50

60 j = j - 1
GO TO 20

70 IF (.NOT.(sum1 > hundrd)) GO TO 80
fn_val = one
fail = .true.
GO TO 90
80 fn_val = 2**(sum1/four)/(pi*axl**2)
90 RETURN
END FUNCTION cfe

END MODULE Nonc_chisq
