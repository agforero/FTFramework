MODULE Solve_Complex_Poly
!     ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-08  Time: 16:02:44

! Latest revision - 10 February 2000
! amiller @ bigpond.net.au

!     ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02, P. 097.

USE quadruple_precision
IMPLICIT NONE

! COMMON AREA
! COMMON /global/ pr, pi, hr, hi, qpr, qpi, qhr, qhi, shr, shi, sr,  &
!    si, tr, ti, pvr, pvi, are, mre, eta, infin, nn
! TYPE (quad) :: sr, si, tr, ti, pvr, pvi, are, mre, eta, infin,  &
!    pr(50), pi(50), hr(50), hi(50), qpr(50), qpi(50), qhr(50),  &
!    qhi(50), shr(50), shi(50)

! N.B. pi has been changed to ci as the constant pi = 3.14159...
!      is available from the quad.precision module.

TYPE (quad), ALLOCATABLE, SAVE  :: pr(:), ci(:), hr(:), hi(:), qpr(:), qpi(:), &
                                   qhr(:), qhi(:), shr(:), shi(:)
TYPE (quad), SAVE               :: sr, si, tr, ti, pvr, pvi, are, mre, eta,  &
                                   infin
INTEGER, SAVE                   :: nn

PRIVATE
PUBLIC  :: cpoly


CONTAINS


SUBROUTINE cpoly(opr, opi, degree, zeror, zeroi, fail)
! FINDS THE ZEROS OF A COMPLEX POLYNOMIAL.
! OPR, OPI  -  TYPE (quad) VECTORS OF REAL AND IMAGINARY PARTS OF THE
!              COEFFICIENTS IN ORDER OF DECREASING POWERS.
! DEGREE    -  INTEGER DEGREE OF POLYNOMIAL.
! ZEROR, ZEROI  -  OUTPUT TYPE (quad) VECTORS OF REAL AND IMAGINARY
!                  PARTS OF THE ZEROS.
! FAIL      -  OUTPUT LOGICAL PARAMETER,  TRUE ONLY IF LEADING COEFFICIENT
!              IS ZERO OR IF CPOLY HAS FOUND FEWER THAN DEGREE ZEROS.
! THE PROGRAM HAS BEEN WRITTEN TO REDUCE THE CHANCE OF OVERFLOW OCCURRING.
! IF IT DOES OCCUR, THERE IS STILL A POSSIBILITY THAT THE ZEROFINDER WILL WORK
! PROVIDED THE OVERFLOWED QUANTITY IS REPLACED BY A LARGE NUMBER.

TYPE (quad), INTENT(IN)   :: opr(:)
TYPE (quad), INTENT(IN)   :: opi(:)
INTEGER, INTENT(IN)       :: degree
TYPE (quad), INTENT(OUT)  :: zeror(:)
TYPE (quad), INTENT(OUT)  :: zeroi(:)
LOGICAL, INTENT(OUT)      :: fail

TYPE (quad) :: xx, yy, cosr, sinr, smalno, base, xxx, zr, zi, bnd
LOGICAL     :: conv
INTEGER     :: cnt1, cnt2, i, idnn2

! INITIALIZATION OF CONSTANTS
CALL mcon(eta, infin, smalno, base)
are = eta
mre = 2.0_dp * SQRT(2.0_dp) * eta
xx = .70710678
yy = -xx
cosr = -.060756474
sinr = .99756405
fail = .false.
nn = degree + 1

! ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
IF (opr(1)%hi == 0.0_dp .AND. opi(1)%hi == 0.0_dp) THEN
  fail = .true.
  RETURN
END IF

! Allocate arrays

IF (ALLOCATED(pr)) DEALLOCATE(pr, ci, hr, hi, qpr, qpi, qhr, qhi, shr, shi)
ALLOCATE(pr(nn), ci(nn), hr(nn), hi(nn), qpr(nn), qpi(nn), qhr(nn),  &
         qhi(nn), shr(nn), shi(nn))

! REMOVE THE ZEROS AT THE ORIGIN IF ANY.
10 IF (opr(nn)%hi == 0.0_dp .AND. opi(nn)%hi == 0.0_dp) THEN
  idnn2 = degree - nn + 2
  zeror(idnn2) = 0.0_dp
  zeroi(idnn2) = 0.0_dp
  nn = nn - 1
  GO TO 10
END IF

! MAKE A COPY OF THE COEFFICIENTS.
DO  i = 1, nn
  pr(i) = opr(i)
  ci(i) = opi(i)
  shr(i) = cmod(pr(i), ci(i))
END DO

! SCALE THE POLYNOMIAL.
bnd = scalepoly(nn, shr, eta, infin, smalno, base)
IF (bnd%hi /= 1.0_dp) THEN
  DO  i = 1, nn
    pr(i) = bnd * pr(i)
    ci(i) = bnd * ci(i)
  END DO
END IF

! START THE ALGORITHM FOR ONE ZERO .
40 IF (nn <= 2) THEN

! CALCULATE THE FINAL ZERO AND RETURN.
  CALL cdivid(-pr(2), -ci(2), pr(1), ci(1), zeror(degree), zeroi(degree))
  RETURN
END IF

! CALCULATE BND, A LOWER BOUND ON THE MODULUS OF THE ZEROS.
DO  i = 1, nn
  shr(i) = cmod(pr(i), ci(i))
END DO
CALL cauchy(nn, shr, shi, bnd)

! OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES OF SHIFTS.
DO  cnt1 = 1, 2

! FIRST STAGE CALCULATION, NO SHIFT.
  CALL noshft(5)

! INNER LOOP TO SELECT A SHIFT.
  DO  cnt2 = 1, 9

! SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY
! 94 DEGREES FROM THE PREVIOUS SHIFT
    xxx = cosr * xx - sinr * yy
    yy = sinr * xx + cosr * yy
    xx = xxx
    sr = bnd * xx
    si = bnd * yy

! SECOND STAGE CALCULATION, FIXED SHIFT.
    CALL fxshft(10*cnt2,zr,zi,conv)
    IF (conv) THEN

! THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.
! IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.
      idnn2 = degree - nn + 2
      zeror(idnn2) = zr
      zeroi(idnn2) = zi
      nn = nn - 1
      DO  i = 1, nn
        pr(i) = qpr(i)
        ci(i) = qpi(i)
      END DO
      GO TO 40
    END IF

! IF THE ITERATION IS UNSUCCESSFUL ANOTHER SHIFT IS CHOSEN.
  END DO

! IF 9 SHIFTS FAIL, THE OUTER LOOP IS REPEATED WITH ANOTHER SEQUENCE OF SHIFTS.
END DO

! THE ZEROFINDER HAS FAILED ON TWO MAJOR PASSES.
! RETURN EMPTY HANDED.
fail = .true.
RETURN
END SUBROUTINE cpoly



SUBROUTINE noshft(l1)
! COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
! POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.

INTEGER, INTENT(IN)  :: l1

TYPE (quad) :: xni, t1, t2
INTEGER     :: i, j, jj, n, nm1

n = nn - 1
nm1 = n - 1
DO  i = 1, n
  xni = nn - i
  hr(i) = xni * pr(i) / n
  hi(i) = xni * ci(i) / n
END DO
DO  jj = 1, l1
  IF (cmod(hr(n), hi(n)) > eta*10.0_dp*cmod(pr(n), ci(n))) THEN
    CALL cdivid(-pr(nn), -ci(nn), hr(n), hi(n), tr, ti)
    DO  i = 1, nm1
      j = nn - i
      t1 = hr(j-1)
      t2 = hi(j-1)
      hr(j) = tr * t1 - ti * t2 + pr(j)
      hi(j) = tr * t2 + ti * t1 + ci(j)
    END DO
    hr(1) = pr(1)
    hi(1) = ci(1)
  ELSE

! IF THE CONSTANT TERM IS ESSENTIALLY ZERO, SHIFT H COEFFICIENTS.
    DO  i = 1, nm1
      j = nn - i
      hr(j) = hr(j-1)
      hi(j) = hi(j-1)
    END DO
    hr(1) = 0.0_dp
    hi(1) = 0.0_dp
  END IF
END DO
RETURN
END SUBROUTINE noshft



SUBROUTINE fxshft(l2, zr, zi, conv)
! COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
! INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
! APPROXIMATE ZERO IF SUCCESSFUL.
! L2 - LIMIT OF FIXED SHIFT STEPS
! ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
! CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION

INTEGER, INTENT(IN)       :: l2
TYPE (quad), INTENT(OUT)  :: zr
TYPE (quad), INTENT(OUT)  :: zi
LOGICAL, INTENT(OUT)      :: conv

TYPE (quad) :: otr, oti, svsr, svsi
LOGICAL     :: test, pasd, bool
INTEGER     :: i, j, n

n = nn - 1

! EVALUATE P AT S.
CALL polyev(nn, sr, si, pr, ci, qpr, qpi, pvr, pvi)
test = .true.
pasd = .false.

! CALCULATE FIRST T = -P(S)/H(S).
CALL calct(bool)

! MAIN LOOP FOR ONE SECOND STAGE STEP.
DO  j = 1, l2
  otr = tr
  oti = ti

! COMPUTE NEXT H POLYNOMIAL AND NEW T.
  CALL nexth(bool)
  CALL calct(bool)
  zr = sr + tr
  zi = si + ti

! TEST FOR CONVERGENCE UNLESS STAGE 3 HAS FAILED ONCE OR THIS
! IS THE LAST H POLYNOMIAL.
  IF (.NOT.(bool.OR..NOT.test.OR.j == l2)) THEN
    IF (cmod(tr-otr,ti-oti) < .5_dp*cmod(zr,zi)) THEN
      IF (pasd) THEN

! THE WEAK CONVERGENCE TEST HAS BEEN PASSED TWICE, START THE THIRD STAGE
! ITERATION, AFTER SAVING THE CURRENT H POLYNOMIAL AND SHIFT.
        DO  i = 1, n
          shr(i) = hr(i)
          shi(i) = hi(i)
        END DO
        svsr = sr
        svsi = si
        CALL vrshft(10,zr,zi,conv)
        IF (conv) RETURN

! THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND RESTORE H,S,PV AND T.
        test = .false.
        DO  i = 1, n
          hr(i) = shr(i)
          hi(i) = shi(i)
        END DO
        sr = svsr
        si = svsi
        CALL polyev(nn, sr, si, pr, ci, qpr, qpi, pvr, pvi)
        CALL calct(bool)
        CYCLE
      END IF
      pasd = .true.
    ELSE
      pasd = .false.
    END IF
  END IF
END DO

! ATTEMPT AN ITERATION WITH FINAL H POLYNOMIAL FROM SECOND STAGE.
CALL vrshft(10, zr, zi, conv)
RETURN
END SUBROUTINE fxshft



SUBROUTINE vrshft(l3, zr, zi, conv)
! CARRIES OUT THE THIRD STAGE ITERATION.
! L3 - LIMIT OF STEPS IN STAGE 3.
! ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
!           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
! CONV    -  .TRUE. IF ITERATION CONVERGES

INTEGER, INTENT(IN)          :: l3
TYPE (quad), INTENT(IN OUT)  :: zr
TYPE (quad), INTENT(IN OUT)  :: zi
LOGICAL, INTENT(OUT)         :: conv

TYPE (quad) :: mp, ms, omp, relstp, r1, r2, tp
LOGICAL     :: b, bool
INTEGER     :: i, j

conv = .false.
b = .false.
sr = zr
si = zi

! MAIN LOOP FOR STAGE THREE
DO  i = 1, l3

! EVALUATE P AT S AND TEST FOR CONVERGENCE.
  CALL polyev(nn, sr, si, pr, ci, qpr, qpi, pvr, pvi)
  mp = cmod(pvr, pvi)
  ms = cmod(sr,si)
  IF (mp <= 20.0_dp*errev(nn, qpr, qpi, ms, mp, are, mre)) THEN

! POLYNOMIAL VALUE IS SMALLER IN VALUE THAN A BOUND ON THE ERROR
! IN EVALUATING P, TERMINATE THE ITERATION.
    conv = .true.
    zr = sr
    zi = si
    RETURN
  END IF
  IF (i /= 1) THEN
    IF (.NOT.(b .OR. mp < omp .OR. relstp >= quad(.05_dp,0.0_dp))) THEN

! ITERATION HAS STALLED. PROBABLY A CLUSTER OF ZEROS.  DO 5 FIXED
! SHIFT STEPS INTO THE CLUSTER TO FORCE ONE ZERO TO DOMINATE.
      tp = relstp
      b = .true.
      IF (relstp < eta) tp = eta
      r1 = SQRT(tp)
      r2 = sr * (1.0_dp+r1) - si * r1
      si = sr * r1 + si * (1.0_dp+r1)
      sr = r2
      CALL polyev(nn, sr, si, pr, ci, qpr, qpi, pvr, pvi)
      DO  j = 1, 5
        CALL calct(bool)
        CALL nexth(bool)
      END DO
      omp = infin
      GO TO 20
    END IF

! EXIT IF POLYNOMIAL VALUE INCREASES SIGNIFICANTLY.
    IF (mp*.1_dp > omp) RETURN
  END IF
  omp = mp

! CALCULATE NEXT ITERATE.
  20 CALL calct(bool)
  CALL nexth(bool)
  CALL calct(bool)
  IF (.NOT.bool) THEN
    relstp = cmod(tr,ti) / cmod(sr,si)
    sr = sr + tr
    si = si + ti
  END IF
END DO
RETURN
END SUBROUTINE vrshft



SUBROUTINE calct(bool)
! COMPUTES  T = -P(S)/H(S).
! BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.

LOGICAL, INTENT(OUT)  :: bool

TYPE (quad) :: hvr, hvi
INTEGER     :: n

n = nn - 1

! EVALUATE H(S).
CALL polyev(n, sr, si, hr, hi, qhr, qhi, hvr, hvi)
bool = cmod(hvr,hvi) <= are * 10.0_dp * cmod(hr(n), hi(n))
IF (.NOT.bool) THEN
  CALL cdivid(-pvr, -pvi, hvr, hvi, tr, ti)
  RETURN
END IF
tr = 0.0_dp
ti = 0.0_dp
RETURN
END SUBROUTINE calct



SUBROUTINE nexth(bool)
! CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
! BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO

LOGICAL, INTENT(IN)  :: bool

TYPE (quad) :: t1, t2
INTEGER     :: j, n

n = nn - 1
IF (.NOT.bool) THEN
  DO  j = 2, n
    t1 = qhr(j-1)
    t2 = qhi(j-1)
    hr(j) = tr * t1 - ti * t2 + qpr(j)
    hi(j) = tr * t2 + ti * t1 + qpi(j)
  END DO
  hr(1) = qpr(1)
  hi(1) = qpi(1)
  RETURN
END IF

! IF H(S) IS ZERO REPLACE H WITH QH.
DO  j = 2, n
  hr(j) = qhr(j-1)
  hi(j) = qhi(j-1)
END DO
hr(1) = 0.0_dp
hi(1) = 0.0_dp
RETURN
END SUBROUTINE nexth



SUBROUTINE polyev(nn, sr, si, pr, ci, qr, qi, pvr, pvi)
! EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
! PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.

INTEGER, INTENT(IN)       :: nn
TYPE (quad), INTENT(IN)   :: sr
TYPE (quad), INTENT(IN)   :: si
TYPE (quad), INTENT(IN)   :: pr(:)
TYPE (quad), INTENT(IN)   :: ci(:)
TYPE (quad), INTENT(OUT)  :: qr(:)
TYPE (quad), INTENT(OUT)  :: qi(:)
TYPE (quad), INTENT(OUT)  :: pvr
TYPE (quad), INTENT(OUT)  :: pvi

TYPE (quad) :: t
INTEGER     :: i

qr(1) = pr(1)
qi(1) = ci(1)
pvr = qr(1)
pvi = qi(1)
DO  i = 2, nn
  t = pvr * sr - pvi * si + pr(i)
  pvi = pvr * si + pvi * sr + ci(i)
  pvr = t
  qr(i) = pvr
  qi(i) = pvi
END DO
RETURN
END SUBROUTINE polyev



FUNCTION errev(nn, qr, qi, ms, mp, are, mre) RESULT(fn_val)
! BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
! QR,QI - THE PARTIAL SUMS
! MS    -MODULUS OF THE POINT
! MP    -MODULUS OF POLYNOMIAL VALUE
! ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION

INTEGER, INTENT(IN)      :: nn
TYPE (quad), INTENT(IN)  :: qr(:)
TYPE (quad), INTENT(IN)  :: qi(:)
TYPE (quad), INTENT(IN)  :: ms
TYPE (quad), INTENT(IN)  :: mp
TYPE (quad), INTENT(IN)  :: are
TYPE (quad), INTENT(IN)  :: mre
TYPE (quad)              :: fn_val

TYPE (quad) :: e
INTEGER     :: i

e = cmod(qr(1), qi(1)) * mre / (are+mre)
DO  i = 1, nn
  e = e * ms + cmod(qr(i), qi(i))
END DO
fn_val = e * (are+mre) - mp * mre
RETURN
END FUNCTION errev



SUBROUTINE cauchy(nn, pt, q, fn_val)
! CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
! POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.

INTEGER, INTENT(IN)       :: nn
TYPE (quad), INTENT(OUT)  :: pt(:), q(:), fn_val

TYPE (quad) :: x, xm, f, dx, df
INTEGER     :: i, n

pt(nn) = -pt(nn)

! COMPUTE UPPER ESTIMATE OF BOUND.
n = nn - 1
x = EXP((LOG(-pt(nn)) - LOG(pt(1))) / n)
IF (pt(n)%hi /= 0.0_dp) THEN

! IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
  xm = -pt(nn) / pt(n)
  IF (xm < x) x = xm
END IF

! CHOP THE INTERVAL (0,X) UNITL F LE 0.
10 xm = x * .1_dp
f = pt(1)
DO  i = 2, nn
  f = f * xm + pt(i)
END DO
IF (f%hi > 0.0_dp) THEN
  x = xm
  GO TO 10
END IF
dx = x

! DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES.
30 IF (ABS(dx/x) > quad(.005_dp,0.0_dp)) THEN
  q(1) = pt(1)
  DO  i = 2, nn
    q(i) = q(i-1) * x + pt(i)
  END DO
  f = q(nn)
  df = q(1)
  DO  i = 2, n
    df = df * x + q(i)
  END DO
  dx = f / df
  x = x - dx
  GO TO 30
END IF
fn_val = x

RETURN
END SUBROUTINE cauchy



FUNCTION scalepoly(nn, pt, eta, infin, smalno, base) RESULT(fn_val)
! RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
! THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
! INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
! BASE.
! PT - MODULUS OF COEFFICIENTS OF P
! ETA, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.

INTEGER, INTENT(IN)      :: nn
TYPE (quad), INTENT(IN)  :: pt(:)
TYPE (quad), INTENT(IN)  :: eta
TYPE (quad), INTENT(IN)  :: infin
TYPE (quad), INTENT(IN)  :: smalno
TYPE (quad), INTENT(IN)  :: base
TYPE (quad)              :: fn_val

TYPE (quad) :: high, low, MAX, MIN, x, sc
INTEGER     :: i, l

! FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
high = SQRT(infin)
low = smalno / eta
MAX = 0.0_dp
MIN = infin
DO  i = 1, nn
  x = pt(i)
  IF (x > MAX) MAX = x
  IF (x%hi /= 0.0_dp .AND. x < MIN) MIN = x
END DO

! SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS.
fn_val = 1.0_dp
IF (MIN >= low .AND. MAX <= high) RETURN
x = low / MIN
IF (x%hi <= 1.0_dp) THEN
  sc = 1.0_dp / (SQRT(MAX)*SQRT(MIN))
ELSE
  sc = x
  IF (infin/sc > MAX) sc = 1.0_dp
END IF
l = LOG(sc) / LOG(base) + .5
fn_val = base ** l

RETURN
END FUNCTION scalepoly



SUBROUTINE cdivid(ar, ai, br, bi, cr, ci)
! COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.

TYPE (quad), INTENT(IN)   :: ar
TYPE (quad), INTENT(IN)   :: ai
TYPE (quad), INTENT(IN)   :: br
TYPE (quad), INTENT(IN)   :: bi
TYPE (quad), INTENT(OUT)  :: cr
TYPE (quad), INTENT(OUT)  :: ci

TYPE (quad) :: r, d, t, infin

IF (br%hi == 0.0_dp .AND. bi%hi == 0.0_dp) THEN

! DIVISION BY ZERO, C = INFINITY.
  CALL mcon(t, infin, t, t)
  cr = infin
  ci = infin
  RETURN
END IF
IF (ABS(br) < ABS(bi)) THEN
  r = br / bi
  d = bi + r * br
  cr = (ar*r+ai) / d
  ci = (ai*r-ar) / d
  RETURN
END IF
r = bi / br
d = br + r * bi
cr = (ar+ai*r) / d
ci = (ai-ar*r) / d
RETURN
END SUBROUTINE cdivid



FUNCTION cmod(r, i) RESULT(fn_val)
! MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.

TYPE (quad), INTENT(IN)  :: r
TYPE (quad), INTENT(IN)  :: i
TYPE (quad)              :: fn_val

TYPE (quad) :: ar, ai

ar = ABS(r)
ai = ABS(i)
IF (ar < ai) THEN
  fn_val = ai * SQRT(1.0_dp + (ar/ai)**2)
  RETURN
END IF
IF (ar > ai) THEN
  fn_val = ar * SQRT(1.0_dp + (ai/ar)**2)
  RETURN
END IF
fn_val = ar * SQRT(2.0_dp)
RETURN
END FUNCTION cmod



SUBROUTINE mcon(eta, infiny, smalno, base)
! MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE PROGRAM.
! THE USER MAY EITHER SET THEM DIRECTLY OR USE THE STATEMENTS BELOW TO
! COMPUTE THEM. THE MEANING OF THE FOUR CONSTANTS ARE -
! ETA       THE MAXIMUM RELATIVE REPRESENTATION ERROR WHICH CAN BE DESCRIBED
!           AS THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!           1.0_dp + ETA > 1.0.
! INFINY    THE LARGEST FLOATING-POINT NUMBER
! SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
! BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED

TYPE (quad), INTENT(OUT)  :: eta, infiny, smalno, base

base   = RADIX(0.0_dp)
eta    = EPSILON(0.0_dp) * 1.0E-15_dp
infiny = HUGE(0.0_dp)
smalno = TINY(0.0_dp)

RETURN
END SUBROUTINE mcon

END MODULE Solve_Complex_Poly



PROGRAM cpolydr
!              DRIVER TO TEST CPOLY

USE quadruple_precision
USE Solve_Complex_Poly
IMPLICIT NONE

LOGICAL     :: fail
TYPE (quad) :: p(50), ci(50), zr(50), zi(50)
INTEGER     :: i

WRITE (6,5000)
p(1) = 1
p(2) = -55
p(3) = 1320
p(4) = -18150
p(5) = 157773
p(6) = -902055
p(7) = 3416930
p(8) = -8409500
p(9) = 12753576
p(10) = -10628640
p(11) = 3628800
ci(1:11) = quad(0.0_dp,0.0_dp)
CALL prtc(11, p, ci)
CALL cpoly(p, ci, 10, zr, zi, fail)
IF (fail) GO TO 60
CALL prtz(10, zr, zi)

20 WRITE (6,5100)
p(1) = 1
p(2) = 0
p(3) = -10001.0001_dp
p(4) = 0
ci(1) = 0
ci(2) = -10001.0001_dp
ci(3) = 0
ci(4) = 1
CALL prtc(4, p, ci)
CALL cpoly(p, ci, 3, zr, zi, fail)
IF (fail) GO TO 70
CALL prtz(3,zr,zi)

30 WRITE (6,5200)
p(1) = 1.0
p(2) = -1.998046875
p(3) = 0.0
p(4) = .7567065954208374_dp
p(5) = -.2002119533717632_dp
p(6) = 1.271507365163416D-2
p(7) = 0
p(8) = -1.154642632172909D-5
p(9) = 1.584803612786345D-7
p(10) = -4.652065399568528D-10
p(11) = 0
ci(1) = 0
ci(2) = p(2)
ci(3) = 2.658859252929688_dp
ci(4) = -7.567065954208374D-1
ci(5) = 0
ci(6) = p(6)
ci(7) = -7.820779428584501D-4
ci(8) = -p(8)
ci(9) = 0
ci(10) = p(10)
ci(11) = 9.094947017729282D-13
CALL prtc(11, p, ci)
CALL cpoly(p, ci, 10, zr, zi, fail)
IF (fail) GO TO 80
CALL prtz(10, zr, zi)

40 WRITE (6,5300)
p(1) = 1
p(2) = -10
p(3) = 3
p(4) = 284
p(5) = -1293
p(6) = 2374
p(7) = -1587
p(8) = -920
p(9) = 2204
p(10) = -1344
p(11) = 288
ci(1) = 0
ci(2) = -10
ci(3) = 100
ci(4) = -334
ci(5) = 200
ci(6) = 1394
ci(7) = -3836
ci(8) = 4334
ci(9) = -2352
ci(10) = 504
ci(11) = 0
CALL prtc(11, p, ci)
CALL cpoly(p, ci, 10, zr, zi, fail)
IF (fail) GO TO 90
CALL prtz(10, zr, zi)

50 WRITE (6,5400)
p(1) = 1
p(2) = 0
p(3) = -264
p(4) = 0
p(5) = 7920
p(6) = 0
p(7) = -59136
p(8) = 0
p(9) = 126720
p(10) = 0
p(11) = -67584
p(12) = 0
p(13) = 4095
ci(1) = 0
ci(2) = -24
ci(3) = 0
ci(4) = 1760
ci(5) = 0
ci(6) = -25344
ci(7) = 0
ci(8) = 101376
ci(9) = 0
ci(10) = -112640
ci(11) = 0
ci(12) = 24576
ci(13) = 0
CALL prtc(13, p, ci)
CALL cpoly(p, ci, 12, zr, zi, fail)
IF (fail) GO TO 100
CALL prtz(12, zr, zi)
STOP

60 WRITE (6,5500)
GO TO 20
70 WRITE (6,5500)
GO TO 30
80 WRITE (6,5500)
GO TO 40
90 WRITE (6,5500)
GO TO 50
100 WRITE (6,5500)
STOP

5000 FORMAT ('1EXAMPLE 1. POLYNOMIAL WITH ZEROS 1,2,...,10.')
5100 FORMAT ('1EXAMPLE 2. ZEROS ON IMAGINARY AXIS DEGREE 3.')
5200 FORMAT ('1EXAMPLE 3. ZEROS AT 1+I,1/2*(1+I)....1/(2**-9)*(1+I)')
5300 FORMAT ('1EXAMPLE 4. MULTIPLE ZEROS')
5400 FORMAT ('1EXAMPLE 5. 12 ZEROS EVENLY DISTRIBUTED ON A CIRCLE OF',  &
             ' RADIUS 1. CENTERED AT 0+2I.')
5500 FORMAT (//' CPOLY HAS FAILED ON THIS EXAMPLE')


CONTAINS


SUBROUTINE prtc(n, p, q)

INTEGER, INTENT(IN)      :: n
TYPE (quad), INTENT(IN)  :: p(50)
TYPE (quad), INTENT(IN)  :: q(50)

WRITE (6,5000) (p(i)%hi,q(i)%hi,i = 1,n)
RETURN
5000 FORMAT (//' COEFFICIENTS'/ 50(2G26.16/))
END SUBROUTINE prtc


SUBROUTINE prtz(n, zr, zi)

INTEGER, INTENT(IN)      :: n
TYPE (quad), INTENT(IN)  :: zr(50)
TYPE (quad), INTENT(IN)  :: zi(50)

WRITE (6,5000) (zr(i)%hi,zi(i)%hi,i = 1,n)
RETURN
5000 FORMAT (//' ZEROS'/50(2G26.16/))
END SUBROUTINE prtz

END PROGRAM cpolydr
