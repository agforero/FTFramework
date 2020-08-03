MODULE toms519

! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-15  Time: 00:33:15

IMPLICIT NONE
INTEGER, PARAMETER             :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), PARAMETER, PRIVATE  :: one = 1.0_dp, zero = 0.0_dp


CONTAINS


FUNCTION pln(n, s2, b2, s1, b1) RESULT(fn_val)
 
! COMPUTES PROB FOR BOUNDARIES OF FORM S*Y + B.  N IS SAMPLE SIZE,
! S2,B2 ARE PARAMETERS FOR UPPER BOUNDARY.
! CALLS SUBROUTINE RAKK.
! TO CALL DURB, CHANGE AS INDICATED BELOW.

IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: s2
REAL (dp), INTENT(IN)  :: b2
REAL (dp), INTENT(IN)  :: s1
REAL (dp), INTENT(IN)  :: b1
REAL (dp)              :: fn_val

REAL (dp)  :: yy(2,202), s(2), b(2), an, ans, c
INTEGER    :: i, k, n1

! MAX SIZE FOR N IS 200
s(1) = s1
s(2) = s2
b(1) = b1
b(2) = b2
IF (n <= 0) THEN
  fn_val = zero
  RETURN
END IF
IF (s(1) < zero) THEN
  s(1) = zero
END IF
IF (s(2) < zero) THEN
  b(2) = s(2) + b(2)
  s(2) = zero
END IF
n1 = 1 + n
an = n
DO  i = 1, 2
  DO  k = 1, n1
    c = REAL(k-1) / an - b(i)
    yy(i,k) = -0.5
    IF (c > zero) THEN
      yy(i,k) = 1.1
      IF (s(i) > c) yy(i,k) = c / s(i)
    END IF
  END DO
END DO
CALL rakk(n, yy, ans)

! TO CALL DURB CHANGE PRECEDING STATEMENT TO
!    CALL DURB (N, YY, ANS)
fn_val = ans

RETURN
END FUNCTION pln



SUBROUTINE rakk(n, yy, ans)
! N IS SAMPLE SIZE.  YY(2,.),YY(1,.) CONTAIN, RESPECTIVELY,
! UPPER AND LOWER BOUNDARY LEVEL POINTS CORRESPONDING TO 0,1/N,2/N,...,1.
! MUST BE N+1 OF EACH.  ROUTINE ASSUMES THEY ARE NONDECREASING.

IMPLICIT NONE

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: yy(2,202)
REAL (dp), INTENT(OUT)     :: ans

REAL (dp)  :: convi, frng(41,2), p, rrngs, sss(41), s, sqrng, sqrrng, t, tf, &
              tmlt, w(201,2), y0, y1
INTEGER    :: i, i1, iext, j, j1, j11, jbot0, jbotp, jtop0, jtop1, jtopp,  &
              k, k0, k01, k0p1, k1, k1p1, kbli, kbloc(41), l, l0, l1, lb,  &
              lb1, ldf, lfg, ll, lt, m1, m2, mbav, mt, n1, nbloc, nnn
REAL (dp), PARAMETER  :: rng = HUGE(one), rrng = TINY(one)
INTEGER, PARAMETER    :: maxnb = 40

! DATA rng, rrng, sqrng, sqrrng, maxnb /1.e292, 1.e-292, 1.e146, 1.e-146, 40/
! RNG AND RRNG ARE MACHINE DEPENDENT CONSTANTS SELECTED
! SO RRNG=1./RNG >= SMALLEST POSITIVE FLOATING
! NUMBER AND RNG <= LARGEST POSITIVE FLOATING
! NUMBER OF MACHINE.  SQRNG AND SQRRNG ARE SQUARE ROOTS OF RNG AND RRNG.
! MAXNB IS THE MAX NUMBER OF BLOCKS INTO WHICH W(.,.) WILL BE PARTITIONED.
! IF MAXNB IS INCREASED REVISE DIMENSION STATEMENT ACCORDINGLY.

sqrng = SQRT(rng)
sqrrng = SQRT(rrng)
n1 = n + 1
tf = one
nnn = 1
iext = 1
m2 = 1
m1 = 2
w(1,m2) = sqrng
nbloc = 1
kbloc(1) = 1
kbloc(2) = 2
frng(1,m2) = one
ans = zero
IF (n1 <= 1) RETURN
yy(2,n1+1) = one
y1 = zero
IF (yy(1,1) <= zero .OR. yy(2,n1) >= one .OR. yy(1,n1) <= one .OR.   &
    yy(2,1) >= zero) RETURN
jbotp = 1
jtopp = 1
jtop1 = 1

! ENTRY POINT FOR SUBSEQUENT ITERATIONS.
10 jbot0 = jbotp
jtop0 = jtop1
mt = m1
m1 = m2
m2 = mt
y0 = y1

! COMPUTE SUMMATION INDICES, JBOTP,JTOPP FOR NEXT ITERATION.
20 IF (yy(2,jtopp+1) <= y0) THEN
  jtopp = jtopp + 1
  GO TO 20
END IF
y1 = MIN(yy(2,jtopp+1), yy(1,jbotp))
IF (y1 < one) THEN
  30 IF (yy(1,jbotp) > y1) GO TO 40
  jbotp = jbotp + 1
  GO TO 30
END IF
y1 = one

! SET EXIT FLAG.  NEXT ITERATION IS FINAL ONE.
iext = 2
jbotp = n1
40 IF (jbotp > jtopp) RETURN
! RETURN WITH ANS=0 SINCE NO PATHS BETWEEN BOUNDARIES.

jtop1 = jtopp
p = y1 - y0
kbloc(nbloc+1) = jtop1 + 1
DO  l = jbot0, jtop1
  w(l,m2) = zero
END DO
l0 = kbloc(nbloc)
l1 = jtop0
DO  l = l0, l1
  IF (w(l,m1) <= sqrrng) THEN
    jtop0 = l - 1
    EXIT
  END IF
END DO
70 IF (jbot0 >= kbloc(2)) THEN

! THE NEXT FEW STATEMENTS ELIMINATE BLOCKS NO LONGER USED.
  DO  i = 1, nbloc
    frng(i,m1) = frng(i+1,m1)
    kbloc(i) = kbloc(i+1)
  END DO
  nbloc = nbloc - 1
  CALL accum(rrng/frng(1,m1),tf,nnn)
  GO TO 70
END IF
kbloc(1) = jbot0
mbav = (maxnb-nbloc+1) / 2
IF (p <= (one/REAL(4*n1))) mbav = 2
mbav = MIN(maxnb, nbloc+mbav)
! MBAV IS UPPER LIMIT ON NUMBER OF BLOCKS DURING NEXT ITERATION.
! MBAV IS SET TO LIMIT NUMBER OF NEW BLOCKS WHEN P IS SMALL.

convi = one
i1 = 1
j1 = 1
90 i = i + 1                 ! Dummy instruction
loop320:  &
DO  i = i1, nbloc
  kbli = kbloc(i)
  IF (i /= i1) THEN
    frng(i,m2) = MIN(frng(i,m1),one/(sqrng*w(kbli-1,m2)))
    convi = (convi/frng(i-j1+1,m1)) * frng(i,m2)
    IF (j1 /= 1) THEN
      k1 = kbloc(i) - kbloc(i-j1+2) + 1
      k0 = kbloc(i-1) - kbloc(i-j1+1) + 1
      IF (k1 < k0) THEN
        k1p1 = k1 + 1
        DO  k = k1p1, k0
          convi = convi * REAL(k) / p
        END DO
      ELSE IF (k1 > k0) THEN
        k0p1 = k0 + 1
        DO  k = k0p1, k1
          convi = convi * p / REAL(k)
        END DO
      END IF
    END IF
  END IF

  IF (convi <= rrng) THEN
    DO  j = j1, j11
! CHANGE UPPER LIMIT IF EARLY EXIT FROM PREVIOUS LOOP
      convi = frng(i,m2) * rng
      k0 = MAX(1, kbloc(i-1) - MIN(jtop0+1, kbloc(i-j+1))+2)
      k1 = kbloc(i) - MIN(jtop0+1, kbloc(i-j+1)) + 1
      DO  k = k0, k1
        convi = convi * p / REAL(k)
      END DO
      convi = sss(j) * convi
      IF (convi > rrng) GO TO 170
    END DO
    SELECT CASE ( iext )
      CASE (    1)
        GO TO 10
      CASE (    2)
        GO TO 360
    END SELECT
    170 j1 = j + 1
  END IF

  DO  j = j1, i
    lfg = 1
! LFG=2 IF NEXT ITERATION ON J BECOMES NECESSARY.
    lb1 = kbloc(i)
    k01 = kbloc(i) - kbloc(i-j+1) + 1
    s = convi
    IF (j /= j1) THEN
      IF (tmlt == zero) THEN
        CYCLE loop320
      ELSE IF (tmlt < zero) THEN
        s = -tmlt * frng(i-j+2,m1)
        IF (s > rrng) GO TO 200
      END IF
      s = tmlt * (frng(i-j+2,m1)*rng)
    END IF

    200 k0 = MAX(kbloc(i)-MIN(kbloc(i-j+2)-1,jtop0)+1, 1)
    k1 = kbloc(i+1) - kbloc(i-j+1)
    tmlt = zero
    rrngs = one
    IF (s < zero) rrngs = -rrng
    sss(j) = s
    DO  k = k0, k1
      ldf = 1 - k
      lb = MAX(kbloc(i),kbloc(i-j+1)-ldf,lb1)
      LT = MIN(kbloc(i+1)-1, MIN(jtop0,kbloc(i-j+2)-1)-ldf)
      IF (lb > lb1) lfg = 2
      IF (lb <= LT) THEN
        DO  l = lb, LT
          IF (w(l,m2) /= (w(l+ldf,m1)*rrngs)*s + w(l,m2)) THEN
            lb1 = l

! WHEN THE SUMMANDS BECOME SO SMALL THAT THE VALUE OF W(L,M2) IS UNCHANGED,
! LB1 IS INCREASED.  THIS PREVENTS FURTHER ADDITIONS TO SUCH W(L,M2).
            GO TO 230
          END IF

          IF (w(l,m2) <= rrng) THEN
            IF (k-k01 > 0) THEN
              GO TO 220
            ELSE
              GO TO 280
            END IF
          END IF
        END DO
      END IF

      lb1 = LT + 1
      IF (lb1 /= kbloc(i+1) .OR. k <= k01) GO TO 280
      220 j11 = j
      IF (w(kbli,m2) <= rrng) lfg = 2
      SELECT CASE ( lfg )
        CASE (    1)
          CYCLE loop320
        CASE (    2)
          GO TO 300
      END SELECT

      230 IF (s < zero) THEN
        DO  l = lb1, LT
          w(l,m2) = (w(l+ldf,m1)*rrngs) * s + w(l,m2)
        END DO
      ELSE IF (s > zero) THEN
        DO  l = lb1, LT
          w(l,m2) = w(l+ldf,m1) * s + w(l,m2)
        END DO
      ELSE
        GO TO 300
      END IF

      280 t = s
      s = s * p / REAL(k)
      IF (s == 0) THEN
        IF (t <= zero) EXIT
        s = -(t*rng) * p / REAL(k)
        rrngs = -rrng
      END IF

! THE NEGATIVE BIT IN S,TMLT IS USED TO
! INDICATE ITS VALUE IS RNG TIMES THE USUAL VALUE.
      IF (k == k01) tmlt = s

! THE ITERATION ON K MUST CONTINUE UNTIL TMLT IS SET.
! TMLT USED IN NEXT ITERATION ON J.
    END DO
    300 j11 = j
  END DO
END DO loop320

! CREATE NEW BLOCK IF NECESSARY.
IF (nbloc >= mbav) THEN
   SELECT CASE ( iext )
    CASE (    1)
      GO TO 10
    CASE (    2)
      GO TO 360
  END SELECT
END IF

l1 = kbloc(nbloc+1) - kbloc(nbloc)
DO  l = 1, l1
  ll = kbloc(nbloc+1) - l
  IF (w(ll,m2) >= sqrrng) EXIT
END DO

IF (ll == jtop1) THEN
  SELECT CASE ( iext )
    CASE (    1)
      GO TO 10
    CASE (    2)
      GO TO 360
  END SELECT
END IF

kbloc(nbloc+2) = kbloc(nbloc+1)
kbloc(nbloc+1) = ll + 1
nbloc = nbloc + 1
l0 = kbloc(nbloc)
DO  l = l0, jtop1
  w(l,m2) = zero
END DO
i1 = nbloc
frng(nbloc,m2) = one
convi = zero
GO TO 90

360 CALL accum(w(n1,m2), tf, nnn)
CALL accum(sqrrng,tf,nnn)
IF (nbloc /= 1) THEN
  DO  i = 2, nbloc
    CALL accum(rrng/frng(i,m2), tf, nnn)
  END DO
END IF
ans = tf
IF (n < nnn) THEN
  DO  j = n1, nnn
    ans = ans / REAL(j)
  END DO
ELSE IF (n > nnn) THEN
  nnn = nnn + 1
  DO  j = nnn, n
    ans = ans * REAL(j)
  END DO
END IF

RETURN
END SUBROUTINE rakk



SUBROUTINE accum(f, tf, nnn)
! THIS SUBROUTINE ANCILLARY TO RAKK AND ACCUMULATES THE
! MULTIPLIERS USED TO KEEP W(.,.) IN RANGE.

IMPLICIT NONE

REAL (dp), INTENT(IN)      :: f
REAL (dp), INTENT(IN OUT)  :: tf
INTEGER, INTENT(IN OUT)    :: nnn

tf = tf * f
DO
  IF (tf > one .OR. tf == zero) EXIT
  nnn = nnn + 1
  tf = tf * REAL(nnn)
END DO

RETURN
END SUBROUTINE accum



SUBROUTINE durb(n,yy,ans)
! DURBINS METHOD.  J.APPL.PROB., 8(1971), PP431-453,
! FORMULAS (36),(37),(40).
! N IS SAMPLE SIZE.  YY(2,.),YY(1,.)
! CONTAIN, RESPECTIVELY, UPPER AND LOWER BOUNDARY
! LEVEL POINTS CORRESPONDING TO 0,1/N,2/N,...,1.  MUST
! BE N+1 OF EACH.  ROUTINE ASSUMES THEY ARE NONDECREASING.

IMPLICIT NONE

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: yy(2,202)
REAL (dp), INTENT(OUT)     :: ans

REAL (dp)  :: alpha(201), beta(201)
! DIMENSION STATEMENT PERMITS MAX SAMPLE SIZE OF 200

INTEGER    :: i, i1, ip, j, k, km1, n1
REAL (dp)  :: comb

n1 = n + 1
ans = zero
IF (yy(1,1) <= zero .OR. yy(2,n1) >= one .OR. yy(1,n1) <= one .OR.   &
    yy(2,1) >= zero) RETURN
ip = 1
DO  i = 1, n1
  IF (yy(2,i) >= zero) THEN
    ip = i
    EXIT
  END IF
  yy(2,i) = zero
END DO
DO  i = 1, n1
  IF (yy(1,i) > one) THEN
    yy(1,i) = one
  END IF
END DO

alpha(1) = zero
beta(1) = one
loop90:  DO  k = 2, n1
  alpha(k) = yy(2,k)**(k-1)
  IF (ip <= (k-1)) THEN
    comb = one
    i1 = k - ip
    DO  i = 1, i1
      comb = comb * REAL(k-i) / REAL(i)
! COMB IS BINOMIAL COEFFICIENT K-1 OVER I-1.
      alpha(k) = alpha(k) - alpha(k-i) * (yy(2,k) - yy(2,k-i))**i * comb
    END DO
  END IF
  comb = one
  DO  j = 1, n1
    IF (yy(2,k) <= yy(1,j)) EXIT
    alpha(k) = alpha(k) - beta(j) * (yy(2,k) - yy(1,j))**(k-j) * comb
! COMB IS BINOMIAL COEFFICIENT K-1 OVER J-1.
    comb = comb * REAL(k-j) / REAL(j)
  END DO
  beta(k) = yy(1,k) ** (k-1)
  IF (ip <= k) THEN
    comb = one
    i1 = k - ip + 1
    DO  i = 1, i1
      beta(k) = beta(k) - alpha(k-i+1) * (yy(1,k) - yy(2,k-i+1))**(i-1) * comb
! COMB IS BINOMIAL COEFFICIENT K-1 OVER I-1.
      comb = comb * REAL(k-i) / REAL(i)
    END DO
    comb = one
  END IF
  km1 = k - 1
  DO  j = 1, km1
    IF (yy(1,j) == one) CYCLE loop90
! COMB IS BINOMIAL COEFFICIENT K-1 OVER J-1.
    beta(k) = beta(k) - beta(j) * (yy(1,k) - yy(1,j))**(k-j) * comb
    comb = comb * REAL(k-j) / REAL(j)
  END DO
END DO loop90

comb = one
i1 = n1 - ip + 1
DO  i = 1, i1
  ans = ans + (one - yy(2,n1-i+1))**(i-1) * alpha(n1-i+1) * comb
  comb = comb * REAL(n1-i) / REAL(i)
END DO
comb = one
DO  j = 1, n1
  IF (yy(1,j) == one) EXIT
  ans = ans + (one - yy(1,j))**(n1-j) * beta(j) * comb
! COMB IS BINOMIAL COEFFICIENT N OVER J-1.
  comb = comb * REAL(n1-j) / REAL(j)
END DO

! FOR PROB OF COMPLEMENTARY EVENT CHANGE NEXT STATEMENT TO
! CONTINUE
ans = one - ans

RETURN
END SUBROUTINE durb

END MODULE toms519
