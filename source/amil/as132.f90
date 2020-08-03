SUBROUTINE simlp(n, x, y, sad, alpha, beta, d, iter, inext, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-19  Time: 09:37:54

!     ALGORITHM AS 132  APPL. STATIST. (1978) VOL.27, NO.3

!     SIMLP:   Fit  Y = ALPHA + BETA.X + error

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(n)
REAL, INTENT(IN)      :: y(n)
REAL, INTENT(OUT)     :: sad
REAL, INTENT(OUT)     :: alpha
REAL, INTENT(OUT)     :: beta
REAL, INTENT(OUT)     :: d(n)
INTEGER, INTENT(OUT)  :: iter
INTEGER, INTENT(OUT)  :: inext(n)
INTEGER, INTENT(OUT)  :: ifault

REAL, PARAMETER  :: acu = 1.0E-06, big = 1.0E19, half = 0.5, zero = 0.0,   &
                    one = 1.0, two = 2.0
INTEGER  :: i, ibas1, ibas2, iflag, iin, iout, isave, j
REAL     :: a1, a2, aaa, aaaa, ahalf, aone, bbb, bbbb, ddd, det,   &
            ratio, rho, subt, sum, t, test, tot1, tot2, y1, y2, zzz

!     Initial settings

ifault = 0
iter = 0
ahalf = half + acu
aone = ahalf + ahalf

!     Determine initial basis

d(1) = zero
y1 = y(1)
ibas1 = 1
a1 = x(1)
DO  i = 2, n
  IF (ABS(a1-x(i)) >= acu) THEN
    a2 = x(i)
    ibas2 = i
    y2 = y(i)
    GO TO 20
  END IF
END DO
ifault = 1
RETURN

!     Calculate initial beta value

20 det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det

!     Calculate initial D-vector

DO  i = 2, n
  ddd = y(i) - (aaaa+bbbb*x(i))
  d(i) = SIGN(one,ddd)
END DO
tot1 = one
tot2 = x(ibas2)
d(ibas2) = -one
DO  i = 2, n
  tot1 = tot1 + d(i)
  tot2 = tot2 + d(i) * x(i)
END DO
t = (a2*tot1-tot2) * det
IF (ABS(t) >= aone) THEN
  det = -det
  GO TO 70
END IF

!     Main iterative loop begins

50 t = (tot2-a1*tot1) * det
IF (ABS(t) < aone) GO TO 130
iflag = 2
iout = ibas2
x(iout) = a1
aaa = a1
bbb = a2
GO TO 80

60 t = (tot2-a2*tot1) * det
IF (ABS(t) < aone) GO TO 130

70 iflag = 1
bbb = a1
aaa = a2
iout = ibas1

80 rho = SIGN(one,t)
t = half * ABS(t)
det = det * rho

!     Perform partial sort of ratios

inext(ibas1) = ibas2
ratio = big
sum = ahalf
DO  i = 1, n
  ddd = (x(i)-aaa) * det
  IF (ddd*d(i) > acu) THEN
    test = (y(i)-aaaa-bbbb*x(i)) / ddd
    IF (test < ratio) THEN
      j = ibas1
      sum = sum + ABS(ddd)
      90 isave = ABS(inext(j))
      IF (test < d(isave)) THEN
        IF (sum >= t) THEN
          subt = ABS((x(isave)-aaa)*det)
          IF (sum-subt >= t) THEN
            sum = sum - subt
            d(isave) = SIGN(1,inext(j))
            inext(j) = inext(isave)
            GO TO 90
          END IF
        END IF
        100 j = isave
        isave = ABS(inext(j))
        IF (test < d(isave)) GO TO 100
      END IF
      inext(i) = inext(j)
      inext(j) = SIGN(i,INT(d(i)))
      d(i) = test
      IF (sum >= t) THEN
        iin = ABS(inext(ibas1))
        ratio = d(iin)
      END IF
    END IF
  END IF
END DO

!     Update basic indicators

iin = ABS(inext(ibas1))
j = iin
120 isave = ABS(inext(j))
IF (isave /= ibas2) THEN
  zzz = SIGN(1,inext(j))
  tot1 = tot1 - zzz - zzz
  tot2 = tot2 - two * zzz * x(isave)
  d(isave) = -zzz
  j = isave
  GO TO 120
END IF
zzz = SIGN(1,inext(ibas1))
tot1 = tot1 - rho - zzz
tot2 = tot2 - rho * bbb - zzz * x(iin)
d(iout) = -rho
iter = iter + 1
IF (iflag /= 1) THEN
  x(ibas2) = a2
  ibas2 = iin
  d(ibas2) = -one
  a2 = x(iin)
  y2 = y(iin)
  det = one / (a1-a2)
  aaaa = (a1*y2-a2*y1) * det
  bbbb = (y1-y2) * det
  GO TO 60
END IF
ibas1 = iin
a1 = x(iin)
d(ibas1) = zero
y1 = y(iin)
det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det
GO TO 50

!     Calculate optimal sum of absolute deviations

130 sad = zero
DO  i = 1, n
  d(i) = y(i) - (aaaa+bbbb*x(i))
  sad = sad + ABS(d(i))
END DO
alpha = aaaa
beta = bbbb

RETURN
END SUBROUTINE simlp
