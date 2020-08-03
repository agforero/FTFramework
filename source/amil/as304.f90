MODULE Fisher_Test
IMPLICIT NONE

CONTAINS


SUBROUTINE fisher (x, m, y, n, total, possib, p, ifault)

!     ALGORITHM AS 304.1 APPL.STATIST. (1996), VOL.45, NO.3

!     Fisher's non-parametric randomization test for two small
!     independent random samples

IMPLICIT NONE
INTEGER, INTENT(IN OUT)  :: m, n
INTEGER, INTENT(OUT)     :: ifault, total, possib
REAL, INTENT(IN OUT)     :: x(*), y(*)
REAL, INTENT(OUT)        :: p

INTEGER, PARAMETER  :: maxsam = 14, maxsiz = 3432

!        Important : set MAXSIZ >= COMB(MAXSAM, MAXSAM / 2)

INTEGER  :: k, size1, size2, ws3(maxsam)
REAL     :: sumx, sumy, ws1(maxsiz), ws2(maxsiz)

! mean(summ, count) = summ / REAL(count)

IF (comb(maxsam, maxsam/2) > maxsiz) THEN
  ifault = 1
ELSE IF (m > maxsam .OR. n > maxsam) THEN
  ifault = 2
ELSE
  ifault = 0
  sumx = SUM( x(1:m) )
  sumy = SUM( y(1:n) )
  IF (mean(sumx, m) > mean(sumy, n)) THEN
    CALL exchng(x, m, y, n, sumx, sumy)
  END IF
  total = 0
  IF (m == n) THEN
    DO k = 1, (m-1)/2
      CALL ktrade(x, m, ws1, size1, ws3, k)
      CALL ktrade(y, n, ws2, size2, ws3, k)
      total = total + trades(ws1, size1, ws2, size2)
      CALL cmplmt(ws1, size1, sumx)
      CALL cmplmt(ws2, size2, sumy)
      total = total + trades(ws1, size1, ws2, size2)
    END DO
    IF (MOD(m, 2) == 0) THEN
      CALL ktrade(x, m, ws1, size1, ws3, k)
      CALL ktrade(y, n, ws2, size2, ws3, k)
      total = total + trades(ws1, size1, ws2, size2)
    END IF
  ELSE
    DO k = 1, MIN(m, n)
      CALL ktrade(x, m, ws1, size1, ws3, k)
      CALL ktrade(y, n, ws2, size2, ws3, k)
      total = total + trades(ws1, size1, ws2, size2)
    END DO
  END IF
  possib = comb(m+n, m)
  p = REAL(total + 1) / REAL(possib)
END IF

RETURN
END SUBROUTINE fisher



FUNCTION mean(summ, count) RESULT(fn_val)
IMPLICIT NONE

REAL, INTENT(IN)     :: summ
INTEGER, INTENT(IN)  :: count
REAL                 :: fn_val

fn_val = summ / REAL(count)

RETURN
END FUNCTION mean


SUBROUTINE exchng (x, m, y, n, sx, sy)

!     ALGORITHM AS 304.2 APPL.STATIST. (1996), VOL.45, NO.3

!     Exchanges the sample data.  Assumes both X and Y have been
!     previously dimensioned to at least max(M, N) elements

IMPLICIT NONE
INTEGER, INTENT(IN OUT)  :: m, n
REAL, INTENT(IN OUT)     :: x(*), y(*), sx, sy

INTEGER  :: c, k
REAL     :: temp

temp = sx
sx = sy
sy = temp

c = MIN(m, n)
DO k = 1, c
  temp = x(k)
  x(k) = y(k)
  y(k) = temp
END DO
IF (m > n) THEN
  DO k = c+1, m
    y(k) = x(k)
  END DO
  n = m
  m = c
ELSE IF (m < n) THEN
  DO k = c+1, n
    x(k) = y(k)
  END DO
  m = n
  n = c
END IF

RETURN
END SUBROUTINE exchng



SUBROUTINE ktrade (w, k, wprime, kprime, ws, r)

!   ALGORITHM AS 304.3 APPL.STATIST. (1996), VOL.45, NO.3

!   Generates and sorts the sums of the R-combinations of the elements of W.

INTEGER, INTENT(IN)      :: k, r
INTEGER, INTENT(IN OUT)  :: ws(*)
INTEGER, INTENT(OUT)     :: kprime
REAL, INTENT(IN OUT)     :: w(*)
REAL, INTENT(OUT)        :: wprime(*)

kprime = comb(k, r)
IF (r <= k - r .OR. r == k) THEN
  CALL gener(w, k, wprime, kprime, ws, r)
  CALL sort(wprime, kprime)
ELSE
  CALL gener(w, k, wprime, kprime, ws, k - r)
  CALL sort(wprime, kprime)
  CALL cmplmt(wprime, kprime, SUM( w(1:k) ))
END IF

RETURN
END SUBROUTINE ktrade



FUNCTION trades (xprime, mprime, yprime, nprime) RESULT(fn_val)

!     ALGORITHM AS 304.4 APPL.STATIST. (1996), VOL.45, NO.3

!     Returns the number of 1-for-1 trades that refutes the null hypothesis.
!     Assumes that XPRIME has the smaller mean and
!     that both arrays are sorted in ascending order.

INTEGER, INTENT(IN)  :: mprime, nprime
REAL, INTENT(IN)     :: xprime(*), yprime(*)
INTEGER              :: fn_val

INTEGER :: i, j

fn_val = 0
i = 1
j = 1
10 IF (j > nprime) GO TO 40
20 IF (xprime(i) >= yprime(j)) GO TO 30
i = i + 1
IF (i <= mprime) GO TO 20
30 fn_val = fn_val + (mprime - i + 1)
j = j + 1
IF (i <= mprime) GO TO 10

40 RETURN
END FUNCTION trades



SUBROUTINE cmplmt (wprime, kprime, sum)

!        ALGORITHM AS 304.5 APPL.STATIST. (1996), VOL.45, NO.3

!        Reverse and complement the data in WPRIME

IMPLICIT NONE
INTEGER, INTENT(IN)   :: kprime
REAL, INTENT(IN)      :: sum
REAL, INTENT(IN OUT)  :: wprime(*)

INTEGER  :: i, j
REAL     :: temp

j = kprime
DO i = 1, kprime / 2 + MOD(kprime, 2)
  temp = wprime(i)
  wprime(i) = REAL(DBLE(sum) - DBLE(wprime(j)))
  wprime(j) = REAL(DBLE(sum) - DBLE(temp))
  j = j - 1
END DO

RETURN
END SUBROUTINE cmplmt



SUBROUTINE gener (w, n, wprime, nprime, INDEX, r)

!     ALGORITHM AS 304.6 APPL.STATIST. (1996), VOL.45, NO.3

!     Computes an array of sums of the various R-combinations of
!     the elements of W

IMPLICIT NONE
INTEGER, INTENT(IN)      :: n, nprime, r
INTEGER, INTENT(IN OUT)  :: INDEX(r)
REAL, INTENT(IN)         :: w(n)
REAL, INTENT(OUT)        :: wprime(nprime)
INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(14, 60)

INTEGER    :: i, j
REAL (dp)  :: sum
LOGICAL    :: init

init = .true.

DO i = 1, nprime
  CALL next(INDEX, r, n, init)
  sum = 0.0D0
  DO j = 1, r
    sum = sum + DBLE(w(INDEX(j)))
  END DO
  wprime(i) = REAL(sum)
END DO

RETURN
END SUBROUTINE gener



SUBROUTINE next (rcombo, r, n, init)

!     ALGORITHM AS 304.7 APPL.STATIST. (1996), VOL.45, NO.3

!     Accepts some R-combination of the first N integers and then
!     computes the next R-combination in the lexicographic
!     ordering of the N! / (R! * (N - R)!) such R-combinations.
!     Returns the first R-combination if the initialization
!     indicator is .true. and then resets the indicator.

IMPLICIT NONE
INTEGER, INTENT(IN)      :: r, n
INTEGER, INTENT(IN OUT)  :: rcombo(r)
LOGICAL, INTENT(IN OUT)  :: init

INTEGER :: i, j, d

IF (init) THEN
  DO i = 1, r
    rcombo(i) = i
  END DO
  init = .false.
ELSE
  d = n - r
  j = r
  
!        The counter J is not prevented from going out of bounds
!        which will happen if there is no next R-combination
  
  20 IF (rcombo(j) < d + j) GO TO 30
  j = j - 1
  GO TO 20

  30 rcombo(j) = rcombo(j) + 1
  DO i = j + 1, r
    rcombo(i) = rcombo(i - 1) + 1
  END DO
END IF

RETURN
END SUBROUTINE next



!        General purpose subroutines

SUBROUTINE sort (x, n)

!     ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3

!     Sorts the N values stored in array X in ascending order

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(n)

INTEGER  :: i, j, incr
REAL     :: temp

incr = 1

!        Loop : calculate the increment

10 incr = 3 * incr + 1
IF (incr <= n) GO TO 10

!        Loop : Shell-Metzner sort

20 incr = incr / 3
i = incr + 1
30 IF (i > n) GO TO 60
temp = x(i)
j = i
40 IF (x(j - incr) < temp) GO TO 50
x(j) = x(j - incr)
j = j - incr
IF (j > incr) GO TO 40
50 x(j) = temp
i = i + 1
GO TO 30
60 IF (incr > 1) GO TO 20

RETURN
END SUBROUTINE sort



FUNCTION comb (n, k) RESULT(fn_val)

!        ALGORITHM AS 304.9 APPL.STATIST. (1996), VOL.45, NO.3

!        Returns the number of combinations of N things taken K at a
!        time or 0 if the parameters are incompatible

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, k
INTEGER              :: fn_val

INTEGER  :: m, i
REAL     :: numer, denom

IF (k < 0 .OR. k > n) THEN
  fn_val = 0
ELSE
  m = n - k
  numer = 1.0
  DO i = n, 1 + MAX(k, m), -1
    numer = numer * REAL(i)
  END DO
  denom = fact(MIN(k, m))
  fn_val = nint(numer / denom)
END IF

RETURN
END FUNCTION comb



FUNCTION fact (n) RESULT(fn_val)

!     ALGORITHM AS 304.10 APPL.STATIST. (1996), VOL.45, NO.3

!     Returns the factorial of N or 0.0 if the parameter is negative.

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL                :: fn_val

INTEGER  :: i

IF (n < 0) THEN
  fn_val = 0.0
ELSE
  fn_val = 1.0
  DO i = 2, n
    fn_val = fn_val * REAL(i)
  END DO
END IF

RETURN
END FUNCTION fact

END MODULE Fisher_Test
