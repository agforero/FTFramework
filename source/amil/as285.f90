MODULE Mult_Normal
! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-15  Time: 09:34:54

!     ALGORITHM AS285 APPL. STATIST. (1993) VOL.42, NO.3

!  FINDS THE PROBABILITY THAT A NORMALLY DISTRIBUTED RANDOM
!  N-VECTOR WITH MEAN 0 AND COVARIANCE COVAR FALLS
!  IN AREA ENCLOSED BY THE EXTERNAL USER-DEFINED FUNCTION F.

IMPLICIT NONE
INTEGER, SAVE  :: ix, iy, iz


CONTAINS


SUBROUTINE mulnor(covar, n, f, sdev, rhigh, init, iseed, t, work, prob,  &
    iter, ifault)

IMPLICIT NONE
REAL, INTENT(IN OUT)     :: covar(*)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: sdev
REAL, INTENT(IN)         :: rhigh
INTEGER, INTENT(IN)      :: init
INTEGER, INTENT(IN OUT)  :: iseed
REAL, INTENT(IN OUT)     :: t(*)
REAL, INTENT(OUT)        :: work(n,*)
REAL, INTENT(OUT)        :: prob
INTEGER, INTENT(OUT)     :: iter
INTEGER, INTENT(OUT)     :: ifault

INTERFACE
  FUNCTION f(n, v) RESULT(fn_val)
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER  :: init1, k, ii, kk, q, nn, nullty, i, j
REAL     :: p, var, cox

REAL, PARAMETER     :: zero = 0.0
INTEGER, PARAMETER  :: maxitr = 10000

ifault = 0
init1 = init
IF (init < 25) init1 = 25
IF (init > 1000) init1 = 1000
cox = 1.0 + 2.0/init1

!  FIND CHOLESKY DECOMPOSITION OF COVARIANCE MATRIX

nn = n*(n+1)/2
CALL chol(covar, n, nn, t, nullty, ifault)
q = n - nullty
ifault = ifault*ifault
IF (ifault /= 0) RETURN

!  TRANSPOSE CHOLESKY FACTOR, OMITTING COLUMNS OF ZEROES.
!  RESULTING VECTOR IS PACKED FORM OF T'J,
!  WRITTEN ROWWISE BEGINNING WITH LAST ROW.

DO  i = 1,n
  work(1:i,i) = zero
END DO
ii = 0
DO  i = 1,n
  DO  j = 1,i
    ii = ii + 1
    work(i,j) = t(ii)
  END DO
END DO
j = 0
DO  i = n,1,-1
  IF (j == nullty) EXIT
  IF (work(i,i) == zero) THEN
    j = j+1
    DO  k = i+1,n
      DO  kk = i, n-j
        work(k,kk) = work(k,kk+1)
      END DO
    END DO
  END IF
END DO

k = 0
DO  j = n,1,-1
  DO  i = 1, MIN(j,q)
    k = k + 1
    t(k) = work(j,i)
  END DO
END DO

!  TAKE PILOT SAMPLE; ESTIMATE VARIANCE OF PROB FROM ORTHP.

prob = zero
var = zero
DO  k=1,init1
  CALL orthp(f, rhigh, n, q, t, work, p, iseed, ifault)
  IF (ifault /= 0) RETURN
  prob = prob + p
  var = var + p*p
END DO
prob = prob/init1
var = var/init1 - prob*prob

iter = INT(cox*var/(sdev*sdev)) + 1
IF (iter > maxitr) ifault = 5

DO  k=init1+1, iter
  CALL orthp(f, rhigh, n, q, t, work, p, iseed, ifault)
  IF (ifault /= 0) RETURN
  prob = prob + (p - prob)/k
END DO

RETURN
END SUBROUTINE mulnor



SUBROUTINE chol(a, n, nn, u, nullty, ifault)

IMPLICIT NONE
INTEGER, INTENT(IN)   :: nn
REAL, INTENT(IN)      :: a(nn)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: u(nn)
INTEGER, INTENT(OUT)  :: nullty
INTEGER, INTENT(OUT)  :: ifault

REAL     :: eta2, x, w
INTEGER  :: j, k, ii, icol, l, kk, m, irow, i

REAL, PARAMETER  :: eta = 1.0E-05, zero = 0.0

! zabs(x) = ABS(x)
! zsqrt(x) = SQRT(x)

ifault = 1
IF (n <= 0) RETURN
ifault = 3
IF (nn /= n*(n+1)/2) RETURN
ifault = 2
nullty = 0
j = 1
k = 0
eta2 = eta*eta
ii = 0
DO  icol = 1,n
  ii = ii + icol
  x = eta2 * a(ii)
  l = 0
  kk = 0
  DO  irow = 1, icol
    kk = kk + irow
    k = k+1
    w = a(k)
    m = j
    DO  i = 1, irow
      l = l + 1
      IF (i == irow) EXIT
      w = w - u(l) * u(m)
      m = m + 1
    END DO
    IF (irow == icol) EXIT
    IF (u(l) == zero) GO TO 30
    u(k) = w / u(l)
    CYCLE
    30 IF (w * w > abs(x * a(kk))) RETURN
    u(k) = zero
  END DO
  IF (abs(w) <= abs(eta * a(k))) GO TO 60
  IF (w < zero) RETURN
  u(k) = sqrt(w)
  GO TO 70
  60 u(k) = zero
  nullty = nullty + 1
  70 j = j + icol
END DO
ifault = 0
RETURN
END SUBROUTINE chol



FUNCTION chiprb(x, df) RESULT(fn_val)

!  FINDS THE PROBABILITY THAT A CHI-SQUARED RANDOM VARIABLE
!  WITH DF DEGREES OF FREEDOM EXCEEDS X.
!  ADAPTED FROM HILL AND PIKE (1967) ALG. 299, CHI-SQUARED INTEGRAL,
!  COMM. ACM, VOL.10, 243-244.

IMPLICIT NONE
REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: df
REAL                 :: fn_val

REAL     :: a, sqrta, y, c, e, z, s
INTEGER  :: odd, n2, i
!     REAL ERF

REAL, PARAMETER  :: rtpirp = 0.564189583547756
REAL, PARAMETER  :: big = 88.03, small = -85.2, zero = 0.0, one = 1.0, &
   half = 0.5

!  NOTE: IF THE FUNCTION "ERF" IS AVAILABLE, THE ABOVE STATEMENT
!  FUNCTION MAY BE REPLACED BY
!     ZERFP1(A) = ERF(A) + ONE

s = zero
IF (x > big) GO TO 50
s = one
IF (x < small) GO TO 50

a = half*x
y = EXP(-a)
odd = MOD(df,2)
IF (odd == 0) THEN
  s = y
  e = one
  z = zero
ELSE
  sqrta = SQRT(a)
  s = zerfp1(-sqrta)
  e = rtpirp/sqrta
  z = -half
END IF
IF (df < 3) GO TO 50

n2 = df/2 - 1 + odd
c = zero
DO  i = 1,n2
  z = z + 1
  e = e*(a/z)
  c = c + e
END DO
s = s + c*y

50 fn_val = s
RETURN
END FUNCTION chiprb



SUBROUTINE orthp(f, rhigh, n, q, t, z, prob, iseed, ifault)

!  FINDS THE PSEUDO-MONTE CARLO ESTIMATE OF PROB USING AN
!  ORTHONORMAL SYSTEM OF VECTORS.

IMPLICIT NONE
REAL, INTENT(IN)         :: rhigh
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: q
REAL, INTENT(IN)         :: t(*)
REAL, INTENT(IN OUT)     :: z(n,*)
REAL, INTENT(OUT)        :: prob
INTEGER, INTENT(IN OUT)  :: iseed
INTEGER, INTENT(OUT)     :: ifault

INTERFACE
  FUNCTION f(n, v) RESULT(fn_val)
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER  :: i, j, k, jj, ii, qp1, qp2
REAL     :: s, a, b, r, fab

! EXTERNAL f
REAL, PARAMETER  :: zero = 0.0, rt2rcp = 0.707106781

prob = zero
qp1 = q+1
qp2 = q+2

!  PUT RANDOMLY ORIENTED ORTHONORMAL SYSTEM OF Q-VECTORS IN Z.

CALL rnortm(n,q, qp1, 1, 1, z(1,qp1), 1, z(1,qp1), iseed, z,  &
    z(1,qp1), fab, z(1,qp2), z(1,qp2), ifault)

!  REPLACE EACH COLUMN OF Z BY TRANS(CHOLESKY FACTOR)*J*COLUMN OF Z.

DO  i=1,q
  jj = 0
  DO  j = n,1,-1
    s = zero
    DO  k=1, MIN(j,q)
      jj = jj + 1
      s = s + t(jj)*z(k,i)
    END DO
    z(j,i) = s
  END DO
END DO

!  ESTIMATE THE PROBABILITY THAT N(0,COVAR) IS WITHIN THE BOUNDARY GIVEN BY F.

DO  i = 1,q
  a = zero
  b = rhigh
  CALL zerovc (f, r, a, b, n, z(1,i), z(1,qp2), ifault)
  prob = prob - chiprb(r*r, q)
END DO

DO  i = 1,q
  a = zero
  b = rhigh
  DO  k = 1,n
    z(k,qp1) = -z(k,i)
  END DO
  CALL zerovc (f, r, a, b, n, z(1,qp1), z(1,qp2), ifault)
  prob = prob - chiprb(r*r, q)
END DO

DO  i = 1,q-1
  DO  j = i+1,q
    DO  jj= -1,1,2
      DO  ii = -1,1,2
        DO  k=1,n
          z(k,qp1) = rt2rcp*(z(k,i)*ii + z(k,j)*jj)
        END DO
        a = zero
        b = rhigh
        CALL zerovc (f,r,a,b,n,z(1,qp1),z(1,qp2),ifault)
        prob = prob - chiprb(r*r,q)
      END DO
    END DO
  END DO
END DO
prob = prob/(2.0*q*q) + 1.0

RETURN
END SUBROUTINE orthp



SUBROUTINE rnortm(lda, n, np1, ibconf, nb, b, ndbi, dbi, iseed, a,  &
    chisq, fab, u, w, ifault)

!      ALGORITHM AS 127  APPL. STATIST.  (1978) VOL.27, NO.2

!      RNORTM GENERATES ORTHOGONAL MATRICES A FROM A DISTRIBUTION WITH
!      DENSITY FUNCTION DEFINED ON THE GROUP OF ORTHOGONAL MATRICES.
!      WHEN THE INPUT PARAMETER IBCONF = 1 THE MATRICES ARE
!      GENERATED FROM THE INVARIANT HAAR MEASURE.

IMPLICIT NONE
INTEGER, INTENT(IN)      :: lda
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: np1
INTEGER, INTENT(IN)      :: ibconf
INTEGER, INTENT(IN)      :: nb
REAL, INTENT(IN OUT)     :: b(nb)
INTEGER, INTENT(IN)      :: ndbi
REAL, INTENT(IN OUT)     :: dbi(ndbi)
INTEGER, INTENT(IN OUT)  :: iseed
REAL, INTENT(OUT)        :: a(lda,n)
REAL, INTENT(IN OUT)     :: chisq(n)
REAL, INTENT(OUT)        :: fab
REAL, INTENT(OUT)        :: u(np1)
REAL, INTENT(OUT)        :: w(n)
INTEGER, INTENT(OUT)     :: ifault

REAL    :: ul, wl, wil, wwil, t2, hv, hw
INTEGER :: i, j, l, np1l
INTEGER :: lp1

REAL, PARAMETER  :: zero = 0.0, one = 1.0

!      STATEMENT FUNCTIONS FOR DOUBLE PRECISION
!      DOUBLE PRECISION SQRT, SIGN, ABS
!      SQRT(UL) = DSQRT(UL)
!      SIGN(UL, WL) = DSIGN(UL, WL)
!      ABS(UL) = DABS(UL)

!      SECTION 1 INITIALIZATION

ifault = 0
GO TO 4

!      ERROR RETURNS

2 ifault = 2
RETURN

!      INITIALIZE DENSITY, H(STORED IN A), AND IBSUB

4 fab = one
DO  i=1,n
  a(i,1:n) = zero
  a(i,i) = one
END DO

!      CONSTRUCT  A  ONE COLUMN AT A TIME

DO  l=1,n
  np1l = np1 - l
  
!      SECTION 2
!      CONSTRUCT U(L...N) UNIFORMLY DISTRIBUTED ON THE SURFACE OF THE
!      (N+1-L)-SPHERE OF RADIUS UL = SQRT(CHISQ(L)).
  
  CALL normal(u(l), np1l+1, iseed, chisq(l))
  ul = SQRT(chisq(l))
  IF (ul == zero) GO TO 2
  IF (l == n) GO TO 301
  
!      SECTION 3
!      CALCULATE  W = B*U  AND STORE ITS LENGTH IN WL.
!      SIGN OF WL PREVENTS LOSS OF SIGNIFICANCE IN WWIL.
  
!      B IS IDENTITY
  
  wl = SIGN(ul,-u(l))
  DO  j=l,n
    w(j) = u(j)
  END DO
  
!      SECTION 4
!      PROJECT  W  ONTO ORTHOGOANL COMPLEMENT OF FIRST L-1
!      COLUMNS OF A, NORMALIZE, AND STORE IN L COLUMN OF A.
!      A(L) = H(L) * W * WIL
!      CALCULATE PROJECTION MATRIX H(L) FOR ORTHOGOANL
!      COMPLEMENT OF FIRST L COLUMNS OF A AND STORE IN LAST
!      N-L COLUMNS OF A.
  
  lp1 = l+1
  wil = one/wl
  wwil = one/(wl - w(l))
  DO  i=1,n
    t2 = zero
    DO  j=l,n
      t2 = t2 + a(i,j) * w(j)
    END DO
    hv = t2 * wil
    hw = (hv - a(i,l)) * wwil
    a(i,l) = hv
    DO  j = lp1, n
      a(i,j) = a(i,j) - hw*w(j)
    END DO
  END DO
  
!      END OF L COLUMN OF MATRIX A
  
!      SECTION 5
!      WHEN L = N ONLY SIGN NEEDS TO BE CHOSEN
!      REMEMBER MATRIX H(N) IS STORED IN N COLUMN OF A.
  
  301 IF (u(l) > zero) CYCLE
  DO  i=1,n
    a(i,l) = -a(i,l)
  END DO
  
!      FAB = ABS(B(LAST) * DBI(N)*FAB)
  
END DO
fab = ABS(fab)

RETURN
END SUBROUTINE rnortm



SUBROUTINE zerovc (f, r, r0, r1, n, v, y, ifault)

!  CALCULATES THE SOLUTION R TO F(R*V) = 0.

IMPLICIT NONE
REAL, INTENT(OUT)     :: r
REAL, INTENT(IN OUT)  :: r0
REAL, INTENT(IN OUT)  :: r1
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: v(n)
REAL, INTENT(OUT)     :: y(n)
INTEGER, INTENT(OUT)  :: ifault

INTERFACE
  FUNCTION f(n, v) RESULT(fn_val)
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER  :: j
REAL     :: f0, f1, fr
REAL, PARAMETER     :: half = 0.5, eta = 1.0E-05
INTEGER, PARAMETER  :: niter = 500

y(1:n) = r0*v(1:n)
f0 = f(n,y)
y(1:n) = r1*v(1:n)
f1 = f(n,y)
IF (SIGN(half,f1) == SIGN(half,f0)) THEN
  ifault = 3
  RETURN
END IF
DO  j = 1,niter
  r = r1 - f1*(r1-r0)/(f1-f0)
  y(1:n) = r*v(1:n)
  fr = f(n, y)
  IF (abs(r-r1) < eta) GO TO 50
  IF (SIGN(half,fr) /= SIGN(half,f1)) THEN
    r0=r1
    f0=f1
  ELSE
    f0 = f0 * half
  END IF
  r1 = r
  f1 = fr
END DO
ifault = 6

50 RETURN
END SUBROUTINE zerovc


FUNCTION zerfp1(a) RESULT(fn_val)
IMPLICIT NONE
REAL, INTENT(IN)  :: a
REAL              :: fn_val

REAL, PARAMETER  :: sqrt2 = 1.41421356237309

fn_val = 2.0*alnorm(sqrt2*a, .false.)
RETURN
END FUNCTION zerfp1



SUBROUTINE normal(u, np1, iseed, chisq)

!  GENERATES AN (NP1-1)-VECTOR OF INDEPENDENT NORMAL VARIATES IN U.
!  ON OUTPUT, CHISQ CONTAINS THE SQUARED LENGTH OF THE VECTOR.

IMPLICIT NONE

INTEGER, INTENT(IN)      :: np1
REAL, INTENT(OUT)        :: u(np1)
INTEGER, INTENT(IN OUT)  :: iseed
REAL, INTENT(OUT)        :: chisq

INTEGER  :: i, ifault

! COMMON /rand/ ix,iy,iz
REAL, PARAMETER  :: zero = 0.0

chisq = zero
DO  i = 1, np1-1
  CALL ppnd(random(), ifault, u(i))
  chisq = chisq + u(i)*u(i)
END DO
RETURN
END SUBROUTINE normal



SUBROUTINE ppnd(p, ifault, fn_val)

!  ALGORITHM AS 111  APPL. STATIST. (1977), VOL.26, NO.1

!  PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA OF P
!  REAL VERSION FOR EPS = 2 **(-31)
!  THE HASH SUMS ARE THE SUMS OF THE MODULI OF THE COEFFICIENTS
!  THEY HAVE NO INHERENT MEANINGS BUT ARE INCLUDED FOR USE IN
!  CHECKING TRANSCRIPTIONS
!  STANDARD FUNCTIONS ABS, ALOG AND SQRT ARE USED

IMPLICIT NONE
REAL, INTENT(IN)      :: p
INTEGER, INTENT(OUT)  :: ifault
REAL, INTENT(OUT)     :: fn_val

REAL :: q, r

REAL, PARAMETER  :: zero = 0.0E0, half = 0.5E0, one = 1.0E0
REAL, PARAMETER  :: split = 0.42E0
REAL, PARAMETER  :: a0 = 2.50662823884E0
REAL, PARAMETER  :: a1 = -18.61500062529E0
REAL, PARAMETER  :: a2 = 41.39119773534E0
REAL, PARAMETER  :: a3 = -25.44106049637E0
REAL, PARAMETER  :: b1 = -8.47351093090E0
REAL, PARAMETER  :: b2 = 23.08336743743E0
REAL, PARAMETER  :: b3 = -21.06224101826E0
REAL, PARAMETER  :: b4 = 3.13082909833E0
REAL, PARAMETER  :: c0 = -2.78718931138E0
REAL, PARAMETER  :: c1 = -2.29796479134E0
REAL, PARAMETER  :: c2 = 4.85014127135E0
REAL, PARAMETER  :: c3 = 2.32121276858E0
REAL, PARAMETER  :: d1 = 3.54388924762E0
REAL, PARAMETER  :: d2 = 1.63706781897E0

ifault = 0
q = p - half
IF (ABS(q) > split) GO TO 1
r = q*q
fn_val = q * (((a3*r + a2)*r + a1) * r + a0) /  &
    ((((b4*r + b3)*r + b2) * r + b1) * r + one)
RETURN

1 r = p
IF (q > zero)r = one - p
IF (r <= zero) GO TO 2
r = SQRT(-LOG(r))
fn_val = (((c3 * r + c2) * r + c1) * r + c0)/ ((d2*r + d1) * r + one)
IF (q < zero) fn_val = -fn_val
RETURN

2 ifault = 1
fn_val = zero
RETURN
END SUBROUTINE ppnd



FUNCTION random() RESULT(fn_val)

!  ALGORITHM AS 183  APPL. STATIST. (1982) VOL. 31, NO. 2

!  RETURNS A PSEUDO-RANDOM NUMBER RECTANGULARLY DISTRIBUTED BETWEEN 0 AND 1.

!  IX, IY, AND IZ SHOULD BE SET TO INTEGER VALUES BETWEEN
!  1 AND 30000 BEFORE FIRST ENTRY

!  INTEGER ARITHMETIC UP TO 30323 IS REQUIRED

REAL  :: fn_val

! INTEGER :: ix, iy, iz
! COMMON /rand/ ix,iy,iz

ix = 171 * MOD(ix, 177) - 2 * (ix / 177)
iy = 172 * MOD(iy, 176) - 35 * (iy / 176)
iz = 170 * MOD(iz, 178) - 63 * (iz / 178)

IF (ix < 0) ix = ix + 30269
IF (iy < 0) iy = iy + 30307
IF (iz < 0) iz = iz + 30323

fn_val = MOD(REAL(ix)/30269.0 + REAL(iy)/30307.0 + REAL(iz)/30323.0, 1.0)
RETURN
END FUNCTION random



FUNCTION alnorm(x, upper) RESULT(fn_val)

!  ALGORITHM AS 66  APPL. STATIST. (1973) VOL. 22, NO. 3

!  EVALUATES THE TAIL AREA OF THE STANDARDIZED NORMAL CURVE
!  FROM X TO INFINITY IF UPPER IS .TRUE. OR
!  FROM MINUS INFINITY TO X IF UPPER IS .FALSE.


REAL, INTENT(IN)     :: x
LOGICAL, INTENT(IN)  :: upper
REAL                 :: fn_val

REAL     :: z, y
LOGICAL  :: up

!  LTONE AND UTZERO MUST BE SET TO SUIT THE PARTICULAR COMPUTER
!  (SEE INTRODUCTORY TEXT)

REAL, PARAMETER  :: ltone = 7.0, utzero = 18.66
REAL, PARAMETER  :: zero = 0.0, half = 0.5, one = 1.0, con = 1.28

up = upper
z = x
IF (z >= zero) GO TO 10
up = .NOT. up
z = -z
10    IF (z <= ltone .OR. up .AND. z <= utzero) GO TO 20
fn_val = zero
GO TO 40
20 y = half *z * z
IF (z > con) GO TO 30

fn_val = half - z*(0.398942280444 - 0.399903438504 * y /  &
    (y + 5.75885480458 - 29.8213557808 / (y + 2.62433121679 + 48.6959930692 /  &
    (y + 5.92885724438))))
GO TO 40

30 fn_val = 0.398942280385 * EXP(-y) / (z - 3.8052E-08 + 1.00000615302 /  &
    (z + 3.98064794E-4 + 1.98615381364 / (z - 0.151679116635 + 5.29330324926 /  &
    (z + 4.8385912808 - 15.1508972451 /  &
    (z + 0.742380924027 + 30.789933034 / (z + 3.99019417011))))))

40 IF (.NOT. up) fn_val = one - fn_val
RETURN
END FUNCTION alnorm

END MODULE Mult_Normal



PROGRAM Driver285
!  DRIVING PROGRAM
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-15  Time: 12:20:56

USE Mult_Normal
IMPLICIT NONE

REAL     :: covar(1275), work(50,52), t(1275)
INTEGER  :: iseed, n, ifault, j, iter, init, nn, nullty
REAL     :: sdev, prob, rhigh
! COMMON /rand/ ix,iy,iz

REAL, PARAMETER  :: eta = 1.0E-5

INTERFACE
  FUNCTION f1(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f1

  FUNCTION f2(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f2

  FUNCTION f3(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f3

  FUNCTION f4(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f4

  FUNCTION f5(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f5

  FUNCTION f6(n,v) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: v(:)
    REAL                 :: fn_val
  END FUNCTION f6
END INTERFACE

! READ IN COVARIANCE AND INITIALIZE RANDOM NORMAL GENERATOR.

OPEN(1, file='FORT.1', status='old')
READ (1, *) n, sdev, rhigh, init, ix, iy, iz
READ (1, *) (covar(j),j=1,n*(n+1)/2)

!  FIND SMALLEST NONZERO EIGENVALUE OF COVARIANCE MATRIX

nn = n*(n+1)/2
CALL chol(covar, n, nn, t, nullty, ifault)
CALL tdiag(n, eta, covar, work(1,1), work(1,2), work(1,3), ifault)
CALL lrvt(n, eta, work(1,1), work(1,2), work(1,3), ifault)
rhigh = rhigh/SQRT(work(nullty+1,1))

! Run through the supplied functions F1 ... F6

DO n = 2, 5
  CALL mulnor(covar, n, f1, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
  WRITE(*, *) 'Function F1, with n = ', n
  WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
  WRITE(*, *) 'ITER = ',iter
END DO

n = 2
CALL mulnor(covar, n, f2, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
WRITE(*, *)
WRITE(*, *) 'Function F2'
WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
WRITE(*, *) 'ITER = ',iter

n = 4
CALL mulnor(covar, n, f3, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
WRITE(*, *)
WRITE(*, *) 'Function F3'
WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
WRITE(*, *) 'ITER = ',iter

n = 3
CALL mulnor(covar, n, f4, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
WRITE(*, *)
WRITE(*, *) 'Function F4'
WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
WRITE(*, *) 'ITER = ',iter

n = 2
CALL mulnor(covar, n, f5, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
WRITE(*, *)
WRITE(*, *) 'Function F5'
WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
WRITE(*, *) 'ITER = ',iter

n = 2
CALL mulnor(covar, n, f6, sdev, rhigh, init, iseed, t, work, prob, iter, ifault)
WRITE(*, *)
WRITE(*, *) 'Function F6'
WRITE(*, *) 'ESTIMATED PROBABILITY IS', prob, 'WITH IFAULT = ', ifault
WRITE(*, *) 'ITER = ',iter

STOP


CONTAINS


SUBROUTINE tdiag(n, tol, a, d, e, z, ifault)

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: tol
REAL, INTENT(IN)      :: a(*)
REAL, INTENT(OUT)     :: d(n)
REAL, INTENT(OUT)     :: e(n)
REAL, INTENT(OUT)     :: z(n,n)
INTEGER, INTENT(OUT)  :: ifault

REAL     :: f, g, h, hh
INTEGER  :: i, j, k, l, i1, j1

REAL, PARAMETER  :: zero = 0.0, one = 1.0

ifault = 1
IF (n <= 1) RETURN
ifault = 0
k = 0
DO  i=1,n
  DO  j=1,i
    k = k+1
    z(i,j) = a(k)
  END DO
END DO
i=n
DO  i1=2,n
  l = i-2
  f = z(i,i-1)
  g=zero
  IF (l < 1) GO TO 25
  DO  k=1,l
    g=g + z(i,k)**2
  END DO
  25 h=g + f*f
  
  IF (g > tol) GO TO 30
  e(i)=f
  d(i)=zero
  GO TO 65

  30 l=l+1
  g=SQRT(h)
  IF (f >= zero) g=-g
  e(i) = g
  h = h-f*g
  z(i,i-1)=f-g
  f=zero
  DO  j=1,l
    z(j,i) = z(i,j)/h
    g = zero
    
    DO  k=1,j
      g = g + z(j,k) *z(i,k)
    END DO
    IF (j >= l)GO TO 47
    j1  = j+1
    DO  k=j1,l
      g=g + z(k,j)*z(i,k)
    END DO
    
    47 e(j)=g/h
    f = f+g*z(j,i)
  END DO
  
  hh = f/(h+h)
  
  DO  j=1,l
    f=z(i,j)
    g=e(j)-hh*f
    e(j)=g
    DO  k=1,j
      z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
    END DO
  END DO
  d(i)=h
  65    i=i-1
END DO
d(1) = zero
e(1) = zero

DO  i=1,n
  l=i-1
  IF (d(i) == zero .OR. l == 0) GO TO 100
  DO  j=1,l
    g=zero
    DO  k=1,l
      g=g + z(i,k)*z(k,i)
    END DO
    DO  k=1,l
      z(k,j) = z(k,j) - g*z(k,i)
    END DO
  END DO
  100   d(i) = z(i,i)
  z(i,i) = one
  IF (l == 0) CYCLE
  DO  j=1,l
    z(i,j) = zero
    z(j,i) = zero
  END DO
END DO
RETURN
END SUBROUTINE tdiag



SUBROUTINE lrvt(n, precis, d, e, z, ifault)

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: precis
REAL, INTENT(IN OUT)  :: d(n)
REAL, INTENT(OUT)     :: e(n)
REAL, INTENT(IN OUT)  :: z(n,n)
INTEGER, INTENT(OUT)  :: ifault

REAL     :: b, c, f, g, h, p, pr, r, s
INTEGER  :: i, j, k, l, n1, jj, m1, m, i1

INTEGER, PARAMETER  :: mits = 30
REAL, PARAMETER     :: zero = 0.0, one = 1.0, two = 2.0

ifault = 2
IF (n <= 1) RETURN
ifault = 1
n1 = n-1
DO  i=2,n
  e(i-1) = e(i)
END DO
e(n) = zero
b = zero
f = zero
DO  l=1,n
  jj=0
  h = precis * (ABS(d(l)) + ABS(e(l)))
  IF (b < h) b = h
  
  DO  m1=l,n
    m=m1
    IF (ABS(e(m)) <= b) EXIT
  END DO
  IF (m == l) GO TO 85
  40 IF (jj == mits) RETURN
  jj = jj+1
  
  p = (d(l+1) - d(l))/(two*e(l))
  r = SQRT(p*p + one)
  pr = p + r
  IF (p < zero) pr = p - r
  h = d(l) - e(l)/pr
  DO  i=l,n
    d(i) = d(i)-h
  END DO
  f=f+h
  
  p=d(m)
  c=one
  s=zero
  m1=m-1
  i=m
  DO  i1=l,m1
    j=i
    i=i-1
    g=c*e(i)
    h=c*p
    IF(ABS(p) >= ABS(e(i))) GO TO 60
    c = p/e(i)
    r = SQRT(c*c+one)
    e(j) = s*e(i)*r
    s = one/r
    c = c/r
    GO TO 70

    60 c=e(i)/p
    r=SQRT(c*c+one)
    e(j) = s*p*r
    s=c/r
    c=one/r
    70 p=c*d(i)-s*g
    d(j)=h + s*(c*g + s*d(i))
    
    DO  k=1,n
      h=z(k,j)
      z(k,j) = s*z(k,i) + c*h
      z(k,i) = c*z(k,i) - s*h
    END DO
  END DO
  e(l)=s*p
  d(l)=c*p
  IF (ABS(e(l)) > b) GO TO 40
  85 d(l) = d(l) + f
END DO

DO  i=1,n1
  k=i
  p=d(i)
  i1=i+1
  DO  j=i1,n
    IF (d(j) >= p) CYCLE
    k=j
    p=d(j)
  END DO
  IF (k == i) CYCLE
  d(k) = d(i)
  d(i) = p
  DO  j=1,n
    p=z(j,i)
    z(j,i)=z(j,k)
    z(j,k)=p
  END DO
END DO
ifault=0
RETURN
END SUBROUTINE lrvt

END PROGRAM Driver285




FUNCTION f1(n,v) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-15  Time: 13:39:21

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

REAL, PARAMETER  :: cubsz = 3.0
INTEGER  :: i

!  BOUNDARY OF REGION IS AN N-DIMENSIONAL CUBE.
!  CUBSZ IS THE DISTANCE FROM THE ORIGIN TO THE SIDE OF THE CUBE

fn_val = 0.0
DO  i = 1,n
  fn_val = MAX(fn_val, ABS(v(i)))
END DO
fn_val = fn_val - cubsz
RETURN
END FUNCTION f1



FUNCTION f2(n,v) RESULT(fn_val)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

REAL     :: b(10)
INTEGER  :: i, j, k

!  THE BOUNDARY IS AN ELLIPSOID OF THE FORM X'AX + BX + C = 0
!  THE UPPER TRIANGLE OF A IS STORED COLUMN-WISE IN PACKED FORM.
!  IN THIS EXAMPLE, A = INVERSE (1  .3), B = 0, C = -1.
!                               (.3  1)

REAL, PARAMETER  :: a(3) = (/ 1.0989011, -.32967033, 1.0989011 /)
REAL, PARAMETER  :: c = -1.0

b(1:2) = 0.0
fn_val = c
k = 0
DO  i = 1,n
  k = k + 1
  fn_val = fn_val + (v(i)**2)*a(k)
  DO  j = i+1,n
    k = k + 1
    fn_val = fn_val + 2.0*v(i)*v(j)*a(k)
  END DO
END DO
fn_val = fn_val + DOT_PRODUCT(b(1:n), v(1:n))
RETURN
END FUNCTION f2



FUNCTION f3(n,v) RESULT(fn_val)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

INTEGER  :: i
REAL     :: temp

!  RECTANGLE IS GIVEN BY A < V < B FOR EACH COORDINATE.
!  THIS SUBROUTINE MAY ALSO BE USED TO CALCULATE THE CDF BY
!  MAKING THE APPROPRIATE SIDES OF THE RECTANGLE SUFFICIENTLY LARGE.
!  THIS EXAMPLE IS THE FOUR-DIMENSIONAL RECTANGLE DEFINED
!  BY  -INFINITY (ALMOST) < V(I) < 2 FOR EACH I.

REAL, PARAMETER  :: a(4) = (/ -20.0, -20.0, -20.0, -20.0 /)
REAL, PARAMETER  :: b(4) = (/ 2.0, 2.0, 2.0, 2.0 /)

fn_val = -1.0
DO  i=1,n
  IF (v(i) > 0.0) THEN
    temp = v(i) - b(i)
  ELSE
    temp = a(i) - v(i)
  END IF
  fn_val  = MAX(fn_val, temp)
END DO
RETURN
END FUNCTION f3



FUNCTION f4(n,v) RESULT(fn_val)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

REAL, PARAMETER  :: radsq = 4.0
INTEGER  :: i

!  BOUNDARY IS SPHERE OF RADIUS SQRT(RADSQ).

fn_val = 0.0
DO  i = 1,n
  fn_val = fn_val + v(i)*v(i)
END DO
fn_val = fn_val - radsq
RETURN
END FUNCTION f4



FUNCTION f5(n,v) RESULT(fn_val)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

INTEGER  :: i, j
REAL     :: temp

!  A CONVEX POLYTOPE WITH P FACES IS DEFINED BY THE INEQUALITIES
!  A(1,I)*V(1) + A(2,I)*V(2) +...+ A(N,I)*V(N) <= 1
!  FOR I=1,...,P.

!  THE FIGURE IN THIS DEMONSTRATION PROGRAM IS THE
!  REGULAR PENTAGON WITH DISTANCE FROM CENTER TO POINT = 1.

INTEGER, PARAMETER  :: p = 5
REAL, PARAMETER     :: a(2,5) = RESHAPE( (/ 0.726542528, -1.0,  &
    1.175570504, 0.381966011, 0.0, 1.236067978, -0.726542528, -1.0,  &
    -1.175570504, 0.381966011 /), (/ 2, 5 /) )

fn_val = -1.0
DO  j=1,n
  fn_val = fn_val + a(j,1)*v(j)
END DO
DO  i=2,p
  temp = -1.0
  DO  j=1,n
    temp = temp + a(j,i)*v(j)
  END DO
  fn_val = MAX(fn_val, temp)
END DO
RETURN
END FUNCTION f5


FUNCTION f6(n,v) RESULT(fn_val)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: v(:)
REAL                 :: fn_val

INTEGER  :: i, j, INDEX
REAL     :: x(2), tangnt

!  THE FIGURE IN THIS DEMONSTRATION PROGRAM IS THE
!  REGULAR PENTAGRAM WITH DISTANCE FROM CENTER TO POINT = 2.618.

INTEGER, PARAMETER  :: p = 5
REAL, PARAMETER     :: tang(5) = (/  &
    -1.376381921, -0.3249196962329, 0.3249196962329, 1.376381921, 0.0 /)
REAL, PARAMETER     :: a(2,5) = RESHAPE( (/ -0.726542528, -1.0,  &
    1.175570504, 0.381966011, 0.726542528, -1.0, 0.0, 1.236067978,  &
    1.175570504, 0.381966011 /), (/ 2, 5 /) )

x(1) = ABS(v(1))
x(2) = v(2)
INDEX = p
IF (x(1) == 0.0) THEN
  INDEX = 3
ELSE
  tangnt = x(2)/x(1)
  DO  i=1,p
    IF (tangnt <= tang(i)) THEN
      INDEX = i
      EXIT
    END IF
  END DO
END IF

fn_val = -1.0
DO  j=1,2
  fn_val = fn_val + a(j,INDEX)*x(j)
END DO

RETURN
END FUNCTION f6
