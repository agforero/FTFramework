MODULE constants_NSWC
! Contains the NSWC functions SPMPAR, DPMPAR, EPSLN, DEPSLN, EXPARG & DXPARG
!-----------------------------------------------------------------------
!     WRITTEN using F90 intrinsics by
!        Alan Miller
!        CSIRO Mathematical & Information Sciences
!        CLAYTON, VICTORIA, AUSTRALIA 3169
!     Latest revision - 1 February 1997
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)

CONTAINS

FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION spmpar



FUNCTION dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg

END MODULE constants_NSWC



SUBROUTINE smplx (a, b0, c, ka, m, n0, ind, ibasis, x, z, iter, mxiter,   &
                  numle, numge, bi, rerr)
!-----------------------------------------------------------------------
!     SIMPLEX PROCEDURE FOR SOLVING LINEAR PROGRAMMING PROBLEMS
!-----------------------------------------------------------------------
! Finds non-negative x's to maximize:
!     c1.x1 + ... + cn.xn
!
! Subject to the constraints:
!     a11.x1 + ... + a1n.xn {<=,=,>=} b1
!              ...
!     am1.x1 + ... + amn.xn {<=,=,>=} bm
!
! The constraints must be ordered so that the <= constraints, if any, come
! first, then any >= constraints, with the equality constraints last.
!
! Arguments:
! a(ka,n0)   INPUT   coefficients in the constraints.
! b0(m)      INPUT   right-hand side of the constraints.
! c(n0)      INPUT   vector of `costs' in the objective function.
! ka         INPUT   first dimension of array a.
! m          INPUT   dimension of array b0.
! n0         INPUT   2nd dimension of array a.
! ind        IN/OUT  If IND = 0 on input, the routine selects its own initial
!                    basis, else if IND = 1 then the indices of the initial
!                    basis should be in IBASIS.
!                    On output:
!                    IND = 0 the problem was solved
!                    IND = 1 the problem has no solution
!                    IND = 2 MXITER iterations were performed; more needed
!                    IND = 3 sufficient accuracy could not be maintained to
!                            solve the problem
!                    IND = 4 the problem has an unbounded solution
!                    IND = 5 input error detected
!                    IND = 6 the solution may have been obtained
! ibasis(m)  IN/OUT  indices of the variables in the basis.
! x()        OUTPUT  Dimension must be >= n + numle + numge.
!                    If IND = 0 or 6, it contains the values of the original,
!                    slack and surplus variables.
! z          OUTPUT  If IND = 0 or 6, contains the value of the objective.
! iter       OUTPUT  number of iterations used.
! mxiter     INPUT   maximum number of iterations.
! numle      INPUT   number of <= constraints.
! numge      INPUT   number of >= constraints.
! bi(m,m)    OUTPUT  the inverse of the basis matrix.
! rerr       OUTPUT  the estimated relative error achieved.
! N.B. The last two arguments in the NSWC routine have been eliminated.
!      These were work spaces.   They are now declared as automatic arrays.
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS JR.
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!------------------------
!     INITIAL VERSION  DEC  1977
!     LAST UPDATE      OCT  1990
!-----------------------------------------------------------------------
!     Converted using F90 intrinsics by
!        Alan Miller
!        CSIRO Mathematical & Information Sciences
!        CLAYTON, VICTORIA, AUSTRALIA 3169
!     Latest revision - 5 February 1997
!-----------------------------------------------------------------------

USE constants_NSWC
IMPLICIT NONE

INTEGER, INTENT(IN)                    :: ka, m, n0, mxiter, numle, numge
INTEGER, INTENT(IN OUT)                :: ind
INTEGER, DIMENSION(:), INTENT(IN OUT)  :: ibasis
INTEGER, INTENT(OUT)                   :: iter
REAL (dp), DIMENSION(:,:), INTENT(IN)  :: a
REAL (dp), DIMENSION(:), INTENT(IN)    :: b0, c
REAL (dp), INTENT(OUT)                 :: z, rerr
REAL (dp), DIMENSION(:), INTENT(OUT)   :: x
REAL (dp), DIMENSION(:,:), INTENT(OUT) :: bi

INTERFACE
  SUBROUTINE smplx1 (a, b0, c, ka, m, n0, ind, ibasis, r, z, iter, mxiter,  &
                     eps0, rerrmn, rerrmx, rerr, numle, numge, bi)
    USE constants_NSWC
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: ka, m, n0, mxiter, numle, numge
    INTEGER, INTENT(OUT)                   :: ind, iter
    REAL (dp), DIMENSION(:,:), INTENT(IN)  :: a
    REAL (dp), DIMENSION(:), INTENT(IN)    :: b0, c
    REAL (dp), DIMENSION(:,:), INTENT(OUT) :: bi
    REAL (dp), DIMENSION(:), INTENT(OUT)   :: r
    REAL (dp), INTENT(OUT)                 :: z, rerr
    REAL (dp), INTENT(IN)                  :: eps0, rerrmn, rerrmx
    INTEGER, DIMENSION(:), INTENT(IN OUT)  :: ibasis
  END SUBROUTINE smplx1
END INTERFACE

!------------------------
!     DIMENSION X(N0+NUMLE+NUMGE)
!------------------------

!     ********** EPS0 IS A MACHINE DEPENDENT PARAMETER. ASSIGN EPS0
!                THE VALUE U WHERE U IS THE SMALLEST POSITIVE FLOATING
!                POINT NUMBER SUCH THAT 1.0 + U .GT. 1.0.

!     Local variables
REAL (dp) :: eps0, rerrmn, rerrmx

eps0 = dpmpar(1)

!------------------------
rerrmn = 10.0_dp*eps0
rerrmx = 1.d-4
IF (eps0 < 1.d-13) rerrmx = 1.d-5

CALL smplx1(a, b0, c, ka, m, n0, ind, ibasis, x, z, iter, mxiter, eps0,   &
        rerrmn, rerrmx, rerr, numle, numge, bi)
RETURN
END SUBROUTINE smplx



SUBROUTINE smplx1 (a, b0, c, ka, m, n0, ind, ibasis, r, z, iter, mxiter,  &
                   eps0, rerrmn, rerrmx, rerr, numle, numge, bi)
!----------------------
!     NSTEP = 1   ELIMINATE THE NEGATIVE VARIABLES
!     NSTEP = 2   PHASE 1 OF THE SIMPLEX ALGORITHM
!     NSTEP = 3   PHASE 2 OF THE SIMPLEX ALGORITHM
!----------------------
!     MXITER = THE MAXIMUM NUMBER OF ITERATIONS PERMITTED
!     ITER = THE NUMBER OF THE CURRENT ITERATION
!     ICOUNT = THE NUMBER OF ITERATIONS SINCE THE LAST INVERSION
!----------------------
!     NUMLE = THE NUMBER OF .LE. CONSTRAINTS
!     NUMGE = THE NUMBER OF .GE. CONSTRAINTS
!----------------------
!     THE ROUTINE ASSUMES THAT THE .LE. CONSTRAINTS PRECEDE THE .GE.
!     CONSTRAINTS AND THAT THE .EQ. CONSTRAINTS COME LAST. THERE ARE
!     M CONSTRAINTS. X(N0+I) IS THE SLACK, SURPLUS, OR ARTIFICIAL
!     VARIABLE FOR THE I-TH CONSTRAINT (I=1, ..., M).
!----------------------
!     N0 = THE NUMBER OF ORGINAL VARIABLES
!     NS = THE NUMBER OF ORGINAL AND SLACK VARIABLES
!     N  = THE NUMBER OF ORGINAL, SLACK, AND SURPLUS VARIABLES
!     NUM = THE TOTAL NUMBER OF VARIABLES
!----------------------
!     RERRMN = THE SMALLEST RELATIVE ERROR TOLERANCE USED
!     RERRMX = THE LARGEST RELATIVE ERROR TOLERACE USED
!     RERR   = THE ESTIMATED CURRENT RELATIVE ERROR
!----------------------
!     ASSUME THAT
!         B0 = (B0(1), ..., B0(M))
!         C  = (C(1), ..., C(N0))
!         Z  = C(1)*X(1)+...+C(N0)*X(N0)
!     THE PROBLEM IS TO MAXIMIZE Z SUBJECT TO
!         AX(LE,EQ,GE)B0
!         X.GE.0
!----------------------
!     ON INPUT IND CAN HAVE THE VALUES
!         IND = 0   NO BEGINNING BASIS IS PROVIDED BY THE USER
!         IND = 1   THE ARRAY IBASIS HAS BEEN SET BY THE USER
!     ON OUTPUT IND IS ASSIGNED ONE OF THE VALUES
!         IND = 0   Z WAS SUCCESSFULLY MAXIMIZED
!         IND = 1   THE PROBLEM HAS NO FEASIBLE SOLUTION
!         IND = 2   MXITER ITERATIONS WERE PERFORMED
!         IND = 3   SUFFICIENT ACCURACY CANNOT BE MAINTAINED
!         IND = 4   THE PROBLEM HAS AN UNBOUNDED SOLUTION
!         IND = 5   THERE IS AN INPUT ERROR
!         IND = 6   Z WAS POSSIBLY MAXIMIZED
!----------------------
!     BASIS IS AN INTEGER ARRAY OF DIMENSION N0+M. FOR J.LE.N
!         BASIS(J) = 1  IF X(J) IS A BASIC VARIABLE
!         BASIS(J) = 0  IF X(J) IS NOT A BASIC VARIABLE
!     IF THE BASIC VARIABLES ARE X(I1), ..., X(IM) THEN
!         IBASIS = (I1, ..., IM)
!     ALSO XB(1), ..., XB(M) ARE THE CORRESPONDING VALUES OF THE
!     BASIC VARIABLES.
!----------------------
!     BI IS AN MXM ARRAY CONTAINING THE INVERSE OF THE BASIS MATRIX.
!----------------------
!     R IS AN ARRAY OF DIMENSION N. ON OUTPUT R CONTAINS THE CURRENT
!     VALUE OF X. DURING COMPUTATION R NORMALLY CONTAINS THE REDUCED
!     COSTS USED FOR THE SELECTION OF THE VARIABLE TO BE MADE BASIC.
!----------------------
USE constants_NSWC
IMPLICIT NONE
INTEGER, INTENT(IN)                    :: ka, m, n0, mxiter, numle, numge
INTEGER, INTENT(OUT)                   :: ind, iter
REAL (dp), DIMENSION(:,:), INTENT(IN)  :: a
REAL (dp), DIMENSION(:), INTENT(IN)    :: b0, c
REAL (dp), DIMENSION(:,:), INTENT(OUT) :: bi
REAL (dp), DIMENSION(:), INTENT(OUT)   :: r
REAL (dp), INTENT(OUT)                 :: z, rerr
REAL (dp), INTENT(IN)                  :: eps0, rerrmn, rerrmx
INTEGER, DIMENSION(:), INTENT(IN OUT)  :: ibasis
!----------------------

INTERFACE
  SUBROUTINE crout1(a, ka, n, iend, indx, temp, ierr)
    USE constants_NSWC
    IMPLICIT NONE
    INTEGER, INTENT(IN)                     :: ka, n, iend
    INTEGER, INTENT(OUT)                    :: ierr
    REAL (dp), DIMENSION(:), INTENT(IN OUT) :: a, temp
    INTEGER, DIMENSION(:), INTENT(IN OUT)   :: indx
  END SUBROUTINE crout1
END INTERFACE

!     Local variables
INTEGER   :: i, ibeg, icount, iend, ierr, ii, il, imin, iout, ip, j, jj,  &
             jmin, jp, k, ki, kj, l, ll, lrow, m0, mcheck, ms, n, npos,   &
             nrow, ns, nstep, num, bflag
REAL (dp) :: xb(m), y(m)
INTEGER   :: basis(m+n0), indx(m)
REAL (dp) :: amax, binorm, bmax, bmin, bnorm, cmin, const, eps, epsi, ratio, &
             rerr1, rmin, rtol, s, sgn, t, tol, total, w, xmax, zero = 0._dp, &
             dsum, dsump, dsumn, dt

!     ****** XMAX IS A MACHINE DEPENDENT CONSTANT. XMAX IS THE
!            LARGEST POSITIVE FLOATING POINT NUMBER.

xmax = dpmpar(3)

!----------------------
iter = 0
icount = 0
mcheck = MIN(5, 1 + m/15)
z = zero

!                CHECK FOR INPUT ERRORS

ms = numle + numge
IF (m < 2 .OR. n0 < 2 .OR. ms > m .OR. ka < m)  &
GO TO 12
DO i = 1, m
  IF (b0(i) < zero) GO TO 12
  xb(i) = zero
END DO
rtol = xmax
DO i = 1, n0
  IF (c(i) /= zero) rtol = MIN(ABS(c(i)), rtol)
END DO
rtol = rerrmx*rtol
GO TO 20

12 ind = 5
RETURN

!     FORMATION OF THE IBASIS AND BASIS ARRAYS.
!     (IF IND = 1  THEN THE IBASIS ARRAY IS DEFINED BY THE USER.)

20 ns = n0 + numle
n = ns + numge
IF (ind == 0) GO TO 30
num = n
DO i = 1, m
  IF (ibasis(i) > n) num = num + 1
END DO
GO TO 32
22 IF (ind == 0) GO TO 590
ind = 0

30 num = n0 + m
DO i = 1, m
  ibasis(i) = n0 + i
END DO
32 bflag = 0
basis(1:n) = 0
DO i = 1, m
  ki = ibasis(i)
  basis(ki) = 1
END DO
IF (ind == 1) GO TO 100

!          CALCULATION OF XB AND BI WHEN IND = 0

rerr = rerrmn
DO j = 1, m
  xb(j) = b0(j)
  bi(1:m,j) = zero
  bi(j,j) = 1.0_dp
END DO
IF (numge == 0) GO TO 630
jmin = numle + 1
DO j = jmin, ms
  xb(j) = -xb(j)
  bi(j,j) = -1.0_dp
END DO
GO TO 601

!                  REORDER THE BASIS

100 ibeg = 1
iend = m
DO i = 1, m
  IF (ibasis(i) > n0) THEN
    indx(ibeg) = ibasis(i)
    ibeg = ibeg + 1
  ELSE
    indx(iend) = ibasis(i)
    iend = iend - 1
  END IF
END DO
IF (iend == m) GO TO 22
ibasis(1:m) = indx(1:m)

!            REINVERSION OF THE BASIS MATRIX

DO j = 1, m
  kj = ibasis(j)
  IF (kj <= n0) GO TO 110
  IF (kj <= ns) GO TO 120
  IF (kj <= n) GO TO 130
  GO TO 120
  
  110 bi(1:m,j) = a(1:m,kj)
  CYCLE
  
  120 l = kj - n0
  bi(1:m,j) = zero
  bi(l,j) = 1.0_dp
  CYCLE
  
  130 l = kj - n0
  bi(1:m,j) = zero
  bi(l,j) = -1.0_dp
END DO

icount = 0
CALL crout1 (bi(1:,1), m, m, iend, indx, y, ierr)
IF (ierr /= 0) GO TO 580

!         CHECK THE ACCURACY OF BI AND RESET RERR

bnorm = zero
DO j = 1, m
  kj = ibasis(j)
  IF (kj <= n0) GO TO 140
  total = 1.0_dp
  GO TO 142
  140 total = SUM( ABS(a(1:m,kj)) )
  142 bnorm = MAX(bnorm, total)
END DO

binorm = zero
DO j = 1, m
  total = SUM( ABS(bi(1:m,j)) )
  binorm = MAX(binorm, total)
END DO
rerr = MAX(rerrmn, eps0*bnorm*binorm)
IF (rerr > 1.d-2) GO TO 580
bflag = 0

!                 RECALCULATION OF XB

DO i = 1, m
  dsump = zero
  dsumn = zero
  DO l = 1, m
    dt = bi(i,l)*b0(l)
    IF (dt > zero) THEN
      dsump = dsump + dt
    ELSE
      dsumn = dsumn + dt
    END IF
  END DO
  xb(i) = dsump + dsumn
  s = dsump
  t = dsumn
  tol = rerrmx*MAX(s, -t)
  IF (ABS(xb(i)) < tol) xb(i) = zero
END DO
GO TO 601

!     FIND THE NEXT VECTOR A(--, JP) TO BE INSERTED INTO THE BASIS

200 jp = 0
rmin = zero
IF (nstep == 3) rmin = -rtol
DO j = 1, n0
  IF (basis(j) /= 0) CYCLE
  IF (r(j) >= rmin) CYCLE
  jp = j
  rmin = r(j)
END DO
IF (n0 == n) GO TO 203
jmin = n0 + 1
rmin = rmin*1.1_dp
DO j = jmin, n
  IF (basis(j) /= 0) CYCLE
  IF (r(j) >= rmin) CYCLE
  jp = j
  rmin = r(j)
END DO
203 IF (jp /= 0) GO TO 300
IF (nstep < 2) THEN
   GO TO 800
ELSE IF (nstep == 2) THEN
  GO TO 230
ELSE
  GO TO 250
END IF

!     INSERT THE VALUES OF THE ORGINAL, SLACK, AND SURPLUS
!             VARIABLES INTO R, THEN TERMINATE.

220 r(1:n) = zero
DO i = 1, m
  ki = ibasis(i)
  IF (ki <= n) r(ki) = xb(i)
END DO
RETURN

!             COMPLETION OF THE NSTEP = 2 CASE

230 DO i = 1, m
  IF (ibasis(i) <= n) CYCLE
  IF (xb(i) > zero) GO TO 800
END DO
GO TO 680

240 IF (icount >= 5) GO TO 100
ind = 1
GO TO 220

!             COMPLETION OF THE NSTEP = 3 CASE

250 IF (rerr > 1.d-2) GO TO 251
ind = 0
GO TO 800
251 IF (icount >= 5) GO TO 100
ind = 6
GO TO 800

!     IF MXITER ITERATIONS HAVE NOT BEEN PERFORMED THEN BEGIN THE
!     NEXT ITERATION. COMPUTE THE JP-TH COLUMN OF BI*A AND STORE IT IN Y.

300 IF (iter < mxiter) GO TO 301
ind = 2
GO TO 220
301 iter = iter + 1
icount = icount + 1
IF (jp > ns) GO TO 330
IF (jp > n0) GO TO 320

nrow = 0
amax = zero
DO i = 1, m
  IF (a(i,jp) == zero) CYCLE
  nrow = nrow + 1
  indx(nrow) = i
  amax = MAX(ABS(a(i,jp)), amax)
END DO
IF (nrow /= 0) GO TO 310
ind = 4
GO TO 220

310 rerr1 = rerrmx*amax
DO i = 1, m
  dsum = zero
  DO ll = 1, nrow
    l = indx(ll)
    dsum = dsum + bi(i,l)*a(l,jp)
  END DO
  y(i) = dsum
  IF (ABS(y(i)) >= 5.d-3) CYCLE
  bmax = zero
  DO l = 1, m
    bmax = MAX(ABS(bi(i,l)), bmax)
  END DO
  tol = rerr1*bmax
  IF (ABS(y(i)) < tol) y(i) = zero
END DO
GO TO 350

320 l = jp - n0
y(1:m) = bi(1:m,l)
GO TO 350

330 l = jp - n0
y(1:m) = -bi(1:m,l)

350 DO i = 1, m
  IF (y(i) /= zero) GO TO 360
END DO
r(jp) = zero
iter = iter - 1
icount = icount - 1
GO TO 200

360 IF (nstep == 2) GO TO 430
IF (nstep > 2) GO TO 440

!     FINDING THE VARIABLE XB(IP) TO BE MADE NONBASIC FOR THE NSTEP = 1 CASE

npos = 0
ip = 0
eps = zero
epsi = xmax
DO i = 1, m
  IF (xb(i) < zero .OR. y(i) <= zero) CYCLE
  ratio = xb(i)/y(i)
  IF (ratio < epsi) THEN
    epsi = ratio
    npos = 1
    indx(1) = i
    CYCLE
  ELSE IF (ratio > epsi) THEN
    CYCLE
  END IF
  npos = npos + 1
  indx(npos) = i
END DO
IF (npos == 0) GO TO 420
IF (epsi == zero) GO TO 460

DO i = 1, m
  IF (xb(i) >= zero .OR. y(i) >= zero) CYCLE
  ratio = xb(i)/y(i)
  IF (ratio > epsi) CYCLE
  IF (ratio < eps) CYCLE
  eps = ratio
  ip = i
END DO
IF (ip /= 0) GO TO 500
GO TO 460

420 DO i = 1, m
  IF (xb(i) >= zero .OR. y(i) >= zero) CYCLE
  ratio = xb(i)/y(i)
  IF (ratio < eps) CYCLE
  eps = ratio
  ip = i
END DO
GO TO 500

!     FINDING THE VARIABLE XB(IP) TO BE MADE NONBASIC FOR THE NSTEP = 2 CASE

430 npos = 0
epsi = xmax
DO i = 1, m
  IF (y(i) <= zero) CYCLE
  ratio = xb(i)/y(i)
  IF (ratio < epsi) THEN
    epsi = ratio
    npos = 1
    indx(1) = i
  ELSE IF (ratio > epsi) THEN
    CYCLE
  END IF
  npos = npos + 1
  indx(npos) = i
END DO
GO TO 450

!     FINDING THE VARIABLE XB(IP) TO BE MADE NONBASIC FOR THE NSTEP = 3 CASE

440 npos = 0
epsi = xmax
DO i = 1, m
  IF (y(i) < zero) THEN
    IF (ibasis(i) <= n) CYCLE
    ip = i
    GO TO 500
  ELSE IF (y(i) > zero) THEN
    ratio = xb(i)/y(i)
    IF (ratio < epsi) THEN
      epsi = ratio
      npos = 1
      indx(1) = i
    ELSE IF (ratio > epsi) THEN
      CYCLE
    END IF
    npos = npos + 1
    indx(npos) = i
  END IF
END DO

450 IF (npos /= 0) GO TO 460
IF (icount >= 5) GO TO 100
ind = 4
GO TO 220

!              TIE BREAKING PROCEDURE

460 ip = indx(1)
IF (npos == 1) GO TO 500
ip = 0
bmin = xmax
cmin = xmax
DO ii = 1, npos
  i = indx(ii)
  l = ibasis(i)
  IF (l > n0) GO TO 461
  IF (c(l) <= zero) cmin = MIN(zero, cmin)
  IF (c(l) > cmin) CYCLE
  imin = i
  cmin = c(l)
  CYCLE
  461 IF (l <= n) GO TO 462
  ip = i
  GO TO 500
  462 lrow = l - n0
  s = b0(lrow)
  IF (lrow <= numle) THEN
    IF (s > bmin) CYCLE
    ip = i
    bmin = s
  ELSE
    s = -s
    bmin = MIN(zero, bmin)
    IF (s > bmin) CYCLE
    ip = i
    bmin = s
  END IF
END DO
IF (cmin <= zero .OR. ip == 0) ip = imin

!               TRANSFORMATION OF XB

500 IF (xb(ip) == zero) GO TO 510
const = xb(ip)/y(ip)
DO i = 1, m
  s = xb(i)
  xb(i) = xb(i) - const*y(i)
  IF (xb(i) >= zero) CYCLE
  IF (s >= zero .OR. xb(i) >= rerrmx*s) xb(i) = zero
END DO
xb(ip) = const

!               TRANSFORMATION OF BI

510 DO j = 1, m
  IF (bi(ip,j) == zero) CYCLE
  const = bi(ip,j)/y(ip)
  bi(1:m,j) = bi(1:m,j) - const*y(1:m)
  bi(ip,j) = const
END DO

!             UPDATING IBASIS AND BASIS

iout = ibasis(ip)
ibasis(ip) = jp
basis(iout) = 0
basis(jp) = 1
IF (iout > n) num = num - 1

!        CHECK THE ACCURACY OF BI AND RESET RERR

IF (rerr > 1.d-2) GO TO 530
k = 0
DO j = 1, m
  kj = ibasis(j)
  IF (kj > n0) CYCLE
  total = DOT_PRODUCT( bi(j,1:m), a(1:m, kj) )
  rerr = MAX(rerr, ABS(1.0_dp - total))
  k = k + 1
  IF (k >= mcheck) GO TO 522
END DO
522 IF (rerr <= 1.d-2) GO TO 600

!        THE ACCURACY CRITERIA ARE NOT SATISFIED

530 IF (icount < 5) GO TO 600
bflag = 1
GO TO 100

580 IF (iter == 0) GO TO 12
IF (bflag == 0) GO TO 590
bflag = 0
DO ip = 1, m
  IF (jp == ibasis(ip)) GO TO 582
END DO
582 ibasis(ip) = iout
basis(jp) = 0
basis(iout) = 1
IF (iout > n) num = num + 1
GO TO 100

590 ind = 3
GO TO 220

!       SET UP THE R ARRAY FOR THE NSTEP = 1 CASE

600 IF (nstep == 2) GO TO 630
IF (nstep > 2) GO TO 700
601 DO j = 1, m
  IF (xb(j) < zero) GO TO 610
END DO
GO TO 630

610 nstep = 1
m0 = 0
DO l = 1, m
  IF (xb(l) >= zero) CYCLE
  m0 = m0 + 1
  indx(m0) = l
END DO

DO j = 1, m
  dsump = zero
  dsumn = zero
  DO ll = 1, m0
    l = indx(ll)
    IF (bi(l,j) < zero) THEN
      dsumn = dsumn + bi(l,j)
    ELSE IF (bi(l,j) > zero) THEN
      dsump = dsump + bi(l,j)
    END IF
  END DO
  y(j) = dsump + dsumn
  s = dsump
  t = dsumn
  tol = rerrmx*MAX(s, -t)
  IF (ABS(y(j)) < tol) y(j) = zero
END DO
GO TO 650

!       SET UP THE R ARRAY FOR THE NSTEP = 2 CASE

630 IF (n == num) GO TO 680
nstep = 2
m0 = 0
DO l = 1, m
  IF (ibasis(l) <= n) CYCLE
  m0 = m0 + 1
  indx(m0) = l
END DO

DO j = 1, m
  dsump = zero
  dsumn = zero
  DO ll = 1, m0
    l = indx(ll)
    IF (bi(l,j) < zero) THEN
      dsumn = dsumn + bi(l,j)
    ELSE IF (bi(l,j) > zero) THEN
      dsump = dsump + bi(l,j)
    END IF
  END DO
  y(j) = -(dsump + dsumn)
  s = dsump
  t = dsumn
  tol = rerrmx*MAX(s, -t)
  IF (ABS(y(j)) < tol) y(j) = zero
END DO

650 DO j = 1, n0
  IF (basis(j) == 0) THEN
    r(j) = DOT_PRODUCT( y(1:m), a(1:m,j) )
  ELSE
    r(j) = zero
  END IF
END DO

660 IF (n0 == ns) GO TO 670
jmin = n0 + 1
DO j = jmin, ns
  r(j) = zero
  IF (basis(j) /= 0) CYCLE
  jj = j - n0
  r(j) = y(jj)
END DO

670 IF (ns == n) GO TO 200
jmin = ns + 1
DO j = jmin, n
  r(j) = zero
  IF (basis(j) /= 0) CYCLE
  jj = j - n0
  r(j) = -y(jj)
END DO
GO TO 200

!      SET UP A NEW R ARRAY FOR THE NSTEP = 3 CASE

680 nstep = 3
DO j = 1, m
  dsum = zero
  DO l = 1, m
    il = ibasis(l)
    IF (il <= n0) dsum = dsum + c(il)*bi(l,j)
  END DO
  y(j) = dsum
END DO

DO j = 1, n0
  r(j) = zero
  IF (basis(j) /= 0) CYCLE
  dsum = -c(j) + DOT_PRODUCT( y(1:m), a(1:m,j) )
  r(j) = dsum
  IF (r(j) >= zero) CYCLE
  tol = rerrmx*ABS(c(j))
  IF (ABS(r(j)) < tol) r(j) = zero
END DO
GO TO 660

!       UPDATE THE R ARRAY FOR THE NSTEP = 3 CASE

700 const = r(jp)
DO j = 1, n0
  IF (basis(j) /= 0) THEN
    r(j) = zero
  ELSE
    total = DOT_PRODUCT( bi(ip,1:m), a(1:m,j) )
    r(j) = r(j) - const*total
    IF (r(j) >= zero) CYCLE
    tol = rerrmx*ABS(c(j))
    IF (ABS(r(j)) < tol) r(j) = zero
  END IF
END DO

IF (n0 == ns) GO TO 720
jmin = n0 + 1
DO j = jmin, ns
  IF (basis(j) /= 0) THEN
    r(j) = zero
  ELSE
    jj = j - n0
    r(j) = r(j) - const*bi(ip,jj)
  END IF
END DO

720 IF (ns == n) GO TO 200
jmin = ns + 1
DO j = jmin, n
  IF (basis(j) /= 0) THEN
    r(j) = zero
  ELSE
    jj = j - n0
    r(j) = r(j) + const*bi(ip,jj)
  END IF
END DO
GO TO 200
!-----------------------------------------------------------------------
!               REFINE XB AND STORE THE RESULT IN Y
!-----------------------------------------------------------------------
800 y(1:m) = zero

m0 = 0
DO j = 1, m
  kj = ibasis(j)
  IF (kj <= n0) GO TO 810
  IF (kj <= ns) GO TO 820
  IF (kj <= n) GO TO 830
  GO TO 820
  
  810 m0 = m0 + 1
  indx(m0) = j
  CYCLE
  
  820 l = kj - n0
  y(l) = xb(j)
  CYCLE
  
  830 l = kj - n0
  y(l) = -xb(j)
END DO

IF (m0 == 0) THEN
  r(1:m) = b0(1:m) - y(1:m)
ELSE
  DO i = 1, m
    dsum = y(i)
    DO jj = 1, m0
      j = indx(jj)
      kj = ibasis(j)
      dsum = dsum + a(i,kj)*xb(j)
    END DO
    r(i) = b0(i) - dsum
  END DO
END IF

rerr1 = MIN(rerrmx, rerr)
DO i = 1, m
  y(i) = zero
  IF (xb(i) < zero) THEN
    sgn = -1.0_dp
    dsump = zero
    dsumn = xb(i)
  ELSE IF (xb(i) > zero) THEN
    sgn = 1.0_dp
    dsump = xb(i)
    dsumn = zero
  ELSE
    CYCLE
  END IF
  DO l = 1, m
    dt = bi(i,l)*r(l)
    IF (dt > zero) THEN
      dsump = dsump + dt
    ELSE
      dsumn = dsumn + dt
    END IF
  END DO
  w = dsump + dsumn
  IF (w == zero) CYCLE
  IF (sgn /= SIGN(1.0_dp, w)) CYCLE
  s = dsump
  t = dsumn
  tol = rerr1*MAX(s, -t)
  IF (ABS(w) > tol) y(i) = w
END DO
IF (nstep == 2) GO TO 870
IF (nstep > 2) GO TO 880

!         CHECK THE REFINEMENT (NSTEP = 1)

DO i = 1, m
  IF (y(i) >= zero) GO TO 861
  IF (y(i) < -rerrmx) GO TO 240
  y(i) = zero
  861 xb(i) = y(i)
END DO
GO TO 630

!         CHECK THE REFINEMENT (NSTEP = 2)

870 DO i = 1, m
  IF (ibasis(i) <= n) GO TO 871
  IF (y(i) > rerrmx) GO TO 240
  y(i) = zero
  871 xb(i) = y(i)
END DO
GO TO 680

!              COMPUTE Z  (NSTEP = 3)

880 dsum = zero
DO i = 1, m
  ki = ibasis(i)
  IF (ki > n0) GO TO 881
  dsum = dsum + c(ki)*y(i)
  881 xb(i) = y(i)
END DO
z = dsum
GO TO 220
END SUBROUTINE smplx1



SUBROUTINE crout1(a, ka, n, iend, indx, temp, ierr)
!     ******************************************************************
!     CROUT PROCEDURE FOR INVERTING MATRICES
!     ******************************************************************
!     A IS A MATRIX OF ORDER N WHERE N IS GREATER THAN OR EQUAL TO 1.
!     THE INVERSE OF A IS COMPUTED AND STORED IN A.

!     KA = LENGTH OF THE COLUMNS OF THE ARRAY A

!     IEND MAY BE 0, 1, ..., N-1.  IT IS ASSUMED THAT EACH OF THE FIRST
!     IEND COLUMNS OF THE MATRIX A CONTAINS ONLY ONE NONZERO ELEMENT
!     AND THAT THE NONZERO ELEMENT IS 1 OR -1.

!     indx IS AN ARRAY OF DIMENSION N-1 OR LARGER THAT IS USED BY THE
!     ROUTINE FOR KEEPING TRACK OF THE ROW INTERCHANGES THAT ARE MADE.

!     TEMP IS A TEMPORARY STORAGE ARRAY.

!     IERR REPORTS THE STATUS OF THE RESULTS. IF A IS NONSINGULAR THEN
!     THE INVERSE OF A IS COMPUTED AND IERR=0. OTHERWISE IF A IS FOUND
!     TO BE SINGULAR THEN IERR=1 AND THE ROUTINE TERMINATES.
!     --------------------
USE constants_NSWC
IMPLICIT NONE
INTEGER, INTENT(IN)                     :: ka, n, iend
INTEGER, INTENT(OUT)                    :: ierr
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: a, temp
INTEGER, DIMENSION(:), INTENT(IN OUT)   :: indx

!     Local variables
INTEGER   :: i, ibeg, ij, ik, il, j, jcol, jj, k, kcol, kcount, kj,  &
             kj0, kk, kl, km1, kp1, l, lj, lj0, lk, lmin, maxdim, mcol,  &
             ncol, nk, nm1, nmj, nmk, nn
REAL (dp) :: zero = 0._dp
REAL (dp) :: dsum, c, pmin, s

maxdim = ka*n
mcol = iend*ka
IF (iend == 0) GO TO 100

!           PROCESS THE FIRST IEND COLUMNS OF A

kcol = 0
DO k = 1, iend
  kk = kcol + k
  nk = kcol + n
  DO lk = kk, nk
    IF (a(lk) < zero) GO TO 20
    IF (a(lk) > zero) GO TO 30
  END DO
  GO TO 300
  
  20 l = lk - kcol
  lj0 = mcol + l
  DO lj = lj0, maxdim, ka
    a(lj) = -a(lj)
  END DO
  
  30 l = lk - kcol
  indx(k) = l
  IF (k == l) GO TO 32
  lj = lk
  DO kj = kk, maxdim, ka
    c = a(kj)
    a(kj) = a(lj)
    a(lj) = c
    lj = lj + ka
  END DO
  32 kcol = kcol + ka
END DO

!           PROCESS THE REMAINING COLUMNS OF A

100 nm1 = n - 1
ierr = 0
pmin = zero
ibeg = iend + 1
IF (ibeg == n) GO TO 190

k = ibeg
km1 = iend
kp1 = k + 1
kcol = mcol
kk = kcol + k
DO kcount = ibeg, nm1
  
!     SEARCH FOR THE K-TH PIVOT ELEMENT (K=IBEG, ..., N-1)
  
  l = k
  s = ABS(a(kk))
  DO i = kp1, n
    ik = kcol + i
    c = ABS(a(ik))
    IF (s >= c) CYCLE
    l = i
    s = c
  END DO
  
  IF (k > ibeg .AND. s >= pmin) GO TO 120
  pmin = s
  IF (s == zero) GO TO 300
  
!              INTERCHANGING ROWS K AND L
  
  120 indx(k) = l
  IF (k == l) GO TO 130
  kj0 = mcol + k
  lj  = mcol + l
  DO kj = kj0, maxdim, ka
    c = a(kj)
    a(kj) = a(lj)
    a(lj) = c
    lj = lj + ka
  END DO
  
!       COMPUTE THE K-TH ROW OF U (K=IBEG, ..., N-1)
  
  130 c = a(kk)
  IF (k > ibeg) GO TO 140
  kj0 = kk + ka
  DO kj = kj0, maxdim, ka
    a(kj) = a(kj)/c
  END DO
  GO TO 160
  
  140 kl = mcol + k
  DO l = ibeg, km1
    temp(l) = a(kl)
    kl = kl + ka
  END DO
  
  kj0 = kk + ka
  DO kj = kj0, maxdim, ka
    jcol = kj - k
    dsum = -a(kj)
    DO l = ibeg, km1
      lj = jcol + l
      dsum = dsum + temp(l)*a(lj)
    END DO
    a(kj) = -dsum/c
  END DO
  
!      COMPUTE THE K-TH COLUMN OF L (K=IBEG+1, ..., N)
  
  160 km1 = k
  k = kp1
  kp1 = k + 1
  kcol = kcol + ka
  kk = kcol + k
  DO l = ibeg, km1
    lk = kcol + l
    temp(l) = a(lk)
  END DO
  
  DO i = k, n
    il = mcol + i
    dsum = zero
    DO l = ibeg, km1
      dsum = dsum + a(il)*temp(l)
      il = il + ka
    END DO
    a(il) = a(il) - dsum
  END DO
END DO

!           CHECK THE N-TH PIVOT ELEMENT

190 ncol = maxdim - ka
nn = ncol + n
c = ABS(a(nn))
IF (c > pmin) GO TO 200
IF (c == zero) GO TO 300

!          REPLACE L WITH THE INVERSE OF L

200 IF (ibeg == n) GO TO 213
jj = mcol + ibeg
i = ka + 1
DO j = ibeg, nm1
  a(jj) = 1.0_dp/a(jj)
  temp(j) = a(jj)
  kj = jj
  DO km1 = j, nm1
    k = km1 + 1
    kj = kj + 1
    dsum = zero
    kl = kj
    DO l = j, km1
      dsum = dsum + a(kl)*temp(l)
      kl = kl + ka
    END DO
    a(kj) = -dsum/a(kl)
    temp(k) = a(kj)
  END DO
  jj = jj + i
END DO
213 a(nn) = 1.0_dp/a(nn)
IF (n == 1) RETURN

!       SOLVE UX = Y WHERE Y IS THE INVERSE OF L

DO nmk = 1, nm1
  k = n - nmk
  lmin = MAX(ibeg, k+1)
  kl = (lmin-1)*ka + k
  DO l = lmin, n
    temp(l) = a(kl)
    a(kl) = zero
    kl = kl + ka
  END DO
  
  kj0 = mcol + k
  DO kj = kj0, maxdim, ka
    dsum = -a(kj)
    lj = (kj - k) + lmin
    DO l = lmin, n
      dsum = dsum + temp(l)*a(lj)
      lj = lj + 1
    END DO
    a(kj) = -dsum
  END DO
END DO

!                 COLUMN INTERCHANGES

jcol = ncol - ka
DO nmj = 1, nm1
  j = n - nmj
  k = indx(j)
  IF (j == k) GO TO 251
  ij = jcol
  ik = (k-1)*ka
  DO i = 1, n
    ij = ij + 1
    ik = ik + 1
    c = a(ij)
    a(ij) = a(ik)
    a(ik) = c
  END DO
  251 jcol = jcol - ka
END DO
RETURN

!                    ERROR RETURN

300 ierr = 1
RETURN
END SUBROUTINE crout1
