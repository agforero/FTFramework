MODULE twod_quad
IMPLICIT NONE

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

CONTAINS


SUBROUTINE twodqd(f, n, x, y, tol, iclose, maxtri, mevals, result,  &
                  error, nu, nd, nevals, iflag)

!  Date written 840518 (YYMMDD)
!  Revision date 840518
!  ELF90 version by Alan Miller 18 August 1997
!  N.B. The last 2 calling arguments have been removed from this version.
!  Category no. D11
!  Keywords: QUADRATURE, TWO DIMENSIONAL, ADAPTIVE, CUBATURE
!  Author: Kahaner, D.K., N.B.S. and Rechard, O.W., Univ. of Denver

!  Purpose: To compute the two-dimensional integral of a function
!           F over a region consisting of N triangles.

!  Description:
!  A total error estimate is obtained and compared with a tolerance - TOL
!  - that is provided as input to the subroutine.   The error tolerance is
!  treated as either relative or absolute depending on the input value of
!  IFLAG.   A 'local quadrature module' is applied to each input triangle and
!  estimates of the total integral and the total error are computed.  The local
!  quadrature module is either subroutine LQM0 or LQM1; the choice between
!  them is determined by the value of the input variable ICLOSE.

!  If the total error estimate exceeds the tolerance, the triangle with the
!  largest absolute error is divided into two triangles by a median to its
!  longest side.   The local quadrature module is then applied to each of
!  the sub-triangles to obtain new estimates.   This process is repeated until
!  either 1) the error tolerance is satisfied, 2) the number of triangles
!  generated exceeds the input parameter MAXTRI, 3) the number of integrand
!  evaluations exceeds the input parameter MEVALS, or 4) the subroutine
!  senses that roundoff error is beginning to contaminate the result.

!  The user must specify MAXTRI.   The user must also specify MEVALS.   This
!  number will be effective in limiting the computation only if it is less than
!  94 * MAXTRI when LQM1 is specified, or 54 * MAXTRI when LQM0 is specified.

!  After the subroutine has returned to the calling program with output values.
!  It can be called again with a smaller value of TOL, and/or a different value
!  of MEVALS.   The tolerance can also be changed from relative to absolute or
!  vice-versa by changing IFLAG.   Unless the parameters NU and ND are reset
!  to zero, the subroutine will restart with the final set of triangles and
!  output values from the previous call.

!  Arguments:

!  F   Function subprogram defining the integrand F(U,V).   An INTERFACE
!      for F must be included in the calling program.

!  N   The number of input triangles.

!  X   A 3 by N array containing the abscissae of the vertices of
!      the N triangles.

!  Y   A 3 by N array containing the ordinates of the vertices of
!      the N triangles.

!  TOL The desired bound on the error.  If IFLAG = 0 on input, TOL
!      is interpreted as a bound on the relative error; if IFLAG = 1
!      the bound is on the absolute error.

!  ICLOSE If ICLOSE = 1 then LQM1 is used, otherwise LQM0 is used.
!      LQM0 uses function values only at interior points of the triangle.
!      LQM1 is usually more accurate but involves evaluating the integrand
!      at more points including some on the boundary of the triangle.
!      It will usually be better to use LQM1 unless the integrand has
!      singularities on the boundary of the triangle.

!  MAXTRI The maximum number of triangles allowed to be generated.

!  MEVALS Maximum number of function evaluations allowed.

!  RESULT Estimate of the integral on output.

!  ERROR Estimate of the absolute total error on output.

!  NU  An integer which must be set to zero for the first call of this routine.
!      Subsequent calls should use the value output from the previous call.

!  ND  An integer which must be set to zero for the first call of this routine.
!      Subsequent calls should use the value output from the previous call.

!  NEVALS The actual number of function evaluations.

!  IFLAG See TOL above for description of input value.
!      On output:
!      IFLAG = 0 means normal termination.
!            = 1 means lack of space to divide more triangles.
!            = 2 means termination because of roundoff noise.
!            = 3 means termination with relative error <= 5 times
!                machine epsilon.
!            = 4 means number of function evaluations > MEVALS.
!            = 9 means IFLAG was not set = 0 or 1 on input.

!  Routines called: HINITD, HINITU, HPACC, HPDEL, HPINS, LQM0, LQM1, TRIDIV

IMPLICIT NONE
INTEGER, INTENT(IN)       :: n, iclose, mevals, maxtri
INTEGER, INTENT(IN OUT)   :: nu, nd, iflag
INTEGER, INTENT(OUT)      :: nevals
REAL (dp), INTENT(IN)     :: x(:,:), y(:,:)
REAL (dp), INTENT(IN OUT) :: tol
REAL (dp), INTENT(OUT)    :: result, error
! In F77 version   REAL    x(3,n), y(3,n)

INTERFACE
  FUNCTION f(x, y) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: x, y
    REAL (dp)             :: fn_val
  END FUNCTION f

  FUNCTION hpfun(a, b) RESULT(yes)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: a, b
    LOGICAL               :: yes
  END FUNCTION hpfun
END INTERFACE

!     Local variables.

INTEGER              :: rndcnt, fadd, i, j, maxp1, iwork(2*maxtri)
LOGICAL              :: full
REAL (dp)            :: a, r, e, u(3), v(3), node(9), node1(9), node2(9),  &
                        epsabs, emach, newres, newerr, small, data(9*maxtri)
REAL (dp), SAVE      :: atot
REAL (dp), PARAMETER :: zero = 0.0_dp, five = 5.0_dp, half = 0.5_dp,  &
                        pt99 = 0.99_dp
!---------------------------------------------+
!                                             |
!     WARNING - Machine dependent constants   |
!                                             |
!---------------------------------------------+
emach = EPSILON(zero)
small = SQRT(emach)

!     If the heaps are empty, apply LQM to each input triangle and
!     place all of the data on the second heap.

maxp1 = maxtri + 1
IF (nu + nd == 0) THEN
  CALL hinitu(maxtri, 9, nu, iwork)
  CALL hinitd(maxtri, 9, nd, iwork(maxp1:))
  atot = zero
  result = zero
  error = zero
  rndcnt = 0
  nevals = 0
  DO i = 1, n
    DO j = 1, 3
      u(j) = x(j,i)
      v(j) = y(j,i)
    END DO
    a = half * ABS(u(1)*v(2) + u(2)*v(3) + u(3)*v(1) - u(1)*v(3)  &
        - u(2)*v(1) - u(3)*v(2))
    atot = atot + a
    IF (iclose == 1) THEN
      CALL lqm1(f, u, v, r, e)
      nevals = nevals + 47
    ELSE
      CALL lqm0(f, u, v, r, e)
      nevals = nevals + 28
    END IF
    result = result + r
    error = error + e
    node(1) = e
    node(2) = r
    node(3) = x(1,i)
    node(4) = y(1,i)
    node(5) = x(2,i)
    node(6) = y(2,i)
    node(7) = x(3,i)
    node(8) = y(3,i)
    node(9) = a
    CALL hpins(maxtri, 9, DATA, nd, iwork(maxp1:), node, greatr)
  END DO
END IF

!     Check that input tolerance is consistent with machine epsilon.

IF (iflag == 0) THEN
  IF (tol <= five * emach) THEN
    tol = five * emach
    fadd = 3
  ELSE
    fadd = 0
  END IF
  epsabs = tol * ABS(result)
ELSE IF (iflag == 1) THEN
  IF (tol <= five * emach * ABS(result)) THEN
    epsabs = five * emach * ABS(result)
  ELSE
    epsabs = tol
  END IF
ELSE
  iflag = 9
  RETURN
END IF

!     Adjust the second heap using the current EPSABS.

2 IF (nd == 0) GO TO 40
j = nd
3 IF (j == 0) GO TO 40
CALL hpacc(maxtri, 9, DATA, nd, iwork(maxp1:), node, j)
IF (node(1) > epsabs * node(9) / atot) THEN
  CALL hpins(maxtri, 9, DATA, nu, iwork, node, greatr)
  CALL hpdel(maxtri, DATA, nd, iwork(maxp1:), greatr, j)
  IF (j > nd) j = j - 1
ELSE
  j = j - 1
END IF
GO TO 3

!     Beginning of main loop, from here to end.

40 IF (nevals >= mevals) THEN
  iflag = 4
  RETURN
END IF
IF (error <= epsabs) THEN
  IF (iflag == 0) THEN
    IF (error <= ABS(result) * tol) THEN
      iflag = fadd
      RETURN
    ELSE
      epsabs = ABS(result) * tol
      GO TO 2
    END IF
  ELSE
    IF (error <= tol) THEN
      iflag = 0
      RETURN
    ELSE IF (error <= five * emach * ABS(result)) THEN
      iflag = 3
      RETURN
    ELSE
      epsabs = five * emach * ABS(result)
      GO TO 2
    END IF
  END IF
END IF

!     If there are too many triangles and the second heap is not empty,
!     remove the bottom triangle from the second heap.   If the second
!     heap is empty, return with IFLAG = 1 or 4.

IF (nu + nd >= maxtri) THEN
  full = .true.
  IF (nd > 0) THEN
    iwork(nu+1) = iwork(maxtri+nd)
    nd = nd - 1
  ELSE
    iflag = 1
    RETURN
  END IF
ELSE
  full = .false.
END IF

!     Find triangle with largest error, divide it into two,
!     and apply LQM to each half.

IF (nd == 0) THEN
  CALL hpacc(maxtri, 9, DATA, nu, iwork, node, 1)
  CALL hpdel(maxtri, DATA, nu, iwork, greatr, 1)
ELSE IF (nu == 0) THEN
  CALL hpacc(maxtri, 9, DATA, nd, iwork(maxp1:), node, 1)
  CALL hpdel(maxtri, DATA, nd, iwork(maxp1:), greatr, 1)
ELSE IF (DATA(iwork(1)) >= DATA(iwork(maxp1))) THEN
  IF (full) iwork(maxtri+nd+2) = iwork(nu)
  CALL hpacc(maxtri, 9, DATA, nu, iwork, node, 1)
  CALL hpdel(maxtri, DATA, nu, iwork, greatr, 1)
ELSE
  IF (full) iwork(nu+2) = iwork(maxtri+nd)
  CALL hpacc(maxtri, 9, DATA, nd, iwork(maxp1:), node, 1)
  CALL hpdel(maxtri, DATA, nd, iwork(maxp1:), greatr, 1)
END IF
CALL tridiv(node, node1, node2, half, 1)
DO j = 1, 3
  u(j) = node1(2*j+1)
  v(j) = node1(2*j+2)
END DO
IF (iclose == 1) THEN
  CALL lqm1(f, u, v, node1(2), node1(1))
  nevals = nevals + 47
ELSE
  CALL lqm0(f, u, v, node1(2), node1(1))
  nevals = nevals + 28
END IF
DO j = 1, 3
  u(j) = node2(2*j+1)
  v(j) = node2(2*j+2)
END DO
IF (iclose == 1) THEN
  CALL lqm1(f, u, v, node2(2), node2(1))
  nevals = nevals + 47
ELSE
  CALL lqm0(f, u, v, node2(2), node2(1))
  nevals = nevals + 28
END IF

newerr = node1(1) + node2(1)
newres = node1(2) + node2(2)
IF (newerr > pt99 * node(1)) THEN
  IF (ABS(node(2) - newres) <= small * ABS(newres)) rndcnt = rndcnt + 1
END IF
result = result - node(2) + newres
error = error - node(1) + newerr

IF (node1(1) > node1(9) * epsabs / atot) THEN
  CALL hpins(maxtri, 9, DATA, nu, iwork, node1, greatr)
ELSE
  CALL hpins(maxtri, 9, DATA, nd, iwork(maxp1:), node1, greatr)
END IF
IF (node2(1) > node2(9) * epsabs / atot) THEN
  CALL hpins(maxtri, 9, DATA, nu, iwork, node2, greatr)
ELSE
  CALL hpins(maxtri, 9, DATA, nd, iwork(maxp1:), node2, greatr)
END IF

IF (rndcnt >= 20) THEN
  iflag = 2
  RETURN
END IF

IF (iflag == 0) THEN
  IF (epsabs < half * tol * ABS(result)) THEN
    epsabs = tol * ABS(result)
    j = nu
    5     IF (j == 0) GO TO 40
    CALL hpacc(maxtri, 9, DATA, nu, iwork, node, j)
    IF (node(1) <= epsabs * node(9) / atot) THEN
      CALL hpins(maxtri, 9, DATA, nd, iwork(maxp1:), node, greatr)
      CALL hpdel(maxtri, DATA, nu, iwork, greatr, j)
      IF (j > nu) j = j - 1
    ELSE
      j = j - 1
    END IF
    GO TO 5
  END IF
END IF
GO TO 40
END SUBROUTINE twodqd



SUBROUTINE tridiv(node, node1, node2, coef, rank)
IMPLICIT NONE
REAL (dp), INTENT(IN)  :: node(:), coef
REAL (dp), INTENT(OUT) :: node1(:), node2(:)
INTEGER, INTENT(IN)    :: rank

!     Local variables.

REAL (dp) :: s(3), coef1, temp, one = 1.0_dp
INTEGER   :: t(3), i, j

coef1 = one - coef
s(1) = (node(3) - node(5))**2 + (node(4) - node(6))**2
s(2) = (node(5) - node(7))**2 + (node(6) - node(8))**2
s(3) = (node(3) - node(7))**2 + (node(4) - node(8))**2
t(1) = 1
t(2) = 2
t(3) = 3
DO i = 1, 2
  DO j = i+1, 3
    IF (s(i) < s(j)) THEN
      temp = t(i)
      t(i) = t(j)
      t(j) = temp
    END IF
  END DO
END DO

IF (t(rank) == 1) THEN
  node1(3) = coef * node(3) + coef1 * node(5)
  node1(4) = coef * node(4) + coef1 * node(6)
  node1(5) = node(5)
  node1(6) = node(6)
  node1(7) = node(7)
  node1(8) = node(8)
  node2(3) = node1(3)
  node2(4) = node1(4)
  node2(5) = node(7)
  node2(6) = node(8)
  node2(7) = node(3)
  node2(8) = node(4)
ELSE IF (t(rank) == 2) THEN
  node1(3) = coef * node(5) + coef1 * node(7)
  node1(4) = coef * node(6) + coef1 * node(8)
  node1(5) = node(7)
  node1(6) = node(8)
  node1(7) = node(3)
  node1(8) = node(4)
  node2(3) = node1(3)
  node2(4) = node1(4)
  node2(5) = node(3)
  node2(6) = node(4)
  node2(7) = node(5)
  node2(8) = node(6)
ELSE
  node1(3) = coef * node(3) + coef1 * node(7)
  node1(4) = coef * node(4) + coef1 * node(8)
  node1(5) = node(3)
  node1(6) = node(4)
  node1(7) = node(5)
  node1(8) = node(6)
  node2(3) = node1(3)
  node2(4) = node1(4)
  node2(5) = node(5)
  node2(6) = node(6)
  node2(7) = node(7)
  node2(8) = node(8)
END IF
node1(9) = coef * node(9)
node2(9) = coef1 * node(9)
RETURN
END SUBROUTINE tridiv



SUBROUTINE lqm1(f, u, v, res11, est)

!     Purpose:
!       Computes the integral of F over the triangle with vertices
!       (U(1),V(1)), (U(2),V(2)), (U(3),V(3)) and estimates the error.

!     Parameters:
!     F     Function to be integrated.   The actual name must be
!           declared EXTERNAL  in the calling program.
!     U(1)..U(3)  Abscissae of vertices
!     V(1)..V(3)  Ordinates of vertices
!     RES11 Approximation to the integral using the Lyness & Jespersen
!           rule of degree 11 and using 28 points
!     EST   Estimate of the absolute error

!     Remarks:
!       Date of last update: 18 JAN 1984 D. Kahaner, NBS

!     Function called: F (User supplied)

IMPLICIT NONE
REAL (dp), INTENT(IN)  :: u(:), v(:)
REAL (dp), INTENT(OUT) :: res11, est

INTERFACE
  FUNCTION f(x, y) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: x, y
    REAL (dp)             :: fn_val
  END FUNCTION f
END INTERFACE


!     Local variables

REAL (dp) :: area, df0, dresc, emach, fv(19), f0, res9, u1, u2, u3, uflow,  &
             v1, v2, v3, x(3), y(3), z1, z2, z3, resab9
INTEGER   :: j, kount, l
REAL (dp), PARAMETER :: half = 0.5_dp, three = 3.0_dp, zero = 0.0_dp,  &
                        one = 1.0_dp, twenty = 20.0_dp, onept5 = 1.5_dp

!     First homogeneous coordinates of degree-9 and degree-11 points
!     taken with multiplicity 3.

REAL (dp), PARAMETER :: zeta1(15) = (/ 0.2063496160252593D-01, 0.1258208170141290D+00, &
        0.6235929287619356D+00, 0.9105409732110941D+00, 0.3683841205473626D-01, &
        0.7411985987844980D+00, 0.9480217181434233D+00, 0.8114249947041546D+00, &
        0.1072644996557060D-01, 0.5853132347709715D+00, 0.1221843885990187D+00, &
        0.4484167758913055D-01, 0.6779376548825902D+00, 0.0D+00,  &
        0.8588702812826364D+00 /)

!     Second homogeneous coordinates.

REAL (dp), PARAMETER :: zeta2(15) = (/ 0.4896825191987370D+00, 0.4370895914929355D+00, &
        0.1882035356190322D+00, 0.4472951339445297D-01, 0.7411985987844980D+00, &
        0.3683841205473626D-01, 0.2598914092828833D-01, 0.9428750264792270D-01, &
        0.4946367750172147D+00, 0.2073433826145142D+00, 0.4389078057004907D+00, &
        0.6779376548825902D+00, 0.4484167758913055D-01, 0.8588702812826364D+00, &
        0.0D+00 /)

!     Weights of mid-point of triangle in 9-point & 11-point formulae.

REAL (dp), PARAMETER :: w90 = 0.9713579628279610D-01, &
                        w110 = 0.8797730116222190D-01

!     Weights in degree-9 & degree-11 rules.

REAL (dp), PARAMETER :: w(15) = (/ 0.3133470022713983D-01, 0.7782754100477543D-01,  &
    0.7964773892720910D-01, 0.2557767565869810D-01, 0.4328353937728940D-01,  &
    0.4328353937728940D-01, 0.8744311553736190D-02, 0.3808157199393533D-01,  &
    0.1885544805613125D-01, 0.7215969754474100D-01, 0.6932913870553720D-01,  &
    0.4105631542928860D-01, 0.4105631542928860D-01, 0.7362383783300573D-02,  &
    0.7362383783300573D-02 /)

!     LIST OF MAJOR VARIABLES
!     -----------------------
!     AREA  Area of the triangle
!     DRESC Approximation to the integral of |F - INTEGRAL/AREA|
!     RESAB9 Approximation to the integral of |F| over the triangle
!     X     Abscissae of integration points
!     Y     Ordinates of the integration points
!     FV    Function values

!-----------------------------------+
!                                   |
!     Machine dependent constants   |
!                                   |
!-----------------------------------+
emach = EPSILON(zero)
uflow = TINY(zero)

!     Estimate degree-9 and degree-11 results for INTEGRAL/AREA and
!     degree-9 result for |F|.

u1 = u(1)
u2 = u(2)
u3 = u(3)
v1 = v(1)
v2 = v(2)
v3 = v(3)
area = ABS(u1*v2 - u2*v1 - u1*v3 + v1*u3 + u2*v3 - v2*u3) * half
f0 = f((u1 + u2 + u3)/three, (v1 + v2 + v3)/three)
res9 = f0 * w90
resab9 = ABS(f0) * w90
fv(1) = f0
kount = 1
res11 = f0 * w110
DO j = 1, 15
  z1 = zeta1(j)
  z2 = zeta2(j)
  z3 = one - z1 - z2
  x(1) = z1*u1 + z2*u2 + z3*u3
  y(1) = z1*v1 + z2*v2 + z3*v3
  x(2) = z2*u1 + z3*u2 + z1*u3
  y(2) = z2*v1 + z3*v2 + z1*v3
  x(3) = z3*u1 + z1*u2 + z2*u3
  y(3) = z3*v1 + z1*v2 + z2*v3
  IF (j <= 6) THEN
    f0 = zero
    df0 = zero
    DO l = 1, 3
      kount = kount + 1
      fv(kount) = f(x(l), y(l))
      f0 = f0 + fv(kount)
      df0 = df0 + ABS(fv(kount))
    END DO
    res9 = res9 + f0*w(j)
    resab9 = resab9 + df0*w(j)
  ELSE
    f0 = f(x(1), y(1)) + f(x(2), y(2)) + f(x(3), y(3))
    res11 = res11 + f0*w(j)
  END IF
END DO

!     Compute degree-9 approximation to the integral of
!     |F - INTEGRAL/AREA|.

dresc = ABS(fv(1) - res9) * w90
kount = 2
DO j = 1, 6
  dresc = dresc + (ABS(fv(kount) - res9) + ABS(fv(kount+1) - res9)  &
  + ABS(fv(kount+2) - res9)) * w(j)
  kount = kount + 3
END DO

!     Compute degree-9 and degree-11 approximations to the integral and
!     its error.

res9 = res9 * area
res11 = res11 * area
resab9 = resab9 * area
dresc = dresc * area
est = ABS(res9 - res11)
IF (dresc /= zero) est = MAX(est, dresc * MIN(one, (twenty * est  &
                                  / dresc)**onept5))
IF (resab9 > uflow) est = MAX(emach * resab9, est)
RETURN
END SUBROUTINE lqm1



SUBROUTINE lqm0(f, u, v, res8, est)

!     Purpose:
!       Computes the integral of F over the triangle with vertices
!       (U(1),V(1)), (U(2),V(2)), (U(3),V(3)) and estimates the error.

!     Parameters:
!     F     Function to be integrated.   The actual name must be
!           declared EXTERNAL  in the calling program.
!     U(1)..U(3)  Abscissae of vertices
!     V(1)..V(3)  Ordinates of vertices
!     RES8 Approximation to the integral using the Lyness & Jespersen
!           rule of degree 8 and using 16 points
!     EST   Estimate of the absolute error

!     Remarks:
!       Date of last update: 10 APRIL 1984 O.W. Rechard, NBS

!     Function called: F (User supplied)

IMPLICIT NONE
REAL (dp), INTENT(IN)  :: u(:), v(:)
REAL (dp), INTENT(OUT) :: res8, est

INTERFACE
  FUNCTION f(x, y) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: x, y
    REAL (dp)             :: fn_val
  END FUNCTION f
END INTERFACE

!     Local variables

REAL (dp) :: area, df0, dresc, emach, fv(19), f0, res6, u1, u2, u3, uflow,  &
             v1, v2, v3, x(3), y(3), z1, z2, z3, resab6
INTEGER   :: j, kount, l
REAL (dp), PARAMETER :: half = 0.5_dp, three = 3.0_dp, zero = 0.0_dp,   &
                        one = 1.0_dp, twenty = 20.0_dp, onept5 = 1.5_dp

!     First homogeneous coordinates of degree-6 and degree-8 points
!     taken with multiplicity 3.

REAL (dp), PARAMETER :: zeta1(9) = (/ 0.5014265096581342D+00,  &
              0.8738219710169965D+00, 0.6365024991213939D+00,  &
              0.5314504984483216D+00, 0.8141482341455413D-01,  &
              0.8989055433659379D+00, 0.6588613844964797D+00,  &
              0.8394777409957211D-02, 0.7284923929554041D+00 /)

!     Second homogeneous coordinates.

REAL (dp), PARAMETER :: zeta2(9) = (/ 0.2492867451709329D+00,  &
              0.6308901449150177D-01, 0.5314504984483216D+00,  &
              0.6365024991213939D+00, 0.4592925882927229D+00,  &
              0.5054722831703103D-01, 0.1705693077517601D+00,  &
              0.7284923929554041D+00, 0.8394777409957211D-00 /)

!     Weights of mid-point of triangle in 6-point & 8-point formulae.

REAL (dp), PARAMETER :: w60 = 0.0_dp, w80 = 0.1443156076777862D+00

!     Weights in degree-6 & degree-8 rules.

REAL (dp), PARAMETER :: w(9) = (/ 0.1167862757263407D+00,  &
          0.5084490637020547D-01, 0.8285107561839291D-01,  &
          0.8285107561839291D-01, 0.9509163426728497D-01,  &
          0.3245849762319813D-01, 0.1032173705347184D+00,  &
          0.2723031417443487D-01, 0.2723031417443487D-01 /)

!     LIST OF MAJOR VARIABLES
!     -----------------------
!     AREA  Area of the triangle
!     DRESC Approximation to the integral of |F - INTEGRAL/AREA|
!     RESAB6 Approximation to the integral of |F| over the triangle
!     X     Abscissae of integration points
!     Y     Ordinates of the integration points
!     FV    Function values

!-----------------------------------+
!                                   |
!     Machine dependent constants   |
!                                   |
!-----------------------------------+
emach = EPSILON(zero)
uflow = TINY(zero)

!     Estimate degree-6 and degree-8 results for INTEGRAL/AREA and
!     degree-6 result for |F|.

u1 = u(1)
u2 = u(2)
u3 = u(3)
v1 = v(1)
v2 = v(2)
v3 = v(3)
area = ABS(u1*v2 - u2*v1 - u1*v3 + v1*u3 + u2*v3 - v2*u3) * half
f0 = f((u1 + u2 + u3)/three, (v1 + v2 + v3)/three)
res6 = f0 * w60
resab6 = ABS(f0) * w60
fv(1) = f0
kount = 1
res8 = f0 * w80
DO j = 1, 9
  z1 = zeta1(j)
  z2 = zeta2(j)
  z3 = one - z1 - z2
  x(1) = z1*u1 + z2*u2 + z3*u3
  y(1) = z1*v1 + z2*v2 + z3*v3
  x(2) = z2*u1 + z3*u2 + z1*u3
  y(2) = z2*v1 + z3*v2 + z1*v3
  x(3) = z3*u1 + z1*u2 + z2*u3
  y(3) = z3*v1 + z1*v2 + z2*v3
  IF (j <= 4) THEN
    f0 = zero
    df0 = zero
    DO l = 1, 3
      kount = kount + 1
      fv(kount) = f(x(l), y(l))
      f0 = f0 + fv(kount)
      df0 = df0 + ABS(fv(kount))
    END DO
    res6 = res6 + f0*w(j)
    resab6 = resab6 + df0*w(j)
  ELSE
    f0 = f(x(1), y(1)) + f(x(2), y(2)) + f(x(3), y(3))
    res8 = res8 + f0*w(j)
  END IF
END DO

!     Compute degree-6 approximation to the integral of
!     |F - INTEGRAL/AREA|.

dresc = ABS(fv(1) - res6) * w60
kount = 2
DO j = 1, 4
  dresc = dresc + (ABS(fv(kount) - res6) + ABS(fv(kount+1) - res6)  &
  + ABS(fv(kount+2) - res6)) * w(j)
  kount = kount + 3
END DO

!     Compute degree-6 and degree-8 approximations to the integral and
!     its error.

res6 = res6 * area
res8 = res8 * area
resab6 = resab6 * area
dresc = dresc * area
est = ABS(res6 - res8)
IF (dresc /= zero) est = MAX(est, dresc * MIN(one, (twenty * est  &
                                  / dresc)**onept5))
IF (resab6 > uflow) est = MAX(emach * resab6, est)
RETURN
END SUBROUTINE lqm0



SUBROUTINE hinitu(nmax, nwds, n, t)

!     Purpose:
!       Initialize the heap programs with T(1) pointing to the top of
!       the heap.   It is called once, at the start.

!     Input:
!     NMAX = Maximum number of nodes.
!     NWDS = Number of words per node.

!     Output:
!     N = 0 (to hold current number of nodes in heap)
!     T = integer array of pointers to the heap nodes.

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nmax, nwds
INTEGER, INTENT(OUT) :: n, t(:)

!     Local variable.

INTEGER :: i

DO i = 1, nmax
  t(i) = (i-1)*nwds + 1
END DO
n = 0
RETURN
END SUBROUTINE hinitu



SUBROUTINE hinitd(nmax, nwds, n, t)

!     Purpose:
!       Initialize the heap programs with T(NMAX) pointing to the bottom
!       of the heap.   It is called once, at the start.

!     Input:
!     NMAX = Maximum number of nodes.
!     NWDS = Number of words per node.

!     Output:
!     N = 0 (to hold current number of nodes in heap)
!     T = integer array of pointers to the heap nodes.

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nmax, nwds
INTEGER, INTENT(OUT) :: n, t(:)

!     Local variable.

INTEGER :: i

DO i = 1, nmax
  t(i) = (nmax-i)*nwds + 1
END DO
n = 0
RETURN
END SUBROUTINE hinitd



SUBROUTINE hpins(nmax, nwds, DATA, n, t, xnode, hpfun)

!     Purpose:
!       Insert a node into a heap.

!     Input:
!     NMAX = Maximum number of nodes.
!     NWDS = Number of words per node.
!     DATA = Work area for storing nodes.
!     N    = Current number of nodes in the tree.
!     T    = Integer array of pointers to the nodes.
!     XNODE = Real array of length NWDS containing the information to
!             be inserted.
!     HPFUN = Name of user-written LOGICAL function to determine the top
!             node. (GREATR in TWODQD)

!     Output:
!     DATA = Work area with the node inserted.
!     N    = Updated number of nodes.
!     T    = Updated pointer array.

IMPLICIT NONE
INTEGER, INTENT(IN)       :: nmax, nwds
INTEGER, INTENT(IN OUT)   :: n, t(:)
REAL (dp), INTENT(IN)     :: xnode(:)
REAL (dp), INTENT(IN OUT) :: data(:)

INTERFACE
  FUNCTION hpfun(a, b) RESULT(yes)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: a, b
    LOGICAL               :: yes
  END FUNCTION hpfun
END INTERFACE

!     Local variables

INTEGER ::   i, j, ipj, jr, j2, jl

IF (n == nmax) RETURN
n = n + 1
ipj = t(n)
DO i = 1, nwds
  DATA(ipj) = xnode(i)
  ipj = ipj + 1
END DO

j = n
20 IF (j == 1) RETURN
jr = t(j)
j2 = j / 2
jl = t(j2)
IF (hpfun(DATA(jl), DATA(jr))) RETURN
t(j2) = t(j)
t(j) = jl
j = j2
GO TO 20
END SUBROUTINE hpins



SUBROUTINE hpbld(nmax, DATA, n, t, hpfun)

!     Purpose:
!       Builds a heap in T from an array of N elements in DATA spaced
!       NWDS words apart.

!     Input:
!     NMAX = Maximum number of nodes.
!     DATA = Work area for storing nodes.
!     N    = Current number of nodes in the tree.
!     T    = Integer array of pointers to the nodes.
!     HPFUN = Name of user-written LOGICAL function to determine the top
!             node. (GREATR in TWODQD)

!     Output:
!     DATA = Work area with the node inserted.
!     T    = Updated pointer array, with T(1) pointing to the top.

IMPLICIT NONE
INTEGER, INTENT(IN)       :: nmax, n
INTEGER, INTENT(IN OUT)   :: t(:)
REAL (dp), INTENT(IN OUT) :: DATA(:)

INTERFACE
  FUNCTION hpfun(a, b) RESULT(yes)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: a, b
    LOGICAL               :: yes
  END FUNCTION hpfun
END INTERFACE

!     Local variable

INTEGER :: indx

IF (nmax < n) RETURN
indx = n / 2
DO
  IF (indx == 0) RETURN
  CALL hpgro(nmax, DATA, n, t, hpfun, indx)
  indx = indx - 1
END DO
RETURN
END SUBROUTINE hpbld



SUBROUTINE hpdel(nmax, DATA, n, t, hpfun, k)

!     Purpose:
!       Delete the K-th element from the heap.

!     Input:
!     NMAX = Maximum number of nodes.
!     DATA = Work area for storing nodes.
!     N    = Current number of nodes in the tree.
!     T    = Integer array of pointers to the nodes.
!     HPFUN = Name of user-written LOGICAL function to determine the top
!             node. (GREATR in TWODQD)
!     K    = Index of the node to be deleted.

!     Output:
!     DATA = Work area with the node inserted.
!     T    = Updated pointer array, with T(1) pointing to the top.

IMPLICIT NONE
INTEGER, INTENT(IN)       :: nmax, k
INTEGER, INTENT(IN OUT)   :: n, t(:)
REAL (dp), INTENT(IN OUT) :: DATA(:)

INTERFACE
  FUNCTION hpfun(a, b) RESULT(yes)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: a, b
    LOGICAL               :: yes
  END FUNCTION hpfun
END INTERFACE

!     Local variables

INTEGER :: kdel, junk, khalve, il, ir

IF (n == 0) RETURN
IF (k == n) THEN
  n = n - 1
  RETURN
END IF

kdel = k
junk = t(kdel)
t(kdel) = t(n)
t(n) = junk
n = n - 1
10 IF (kdel == 1) THEN
  CALL hpgro(nmax, DATA, n, t, hpfun, kdel)
  RETURN
ELSE
  khalve = kdel / 2
  il = t(khalve)
  ir = t(kdel)
  IF (hpfun(DATA(il), DATA(ir))) THEN
    CALL hpgro(nmax, DATA, n, t, hpfun, kdel)
    RETURN
  ELSE
    t(khalve) = ir
    t(kdel) = il
    kdel = khalve
  END IF
END IF
GO TO 10
END SUBROUTINE hpdel



SUBROUTINE hpgro(nmax, DATA, n, t, hpfun, i)

!     Purpose:
!       Forms a heap out of a tree used privately by HPBLD.
!       The top of the tree is stored in location T(I), the first son
!       is in location T(2I), next son in T(2I+1), etc.   Each branch
!       is assumed to be a tree.

!     Input:
!     NMAX = Maximum number of nodes.
!     DATA = Work area for storing nodes.
!     N    = Current number of nodes in the tree.
!     T    = Integer array of pointers to the nodes.
!     HPFUN = Name of user-written LOGICAL function to determine the top
!             node. (GREATR in TWODQD)
!     I    = Pointer to top of current tree.

!     Output:
!     DATA = Work area with the node inserted.
!     T    = Updated pointer array, with T(1) pointing to the top.

IMPLICIT NONE
INTEGER, INTENT(IN)       :: nmax, n, i
INTEGER, INTENT(IN OUT)   :: t(:)
REAL (dp), INTENT(IN OUT) :: DATA(:)

INTERFACE
  FUNCTION hpfun(a, b) RESULT(yes)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: a, b
    LOGICAL               :: yes
  END FUNCTION hpfun
END INTERFACE

!     Local variables

INTEGER :: j, k, ir, il, itemp

IF (n > nmax) RETURN

k = i
10 j = 2*k

!     Test whether J-th element is a leaf.

IF (j > n) RETURN

!     If there is more than one son, find the smallest.

IF (j == n) GO TO 20
ir = t(j)
il = t(j+1)
IF (hpfun(DATA(il), DATA(ir))) j = j + 1

!     If a son is larger than father, interchange.

20 il = t(k)
ir = t(j)
IF (hpfun(DATA(il), DATA(ir))) RETURN
itemp = t(j)
t(j) = t(k)
t(k) = itemp
k = j
GO TO 10

END SUBROUTINE hpgro



SUBROUTINE hpacc(nmax, nwds, DATA, n, t, xnode, k)

!     Purpose:
!       Access information for the K-th node of the heap.

!     Input:
!     NMAX = Maximum number of nodes.
!     NWDS = Number of words per node.
!     DATA = Work area for storing nodes.
!     N    = Current number of nodes in the tree.
!     T    = Integer array of pointers to the nodes.
!     K    = The node for which information is required.

!     Output:
!     XNODE = Real array of length NWDS in which the information will
!             be returned.

IMPLICIT NONE
INTEGER, INTENT(IN)    :: nmax, nwds, n, t(:), k
REAL (dp), INTENT(IN)  :: DATA(:)
REAL (dp), INTENT(OUT) :: xnode(:)

!     Local variables

INTEGER :: ipj, i

IF (k < 1 .OR. k > n .OR. n > nmax) RETURN
ipj = t(k)
DO i = 1, nwds
  xnode(i) = DATA(ipj)
  ipj = ipj + 1
END DO
RETURN
END SUBROUTINE hpacc



FUNCTION greatr(a, b) RESULT(yes)
IMPLICIT NONE
REAL (dp), INTENT(IN) :: a, b
LOGICAL               :: yes

yes = (a > b)
RETURN
END FUNCTION greatr

END MODULE twod_quad



PROGRAM test_twodqd
! Test twodqd by evaluating:
!        1   1
!       Int Int exp(-x^2.y^2) dy dx = 0.9059444..
!        0   0

USE twod_quad
IMPLICIT NONE

INTEGER, PARAMETER   :: n = 2, maxtri = 500, mevals = 20000
INTEGER              :: iclose, nu, nd, iflag, nevals
REAL (dp)            :: x(3,n), y(3,n), tol = 1.D-04, result, error
REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp

INTERFACE
  FUNCTION f(x, y) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: x, y
    REAL (dp)             :: fn_val
  END FUNCTION f
END INTERFACE

! The square is divided into 2 triangles.
! 2nd subscript is the triangle number.
                                        !      (2,2)  (1,2)
x(1,1) = zero                           !   (3,1)+-----+
x(2,1) = one                            !        |\    |
x(3,1) = zero                           !        | \ 2 |
y(1,1) = zero                           !        |  \  |
y(2,1) = zero                           !        | 1 \ |
y(3,1) = one                            !        |    \|
x(:,2) = one - x(:,1)                   !   (1,1)+-----+(3,2)
y(:,2) = one - y(:,1)                   !            (2,1)

WRITE(*, *) 'Integral should be 0.905940476..'
WRITE(*, *)

DO iclose = 0, 1
  WRITE(*, '(" Using ICLOSE = ", i2)') iclose
  iflag = 0
  nu = 0
  nd = 0
  CALL twodqd(f, n, x, y, tol, iclose, maxtri, mevals, result,  &
                  error, nu, nd, nevals, iflag)
  WRITE(*, '(a, i2)') ' IFLAG = ', iflag
  WRITE(*, '(" Result = ", f11.7, "  Estimated error = ", g11.3)')  &
                                                         result, error
  WRITE(*, '(" No. of function evaluations = ", i5)') nevals
  WRITE(*, '(" NU, ND = ", 2i6)') nu, nd
  WRITE(*, *)
END DO

STOP
END PROGRAM test_twodqd



FUNCTION f(x, y) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

fn_val = EXP(- (x*y)**2 )
RETURN
END FUNCTION f
