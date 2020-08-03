MODULE initialize
! Module to initialize SPNRNT & MVNFNC; replaces the use of ENTRY statements
IMPLICIT NONE

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, PARAMETER  :: nl = 64, nd = 3, nl2 = 50
REAL (dp), SAVE     :: aa(nl), bb(nl), u(nl,nl), cov((nl2*(nl2+1))/2), d1, e1
INTEGER, SAVE       :: infi(nl)

! Limit beyond which the tail area will be taken as zero,
! for the univariate normal.
REAL (dp), SAVE     :: xinfin

END MODULE initialize



MODULE mvnpack

! This is Alan Genz's MVNPACK package of routines for multivariate normal
! integrals, translated to Fortran 90 style by Alan Miller

! 5 March 2001 Changed NCVSRT & LIMITS to handle differences in probabilities
!              in upper tails to be handled much more accurately, and hence
!              for occasional divisions by zero to be eliminated.
! Latest revision - 10 March 2001

USE initialize
IMPLICIT NONE
REAL (dp), PARAMETER, PRIVATE  :: zero = 0._dp, one = 1._dp

CONTAINS


SUBROUTINE ranmvn(n, lower, upper, infin, correl, maxpts,  &
                  abseps, releps, error, value, inform)

!     A subroutine for computing multivariate normal probabilities.
!     This subroutine uses the Monte-Carlo algorithm given in the paper
!     "Numerical Computation of Multivariate Normal Probabilities", in
!     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113
!          Email : alangenz@wsu.edu

!  Parameters

!     N      INTEGER, the number of variables.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, array of correlation coefficients; the correlation
!            coefficient in row I column J of the correlation matrix
!            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
!     MAXPTS INTEGER, maximum number of function values allowed. This
!            parameter can be used to limit the time taken. A
!            sensible strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS
!                           function vaules used; increase MAXPTS to
!                           decrease ERROR;
!            if INFORM = 2, N > 100 or N < 1.

INTEGER, INTENT(IN)    :: n, infin(:), maxpts
INTEGER, INTENT(OUT)   :: inform
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:), abseps, releps
REAL (dp), INTENT(OUT) :: error, value

! Local variables
INTEGER   :: mpt, infis, ivls
REAL (dp) :: d, e, eps

IF ( n > 100 .OR. n < 1 ) THEN
  inform = 2
  value = zero
  error = one
  RETURN
END IF
CALL mvnnit(n, correl, lower, upper, infin, infis, d, e)
inform = 0
IF ( n-infis == 0 ) THEN
  value = one
  error = zero
ELSE IF ( n-infis == 1 ) THEN
  value = e - d
  error = 2E-16
ELSE
  
!        Call then Monte-Carlo integration subroutine
  
  mpt = 25 + 10*n
  CALL rcrude(n-infis-1, mpt, mvnfnc, error, value, 0)
  ivls = mpt
  10 eps = MAX( abseps, releps*ABS(value) )
  IF ( error > eps .AND. ivls < maxpts ) THEN
    mpt = MAX( MIN( INT(mpt*(error/(eps))**2),maxpts-ivls ), 10 )
    CALL rcrude(n-infis-1, mpt, mvnfnc, error, value, 1)
    ivls = ivls + mpt
    GO TO 10
  END IF
  IF ( error > eps .AND. ivls >= maxpts ) inform = 1
END IF
RETURN
END SUBROUTINE ranmvn


SUBROUTINE rcrude(ndim, maxpts, functn, absest, finest, ir)

!     Crude Monte-Carlo Algorithm with simple antithetic variates
!      and weighted results on restart

INTEGER, INTENT(IN)    :: ndim, maxpts, ir
REAL (dp), INTENT(OUT) :: finest, absest

! Local variables
INTEGER         :: m, k, npts
REAL (dp)       :: x(100), fun, varsqr, varprd, findif, finval
REAL (dp), SAVE :: varest

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

IF ( ir <= 0 ) THEN
  varest = zero
  finest = zero
END IF
finval = zero
varsqr = zero
npts = maxpts/2
DO m = 1,npts
  DO k = 1,ndim
    x(k) = uni()
  END DO
  fun = functn(ndim, x)
  x(1:ndim) = one - x(1:ndim)
  fun = ( functn(ndim, x) + fun )/2
  findif = ( fun - finval )/m
  varsqr = ( m - 2 )*varsqr/m + findif**2
  finval = finval + findif
END DO
varprd = varest*varsqr
finest = finest + ( finval - finest )/(one + varprd)
IF ( varsqr > zero ) varest = (one + varprd)/varsqr
absest = 3*SQRT( varsqr/(one + varprd ) )

RETURN
END SUBROUTINE rcrude


SUBROUTINE sphmvn(n, lower, upper, infin, correl, maxpts,  &
                  abseps, releps, error, value, inform)

!     A subroutine for computing multivariate normal probabilities.
!     This subroutine uses a Monte-Carlo algorithm given in the paper
!       "Three Digit Accurate Multiple Normal Probabilities",
!          pp. 369-380, Numer. Math. 35(1980), by I. Deak


!  Parameters

!     N      INTEGER, the number of variables.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, array of correlation coefficients; the correlation
!            coefficient in row I column J of the correlation matrix
!            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
!     MAXPTS INTEGER, maximum number of function values allowed. This
!            parameter can be used to limit the time. A sensible
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL, estimated absolute error, with 99% confidence level.
!     VALUE  REAL, estimated value for the integral
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS
!                           function vaules used; increase MAXPTS to
!                           decrease ERROR;
!            if INFORM = 2, N > 50.

INTEGER, INTENT(IN)    :: n, infin(:), maxpts
INTEGER, INTENT(OUT)   :: inform
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:), abseps, releps
REAL (dp), INTENT(OUT) :: error, value

! Local variables
INTEGER   :: infis, mpt, ns, ivls
REAL (dp) :: d, e, eps

IF ( n > 64 ) THEN
  inform = 2
  value = zero
  error = one
  RETURN
END IF

CALL spnrnt(n, correl, lower, upper, infin, infis, d, e, ns)
inform = 0
IF ( n-infis == 0 ) THEN
  value = one
  error = zero
ELSE IF ( n-infis == 1 ) THEN
  value = e - d
  error = 2E-16
ELSE
  
!        Call then Monte-Carlo integration subroutine
  
  mpt = 25 + ns/n**3
  CALL scrude( n-infis, mpt, error, value, 0 )
  ivls = mpt*ns
  10 eps = MAX( abseps, releps*ABS(value) )
  IF ( error > eps .AND. ivls < maxpts ) THEN
    mpt = MAX( MIN( INT(mpt*(error/(eps))**2),( maxpts - ivls )/ns ), 10 )
    CALL scrude( n-infis, mpt, error, value, 1 )
    ivls = ivls + mpt*ns
    GO TO 10
  END IF
  IF ( error > eps .AND. ivls >= maxpts ) inform = 1
END IF

RETURN
END SUBROUTINE sphmvn


FUNCTION spnrml(n) RESULT(fn_val)

!     Integrand subroutine

USE initialize
INTEGER, INTENT(IN) :: n
REAL (dp)           :: fn_val

! Local variables

INTEGER    :: i, j, k, ns
INTEGER    :: is(nl), ic(nl)
REAL (dp)  :: rs, tmp, bt, y(n)

!    First generate U = COV*(random orthogonal matrix)

DO k = n-1, 1, -1
  tmp = zero
  DO j = k, n
    y(j) = rnor()
    tmp = tmp + y(j)**2
  END DO
  tmp = -SQRT(tmp)
  bt = 1/( tmp*( y(k) + tmp ) )
  y(k) = y(k) + tmp
  DO i = 1, n
    tmp = bt * DOT_PRODUCT( u(i,k:n), y(k:n) )
    u(i,k:n) = u(i,k:n) - tmp * y(k:n)
  END DO
END DO

!     Compute integrand average

rs = SQRT( DBLE(nd) )
DO i = 1,nd
  ic(i) = i
END DO
ic(nd+1) = n+1
fn_val = zero
ns = 0

10 is(1:nd) = -1
20 DO i = 1, n
  tmp = zero
  DO j = 1,nd
    tmp = tmp + is(j)*u(i,ic(j))
  END DO
  y(i) = tmp/rs
END DO

ns = ns + 1
fn_val = fn_val + ( sphlim( n, aa, bb, infi, y ) - fn_val )/ns
DO i = 1, nd
  is(i) = is(i) + 2
  IF ( is(i) < 2 ) GO TO 20
  is(i) = -1
END DO
DO i = 1, nd
  ic(i) = ic(i) + 1
  IF ( ic(i) < ic(i+1)  ) GO TO 10
  ic(i) = i
END DO
fn_val = fn_val/2

RETURN
END FUNCTION spnrml


SUBROUTINE spnrnt( n, correl, lower, upper, infin, infis, d, e, nso )

!     Initialisation

USE initialize
INTEGER, INTENT(IN)    :: n, infin(:)
INTEGER, INTENT(OUT)   :: infis, nso
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:)
REAL (dp), INTENT(OUT) :: d, e

! Local variables
INTEGER   :: i, ij, j, k
REAL (dp) :: tmp

nso = 1
infis = 0
DO i = 1, nd
  nso = 2*nso*( n - infis - i + 1 )/i
END DO

ij = 0
DO i = 1, n
  infi(i) = infin(i)
  IF ( infi(i) < 0 ) THEN
    infis = infis + 1
  ELSE
    aa(i) = zero
    bb(i) = zero
    IF ( infi(i) /= 0 ) aa(i) = lower(i)
    IF ( infi(i) /= 1 ) bb(i) = upper(i)
  END IF
  DO j = 1, i-1
    ij = ij + 1
    u(i,j) = correl(ij)
    u(j,i) = zero
  END DO
  u(i,i) = one
END DO

!     First move any doubly infinite limits to innermost positions

IF ( infis < n ) THEN
  DO i = n, n-infis+1, -1
    IF ( infi(i) >= 0 ) THEN
      DO j = 1,i-1
        IF ( infi(j) < 0 ) THEN
          DO k = 1, j-1
            tmp = u(j,k)
            u(j,k) = u(i,k)
            u(i,k) = tmp
          END DO
          DO k = j+1, i-1
            tmp = u(i,k)
            u(i,k) = u(k,j)
            u(k,j) = tmp
          END DO
          DO k = i+1, n
            tmp = u(k,j)
            u(k,j) = u(k,i)
            u(k,i) = tmp
          END DO
          tmp = aa(j)
          aa(j) = aa(i)
          aa(i) = tmp
          tmp = bb(j)
          bb(j) = bb(i)
          bb(i) = tmp
          tmp = infi(j)
          infi(j) = infi(i)
          infi(i) = tmp
          GO TO 30
        END IF
      END DO
    END IF
  30 END DO
END IF

!     Determine Cholesky decomposition

DO j = 1, n-infis
  DO i = j, n-infis
    tmp = u(i,j)
    DO k = 1, j-1
      tmp = tmp - u(i,k)*u(j,k)
    END DO
    IF ( i == j ) THEN
      u(j,j) = SQRT( MAX( tmp, zero ) )
    ELSE IF ( u(i,i) > zero ) THEN
      u(i,j) = tmp/u(j,j)
    ELSE
      u(i,j) = zero
    END IF
  END DO
END DO

DO i = 1, n-infis
  IF ( u(i,i) > zero ) THEN
    IF ( infi(i) /= 0 ) aa(i) = aa(i)/u(i,i)
    IF ( infi(i) /= 1 ) bb(i) = bb(i)/u(i,i)
    DO j = 1,i
      u(i,j) = u(i,j)/u(i,i)
    END DO
  END IF
END DO
CALL limits( aa(1), bb(1), infi(1), d, e )

RETURN
END SUBROUTINE spnrnt



FUNCTION sphlim( n, a, b, infi, y ) RESULT(fn_val)

REAL (dp), INTENT(IN) :: a(:), b(:), y(:)
INTEGER, INTENT(IN)   :: infi(:), n
REAL (dp)             :: fn_val

! Local variables
INTEGER   :: i
REAL (dp) :: cmn, cmx

cmn = -10*n
cmx =  10*n
DO i = 1,n
  IF ( y(i) > zero ) THEN
    IF ( infi(i) /= 1 ) cmx = MIN( cmx, b(i)/y(i) )
    IF ( infi(i) /= 0 ) cmn = MAX( cmn, a(i)/y(i) )
  ELSE
    IF ( infi(i) /= 1 ) cmn = MAX( cmn, b(i)/y(i) )
    IF ( infi(i) /= 0 ) cmx = MIN( cmx, a(i)/y(i) )
  END IF
END DO
IF ( cmn < cmx ) THEN
  IF ( cmn >= zero .AND. cmx >= zero ) THEN
    fn_val = sphinc( n,  cmx ) - sphinc( n,  cmn )
  ELSE If ( cmn < zero .AND. cmx >= zero ) THEN
    fn_val = sphinc( n, -cmn ) + sphinc( n,  cmx )
  ELSE
    fn_val = sphinc( n, -cmn ) - sphinc( n, -cmx )
  END IF
ELSE
  fn_val = zero
END IF

RETURN
END FUNCTION sphlim


SUBROUTINE scrude( ndim, maxpts, absest, finest, ir )

!     Crude Monte-Carlo Algorithm for Deak method with
!      weighted results on restart

INTEGER, INTENT(IN)    :: ndim, maxpts, ir
REAL (dp), INTENT(OUT) :: finest, absest

! Local variables
INTEGER         :: m
REAL (dp)       :: varsqr, varprd, findif, finval
REAL (dp), SAVE :: varest

IF ( ir <= 0 ) THEN
  varest = zero
  finest = zero
END IF
finval = zero
varsqr = zero
DO m = 1,maxpts
  findif = ( spnrml(ndim) - finval )/m
  finval = finval + findif
  varsqr = ( m - 2 )*varsqr/m + findif**2
END DO
varprd = varest*varsqr
finest = finest + ( finval - finest )/(one + varprd)
IF ( varsqr > zero ) varest = (one + varprd)/varsqr
absest = 3._dp*SQRT( varsqr/( one + varprd ) )

RETURN
END SUBROUTINE scrude


FUNCTION sphinc( n, r ) RESULT(fn_val)

!                   R
!     SPHINC =  K  I  exp(-t*t/2) t**(N-1) dt, for N > 1.
!                N  0

INTEGER, INTENT(IN)   :: n
REAL (dp), INTENT(IN) :: r
REAL (dp)             :: fn_val

! Local variables

INTEGER              :: i
REAL (dp)            :: rr, pf
REAL (dp), PARAMETER :: rp = 2.5066282746310004D0

IF ( r > zero ) THEN
  IF ( r < 5*n ) THEN
    rr = r*r
    pf = one
    DO i = n-2, 2, -2
      pf = one + rr*pf/i
    END DO
    IF ( MOD( n, 2 ) == 0 ) THEN
      fn_val = one - EXP( LOG(pf) - rr/2._dp )
    ELSE
      fn_val = one  - 2._dp*( phi(-r) + EXP( LOG(r*pf) - rr/2._dp )/rp )
    END IF
  ELSE
    fn_val = one
  END IF
ELSE
  fn_val = zero
END IF

RETURN
END FUNCTION sphinc


FUNCTION rnor() RESULT(fn_val)

!       RNOR generates normal random numbers with zero mean and
!       unit standard deviation, often denoted N(0,1).
!       Adapted from RNOR in "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988

REAL  :: fn_val

! Local variables

REAL            :: x, y, s, vt
REAL, PARAMETER :: aa = 12.37586, b = 0.4878992, c = 12.67706, &
                        c1 = 0.9689279, c2 = 1.301198,  &
                        pc = 0.1958303E-1, xn = 2.7769943
INTEGER         :: j
REAL, PARAMETER :: v(65) = (/   &
       0.3409450, 0.4573146, 0.5397793, 0.6062427, 0.6631691,  &
       0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,  &
       0.9490791, 0.9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917,  &
       1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635,  &
       1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929,  &
       1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454,  &
       1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422,  &
       1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166,  &
       1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713,  &
       2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117,  &
       2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943 /)

y = uni()
j = MOD( INT( uni()*128 ), 64 ) + 1

!     Pick sign as Y+Y-1 is positive or negative

vt = v(j+1)
fn_val = (y+y-1)*vt
IF ( ABS(fn_val) > v(j) ) THEN
  x = ( ABS(fn_val)-v(j) )/( vt-v(j) )
  y = uni()
  s = x + y
  IF ( s > c1 ) THEN
    IF ( s > c2 .OR. y > c - aa*EXP(-(b-b*x)**2/2) ) THEN
      fn_val = SIGN( b-b*x, fn_val )
    ELSE
      IF( EXP(-vt**2/2) + y*pc/vt > EXP(-fn_val**2/2) ) THEN
        DO
          x = 0.3601016*LOG( uni() )
          IF ( -2*LOG( uni() ) > x**2 ) EXIT
        END DO
        fn_val = SIGN( xn-x, fn_val )
      END IF
    END IF
  END IF
END IF

RETURN
END FUNCTION rnor


SUBROUTINE sadmvn(n, lower, upper, infin, correl, maxpts,  &
                  abseps, releps, error, value, inform)

!     A subroutine for computing multivariate normal probabilities.
!     This subroutine uses an algorithm given in the paper
!     "Numerical Computation of Multivariate Normal Probabilities", in
!     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113
!          Email : alangenz@wsu.edu

!  Parameters

!     N      INTEGER, the number of variables.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, array of correlation coefficients; the correlation
!            coefficient in row I column J of the correlation matrix
!            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
!     MAXPTS INTEGER, maximum number of function values allowed. This
!            parameter can be used to limit the time taken. A
!            sensible strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS
!                           function vaules used; increase MAXPTS to
!                           decrease ERROR;
!            if INFORM = 2, N > 20 or N < 1.

INTEGER, INTENT(IN)    :: n, infin(:), maxpts
INTEGER, INTENT(OUT)   :: inform
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:), abseps, releps
REAL (dp), INTENT(OUT) :: error, value

! Local variables
INTEGER   :: infis, m, rulcls, totcls, newcls, maxcls
REAL (dp) :: oldval, d, e

IF ( n > 20 .OR. n < 1 ) THEN
  inform = 2
  value = zero
  error = one
  RETURN
END IF

CALL mvnnit(n, correl, lower, upper, infin, infis, d, e)
inform = 0
m = n - infis
IF ( m == 0 ) THEN
  value = one
  error = zero
ELSE IF ( m == 1 ) THEN
  value = e - d
  error = 2E-16
ELSE
  
!        Call the subregion adaptive integration subroutine
  
  m = m - 1
  rulcls = 1
  CALL adapt(m, rulcls, 0, mvnfnc, abseps, releps, error, value, inform)
  maxcls = MIN( 10*rulcls, maxpts )
  totcls = 0
  CALL adapt(m, totcls, maxcls, mvnfnc, abseps, releps, error, value, inform)
  IF ( error > MAX( abseps, releps*ABS(value) ) ) THEN
    10 oldval = value
    maxcls = MAX( 2*rulcls, MIN( 3*maxcls/2, maxpts - totcls ) )
    newcls = -1
    CALL adapt(m, newcls, maxcls, mvnfnc, abseps, releps, error, value, inform)
    totcls = totcls + newcls
    error = ABS(value-oldval) + SQRT(rulcls*error**2/totcls)
    IF ( error > MAX( abseps, releps*ABS(value) ) ) THEN
      IF ( maxpts - totcls > 2*rulcls ) GO TO 10
    ELSE
      inform = 0
    END IF
  END IF
END IF

RETURN
END SUBROUTINE sadmvn


SUBROUTINE kromvn(n, lower, upper, infin, correl, maxpts,  &
                  abseps, releps, error, value, inform)

!     A subroutine for computing multivariate normal probabilities.
!     This subroutine uses an algorithm given in the paper
!     "Numerical Computation of Multivariate Normal Probabilities", in
!     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113
!          Email : alangenz@wsu.edu

!  Parameters

!     N      INTEGER, the number of variables.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, array of correlation coefficients; the correlation
!            coefficient in row I column J of the correlation matrix
!            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
!     MAXPTS INTEGER, maximum number of function values allowed. This
!            parameter can be used to limit the time. A sensible
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS    REAL absolute error tolerance.
!     RELEPS    REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS
!                           function values used; increase MAXPTS to
!                           decrease ERROR;
!            if INFORM = 2, N > 40 or N < 1.

INTEGER, INTENT(IN)    :: n, infin(:), maxpts
INTEGER, INTENT(OUT)   :: inform
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:), releps, abseps
REAL (dp), INTENT(OUT) :: error, value

! Local variables
INTEGER   :: infis, ivls
REAL (dp) :: e, d

IF ( n > 40 .OR. n < 1 ) THEN
  inform = 2
  value = zero
  error = one
  RETURN
END IF
CALL mvnnit(n, correl, lower, upper, infin, infis, d, e)
inform = 0
IF ( n-infis == 0 ) THEN
  value = one
  error = zero
ELSE IF ( n-infis == 1 ) THEN
  value = e - d
  error = 2E-16
ELSE
  
!        Call the lattice rule integration integration subroutine
  
  ivls = 0
  CALL krobov(n-infis-1, ivls, maxpts, mvnfnc, abseps, releps,  &
              error, value, inform)
END IF

RETURN
END SUBROUTINE kromvn


SUBROUTINE krobov(ndim, minvls, maxvls, functn, abseps, releps,  &
                  abserr, finest, ifail)

!  Automatic Multidimensional Integration Subroutine

!         AUTHOR: Alan Genz
!                 Department of Mathematics
!                 Washington State University
!                 Pulman, WA 99164-3113

!         Last Change: 5/15/93

!  KROBOV computes an approximation to the integral

!      1  1     1
!     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
!      0  0     0


!  KROBOV uses randomized Korobov rules. The primary references are
!  "Randomization of Number Theoretic Methods for Multiple Integration"
!   R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
!  and
!   "Optimal Parameters for Multidimensional Integration",
!    P. Keast, SIAM J Numer Anal, 10, pp.831-838.

!**************  Parameters for KROBOV  *******************************
!***** Input parameters
!  NDIM    Number of variables, must exceed 1, but not exceed 40
!  MINVLS  Integer minimum number of function evaluations allowed.
!          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
!          routine assumes a previous call of KROBOV has been made with
!          the same integrand and continues that calculation.
!  MAXVLS  Integer maximum number of function evaluations allowed.
!  FUNCTN  EXTERNALly declared user defined function to be integrated.
!          It must have parameters (NDIM,Z), where Z is a real array
!          of dimension NDIM.
!  ABSEPS  Required absolute accuracy.
!  RELEPS  Required relative accuracy.
!***** Output parameters
!  MINVLS  Actual number of function evaluations used by KROBOV.
!  ABSERR  Estimated absolute accuracy of FINEST.
!  FINEST  Estimated value of integral.
!  IFAIL   IFAIL = 0 for normal exit, when
!                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
!                  and
!                     INTVLS <= MAXCLS.
!          IFAIL = 1 If MAXVLS was too small for KROBOV to obtain the
!                  required accuracy. In this case KROBOV returns a
!                  value FINEST with estimated absolute accuracy ABSERR.
!***********************************************************************

INTEGER, INTENT(IN)     :: ndim, maxvls
INTEGER, INTENT(IN OUT) :: minvls
INTEGER, INTENT(OUT)    :: ifail
REAL (dp), INTENT(IN)   :: abseps, releps
REAL (dp), INTENT(OUT)  :: finest, abserr

! Local variables
INTEGER            :: lpsamp, i, intvls
INTEGER, SAVE      :: sampls, np
INTEGER, PARAMETER :: plim = 17, nlim = 40, minsmp = 8
REAL (dp)          :: vk(nlim), difint, finval, varsqr, varest, varprd, value

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

INTEGER, PARAMETER :: p(17) = (/ 173, 263, 397, 593, 887, 1327, 1997, 2999,  &
                                 4493, 6737, 10111, 15161, 22751, 34127,  &
                                 51193, 76801, 115183 /)
INTEGER, PARAMETER :: c(17,40) = RESHAPE( (/         &
  73,111,163,229,192,513,839,1148,1360,2602,3071,6280,9438,12543,19397,22500,33681, &
  34,106,168,40,424,195,146,1406,383,2818,3114,2350,3557,7726,10691,16187,26530, &
  57,51,164,268,55,599,860,1192,842,3240,1170,1452,817,6142,20935,33739,42347, &
  9,36,133,240,221,661,183,1094,2157,2528,3432,7009,9167,11635,12665,23268,15491, &
  12,48,23,31,179,443,121,1290,30,2260,2726,4273,7379,6071,5132,23268,15491, &
  2,110,64,119,242,632,11,632,959,3141,1098,4273,811,6071,5132,23268,15491, &
  16,2,2,71,242,251,11,341,3,2857,3371,2538,811,6071,5132,23268,15491, &
  30,2,2,296,2,603,793,785,717,1484,185,3976,6077,6071,5132,23268,15491, &
  30,2,106,130,2,663,998,393,1107,2113,4,4273,6205,1886,5132,23268,24490, &
  42,2,80,199,11,2,2,1499,2,2265,3143,2716,2505,1993,5132,20841,5708, &
  70,70,80,149,11,425,2,2,2,2,5055,2716,2896,634,4347,14733,24490, &
  86,70,126,149,11,425,110,798,2,2,2,6143,4055,836,1087,5916,34535, &
  2,48,16,149,394,603,110,808,836,2207,2,6143,4055,1993,2451,4270,204, &
  53,2,16,149,394,425,236,798,836,2207,2,91,4055,2095,2134,13443,23709, &
  53,2,16,149,439,425,110,918,1134,2207,2,107,900,142,10539,13443,5474, &
  30,70,16,31,394,525,236,393,836,542,2,2716,1333,6542,3595,21793,14302, &
  30,124,16,130,394,412,147,924,836,132,2,2,1154,13195,3892,21793,38969, &
  5,124,16,149,394,412,147,924,426,934,334,2,2,4,6301,871,11298, &
  42,70,107,149,394,412,110,2,898,378,1254,2,2,4,3595,18819,2, &
  42,48,80,79,439,412,190,2,898,378,4146,2,2,2,2,2,22095, &
  70,48,2,119,394,412,147,2,65,2099,617,2,2,2,2,2,7317, &
  42,48,2,119,394,82,147,1499,836,934,1879,2,2,2,2,2,2, &
  53,48,2,31,394,82,147,2,836,225,2,2,2,2,2,2,2, &
  42,108,32,82,101,82,147,1016,836,225,2,2,2,2,2,2,2, &
  42,65,32,130,378,603,147,798,836,225,1146,2,2,2,2,2,2, &
  53,48,32,122,394,580,147,798,216,169,475,2,2,2,2,2,2, &
  42,48,31,122,394,580,236,798,104,378,4725,2,2,2,2,2,2, &
  53,70,64,122,394,444,110,808,300,2257,2,2,2,2,2,2,2, &
  53,2,31,122,394,82,110,270,836,2257,2,2,2,2,2,2,2, &
  2,20,31,2,394,82,147,1344,1022,2257,475,2,2,2,2,2,2, &
  86,2,4,130,202,276,110,798,1022,2257,475,2,2,2,2,2,2, &
  2,2,4,130,279,601,110,798,1022,934,475,2,2,2,2,2,2, &
  2,2,4,130,394,276,632,798,1022,2576,475,2,2,2,2,2,2, &
  2,2,126,130,279,276,147,798,1420,934,475,2,2,2,2,2,2, &
  2,2,16,2,2,276,147,798,1478,934,638,2,2,2,2,2,2, &
  2,2,16,2,2,276,148,798,1478,934,638,2,2,2,2,2,2, &
  2,2,16,82,2,112,2,798,1478,934,638,2,2,2,2,2,2, &
  2,2,16,82,2,112,147,2,283,934,2,2,2,2,2,2,2, &
  2,2,16,2,2,112,147,2,2246,2257,3107,2,2,2,2,2,2, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /), (/ 17, 40 /) )

ifail = 1
intvls = 0
finest = zero
varest = zero
IF ( minvls >= 0 ) THEN
  finest = zero
  sampls = minsmp
  DO i = 1,plim
    np = i
    IF ( minvls < sampls*p(i) ) GO TO 10
  END DO
  sampls = MAX( minsmp, minvls/p(np) )
END IF
10 vk(1) = one/p(np)
DO i = 2,ndim
  vk(i) = MOD( c(np,ndim-1)*vk(i-1), one )
END DO
finval = zero
varsqr = zero
lpsamp = sampls/2
DO i = 1, lpsamp
  CALL krosum( ndim, value, p(np), vk, functn)
  difint = ( value - finval )/i
  finval = finval + difint
  varsqr = ( i - 2 )*varsqr/i + difint**2
END DO
intvls = intvls + 2*lpsamp*p(np)
varprd = varest*varsqr
finest = finest + ( finval - finest )/(one + varprd )
IF ( varsqr > zero ) varest = (one + varprd )/varsqr
abserr = 3.*SQRT( varsqr/(one + varprd ) )
IF ( abserr > MAX( abseps, ABS(finest)*releps ) ) THEN
  IF ( np < plim ) THEN
    np = np + 1
  ELSE
    sampls = MAX( minsmp,MIN(3*sampls/2,(maxvls-intvls)/p(np)) )
  END IF
  IF ( intvls+sampls*p(np) <= maxvls ) GO TO 10
ELSE
  ifail = 0
END IF

RETURN
END SUBROUTINE krobov


SUBROUTINE krosum( ndim, sumkro, npts, vk, functn)

REAL (dp), INTENT(IN)  :: vk(:)
INTEGER, INTENT(IN)    :: ndim, npts
REAL (dp), INTENT(OUT) :: sumkro

! Local variables
INTEGER   :: k, j
REAL (dp) :: sumfun, alpha(ndim), x(ndim)

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

sumfun = zero
DO j = 1,ndim
  alpha(j) = uni()
END DO
DO k = 1,npts
  DO j = 1,ndim
    x(j) = ABS( 2*MOD( alpha(j) + vk(j)*k, one ) - 1 )
  END DO
  sumfun = sumfun + functn(ndim,x)/2
  DO j = 1,ndim
    x(j) = 1 - x(j)
  END DO
  sumfun = sumfun + functn(ndim,x)/2
END DO
sumkro = sumfun/npts

RETURN
END SUBROUTINE krosum


FUNCTION mvnfnc(n, w) RESULT(fn_val)

!     Integrand subroutine

USE initialize
INTEGER, INTENT(IN)   :: n
REAL (dp), INTENT(IN) :: w(:)
REAL (dp)             :: fn_val

! Local variables
INTEGER    :: ij, i, j
REAL (dp)  :: prod, di, ei, emd, sum, y(n)

di = d1
ei = e1
prod = e1 - d1
emd = prod
ij = 1
DO i = 1,n
  y(i) = phinv( di + w(i)*emd )
  sum = zero
  DO j = 1,i
    ij = ij + 1
    sum = sum + cov(ij)*y(j)
  END DO
  ij = ij + 1
  IF ( cov(ij) > zero ) THEN
    CALL limits( aa(i+1)-sum, bb(i+1)-sum, infi(i+1), di, ei, emd )
  ELSE
    di = ( one + SIGN( one, aa(i+1)-sum ) )/2
    ei = ( one + SIGN( one, bb(i+1)-sum ) )/2
  END IF
  prod = prod*emd
END DO
fn_val = prod

RETURN
END FUNCTION mvnfnc


SUBROUTINE mvnnit(n, correl, lower, upper, infin, infis, d, e)

!     Initialization and computation of covariance Cholesky factor.

USE initialize
INTEGER, INTENT(IN)    :: n, infin(:)
INTEGER, INTENT(OUT)   :: infis
REAL (dp), INTENT(IN)  :: correl(:), lower(:), upper(:)
REAL (dp), INTENT(OUT) :: d, e

! Set xinfin (infinity).
! Such that EXP(-2.x^2) > 10^(12) times TINY
xinfin = SQRT( -2.0_dp * LOG(1.e+12 * TINY(one)) )

CALL ncvsrt(n, lower, upper, correl, infin, infis, aa, bb, infi, d, e)
d1 = d
e1 = e
IF ( n - infis == 2 ) THEN
  d = SQRT( one + cov(2)**2 )
  aa(2) = aa(2)/d
  bb(2) = bb(2)/d
  e = bvn( aa, bb, infi, cov(2)/d )
  d = zero
  infis = infis + 1
END IF

RETURN
END SUBROUTINE mvnnit



SUBROUTINE limits( a, b, infin, lower, upper, up_low )
REAL (dp), INTENT(IN)             :: a, b
INTEGER, INTENT(IN)               :: infin
REAL (dp), INTENT(OUT), OPTIONAL  :: lower, upper, up_low

! Local variables
REAL (dp)  :: ap, aq, bp, bq

ap = zero
aq = one
bp = one
bq = zero
IF ( infin >= 0 ) THEN
  IF ( infin /= 0 ) CALL normp(a, ap, aq)
  IF ( infin /= 1 ) CALL normp(b, bp, bq)
END IF

IF (PRESENT(lower)) lower = ap
IF (PRESENT(upper)) upper = bp

! N.B. When the difference is required, use the small tail probabilities.
IF (PRESENT(up_low)) THEN
  SELECT CASE (infin)
    CASE (:-1)
      up_low = one
    CASE (0)
      up_low = bp
    CASE (1)
      up_low = aq
    CASE (2:)
      IF (ap+bp < one) THEN
        up_low = bp - ap
      ELSE
        up_low = aq - bq
      END IF
  END SELECT
END IF

RETURN
END SUBROUTINE limits


SUBROUTINE ncvsrt(n, lower, upper, correl, infin, infis, a, b, infi, d, e)

! N.B. Argument COV has been removed.

!     Subroutine to sort integration limits.

INTEGER, INTENT(IN)     :: n, infin(:)
INTEGER, INTENT(OUT)    :: infi(:)
INTEGER, INTENT(OUT)    :: infis
REAL (dp), INTENT(IN)   :: lower(:), upper(:), correl(:)
REAL (dp), INTENT(OUT)  :: a(:), b(:), d, e

! Local variables
INTEGER   :: i, j, k, ij, ii, jmin
REAL (dp) :: ai, bi, sum, sumsq
REAL (dp) :: cvdiag, amin, bmin, emd, emd_min, yl, yu, y(n)
REAL (dp), PARAMETER :: sqtwpi = 2.50662827463100050240_dp

! Reset INFIN if the upper or lower limits are too extreme.
DO i = 1, n
  infi(i) = infin(i)
  SELECT CASE (infin(i))
    CASE (0)
      IF (upper(i) > xinfin) infi(i) = -1
    CASE (1)
      IF (lower(i) < -xinfin) infi(i) = -1
    CASE (2:)
      IF (lower(i) < -xinfin) infi(i) = 0
      IF (upper(i) > xinfin) infi(i) = MIN( infi(i)-1, 1)
  END SELECT
END DO

ij = 0
ii = 0
infis = 0
DO i = 1,n
  IF ( infi(i) < 0 ) THEN
    infis = infis + 1
  ELSE
    a(i) = zero
    b(i) = zero
    IF ( infi(i) /= 0 ) a(i) = lower(i)
    IF ( infi(i) /= 1 ) b(i) = upper(i)
  END IF
  DO j = 1,i-1
    ij = ij + 1
    ii = ii + 1
    cov(ij) = correl(ii)
  END DO
  ij = ij + 1
  cov(ij) = 1
END DO

!     First move any doubly infinite limits to innermost positions

IF ( infis < n ) THEN
  DO i = n, n-infis+1, -1
    IF ( infi(i) >= 0 ) THEN
      DO j = 1,i-1
        IF ( infi(j) < 0 ) THEN
          CALL rcswap(j, i, a, b, infi, n)
          GO TO 10
        END IF
      END DO
    END IF
  10 END DO

!     Sort remaining limits and determine Cholesky decomposition
  
  ii = 0
  DO i = 1,n-infis
    
!     Determine the integration limits for variable with minimum
!      expected probability and interchange that variable with Ith.
    
    emd_min = one
    jmin = i
    cvdiag = zero
    ij = ii
    DO j = i, n-infis
      sum = zero
      sumsq = zero
      DO k = 1, i-1
        sum = sum + cov(ij+k)*y(k)
        sumsq = sumsq + cov(ij+k)**2
      END DO
      ij = ij + j
      sumsq = SQRT( MAX( cov(ij)-sumsq, zero ) )
      IF ( sumsq > zero ) THEN
        ai = (a(j) - sum)/sumsq
        IF (ABS(ai) > xinfin) CALL adjust_infi(ai, infi(j))
        bi = (b(j) - sum)/sumsq
        IF (ABS(bi) > xinfin) CALL adjust_infi(bi, infi(j))
        CALL limits( ai, bi, infi(j), UP_LOW=emd )
        IF ( emd_min >= emd ) THEN
          jmin = j
          amin = ai
          bmin = bi
          emd_min = emd
          cvdiag = sumsq
        END IF
      END IF
    END DO
    IF ( jmin /= i) CALL rcswap(i, jmin, a, b, infi, n)

!     Compute Ith column of Cholesky factor.
    
    ij = ii + i
    cov(ij) = cvdiag
    DO j = i+1, n-infis
      IF ( cvdiag > zero ) THEN
        sum = cov(ij+i)
        DO k = 1, i-1
          sum = sum - cov(ii+k)*cov(ij+k)
        END DO
        cov(ij+i) = sum/cvdiag
      ELSE
        cov(ij+i) = zero
      END IF
      ij = ij + j
    END DO
    
!     Compute expected value for Ith integration variable and
!     scale Ith covariance matrix row and limits.
    
    IF ( cvdiag > zero ) THEN
      yl = zero
      yu = zero
      IF ( infi(i) /= 0 ) yl = -EXP(-amin**2/2) / sqtwpi
      IF ( infi(i) /= 1 ) yu = -EXP(-bmin**2/2) / sqtwpi
      IF (yu - yl > EPSILON(zero)) THEN
        y(i) = (yu - yl) / emd_min
      ELSE
        y(i) = zero
      END IF
      DO j = 1,i
        ii = ii + 1
        cov(ii) = cov(ii) / cvdiag
      END DO
      IF ( infi(i) /= 0 ) a(i) = a(i) / cvdiag
      IF ( infi(i) /= 1 ) b(i) = b(i) / cvdiag
    ELSE
      y(i) = zero
      ii = ii + i
    END IF
  END DO
  CALL limits( a(1), b(1), infi(1), d, e)
END IF

RETURN
END SUBROUTINE ncvsrt


FUNCTION condit( n, symin ) RESULT(fn_val)

!     Computes condition number of symmetric matix in situ

INTEGER, INTENT(IN)   :: n
REAL (dp), INTENT(IN) :: symin(:)
REAL (dp)             :: fn_val

! Local variables
INTEGER, PARAMETER :: nl = 50
REAL (dp)          :: sum, rowmx, rowmxi, sym(nl*(nl+1)/2)
INTEGER            :: ii, ij, i, j, im

rowmx = zero
ij = 0
DO i = 1,n
  sum = zero
  im = (i-2)*(i-1)/2
  DO j = 1,i-1
    im = im + 1
    sum = sum + ABS(symin(im))
    ij = ij + 1
    sym(ij) = symin(im)
  END DO
  sum = sum + 1
  ij = ij + 1
  sym(ij) = 1
  im = im + i
  DO j = i,n-1
    sum = sum + ABS(symin(im))
    im = im + j
  END DO
  rowmx = MAX( sum, rowmx )
END DO
CALL syminv(n, sym)
rowmxi = zero
ii = 0
DO i = 1,n
  sum = zero
  ij = ii
  DO j = 1,i
    ij = ij + 1
    sum = sum + ABS(sym(ij))
  END DO
  DO j = i,n-1
    ij = ij + j
    sum = sum + ABS(sym(ij))
  END DO
  rowmxi = MAX( sum, rowmxi )
  ii = ii + i
END DO
fn_val = rowmx*rowmxi

RETURN
END FUNCTION condit


SUBROUTINE syminv(n, lowinv, det)

!     Computes lower symmetric inverse and determinant in situ

INTEGER, INTENT(IN)              :: n
REAL (dp), INTENT(IN OUT)        :: lowinv(:)
REAL (dp), INTENT(OUT), OPTIONAL :: det

! Local variables
INTEGER   :: i, ii

CALL cholsk(n, lowinv)
IF (PRESENT(det)) THEN
  det = one
  ii = 0
  DO i = 1,n
    ii = ii + i
    det = det*lowinv(ii)
  END DO
END IF

CALL cholnv(n, lowinv)
CALL cholpi(n, lowinv)

RETURN
END SUBROUTINE syminv


SUBROUTINE cholsk(n, chofac)

!     Computes Choleski factor in situ

INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: chofac(:)

! Local variables
INTEGER   :: i, ii, j, jj, k
REAL (dp) :: s, t

jj = 0
DO j = 1,n
  ii = jj
  DO i = j,n
    s = chofac(ii+j)
    DO k = 1,j-1
      s = s - chofac(ii+k)*chofac(jj+k)
    END DO
    IF ( i == j ) THEN
      t = SQRT( MAX( s, zero ) )
      chofac(ii+j) = t
    ELSE
      chofac(ii+j) = s/t
    END IF
    ii = ii + i
  END DO
  jj = jj + j
END DO

RETURN
END SUBROUTINE cholsk


SUBROUTINE cholnv(n, choinv)

!     Inverts a lower triangular matrix in situ

INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: choinv(:)

! Local variables
INTEGER   :: i, ii, j, jj, k, kk
REAL (dp) :: s, t


ii = 0
DO i = 1,n
  t = 1/choinv(ii+i)
  jj = 0
  DO j = 1,i-1
    s = zero
    jj = jj + j
    kk = jj
    DO k = j,i-1
      s = s + choinv(ii+k)*choinv(kk)
      kk = kk + k
    END DO
    choinv(ii+j) = -s*t
  END DO
  ii = ii + i
  choinv(ii) = t
END DO

RETURN
END SUBROUTINE cholnv


SUBROUTINE cholpi(n, chopdi)

!     Multiplies Choleski inverse factors in situ

INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: chopdi(:)

! Local variables
INTEGER   :: i, ii, j, jj, k, kk
REAL (dp) :: s

ii = 0
DO i = 1,n
  DO j = 1,i
    s = zero
    jj = ii + i
    kk = ii + j
    DO k = i,n
      s = s + chopdi(kk)*chopdi(jj)
      jj = jj + k
      kk = kk + k
    END DO
    chopdi(ii+j) = s
  END DO
  ii = ii + i
END DO

RETURN
END SUBROUTINE cholpi


SUBROUTINE rcswap(p, q, a, b, infin, n)

! N.B. Arguments C has been removed.

!     Swaps rows and columns P and Q in situ.

REAL (dp), INTENT(IN OUT) :: a(:), b(:)
INTEGER, INTENT(IN)       :: p, q, n
INTEGER, INTENT(IN OUT)   :: infin(:)

! Local variables
INTEGER   :: i, j, ii, jj
REAL (dp) :: t

t = a(p)
a(p) = a(q)
a(q) = t
t = b(p)
b(p) = b(q)
b(q) = t
j = infin(p)
infin(p) = infin(q)
infin(q) = j

jj = (p*(p-1))/2
ii = (q*(q-1))/2
t = cov(jj+p)
cov(jj+p) = cov(ii+q)
cov(ii+q) = t
DO j = 1, p-1
  t = cov(jj+j)
  cov(jj+j) = cov(ii+j)
  cov(ii+j) = t
END DO
jj = jj + p
DO i = p+1, q-1
  t = cov(jj+p)
  cov(jj+p) = cov(ii+i)
  cov(ii+i) = t
  jj = jj + i
END DO
ii = ii + q
DO i = q+1, n
  t = cov(ii+p)
  cov(ii+p) = cov(ii+q)
  cov(ii+q) = t
  ii = ii + i
END DO

RETURN
END SUBROUTINE rcswap


SUBROUTINE normp(z, p, q, pdf)

! Normal distribution probabilities accurate to 1.e-15.
! Z = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
! PDF = the probability density.

! Based upon algorithm 5666 for the error function, from:
! Hart, J.F. et al, 'Computer Approximations', Wiley 1968

! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 30 March 1986
! Latest revision of Fortran 90 version - 12 August 1997

IMPLICIT NONE
REAL (dp), INTENT(IN)            :: z
REAL (dp), INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL (dp), PARAMETER :: p0 = 220.2068679123761D0, p1 = 221.2135961699311D0,  &
                        p2 = 112.0792914978709D0, p3 = 33.91286607838300D0,  &
                        p4 = 6.373962203531650D0, p5 = .7003830644436881D0,  &
                        p6 = .3526249659989109D-01,  &
                        q0 = 440.4137358247522D0, q1 = 793.8265125199484D0,  &
                        q2 = 637.3336333788311D0, q3 = 296.5642487796737D0,  &
                        q4 = 86.78073220294608D0, q5 = 16.06417757920695D0,  &
                        q6 = 1.755667163182642D0, q7 = .8838834764831844D-1, &
                        cutoff = 7.071D0, root2pi = 2.506628274631001D0
REAL (dp)            :: zabs, expntl, pp, qq, ppdf

zabs = ABS(z)

! |Z| > 37.

IF (zabs > 37.d0) THEN
  IF (PRESENT(pdf)) pdf = zero
  IF (z > zero) THEN
    IF (PRESENT(p)) p = one
    IF (PRESENT(q)) q = zero
  ELSE
    IF (PRESENT(p)) p = zero
    IF (PRESENT(q)) q = one
  END IF
  RETURN
END IF

! |Z| <= 37.

expntl = EXP(-0.5D0*zabs**2)
ppdf = expntl/root2pi
IF (PRESENT(pdf)) pdf = ppdf

! |Z| < CUTOFF = 10/sqrt(2).

IF (zabs < cutoff) THEN
  pp = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + p2)*zabs     &
                   + p1)*zabs + p0) / (((((((q7*zabs + q6)*zabs + q5)*zabs &
                   + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +q0)
  
! |Z| >= CUTOFF.
  
ELSE
  pp = ppdf/(zabs + one/(zabs + 2.d0/(zabs + 3.d0/(zabs + 4.d0/(zabs + 0.65D0)))))
END IF

IF (z < zero) THEN
  qq = one - pp
ELSE
  qq = pp
  pp = one - qq
END IF

IF (PRESENT(p)) p = pp
IF (PRESENT(q)) q = qq

RETURN
END SUBROUTINE normp



FUNCTION phi(z) RESULT(p)
REAL (dp), INTENT(IN) :: z
REAL (dp)             :: p

CALL normp(z, p)

RETURN
END FUNCTION phi



SUBROUTINE ppnd16(p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE

REAL (dp), INTENT(IN)  :: p
INTEGER, INTENT(OUT)   :: ifault
REAL (dp), INTENT(OUT) :: normal_dev

! Local variables

REAL (dp) :: half = 0.5d0, split1 = 0.425d0, split2 = 5.d0,  &
             const1 = 0.180625d0, const2 = 1.6d0, q, r

! Coefficients for P close to 0.5

REAL (dp) :: a0 = 3.3871328727963666080D0, &
             a1 = 1.3314166789178437745D+2, &
             a2 = 1.9715909503065514427D+3, &
             a3 = 1.3731693765509461125D+4, &
             a4 = 4.5921953931549871457D+4, &
             a5 = 6.7265770927008700853D+4, &
             a6 = 3.3430575583588128105D+4, &
             a7 = 2.5090809287301226727D+3, &
             b1 = 4.2313330701600911252D+1, &
             b2 = 6.8718700749205790830D+2, &
             b3 = 5.3941960214247511077D+3, &
             b4 = 2.1213794301586595867D+4, &
             b5 = 3.9307895800092710610D+4, &
             b6 = 2.8729085735721942674D+4, &
             b7 = 5.2264952788528545610D+3
! HASH SUM AB    55.8831928806149014439

! Coefficients for P not close to 0, 0.5 or 1.

REAL (dp) :: c0 = 1.42343711074968357734D0, &
             c1 = 4.63033784615654529590D0, &
             c2 = 5.76949722146069140550D0, &
             c3 = 3.64784832476320460504D0, &
             c4 = 1.27045825245236838258D0, &
             c5 = 2.41780725177450611770D-1, &
             c6 = 2.27238449892691845833D-2, &
             c7 = 7.74545014278341407640D-4, &
             d1 = 2.05319162663775882187D0, &
             d2 = 1.67638483018380384940D0, &
             d3 = 6.89767334985100004550D-1, &
             d4 = 1.48103976427480074590D-1, &
             d5 = 1.51986665636164571966D-2, &
             d6 = 5.47593808499534494600D-4, &
             d7 = 1.05075007164441684324D-9
! HASH SUM CD    49.33206503301610289036

! Coefficients for P near 0 or 1.

REAL (dp) :: e0 = 6.65790464350110377720D0, &
             e1 = 5.46378491116411436990D0, &
             e2 = 1.78482653991729133580D0, &
             e3 = 2.96560571828504891230D-1, &
             e4 = 2.65321895265761230930D-2, &
             e5 = 1.24266094738807843860D-3, &
             e6 = 2.71155556874348757815D-5, &
             e7 = 2.01033439929228813265D-7, &
             f1 = 5.99832206555887937690D-1, &
             f2 = 1.36929880922735805310D-1, &
             f3 = 1.48753612908506148525D-2, &
             f4 = 7.86869131145613259100D-4, &
             f5 = 1.84631831751005468180D-5, &
             f6 = 1.42151175831644588870D-7, &
             f7 = 2.04426310338993978564D-15
! HASH SUM EF    47.52583317549289671629

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
  ELSE
    r = r - split2
    normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF

RETURN
END SUBROUTINE ppnd16


FUNCTION phinv(p) RESULT(normal_dev)
IMPLICIT NONE
REAL (dp), INTENT(IN) :: p
REAL (dp)             :: normal_dev

! Local variable
INTEGER :: ifault

CALL ppnd16(p, normal_dev, ifault)
IF (ifault /= 0) THEN
  IF (p <= zero) THEN
    normal_dev = -10._dp
  ELSE
    normal_dev = 10._dp
  END IF
END IF

RETURN
END FUNCTION phinv


FUNCTION bvn ( lower, upper, infin, correl ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!  Parameters

!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, correlation coefficient.

REAL (dp), INTENT(IN) :: lower(:), upper(:), correl
INTEGER, INTENT(IN)   :: infin(:)
REAL (dp)             :: fn_val

IF ( infin(1) == 2  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )    &
         - bvnu ( upper(1), lower(2), correl )  &
         - bvnu ( lower(1), upper(2), correl )  &
         + bvnu ( upper(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
         - bvnu ( upper(1), lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
         - bvnu ( lower(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
         - bvnu ( -lower(1), -upper(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
         - bvnu ( -upper(1), -lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( lower(1), -upper(2), -correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( -upper(1), lower(2), -correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )
END IF

RETURN
END FUNCTION bvn


FUNCTION bvnu( sh, sk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!       Yihong Ge
!       Department of Computer Science and Electrical Engineering
!       Washington State University
!       Pullman, WA 99164-2752
!       Email : yge@eecs.wsu.edu
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! BVN - calculate the probability that X is larger than SH and Y is
!       larger than SK.

! Parameters

!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights

REAL (dp), INTENT(IN) :: sh, sk, r
REAL (dp)             :: fn_val

! Local variables
INTEGER              :: i, lg, ng
REAL (dp), PARAMETER :: twopi = 6.283185307179586
REAL (dp)            :: as, a, b, c, d, rs, xs
REAL (dp)            :: bvn, sn, asr, h, k, bs, hs, hk
!     Gauss Legendre Points and Weights, N =  6
! DATA ( w(i,1), x(i,1), i = 1,3) /  &
! 0.1713244923791705D+00,-0.9324695142031522D+00,  &
! 0.3607615730481384D+00,-0.6612093864662647D+00,  &
! 0.4679139345726904D+00,-0.2386191860831970D+00/
!     Gauss Legendre Points and Weights, N = 12
! DATA ( w(i,2), x(i,2), i = 1,6) /  &
! 0.4717533638651177D-01,-0.9815606342467191D+00,  &
! 0.1069393259953183D+00,-0.9041172563704750D+00,  &
! 0.1600783285433464D+00,-0.7699026741943050D+00,  &
! 0.2031674267230659D+00,-0.5873179542866171D+00,  &
! 0.2334925365383547D+00,-0.3678314989981802D+00,  &
! 0.2491470458134029D+00,-0.1252334085114692D+00/
!     Gauss Legendre Points and Weights, N = 20
! DATA ( w(i,3), x(i,3), i = 1,10) /  &
! 0.1761400713915212D-01,-0.9931285991850949D+00,  &
! 0.4060142980038694D-01,-0.9639719272779138D+00,  &
! 0.6267204833410906D-01,-0.9122344282513259D+00,  &
! 0.8327674157670475D-01,-0.8391169718222188D+00,  &
! 0.1019301198172404D+00,-0.7463319064601508D+00,  &
! 0.1181945319615184D+00,-0.6360536807265150D+00,  &
! 0.1316886384491766D+00,-0.5108670019508271D+00,  &
! 0.1420961093183821D+00,-0.3737060887154196D+00,  &
! 0.1491729864726037D+00,-0.2277858511416451D+00,  &
! 0.1527533871307259D+00,-0.7652652113349733D-01/
REAL (dp), PARAMETER :: w(10,3) = RESHAPE( (/      &
      0.1713244923791705D+00, 0.3607615730481384D+00, 0.4679139345726904D+00, &
        0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                  &
      0.4717533638651177D-01, 0.1069393259953183D+00, 0.1600783285433464D+00, &
      0.2031674267230659D+00, 0.2334925365383547D+00, 0.2491470458134029D+00, &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      0.1761400713915212D-01, 0.4060142980038694D-01, 0.6267204833410906D-01, &
      0.8327674157670475D-01, 0.1019301198172404D+00, 0.1181945319615184D+00, &
      0.1316886384491766D+00, 0.1420961093183821D+00, 0.1491729864726037D+00, &
      0.1527533871307259D+00 /), (/ 10, 3 /) )
REAL (dp), PARAMETER :: x(10,3) = RESHAPE( (/      &
      -0.9324695142031522D+00, -0.6612093864662647D+00,   &
      -0.2386191860831970D+00,  0.0D0, 0.0D0, 0.0D0,      &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9815606342467191D+00, -0.9041172563704750D+00,   &
      -0.7699026741943050D+00, -0.5873179542866171D+00,   &
      -0.3678314989981802D+00, -0.1252334085114692D+00,   &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9931285991850949D+00, -0.9639719272779138D+00,   &
      -0.9122344282513259D+00, -0.8391169718222188D+00,   &
      -0.7463319064601508D+00, -0.6360536807265150D+00,   &
      -0.5108670019508271D+00, -0.3737060887154196D+00,   &
      -0.2277858511416451D+00, -0.7652652113349733D-01 /), (/ 10, 3 /) )

IF ( ABS(r) < 0.3 ) THEN
  ng = 1
  lg = 3
ELSE IF ( ABS(r) < 0.75 ) THEN
  ng = 2
  lg = 6
ELSE
  ng = 3
  lg = 10
END IF
h = sh
k = sk
hk = h*k
bvn = zero
IF ( ABS(r) < 0.925 ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN(r)
  DO  i = 1, lg
    sn = SIN(asr*( x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
    sn = SIN(asr*(-x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
  END DO
  bvn = bvn*asr/(2*twopi) + phi(-h)*phi(-k)
ELSE
  IF ( r < zero ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < one ) THEN
    as = ( one - r )*( one + r )
    a = SQRT(as)
    bs = ( h - k )**2
    c = ( 4. - hk )/8
    d = ( 12. - hk )/16.
    bvn = a*EXP( -(bs/as + hk)/2. )  &
    *( one - c*(bs - as)*(one - d*bs/5.)/3. + c*d*as*as/5. )
    IF ( hk > -160. ) THEN
      b = SQRT(bs)
      bvn = bvn - EXP(-hk/2)*SQRT(twopi)*phi(-b/a)*b  &
                      *( one - c*bs*( one - d*bs/5. )/3. )
    END IF
    a = a/2
    DO i = 1, lg
      xs = ( a*(x(i,ng) + one) )**2
      rs = SQRT( one - xs )
      bvn = bvn + a*w(i,ng)*( EXP( -bs/(2*xs) - hk/(1+rs) )/rs  &
                - EXP( -(bs/xs+hk)/2. )*( one + c*xs*( one + d*xs ) ) )
      xs = as*(-x(i,ng) + one)**2/4.
      rs = SQRT( 1 - xs )
      bvn = bvn + a*w(i,ng)*EXP( -(bs/xs + hk)/2 )  &
                * ( EXP( -hk*(one - rs)/(2*(one + rs)) )/rs - &
                       ( one + c*xs*( one + d*xs ) ) )
    END DO
    bvn = -bvn/twopi
  END IF
  IF ( r > 0 ) bvn =  bvn + phi( -MAX( h, k ) )
  IF ( r < 0 ) bvn = -bvn + MAX( zero, phi(-h) - phi(-k) )
END IF
fn_val = bvn

RETURN
END FUNCTION bvnu


FUNCTION uni() RESULT(fn_val)

!     Random number generator, adapted from F. James
!     "A Review of Random Number Generators"
!      Comp. Phys. Comm. 60(1990), pp. 329-344.

REAL :: fn_val

! Local variables
REAL, SAVE      :: seeds(24) = (/   &
        0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,  &
        0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043,  &
        0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127,  &
        0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357 /)
REAL, PARAMETER :: twom24 = 1.0/16777216
INTEGER, SAVE   :: i = 24, j = 10
REAL, SAVE      :: carry = 0.0

fn_val = seeds(i) - seeds(j) - carry
IF ( fn_val < zero ) THEN
  fn_val = fn_val + one
  carry = twom24
ELSE
  carry = zero
END IF

seeds(i) = fn_val
i = 24 - MOD( 25-i, 24 )
j = 24 - MOD( 25-j, 24 )

RETURN
END FUNCTION uni


SUBROUTINE adapt(ndim, mincls, maxcls, functn,  &
                 absreq, relreq, absest, finest, inform)

!   Adaptive Multidimensional Integration Subroutine

!   Author: Alan Genz
!           Department of Mathematics
!           Washington State University
!           Pullman, WA 99164-3113 USA

!  This subroutine computes an approximation to the integral

!      1 1     1
!     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
!      0 0     0

!**************  Parameters for ADAPT  ********************************

!***** Input Parameters

!  NDIM    Integer number of integration variables.
!  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS
!          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes
!          that a previous call of ADAPT has been made with the same
!          integrand and continues that calculation.
!  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS
!          must be >= RULCLS, the number of function calls required for
!          one application of the basic integration rule.
!           IF ( NDIM .EQ. 1 ) THEN
!              RULCLS = 11
!           ELSE IF ( NDIM .LT. 15 ) THEN
!              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1
!           ELSE
!              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3
!           ENDIF
!  FUNCTN  Externally declared real user defined integrand. Its
!          parameters must be (NDIM, Z), where Z is a real array of
!          length NDIM.
!  ABSREQ  Real required absolute accuracy.
!  RELREQ  Real required relative accuracy.

!***** Output Parameters

!  MINCLS  Actual number of FUNCTN calls used by ADAPT.
!  ABSEST  Real estimated absolute accuracy.
!  FINEST  Real estimated value of integral.
!  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
!                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
!          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
!                     result FINEST to within the requested accuracy.
!          INFORM = 2 if MINCLS > MAXCLS or RULCLS > MAXCLS.

!***********************************************************************

!     Begin driver routine. This routine partitions the working storage
!      array and then calls the main subroutine ADBASE.

INTEGER, INTENT(IN)     :: ndim, maxcls
INTEGER, INTENT(IN OUT) :: mincls
INTEGER, INTENT(OUT)    :: inform
REAL (dp), INTENT(IN)   :: absreq, relreq
REAL (dp), INTENT(OUT)  :: absest, finest

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

! Local variables
INTEGER       :: mxrgns, rulcls, lenrul
INTEGER, SAVE :: sbrgns

IF ( ndim == 1 ) THEN
  lenrul = 5
  rulcls = 9
ELSE IF ( ndim < 12 ) THEN
  lenrul = 6
  rulcls = 2**ndim + 2*ndim*(ndim+2) + 1
ELSE
  lenrul = 6
  rulcls = 1 + 2*ndim*(1+2*ndim)
END IF

IF ( rulcls <= maxcls .AND. mincls <= maxcls ) THEN
  mxrgns = 500
  CALL adbase(ndim, mincls, maxcls, functn, absreq, relreq,  &
              absest, finest, sbrgns, mxrgns, rulcls, lenrul, inform)
ELSE
  inform = 2
  mincls = rulcls
END IF

RETURN
END SUBROUTINE adapt


SUBROUTINE bsinit(ndim, w, lenrul, g)

!     For initializing basic rule weights and symmetric sum parameters.

INTEGER, INTENT(IN)    :: ndim, lenrul
! REAL (dp) :: w(lenrul,4), g(ndim,lenrul)
REAL (dp), INTENT(OUT) :: w(:,:), g(:,:)

! Local variables
INTEGER, PARAMETER :: numnul = 4, sdim = 12
INTEGER            :: rulpts(6), i, j
REAL (dp)          :: lam1, lam2, lam3, lamp, rulcon

!     The following code determines rule parameters and weights for a
!      degree 7 rule (W(1,1),...,W(5,1)), two degree 5 comparison rules
!      (W(1,2),...,W(5,2) and W(1,3),...,W(5,3)) and a degree 3
!      comparison rule (W(1,4),...W(5,4)).

!       If NDIM = 1, then LENRUL = 5 and total points = 9.
!       If NDIM < SDIM, then LENRUL = 6 and
!                      total points = 1+2*NDIM*(NDIM+2)+2**NDIM.
!       If NDIM > = SDIM, then LENRUL = 6 and
!                      total points = 1+2*NDIM*(1+2*NDIM).

DO i = 1,lenrul
  g(:ndim,i) = zero
  w(i,:ndim) = zero
END DO
rulpts(5) = 2*ndim*(ndim-1)
rulpts(4) = 2*ndim
rulpts(3) = 2*ndim
rulpts(2) = 2*ndim
rulpts(1) = 1
lamp = 0.85
lam3 = 0.4707
lam2 = 4/(15 - 5/lam3)
w(5,1) = ( 3 - 5*lam3 )/( 180*(lam2-lam3)*lam2**2 )
IF ( ndim < sdim ) THEN
  lam1 = 8*lam3*(31*lam3-15)/( (3*lam3-1)*(5*lam3-3)*35 )
  w(lenrul,1) = 1/(3*lam3)**3/2**ndim
ELSE
  lam1 = ( lam3*(15 - 21*lam2) + 35*(ndim-1)*(lam2-lam3)/9 )  &
         / ( lam3*(21 - 35*lam2) + 35*(ndim-1)*(lam2/lam3-1)/9 )
  w(6,1) = 1/(4*(3*lam3)**3)
END IF
w(3,1) = ( 15 - 21*(lam3+lam1) + 35*lam3*lam1 )  &
         / ( 210*lam2*(lam2-lam3)*(lam2-lam1) ) - 2*(ndim-1)*w(5,1)
w(2,1) = ( 15 - 21*(lam3+lam2) + 35*lam3*lam2 )  &
         / ( 210*lam1*(lam1-lam3)*(lam1-lam2) )
IF ( ndim < sdim ) THEN
  rulpts(lenrul) = 2**ndim
  lam3 = SQRT(lam3)
  g(1:ndim,lenrul) = lam3
ELSE
  w(6,1) = 1/(4*(3*lam3)**3)
  rulpts(6) = 2*ndim*(ndim-1)
  lam3 = SQRT(lam3)
  g(1:2,6) = lam3
END IF
IF ( ndim > 1 ) THEN
  w(5,2) = 1/(6*lam2)**2
  w(5,3) = 1/(6*lam2)**2
END IF
w(3,2) = ( 3 - 5*lam1 )/( 30*lam2*(lam2-lam1) )- 2*(ndim-1)*w(5,2)
w(2,2) = ( 3 - 5*lam2 )/( 30*lam1*(lam1-lam2) )
w(4,3) = ( 3 - 5*lam2 )/( 30*lamp*(lamp-lam2) )
w(3,3) = ( 3 - 5*lamp )/( 30*lam2*(lam2-lamp) )- 2*(ndim-1)*w(5,3)
w(2,4) = 1/(6*lam1)
lamp = SQRT(lamp)
lam2 = SQRT(lam2)
lam1 = SQRT(lam1)
g(1,2) = lam1
g(1,3) = lam2
g(1,4) = lamp
IF ( ndim > 1 ) THEN
  g(1,5) = lam2
  g(2,5) = lam2
END IF
DO j = 1, numnul
  w(1,j) = one - DOT_PRODUCT( rulpts(2:lenrul), w(2:lenrul,j) )
END DO
rulcon = 2
CALL rulnrm( lenrul, numnul, rulpts, w, rulcon )

RETURN
END SUBROUTINE bsinit


SUBROUTINE rulnrm( lenrul, numnul, rulpts, w, rulcon )

INTEGER, INTENT(IN)       :: lenrul, numnul, rulpts(:)
REAL (dp), INTENT(IN)     :: rulcon
REAL (dp), INTENT(IN OUT) :: w(:,:)

! Local variables
INTEGER   :: i, j, k
REAL (dp) :: alpha, normcf, normnl

!     Compute orthonormalized null rules.

normcf = DOT_PRODUCT( rulpts(1:lenrul), w(1:lenrul,1)**2 )
DO k = 2,numnul
  w(1:lenrul,k) = w(1:lenrul,k) - w(1:lenrul,1)
  DO j = 2,k-1
    alpha = zero
    DO i = 1,lenrul
      alpha = alpha + rulpts(i)*w(i,j)*w(i,k)
    END DO
    alpha = -alpha/normcf
    w(1:lenrul,k) = w(1:lenrul,k) + alpha*w(1:lenrul,j)
  END DO
  normnl = zero
  DO i = 1,lenrul
    normnl = normnl + rulpts(i)*w(i,k)*w(i,k)
  END DO
  alpha = SQRT(normcf/normnl)
  w(1:lenrul,k) = alpha*w(1:lenrul,k)
END DO
DO j = 2, numnul
  w(1:lenrul,j) = w(1:lenrul,j) / rulcon
END DO

RETURN
END SUBROUTINE rulnrm


SUBROUTINE adbase(ndim, mincls, maxcls, functn, absreq, relreq,  &
                  absest, finest, sbrgns, mxrgns, rulcls, lenrul, inform)

!        Main adaptive integration subroutine

INTEGER, INTENT(IN)     :: ndim, maxcls, mxrgns, lenrul, rulcls
INTEGER, INTENT(IN OUT) :: sbrgns, mincls
INTEGER, INTENT(OUT)    :: inform
REAL (dp), INTENT(IN)   :: absreq, relreq
REAL (dp), INTENT(OUT)  :: absest, finest

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

! Local variables
REAL (dp), ALLOCATABLE :: errors(:), values(:), lowers(:,:), uppers(:,:),   &
                          meshes(:,:), weghts(:,:), points(:,:), lower(:),  &
                          upper(:), width(:), mesh(:)

! Local variables
INTEGER   :: i, j, pontrs(mxrgns), nwrgns, divaxn, top, rgncls, funcls,  &
             difcls

! Allocation of working arrays
ALLOCATE( errors(mxrgns), values(mxrgns), lowers(ndim,mxrgns),  &
          uppers(ndim,mxrgns), meshes(ndim,mxrgns), weghts(lenrul,4),  &
          points(ndim,lenrul), lower(ndim), upper(ndim), width(ndim),  &
          mesh(ndim) )

!     Initialization of subroutine

inform = 2
funcls = 0
CALL bsinit(ndim, weghts, lenrul, points)
IF ( mincls >= 0) THEN
  
!       When MINCLS >= 0 determine initial subdivision of the
!       integration region and apply basic rule to each subregion.
  
  sbrgns = 0
  DO i = 1,ndim
    lower(i) = zero
    mesh(i) = one
    width(i) = one/(2._dp*mesh(i))
    upper(i) = one
  END DO
  divaxn = 0
  rgncls = rulcls
  nwrgns = 1
  10 CALL differ(ndim, lower, upper, width, functn, divaxn, difcls)
  IF ( rgncls*(mesh(divaxn)+1)/mesh(divaxn) <= mincls ) THEN
    rgncls = rgncls*(mesh(divaxn)+1)/mesh(divaxn)
    nwrgns = nwrgns*(mesh(divaxn)+1)/mesh(divaxn)
    width(divaxn) = width(divaxn)*mesh(divaxn)/(mesh(divaxn)+1)
    mesh(divaxn) = mesh(divaxn) + 1
    GO TO 10
  END IF
  IF ( nwrgns <= mxrgns ) THEN
    DO i = 1,ndim
      upper(i) = lower(i) + 2*width(i)
      mesh(i) = one
    END DO
  END IF
  
!     Apply basic rule to subregions and store results in heap.
  
  20 sbrgns = sbrgns + 1
  CALL basrul(ndim, lower, upper, width, functn, weghts, lenrul, points,  &
              errors(sbrgns), values(sbrgns))
  CALL trestr(sbrgns, sbrgns, pontrs, errors)
  DO i = 1,ndim
    lowers(i,sbrgns) = lower(i)
    uppers(i,sbrgns) = upper(i)
    meshes(i,sbrgns) = mesh(i)
  END DO
  DO i = 1,ndim
    lower(i) = upper(i)
    upper(i) = lower(i) + 2*width(i)
    IF ( lower(i)+width(i) < 1 )  GO TO 20
    lower(i) = zero
    upper(i) = lower(i) + 2*width(i)
  END DO
  funcls = sbrgns*(difcls+rulcls)
END IF

!     Check for termination

30 finest = zero
absest = zero
DO i = 1, sbrgns
  finest = finest + values(i)
  absest = absest + errors(i)
END DO
IF ( absest > MAX( absreq, relreq*ABS(finest) )  &
     .OR. funcls < mincls ) THEN
  
!     Prepare to apply basic rule in (parts of) subregion with
!     largest error.
  
  top = pontrs(1)
  rgncls = rulcls
  DO i = 1,ndim
    lower(i) = lowers(i,top)
    upper(i) = uppers(i,top)
    mesh(i) = meshes(i,top)
    width(i) = (upper(i)-lower(i))/(2*mesh(i))
    rgncls = rgncls*mesh(i)
  END DO
  CALL differ(ndim, lower, upper, width, functn, divaxn, difcls)
  funcls = funcls + difcls
  rgncls = rgncls*(mesh(divaxn)+1)/mesh(divaxn)
  IF ( funcls + rgncls <= maxcls ) THEN
    IF ( sbrgns + 1 <= mxrgns ) THEN
      
!     Prepare to subdivide into two pieces.
      
      nwrgns = 1
      width(divaxn) = width(divaxn)/2
    ELSE
      nwrgns = 0
      width(divaxn) = width(divaxn)*mesh(divaxn)/( mesh(divaxn) + 1 )
      meshes(divaxn,top) = mesh(divaxn) + 1
    END IF
    IF ( nwrgns > 0 ) THEN
      
!     Only allow local subdivision when space is available.
      
      DO j = sbrgns+1,sbrgns+nwrgns
        DO i = 1,ndim
          lowers(i,j) = lower(i)
          uppers(i,j) = upper(i)
          meshes(i,j) = mesh(i)
        END DO
      END DO
      uppers(divaxn,top) = lower(divaxn) + 2*width(divaxn)
      lowers(divaxn,sbrgns+1) = uppers(divaxn,top)
    END IF
    funcls = funcls + rgncls
    CALL basrul(ndim, lowers(1:,top), uppers(1:,top), width,  &
                functn, weghts, lenrul, points, errors(top), values(top))
    CALL trestr(top, sbrgns, pontrs, errors)
    DO i = sbrgns+1, sbrgns+nwrgns
      
!     Apply basic rule and store results in heap.
      
      CALL basrul(ndim, lowers(1:,i), uppers(1:,i), width, functn, weghts,  &
                  lenrul, points, errors(i), values(i))
      CALL trestr(i, i, pontrs, errors)
    END DO
    sbrgns = sbrgns + nwrgns
    GO TO 30
  ELSE
    inform = 1
  END IF
ELSE
  inform = 0
END IF
mincls = funcls
DEALLOCATE( errors, values, lowers, uppers, meshes, weghts, points, lower,  &
            upper, width, mesh )

RETURN
END SUBROUTINE adbase


SUBROUTINE basrul( ndim, a, b, width, functn, w, lenrul, g, rgnert, basest )

!     For application of basic integration rule

INTEGER, INTENT(IN)       :: lenrul, ndim
REAL (dp), INTENT(IN)     :: a(:), b(:), width(:), w(:,:)
REAL (dp), INTENT(IN OUT) :: g(:,:)
REAL (dp), INTENT(OUT)    :: rgnert, basest

! Local variables
INTEGER   :: i
REAL (dp) :: fsymsm, rgncmp, rgnval, rgnvol, rgncpt, rgnerr, center(ndim)

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

!     Compute Volume and Center of Subregion

rgnvol = 1
DO i = 1,ndim
  rgnvol = 2*rgnvol*width(i)
  center(i) = a(i) + width(i)
END DO
basest = zero
rgnert = zero

!     Compute basic rule and error

10 rgnval = zero
rgnerr = zero
rgncmp = zero
rgncpt = zero
DO i = 1,lenrul
  CALL fulsum(ndim, center, width, g(1:,i), functn, fsymsm)
!     Basic Rule
  rgnval = rgnval + w(i,1)*fsymsm
!     First comparison rule
  rgnerr = rgnerr + w(i,2)*fsymsm
!     Second comparison rule
  rgncmp = rgncmp + w(i,3)*fsymsm
!     Third Comparison rule
  rgncpt = rgncpt + w(i,4)*fsymsm
END DO

!     Error estimation

rgnerr = SQRT(rgncmp**2 + rgnerr**2)
rgncmp = SQRT(rgncpt**2 + rgncmp**2)
IF ( 4*rgnerr < rgncmp ) rgnerr = rgnerr/2
IF ( 2*rgnerr > rgncmp ) rgnerr = MAX( rgnerr, rgncmp )
rgnert = rgnert +  rgnvol*rgnerr
basest = basest +  rgnvol*rgnval

!     When subregion has more than one piece, determine next piece and
!      loop back to apply basic rule.

DO i = 1,ndim
  center(i) = center(i) + 2*width(i)
  IF ( center(i) < b(i) ) GO TO 10
  center(i) = a(i) + width(i)
END DO

RETURN
END SUBROUTINE basrul


SUBROUTINE fulsum(s, center, hwidth, g, f, fn_val)

!***  To compute fully symmetric basic rule sum

INTERFACE
  FUNCTION f(ndim, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER, INTENT(IN)       :: s
REAL (dp), INTENT(IN)     :: center(:), hwidth(:)
REAL (dp), INTENT(IN OUT) :: g(:)
REAL (dp), INTENT(OUT)    :: fn_val

! Local variables
INTEGER   :: ixchng, lxchng, i, l
REAL (dp) :: intsum, gl, gi, x(s)

fn_val = zero

!     Compute centrally symmetric sum for permutation of G

10 intsum = zero
x(1:s) = center(1:s) + g(1:s)*hwidth(1:s)
20 intsum = intsum + f(s,x)
DO i = 1,s
  g(i) = -g(i)
  x(i) = center(i) + g(i)*hwidth(i)
  IF ( g(i) < zero ) GO TO 20
END DO
fn_val = fn_val + intsum

!     Find next distinct permutation of G and loop back for next sum

DO i = 2,s
  IF ( g(i-1) > g(i) ) THEN
    gi = g(i)
    ixchng = i - 1
    DO l = 1,(i-1)/2
      gl = g(l)
      g(l) = g(i-l)
      g(i-l) = gl
      IF ( gl <= gi ) ixchng = ixchng - 1
      IF ( g(l) > gi ) lxchng = l
    END DO
    IF ( g(ixchng) <= gi ) ixchng = lxchng
    g(i) = g(ixchng)
    g(ixchng) = gi
    GO TO 10
  END IF
END DO

!     End loop for permutations of G and associated sums

!     Restore original order to G's

DO i = 1,s/2
  gi = g(i)
  g(i) = g(s+1-i)
  g(s+1-i) = gi
END DO

RETURN
END SUBROUTINE fulsum


SUBROUTINE differ(ndim, a, b, width, functn, divaxn, difcls)

!     Compute fourth differences and subdivision axes

INTEGER, INTENT(IN)     :: ndim
INTEGER, INTENT(IN OUT) :: divaxn
INTEGER, INTENT(OUT)    :: difcls
REAL (dp), INTENT(IN)   :: a(:), b(:), width(:)

! Local variables
INTEGER   :: i
REAL (dp) :: frthdf, funcen, widthi, z(ndim), dif(ndim)

INTERFACE
  FUNCTION functn(ndim, x) RESULT(fn_val)
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)   :: ndim
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val
  END FUNCTION functn
END INTERFACE

difcls = 0
divaxn = MOD( divaxn, ndim ) + 1
IF ( ndim > 1 ) THEN
  DO i = 1,ndim
    dif(i) = zero
    z(i) = a(i) + width(i)
  END DO
  10 funcen = functn(ndim, z)
  DO i = 1,ndim
    widthi = width(i)/5
    frthdf = 6*funcen
    z(i) = z(i) - 4*widthi
    frthdf = frthdf + functn(ndim,z)
    z(i) = z(i) + 2*widthi
    frthdf = frthdf - 4*functn(ndim,z)
    z(i) = z(i) + 4*widthi
    frthdf = frthdf - 4*functn(ndim,z)
    z(i) = z(i) + 2*widthi
    frthdf = frthdf + functn(ndim,z)
!     Do not include differences below roundoff
    IF ( funcen + frthdf/8 /= funcen )  &
    dif(i) = dif(i) + ABS(frthdf)*width(i)
    z(i) = z(i) - 4*widthi
  END DO
  difcls = difcls + 4*ndim + 1
  DO i = 1,ndim
    z(i) = z(i) + 2*width(i)
    IF ( z(i) < b(i) ) GO TO 10
    z(i) = a(i) + width(i)
  END DO
  DO i = 1,ndim
    IF ( dif(divaxn) < dif(i) ) divaxn = i
  END DO
END IF

RETURN
END SUBROUTINE differ


SUBROUTINE trestr(pointr, sbrgns, pontrs, rgners)
!***BEGIN PROLOGUE TRESTR
!***PURPOSE TRESTR maintains a heap for subregions.
!***DESCRIPTION TRESTR maintains a heap for subregions.
!            The subregions are ordered according to the size of the
!            greatest error estimates of each subregion (RGNERS).

!   PARAMETERS

!     POINTR Integer.
!            The index for the subregion to be inserted in the heap.
!     SBRGNS Integer.
!            Number of subregions in the heap.
!     PONTRS Integer array of dimension SBRGNS.
!            Used to store the indices for the greatest estimated errors
!            for each subregion.
!     RGNERS Real array of dimension SBRGNS.
!            Used to store the greatest estimated errors for each subregion.

!***ROUTINES CALLED NONE
!***END PROLOGUE TRESTR

!   Global variables.

INTEGER, INTENT(IN)       :: pointr, sbrgns
REAL (dp), INTENT(IN OUT) :: rgners(:)
INTEGER, INTENT(IN OUT)   :: pontrs(:)

!   Local variables.

!   RGNERR Intermediate storage for the greatest error of a subregion.
!   SUBRGN Position of child/parent subregion in the heap.
!   SUBTMP Position of parent/child subregion in the heap.

INTEGER   :: subrgn, subtmp
REAL (dp) :: rgnerr

!***FIRST PROCESSING STATEMENT TRESTR

rgnerr = rgners(pointr)
IF ( pointr == pontrs(1)) THEN
  
!        Move the new subregion inserted at the top of the heap
!        to its correct position in the heap.
  
  subrgn = 1
  10 subtmp = 2*subrgn
  IF ( subtmp <= sbrgns ) THEN
    IF ( subtmp /= sbrgns ) THEN
      
!              Find maximum of left and right child.
      
      IF ( rgners(pontrs(subtmp)) <  &
           rgners(pontrs(subtmp+1)) ) subtmp = subtmp + 1
    END IF
    
!           Compare maximum child with parent.
!           If parent is maximum, then done.
    
    IF ( rgnerr < rgners(pontrs(subtmp)) ) THEN
      
!              Move the pointer at position subtmp up the heap.
      
      pontrs(subrgn) = pontrs(subtmp)
      subrgn = subtmp
      GO TO 10
    END IF
  END IF
ELSE
  
!        Insert new subregion in the heap.
  
  subrgn = sbrgns
  20 subtmp = subrgn/2
  IF ( subtmp >= 1 ) THEN
    
!           Compare child with parent. If parent is maximum, then done.
    
    IF ( rgnerr > rgners(pontrs(subtmp)) ) THEN
      
!              Move the pointer at position subtmp down the heap.
      
      pontrs(subrgn) = pontrs(subtmp)
      subrgn = subtmp
      GO TO 20
    END IF
  END IF
END IF
pontrs(subrgn) = pointr

!***END TRESTR
RETURN
END SUBROUTINE trestr



SUBROUTINE adjust_infi(a, infi)

! Adjust INFI when a limit is too far out in the tail

REAL (dp), INTENT(IN)    :: a
INTEGER, INTENT(IN OUT)  :: infi

IF (a < - xinfin) THEN
  IF (infi == 1) THEN
    infi = -1
  ELSE IF (infi > 1) THEN
    infi = 0
  END IF
ELSE IF (a > xinfin) THEN
  IF (infi == 0) THEN
    infi = -1
  ELSE IF (infi > 1) THEN
    infi = 1
  END IF
END IF

RETURN
END SUBROUTINE adjust_infi

END MODULE mvnpack
