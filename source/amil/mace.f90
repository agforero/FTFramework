MODULE Ace
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

! BLOCK DATA
! COMMON /parms/ itape, maxit, nterm, span, alpha, big

!     block data
!     common /parms/ itape,maxit,nterm,span,alpha,big

!------------------------------------------------------------------

! these procedure parameters can be changed in the calling routine
! by defining the above labeled common and resetting the values with
! executable statements.

! itape : fortran file number for printer output.
!         (itape.le.0 => no printer output.)
! maxit : maximum number of iterations.
! nterm : number of consecutive iterations for which
!         rsq must change less than delcor for convergence.
! span, alpha : super smoother parameters (see below).
! big : a large representable floating point number.

!------------------------------------------------------------------

! DATA itape, maxit, nterm, span, alpha, big /6, 20, 3, 0.0, 0.0, 1.0E20/
! END

! COMMON /parms/ itape, maxit, nterm, span, alpha, big
INTEGER, SAVE  :: itape = 6, maxit = 20, nterm = 3
REAL, SAVE     :: span = 0.0, alpha = 0.0, big = 1.E20


! BLOCK DATA
! COMMON /spans/ spans(3) /consts/ big, sml, eps

!---------------------------------------------------------------

! this sets the compile time (default) values for various
! internal parameters :

! spans : span values for the three running linear smoothers.
! spans(1) : tweeter span.
! spans(2) : midrange span.
! spans(3) : woofer span.
! (these span values should be changed only with care.)
! big : a large representable floating point number.
! sml : a small number. should be set so that (sml)**(10.0) does
!       not cause floating point underflow.
! eps : used to numerically stabilize slope calculations for
!       running linear fits.

! these parameter values can be changed by declaring the
! relevant labeled common in the main program and resetting
! them with executable statements.

!-----------------------------------------------------------------

! DATA spans, big, sml, eps /0.05, 0.2, 0.5, 1.0E20, 1.0E-7, 1.0E-3/
! END

REAL, SAVE  :: spans(3) = 0.05, sml = 1.0E-7, eps = 1.0E-3


CONTAINS


SUBROUTINE mace(p, n, x, y, w, l, delrsq, ns, tx, ty, rsq, ierr, m, z)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-10  Time: 18:07:25

!------------------------------------------------------------------

! estimate multiple optimal transformations for regression and
! correlation by alternating conditional expectation estimates.

! version 3/28/85.

! Breiman and Friedman `Estimating optimal transformations for multiple
! regression and correlation',
! Journal of the American Statistical Association, Vol.80, pp.580-614 (1985).

! coded  and copywrite (c) 1985 by:

!                        jerome h. friedman
!                     department of statistics
!                               and
!                stanford linear accelerator center
!                        stanford university

! all rights reserved.

! N.B. Jerry Friedman has given permission for this version to be made
!      available on this web site.

! input:

!    n : number of observations.
!    p : number of predictor variables for each observation.
!    x(p,n) : predictor data matrix.
!    y(n) : response values for the observations.
!       missing values are signified by a value (response or
!       predictor) greater than or equal to big.
!       (see below - default, big = 1.0e20)
!    w(n) : weights for the observations.
!    l(p+1) : flag for each variable.
!       l(1) through l(p) : predictor variables.
!       l(p+1) : response variable.
!       l(i)=0 => ith variable not to be used.
!       l(i)=1 => ith variable assumes orderable values.
!       l(i)=2 => ith variable assumes circular (periodic) values
!                 in the range (0.0,1.0) with period 1.0.
!       l(i)=3 => ith variable transformation is to be monotone.
!       l(i)=4 => ith variable transformation is to be linear.
!       l(i)=5 => ith variable assumes categorical (unorderable) values.
!   delrsq : termination threshold. iteration stops when
!       rsq changes less than delrsq in nterm
!       consecutive iterations (see below - default, nterm=3).
!   ns : number of eigensolutions (sets of transformations).

! output:

!   tx(n,p,ns) : predictor transformations.
!      tx(j,i,k) = transformed value of ith predictor for jth obs
!                  for kth eigensolution.
!   ty(n,ns) = response transformations.
!      ty(j,k) = transformed response value for jth observation
!                for kth eigensolution.
!   rsq(ns) = fraction of variance(ty<y>)
!                       p
!         explained by sum tx(i)<x(i)>  for each eigensolution.
!                      i=1
!   ierr : error flag.
!      ierr = 0 : no errors detected.
!      ierr > 0 : error detected - see format statements below.

! scratch:

!    m(n,p+1), z(n,12) : internal working storage.

! note: mace uses an iterative procedure for solving the optimization
!    problem. default starting transformations are ty(j,k)=y(j),
!    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns.  Other starting transformations
!    can be specified (if desired) for either the response and/or any of
!    the predictor variables.  This is signaled by negating the
!    corresponding l(i) value and storing the starting transformed values
!    in the corresponding array (ty(j,k), tx(j,i,k)) before calling mace.

!------------------------------------------------------------------

INTEGER, INTENT(IN)      :: p
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: x(p,n)
REAL, INTENT(IN)         :: y(n)
REAL, INTENT(IN)         :: w(n)
INTEGER, INTENT(IN OUT)  :: l(1)
REAL, INTENT(IN OUT)     :: delrsq
INTEGER, INTENT(IN)      :: ns
REAL, INTENT(OUT)        :: tx(n,p,ns)
REAL, INTENT(OUT)        :: ty(n,ns)
REAL, INTENT(OUT)        :: rsq(ns)
INTEGER, INTENT(OUT)     :: ierr
INTEGER, INTENT(OUT)     :: m(n,1)
REAL, INTENT(OUT)        :: z(n,12)

INTEGER   :: i, is, ism1, iter, j, js, k, np, nit, nt, pp1
REAL      :: cmn, cmx, ct(10), rsqi
REAL (dp) :: sm, sv, sw, sw1

ierr = 0
pp1 = p + 1
sm = 0.0
sv = sm
sw = sv
sw1 = sw
DO  i = 1, pp1
  IF (l(i) < -5.OR.l(i) > 5) THEN
    ierr = 6
    IF (itape > 0) WRITE (itape,5800) i, l(i)
  END IF
END DO
IF (ierr /= 0) RETURN
IF (l(pp1) == 0) THEN
  ierr = 4
  IF (itape > 0) WRITE (itape,5600) pp1
  RETURN
END IF
np = 0
DO  i = 1, p
  IF (l(i) /= 0) np = np + 1
END DO
IF (np <= 0) THEN
  ierr = 5
  IF (itape > 0) WRITE (itape,5700) p
  RETURN
END IF
DO  j = 1, n
  sw = sw + w(j)
END DO
IF (sw <= 0.0) THEN
  ierr = 1
  IF (itape > 0) WRITE (itape,5300)
  RETURN
END IF
DO  is = 1, ns
  IF (itape > 0) WRITE (itape,5000) is
  DO  j = 1, n
    IF (l(pp1) > 0) ty(j,is) = y(j)
  END DO
  DO  i = 1, p
    IF (l(i) == 0) THEN
      DO  j = 1, n
        tx(j,i,is) = 0.0
      END DO
    ELSE
      IF (l(i) > 0) THEN
        DO  j = 1, n
          tx(j,i,is) = x(i,j)
        END DO
      END IF
      DO  j = 1, n
        IF (tx(j,i,is) < big) THEN
          sm = sm + w(j) * tx(j,i,is)
          sw1 = sw1 + w(j)
        END IF
      END DO
      IF (sw1 <= 0.0) THEN
        DO  j = 1, n
          tx(j,i,is) = 0.0
        END DO
        sm = 0.0
        sw1 = sm
      ELSE
        sm = sm / sw1
        DO  j = 1, n
          IF (tx(j,i,is) < big) THEN
            tx(j,i,is) = tx(j,i,is) - sm
          ELSE
            tx(j,i,is) = 0.0
          END IF
        END DO
        sm = 0.0
        sw1 = sm
      END IF
    END IF
  END DO
  DO  j = 1, n
    IF (ty(j,is) < big) THEN
      sm = sm + w(j) * ty(j,is)
      sw1 = sw1 + w(j)
    END IF
  END DO
  IF (sw1 <= 0.0) THEN
    ierr = 1
    IF (itape > 0) WRITE (itape,5300)
    RETURN
  END IF
  sm = sm / sw1
  DO  j = 1, n
    IF (ty(j,is) < big) THEN
      ty(j,is) = ty(j,is) - sm
    ELSE
      ty(j,is) = 0.0
    END IF
  END DO
  DO  j = 1, n
    sv = sv + w(j) * ty(j,is) ** 2
  END DO
  sv = sv / sw
  IF (sv > 0.0) THEN
    sv = 1.0 / SQRT(sv)
  ELSE
    IF (l(pp1) > 0) THEN
      ierr = 2
      IF (itape > 0) WRITE (itape,5400)
    ELSE
      ierr = 3
      IF (itape > 0) WRITE (itape,5500) is
    END IF
    RETURN
  END IF
  DO  j = 1, n
    ty(j,is) = ty(j,is) * sv
  END DO
  IF (is == 1) THEN
    DO  j = 1, n
      m(j,pp1) = j
      z(j,2) = y(j)
    END DO
    CALL sort(z(1:,2), m(1:,pp1), 1, n)
    DO  i = 1, p
      IF (l(i) /= 0) THEN
        DO  j = 1, n
          m(j,i) = j
          z(j,2) = x(i,j)
        END DO
        CALL sort(z(1:,2), m(1:,i), 1, n)
      END IF
    END DO
  END IF
  CALL scale(p, n, w, sw, ty(1:,is), tx(1:,1:,is), delrsq, p, z(1:,5), z(1:,6))
  rsq(is) = 0.0
  iter = 0
  nterm = MIN(nterm,10)
  nt = 0
  ct(1:nterm) = 100.0

  190   iter = iter + 1
  nit = 0
  200   rsqi = rsq(is)
  nit = nit + 1
  DO  j = 1, n
    z(j,5) = ty(j,is)
    DO  i = 1, p
      IF (l(i) /= 0) z(j,5) = z(j,5) - tx(j,i,is)
    END DO
  END DO
  DO  i = 1, p
    IF (l(i) /= 0) THEN
      DO  j = 1, n
        k = m(j,i)
        z(j,1) = z(k,5) + tx(k,i,is)
        z(j,2) = x(i,k)
        z(j,4) = w(k)
      END DO
      CALL smothr(ABS(l(i)), n, z(1,2), z, z(1,4), z(1,3), z(1,6))
      sm = 0.0
      DO  j = 1, n
        sm = sm + z(j,4) * z(j,3)
      END DO
      sm = sm / sw
      DO  j = 1, n
        z(j,3) = z(j,3) - sm
      END DO
      sv = 0.0
      DO  j = 1, n
        sv = sv + z(j,4) * (z(j,1)-z(j,3)) ** 2
      END DO
      sv = 1.0 - sv / sw
      IF (sv > rsq(is)) THEN
        rsq(is) = sv
        DO  j = 1, n
          k = m(j,i)
          tx(k,i,is) = z(j,3)
          z(k,5) = z(j,1) - z(j,3)
        END DO
      END IF
    END IF
  END DO
  IF (np /= 1 .AND. rsq(is)-rsqi > delrsq .AND. nit < maxit) THEN
    GO TO 200
  END IF
  DO  j = 1, n
    k = m(j,pp1)
    z(j,2) = y(k)
    z(j,4) = w(k)
    z(j,1) = 0.0
    DO  i = 1, p
      IF (l(i) /= 0) z(j,1) = z(j,1) + tx(k,i,is)
    END DO
  END DO
  CALL smothr(ABS(l(pp1)), n, z(1,2), z, z(1,4), z(1,3), z(1,6))
  IF (is > 1) THEN
    ism1 = is - 1
    DO  js = 1, ism1
      sm = 0.0
      DO  j = 1, n
        k = m(j,pp1)
        sm = sm + w(k) * z(j,3) * ty(k,js)
      END DO
      sm = sm / sw
      DO  j = 1, n
        k = m(j,pp1)
        z(j,3) = z(j,3) - sm * ty(k,js)
      END DO
    END DO
  END IF
  sm = 0.0
  sv = sm
  DO  j = 1, n
    k = m(j,pp1)
    sm = sm + w(k) * z(j,3)
    z(k,2) = z(j,1)
  END DO
  sm = sm / sw
  DO  j = 1, n
    z(j,3) = z(j,3) - sm
    sv = sv + z(j,4) * z(j,3) ** 2
  END DO
  sv = sv / sw
  IF (sv > 0.0) THEN
    sv = 1.0 / SQRT(sv)
  ELSE
    ierr = 3
    IF (itape > 0) WRITE (itape,5500) is
    RETURN
  END IF
  DO  j = 1, n
    k = m(j,pp1)
    ty(k,is) = z(j,3) * sv
  END DO
  sv = 0.0
  DO  j = 1, n
    sv = sv + w(j) * (ty(j,is)-z(j,2)) ** 2
  END DO
  rsq(is) = 1.0 - sv / sw
  IF (itape > 0) WRITE (itape,5200) iter, rsq(is)
  nt = MOD(nt,nterm) + 1
  ct(nt) = rsq(is)
  cmn = 100.0
  cmx = -100.0
  DO  i = 1, nterm
    cmn = MIN(cmn,ct(i))
    cmx = MAX(cmx,ct(i))
  END DO
  IF (cmx-cmn > delrsq .AND. iter < maxit) THEN
    GO TO 190
  END IF
  IF (itape > 0) WRITE (itape,5100) is, rsq(is)
END DO
RETURN
5000 FORMAT ('0eigensolution ',i2,':')
5100 FORMAT (' eigensolution ',i2,'   r**2  =  1 - e**2  =',g12.4)
5200 FORMAT ('     iteration ',i2,'   r**2  =  1 - e**2  =',g12.4)
5300 FORMAT (' ierr=1: sum of weights (w) not positive.')
5400 FORMAT (' ierr=2: y has zero variance.')
5500 FORMAT (' ierr=3: ty(.,',i2,') has zero variance.')
5600 FORMAT (' ierr=4: l(',i2,') must be nonzero.')
5700 FORMAT (' ierr=5: at least one l(1)-l(',i2,') must be nonzero.')
5800 FORMAT (' ierr=6: l(',i2,') =',g12.4, ' must be in the range (-5, 5).')
END SUBROUTINE mace


SUBROUTINE model(p, n, y, w, l, tx, ty, f, t, m, z)

!--------------------------------------------------------------------

! computes response predictive  function f for the model yhat = f(t),
! where
!                                        p
!            f(t) = e(y : t),     t =   sum  tx<i> ( x<i> )
!                                       i=1
! using the x transformations tx constructed by subroutine ace.
! if y is a categorical variable (classification) then
!                                -1
!                       f(t) = ty  (t).
! input:

!    p,n,y,w,l : same input as for subroutine ace.
!    tx,ty,m,z : output from subroutine ace.

! output:

!    f(n),t(n) : input for subroutine acemod.

! note: this subroutine must be called before subroutine acemod.

!-------------------------------------------------------------------

INTEGER, INTENT(IN)      :: p
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: y(n)
REAL, INTENT(IN)         :: w(n)
INTEGER, INTENT(IN OUT)  :: l(1)
REAL, INTENT(IN)         :: tx(n,p)
REAL, INTENT(IN)         :: ty(n)
REAL, INTENT(OUT)        :: f(n)
REAL, INTENT(OUT)        :: t(n)
INTEGER, INTENT(OUT)     :: m(n,1)
REAL, INTENT(OUT)        :: z(n,12)

INTEGER  :: i, j, j1, j2, k, pp1
REAL     :: s

pp1 = p + 1
IF (ABS(l(pp1)) == 5) THEN
  DO  j = 1, n
    t(j) = ty(j)
    m(j,pp1) = j
  END DO
ELSE
  DO  j = 1, n
    s = 0.0
    DO  i = 1, p
      s = s + tx(j,i)
    END DO
    t(j) = s
    m(j,pp1) = j
  END DO
END IF
CALL sort(t, m(1,pp1), 1, n)
DO  j = 1, n
  k = m(j,pp1)
  z(j,2) = w(k)
  IF (y(k) < big) THEN
    z(j,1) = y(k)
  ELSE
    j1 = j
    j2 = j1
    40     IF (y(m(j1,pp1)) >= big) THEN
      j1 = j1 - 1
      IF (j1 >= 1) THEN
        GO TO 40
      END IF
    END IF
    50     IF (y(m(j2,pp1)) >= big) THEN
      j2 = j2 + 1
      IF (j2 <= n) THEN
        GO TO 50
      END IF
    END IF
    IF (j1 < 1) THEN
      k = j2
    ELSE
      IF (j2 > n) THEN
        k = j1
      ELSE
        IF (t(j)-t(j1) < t(j2)-t(j)) THEN
          k = j1
        ELSE
          k = j2
        END IF
      END IF
    END IF
    z(j,1) = y(m(k,pp1))
    t(j) = t(k)
  END IF
END DO
IF (ABS(l(pp1)) == 5) THEN
  DO  j = 1, n
    f(j) = z(j,1)
  END DO
ELSE
  CALL smothr(1, n, t, z, z(1,2), f, z(1,6))
END IF
RETURN
END SUBROUTINE model


SUBROUTINE acemod(v, p, n, x, l, tx, f, t, m, yhat)
!--------------------------------------------------------------------

! computes response y estimates from the model

!                yhat =  f ( t( v ) )

! using the x transformations tx constructed by subroutine ace and
! the predictor function (f,t) constructed by subroutine model.

! input:

!       v(p) : vector of predictor values.
!    p,n,x,l : same input as for subroutine ace.
!       tx,m : output from subroutine ace.
!        f,t : output from subroutine model.

! output:

!    yhat : estimated response value for v.

! note: this subroutine must not be called before subroutine model.

!-------------------------------------------------------------------

INTEGER, INTENT(IN)      :: p
REAL, INTENT(IN)         :: v(p)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: x(p,n)
INTEGER, INTENT(IN OUT)  :: l(1)
REAL, INTENT(IN)         :: tx(n,p)
REAL, INTENT(IN)         :: f(n)
REAL, INTENT(IN)         :: t(n)
INTEGER, INTENT(IN)      :: m(n,1)
REAL, INTENT(OUT)        :: yhat

INTEGER  :: i, jh, jl, low, high, place
REAL     :: th, vi, xt

th = 0.0
DO  i = 1, p
  IF (l(i) /= 0) THEN
    vi = v(i)
    IF (vi >= big) THEN
      IF (x(i,m(n,i)) >= big) th = th + tx(m(n,i),i)
    ELSE
      IF (vi <= x(i,m(1,i))) THEN
        place = 1
      ELSE
        IF (vi >= x(i,m(n,i))) THEN
          place = n
        ELSE
          low = 0
          high = n + 1
          10   IF (low+1 < high) THEN
            place = (low+high) / 2
            xt = x(i,m(place,i))
            IF (vi == xt) GO TO 20
            IF (vi < xt) THEN
              high = place
              GO TO 10
            END IF
            low = place
            GO TO 10
          END IF
          IF (ABS(l(i)) == 5) CYCLE
          jl = m(low,i)
          jh = m(high,i)
          IF (x(i,jh) >= big) THEN
            th = th + tx(jl,i)
            CYCLE
          END IF
          th = th + tx(jl,i) + (tx(jh,i)-tx(jl,i)) * (vi-x(i,jl))  &
              / (x(i,jh)-x(i,jl))
          CYCLE
        END IF
      END IF
      20  th = th + tx(m(place,i),i)
    END IF
  END IF
END DO
IF (th <= t(1)) THEN
  yhat = f(1)
  RETURN
END IF
IF (th >= t(n)) THEN
  yhat = f(n)
  RETURN
END IF
low = 0
high = n + 1
40 IF (low+1 < high) THEN
  place = (low+high) / 2
  xt = t(place)
  IF (th == xt) THEN
    yhat = f(place)
    RETURN
  END IF
  IF (th < xt) THEN
    high = place
    GO TO 40
  END IF
  low = place
  GO TO 40
END IF
IF (ABS(l(p+1)) == 5) THEN
  IF (th-t(low) <= t(high)-th) THEN
    yhat = f(low)
    GO TO 50
  END IF
  yhat = f(high)
ELSE
  yhat = f(low) + (f(high)-f(low)) * (th-t(low)) / (t(high)- t(low))
END IF
50 RETURN
END SUBROUTINE acemod



SUBROUTINE smothr(l, nn, x, y, w, smo, scr)

INTEGER, INTENT(IN)  :: l
INTEGER, INTENT(IN)  :: nn
REAL, INTENT(IN)     :: x(nn)
REAL, INTENT(IN)     :: y(nn)
REAL, INTENT(IN)     :: w(nn)
REAL, INTENT(OUT)    :: smo(nn)
REAL, INTENT(OUT)    :: scr(nn,7)

REAL (dp) :: sm, sw, a, b, d
INTEGER   :: i, j, j0, n, np1

n = nn
sm = 0.0
sw = sm
10 IF (x(n) >= big) THEN
  sm = sm + w(n) * y(n)
  sw = sw + w(n)
  n = n - 1
  IF (n >= 1) THEN
    GO TO 10
  END IF
END IF
IF (n < nn) THEN
  np1 = n + 1
  sm = sm / sw
  DO  j = np1, nn
    smo(j) = sm
  END DO
END IF
IF (n < 1) RETURN
IF (l >= 5) THEN
  j = 1
  30   j0 = j
  sm = w(j) * y(j)
  sw = w(j)
  IF (j < n) THEN
    40  IF (x(j+1) <= x(j)) THEN
      j = j + 1
      sm = sm + w(j) * y(j)
      sw = sw + w(j)
      IF (j < n) THEN
        GO TO 40
      END IF
    END IF
  END IF
  sm = sm / sw
  DO  i = j0, j
    smo(i) = sm
  END DO
  j = j + 1
  IF (j > n) GO TO 120
  GO TO 30
END IF
IF (l == 4) THEN
  sm = 0.0
  sw = sm
  b = sw
  d = b
  DO  j = 1, n
    sm = sm + w(j) * x(j) * y(j)
    sw = sw + w(j) * x(j) ** 2
    b = b + w(j) * x(j)
    d = d + w(j)
  END DO
  a = sw - (b**2) / d
  IF (a <= 0.0) THEN
    a = 0.0
  ELSE
    a = sm / a
  END IF
  b = b / d
  DO  j = 1, n
    smo(j) = a * (x(j)-b)
  END DO
ELSE
  CALL supsmu(n, x, y, w, l, span, alpha, smo, scr)
  IF (l == 3) THEN
    DO  j = 1, n
      scr(j,1) = smo(j)
      scr(n-j+1,2) = scr(j,1)
    END DO
    CALL montne(scr, n)
    CALL montne(scr(1,2), n)
    sm = 0.0
    sw = sm
    DO  j = 1, n
      sm = sm + (smo(j)-scr(j,1)) ** 2
      sw = sw + (smo(j)-scr(n-j+1,2)) ** 2
    END DO
    IF (sm < sw) THEN
      DO  j = 1, n
        smo(j) = scr(j,1)
      END DO
    ELSE
      DO  j = 1, n
        smo(j) = scr(n-j+1,2)
      END DO
    END IF
  END IF
END IF
120 RETURN
END SUBROUTINE smothr


SUBROUTINE montne(x, n)

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(n)

INTEGER  :: bb, eb, br, er, bl, el, i
REAL     :: pmn

bb = 0
eb = bb
10 IF (eb < n) THEN
  bb = eb + 1
  eb = bb
  20   IF (eb < n) THEN
    IF (x(bb) == x(eb+1)) THEN
      eb = eb + 1
      GO TO 20
    END IF
  END IF
  30   IF (eb < n) THEN
    IF (x(eb) > x(eb+1)) THEN
      br = eb + 1
      er = br
      40  IF (er < n) THEN
        IF (x(er+1) == x(br)) THEN
          er = er + 1
          GO TO 40
        END IF
      END IF
      pmn = (x(bb)*(eb-bb+1) + x(br)*(er-br+1)) / (er-bb+1)
      eb = er
      DO  i = bb, eb
        x(i) = pmn
      END DO
    END IF
  END IF
  IF (bb <= 1) GO TO 10
  IF (x(bb-1) <= x(bb)) GO TO 10
  bl = bb - 1
  el = bl
  60   IF (bl > 1) THEN
    IF (x(bl-1) == x(el)) THEN
      bl = bl - 1
      GO TO 60
    END IF
  END IF
  pmn = (x(bb)*(eb-bb+1)+x(bl)*(el-bl+1)) / (eb-bl+1)
  bb = bl
  DO  i = bb, eb
    x(i) = pmn
  END DO
  GO TO 30
END IF
RETURN
END SUBROUTINE montne


SUBROUTINE scale(p, n, w, sw, ty, tx, eps, maxit, r, sc)

INTEGER, INTENT(IN)    :: p
INTEGER, INTENT(IN)    :: n
REAL, INTENT(IN)       :: w(n)
REAL (dp), INTENT(IN)  :: sw
REAL, INTENT(IN OUT)   :: ty(n)
REAL, INTENT(IN OUT)   :: tx(n,p)
REAL, INTENT(IN OUT)   :: eps
INTEGER, INTENT(IN)    :: maxit
REAL, INTENT(OUT)      :: r(n)
REAL, INTENT(OUT)      :: sc(p,5)

REAL (dp) :: s, h, t, u, gama, delta, v
INTEGER   :: i, iter, j, nit

DO  i = 1, p
  sc(i,1) = 0.0
END DO
nit = 0
20 nit = nit + 1
DO  i = 1, p
  sc(i,5) = sc(i,1)
END DO
DO  iter = 1, p
  DO  j = 1, n
    s = 0.0
    DO  i = 1, p
      s = s + sc(i,1) * tx(j,i)
    END DO
    r(j) = (ty(j)-s) * w(j)
  END DO
  DO  i = 1, p
    s = 0.0
    DO  j = 1, n
      s = s + r(j) * tx(j,i)
    END DO
    sc(i,2) = -2.0 * s / sw
  END DO
  s = 0.0
  DO  i = 1, p
    s = s + sc(i,2) ** 2
  END DO
  IF (s <= 0.0) EXIT
  IF (iter == 1) THEN
    DO  i = 1, p
      sc(i,3) = -sc(i,2)
    END DO
    h = s
  ELSE
    gama = s / h
    h = s
    DO  i = 1, p
      sc(i,3) = -sc(i,2) + gama * sc(i,4)
    END DO
  END IF
  s = 0.0
  t = s
  DO  j = 1, n
    u = 0.0
    DO  i = 1, p
      u = u + sc(i,3) * tx(j,i)
    END DO
    s = s + u * r(j)
    t = t + w(j) * u ** 2
  END DO
  delta = s / t
  DO  i = 1, p
    sc(i,1) = sc(i,1) + delta * sc(i,3)
    sc(i,4) = sc(i,3)
  END DO
END DO

v = 0.0
DO  i = 1, p
  v = MAX(v, ABS(sc(i,1)-sc(i,5)))
END DO
IF (v >= eps.AND.nit < maxit) THEN
  GO TO 20
END IF
DO  i = 1, p
  DO  j = 1, n
    tx(j,i) = sc(i,1) * tx(j,i)
  END DO
END DO
RETURN
END SUBROUTINE scale



SUBROUTINE sort(v, a, ii, jj)
 
!  Puts into a the permutation vector which sorts v into
!  increasing order.  Only elements from ii to jj are considered.
!  arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements

!  This is a modification of CACM algorithm #347 by R. C. Singleton,
!  which is a modified Hoare quicksort.

REAL, INTENT(IN OUT)     :: v(:)
INTEGER, INTENT(IN)      :: jj
INTEGER, INTENT(IN OUT)  :: a(jj)
INTEGER, INTENT(IN)      :: ii

INTEGER  :: iu(20), il(20)
INTEGER  :: i, ij, j, k, l, m, t, tt
REAL     :: vt, vtt

m = 1
i = ii
j = jj
10 IF (i >= j) GO TO 60
20 k = i
ij = (j+i) / 2
t = a(ij)
vt = v(ij)
IF (v(i) > vt) THEN
  a(ij) = a(i)
  a(i) = t
  t = a(ij)
  v(ij) = v(i)
  v(i) = vt
  vt = v(ij)
END IF
l = j
IF (v(j) >= vt) GO TO 40
a(ij) = a(j)
a(j) = t
t = a(ij)
v(ij) = v(j)
v(j) = vt
vt = v(ij)
IF (v(i) <= vt) GO TO 40
a(ij) = a(i)
a(i) = t
t = a(ij)
v(ij) = v(i)
v(i) = vt
vt = v(ij)
GO TO 40

30 a(l) = a(k)
a(k) = tt
v(l) = v(k)
v(k) = vtt

40 l = l - 1
IF (v(l) > vt) GO TO 40
tt = a(l)
vtt = v(l)
50 k = k + 1
IF (v(k) < vt) GO TO 50
IF (k <= l) GO TO 30
IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 70
END IF
il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 70
60 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)
70 IF (j-i > 10) GO TO 20
IF (i == ii) GO TO 10
i = i - 1

80 i = i + 1
IF (i == j) GO TO 60
t = a(i+1)
vt = v(i+1)
IF (v(i) <= vt) GO TO 80
k = i

90 a(k+1) = a(k)
v(k+1) = v(k)
k = k - 1
IF (vt < v(k)) GO TO 90
a(k+1) = t
v(k+1) = vt
GO TO 80
END SUBROUTINE sort



SUBROUTINE supsmu(n, x, y, w, iper, span, alpha, smo, sc)
 
!------------------------------------------------------------------

! super smoother (Friedman, 1984).

! version 10/10/84

! coded  and copywrite <c> 1984 by:

!                        jerome h. friedman
!                     department of statistics
!                               and
!                stanford linear accelerator center
!                        stanford university

! all rights reserved.


! input:
!    n : number of observations (x,y - pairs).
!    x(n) : ordered abscissa values.
!    y(n) : corresponding ordinate (response) values.
!    w(n) : weight for each (x,y) observation.
!    iper : periodic variable flag.
!       iper=1 => x is ordered interval variable.
!       iper=2 => x is a periodic variable with values
!                 in the range (0.0,1.0) and peroid 1.0.
!    span : smoother span (fraction of observations in window).
!           span=0.0 => automatic (variable) span selection.
!    alpha : controles high frequency (small span) penality
!            used with automatic span selection (bass tone control).
!            (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
! output:
!   smo(n) : smoothed ordinate (response) values.
! scratch:
!   sc(n,7) : internal working storage.

! note:
!    for small samples (n < 40) or if there are substantial serial
!    correlations between obserations close in x - value, then
!    a prespecified fixed span smoother (span > 0) should be
!    used. reasonable span values are 0.3 to 0.5.

!------------------------------------------------------------------

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: x(n)
REAL, INTENT(IN)      :: y(n)
REAL, INTENT(IN)      :: w(n)
INTEGER, INTENT(IN)   :: iper
REAL, INTENT(IN OUT)  :: span
REAL, INTENT(IN OUT)  :: alpha
REAL, INTENT(OUT)     :: smo(n)
REAL, INTENT(IN OUT)  :: sc(n,7)

! COMMON /spans/ spans(3) /consts/ big, sml, eps
INTEGER  :: i, j, jper
REAL     :: a, f, h(n), resmin, scale, sw, sy, vsmlsq

IF (x(n) <= x(1)) THEN
  sy = 0.0
  sw = sy
  DO  j = 1, n
    sy = sy + w(j) * y(j)
    sw = sw + w(j)
  END DO
  a = 0.0
  IF (sw > 0.0) a = sy / sw
  DO  j = 1, n
    smo(j) = a
  END DO
  RETURN
END IF
i = n / 4
j = 3 * i
scale = x(j) - x(i)
30 IF (scale <= 0.0) THEN
  IF (j < n) j = j + 1
  IF (i > 1) i = i - 1
  scale = x(j) - x(i)
  GO TO 30
END IF
vsmlsq = (eps*scale) ** 2
jper = iper
IF (iper == 2 .AND. (x(1) < 0.0 .OR. x(n) > 1.0)) jper = 1
IF (jper < 1 .OR. jper > 2) jper = 1
IF (span > 0.0) THEN
  CALL smooth(n, x, y, w, span, jper, vsmlsq, smo, sc)
  RETURN
END IF
DO  i = 1, 3
  CALL smooth(n, x, y, w, spans(i), jper, vsmlsq, sc(1:,2*i-1), sc(1:,7))
  CALL smooth(n, x, sc(1:,7), w, spans(2), -jper, vsmlsq, sc(1:,2*i), h)
END DO
DO  j = 1, n
  resmin = big
  DO  i = 1, 3
    IF (sc(j,2*i) < resmin) THEN
      resmin = sc(j,2*i)
      sc(j,7) = spans(i)
    END IF
  END DO
  IF (alpha > 0.0 .AND. alpha <= 10.0 .AND. resmin < sc(j,6) .AND.   &
      resmin > 0.0) sc(j,7) = sc(j,7) + (spans(3)-sc(j,7)) *  &
      MAX(sml,resmin/sc(j,6)) ** (10.0-alpha)
END DO
CALL smooth(n, x, sc(1,7), w, spans(2), -jper, vsmlsq, sc(1,2), h)
DO  j = 1, n
  IF (sc(j,2) <= spans(1)) sc(j,2) = spans(1)
  IF (sc(j,2) >= spans(3)) sc(j,2) = spans(3)
  f = sc(j,2) - spans(2)
  IF (f < 0.0) THEN
    f = -f / (spans(2)-spans(1))
    sc(j,4) = (1.0-f) * sc(j,3) + f * sc(j,1)
  ELSE
    f = f / (spans(3)-spans(2))
    sc(j,4) = (1.0-f) * sc(j,3) + f * sc(j,5)
  END IF
END DO
CALL smooth(n, x, sc(1,4), w, spans(1), -jper, vsmlsq, smo, h)
RETURN
END SUBROUTINE supsmu


SUBROUTINE smooth(n, x, y, w, span, iper, vsmlsq, smo, acvr)

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: x(n)
REAL, INTENT(IN)      :: y(n)
REAL, INTENT(IN)      :: w(n)
REAL, INTENT(IN)      :: span
INTEGER, INTENT(IN)   :: iper
REAL, INTENT(IN OUT)  :: vsmlsq
REAL, INTENT(OUT)     :: smo(n)
REAL, INTENT(OUT)     :: acvr(n)

INTEGER    :: i, ibw, in, it, j, j0, jper, out
REAL (dp)  :: wt, fbo, fbw, xm, ym, tmp, var, cvar, a, h, sy, xti, xto

xm = 0.0
ym = xm
var = ym
cvar = var
fbw = cvar
jper = ABS(iper)
ibw = 0.5 * span * n + 0.5
IF (ibw < 2) ibw = 2
it = 2 * ibw + 1
DO  i = 1, it
  j = i
  IF (jper == 2) j = i - ibw - 1
  xti = x(j)
  IF (j < 1) THEN
    j = n + j
    xti = x(j) - 1.0
  END IF
  wt = w(j)
  fbo = fbw
  fbw = fbw + wt
  IF (fbw > 0.0) xm = (fbo*xm+wt*xti) / fbw
  IF (fbw > 0.0) ym = (fbo*ym+wt*y(j)) / fbw
  tmp = 0.0
  IF (fbo > 0.0) tmp = fbw * wt * (xti-xm) / fbo
  var = var + tmp * (xti-xm)
  cvar = cvar + tmp * (y(j)-ym)
END DO
DO  j = 1, n
  out = j - ibw - 1
  in = j + ibw
  IF (.NOT.((jper /= 2) .AND. (out < 1 .OR. in > n))) THEN
    IF (out < 1) THEN
      out = n + out
      xto = x(out) - 1.0
      xti = x(in)
    ELSE
      IF (in > n) THEN
        in = in - n
        xti = x(in) + 1.0
        xto = x(out)
      ELSE
        xto = x(out)
        xti = x(in)
      END IF
    END IF
    wt = w(out)
    fbo = fbw
    fbw = fbw - wt
    tmp = 0.0
    IF (fbw > 0.0) tmp = fbo * wt * (xto-xm) / fbw
    var = var - tmp * (xto-xm)
    cvar = cvar - tmp * (y(out)-ym)
    IF (fbw > 0.0) xm = (fbo*xm-wt*xto) / fbw
    IF (fbw > 0.0) ym = (fbo*ym-wt*y(out)) / fbw
    wt = w(in)
    fbo = fbw
    fbw = fbw + wt
    IF (fbw > 0.0) xm = (fbo*xm+wt*xti) / fbw
    IF (fbw > 0.0) ym = (fbo*ym+wt*y(in)) / fbw
    tmp = 0.0
    IF (fbo > 0.0) tmp = fbw * wt * (xti-xm) / fbo
    var = var + tmp * (xti-xm)
    cvar = cvar + tmp * (y(in)-ym)
  END IF
  a = 0.0
  IF (var > vsmlsq) a = cvar / var
  smo(j) = a * (x(j)-xm) + ym
  IF (iper > 0) THEN
    h = 0.0
    IF (fbw > 0.0) h = 1.0 / fbw
    IF (var > vsmlsq) h = h + (x(j)-xm) ** 2 / var
    acvr(j) = 0.0
    a = 1.0 - w(j) * h
    IF (a > 0.0) THEN
      acvr(j) = ABS(y(j)-smo(j)) / a
    ELSE
      IF (j > 1) THEN
        acvr(j) = acvr(j-1)
      END IF
    END IF
  END IF
END DO
j = 1
30 j0 = j
sy = smo(j) * w(j)
fbw = w(j)
IF (j < n) THEN
  40   IF (x(j+1) <= x(j)) THEN
    j = j + 1
    sy = sy + w(j) * smo(j)
    fbw = fbw + w(j)
    IF (j < n) THEN
      GO TO 40
    END IF
  END IF
END IF
IF (j > j0) THEN
  a = 0.0
  IF (fbw > 0.0) a = sy / fbw
  DO  i = j0, j
    smo(i) = a
  END DO
END IF
j = j + 1
IF (j <= n) THEN
  GO TO 30
END IF
RETURN
END SUBROUTINE smooth

END MODULE Ace
