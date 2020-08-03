SUBROUTINE chol (a, n, nn, u, nullty, ifault)

!    Algorithm AS6, Applied Statistics, vol.17, (1968)

!    Given a symmetric matrix order n as lower triangle in a( ),
!    calculates an upper triangle, u( ), such that uprime * u = a.
!    a must be positive semi-definite.  eta is set to multiplying
!    factor determining effective zero for pivot.

!    arguments:-
!    a()     = input, a +ve definite matrix stored in lower-triangular form.
!    n       = input, the order of a
!    nn      = input, the size of the a and u arrays >= n*(n+1)/2
!    u()     = output, a lower triangular matrix such that u*u' = a.
!              a & u may occupy the same locations.
!    nullty  = output, the rank deficiency of a.
!    ifault  = output, error indicator
!                    = 1 if n < 1
!                    = 2 if a is not +ve semi-definite
!                    = 3 if nn < n*(n+1)/2
!                    = 0 otherwise

!***********************************************************************

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)  :: a(:)
INTEGER, INTENT(IN)    :: n, nn
REAL (dp), INTENT(OUT) :: u(:)
INTEGER, INTENT(OUT)   :: nullty, ifault

! Local variables
REAL (dp)              :: eta2, x, w
INTEGER                :: i, icol, ii, irow, j, k, kk, l, m

!       The value of eta will depend on the word-length of the
!       computer being used.  See introductory text.

REAL (dp), PARAMETER   :: eta = 1.E-09_dp, zero = 0.0_dp

ifault = 1
IF (n <=  0) RETURN
ifault = 3
IF (nn < n*(n+1)/2) RETURN
ifault = 2
nullty = 0
j = 1
k = 0
eta2 = eta*eta
ii = 0

!       Factorize column by column, icol = column no.

DO icol = 1, n
  ii = ii + icol
  x = eta2*a(ii)
  l = 0
  kk = 0
  
!       IROW = row number within column ICOL
  
  DO irow = 1, icol
    kk = kk + irow
    k = k + 1
    w = a(k)
    m = j
    DO i = 1, irow
      l = l + 1
      IF (i == irow) EXIT
      w = w - u(l)*u(m)
      m = m + 1
    END DO

    IF (irow == icol) EXIT
    IF (u(l) /= zero) THEN
      u(k) = w / u(l)
    ELSE
      IF (w*w > ABS(x*a(kk))) RETURN
      u(k) = zero
    END IF
  END DO

  IF (ABS(w) > ABS(eta*a(k))) THEN
    IF (w < zero) RETURN
    u(k) = SQRT(w)
  ELSE
    u(k) = zero
    nullty = nullty + 1
  END IF
  j = j + icol
END DO

ifault = 0
RETURN
END SUBROUTINE chol




SUBROUTINE subchl (a, b, n, u, nullty, ifault, det)

!   REMARK ASR 44  APPL. STATIST. (1982) VOL. 31, NO. 3

!   A revised and enhanced version of
!   ALGORITHM AS 6  APPL. STATIST. (1968) VOL. 17, NO. 2

! N.B. Argument NDIM has been removed, and DET is an OPTIONAL argument

!   Given a symmetric matrix of order N as lower triangle in A(),
!   calculates an upper triangle, U(), such that U'U = the sub-matrix
!   of A whose rows and columns are specified in the integer array B().
!   U() may coincide with A().   A() must be +ve semi-definite.
!   ETA is set to multiplying factor determining effective zero for a pivot.
!   NULLTY is returned as number of effective zero pivots.
!   IFAULT is returned as 1 if N <= 0, 2 if A() is not +ve semi-
!   definite, otherwise 0 is returned.

IMPLICIT NONE
INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)            :: a(:)
REAL (dp), INTENT(OUT)           :: u(:)
REAL (dp), INTENT(OUT), OPTIONAL :: det
INTEGER, INTENT(IN)              :: b(:), n
INTEGER, INTENT(OUT)             :: nullty, ifault

!     Local variables

REAL (dp) :: w, eta2, x
INTEGER   :: i, icol, ii, ij, irow, j, jj, k, kk, l, m

!     The value of ETA below will depend upon the word length of the
!     computer being used.

REAL (dp), PARAMETER :: eta = 1.E-14_dp, zero = 0.0_dp, one = 1.0_dp

ifault = 1
IF (n <= 0) GO TO 90
ifault = 2
nullty = 0
IF (PRESENT(det)) det = one
j = 1
k = 0
eta2 = eta*eta
DO icol = 1, n
  ij = b(icol)*(b(icol) - 1)/2
  ii = ij + b(icol)
  x = eta2*a(ii)
  l = 0
  DO irow = 1, icol
    kk = b(irow)*(b(irow) + 1)/2
    k = k + 1
    jj = ij + b(irow)
    w = a(jj)
    m = j
    DO i = 1, irow
      l = l + 1
      IF (i == irow) GO TO 20
      w = w - u(l)*u(m)
      m = m + 1
    END DO
    20 IF (irow == icol) GO TO 50
    IF (u(l) == zero) GO TO 30
    u(k) = w/u(l)
    CYCLE
    30 IF (w*w > ABS(x*a(kk))) GO TO 90
    u(k) = zero
  END DO

  50 IF (ABS(w) <= ABS(eta*a(kk))) GO TO 60
  IF (w < zero) GO TO 90
  u(k) = SQRT(w)
  GO TO 70
  60 u(k) = zero
  nullty = nullty + 1
  70 j = j + icol
  IF (PRESENT(det)) det = det*u(k)*u(k)
END DO

ifault = 0
90 RETURN
END SUBROUTINE subchl
