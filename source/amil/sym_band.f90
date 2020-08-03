MODULE Eigen_symmetric_banded
! Calculate eigenvalues & vectors of symmetric banded matrices

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-18  Time: 23:59:29

IMPLICIT NONE
INTEGER, PARAMETER             :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), PARAMETER, PRIVATE  :: one = 1.0_dp, zero = 0.0_dp

CONTAINS


SUBROUTINE rsb(nm, n, mb, a, w, matz, z, ierr)

! N.B. Arguments FV1 & FV2 have been removed.
 
INTEGER, INTENT(IN)        :: nm
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mb
REAL (dp), INTENT(IN OUT)  :: a(:,:)   ! a(nm,mb)
REAL (dp), INTENT(IN OUT)  :: w(:)
INTEGER, INTENT(IN)        :: matz
REAL (dp), INTENT(IN OUT)  :: z(:,:)   ! z(nm,n)
INTEGER, INTENT(OUT)       :: ierr

! Local variables

LOGICAL    :: tf
REAL (dp)  :: fv1(n), fv2(n)

!  This subroutine calls the recommended sequence of subroutines from the
!  eigensystem subroutine package (eispack) to find the eigenvalues and
!  eigenvectors (if desired) of a real symmetric band matrix.

!  on input

!     nm  must be set to the row dimension of the two-dimensional array
!     parameters as declared in the calling program dimension statement.

!     n  is the order of the matrix  a.

!     mb  is the half band width of the matrix, defined as the number of
!     adjacent diagonals, including the principal diagonal, required to
!     specify the non-zero portion of the lower triangle of the matrix.

!     a  contains the lower triangle of the real symmetric band matrix.
!     Its lowest subdiagonal is stored in the last  n+1-mb  positions of the
!     first column, its next subdiagonal in the last  n+2-mb  positions of the
!     second column, further subdiagonals similarly, and finally its principal
!     diagonal in the  n  positions of the last column.
!     Contents of storages not part of the matrix are arbitrary.

!     matz  is an integer variable set equal to zero if only eigenvalues are
!     desired.   Otherwise it is set to any non-zero integer for both
!     eigenvalues and eigenvectors.

!  on output

!     w  contains the eigenvalues in ascending order.

!     z  contains the eigenvectors if matz is not zero.

!     ierr  is an integer output variable set equal to an error completion
!        code described in the documentation for tqlrat and tql2.
!        The normal completion code is zero.

!     fv1  and  fv2  are temporary storage arrays.

!  questions and comments should be directed to burton s. garbow,
!  mathematics and computer science div, argonne national laboratory

!  this version dated august 1983.

!  ------------------------------------------------------------------

IF (n > nm) THEN
  ierr = 10 * n
ELSE
  IF (mb <= 0) THEN
    ierr = 12 * n
  ELSE
    IF (mb > n) THEN
      ierr = 12 * n
    ELSE
      
      IF (matz == 0) THEN
!     .......... find eigenvalues only ..........
        tf = .false.
        CALL bandr(n, mb, a, w, fv1, fv2, tf, z)
        CALL tqlrat(n, w, fv2, ierr)
      ELSE
!     .......... find both eigenvalues and eigenvectors ..........
        tf = .true.
        CALL bandr(n, mb, a, w, fv1, fv1, tf, z)
        CALL tql2(n, w, fv1, z, ierr)
      END IF
    END IF
  END IF
END IF

RETURN
END SUBROUTINE rsb



SUBROUTINE bandr(n, mb, a, d, e, e2, matz, z)

! N.B. Argument NM has been removed.
 
!  REFORMULATED S2 IN LOOP 500 TO AVOID OVERFLOW. (9/29/89 BSG)

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mb
REAL (dp), INTENT(IN OUT)  :: a(:,:)   ! a(nm,mb)
REAL (dp), INTENT(OUT)     :: d(:)
REAL (dp), INTENT(OUT)     :: e(:)
REAL (dp), INTENT(OUT)     :: e2(:)
LOGICAL, INTENT(IN)        :: matz
REAL (dp), INTENT(OUT)     :: z(:,:)   ! z(nm,n)

! Local variables

INTEGER    :: j, k, l, r, i1, i2, j1, j2, kr, mr, m1, n2, r1, ugl, maxl, maxr
REAL (dp)  :: g, u, b1, b2, c2, f1, f2, s2, dmin, dminrt


!  this subroutine is a translation of the algol procedure bandrd,
!  num. math. 12, 231-241(1968) by schwarz.
!  handbook for auto. comp., vol.ii-linear algebra, 273-283(1971).

!  this subroutine reduces a real symmetric band matrix
!  to a symmetric tridiagonal matrix using and optionally
!  accumulating orthogonal similarity transformations.

!  on input

!     nm must be set to the row dimension of two-dimensional array
!       parameters as declared in the calling program dimension statement.

!     n is the order of the matrix.

!     mb is the (half) band width of the matrix, defined as the number of
!       adjacent diagonals, including the principal diagonal, required to
!       specify the non-zero portion of the lower triangle of the matrix.

!     a contains the lower triangle of the symmetric band input
!       matrix stored as an n by mb array.  its lowest subdiagonal
!       is stored in the last n+1-mb positions of the first column,
!       its next subdiagonal in the last n+2-mb positions of the
!       second column, further subdiagonals similarly, and finally
!       its principal diagonal in the n positions of the last column.
!       contents of storages not part of the matrix are arbitrary.

!     matz should be set to .true. if the transformation matrix is
!       to be accumulated, and to .false. otherwise.

!  on output

!     a has been destroyed, except for its last two columns which
!       contain a copy of the tridiagonal matrix.

!     d contains the diagonal elements of the tridiagonal matrix.

!     e contains the subdiagonal elements of the tridiagonal
!       matrix in its last n-1 positions.  e(1) is set to zero.

!     e2 contains the squares of the corresponding elements of e.
!       e2 may coincide with e if the squares are not needed.

!     z contains the orthogonal transformation matrix produced in the
!       reduction if matz has been set to .true.  Otherwise, z is not
!       referenced.

!  questions and comments should be directed to burton s. garbow,
!  mathematics and computer science div, argonne national laboratory

!  this version dated september 1989.

!  ------------------------------------------------------------------

dmin = EPSILON(one)
dminrt = SQRT(dmin)
!     .......... initialize diagonal scaling matrix ..........
d(1:n) = one

IF (matz) THEN
  
  DO  j = 1, n
    
    z(j,1:n) = zero
    
    z(j,j) = one
  END DO
END IF

m1 = mb - 1
IF (m1 < 1) THEN
  GO TO 230
ELSE IF (m1 == 1) THEN
  GO TO 180
END IF
n2 = n - 2

DO  k = 1, n2
  maxr = MIN(m1,n-k)
!     .......... for r=maxr step -1 until 2 do -- ..........
  loop120:  DO  r1 = 2, maxr
    r = maxr + 2 - r1
    kr = k + r
    mr = mb - r
    g = a(kr,mr)
    a(kr-1,1) = a(kr-1,mr+1)
    ugl = k
    
    DO  j = kr, n, m1
      j1 = j - 1
      j2 = j1 - 1
      IF (g == zero) CYCLE loop120
      b1 = a(j1,1) / g
      b2 = b1 * d(j1) / d(j)
      IF (ABS(b1) > one) THEN
        u = one / b1
        s2 = u / (u+b2)
      ELSE
        s2 = one / (one+b1*b2)
      END IF
      
      IF (s2 < 0.5D0) THEN
        b1 = g / a(j1,1)
        b2 = b1 * d(j) / d(j1)
        c2 = one - s2
        d(j1) = c2 * d(j1)
        d(j) = c2 * d(j)
        f1 = 2.0D0 * a(j,m1)
        f2 = b1 * a(j1,mb)
        a(j,m1) = -b2 * (b1*a(j,m1) - a(j,mb)) - f2 + a(j,m1)
        a(j1,mb) = b2 * (b2*a(j,mb) + f1) + a(j1,mb)
        a(j,mb) = b1 * (f2-f1) + a(j,mb)
        
        DO  l = ugl, j2
          i2 = mb - j + l
          u = a(j1,i2+1) + b2 * a(j,i2)
          a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
          a(j1,i2+1) = u
        END DO
        
        ugl = j
        a(j1,1) = a(j1,1) + b2 * g
        IF (j /= n) THEN
          maxl = MIN(m1,n-j1)
          
          DO  l = 2, maxl
            i1 = j1 + l
            i2 = mb - l
            u = a(i1,i2) + b2 * a(i1,i2+1)
            a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
            a(i1,i2) = u
          END DO
          
          i1 = j + m1
          IF (i1 <= n) THEN
            g = b2 * a(i1,1)
          END IF
        END IF
        IF (.NOT.matz) CYCLE
        
        DO  l = 1, n
          u = z(l,j1) + b2 * z(l,j)
          z(l,j) = -b1 * z(l,j1) + z(l,j)
          z(l,j1) = u
        END DO
        
      ELSE
        
        u = d(j1)
        d(j1) = s2 * d(j)
        d(j) = s2 * u
        f1 = 2.0D0 * a(j,m1)
        f2 = b1 * a(j,mb)
        u = b1 * (f2-f1) + a(j1,mb)
        a(j,m1) = b2 * (b1*a(j,m1) - a(j1,mb)) + f2 - a(j,m1)
        a(j1,mb) = b2 * (b2*a(j1,mb) + f1) + a(j,mb)
        a(j,mb) = u
        
        DO  l = ugl, j2
          i2 = mb - j + l
          u = b2 * a(j1,i2+1) + a(j,i2)
          a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
          a(j1,i2+1) = u
        END DO
        
        ugl = j
        a(j1,1) = b2 * a(j1,1) + g
        IF (j /= n) THEN
          maxl = MIN(m1,n-j1)
          
          DO  l = 2, maxl
            i1 = j1 + l
            i2 = mb - l
            u = b2 * a(i1,i2) + a(i1,i2+1)
            a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
            a(i1,i2) = u
          END DO
          
          i1 = j + m1
          IF (i1 <= n) THEN
            g = a(i1,1)
            a(i1,1) = b1 * a(i1,1)
          END IF
        END IF
        IF (matz) THEN
          
          DO  l = 1, n
            u = b2 * z(l,j1) + z(l,j)
            z(l,j) = -z(l,j1) + b1 * z(l,j)
            z(l,j1) = u
          END DO
        END IF
      END IF
      
    END DO
    
  END DO loop120
  
  IF (MOD(k,64) == 0) THEN
!     .......... rescale to avoid underflow or overflow ..........
    DO  j = k, n
      IF (d(j) < dmin) THEN
        maxl = MAX(1,mb+1-j)
        
        DO  l = maxl, m1
          a(j,l) = dminrt * a(j,l)
        END DO
        
        IF (j /= n) THEN
          maxl = MIN(m1,n-j)
          
          DO  l = 1, maxl
            i1 = j + l
            i2 = mb - l
            a(i1,i2) = dminrt * a(i1,i2)
          END DO
        END IF
        
        IF (matz) THEN
          
          DO  l = 1, n
            z(l,j) = dminrt * z(l,j)
          END DO
        END IF
        
        a(j,mb) = dmin * a(j,mb)
        d(j) = d(j) / dmin
      END IF
    END DO
  END IF
  
END DO
!     .......... form square root of scaling matrix ..........
180 DO  j = 2, n
  e(j) = SQRT(d(j))
END DO

IF (matz) THEN
  
  DO  j = 1, n
    z(j,2:n) = e(2:n) * z(j,2:n)
  END DO
END IF

u = one

DO  j = 2, n
  a(j,m1) = u * e(j) * a(j,m1)
  u = e(j)
  e2(j) = a(j,m1) ** 2
  a(j,mb) = d(j) * a(j,mb)
  d(j) = a(j,mb)
  e(j) = a(j,m1)
END DO

d(1) = a(1,mb)
e(1) = zero
e2(1) = zero
GO TO 250

230 DO  j = 1, n
  d(j) = a(j,mb)
  e(j) = zero
  e2(j) = zero
END DO

250 RETURN
END SUBROUTINE bandr



SUBROUTINE tqlrat(n, d, e2, ierr)
 
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: d(:)
REAL (dp), INTENT(OUT)     :: e2(:)
INTEGER, INTENT(OUT)       :: ierr

! Local variables

INTEGER    :: i, j, l, m, ii, l1, mml
REAL (dp)  :: b, c, f, g, h, p, r, s, t

!  this subroutine is a translation of the algol procedure tqlrat,
!  algorithm 464, comm. acm 16, 689(1973) by reinsch.

!  this subroutine finds the eigenvalues of a symmetric
!  tridiagonal matrix by the rational ql method.

!  on input

!     n is the order of the matrix.

!     d contains the diagonal elements of the input matrix.

!     e2 contains the squares of the subdiagonal elements of the
!       input matrix in its last n-1 positions.  e2(1) is arbitrary.

!   on output

!     d contains the eigenvalues in ascending order.  if an
!       error exit is made, the eigenvalues are correct and
!       ordered for indices 1,2,...ierr-1, but may not be
!       the smallest eigenvalues.

!     e2 has been destroyed.

!     ierr is set to
!       zero       for normal return,
!       j          if the j-th eigenvalue has not been
!                  determined after 30 iterations.

!  calls pythag for  SQRT(a*a + b*b) .

!  questions and comments should be directed to burton s. garbow,
!  mathematics and computer science div, argonne national laboratory

!  this version dated august 1983.

!  ------------------------------------------------------------------

ierr = 0
IF (n /= 1) THEN
  
  DO  i = 2, n
    e2(i-1) = e2(i)
  END DO
  
  f = zero
  t = zero
  e2(n) = zero
  
  DO  l = 1, n
    j = 0
    h = ABS(d(l)) + SQRT(e2(l))
    IF (t <= h) THEN
      t = h
      b = ABS(t) * EPSILON(one)
      c = b * b
    END IF
!     .......... look for small squared sub-diagonal element ..........
    DO  m = l, n
      IF (e2(m) <= c) GO TO 30
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
    END DO
    
    30 IF (m /= l) THEN
      40 IF (j == 30) GO TO 100
      j = j + 1
!     .......... form shift ..........
      l1 = l + 1
      s = SQRT(e2(l))
      g = d(l)
      p = (d(l1)-g) / (2.0D0*s)
      r = pythag(p, one)
      d(l) = s / (p + SIGN(r,p))
      h = g - d(l)
      
      DO  i = l1, n
        d(i) = d(i) - h
      END DO
      
      f = f + h
!     .......... rational ql transformation ..........
      g = d(m)
      IF (g == zero) g = b
      h = g
      s = zero
      mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
      DO  ii = 1, mml
        i = m - ii
        p = g * h
        r = p + e2(i)
        e2(i+1) = s * r
        s = e2(i) / r
        d(i+1) = h + s * (h+d(i))
        g = d(i) - e2(i) / g
        IF (g == zero) g = b
        h = g * p / r
      END DO
      
      e2(l) = s * g
      d(l) = h
!     .......... guard against underflow in convergence test ..........
      IF (h /= zero) THEN
        IF (ABS(e2(l)) > ABS(c/h)) THEN
          e2(l) = h * e2(l)
          IF (e2(l) /= zero) GO TO 40
        END IF
      END IF
    END IF
    p = d(l) + f
!     .......... order eigenvalues ..........
    IF (l /= 1) THEN
!     .......... for i=l step -1 until 2 do -- ..........
      DO  ii = 2, l
        i = l + 2 - ii
        IF (p >= d(i-1)) GO TO 80
        d(i) = d(i-1)
      END DO
    END IF
    
    i = 1
    80 d(i) = p
  END DO
  
  GO TO 110
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
  100 ierr = l
END IF

110 RETURN
END SUBROUTINE tqlrat



FUNCTION pythag(a, b) RESULT(fn_val)
 
REAL (dp), INTENT(IN)  :: a, b
REAL (dp)              :: fn_val

!     finds sqrt(a**2 + b**2) without overflow or destructive underflow

REAL (dp) :: r, s, t, u

fn_val = MAX(ABS(a), ABS(b))
IF (fn_val /= zero) THEN
  r = (MIN(ABS(a), ABS(b))/fn_val) ** 2
  DO
    t = 4.0_dp + r
    IF (t /= 4.0_dp) THEN
      s = r / t
      u = one + 2.0_dp * s
      fn_val = u * fn_val
      r = (s/u) ** 2 * r
      CYCLE
    END IF
    EXIT
  END DO
END IF

RETURN
END FUNCTION pythag



SUBROUTINE tql2(n, d, e, z, ierr)

! N.B. Argument NM has been removed.
 
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: d(:)
REAL (dp), INTENT(OUT)     :: e(:)
REAL (dp), INTENT(IN OUT)  :: z(:,:)   ! z(nm,n)
INTEGER, INTENT(OUT)       :: ierr

! Local variables

INTEGER    :: i, j, k, l, m, ii, l1, l2, mml
REAL (dp)  :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2

!  this subroutine is a translation of the algol procedure tql2,
!  num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and wilkinson.
!  handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).

!  this subroutine finds the eigenvalues and eigenvectors of a symmetric
!  tridiagonal matrix by the ql method.  The eigenvectors of a full symmetric
!  matrix can also be found if  tred2  has been used to reduce this full
!  matrix to tridiagonal form.

!  on input

!     nm must be set to the row dimension of two-dimensional array parameters
!       as declared in the calling program dimension statement.

!     n is the order of the matrix.

!     d contains the diagonal elements of the input matrix.

!     e contains the subdiagonal elements of the input matrix
!       in its last n-1 positions.  e(1) is arbitrary.

!     z contains the transformation matrix produced in the reduction by
!       tred2, if performed.  If the eigenvectors of the tridiagonal matrix
!       are desired, z must contain the identity matrix.

!   on output

!     d contains the eigenvalues in ascending order.  If an error exit is made,
!       the eigenvalues are correct but unordered for indices 1,2,...,ierr-1.

!     e has been destroyed.

!     z contains orthonormal eigenvectors of the symmetric
!       tridiagonal (or full) matrix.  if an error exit is made,
!       z contains the eigenvectors associated with the stored eigenvalues.

!     ierr is set to
!       zero       for normal return,
!       j          if the j-th eigenvalue has not been
!                  determined after 30 iterations.

!  calls pythag for  SQRT(a*a + b*b) .

!  questions and comments should be directed to burton s. garbow,
!  mathematics and computer science div, argonne national laboratory

!  this version dated august 1983.

!  ------------------------------------------------------------------

ierr = 0
IF (n /= 1) THEN
  
  DO  i = 2, n
    e(i-1) = e(i)
  END DO
  
  f = zero
  tst1 = zero
  e(n) = zero
  
  DO  l = 1, n
    j = 0
    h = ABS(d(l)) + ABS(e(l))
    IF (tst1 < h) tst1 = h
!     .......... look for small sub-diagonal element ..........
    DO  m = l, n
      tst2 = tst1 + ABS(e(m))
      IF (tst2 == tst1) GO TO 30
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
    END DO
    
    30 IF (m /= l) THEN
      40 IF (j == 30) GO TO 120
      j = j + 1
!     .......... form shift ..........
      l1 = l + 1
      l2 = l1 + 1
      g = d(l)
      p = (d(l1)-g) / (2.0_dp*e(l))
      r = pythag(p,one)
      d(l) = e(l) / (p + SIGN(r,p))
      d(l1) = e(l) * (p + SIGN(r,p))
      dl1 = d(l1)
      h = g - d(l)
      IF (l2 <= n) THEN
        
        DO  i = l2, n
          d(i) = d(i) - h
        END DO
      END IF
      
      f = f + h
!     .......... ql transformation ..........
      p = d(m)
      c = one
      c2 = c
      el1 = e(l1)
      s = zero
      mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
      DO  ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c*g+s*d(i))
!     .......... form vector ..........
        DO  k = 1, n
          h = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * h
          z(k,i) = c * z(k,i) - s * h
        END DO
        
      END DO
      
      p = -s * s2 * c3 * el1 * e(l) / dl1
      e(l) = s * p
      d(l) = c * p
      tst2 = tst1 + ABS(e(l))
      IF (tst2 > tst1) GO TO 40
    END IF
    d(l) = d(l) + f
  END DO
!     .......... order eigenvalues and eigenvectors ..........
  DO  ii = 2, n
    i = ii - 1
    k = i
    p = d(i)
    
    DO  j = ii, n
      IF (d(j) < p) THEN
        k = j
        p = d(j)
      END IF
    END DO
    
    IF (k /= i) THEN
      d(k) = d(i)
      d(i) = p
      
      DO  j = 1, n
        p = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = p
      END DO
    END IF
    
  END DO
  
  GO TO 130
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
  120 ierr = l
END IF

130 RETURN
END SUBROUTINE tql2

END MODULE Eigen_symmetric_banded



PROGRAM Test_Eigen_symmetric_banded
USE Eigen_symmetric_banded
IMPLICIT NONE

INTEGER    :: i, ierr, mb, nrows
REAL (dp)  :: sumeig
REAL (dp), ALLOCATABLE  :: a(:,:), eigenval(:), eigenvector(:,:)

WRITE(*, '(a)', ADVANCE='NO') ' Enter number of rows in matrix: '
READ(*, *) nrows
WRITE(*, '(a)', ADVANCE='NO') ' Enter number of bands on and below diagonal: '
READ(*, *) mb
IF (mb >= nrows/2) THEN
  WRITE(*, *) ' ** Number of bands must be less than half the number of rows'
  STOP
END IF

ALLOCATE( a(nrows,mb), eigenval(nrows), eigenvector(nrows,nrows) )

DO i = 1, mb
  a(1:nrows,i) = i
END DO

CALL rsb(nrows, nrows, mb, a, eigenval, 1, eigenvector, ierr)
WRITE(*, *) 'IERR = ', ierr
sumeig = SUM(eigenval)
WRITE(*, '(a, g16.8)') ' Sum of eigenvalues = ', sumeig

DO i = 1, nrows
  WRITE(*, '(a, i3, a, g16.8)') ' Eigenvalue', i, ' = ', eigenval(i)
  IF (ABS(eigenval(i)) < 1.0e-10_dp) CYCLE
  WRITE(*, *) 'Eigenvector'
  WRITE(*, '(" ", 6g13.5)') eigenvector(:,i)
  WRITE(*, *)
END DO

STOP
END PROGRAM Test_Eigen_symmetric_banded
