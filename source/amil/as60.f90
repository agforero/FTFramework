MODULE AS60

! Calculate the eigenvalues and eigenvectors of a real symmetric matrix.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

CONTAINS


SUBROUTINE tdiag (n, tol, a, d, e, z, maxdim, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-03  Time: 12:16:20

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(OUT)   :: tol
INTEGER, INTENT(IN OUT)  :: maxdim
REAL (dp), INTENT(IN)    :: a(maxdim,maxdim)
REAL (dp), INTENT(OUT)   :: d(maxdim)
REAL (dp), INTENT(OUT)   :: e(maxdim)
REAL (dp), INTENT(OUT)   :: z(maxdim,maxdim)
INTEGER, INTENT(OUT)     :: ifault

REAL (dp), PARAMETER  :: eta = 1.0D-37, eps = 1.0D-14

!    Algorithm as 60.1 appl.statist. (1973) vol.22 no.2

!    Reduces real symmetric matrix to tridiagonal form

!    tol is a machine dependent constant , tol = eta/eps , where
!    eta is the smallest positive number representable in the
!    computer and eps is the smallest positive number for which
!    1+eps.ne.1.

!      eta=eps*tol
!      eps=0.7105427358e-14
!      tol=0.3131513063e-293
!      precis=1.0e-14

!      NB
!        Real constants must be <= 15 decimal digits
!        the range of a real constant is from 1.0e-293 to 1.0e+322

INTEGER    :: i, i1, j, j1, k, l
REAL (dp)  :: f, g, h, hh

tol=eta/eps
ifault=1
IF (n <= 1) RETURN
ifault=0
DO  i=1,n
  z(i,1:i)=a(i,1:i)
END DO
i=n
DO  i1=2,n
  l=i-2
  f=z(i,i-1)
  g=zero
  IF (l < 1) GO TO 30
  g = SUM( z(i,1:l)**2 )
  30 h=g + f*f
  
!       If g is too small for orthogonality to be guaranteed, the
!       transformation is skipped
  
  IF (g > tol) GO TO 40
  e(i)=f
  d(i)=zero
  GO TO 100

  40 l=l+1
  g=SQRT(h)
  IF (f >= zero) g=-g
  e(i)=g
  h=h - f*g
  z(i,i-1)=f-g
  f=zero
  DO  j=1,l
    z(j,i)=z(i,j)/h
    g=zero
    
!       Form element of a * u
    
    DO  k=1,j
      g=g + z(j,k)*z(i,k)
    END DO
    IF (j >= l) GO TO 70
    j1=j+1
    DO  k=j1,l
      g=g + z(k,j)*z(i,k)
    END DO
    
!       Form element of p
    
    70 e(j)=g/h
    f=f + g*z(j,i)
  END DO
  
!       Form k
  
  hh=f/(h+h)
  
!       Form reduced a
  
  DO  j=1,l
    f=z(i,j)
    g=e(j)-hh*f
    e(j)=g
    DO  k=1,j
      z(j,k)=z(j,k) - f*e(k) - g*z(i,k)
    END DO
  END DO
  d(i)=h
  100 i=i-1
END DO
d(1)=zero
e(1)=zero

!       Accumulation of transformation matrices

DO  i=1,n
  l=i-1
  IF (d(i) == zero .OR. l == 0) GO TO 140
  DO  j=1,l
    g=zero
    DO  k=1,l
      g=g + z(i,k)*z(k,j)
    END DO
    DO  k=1,l
      z(k,j)=z(k,j) - g*z(k,i)
    END DO
  END DO
  140 d(i)=z(i,i)
  z(i,i)=one
  IF (l == 0) CYCLE
  DO  j=1,l
    z(i,j)=zero
    z(j,i)=zero
  END DO
END DO
RETURN
END SUBROUTINE tdiag


SUBROUTINE lrvt (n, precis, d, e, z, ifault, maxdim)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: precis
INTEGER, INTENT(IN)        :: maxdim
REAL (dp), INTENT(IN OUT)  :: d(maxdim)
REAL (dp), INTENT(OUT)     :: e(maxdim)
REAL (dp), INTENT(IN OUT)  :: z(maxdim,maxdim)
INTEGER, INTENT(OUT)       :: ifault

!    Algorithm AS 60.2 appl.statist. (1973) vol.22, no.2

!    Finds latent roots and vectors of tridiagonal matrix

INTEGER, PARAMETER    :: mits = 30
REAL (dp), PARAMETER  :: two = 2.0D0

INTEGER    :: i, i1, j, jj, k, l, m, m1, n1
REAL (dp)  :: b, c, f, g, h, p, pr, r, s

precis=1.0D-14
ifault=2
IF (n <= 1) RETURN
ifault=1
n1=n-1
DO  i=2,n
  e(i-1)=e(i)
END DO
e(n)=zero
b=zero
f=zero
DO  l=1,n
  jj=0
  h=precis*(ABS(d(l))+ABS(e(l)))
  IF (b < h) b=h
  
!       Look for small sub-diagonal element
  
  DO  m1=l,n
    m=m1
    IF (ABS(e(m)) <= b) EXIT
  END DO
  IF (m == l) GO TO 90
  40 IF (jj == mits) RETURN
  jj=jj+1
  
!       Form shift
  
  p=(d(l+1)-d(l))/(two*e(l))
  r=SQRT(p*p + one)
  pr=p+r
  IF (p < zero) pr=p-r
  h=d(l)-e(l)/pr
  DO  i=l,n
    d(i)=d(i)-h
  END DO
  f=f+h
  
!       QL transformation
  
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
    IF (ABS(p) >= ABS(e(i))) GO TO 60
    c=p/e(i)
    r=SQRT(c*c + one)
    e(j)=s*e(i)*r
    s=one/r
    c=c/r
    GO TO 70
    60 c=e(i)/p
    r=SQRT(c*c + one)
    e(j)=s*p*r
    s=c/r
    c=one/r
    70 p=c*d(i) - s*g
    d(j)=h + s*(c*g + s*d(i))
    
!       Form vector
    
    DO  k=1,n
      h=z(k,j)
      z(k,j)=s*z(k,i) + c*h
      z(k,i)=c*z(k,i) - s*h
    END DO
  END DO
  e(l)=s*p
  d(l)=c*p
  IF (ABS(e(l)) > b) GO TO 40
  90 d(l)=d(l) + f
END DO

!       Order latent roots and vectors

DO  i=1,n1
  k=i
  p=d(i)
  i1=i+1
  DO  j=i1,n
    IF (d(j) <= p) CYCLE
    k=j
    p=d(j)
  END DO
  IF (k == i) CYCLE
  d(k)=d(i)
  d(i)=p
  DO  j=1,n
    p=z(j,i)
    z(j,i)=z(j,k)
    z(j,k)=p
  END DO
END DO
ifault=0
RETURN
END SUBROUTINE lrvt

END MODULE AS60
