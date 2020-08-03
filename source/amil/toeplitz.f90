SUBROUTINE dtoplx (a, b, x, n, g, h, ierr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-17  Time: 09:28:46
 
!-----------------------------------------------------------------------
!            SOLUTION OF THE TOEPLITZ SYSTEM OF EQUATIONS

!               SUM(J = 1,...,N) A(N+I-J)*X(J) = B(I)

!     FOR I = 1,...,N.
!-----------------------------------------------------------------------
!     DOUBLE PRECISION A(2*N - 1)
!----------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

REAL (dp), INTENT(IN)   :: a(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: b(n)
REAL (dp), INTENT(OUT)  :: x(n)
REAL (dp), INTENT(OUT)  :: g(n)
REAL (dp), INTENT(OUT)  :: h(n)
INTEGER, INTENT(OUT)    :: ierr

INTEGER    :: j, k, l, m, MAX, mp1, nml, npl
REAL (dp)  :: c, c1, c2, gd, gj, gk, gn, hj, hk, hn, xd, xn

IF (a(n) == 0.d0) GO TO 100
ierr = 0
x(1) = b(1)/a(n)
IF (n == 1) RETURN
g(1) = a(n - 1)/a(n)
h(1) = a(n + 1)/a(n)
mp1 = 1

!     COMPUTE NUMERATOR AND DENOMINATOR OF X(M+1)

10 m = mp1
mp1 = m + 1
xn = -b(mp1)
xd = -a(n)
DO  j = 1,m
  l = mp1 - j
  npl = n + l
  xn = xn + a(npl)*x(j)
  xd = xd + a(npl)*g(l)
END DO
IF (xd == 0.d0) GO TO 100
x(mp1) = xn/xd

!        COMPUTE X

c = x(mp1)
DO  j = 1,m
  l = mp1 - j
  x(j) = x(j) - c*g(l)
END DO
IF (mp1 == n) RETURN

!     COMPUTE NUMERATOR AND DENOMINATOR OF G(M+1) AND H(M+1)

l = n - mp1
gn = -a(l)
gd = -a(n)
l = n + mp1
hn = -a(l)
DO  j = 1,m
  l = mp1 - j
  nml = n - l
  npl = n + l
  gn = gn + a(nml)*g(j)
  gd = gd + a(nml)*h(l)
  hn = hn + a(npl)*h(j)
END DO
IF (gd == 0.d0) GO TO 100
g(mp1) = gn/gd
h(mp1) = hn/xd

!     COMPUTE G AND H

c1 = g(mp1)
c2 = h(mp1)
MAX = mp1/2
k = m
DO  j = 1,MAX
  gj = g(j)
  gk = g(k)
  hj = h(j)
  hk = h(k)
  g(j) = gj - c1*hk
  g(k) = gk - c1*hj
  h(j) = hj - c2*gk
  h(k) = hk - c2*gj
  k = k - 1
END DO
GO TO 10

!     ERROR RETURN

100 ierr = 1
RETURN
END SUBROUTINE dtoplx
