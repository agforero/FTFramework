MODULE Gauss_Weights

!      ALGORITHM 125, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN COMMUNICATIONS OF THE ACM,
!      VOL. 5, NO. 10, October, 1962, P.  510--511

! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-03  Time: 14:03:32

! Original Fortran translation donated by

! M.Dow@anu.edu.au
! ANUSF,  Australian National University
! Canberra, Australia

! Tidied up to use workspace arrays, dimension arrays to
! required lengths, use 0: for dimensions of A, add comments,
! and NAG tools to layout source

! trh (20/07/97)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

CONTAINS


SUBROUTINE weightc(n, q, e, eps, w, x)

! Computes the abscissae x(i) and the weight coefficients w(i) for a
! Gaussian quadrature method
!      \int_0^b w(x)f(x) dx \approx \sum_{i=1}^{n}w_if(x_i) where
! \int_0^b w(x) dx = 1 and w(x)>= 0. The method requires the order n, a
! tolerance eps and the 2n-1 first coefficients of the continued fraction
!      \int_0^b {w(x) \over z-x} = { 1| \over |z} - {q_1 | \over |1} -
!      {e_1 | \over |z} - {q_2 | \over |1} -
!      {e_2 | \over |z} - \cdots
! to be given, the latter in the two arrays q(n) and e(n-1) all components
! of which are automatically positive by virtue of the condition w(x)>= 0.
! The method works as well if the upper bound b is actually infinity
! (note that b does not appear directly as an argument!) or if the density
! function w(x) dx is replaced by da(x) with a monotonically increasing a(x)
! with at least n points of of variation. The tolerance eps should be given
! in accordance to the machines accuracy (preferably by using the value of
! d1mach(4)).  The result is delivered as two arrays w(n) (the weight
! coefficients) and x(n) (the abscissae).  For a description of the method see
! H Rutishauser, ``On a modification of the QD-algorithm with Graeffe-type
! convergence'' [Proceedings of the IFIPS Congress, Munich, 1962].

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: q(n)
REAL (dp), INTENT(IN)      :: e(n-1)
REAL (dp), INTENT(IN OUT)  :: eps
REAL (dp), INTENT(OUT)     :: w(n)
REAL (dp), INTENT(OUT)     :: x(n)

! Allocate workspaces G & A
REAL (dp)  :: g(n), a(0:n,0:7)
!     ..
!     .. Local Scalars ..
REAL (dp)  :: m, p
INTEGER    :: k
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC ABS,EXP,LOG
!     ..

x(1) = q(1) + e(1)
DO k = 2,n
  g(k-1) = e(k-1)*q(k)/x(k-1)
  IF (k == n) THEN
    x(k) = q(k) - g(k-1)
    
  ELSE
    x(k) = q(k) + e(k) - g(k-1)
  END IF
  
  g(k-1) = g(k-1)/x(k)
  w(k-1) = x(k)/x(k-1)
  x(k-1) = LOG(x(k-1))
END DO
x(n) = LOG(x(n))
p = 1
30  DO k = 1,n - 1
  IF (ABS(g(k)*w(k)) > eps) GO TO 40
END DO
GO TO 50

40 CALL qdgraeffe(n, x, g, w, a)
p = 2*p
GO TO 30

! What follows is a peculiar method to compute the w(k) from
! the given ratios g_k = w_{k+1}/w_k such that
!  \sum_{k=1}^n w_k = 1, but the straightforward formulae to do
! this might well produce overflow of exponent

50 w(1) = 1
m = 0
DO k = 1,n - 1
  w(k+1) = w(k)*g(k)
  IF (w(k) > m) m = w(k)
END DO
! WRITE (6,FMT=*) 'm=',m
m = SUM( w(1:n) )
DO k = 1,n
  w(k) = w(k)/m
  x(k) = EXP(x(k)/p)
END DO
RETURN

END SUBROUTINE weightc



SUBROUTINE red(a, f, n)

!  Subroutine RED reduces a heptadiagonal matrix a to tridiagonal form as
!  described in the paper referenced above.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(0:n,0:7)
REAL (dp), INTENT(IN)      :: f(n)

!     ..
REAL (dp)  :: c
INTEGER    :: j, k
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC ABS
!     ..
DO k = 1,n - 1
  DO j = k,n - 1
    
    c = -f(j)*a(j,7)/a(j,2)
    a(j,7) = 0
    a(j+1,2) = a(j+1,2) + c*a(j,5)
    a(j,1) = a(j,1) - c*f(j)*a(j,4)
    a(j,6) = a(j,6) - c*a(j+1,1)
    a(j+1,3) = a(j+1,3) - c*a(j+1,6)
  END DO
  DO j = k,n - 1
    c = -f(j)*a(j,4)/a(j,1)
    a(j,4) = 0
    a(j+1,1) = a(j+1,1) + c*a(j,6)
    a(j+1,6) = a(j+1,6) + c*a(j+1,3)
    a(j,5) = a(j,5) - c*a(j+1,2)
    a(j+1,0) = a(j+1,0) - c*a(j+1,5)
  END DO
  DO j = k + 1,n - 1
    c = -a(j,3)/a(j-1,6)
    a(j,3) = 0
    a(j,6) = a(j,6) + c*a(j,1)
    a(j-1,5) = a(j-1,5) - c*f(j)*f(j)*a(j,0)
    a(j,2) = a(j,2) - c*f(j)*f(j)*a(j,5)
    a(j,7) = a(j,7) - c*f(j)*a(j+1,2)
  END DO
  DO j = k + 1,n - 1
    c = -a(j,0)/a(j-1,5)
    a(j,0) = 0
    a(j+1,2) = a(j+1,2) + c*f(j)*a(j,7)
    a(j,5) = a(j,5) + c*a(j,2)
    a(j,1) = a(j,1) - c*f(j)*f(j)*a(j,6)
    a(j,4) = a(j,4) - c*f(j)*a(j+1,1)
  END DO
END DO
RETURN

END SUBROUTINE red



SUBROUTINE qdgraeffe(n, h, g, f, a)

! Subroutine QDGRAEFFE computes for a given continued fraction
!     f(z) = { 1| \over |z} - {q_1 | \over |1} -
!      {e_1 | \over |z} - {q_2 | \over |1} -
!      {e_2 | \over |z} - \cdots - {q_n | \over |1}
! another one, the poles of which are the squares of the poles of f(z)
! However QDGRAEFFE uses not the coefficients q_1 ... q_n and
! e_1 ... e_{n-1} of f(z) but the quotients f_k = q_{k+1}/q_k and
! g_k = e_k/q_{k+1} for k=1,2,...,n-1 and the h_k = ln(abs(q_k)) for
! k=1,2,...,n, and the results are delivered in the same form.  Routine
! QDGRAEFFE can be used independently, but requires subroutine RED

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: h(n)
REAL (dp), INTENT(OUT)  :: g(n)
REAL (dp), INTENT(OUT)  :: f(n)
REAL (dp), INTENT(OUT)  :: a(0:n,0:7)

!     ..
!     .. Local Scalars ..
INTEGER  :: k
!     ..
!     .. External Subroutines ..
! EXTERNAL red
!     ..

g(n) = 0
f(n) = 0
DO k = 1,n
  a(k-1,4) = 1.0_dp
  a(k-1,5) = 1.0_dp
  a(k,2) = 1 + g(k)*f(k)
  a(k,1) = 1 + g(k)*f(k)
  a(k,6) = g(k)
  a(k,7) = g(k)
  a(k,0) = 0.0_dp
  a(k,3) = 0.0_dp
END DO
a(n,5) = 0.0_dp

! The array a represents the heptadiagonal matrix Q of the paper
! cited above, but with the modifications needed to avoid large numbers
! and with a peculiar arrangement.

CALL red(a, f, n)
DO k = 1,n
  h(k) = 2*h(k) + LOG(ABS(a(k,1)*a(k,2)))
END DO
DO k = 1,n - 1
  f(k) = f(k)*f(k)*a(k+1,2)*a(k+1,1)/ (a(k,1)*a(k,2))
  g(k) = a(k,5)*a(k,6)/ (a(k+1,1)*a(k+1,2))
END DO
RETURN

END SUBROUTINE qdgraeffe

END MODULE Gauss_Weights



PROGRAM main
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-03  Time: 14:01:08
 
! Driver for  CACM Alg 125 Rutishauser
! Author M.Dow@anu.edu.au
!        ANUSF,  Australian National University
!        Canberra, Australia

! Tidied up to use workspace arrays, real parameter values
! and put through nag tools

! trh (20/07/97)

USE Gauss_Weights
IMPLICIT NONE

!     ..
!     .. Parameters ..
INTEGER, PARAMETER    :: nm = 100
REAL (dp), PARAMETER  :: zero=0.0D0, one=1.0D0, two=2.0D0, three=3.0D0
!     .. Local Scalars ..
REAL (dp)  :: a, b, eps, exact, s
INTEGER    :: i, n, p
!     ..
!     .. Local Arrays ..
REAL (dp) :: e(nm), q(nm), w(nm), x(nm)
!     ..
!     .. External Subroutines ..
! EXTERNAL weightcoeff
!     ..

10 WRITE (6, FMT=*) ' Enter n, eps: '
READ (5, FMT=*, END=20) n, eps
IF (eps <= zero) eps = 1.D-15
IF (n >= nm-1) STOP
WRITE (6, FMT=*) 'n=', n, ' eps=', eps

!  q,e for w=1/2 interval (0,2)
! Ref: W.Jones & W.Thron Continued Fractions ...
!      Encyclopedia of maths...vol 11 p24
!      Continued Fraction expansion of log(1+x),
! transformed to Rutishauser form p33

DO i = 2,nm
  q(i) = two*i*i/ (two*i* (two*i-one))
  e(i) = two*i*i/ (two*i* (two*i+one))
END DO
q(1) = one
e(1) = one/three
CALL weightc(n, q, e, eps, w, x)

! Adjust weights, zeros for w=1, interval (-1,1)
DO i = 1, n
  w(i) = two*w(i)
  x(i) = x(i) - one
END DO

DO i = 1, MIN(n,10)
  WRITE (6, FMT=9000) w(i), x(i)
END DO

! Check for x^4
p = 4
a = -one
b = one
exact = (b** (p+1)-a** (p+1))/ (p+1)
s = zero
DO i = 1, n
  s = s + w(i)*x(i)**p
END DO
WRITE (6, FMT=9010) p, exact, s, exact - s
GO TO 10

20 STOP

9000 FORMAT (f20.16, '  ', f20.16)
9010 FORMAT (i3, ' Exact=', f20.16, ' Quadrature=', f20.16, ' Error=', e14.7)
END PROGRAM main
