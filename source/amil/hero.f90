PROGRAM hero

! Nick Stoke's program to find the roots and weights for half-Hermite
! integration.   Converted to quadruple-precision by Alan Miller
! Latest revision - 15 November 1997

USE quadruple_precision
IMPLICIT NONE
TYPE (quad)          :: f(60), g(60), h(3, 3), hi, xx(40), x, x1, x2, xd, &
                        xw(40)
REAL (dp), PARAMETER :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp
INTEGER              :: i, i1, i1j, il1, j, j1, k, l, m, n1

!     OPEN files for output

OPEN(8, FILE='absciss.dat')
OPEN(9, FILE='weights.dat')

!     Set I1 = 60, J1 = no. of iterations = 40, N1 = no. of points.

i1 = 60
j1 = 40
n1 = 20
DO i = 1,60
  f(i) = quad(zero, zero)
  g(i) = quad(zero, zero)
END DO

!     Iterative formation of recurrence relation coefficients.

DO j = 1, j1
  i1j = i1-j
  DO i = 2,i1j
    hi = quad(DBLE(i-1), zero)
    f(i) = hi/DBLE(6) + f(i)*(g(i)*(g(i-1)-g(i)) + half +  &
           f(i)-f(i+1)) / (3.d0*f(i) - 1.5D0*hi)
    g(i-1) = SQRT(hi - half - f(i) - f(i-1))
  END DO
END DO

il1 = i1-1
DO i = 2, il1
  f(i) = SQRT(f(i))
END DO

!     Find the zeroes by Newton's method, evaluating the polynomials
!     by recurrence.

DO j = 1, n1
  x = quad(zero, zero)
  x1 = quad(zero, zero)
  DO k = 1, j
    l = 0
    DO
      h(2,1) = quad(one, zero)
      h(2,2) = quad(zero, zero)
      h(3,2) = quad(zero, zero)
      h(3,1) = quad(zero, zero)
      DO m = 1, j
        DO i = 2, 3
          h(i,3) = h(i,2)
          h(i,2) = h(i,1)
          h(i,1) = (h(i,2)*(g(m)-x) - h(i,3)*f(m) - (i-2)*h(i-1,2))/f(m+1)
        END DO
      END DO

      xd = -h(2,1)/h(3,1)
      x = x+xd
      IF(ABS(xd%hi) > 1.d-12) CYCLE
      l = l+1
      IF(l > 1) EXIT
    END DO
    xw(k) = - quad(half, zero) / (h(2,2)*h(3,1)*g(1)*f(j+1))

!   XW = Weights.

    xx(k) = x
    x2 = x1
    x1 = SQRT(x)
    x = x1 + x1 - x2
    x = x*x
  END DO
  WRITE(8,10) j, xx(1:j)%hi
  WRITE(9,10) j, xw(1:j)%hi
END DO
STOP

10 FORMAT(i3, "  ", 3(e21.15, ", ") / ("     ", 3(e21.15, ", ")) )
END PROGRAM hero
