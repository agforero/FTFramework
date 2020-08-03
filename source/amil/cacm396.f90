FUNCTION t_quantile(p, ndf, normdev) RESULT(fn_val)

!  Calculates the two-tail quantiles of Student's t-distribution.
!  The user must supply a function `normdev' to return quantiles of the
!  the normal distribution, such as Applied Statistics algorithm AS 241.

!  Translation from Algol by Alan Miller, CSIRO Division of Mathematics
!  & Statistics, Clayton, Victoria 3169, Australia of:

!  Algorithm 396: Student's t-quantiles by G.W. Hill
!  Comm. A.C.M., vol.13(10), 619-620, October 1970

!  N.B. The number of degrees of freedom is not necessarily an integer

!  Arguments:
!  P       double precision    Area under both tails
!  ndf     double precision    Number of degrees of freedom
!  normal  double precision    External function to return quantiles
!                              of the standard normal curve
!  ier     integer             Error indicator = 0 for normal exit
!                                              = 1 if ndf < 1

!  ELF90-compatible version by Alan Miller
!  Latest revision - 7 August 1997

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN) :: p, ndf
REAL (dp)             :: fn_val

INTERFACE
  FUNCTION normdev(p) RESULT(x)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: p
    REAL (dp)             :: x
  END FUNCTION normdev
END INTERFACE

!     Local variables

REAL (dp)            :: a, b, c, d, prob, x, y
REAL (dp), PARAMETER :: half_pi = 1.5707963267949D0, eps = 1.d-12, one = 1.d0

IF (ndf < one .OR. p > one .OR. p <= 0.d0) THEN
  WRITE(*, *) 'Error in t_quantile: ndf < 1 or p not between 0 and 1'
  RETURN
END IF

IF ( ABS(ndf - 2.d0) < eps ) THEN
  fn_val = SQRT(2.d0 / ( p * (2.d0 - p) ) - 2.d0)
  RETURN
  
ELSE IF (ndf < one+eps) THEN
  prob = p * half_pi
  fn_val = COS(prob) / SIN(prob)
  RETURN
  
ELSE
  a = one / (ndf - 0.5D0)
  b = 48.d0 / a**2
  c = ((20700.d0 * a / b - 98.d0) * a - 16.d0) * a + 96.36D0
  d = ((94.5D0 / (b + c) - 3.d0) / b + one) * SQRT(a * half_pi)* ndf
  x = d * p
  y = x ** (2.d0 / ndf)
  
  IF (y > 0.05D0 + a) THEN
    
!     Asymptotic inverse expansion about normal
    
    x = normdev(0.5 * p)
    y = x**2
    IF (ndf < 5.d0) c = c + 0.3D0 * (ndf - 4.5D0) * (x + 0.6D0)
    c = (((0.05D0 * d * x - 5.d0) * x - 7.d0) * x - 2.d0) * x + b+ c
    y = (((((0.4D0*y + 6.3D0) * y + 36.d0) * y + 94.5D0) / c - y  &
         - 3.d0) / b + one) * x
    y = a * y**2
    IF (y > 0.002D0) THEN
      y = EXP(y) - one
    ELSE
      y = 0.5D0 * y**2 + y
    END IF
  ELSE
    
    y = ((one / (((ndf + 6.d0) / (ndf * y) - 0.089D0 * d -  &
        0.822D0) * (ndf + 2.d0) * 3.d0) + 0.5D0 / (ndf + 4.d0))  &
        * y - one) * (ndf + one) / (ndf + 2.d0) + one / y
  END IF
END IF
fn_val = SQRT(ndf * y)

RETURN
END FUNCTION t_quantile
