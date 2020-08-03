FUNCTION zeroin(f, ax, bx, aerr, rerr) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-14  Time: 12:32:54
 
!-----------------------------------------------------------------------

!         FINDING A ZERO OF THE FUNCTION F(X) IN THE INTERVAL (AX,BX)

!                       ------------------------

!  INPUT...

!  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN THE
!         CLOSED INTERVAL (AX,BX).  IT IS ASSUMED THAT F IS CONTINUOUS,
!         AND THAT F(AX) AND F(BX) HAVE DIFFERENT SIGNS.
!  AX     LEFT ENDPOINT OF THE INTERVAL
!  BX     RIGHT ENDPOINT OF THE INTERVAL
!  AERR   THE ABSOLUTE ERROR TOLERANCE TO BE SATISFIED
!  RERR   THE RELATIVE ERROR TOLERANCE TO BE SATISFIED

!  OUTPUT...

!         ABCISSA APPROXIMATING A ZERO OF F IN THE INTERVAL (AX,BX)

!-----------------------------------------------------------------------
!  ZEROIN IS A SLIGHTLY MODIFIED TRANSLATION OF THE ALGOL PROCEDURE
!  ZERO GIVEN BY RICHARD BRENT IN ALGORITHMS FOR MINIMIZATION WITHOUT
!  DERIVATIVES, PRENTICE-HALL, INC. (1973).
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

REAL (dp), INTENT(IN)  :: ax
REAL (dp), INTENT(IN)  :: bx
REAL (dp), INTENT(IN)  :: aerr
REAL (dp), INTENT(IN)  :: rerr
REAL (dp)              :: fn_val

! EXTERNAL f
INTERFACE
  FUNCTION f(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp)  :: a, b, c, d, e, eps, fa, fb, fc, tol, xm, p, q, r, s, atol, rtol

!  COMPUTE EPS, THE RELATIVE MACHINE PRECISION

eps = EPSILON(0.0_dp)

! INITIALIZATION

a = ax
b = bx
fa = f(a)
fb = f(b)
atol = 0.5 * aerr
rtol = MAX(0.5_dp*rerr, 2.0_dp*eps)

! BEGIN STEP

10 c = a
fc = fa
d = b - a
e = d
20 IF (ABS(fc) < ABS(fb)) THEN
  a = b
  b = c
  c = a
  fa = fb
  fb = fc
  fc = fa
END IF

! CONVERGENCE TEST

tol = rtol * MAX(ABS(b),ABS(c)) + atol
xm = 0.5 * (c-b)
IF (ABS(xm) > tol) THEN
  IF (fb /= 0.0) THEN
    
! IS BISECTION NECESSARY
    
    IF (ABS(e) >= tol) THEN
      IF (ABS(fa) > ABS(fb)) THEN
        
! IS QUADRATIC INTERPOLATION POSSIBLE
        
        IF (a == c) THEN
          
! LINEAR INTERPOLATION
          
          s = fb / fc
          p = (c-b) * s
          q = 1.0 - s
        ELSE
          
! INVERSE QUADRATIC INTERPOLATION
          
          q = fa / fc
          r = fb / fc
          s = fb / fa
          p = s * ((c-b)*q*(q-r)-(b-a)*(r-1.0))
          q = (q-1.0) * (r-1.0) * (s-1.0)
        END IF
        
! ADJUST SIGNS
        
        IF (p > 0.0) q = -q
        p = ABS(p)
        
! IS INTERPOLATION ACCEPTABLE
        
        IF (2.0*p < (3.0*xm*q-ABS(tol*q))) THEN
          IF (p < ABS(0.5*e*q)) THEN
            e = d
            d = p / q
            GO TO 30
          END IF
        END IF
      END IF
    END IF
    
! BISECTION
    
    d = xm
    e = d
    
! COMPLETE STEP
    
    30 a = b
    fa = fb
    IF (ABS(d) > tol) b = b + d
    IF (ABS(d) <= tol) b = b + SIGN(tol,xm)
    fb = f(b)
    IF (fb*(fc/ABS(fc)) > 0.0) GO TO 10
    GO TO 20
  END IF
END IF

! DONE

fn_val = b
RETURN
END FUNCTION zeroin
