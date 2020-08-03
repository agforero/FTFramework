FUNCTION polyarea(x, y, nb) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:06

IMPLICIT NONE

REAL, INTENT(IN)     :: x(:)
REAL, INTENT(IN)     :: y(:)
INTEGER, INTENT(IN)  :: nb
REAL                 :: fn_val

!*****************************************************************

!   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
! itself.

!*****************************************************************

INTEGER  :: i, n, nm1
REAL     :: a

n = nb
IF (x(1) == x(n) .AND. y(1) == y(n)) n = n - 1

SELECT CASE (n)
  CASE (:2)
     fn_val = 0.0

  CASE (3)
     fn_val= 0.5*((x(2) - x(1))*(y(3) - y(1)) - (x(3) - x(1))*(y(2) - y(1)))

  CASE DEFAULT
    nm1 = n - 1
    a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
    DO  i = 2, nm1
      a = a + x(i)*(y(i+1) - y(i-1))
    END DO
    fn_val = 0.5*a
END SELECT

RETURN
END FUNCTION polyarea
