SUBROUTINE hybrid(nparam, param, n, t0, rho, nt, bl, bu, np, points, low,  &
                  ifault)

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

!     .. Scalar Arguments ..

INTEGER, INTENT(IN)        :: nparam
REAL (dp), INTENT(OUT)     :: param(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: t0
REAL (dp), INTENT(IN)      :: rho
INTEGER, INTENT(IN)        :: nt
REAL (dp), INTENT(IN OUT)  :: bl(:)
REAL (dp), INTENT(IN)      :: bu(:)
INTEGER, INTENT(OUT)       :: np
REAL (dp), INTENT(OUT)     :: points(:,:)    ! points(nparam,n)
REAL (dp), INTENT(OUT)     :: low(:)
INTEGER, INTENT(OUT)       :: ifault

!     ..
!     .. Local Scalars ..
REAL (dp) :: bst, diff, newobj, objf, p, t, u
INTEGER   :: i,i1,i2
!     ..
!     .. External Functions ..
INTERFACE
  FUNCTION random (bl, bu) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN) :: bl, bu
    REAL (dp)             :: fn_val
  END FUNCTION random
END INTERFACE
!     ..
!     .. External Subroutines ..
INTERFACE
  SUBROUTINE funct1 (np, param, objf)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: np
    REAL (dp), INTENT(IN)  :: param(:)
    REAL (dp), INTENT(OUT) :: objf
  END SUBROUTINE funct1

  SUBROUTINE tradnl (np, param, bl, bu, objf)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)        :: np
    REAL (dp), INTENT(IN OUT)  :: param(:)
    REAL (dp), INTENT(IN)      :: bl(:), bu(:)
    REAL (dp), INTENT(OUT)     :: objf
  END SUBROUTINE tradnl
END INTERFACE

!      Define the constants ZERO and ONE and set
!      initial values for BST and IFAULT

!     .. Data statements ..
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp
!     ..

ifault = 0
bst = 1.0E+10_dp

!      ******************************************
!      *** Start of Annealing component to    ***
!      *** generate a set of starting points. ***
!      ******************************************

!      Initialise the temperature

t = t0

!      Choose an initial point

DO i = 1,nparam
  param(i) = random(bl(i), bu(i))
END DO

!      Calculate the function value for the initial point

CALL funct1(nparam, param, objf)

!      Begin the loop over I1 until Nt Temp. reductions.

DO i1 = 1,nt
  
!      Set the counter for the No of points accepted to one
  
  np = 1
  
!      Begin the loop over I2 until N points have been considered.
  
  DO i2 = 1,n
    
!      Record the present point
    
    DO i = 1,nparam
      points(i,np+1) = param(i)
    END DO
    
!      Select a new point
    
    DO i = 1,nparam
      param(i) = random(bl(i), bu(i))
    END DO
    
!      Calculate the function value of this new point
    
    CALL funct1(nparam, param, newobj)
    
!      Calculate the difference
    
    diff = newobj - objf
    
!      If the new point is better than the old one then increase the
!      Np counter, and record it as a possible starting point. Move
!      to this new point.
    
    IF (diff < zero) THEN
      np = np + 1
      DO i = 1,nparam
        points(i,np) = param(i)
      END DO
      objf = newobj
      IF (newobj < bst) THEN
        bst = newobj
        DO i = 1,nparam
          points(i,1) = param(i)
        END DO
      END IF
      
!      Otherwise use the Metropolis Criterion
      
    ELSE
      
!      P is the comparison statistic. U is a standard Uniform variate.
      
      p = EXP(-diff/t)
      u = random(zero, one)
      
!      If U < P then accept this point.  Increase the Np counter.
!      Record this point as a possible starting point, and move
!      to this new point.
      
      IF (u < p) THEN
        np = np + 1
        objf = newobj
        DO i = 1,nparam
          points(i,np) = param(i)
        END DO
        
!      If U >= P then reject this point and return to the old point.
        
      ELSE
        DO i = 1,nparam
          param(i) = points(i,np+1)
        END DO
      END IF
      
!      Both alternatives for the comparison of the new and old point
!      Have been considered.
      
    END IF
    
!      Continue with the loop until N new points have been considered
    
  END DO
  
!      Set the error indicator to 1 if the algorithm appears to have converged.
  
  IF (np == 1) THEN
    ifault = i1
  END IF
  
!      Once N points have been considered, lower the Temperature.
  
  t = t*rho
  
!      Continue with the loop until the Temp. has been reduced Nt times
  
END DO

!       *****************************************
!       *** End of Annealing component, NP    ***
!       *** starting points are now stored in ***
!       *** the array POINTS.                 ***
!       *****************************************

!       Set an initially high value of LOW(NPARAM+1)

low(nparam+1) = 1.0E+10_dp

!       *******************************************
!       *** Start of the traditional component, ***
!       *** which loops over the NP starting    ***
!       *** points calling the traditional      ***
!       *** minimisation routine for each.      ***
!       *******************************************

!       Begin the loop over the Np starting points.

DO i1 = 1,np
  
!       Set the PARAM array to be the values of the I1st start point
  
  DO i2 = 1,nparam
    param(i2) = points(i2,i1)
  END DO
  
!       Call traditional minimisation routine
  
  CALL tradnl(nparam, param, bl, bu, objf)
  
!       Check to see if this starting point gives the
!       lowest function value to date
  
  IF (objf < low(nparam+1)) THEN
    DO i = 1,nparam
      low(i) = param(i)
    END DO
    low(nparam+1) = objf
  END IF
  
!       End the Np loop
  
END DO

!       **************************************
!       *** End of the second, traditional ***
!       *** component.                     ***
!       **************************************

!       Return to calling routine

RETURN

!       End of subroutine HYBRID

END SUBROUTINE hybrid



FUNCTION random (bl, bu) RESULT(fn_val)

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), INTENT(IN) :: bl, bu
REAL (dp)             :: fn_val

CALL RANDOM_NUMBER(fn_val)
fn_val = bl + (bu - bl)*fn_val

RETURN
END FUNCTION random
