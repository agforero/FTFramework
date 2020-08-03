!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 2001

FUNCTION alnorm( x, upper ) RESULT( fn_val )
   IMPLICIT NONE
   INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 100)
   REAL(DP), INTENT(IN)   ::  x
   LOGICAL,   INTENT(IN)  ::  upper
   REAL(DP)               ::  fn_val

   !  Local variables
   REAL(DP), PARAMETER   ::  zero=0.0_DP, one=1.0_DP, half=0.5_DP, con=1.28_DP
   REAL(DP)              ::  z, y
   LOGICAL               ::  up

   !  Machine dependent constants
   REAL(DP), PARAMETER  ::  ltone = 7.0_DP, utzero = 18.66_DP
   REAL(DP), PARAMETER  ::  p = 0.398942280444_DP, q = 0.39990348504_DP,   &
                            r = 0.398942280385_DP, a1 = 5.75885480458_DP,  &
                            a2 = 2.62433121679_DP, a3 = 5.92885724438_DP,  &
                            b1 = -29.8213557807_DP, b2 = 48.6959930692_DP, &
                            c1 = -3.8052E-8_DP, c2 = 3.98064794E-4_DP,     &
                            c3 = -0.151679116635_DP, c4 = 4.8385912808_DP, &
                            c5 = 0.742380924027_DP, c6 = 3.99019417011_DP, &
                            d1 = 1.00000615302_DP, d2 = 1.98615381364_DP,  &
                            d3 = 5.29330324926_DP, d4 = -15.1508972451_DP, &
                            d5 = 30.789933034_DP

   up = upper
   z = x
   IF( z < zero ) THEN
      up = .NOT. up
      z = -z
   END IF
   IF( z <= ltone  .OR.  (up  .AND.  z <= utzero) ) THEN
      y = half*z*z
      IF( z > con ) THEN
         fn_val = r*EXP( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
      ELSE
         fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
      END IF
   ELSE
      fn_val = zero
   END IF

   IF( .NOT. up ) fn_val = one - fn_val
   RETURN
END FUNCTION alnorm
