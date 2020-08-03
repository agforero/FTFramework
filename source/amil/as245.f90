FUNCTION alngam(xvalue) RESULT(fn_val)

!     ALGORITHM AS245  APPL. STATIST. (1989) VOL. 38, NO. 2

!     Calculation of the logarithm of the gamma function

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: xvalue
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x, x1, x2, y

!     Coefficients of rational functions

REAL (dp), PARAMETER :: r1(9) = (/ -2.66685511495_dp, -24.4387534237_dp,  &
                                   -21.9698958928_dp,  11.1667541262_dp,  &
                                    3.13060547623_dp,  0.607771387771_dp, &
                                    11.9400905721_dp,  31.4690115749_dp,  &
                                    15.2346874070_dp /)
REAL (dp), PARAMETER :: r2(9) = (/ -78.3359299449_dp, -142.046296688_dp,  &
                                    137.519416416_dp,  78.6994924154_dp,  &
                                    4.16438922228_dp,  47.0668766060_dp,  &
                                    313.399215894_dp,  263.505074721_dp,  &
                                    43.3400022514_dp /)
REAL (dp), PARAMETER :: r3(9) = (/ -2.12159572323E5_dp,  2.30661510616E5_dp,  &
                                    2.74647644705E4_dp, -4.02621119975E4_dp,  &
                                   -2.29660729780E3_dp, -1.16328495004E5_dp,  &
                                   -1.46025937511E5_dp, -2.42357409629E4_dp,  &
                                   -5.70691009324E2_dp /)
REAL (dp), PARAMETER :: r4(5) = (/ 0.279195317918525_dp, 0.4917317610505968_dp, &
                                   0.0692910599291889_dp, 3.350343815022304_dp, &
                                   6.012459259764103_dp /)

!     Fixed constants

REAL (dp), PARAMETER :: alr2pi = 0.918938533204673_dp, four = 4._dp,  &
                        half = 0.5_dp, one = 1._dp, onep5 = 1.5_dp,   &
                        twelve = 12._dp, zero = 0._dp

!     Machine-dependant constants.
!     A table of values is given at the top of page 399 of the paper.
!     These values are for the IEEE double-precision format for which
!     B = 2, t = 53 and U = 1023 in the notation of the paper.

REAL (dp), PARAMETER :: xlge = 5.10E6_dp, xlgst = HUGE(1.0_dp)

x = xvalue
fn_val = zero

!     Test for valid function argument

IF (x >= xlgst) THEN
  WRITE(*, *) 'AS 245: Argument x too large'
  RETURN
END IF
IF (x <= zero) THEN
  WRITE(*, *) 'AS 245: Argument x <= 0'
  RETURN
END IF

!     Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined

IF (x < onep5) THEN
  IF (x < half) THEN
    fn_val = -LOG(x)
    y = x + one
    
!     Test whether X < machine epsilon
    
    IF (y == one) RETURN
  ELSE
    fn_val = zero
    y = x
    x = (x - half) - half
  END IF
  fn_val = fn_val + x * ((((r1(5)*y + r1(4))*y + r1(3))*y + r1(2))*y + r1(1)) / &
                    ((((y + r1(9))*y + r1(8))*y+ r1(7))*y + r1(6))
  RETURN
END IF

!     Calculation for 1.5 <= X < 4.0

IF (x < four) THEN
  y = (x - one) - one
  fn_val = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x + r2(1)) /  &
               ((((x + r2(9))*x + r2(8))*x + r2(7))*x+ r2(6))
  RETURN
END IF

!     Calculation for 4.0 <= X < 12.0

IF (x < twelve) THEN
  fn_val = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) /  &
           ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
  RETURN
END IF

!     Calculation for X >= 12.0

y = LOG(x)
fn_val = x * (y - one) - half * y + alr2pi
IF (x > xlge) RETURN
x1 = one / x
x2 = x1 * x1
fn_val = fn_val + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) /  &
         ((x2 + r4(5))*x2 + r4(4))
RETURN
END FUNCTION alngam
