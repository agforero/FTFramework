FUNCTION sqmcor(x, ip, n, rho2) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-24  Time: 14:51:12

!    ALGORITHM AS 260  APPL. STATIST. (1991) VOL. 40, NO. 1

!    Computes the C.D.F. for the distribution of the square of the
!    multiple correlation coefficient with parameters X, IP, N, and RHO2.
!    X is the sample value of R**2.
!    IP is the number of predictors, including 1 for the constant if one is
!    being fitted.
!    N is the number of cases.
!    RHO2 is the population value of the squared multiple correlation
!    coefficient (often set = 0).

!    The following auxiliary algorithms are required:
!    ALNGAM - Log-gamma function  (ALNGAM  AS 245)
!    BETAIN - Incomplete beta function (AS 63)

! N.B. Argument IFAULT has been removed.

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: ip
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: rho2
REAL (dp)              :: fn_val

REAL (dp)  :: a, b, beta, errbd, gx, q, sumq, temp, term, xj
! REAL :: ALNGAM, betain
! EXTERNAL ALNGAM, betain
INTERFACE
  FUNCTION betain(x, p, q, beta) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 100)
    REAL (dp), INTENT(IN)  :: x, p, q, beta
    REAL (dp)              :: fn_val
  END FUNCTION betain

  FUNCTION alngam(xvalue) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 100)
    REAL (dp), INTENT(IN)  :: xvalue
    REAL (dp)              :: fn_val
  END FUNCTION alngam
END INTERFACE

REAL (dp), PARAMETER  :: errmax = 1.0D-8, zero = 0.0_dp, half = 0.5_dp,   &
                         one = 1.0_dp
INTEGER, PARAMETER    :: itrmax = 100

fn_val = x
IF (rho2 < zero .OR. rho2 > one .OR. ip < 2 .OR. n <= ip) THEN
  WRITE(*, *) 'Error in Function sqmcor, rho2 must be in the range (0,1)'
  RETURN
END IF
IF (x < zero .OR. x > one) THEN
  WRITE(*, *) 'Error in Function sqmcor, x must be in the range (0,1)'
  RETURN
END IF

IF (x == zero .OR. x == one) RETURN

a = half * (ip - 1)
b = half * (n - ip)

!        Initialize the series

beta = EXP(ALNGAM(a) + ALNGAM(b) - ALNGAM(a + b))
temp = betain(x, a, b, beta)

!        There is no need to test IFAULT since all of the
!        parameter values have already been checked

gx = EXP(a * LOG(x) + b * LOG(one - x) - LOG(a)) / beta
q = (one - rho2) ** (a + b)
xj = zero
term = q * temp
sumq = one - q
fn_val = term

!        Perform recurrence until convergence is achieved

10 xj = xj + one
temp = temp - gx
gx = gx * (a + b + xj - one) * x / (a + xj)
q = q * (a + b + xj - one) * rho2 / xj
sumq = sumq - q
term = temp * q
fn_val = fn_val + term

!        Check for convergence and act accordingly

errbd = (temp - gx) * sumq
IF (INT(xj) < itrmax .AND. errbd > errmax) GO TO 10
IF (errbd > errmax) THEN
  WRITE(*, *) 'Error in SQMCOR, errbd > errmax'
END IF
RETURN
END FUNCTION sqmcor



FUNCTION betain(x, p, q, beta) RESULT(fn_val)

!     Algorithm AS 63  Applied Statistics (1973), vol.22, no.3

!     Computes incomplete beta function ratio for arguments
!     x between zero and one, p and q positive.
!     log of complete beta function, beta, is assumed to be known

! ELF90-compatible version by Alan Miller
! N.B. Argument IFAULT has been removed

! Latest revision - 29 November 1997

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: x, p, q, beta
REAL (dp)             :: fn_val

! Local variables
LOGICAL    :: indx
INTEGER    :: ns
REAL (dp)  :: psq, cx, xx, pp, qq, term, ai, rx, temp

!     define accuracy and initialise

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, acu = 1.0E-14_dp

fn_val = x

!     test for admissibility of arguments

IF(p <= zero .OR. q <= zero) THEN
  WRITE(*, *) 'AS63: Either p or q <= 0'
  RETURN
END IF
IF(x < zero .OR. x > one) THEN
  WRITE(*, *) 'AS63: Argument x outside range (0, 1)'
  RETURN
END IF
IF(x == zero .OR. x == one) RETURN

!     change tail if necessary and determine s

psq = p + q
cx = one - x
IF(p < psq*x) THEN
  xx = cx
  cx = x
  pp = q
  qq = p
  indx = .true.
ELSE
  xx = x
  pp = p
  qq = q
  indx = .false.
END IF
term = one
ai = one
fn_val = one
ns = qq + cx*psq

!     Use Soper's reduction formulae.

rx = xx/cx
3 temp = qq - ai
IF(ns == 0) rx = xx
4 term = term*temp*rx/(pp+ai)
fn_val = fn_val + term
temp = ABS(term)
IF(temp <= acu .AND. temp <= acu*fn_val) GO TO 5
ai = ai + one
ns = ns - 1
IF(ns >= 0) GO TO 3
temp = psq
psq = psq + one
GO TO 4

!     calculate result

5 fn_val = fn_val*EXP(pp*LOG(xx) + (qq-one)*LOG(cx) - beta)/pp
IF(indx) fn_val = one - fn_val

RETURN
END FUNCTION betain



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
