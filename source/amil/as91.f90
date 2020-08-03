FUNCTION ppchi2(p, v, g) RESULT(fn_val)

! N.B. Argument IFAULT has been removed.

! This version by Alan Miller
! amiller @ bigpond.net.au
! Latest revision - 27 October 2000

!  Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35

!  To evaluate the percentage points of the chi-squared
!  probability distribution function.

!  p must lie in the range 0.000002 to 0.999998,
!  v must be positive,
!  g must be supplied and should be equal to ln(gamma(v/2.0))

!  Incorporates the suggested changes in AS R85 (vol.40(1), pp.233-5, 1991)
!  which should eliminate the need for the limited range for p above,
!  though these limits have not been removed from the routine.

!  If IFAULT = 4 is returned, the result is probably as accurate as
!  the machine will allow.

!  Auxiliary routines required: PPND = AS 111 (or AS 241) and GAMMAD = AS 239.

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: p
REAL (dp), INTENT(IN)  :: v
REAL (dp), INTENT(IN)  :: g
REAL (dp)              :: fn_val

INTERFACE
  FUNCTION gammad(x, p) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: x, p
    REAL (dp)             :: fn_val
  END FUNCTION gammad

  SUBROUTINE ppnd16 (p, normal_dev, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: p
    INTEGER, INTENT(OUT)    :: ifault
    REAL (dp), INTENT(OUT)  :: normal_dev
  END SUBROUTINE ppnd16
END INTERFACE

! Local variables

REAL (dp)  :: a, b, c, p1, p2, q, s1, s2, s3, s4, s5, s6, t, x, xx
INTEGER    :: i, if1

INTEGER, PARAMETER    :: maxit = 20
REAL (dp), PARAMETER  :: aa = 0.6931471806_dp, e = 0.5e-06_dp,         &
                         pmin = 0.000002_dp, pmax = 0.999998_dp,       &
                         zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,   &
                         two = 2.0_dp, three = 3.0_dp, six = 6.0_dp,   &
                         c1 = 0.01_dp, c2 = 0.222222_dp, c3 = 0.32_dp, &
                         c4 = 0.4_dp, c5 = 1.24_dp, c6 = 2.2_dp,       &
                         c7 = 4.67_dp, c8 = 6.66_dp, c9 = 6.73_dp,     &
                         c10 = 13.32_dp, c11 = 60.0_dp, c12 = 70.0_dp, &
                         c13 = 84.0_dp, c14 = 105.0_dp, c15 = 120.0_dp, &
                         c16 = 127.0_dp, c17 = 140.0_dp, c18 = 175.0_dp, &
                         c19 = 210.0_dp, c20 = 252.0_dp, c21 = 264.0_dp, &
                         c22 = 294.0_dp, c23 = 346.0_dp, c24 = 420.0_dp, &
                         c25 = 462.0_dp, c26 = 606.0_dp, c27 = 672.0_dp, &
                         c28 = 707.0_dp, c29 = 735.0_dp, c30 = 889.0_dp, &
                         c31 = 932.0_dp, c32 = 966.0_dp, c33 = 1141.0_dp, &
                         c34 = 1182.0_dp, c35 = 1278.0_dp, c36 = 1740.0_dp, &
                         c37 = 2520.0_dp, c38 = 5040.0_dp

!       Test arguments and initialise

fn_val = -one
IF (p < pmin .OR. p > pmax) THEN
  WRITE(*, *) 'Error in PPCHI2: p must be between 0.000002 & 0.999998'
  RETURN
END IF
IF (v <= zero) THEN
  WRITE(*, *) 'Error in PPCHI2: Number of deg. of freedom <= 0'
  RETURN
END IF

xx = half * v
c = xx - one

!       Starting approximation for small chi-squared

IF (v < -c5 * LOG(p)) THEN
  fn_val = (p * xx * EXP(g + xx * aa)) ** (one/xx)
  IF (fn_val < e) GO TO 6
  GO TO 4
END IF

!       Starting approximation for v less than or equal to 0.32

IF (v > c3) GO TO 3
fn_val = c4
a = LOG(one-p)

2 q = fn_val
p1 = one + fn_val * (c7+fn_val)
p2 = fn_val * (c9 + fn_val * (c8 + fn_val))
t = -half + (c7 + two * fn_val) / p1 - (c9 + fn_val * (c10 + three * fn_val)) / p2
fn_val = fn_val - (one - EXP(a + g + half * fn_val + c * aa) * p2 / p1) / t
IF (ABS(q / fn_val - one) > c1) GO TO 2
GO TO 4

!       Call to algorithm AS 241 - note that p has been tested above.

3 CALL ppnd16(p, x, if1)

!       Starting approximation using Wilson and Hilferty estimate

p1 = c2 / v
fn_val = v * (x * SQRT(p1) + one - p1) ** 3

!       Starting approximation for p tending to 1

IF (fn_val > c6 * v + six) fn_val = -two * (LOG(one-p) - c * LOG(half * fn_val) + g)

!       Call to algorithm AS 239 and calculation of seven term Taylor series

4 DO i = 1, maxit
  q = fn_val
  p1 = half * fn_val
  p2 = p - gammad(p1, xx)

  t = p2 * EXP(xx * aa + g + p1 - c * LOG(fn_val))
  b = t / fn_val
  a = half * t - b * c
  s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 + c11 * a))))) / c24
  s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 * a)))) / c37
  s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37
  s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 + c36 * a))) / c38
  s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37
  s6 = (c15 + c * (c23 + c16 * c)) / c38
  fn_val = fn_val + t * (one + half * t * s1 - b * c * (s1 - b *   &
           (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))
  IF (ABS(q / fn_val - one) > e) RETURN
END DO

WRITE(*, *) 'Error in PPCHI2: Max. number of iterations exceeded'

6 RETURN
END FUNCTION ppchi2
