FUNCTION q_lngm(x) RESULT(b)

! Extended arithmetic calculation of the logarithm of the gamma
! function (N.B. gamma(x) = (x-1)! in factorial notation):
!  b = ln (gamma (x))
! where all quantities are in quadruple-precision.
! The result (b) may occupy the same location as the input value (x).
! x must not equal a negative integer; in this case a warning is given
! and no result is calculated.

! Algorithm:  The true Stirling's approximation (not the version derived by
! de Moivre which is usually quoted in textbooks) is used for x >= 16.
! For smaller arguments, a Lanczos-type approximation is used.

! Programmer: Alan Miller
!
! Latest Fortran 77 revision - 11 May 1988
! Fortran 90 version - 21 August 1997

USE quadruple_precision
IMPLICIT NONE
TYPE (quad), INTENT(IN) :: x
TYPE (quad)             :: b

! Local variables

LOGICAL            :: large
INTEGER            :: i, j
TYPE (quad)        :: z, total, temp, zz, term
REAL (dp)          :: zero = 0._dp, half = 0.5_dp, one = 1.0_dp
INTEGER, PARAMETER :: g = 14

! Table of values of
!   a(2j) = B(2j).[2**(2j-1) - 1]/[2j.(2j-1).2**(2j-1)]
! where the B(2j)'s are the Bernouilli numbers.
! Stirling's approximation for log(gamma(x)) is then
!   log(gamma(x)) = 0.5*log(2.pi) + z.log(z) - z - sum[a(2j)/z**(2j-1)]
! where z = x - 0.5.

TYPE (quad), PARAMETER :: st_coeff(16) = (/  &
    quad( 0.4166666666666667D-01, -.4625929269271486D-17),  &
    quad( -.2430555555555556D-02, 0.4722302795714642D-18),  &
    quad( 0.7688492063492064D-03, -.8742455613057720D-19),  &
    quad( -.5905877976190476D-03, -.3820521941139396D-19),  &
    quad( 0.8401067971380472D-03, -.2482348408384320D-19),  &
    quad( -.1916590625086719D-02, 0.2283600073631585D-18),  &
    quad( 0.6409473908253206D-02, -.7561615151693776D-18),  &
    quad( -.2954975178039152D-01, 0.2751974392739155D-17),  &
    quad( 0.1796430017910384D+00, 0.5451863036096252D-16),  &
    quad( -.1392429561052216D+01, 0.1922261747261675D-16),  &
    quad( 0.1340285765318479D+02, -.1330122271048497D-14),  &
    quad( -.1568482659282295D+03, 0.5486595470093716D-13),  &
    quad( 0.2193103267973761D+04, -.1576457483073076D-12),  &
    quad( -.3610877098469368D+05, -.1128718158624283D-11),  &
    quad( 0.6914722675633456D+06, -.9528010535925612D-11),  &
    quad( -.1523822153940742D+08, 0.4711160925202246D-08) /)

!     Lanczos-type coeffs. for g = 14

TYPE (quad), PARAMETER :: a(16) = (/  &
    quad( 0.1000000000000000D+01, 0.4846757479531194D-21),  &
    quad( 0.1069184011592089D+07, -.5763228334420908D-09),  &
    quad( -.4731034435392098D+07, 0.1831263028328363D-09),  &
    quad( 0.8831027007071320D+07, -.3562876555263793D-09),  &
    quad( -.9057605936125744D+07, -.8123034251454294D-09),  &
    quad( 0.5575099541767454D+07, 0.4521232468551490D-09),  &
    quad( -.2113664276594438D+07, 0.2261448706472902D-09),  &
    quad( 0.4882973818694724D+06, -.2584262485630417D-10),  &
    quad( -.6580271859490060D+05, 0.6829382319982582D-11),  &
    quad( 0.4754359892166024D+04, -.6321699947095538D-12),  &
    quad( -.1588514130362740D+03, 0.1630647068074753D-13),  &
    quad( 0.1878850329497677D+01, -.5694435459399838D-16),  &
    quad( -.4589872822043718D-02, 0.3894567587921014D-18),  &
    quad( 0.5931363621107230D-06, -.1388883823219028D-21),  &
    quad( 0.1024789274886151D-11, 0.2637600537816092D-27),  &
    quad( -.3355721053959592D-12, 0.4593372363497046D-29) /)

! Check negative arguments.

IF (x%hi <= zero) THEN
  WRITE(*, *)' *** Negative argument in q_lngm ***'
  RETURN
END IF

IF (x%hi == one .OR. x%hi == 2.d0) THEN
  IF (x%lo == zero) THEN
    b%hi = zero
    b%lo = zero
    RETURN
  END IF
END IF

IF (x%hi >= 16._dp) GO TO 20

!  Fit a Lanczos-type approximation if x%hi < 16.
!  i.e. ln(gamma(x)) = lnsqrt(2.pi) + (x-.5)ln(x+g-.5) - (x+g-.5) + ln(A(x))
!  where A(x) = a0 + a1/x + a2/(x+1) + ... + ak/(x+g).
!  g = 14 here to achieve about 29 significant digits accuracy, except in the
!  vicinity of 1 and 2, where ln(gamma) = 0, where the absolute accuracy is
!  about 30 decimal digits.

z = a(1)
DO i = 2, g + 2
  z = z + a(i) / (x + DBLE(i-2))
END DO
temp = x - half
zz = temp + DBLE(g)
b = lnsqrt2pi + temp*log(zz) - zz + log(z)

RETURN

! Stirling's approximation.

20 total = lnsqrt2pi
z = x - half
total = total + z * log(z) - z
zz = z * z
term = z
large = .true.
DO j = 1, 16
  IF (large) THEN
    temp = st_coeff(j) / term
    total = total - temp
    term = term * zz
    IF (ABS(temp%hi) < ABS(total%lo)) large = .false.
  ELSE
    temp%lo = st_coeff(j)%hi / term%hi
    total%lo = total%lo - temp%lo
    IF (ABS(temp%lo) < 1.d-30) CYCLE
    term%hi = term%hi * zz%hi
  END IF
END DO

! The smallest bit of sum%hi & the largest bit of sum%lo may be out
! of alignment, so add zero to re-align.

b = total + zero

RETURN
END FUNCTION q_lngm



PROGRAM test_q_lngm
! Test the quadruple-precision lngamma function using
!    lngamma(m+x) - lngamma(n+x) = log[ (m+x-1).(m+x-2) .. (n+x) ]

USE quadruple_precision
IMPLICIT NONE

INTEGER                :: m, n, i
TYPE (quad)            :: x, lngmm, lngmn, diff, log_prod, error, temp
TYPE (quad), PARAMETER :: qone = quad( 1.0_dp, 0.0_dp)

INTERFACE
  FUNCTION q_lngm(x) RESULT(b)
    USE quadruple_precision
    IMPLICIT NONE
    TYPE (quad), INTENT(IN) :: x
    TYPE (quad)             :: b
  END FUNCTION q_lngm
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter 2 different positive integers: '
  READ(*, *) m, n
  IF (m == n) THEN
    WRITE(*, *) 'The integers must be different - try again!'
    CYCLE
  END IF
  IF (n > m) THEN
    i = m
    m = n
    n = i
  END IF
  CALL RANDOM_NUMBER(x%hi)
  x%lo = 0._dp
  lngmm = q_lngm(x+DBLE(m))
  WRITE(*, '(" lngamma(", f20.15, ") = ", f22.16)') x%hi+DBLE(m), lngmm%hi
  lngmn = q_lngm(x+DBLE(n))
  WRITE(*, '(" lngamma(", f20.15, ") = ", f22.16)') x%hi+DBLE(n), lngmn%hi
  diff = lngmm - lngmn
  temp = x + DBLE(m-1)
  log_prod = log(temp)
  DO i = 1, m-n-1
    temp = temp - qone
    log_prod = log_prod + log(temp)
  END DO
  error = diff - log_prod
  WRITE(*, '(" Diff =", f20.15, "  log_prod =", f20.15, "  Error =", g11.3)')  &
                diff%hi, log_prod%hi, error%hi
END DO

STOP
END PROGRAM test_q_lngm
