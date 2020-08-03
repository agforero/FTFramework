PROGRAM test_kemp_binomial
IMPLICIT NONE

INTEGER :: n, nrep, i, ix
REAL    :: p, mean = 0.0, sumsq = 0.0, stdev = 0.0, dev
LOGICAL :: first

INTERFACE
  FUNCTION random_binomial(n, p, first) RESULT(ran_binom)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN)    :: p
    LOGICAL, INTENT(IN) :: first
    INTEGER             :: ran_binom
  END FUNCTION random_binomial
END INTERFACE

WRITE(*, *) 'Enter n & p: '
READ(*, *) n, p
WRITE(*, *) 'How many random deviates?: '
READ(*, *) nrep

first = .TRUE.
DO i = 1, nrep
  ix = random_binomial(n, p, first)
  dev = ix - mean
  mean = mean + dev/i
  sumsq = sumsq + dev*(ix - mean)
  first = .FALSE.
END DO

IF (nrep > 1) stdev = SQRT(sumsq / (nrep-1))
WRITE(*, '(1x, "Mean = ", f8.2, 5x, "Std. devn. = ", f8.2)') mean, stdev

STOP
END PROGRAM test_kemp_binomial



FUNCTION random_binomial(n, p, first) RESULT(ran_binom)

! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN)    :: p
LOGICAL, INTENT(IN) :: first
INTEGER             :: ran_binom

INTERFACE
  FUNCTION bin_prob(n, p, r) RESULT(prob)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, r
    REAL, INTENT(IN)    :: p
    REAL                :: prob
  END FUNCTION bin_prob
END INTERFACE

!     Local variables

INTEGER             :: ru, rd
INTEGER, SAVE       :: r0
REAL                :: u, pd, pu, zero = 0.0, one = 1.0
REAL, SAVE          :: odds_ratio, p_r

IF (first) THEN
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
END IF

CALL RANDOM_NUMBER(u)
u = u - p_r
IF (u < zero) THEN
  ran_binom = r0
  RETURN
END IF

pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * REAL(n-rd))
    u = u - pd
    IF (u < zero) THEN
      ran_binom = rd
      RETURN
    END IF
  END IF

  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / REAL(ru)
    u = u - pu
    IF (u < zero) THEN
      ran_binom = ru
      RETURN
    END IF
  END IF
END DO

!     This point should not be reached, but just in case:

ran_binom = r0
RETURN

END FUNCTION random_binomial



FUNCTION bin_prob(n, p, r) RESULT(prob)
!     Calculate a binomial probability

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, r
REAL, INTENT(IN)    :: p
REAL                :: prob

!     Local variable
REAL                :: one = 1.0

prob = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
                + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN


CONTAINS


FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

IMPLICIT NONE
INTEGER, PARAMETER          :: doub_prec = SELECTED_REAL_KIND(15, 60)
REAL(doub_prec), INTENT(IN) :: z
REAL(doub_prec)             :: lanczos

! Local variables

REAL(doub_prec)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                              -1259.139216722289D0, 771.3234287757674D0, &
                              -176.6150291498386D0, 12.50734324009056D0, &
                              -0.1385710331296526D0, 0.9934937113930748D-05, &
                               0.1659470187408462D-06 /), zero = 0.D0,   &
                               one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
                               half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER          :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma

END FUNCTION bin_prob




