SUBROUTINE random_Zipf(a, r, setup)
! Generate a random deviate (r) from the Zipf (or Zeta) distribution.
! If setup = .TRUE. (or local variable b < 0), two local variables (b and
! const) are set/reset for future calls, otherwise it is assumed that
! the value of a is unchanged since the last call.

! The value of a must be > 1.0 (NOT CHECKED)

! Algorithm from page 551 of:
! Devroye, Luc (1986) `Non-uniform random variate generation',
! Springer-Verlag: Berlin.   ISBN 3-540-96305-7 (also 0-387-96305-7)

! Programmer: Alan Miller (amiller @ users,bigpond.net.au)
! Latest revision - 3 February 1998

IMPLICIT NONE
REAL, INTENT(IN)        :: a
INTEGER, INTENT(OUT)    :: r
LOGICAL, INTENT(IN OUT) :: setup

! Local variables
REAL, SAVE      :: b = -1.0, const
REAL            :: t, u, v
REAL, PARAMETER :: one = 1.0

IF (b < 0.0 .OR. setup) THEN
  b = 2.0**(a - one)
  const = -one/(a - one)
  setup = .FALSE.
END IF

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  r = FLOOR(u**const)
  t = (one + one/r)**(a - one)
  IF (v*r*(t - one)/(b - one) <= t/b) EXIT
END DO

RETURN
END SUBROUTINE random_Zipf



PROGRAM t_random_Zipf
IMPLICIT NONE
INTEGER              :: r, i, nfreq, count, i1
INTEGER, ALLOCATABLE :: freq(:)
REAL, ALLOCATABLE    :: expctd(:)
REAL                 :: a, total
REAL, PARAMETER      :: one = 1.0
LOGICAL              :: setup

INTERFACE
  SUBROUTINE random_Zipf(a, r, setup)
    IMPLICIT NONE
    REAL, INTENT(IN)        :: a
    INTEGER, INTENT(OUT)    :: r
    LOGICAL, INTENT(IN OUT) :: setup
  END SUBROUTINE random_Zipf
END INTERFACE

WRITE(*, *) 'Enter a (> 1.0): '
READ(*, *) a
setup = .TRUE.
! 1/20 = 5%   Store individual frequencies up to last 5%
nfreq = NINT( (20.0/(a - one))**(one/(a-one)) )
ALLOCATE( freq(nfreq), expctd(nfreq) )
freq = 0

! Generate 1000 random deviates from the Zipf distribution

DO i = 1, 1000
  CALL random_Zipf(a, r, setup)
  r = MIN(r, nfreq)
  freq(r) = freq(r) + 1
END DO

! Calculate approximate expected frequencies
total = 0.0
DO i = 1, nfreq-1
  expctd(i) = REAL(i)**(-a)
  total = total + expctd(i)
END DO
expctd(nfreq) = one / ( (a-one)*(nfreq - 0.5)**(a-one) )
total = total + expctd(nfreq)
expctd = expctd * (1000. / total)

WRITE(*, *) '    Range      Obs.freq.  Expctd.freq.'
i = 1
DO
  count = 0
  total = 0
  i1 = i
  DO
    count = count + freq(i)
    total = total + expctd(i)
    IF (total > 25. .OR. i >= nfreq) EXIT
    i = i + 1
  END DO
  IF (i == nfreq) i = 999999
  WRITE(*, '(i6, " - ", i6, i7, "    ", f9.2)') i1, i, count, total
  IF (i == 999999) EXIT
  i = i + 1
END DO

STOP
END PROGRAM t_random_Zipf
