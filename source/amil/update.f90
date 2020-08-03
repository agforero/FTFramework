! This is a collection of 3 routines to calculate the new mean
! and sum of squares of deviations about the mean (and hence the new
! sample variance or standard deviation) after:

! 1. A new observation becomes available (routine update).
! 2. An observation is dropped (routine downdate).
! 3. An observation is replaced with a new one (routine replace).

! Latest revision - 24 November 2001
! Alan Miller (amiller @ bigpond.net.au)


SUBROUTINE update(n, mean, sumsq, x)
! x contains the value of the new case.

IMPLICIT NONE

INTEGER, INTENT(IN OUT)  :: n
REAL, INTENT(IN OUT)     :: mean, sumsq
REAL, INTENT(IN)         :: x

! Local variable
REAL  :: dev

n = n + 1
dev = x - mean
mean = mean + dev/n
sumsq = sumsq + dev*(x - mean)

RETURN
END SUBROUTINE update



SUBROUTINE downdate(n, mean, sumsq, x)
! x contained the value to be removed.

IMPLICIT NONE

INTEGER, INTENT(IN OUT)  :: n
REAL, INTENT(IN OUT)     :: mean, sumsq
REAL, INTENT(IN)         :: x

! Local variable
REAL  :: dev

n = n - 1
dev = x - mean
mean = mean - dev/n
sumsq = sumsq - dev*(x - mean)

RETURN
END SUBROUTINE downdate



SUBROUTINE replace(n, mean, sumsq, x, y)
! x is the value to be removed and replaced with the value y.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: mean, sumsq
REAL, INTENT(IN)      :: x, y

! Local variables
REAL  :: diff, devx

diff = y - x
devx = x - mean
mean = mean + diff/n
sumsq = sumsq + diff*(devx + y - mean)

RETURN
END SUBROUTINE replace
