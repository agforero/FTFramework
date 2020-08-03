MODULE randdp
! Replaces COMMON /randpp/ and defines dp

IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(15, 60)

REAL (dp), SAVE, PUBLIC    :: poly(101), other, offset
INTEGER, SAVE, PUBLIC      :: index

CONTAINS

!   Nick Maclaren's double precision random number generator.
!   This version, which is compatible with Lahey's ELF90 compiler,
!   is by Alan Miller ( alan @ vic.cmis.csiro.au, www.vic.cmis.csiro.au/~alan )

!   Latest revision - 18 December 1997

!   Copyright (C) 1992  N.M. Maclaren
!   Copyright (C) 1992  The University of Cambridge

!   This software may be reproduced and used freely, provided that all
!   users of it agree that the copyright holders are not liable for any
!   damage or injury caused by use of this software and that this
!   condition is passed onto all subsequent recipients of the software,
!   whether modified or not.



SUBROUTINE sdprnd (iseed)

INTEGER, INTENT(IN)  :: iseed

! Local variables
REAL (dp)            :: x
REAL (dp), PARAMETER :: xmod = 1000009711.0_dp, ymod = 33554432.0_dp
INTEGER              :: ix, iy, iz, i
LOGICAL, SAVE        :: inital = .TRUE.

!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.

IF (inital .OR. iseed /= 0) THEN
  inital = .false.
ELSE
  RETURN
END IF

!   INDEX must be initialised to an integer between 1 and 101 inclusive,
!   POLY(1...N) to integers between 0 and 1000009710 inclusive (not all 0),
!   and OTHER to a non-negative proper fraction with denominator 33554432.
!   It uses the Wichmann-Hill generator to do this.

ix = MOD(ABS(iseed), 10000) + 1
iy = 2*ix + 1
iz = 3*ix + 1
DO i = -10,101
  IF (i >= 1) poly(i) = AINT(xmod*x)
  ix = MOD(171*ix, 30269)
  iy = MOD(172*iy, 30307)
  iz = MOD(170*iz, 30323)
  x = MOD(DBLE(ix)/30269.0_dp + DBLE(iy)/30307.0_dp + DBLE(iz)/30323.0_dp, 1.0_dp)
END DO
other = AINT(ymod*x)/ymod
offset = 1.0_dp/ymod
INDEX = 1

RETURN
END SUBROUTINE sdprnd



FUNCTION dprand() RESULT(fn_val)

REAL (dp)            :: fn_val

! Local variables
REAL (dp)            :: x, y

! N.B. ymod has been removed from the previous DATA statement; it caused a
!      fatal error as it is not used.
REAL (dp), PARAMETER :: xmod = 1000009711.0_dp, xmod2 = 2000019422.0_dp, &
                        xmod4 = 4000038844.0_dp, tiny = 1.0E-17_dp,  &
                        zero = 0.0_dp, one = 1.0_dp
INTEGER              :: n
LOGICAL, SAVE        :: inital = .TRUE.

!   This returns a uniform (0,1) random number, with extremely good
!   uniformity properties.  It assumes that REAL (dp) provides
!   at least 33 bits of accuracy, and uses a power of two base.

IF (inital) THEN
  CALL sdprnd (0)
  inital = .false.
END IF

!   See [Knuth] for why this implements the algorithm described in the paper.
!   Note that this code is tuned for machines with fast REAL (dp), but
!   slow multiply and divide; many, many other options are possible.

n = INDEX - 64
IF (n <= 0) n = n + 101
x = poly(INDEX) + poly(INDEX)
x = xmod4 - poly(n) - poly(n) - x - x - poly(INDEX)
IF (x < zero) THEN
  IF (x < -xmod) x = x + xmod2
  IF (x < zero) x = x + xmod
ELSE
  IF (x >= xmod2) THEN
    x = x - xmod2
    IF (x >= xmod) x = x - xmod
  END IF
  IF (x >= xmod) x = x - xmod
END IF
poly(INDEX) = x
INDEX = INDEX + 1
IF (INDEX > 101) INDEX = INDEX - 101

!   Add in the second generator modulo 1, and force to be non-zero.
!   The restricted ranges largely cancel themselves out.

DO
  y = 37.0_dp*other + offset
  other = y - AINT(y)
  IF (other /= zero) EXIT
END DO

x = x/xmod + other
IF (x >= one) x = x - one
fn_val = x + tiny

RETURN
END FUNCTION dprand

END MODULE randdp
