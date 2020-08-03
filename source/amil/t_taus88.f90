PROGRAM write_random
! Write a large file of random integers for use by DIEHARD.
! This version uses Pierre L'Ecuyer's TUAS88 random number generator.

IMPLICIT NONE
CHARACTER (LEN=25)  :: fname
INTEGER             :: b(4096), i, i1, i2, i3, j

! These are unsigned integers in the C version of TAUS88
INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890

WRITE(*, '(a)', ADVANCE='NO') ' Enter name for your binary file: '
READ(*, *) fname
OPEN (1, FILE=fname, FORM='unformatted', ACCESS='direct', RECL=16384)

WRITE(*, '(a)', ADVANCE='NO') ' Enter 3 integers as random number seeds: '
READ(*, *) i1, i2, i3
CALL init_seeds(i1, i2, i3)

DO i = 1, 700
  DO j = 1, 4096
    b(j) = taus88()
  END DO
  WRITE(1, REC=i) b
END DO

STOP

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3)

INTEGER, INTENT(IN) :: i1, i2, i3

s1 = i1
s2 = i2
s3 = i3
IF (IAND(s1,-2) == 0) s1 = i1 - 1023
IF (IAND(s2,-8) == 0) s2 = i2 - 1023
IF (IAND(s3,-16) == 0) s3 = i3 - 1023

RETURN
END SUBROUTINE init_seeds



FUNCTION taus88() RESULT(irand)
! Translated from C function in:
! Reference:
! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
! generators', Math. of Comput., 65, 203-213.

! Modified to return a random integer.

! The cycle length is claimed to be about 2^(88) or about 3E+26.
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

INTEGER   :: irand

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
irand = IEOR( IEOR(s1,s2), s3)

RETURN
END FUNCTION taus88

END PROGRAM write_random
