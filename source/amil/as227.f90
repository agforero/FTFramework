SUBROUTINE gcount(n, apply, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-04-26  Time: 20:20:16

! ALGORITHM AS227  APPL. STATIST. (1987) VOL. 36, NO. 2, pp. 245-9.

! Generates all possible N-bit binary codes, and applies a users
! procedure for each code generated.

! Translated from Algol 60.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(OUT)  :: ifault

! EXTERNAL apply

INTERFACE
  SUBROUTINE apply(n, change, STATUS)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n, change
    LOGICAL, INTENT(IN)  :: STATUS(n)
  END SUBROUTINE apply
END INTERFACE

! Local variables

INTEGER  :: change, i, tpoint(n)
LOGICAL  :: STATUS(n)

IF (n < 1) THEN
  ifault = 1
  RETURN
END IF
ifault = 0

! Initialize and make first call to user's routine.

DO  i = 1, n
  STATUS(i) = .false.
  tpoint(i) = i + 1
END DO
CALL apply(n, n, STATUS)

! Generate a new code.   The user's routine is called twice each
! cycle; the first time the bit which changes is bit 1.

20 IF (STATUS(1)) THEN
  STATUS(1) = .false.
  change = tpoint(2)
ELSE
  STATUS(1) = .true.
  change = 2
END IF
CALL apply(n, 1, STATUS)

! Check if count exhausted.

IF (change > n) RETURN

IF (STATUS(change)) THEN
  STATUS(change) = .false.
  tpoint(change) = tpoint(change+1)
ELSE
  STATUS(change) = .true.
  tpoint(change) = change + 1
END IF
CALL apply(n, change, STATUS)

GO TO 20
END SUBROUTINE gcount



SUBROUTINE print_comb(n, change, STATUS)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, change
LOGICAL, INTENT(IN)  :: STATUS(n)

WRITE(*, '(" ", 4L2)') STATUS
RETURN
END SUBROUTINE print_comb



PROGRAM t_as227
! A simple program to print out the combinations when n = 4.

IMPLICIT NONE

INTEGER  :: n = 4, ifault

INTERFACE
  SUBROUTINE gcount(n, apply, ifault)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: n
    INTEGER, INTENT(OUT)  :: ifault
    INTERFACE
      SUBROUTINE apply(n, change, STATUS)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: n, change
        LOGICAL, INTENT(IN)  :: STATUS(n)
      END SUBROUTINE apply
    END INTERFACE
  END SUBROUTINE gcount

  SUBROUTINE print_comb(n, change, STATUS)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n, change
    LOGICAL, INTENT(IN)  :: STATUS(n)
  END SUBROUTINE print_comb
END INTERFACE

CALL gcount(n, print_comb, ifault)

STOP
END PROGRAM t_as227
