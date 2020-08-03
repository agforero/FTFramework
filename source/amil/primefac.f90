!-----------------------------------------------------------------------
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-14  Time: 11:42:45
 
!  PRMFAC

!  DECOMPOSE A NUMBER INTO ITS PRIME FACTORS.

!  CODED AT MADISON ACADEMIC COMPUTING CENTER,
!  UNIVERSITY OF WISCONSIN, MADISON

!  FORTRAN 77 VERSION 88.09

!-----------------------------------------------------------------------

!--------CALLING SEQUENCE

!  CALL PRMFAC (NUMBER, IOPT, NPRM, IPRM, IEXP, *ERR)

!  NUMBER - INTEGER CONSTANT OR VARIABLE, NUMBER TO BE DECOMPOSED
!           INTO PRIME FACTORS.  NUMBER .GE. 2.

!  IOPT   - INTEGER CONSTANT OR VARIABLE, CONTROLS PRINTING OF RESULTS.
!           IOPT .EQ. 0 - RESULTS ARE NOT PRINTED.
!           IOPT .NE. 0 - RESULTS ARE PRINTED.

!  NPRM   - INTEGER VARIABLE, WILL CONTAIN THE NO. OF DISTINCT PRIME
!           FACTORS OF THE NUMBER.

!  IPRM   - INTEGER ARRAY OF SIZE AT LEAST 9, WILL CONTAIN THE PRIME
!           FACTORS OF THE NUMBER.

!  IEXP   - INTEGER ARRAY OF SIZE AT LEAST 9, WILL CONTAIN THE
!           EXPONENTS OF THE CORRESPONDING PRIME FACTORS.

!  *ERR   - STATEMENT NUMBER PRECEDED WITH AN ASTERISK (*).
!           IF AN ERROR CONDITION IS DETECTED, CONTROL IS RETURNED
!           TO THAT STATEMENT IN THE CALLING PROGRAM.
!           This has been deleted below.

!--------NOTES

!  (1)  UPON RETURN FROM PRMFAC,
!       NUMBER = IPRM(1)**IEXP(1) * IPRM(2)**IEXP(2) * ... *
!               IPRM(NPRM)**IEXP(NPRM)

!  (2)  A NUMBER REPRESENTED BY A (SINGLE-PRECISION) INTEGER
!       VALUE ON THE VMS VAX CLUSTER CAN HAVE AT MOST 9 DISTINCT
!       PRIME FACTORS.  ON MACHINES WHERE THE MAXIMUM INTEGER IS
!       LARGER THAN 2**31 - 1, IPRM AND IEXP WOULD, IN GENERAL,
!       HAVE TO BE DIMENSIONED LARGER SINCE LARGER NUMBERS MAY
!       HAVE MORE THAN 9 DISTINCT PRIME FACTORS.

!-----------------------------------------------------------------------

SUBROUTINE prmfac (NUMBER, iopt, nprm, iprm, iexp)

!--------PARAMETERS

IMPLICIT NONE

INTEGER, INTENT(IN)   :: NUMBER
INTEGER, INTENT(IN)   :: iopt
INTEGER, INTENT(OUT)  :: nprm
INTEGER, INTENT(OUT)  :: iprm(:)
INTEGER, INTENT(OUT)  :: iexp(:)

!--------LOCAL VARIABLES
INTEGER            :: div, indx, isize, j, n, offset, olddiv, quo, rem

!---------DATA TO OBTAIN TRIAL DIVISORS 2, 3, 5, 7 AND ALL
!         HIGHER NUMBERS NOT DIVISIBLE BY 2, 3, 5, 7.
INTEGER, PARAMETER  :: base(52) = (/   &
    211, 209, 199, 197, 193, 191, 187, 181, 179, 173, 169, 167, 163, &
    157, 151, 149, 143, 139, 137, 131, 127, 121, 113, 109, 107, 103, &
    101,  97,  89,  83,  79,  73,  71,  67,  61,  59,  53,  47,  43, &
     41,  37,  31,  29,  23,  19,  17,  13,  11,   7,   5,   3,   2 /)

!--------FORMAT TEMPLATE
CHARACTER (LEN=16) :: i4mat = '(1X,A8,1X,--I--)'
CHARACTER (LEN=8)  :: txtlst(2) = (/ 'PRIME   ', 'EXPONENT' /)

!--------LOGICAL UNIT FOR STANDARD OUTPUT
INTEGER, PARAMETER :: prn = 6


!--------CHECK NUMBER.  MUST BE .GE. 2
IF (NUMBER < 2) GO TO 140

!--------INITIALIZATIONS.
j = 0
n = NUMBER
olddiv = 0
offset = 0
indx = 53

!--------GET NEXT TRIAL DIVISOR.
DO
  indx = indx - 1
  IF (indx <= 0) THEN
    indx = 48
    offset = offset + 210
  END IF
  div = offset + base(indx)

!--------TEST TRIAL DIVISOR.
  DO
    quo = n / div
    rem = n - quo*div
    IF (rem /= 0) EXIT

!--------FACTOR FOUND, ZERO REMAINDER.
    n = quo
    IF (div <= olddiv) THEN

!--------MULTIPLE FACTOR.
      iexp(j) = iexp(j) + 1
      CYCLE
    END IF

!--------NEW FACTOR.
    j = j + 1
    iprm(j) = div
    iexp(j) = 1
    olddiv = div
  END DO

!--------NOT A FACTOR, POSITIVE REMAINDER.  CHECK DIVISOR SIZE.
  IF (div >= quo) EXIT
END DO

!--------FINISHED, WHAT ISN'T FACTORED IS A PRIME (OR 1).
IF (n > 1) THEN
  j = j + 1
  iexp(j) = 1
  iprm(j) = n
END IF
nprm = j

!--------PRINT RESULTS IF REQUESTED.
IF (iopt == 0) GO TO 130
WRITE (prn,70)
70 FORMAT ( / ' ', 7('----------') )
IF (nprm /= 1 .OR. iexp(1) /= 1) GO TO 90

!--------NUMBER IS PRIME
WRITE (prn,80) NUMBER
80 FORMAT ( / t3, i12, ' IS PRIME' )
GO TO 120

!--------NUMBER IS COMPOSITE
90 WRITE (prn,100) NUMBER
100 FORMAT ( / t3, i12, ' FACTORS AS FOLLOWS' / )

!--------DETERMINE SUITABLE FORMAT FOR FACTORS AND EXPONENTS AND PRINT.
isize = LOG10 (REAL(iprm(nprm)) + 0.5D0)
isize = MAX (6,isize+3)
WRITE (i4mat(11:12), 110) nprm
110 FORMAT (i2.2)
WRITE (i4mat(14:15), 110) isize
WRITE (prn, i4mat) txtlst(1), iprm(1:nprm)
WRITE (prn, i4mat) txtlst(2), iexp(1:nprm)

120 WRITE (prn,70)

!--------NORMAL EXIT
130 RETURN

!--------ERROR, ISSUE MESSAGE AND TAKE ERROR EXIT.
140 WRITE (prn,150) NUMBER
150 FORMAT (' *** ERROR IN PRMFAC, NUMBER =', i12, ', NUMBER MUST BE >= 2')

!--------ERROR EXIT
RETURN

!--------END OF PRMFAC
END SUBROUTINE prmfac



PROGRAM find_prime_factors
IMPLICIT NONE

INTERFACE
  SUBROUTINE prmfac (NUMBER, iopt, nprm, iprm, iexp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: NUMBER
    INTEGER, INTENT(IN)   :: iopt
    INTEGER, INTENT(OUT)  :: nprm
    INTEGER, INTENT(OUT)  :: iprm(:)
    INTEGER, INTENT(OUT)  :: iexp(:)
  END SUBROUTINE prmfac
END INTERFACE

INTEGER  :: iexp(10), iopt, iprm(10), nprm, NUMBER

DO
  iopt = 1
  WRITE(*, '(a)', ADVANCE='NO') ' Enter number to be factored: '
  READ(*, *) NUMBER
  CALL prmfac (NUMBER, iopt, nprm, iprm, iexp)
  WRITE(*, *)
END DO

STOP
END PROGRAM find_prime_factors
