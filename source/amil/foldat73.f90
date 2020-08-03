PROGRAM Fold_at_column73
! The fixed format for Fortran reserves columns 73-80 for sequence numbers.
! There is a Fortran 77 compiler which allows the code to flow beyond
! column 72.   It is madness to use such a non-standard extension as the
! code is not portable.

! This program takes such code and folds long lines at column 73.
! It is assumed that there are no sequence numbers in columns 73-80;
! if there are any, then they are folded as well!

! Programmer:  Alan Miller
! http://users.bigpond.net.au/amiller/
! http://www.ozemail.com.au/~milleraj
! Latest revision - 13 January 2001


IMPLICIT NONE
CHARACTER (LEN=50)  :: infile, outfile
INTEGER             :: iostatus, pos
CHARACTER (LEN=132) :: text


WRITE(*, '(a)', ADVANCE='NO') ' Enter name of input file: '
READ(*, *) infile
OPEN(8, FILE=infile, STATUS='OLD', IOSTAT=iostatus)
IF (iostatus /= 0) THEN
  WRITE(*, '(a)') ' ** Unable to open file: ' // TRIM(infile)
  STOP
END IF

! Give the extension .STD to the output file
! User can change the name later!
pos = INDEX(TRIM(infile), '.')
IF (pos == 0) THEN
  outfile = TRIM(infile) // '.STD'
ELSE
  outfile = infile(:pos-1) // '.STD'
END IF
OPEN(9, FILE=outfile)

! Now read the infile and write to the outfile, folding long lines
DO
  READ(8, '(a)', IOSTAT=iostatus) text
  IF (iostatus /= 0) EXIT
  IF (LEN_TRIM(text) <= 72) THEN
    WRITE(9, '(a)') TRIM(text)
  ELSE

! If this is not a FORMAT statement, try to fold after a blank or comma
    IF (INDEX(text(7:20), 'FORMAT') > 0 .OR.  &
        INDEX(text(7:20), 'format') > 0) THEN
      WRITE(9, '(a)') text(1:72)
      WRITE(9, '(a)') '     +' // TRIM(text(73:))
    ELSE
      pos = SCAN(text(65:72), ' ,', BACK=.TRUE.)
      IF (pos == 0) THEN
        pos = 73
      ELSE
        pos = pos + 65
      END IF
        WRITE(9, '(a)') text(1:pos-1)
        WRITE(9, '(a)') '     +' // TRIM(text(pos:))
    END IF
  END IF
END DO

STOP
END PROGRAM Fold_at_column73
