SUBROUTINE udruns(x, n, uv, dv, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-04  Time: 10:53:11

!      ALGORITHM AS 157 APPL. STATIST. (1981), VOL. 30, NO. 1

!      THE RUNS-UP AND RUNS-DOWN TEST

! N.B. The sample size must be at least 4000.

IMPLICIT NONE

REAL, INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: uv
REAL, INTENT(OUT)     :: dv
INTEGER, INTENT(OUT)  :: ifault

! Local variables

INTEGER :: i, j, ucount(6), dcount(6), ru, rd
REAL    :: b(6)

!      SET-UP THE A AND B MATRICES

REAL :: a( 6, 6 ) = RESHAPE( (/  &
         4529.4,  9044.9, 13568.0,  18091.0,  22615.0,  27892.0,    &
         9044.9, 18097.0, 27139.0,  36187.0,  45234.0,  55789.0,    &
        13568.0, 27139.0, 40721.0,  54281.0,  67852.0,  83685.0,    &
        18091.0, 36187.0, 54281.0,  72414.0,  90470.0, 111580.0,    &
        22615.0, 45234.0, 67852.0,  90470.0, 113262.0, 139476.0,    &
        27892.0, 55789.0, 83685.0, 111580.0, 139476.0, 172860.0 /), &
        (/ 6, 6 /) )

ifault = 0
IF (n < 4000) GO TO 500

b(1) = 1.0 / 6.0
b(2) = 5.0 / 24.0
b(3) = 11.0 / 120.0
b(4) = 19.0 / 720.0
b(5) = 29.0 / 5040.0
b(6) = 1.0 / 840.0

ucount(1:6) = 0
dcount(1:6) = 0

!    THE LOOP THAT ENDS AT LINE 300 DETERMINES THE NUMBER OF
!    RUNS-UP AND RUNS-DOWN OF LENGTH I FOR I=1(1)5 AND THE NUMBER
!    OF RUNS-UP AND RUNS-DOWN OF LENGTH GREATER OR EQUAL TO 6

ru = 1
rd = 1
DO  j = 2, n
  
!    THE FOLLOWING STATEMENT TESTS FOR BOTH RUNS-UP AND RUNS-DOWN BREAK-POINTS.
!    IF A RUN-DOWN BREAK-POINT IS DETECTED - GO TO 200, OTHERWISE - GO TO 150.
!    ( RU AND RD ACT AS LOCAL COUNTERS FOR THE NUMBER OF RUNS-UP AND
!    RUNS-DOWN RESPECTIVELY.)
!    A TEST IS ALSO MADE FOR DATA TIES BETWEEN ADJACENT ELEMENTS.
!    IF A TIE IS DETECTED - GO TO 600.
  
  IF (x(j) - x(j - 1) < 0.0) THEN
    ucount(ru) = ucount(ru) + 1
    ru = 1
    IF (rd < 6) rd = rd + 1
  ELSE IF (x(j) - x(j - 1) == 0.0) THEN
    GO TO   600
  ELSE
    dcount(rd) = dcount(rd) + 1
    rd = 1
    IF (ru < 6) ru = ru + 1
  END IF
END DO
ucount(ru) = ucount(ru) + 1
dcount(rd) = dcount(rd) + 1

!       CALCULATE THE TEST STATISTICS UV AND DV.

uv = 0.0
dv = 0.0
DO  i = 1, 6
  DO  j = 1, 6
    uv = uv + (ucount(i) - n * b(i)) * (ucount(j) - n * b(j)) * a(i, j)
    dv = dv + (dcount(i) - n * b(i)) * (dcount(j) - n * b(j)) * a(i, j)
  END DO
END DO
uv = uv / n
dv = dv / n
RETURN

500 ifault = n
RETURN

600 ifault = 1
RETURN
END SUBROUTINE udruns
