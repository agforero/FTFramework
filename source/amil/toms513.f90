SUBROUTINE trans(a, m, n, mn, iok)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-03  Time: 17:05:40
 
! *****
!  ALGORITHM 380 - REVISED
! *****
!  A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH CONTAINS THE MXN
!  MATRIX TO BE TRANSPOSED (STORED COLUMNWISE).
!  IOK INDICATES THE SUCCESS OR FAILURE OF THE ROUTINE.
!  NORMAL RETURN  IOK=0
!  ERRORS         IOK=-1 ,MN NOT EQUAL TO M*N
!                 IOK > 0, (SHOULD NEVER OCCUR), IN THIS CASE WE SET IOK
!  EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH IS COMPLETED BUT SOME LOOPS
!  HAVE NOT BEEN MOVED.
!  NOTE * MOVE(I) WILL STAY ZERO FOR FIXED POINTS

IMPLICIT NONE

REAL, INTENT(IN OUT)  :: a(:)
INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: mn
INTEGER, INTENT(OUT)  :: iok

! Local variables
INTEGER  :: i, i1, i2, i1c, i2c, im, ir0, ir1, ir2, j, j1, k, kmi,  &
            MAX, n1, ncount
REAL     :: b, c, d

!  MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK USED TO STORE INFORMATION
!  TO SPEED UP THE PROCESS.  THE VALUE IWRK=(M+N)/2 IS RECOMMENDED.
INTEGER  :: iwrk
INTEGER  :: move( (m+n)/2 )

! CHECK ARGUMENTS AND INITIALIZE.
IF (m >= 2 .AND. n >= 2) THEN
  IF (mn /= m*n) GO TO 140
  IF (m == n) GO TO 90
  ncount = 2
  k = mn - 1
  iwrk = (m+n)/2
  move(1:iwrk) = 0
  IF (m >= 3 .AND. n >= 3) THEN

! CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM FOR GCD(M-1,N-1).
    ir2 = m - 1
    ir1 = n - 1
    20     ir0 = MOD(ir2,ir1)
    ir2 = ir1
    ir1 = ir0
    IF (ir0 /= 0) GO TO 20
    ncount = ncount + ir2 - 1
  END IF

! SET INITIAL VALUES FOR SEARCH
  i = 1
  im = m
! AT LEAST ONE LOOP MUST BE RE-ARRANGED
  GO TO 60

! SEARCH FOR LOOPS TO REARRANGE
  30   MAX = k - i
  i = i + 1
  IF (i > MAX) GO TO 120
  im = im + m
  IF (im > k) im = im - k
  i2 = im
  IF (i == i2) GO TO 30
  IF (i > iwrk) GO TO 50
  IF (move(i) == 0) GO TO 60
  GO TO 30

  40   i2 = m * i1 - k * (i1/n)
  50   IF (i2 > i.AND.i2 < MAX) THEN
    i1 = i2
    GO TO 40
  END IF
  IF (i2 /= i) GO TO 30

! REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
  60   i1 = i
  kmi = k - i
  b = a(i1+1)
  i1c = kmi
  c = a(i1c+1)
  70   i2 = m * i1 - k * (i1/n)
  i2c = k - i2
  IF (i1 <= iwrk) move(i1) = 2
  IF (i1c <= iwrk) move(i1c) = 2
  ncount = ncount + 2
  IF (i2 /= i) THEN
    IF (i2 /= kmi) THEN
      a(i1+1) = a(i2+1)
      a(i1c+1) = a(i2c+1)
      i1 = i2
      i1c = i2c
      GO TO 70
    END IF

! FINAL STORE AND TEST FOR FINISHED
    d = b
    b = c
    c = d
  END IF
  a(i1+1) = b
  a(i1c+1) = c
  IF (ncount < mn) GO TO 30
END IF

! NORMAL RETURN
80 iok = 0
RETURN

! IF MATRIX IS SQUARE, EXCHANGE ELEMENTS A(I,J) AND A(J,I).
90 n1 = n - 1
DO  i = 1, n1
  j1 = i + 1
  DO  j = j1, n
    i1 = i + (j-1) * n
    i2 = j + (i-1) * m
    b = a(i1)
    a(i1) = a(i2)
    a(i2) = b
  END DO
END DO
GO TO 80

! ERROR RETURNS.
120 iok = i
130 RETURN
140 iok = -1
GO TO 130
END SUBROUTINE trans


