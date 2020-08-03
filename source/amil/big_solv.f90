SUBROUTINE dple(rowk, n, b, c, ierr)
! Code from the Naval Surface Warfare Center Mathematics Library.
! Solves A.c = b where the n x n matrix is not stored but read in,
! one row at a time.
! The amount of storage required is slightly over (n/2)^2 for large n,
! or just over a quarter of that which would be required to store A.

! The user must provide the routine ROWK to read in row K of the matrix.
! The rows are read sequentially; each row is read once only.

! N.B. Arguments D & IP have been removed.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-10-20  Time: 22:49:06
 
!     ******************************************************************
!     SOLUTION OF LINEAR EQUATIONS WITH REDUCED STORAGE
!     ******************************************************************

! The partial pivot Henderson-Wassyng algorithm is used.
! Reference:
! Wassyng, A. "Solving Ax = b: A method with reduced storage requirements",
! SIAM J. Numerical Analysis, vol.19 (1982), pp.197-204.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: b(:)
REAL (dp), INTENT(OUT)  :: c(:)
INTEGER, INTENT(OUT)    :: ierr

INTERFACE
  SUBROUTINE rowk(n, k, row)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n, k
    REAL (dp), INTENT(OUT)  :: row(:)
  END SUBROUTINE rowk
END INTERFACE

! Local variables

REAL (dp)  :: bk, cj, ck, c1, dkj
INTEGER    :: i, iflag, ij, ijold, ik, ip(n-1), j, k, kjold, km1, kp1,  &
              last, lastm1, lcol, lcolp1, m, MAX, mjold, nm1, np1
REAL (dp), PARAMETER    :: zero = 0.0_dp
REAL (dp), ALLOCATABLE  :: d(:)

!     SET THE NECESSARY CONSTANTS

ierr = 0
np1 = n + 1
MAX = n * n / 4 + n + 3
ALLOCATE( d(MAX) )
k = 1
iflag = -1

!     GET THE FIRST COLUMN OF THE TRANSPOSED SYSTEM

CALL rowk(n, 1, c)
bk = b(1)

IF (n <= 1) THEN
  IF (c(1) == zero) GO TO 130
  c(1) = bk / c(1)
  RETURN
END IF

!     FIND THE PIVOT FOR COLUMN 1

m = 1
DO  i = 2, n
  IF (ABS(c(m)) < ABS(c(i))) m = i
END DO

ip(1) = m
c1 = c(m)
c(m) = c(1)
c(1) = c1
IF (c(1) /= zero) THEN
  
!     FIND THE FIRST ELEMENTARY MATRIX AND STORE IT IN D
  
  DO  i = 2, n
    d(i-1) = -c(i) / c(1)
  END DO
  d(n) = bk / c(1)
  
!     K LOOP - EACH K FOR A NEW COLUMN OF THE TRANSPOSED SYSTEM
  
  DO  k = 2, n
    kp1 = k + 1
    km1 = k - 1
    
!       GET COLUMN K
    
    CALL rowk(n, k, c)
    DO  j = 1, km1
      m = ip(j)
      cj = c(j)
      c(j) = c(m)
      c(m) = cj
    END DO
    bk = b(k)
    
    iflag = -iflag
    lcol = np1 - k
    lcolp1 = lcol + 1
    lastm1 = 1
    last = MAX - n + k
    IF (k /= 2) THEN
      
      lastm1 = MAX - n + km1
      IF (iflag < 0) last = last - n + k - 2
      IF (iflag > 0) lastm1 = lastm1 - n + k - 3
    END IF
    
!     J LOOP - EFFECT OF COLUMNS 1 TO K-1 OF L-INVERSE
    
    DO  j = 1, km1
      cj = c(j)
      ij = (j-1) * lcolp1
      IF (j == km1) ij = lastm1 - 1
      
!     I LOOP - EFFECT OF L-INVERSE ON ROWS K TO N+1
      
      DO  i = k, n
        ij = ij + 1
        c(i) = c(i) + d(ij) * cj
      END DO
      bk = bk - d(ij+1) * cj
    END DO
    
!       K=N CASE
    
    m = k
    IF (k >= n) THEN
      IF (c(k) == zero) GO TO 130
      d(last) = bk / c(k)
    ELSE
      
!     FIND THE PIVOT
      
      DO  i = kp1, n
        IF (ABS(c(m)) < ABS(c(i))) m = i
      END DO
      
      ip(k) = m
      ck = c(m)
      c(m) = c(k)
      c(k) = ck
      IF (c(k) == zero) GO TO 130
      
!     FIND THE K-TH ELEMENTARY MATRIX
      
      ik = last
      DO  i = kp1, n
        d(ik) = -c(i) / c(k)
        ik = ik + 1
      END DO
      d(ik) = bk / c(k)
    END IF
    
!     FORM THE PRODUCT OF THE ELEMENTARY MATRICES
    
    DO  j = 1, km1
      kjold = j * lcolp1 + k - np1
      mjold = kjold + m - k
      ij = (j-1) * lcol
      ijold = ij + j
      IF (j == km1) THEN
        
        kjold = lastm1
        mjold = lastm1 + m - k
        ijold = lastm1
      END IF
      
      ik = last - 1
      dkj = d(mjold)
      d(mjold) = d(kjold)
      DO  i = kp1, np1
        ij = ij + 1
        ijold = ijold + 1
        ik = ik + 1
        d(ij) = d(ijold) + d(ik) * dkj
      END DO
    END DO
  END DO
  
  last = MAX
  IF (iflag < 0) last = MAX - 2
  d(n) = d(last)
  
!     INSERT THE SOLUTION IN C
  
  c(1:n) = d(1:n)
  
  nm1 = n - 1
  DO  i = 1, nm1
    k = n - i
    m = ip(k)
    ck = c(k)
    c(k) = c(m)
    c(m) = ck
  END DO
  RETURN
END IF

!     THE SYSTEM IS SINGULAR

130 ierr = k
RETURN
END SUBROUTINE dple



SUBROUTINE rowk(n, k, row)
IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)     :: n, k
REAL (dp), INTENT(OUT)  :: row(:)

READ(UNIT=9) row(1:n)

RETURN
END SUBROUTINE rowk



PROGRAM test_big_solve
IMPLICIT NONE

INTERFACE
  SUBROUTINE dple(rowk, n, b, c, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: b(:)
    REAL (dp), INTENT(OUT)  :: c(:)
    INTEGER, INTENT(OUT)    :: ierr

    INTERFACE
      SUBROUTINE rowk(n, k, row)
        IMPLICIT NONE
        INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
        INTEGER, INTENT(IN)     :: n, k
        REAL (dp), INTENT(OUT)  :: row(:)
      END SUBROUTINE rowk
    END INTERFACE
  END SUBROUTINE dple

  SUBROUTINE rowk(n, k, row)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n, k
    REAL (dp), INTENT(OUT)  :: row(:)
  END SUBROUTINE rowk
END INTERFACE

INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER                 :: i, ierr, n, row
REAL (dp), ALLOCATABLE  :: b(:), c(:)
REAL (dp)               :: errmax, total
REAL                    :: finish, start

WRITE(*, '(a)', ADVANCE='NO') ' Enter no. of rows (n): '
READ(*, *) n
ALLOCATE( b(n), c(n) )

! The solution will be:
! c = 1, 2, 3, ..., n

! Generate a matrix of random numbers which is written to a scratch file

OPEN(UNIT=9, STATUS='SCRATCH', FORM='UNFORMATTED')
DO row = 1, n
  CALL RANDOM_NUMBER( c )
  total = 0
  DO i = 1, n
    total = total + i * c(i)
  END DO
  b(row) = total
  WRITE(9) c
END DO
REWIND (9)

CALL CPU_TIME(start)
CALL dple(rowk, n, b, c, ierr)
CALL CPU_TIME(finish)

WRITE(*, '(a, F9.2, a)') ' Time taken = ', finish-start, 'sec.'
errmax = 0.0_dp
DO i = 1, n
  errmax = MAX( errmax, ABS(c(i)-i) )
END DO
WRITE(*, '(a, g12.3)') ' Max. error = ', errmax

STOP

END PROGRAM test_big_solve
