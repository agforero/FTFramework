SUBROUTINE dple(rowk, n, b, c, ierr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-16  Time: 12:26:32
 
!     ******************************************************************
!     SOLUTION OF LINEAR EQUATIONS WITH REDUCED STORAGE
!     ******************************************************************

! Uses the Henderson-Wassyng partial pivot algorithm.
! Wassyng, A. 'Solving Ax = b: A method with reduced storage requirements',
! SIAM J. Numerical Analysis, vol.19 (1982), pp. 197-204.

! The user must provide a routine ROWK to return the requested row of the
! matrix A.

! N.B. Arguments D and IP have been removed.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: b(n)
REAL (dp), INTENT(OUT)  :: c(n)
INTEGER, INTENT(OUT)    :: ierr

! EXTERNAL rowk
INTERFACE
  SUBROUTINE rowk(n, k, r)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)     :: n, k
    REAL (dp), INTENT(OUT)  :: r(:)
  END SUBROUTINE rowk
END INTERFACE

! Local variables

REAL (dp)  :: bk, cj, ck, c1, dkj
REAL (dp), PARAMETER  :: zero = 0.0_dp
REAL (dp)  :: wk(n*n/4 + n + 3)
INTEGER    :: i, iflag, ij, ijold, ik, iwk(n), j, k, kjold, km1, kp1,   &
              last, lastm1, lcol, lcolp1, m, maxwk, mjold, nm1, np1

!     SET THE NECESSARY CONSTANTS

ierr = 0
maxwk = n * n / 4 + n + 3
np1 = n + 1
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

iwk(1) = m
c1 = c(m)
c(m) = c(1)
c(1) = c1
IF (c(1) /= zero) THEN
  
!     FIND THE FIRST ELEMENTARY MATRIX AND STORE IT IN D
  
  DO  i = 2, n
    wk(i-1) = -c(i) / c(1)
  END DO
  wk(n) = bk / c(1)
  
!     K LOOP - EACH K FOR A NEW COLUMN OF THE TRANSPOSED SYSTEM
  
  DO  k = 2, n
    kp1 = k + 1
    km1 = k - 1
    
!       GET COLUMN K
    
    CALL rowk(n, k, c)
    DO  j = 1, km1
      m = iwk(j)
      cj = c(j)
      c(j) = c(m)
      c(m) = cj
    END DO
    bk = b(k)
    
    iflag = -iflag
    lcol = np1 - k
    lcolp1 = lcol + 1
    lastm1 = 1
    last = maxwk - n + k
    IF (k /= 2) THEN
      
      lastm1 = maxwk - n + km1
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
        c(i) = c(i) + wk(ij) * cj
      END DO
      bk = bk - wk(ij+1) * cj
    END DO
    
!       K=N CASE
    
    m = k
    IF (k >= n) THEN
      IF (c(k) == zero) GO TO 130
      wk(last) = bk / c(k)
    ELSE
      
!     FIND THE PIVOT
      
      DO  i = kp1, n
        IF (ABS(c(m)) < ABS(c(i))) m = i
      END DO
      
      iwk(k) = m
      ck = c(m)
      c(m) = c(k)
      c(k) = ck
      IF (c(k) == zero) GO TO 130
      
!     FIND THE K-TH ELEMENTARY MATRIX
      
      ik = last
      DO  i = kp1, n
        wk(ik) = -c(i) / c(k)
        ik = ik + 1
      END DO
      wk(ik) = bk / c(k)
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
      dkj = wk(mjold)
      wk(mjold) = wk(kjold)
      DO  i = kp1, np1
        ij = ij + 1
        ijold = ijold + 1
        ik = ik + 1
        wk(ij) = wk(ijold) + wk(ik) * dkj
      END DO
    END DO
  END DO
  
  last = maxwk
  IF (iflag < 0) last = maxwk - 2
  wk(n) = wk(last)
  
!     INSERT THE SOLUTION IN C
  
  c(1:n) = wk(1:n)
  
  nm1 = n - 1
  DO  i = 1, nm1
    k = n - i
    m = iwk(k)
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



SUBROUTINE rowk(n, k, r)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
INTEGER, INTENT(IN)     :: n, k
REAL (dp), INTENT(OUT)  :: r(:)

INTEGER    :: i
REAL (dp)  :: temp

r(k) = 1.0_dp
IF (k > 1) THEN
  temp = 1.0_dp
  DO i = 1, k-1
    temp = 0.99_dp * temp
    r(k-i) = temp
  END DO
END IF

IF (k < n) THEN
  temp = 1.0_dp
  DO i = 1, n-k
    temp = 0.99_dp * temp
    r(k+i) = temp
  END DO
END IF

RETURN
END SUBROUTINE rowk



PROGRAM Test_DPLE
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

REAL (dp)  :: a(100), x(100), b(100), error, temp, soln(100)
INTEGER    :: i, ierr, row

INTERFACE
  SUBROUTINE dple(rowk, n, b, c, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: b(n)
    REAL (dp), INTENT(OUT)  :: c(n)
    INTEGER, INTENT(OUT)    :: ierr
    INTERFACE
      SUBROUTINE rowk(n, k, r)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
        INTEGER, INTENT(IN)     :: n, k
        REAL (dp), INTENT(OUT)  :: r(:)
      END SUBROUTINE rowk
    END INTERFACE
  END SUBROUTINE dple

  SUBROUTINE rowk(n, k, r)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)     :: n, k
    REAL (dp), INTENT(OUT)  :: r(:)
  END SUBROUTINE rowk
END INTERFACE

! Generate the solution using the random number generator
CALL RANDOM_NUMBER(x)
! Calculate the right-hand side
DO row = 1, 100
  CALL rowk(100, row, a)
  b(row) = DOT_PRODUCT(a, x)
END DO

CALL dple(rowk, 100, b, soln, ierr)
IF (ierr == 0) THEN
  temp = 0.0
  DO i = 1, 100
    error = ABS(soln(i) - x(i))
    temp = MAX(temp, error)
  END DO
  WRITE(*, '(a, g13.5)') ' Max. error = ', temp
ELSE
  WRITE(*, *) ' Error = ', ierr
END IF

STOP
END PROGRAM Test_DPLE
