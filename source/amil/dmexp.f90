MODULE Matrix_Exponential
!     DMEXP COMPUTES EXP(A) AND STORES IT IN Z WHERE A IS A MATRIX
!     OF ORDER N.  A IS DESTROYED BY THE ROUTINE.

! Reference:
! Ward, Robert C. (1977) 'Numerical computation of the matrix exponential with
! accuracy estimate', SIAM J. Numerical Analysis, vol.14, pp. 600-610.

! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-19  Time: 10:46:30

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


SUBROUTINE dmexp(a, n, z, ierr)
! Arguments KA, KZ & WK has been removed.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: a(:,:)
REAL (dp), INTENT(OUT)  :: z(:,:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp) :: anorm, anorm1, factor, p, q, s, s1, wk(n,n+12)
INTEGER   :: i, igh, j, k, kp1, l, ll, low, m, np1, np10

! ----------------------------------------------------------------------
!     DMEXP COMPUTES EXP(A) AND STORES IT IN Z WHERE A IS A MATRIX
!     OF ORDER N.  A IS DESTROYED BY THE ROUTINE.

!     WK IS AN ARRAY OF DIMENSION (N,N+12).  WK IS A WORK SPACE
!     FOR THE ROUTINE.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!        IERR = 0  EXP(A) WAS SUCCESSFULLY COMPUTED.
!        IERR = 1  THE NORM OF A IS TOO LARGE.
!        IERR = 2  THE PADE DENOMINATOR MATRIX IS SINGULAR.

!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN VIRGINIA
!     ---------------
!     COEFFICIENTS FOR (12,12) PADE TABLE ENTRY
!     ---------------

REAL (dp), PARAMETER  :: c(12) = (/  &
  .500000000000000000000000000000D+00, .119565217391304347826086956522D+00,  &
  .181159420289855072463768115942D-01, .194099378881987577639751552795D-02,  &
  .155279503105590062111801242236D-03, .953470633104500381388253241800D-05,  &
  .454033634811666848280120591333D-06, .166924130445465753044161982108D-07,  &
  .463678140126293758456005505855D-09, .927356280252587516912011011710D-11,  &
  .120435880552284093105455975547D-12, .772024875335154442983692150941D-15 /)
! ----------------------------------------------------------------------

ierr = 0
IF (n <= 1) THEN
  z(1,1) = EXP(a(1,1))
  RETURN
END IF

!       BALANCE A AND SELECT THE SMALLER OF THE 1-NORM
!              AND INFINITY-NORM OF THE RESULT

CALL dbal(n, a, low, igh, wk(1:n,n+12))
anorm = 0.d0
anorm1 = 0.d0
DO  j = 1, n
  s = 0.d0
  s1 = 0.d0
  DO  i = 1, n
    s = s + ABS(a(j,i))
    s1 = s1 + ABS(a(i,j))
  END DO
  anorm = MAX(s,anorm)
  anorm1 = MAX(s1,anorm1)
END DO

anorm = MIN(anorm,anorm1)
s = anorm + 0.1D0
IF (s /= anorm) THEN
  
!              SELECT THE NORMALIZATION FACTOR
  
  m = 0
  IF (anorm > 1.d0) THEN
    factor = 1.d0
    30 m = m + 1
    factor = 2.d0 * factor
    IF (anorm > factor) GO TO 30
    
!                NORMALIZE THE MATRIX A
    
    a(1:n,1:n) = a(1:n,1:n) / factor
  END IF
  
  np1 = n + 1
  np10 = n + 10
  DO  j = 1, n
    
!     COMPUTE THE J-TH COLUMN OF THE FIRST 12 POWERS OF A
    
    DO  i = 1, n
      s = 0.d0
      DO  l = 1, n
        s = s + a(i,l) * a(l,j)
      END DO
      wk(i,np1) = s
    END DO
    
    DO  k = np1, np10
      kp1 = k + 1
      DO  i = 1, n
        s = 0.d0
        DO  l = 1, n
          s = s + a(i,l) * wk(l,k)
        END DO
        wk(i,kp1) = s
      END DO
    END DO
    
!     COMPUTE THE J-TH COLUMN OF THE NUMERATOR AND DENOMINATOR
!                  OF THE PADE APPROXIMATION
    
    DO  i = 1, n
      p = 0.d0
      q = 0.d0
      k = 12
      l = n + 11
      DO  ll = 1, 11
        s = c(k) * wk(i,l)
        p = s + p
        q = s - q
        k = k - 1
        l = l - 1
      END DO
      s = c(1) * a(i,j)
      z(i,j) = p + s
      wk(i,j) = q - s
      IF (i == j) THEN
        z(i,j) = z(i,j) + 1.d0
        wk(i,j) = wk(i,j) + 1.d0
      END IF
    END DO
  END DO
  
!        CALCULATE EXP(A) BY SOLVING  WK * EXP(A) = Z
  
  CALL dpslv(n, n, wk, z, ierr)
  IF (ierr /= 0) GO TO 200
  IF (m /= 0) THEN
    
!          TAKE OUT THE EFFECT OF THE NORMALIZATION
!                   OPERATION ON EXP(A)
    
    DO  k = 1, m
      DO  j = 1, n
        DO  i = 1, n
          s = 0.d0
          DO  l = 1, n
            s = s + z(i,l) * z(l,j)
          END DO
          wk(i,j) = s
        END DO
      END DO

      z(1:n,1:n) = wk(1:n,1:n)
    END DO
  END IF
  
!            TAKE OUT THE EFFECT OF THE BALANCING
!                   OPERATION ON EXP(A)
  
  CALL dbalnv(n, z, low, igh, wk(1:n,n+12))
  RETURN
END IF

!                     ERROR RETURN

ierr = 1
RETURN
200 ierr = 2
RETURN
END SUBROUTINE dmexp



SUBROUTINE dbal(n, a, low, igh, scale)
! Argument NM has been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(:,:)
INTEGER, INTENT(OUT)       :: low
INTEGER, INTENT(OUT)       :: igh
REAL (dp), INTENT(OUT)     :: scale(:)

INTEGER :: i, j, k, l, m, jj, iexc, radix

REAL (dp) :: c, f, g, r, s, b2
LOGICAL   :: noconv, skip
!-----------------------------------------------------------------------

!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).

!     DBAL BALANCES A REAL (dp) REAL MATRIX AND ISOLATES
!     EIGENVALUES WHENEVER POSSIBLE.

!     ON INPUT-

!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL ARRAY
!          PARAMETERS AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,

!        N IS THE ORDER OF THE MATRIX,

!        A CONTAINS THE INPUT MATRIX TO BE BALANCED.

!     ON OUTPUT-

!        A CONTAINS THE BALANCED MATRIX,

!        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J) IS EQUAL TO ZERO IF
!           (1) I IS GREATER THAN J AND
!           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,

!        SCALE CONTAINS INFORMATION DETERMINING THE
!           PERMUTATIONS AND SCALING FACTORS USED.

!     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
!     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
!     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
!     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
!        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
!                 = D(J,J),      J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
!     THEN 1 TO LOW-1.

!     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.

!     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
!     DBAL IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
!     K,L HAVE BEEN REVERSED.)

!-----------------------------------------------------------------------

!     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.

radix = 2
b2 = radix ** 2
k = 1
l = n
GO TO 50

!     ********** IN-LINE PROCEDURE FOR ROW AND COLUMN EXCHANGE **********
10 scale(m) = j
IF (j /= m) THEN
  
  DO  i = 1, l
    f = a(i,j)
    a(i,j) = a(i,m)
    a(i,m) = f
  END DO
  
  DO  i = k, n
    f = a(j,i)
    a(j,i) = a(m,i)
    a(m,i) = f
  END DO
END IF

SELECT CASE ( iexc )
  CASE (    1)
    GO TO 40
  CASE (    2)
    GO TO 80
END SELECT

!     ***** SEARCH FOR ROWS ISOLATING AN EIGENVALUE AND PUSH THEM DOWN *****
40 IF (l == 1) GO TO 200
l = l - 1
!     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
skip = .FALSE.
50 DO  jj = 1, l
  j = l + 1 - jj
  
  DO  i = 1, l
    IF (i /= j) THEN
      IF (a(j,i) /= 0.d0) THEN
        skip = .TRUE.
        EXIT
      END IF
    END IF
  END DO
  IF (skip) CYCLE
  
  m = l
  iexc = 1
  GO TO 10
END DO

GO TO 90
!     **** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE AND PUSH THEM LEFT ****
80 k = k + 1

90 skip = .FALSE.
DO  j = k, l
  
  DO  i = k, l
    IF (i /= j) THEN
      IF (a(i,j) /= 0.d0) THEN
        skip = .TRUE.
        EXIT
      END IF
    END IF
  END DO
  IF (skip) CYCLE
  
  m = k
  iexc = 2
  GO TO 10
END DO

!     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
scale(k:l) = 1.d0
!     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
130 noconv = .false.

DO  i = k, l
  c = 0.d0
  r = 0.d0
  
  DO  j = k, l
    IF (j /= i) THEN
      c = c + ABS(a(j,i))
      r = r + ABS(a(i,j))
    END IF
  END DO
!     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********
  IF (c /= 0.d0 .AND. r /= 0.d0) THEN
    g = r / radix
    f = 1.d0
    s = c + r
    150 IF (c < g) THEN
      f = f * radix
      c = c * b2
      GO TO 150
    END IF
    g = r * radix
    160 IF (c >= g) THEN
      f = f / radix
      c = c / b2
      GO TO 160
    END IF
!     ********** NOW BALANCE **********
    IF ((c+r)/f < 0.95D0*s) THEN
      g = 1.d0 / f
      scale(i) = scale(i) * f
      noconv = .true.
      
      DO  j = k, n
        a(i,j) = a(i,j) * g
      END DO
      
      DO  j = 1, l
        a(j,i) = a(j,i) * f
      END DO
    END IF
  END IF
  
END DO

IF (noconv) GO TO 130

200 low = k
igh = l
RETURN
END SUBROUTINE dbal



SUBROUTINE dbalnv(n, z, low, igh, scale)
! Argument NZ has been removed.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: z(:,:)
INTEGER, INTENT(IN)     :: low
INTEGER, INTENT(IN)     :: igh
REAL (dp), INTENT(IN)   :: scale(:)

INTEGER   :: i, j, k, ii
REAL (dp) :: s

!-----------------------------------------------------------------------
!        GIVEN A MATRIX A OF ORDER N. DBAL TRANSFORMS A INTO THE MATRIX B
!     BY THE SIMILARITY TRANSFORMATION
!             B = D**(-1)*TRANSPOSE(P)*A*P*D
!     WHERE D IS A DIAGONAL MATRIX AND P A PERMUTATION MATRIX.
!     THE INFORMATION CONCERNING D AND P IS STORED IN IGH, LOW, AND SCALE.
!     THE ORDER IN WHICH THE INTERCHANGES WERE MADE IS N TO IGH + 1,
!     AND THEN 1 TO LOW - 1.

!        Z IS A MATRIX OF ORDER N. DBALNV TRANSFORMS Z INTO THE
!     MATRIX W USING THE INVERSE SIMILARITY TRANSFORM
!             W = P*D*Z*D**(-1)*TRANSPOSE(P)

!     ON INPUT-

!        NZ IS THE ROW DIMENSION OF THE MATRIX Z IN THE CALLING PROGRAM,

!        N IS THE ORDER OF THE MATRIX,

!        LOW AND IGH ARE INTEGERS DETERMINED BY  DBAL,

!        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
!          AND SCALING FACTORS USED BY  DBAL,

!     ON OUTPUT-

!        Z CONTAINS THE TRANSFORMED MATRIX W

!-----------------------------------------------------------------------

IF (igh /= low) THEN
  
  DO  i = low, igh
    s = scale(i)
    z(i,1:n) = z(i,1:n) * s
  END DO
  
  DO  j = low, igh
    s = 1.d0 / scale(j)
    z(1:n,j) = z(1:n,j) * s
  END DO
END IF

!     ********- FOR I=LOW-1 STEP -1 UNTIL 1,
!               IGH+1 STEP 1 UNTIL N DO -- **********

DO  ii = 1, n
  i = ii
  IF (i < low.OR.i > igh) THEN
    IF (i < low) i = low - ii
    k = scale(i)
    IF (k /= i) THEN
      
      DO  j = 1, n
        s = z(i,j)
        z(i,j) = z(k,j)
        z(k,j) = s
      END DO
      
      DO  j = 1, n
        s = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = s
      END DO
    END IF
  END IF
END DO
RETURN
END SUBROUTINE dbalnv



SUBROUTINE dpslv(n, m, a, b, ierr)
! Arguments KA & KB have been removed.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: a(:,:)
REAL (dp), INTENT(IN OUT)  :: b(:,:)
INTEGER, INTENT(OUT)       :: ierr

REAL (dp) :: p, t
INTEGER   :: i, j, k, km1, kp1, l, nm1

!     ------------------------------------------------------------------
!     PARTIAL PIVOT SOLUTION OF A*X = B WHERE A IS A MATRIX OF
!     ORDER N AND B IS A MATRIX HAVING N ROWS AND M COLUMNS.
!     THE SOLUTION MATRIX X IS STORED IN B.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!        IERR = 0   THE EQUATIONS HAVE BEEN SOLVED.
!        IERR = J   THE J-TH PIVOT ELEMENT WAS FOUND TO BE 0.
!     ------------------------------------------------------------------
ierr = 0
nm1 = n - 1
IF (nm1 /= 0) THEN
  DO  k = 1, nm1
    
!     SEARCH FOR THE K-TH PIVOT ELEMENT
    
    p = 0.d0
    DO  i = k, n
      t = ABS(a(i,k))
      IF (p < t) THEN
        p = t
        l = i
      END IF
    END DO
    IF (p == 0.d0) GO TO 120
    IF (k /= l) THEN
      
!     INTERCHANGE ROWS K AND L
      
      DO  j = k, n
        t = a(k,j)
        a(k,j) = a(l,j)
        a(l,j) = t
      END DO
      DO  j = 1, m
        t = b(k,j)
        b(k,j) = b(l,j)
        b(l,j) = t
      END DO
    END IF
    
!     ELIMINATE THE COEFFICIENTS OF X(K) IN ROWS I = K+1,...,N
    
    p = a(k,k)
    kp1 = k + 1
    DO  i = kp1, n
      t = a(i,k) / p
      DO  j = kp1, n
        a(i,j) = a(i,j) - t * a(k,j)
      END DO
      DO  j = 1, m
        b(i,j) = b(i,j) - t * b(k,j)
      END DO
    END DO
  END DO
  IF (a(n,n) == 0.d0) GO TO 130
  
!     BACKSOLVE THE TRIANGULAR SET OF EQUATIONS
  
  DO  j = 1, m
    k = n
    km1 = nm1
    DO  l = 2, n
      b(k,j) = b(k,j) / a(k,k)
      t = b(k,j)
      DO  i = 1, km1
        b(i,j) = b(i,j) - t * a(i,k)
      END DO
      k = km1
      km1 = k - 1
    END DO
    b(1,j) = b(1,j) / a(1,1)
  END DO
  RETURN
END IF

!     CASE WHEN N = 1

IF (a(1,1) /= 0.d0) THEN
  b(1,1:m) = b(1,1:m) / a(1,1)
  RETURN
END IF

!     ERROR RETURN

ierr = 1
RETURN
120 ierr = k
RETURN
130 ierr = n
RETURN
END SUBROUTINE dpslv

END MODULE Matrix_Exponential
