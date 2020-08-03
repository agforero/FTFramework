MODULE skew_symmetric

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS


SUBROUTINE trizd(n, a, e)

! ****

!  FUNCTION  -  REDUCES A REAL SKEW-SYMMETRIC MATRIX TO A SKEW-SYMMETRIC
!                 TRIDIAGONAL MATRIX USING ORTHOGONAL SIMILARITY
!                 TRANSFORMATIONS

!  PARAMETERS

!     N        - INPUT INTEGER SPECIFYING THE ORDER OF A

!     A(NA,N)  - ON INPUT, A CONTAINS THE REAL SKEW-SYMMETRIC MATRIX.
!                  ONLY THE STRICT LOWER TRIANGLE OF THE MATRIX NEED
!                  BE SUPPLIED.
!                ON OUTPUT, A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
!                  TRANSFORMATIONS USED IN THE REDUCTION IN ITS FULL LOWER
!                  TRIANGLE.  THE STRICT UPPER TRIANGLE OF A IS UNALTERED.

!     E(N)     - OUTPUT ARRAY CONTAINING THE LOWER SUBDIAGONAL ELEMENTS
!                  OF THE TRIDIAGONAL MATRIX IN ITS LAST N-1 POSITIONS.
!                  E(1) IS SET TO ZERO.

!  REQUIRED FUNCTIONS - ABS,SIGN,SQRT

! ****

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: a(:,:)
REAL (dp), INTENT(OUT)  :: e(:)

REAL (dp)  :: f, g, h, scale
INTEGER    :: i, j, k, l, ii

IF (n /= 1) THEN
  
! *** MAIN DO LOOP  I=N STEP -1 UNTIL 2
  
  DO  ii = 2, n
    i = n + 2 - ii
    l = i - 1
    h = 0.
    scale = SUM( ABS(a(i,1:l)) )
    
! *** NORMALIZE ROW
    
    IF (scale == 0.) THEN
      e(i) = 0.
    ELSE
      
! *** COMPUTE ELEMENTS OF U VECTOR
      
      DO  k = 1, l
        a(i,k) = a(i,k) / scale
        h = h + a(i,k) * a(i,k)
      END DO
      
      f = a(i,l)
      g = -SIGN(SQRT(h),f)
      e(i) = scale * g
      h = h - f * g
      a(i,l) = f - g
      IF (l /= 1) THEN
        
! *** COMPUTE ELEMENTS OF A*U/H
        
        DO  j = 1, l
          g = DOT_PRODUCT( a(j,1:j-1), a(i,1:j-1) ) -  &
              DOT_PRODUCT( a(j+1:l,j), a(i,j+1:l) )

          e(j) = g / h
        END DO
        
! *** COMPUTE REDUCED A
        
        DO  j = 2, l
          f = a(i,j)
          g = e(j)
          
          DO  k = 1, j-1
            a(j,k) = a(j,k) + f * e(k) - g * a(i,k)
          END DO
        END DO
      END IF
      
      a(i,1:l) = scale * a(i,1:l)
    END IF
    
    a(i,i) = scale * SQRT(h)
  END DO
END IF

e(1) = 0.
RETURN
END SUBROUTINE trizd



SUBROUTINE imzd(n, e, matz, skew, z, ierr)

! ****

!  FUNCTION  -  COMPUTE THE EIGENVALUES AND OPTIONALLY THE EIGENVECTORS OF A
!                 SYMMETRIC TRIDIAGONAL MATRIX WITH ZERO DIAGONALS OR A SKEW-
!                 SYMMETRIC TRIDIAGONAL MATRIX USING AN IMPLICIT QR-TYPE
!                 ITERATION

!  PARAMETERS

!     N        - INPUT INTEGER SPECIFYING THE ORDER OF THE TRIDIAGONAL MATRIX

!     E(N)     - ON INPUT, ARRAY CONTAINING THE LOWER SUBDIAGONAL ELEMENTS OF
!                  THE TRIDIAGONAL MATRIX IN ITS LAST N-1 POSITIONS.
!                  E(1) IS ARBITRARY.
!                ON OUTPUT, ARRAY CONTAINS THE EIGENVALUES. THE NON-ZERO
!                  EIGENVALUES OCCUR IN PAIRS WITH OPPOSITE SIGNS AND ARE FOUND
!                  IN ADJACENT LOCATIONS IN E. THE EIGENVALUES OF SYMMETRIC
!                  MATRICES ARE REAL AND THE EIGENVALUES OF SKEW-SYMMETRIC
!                  MATRICES ARE PURELY IMAGINARY COMPLEX NUMBERS.  IF AN ERROR
!                  EXIT IS MADE, THE EIGENVALUES ARE CORRECT FOR INDICES
!                  IERR+1, IERR+2...N

!     MATZ     - INPUT LOGICAL VARIABLE SPECIFYING THE EIGENVECTOR
!                  OPTION
!                  = .TRUE.   EIGENVECTORS ARE TO BE COMPUTED
!                  = .FALSE.  EIGENVECTORS ARE NOT TO BE COMPUTED

!     SKEW     - INPUT LOGICAL VARIABLE SPECIFYING TYPE OF INPUT MATRIX
!                  = .TRUE.   INPUT TRIDIAGONAL MATRIX IS SKEW-SYMMETRIC
!                  = .FALSE.  INPUT TRIDIAGONAL MATRIX IS SYMMETRIC WITH ZERO
!                               DIAGONALS
!                  SKEW IS NOT REFERENCED IF MATZ = .FALSE.

!     Z(NZ,N)  - OUTPUT ARRAY CONTAINING THE ORTHOGONAL EIGENVECTORS OF THE
!                  INPUT TRIDIAGONAL MATRIX. EIGENVECTORS CORRESPONDING TO ZERO
!                  EIGENVALUES ARE NORMALIZE TO UNIT 2-NORM (LENGTH) AND THOSE
!                  CORRESPONDING TO NON-ZERO EIGENVALUES HAVE 2-NORM OF SQUARE
!                  ROOT 2.  IF THE J-TH EIGENVALUE IS ZERO OR REAL (I.E. E(J)),
!                  ITS EIGENVECTOR IS FOUND IN THE J-TH COLUMN OF Z.
!                  IF THE J-TH EIGENVALUE IS IMAGINARY (I.E. E(J)*I) WITH
!                  E(J+1) = -E(J), THE REAL PART OF ITS EIGENVECTOR IS FOUND IN
!                  THE J-TH COLUMN OF Z AND ITS IMAGINARY PART FOUND IN THE
!                  (J+1)-TH COLUMN.  IF AN ERROR EXIT IS MADE, Z CONTAINS THE
!                  EIGENVECTORS ASSOCIATED WITH THE STORED EIGENVALUES.
!                  Z IS NOT REFERENCED IF MATZ = .FALSE.

!     IERR     - OUTPUT ERROR CODE
!                  = 0   NORMAL RETURN (ALL EIGENVALUES/VECTORS FOUND)
!                  = J   IF THE J-TH EIGENVALUE HAS NOT BEEN DETERMINED AFTER
!                          30 ITERATIONS

!  REQUIRED FUNCTIONS - ABS,SIGN,SQRT,MOD

! ****

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: e(:)
LOGICAL, INTENT(IN)     :: matz
LOGICAL, INTENT(IN)     :: skew
REAL (dp), INTENT(OUT)  :: z(:,:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp)  :: f, g, q, c, s, r, p, test, tmag
INTEGER    :: i, j, k, l, m, l0, l0m1, m0, ls, im1, jp1, km1, kp1, lm1,  &
              mm1, ip3, ieo, its, ip1
REAL (dp), PARAMETER :: eps = EPSILON( 0.0_dp )

IF (matz) THEN
  
! *** PLACE IDENTITY MATRIX IN Z
  
  DO  i = 1, n
    z(i,1:n) = 0.
    z(i,i) = 1.
  END DO
END IF

ierr = 0
m = n
mm1 = m - 1
e(1) = 0.
its = 0

30 IF (m >= 2) THEN
  m0 = m
  
! *** SEARCH FOR NEXT SUBMATRIX TO SOLVE  (MATRIX SPLITTING)
  
  f = 0.
  DO  i = 1, mm1
    j = m - i
    jp1 = j + 1
    g = ABS(e(jp1))
    tmag = ABS(e(j)) + f
    test = tmag + g
    IF (ABS(test - tmag) < eps*tmag) GO TO 50
    f = g
  END DO
  jp1 = 1
  
  50 l0 = jp1 + 1
  l0m1 = jp1
  IF (l0m1 == m) GO TO 160
  IF (matz) THEN
    IF (skew) THEN
      
! *** PLACE CORRECT SIGN ON IDENTITY DIAGONALS
      
      DO  i = l0m1, m, 4
        z(i,i) = -z(i,i)
        ip3 = i + 3
        IF (ip3 > m) GO TO 70
        z(ip3,ip3) = -z(ip3,ip3)
      END DO
    END IF
  END IF
  
  70 IF (l0 == m) GO TO 170
  ieo = m - l0
  ieo = MOD(ieo,2)
  l = l0
  IF (ieo == 0) GO TO 130
  
! *** FIND ZERO EIGENVALUE OF ODD ORDERED SUBMATRICES
  
  c = 0.
  s = -1.
  DO  i = l0, mm1, 2
    k = mm1 + l0 - i
    kp1 = k + 1
    q = -s * e(kp1)
    e(kp1) = c * e(kp1)
    IF (ABS(e(k)) <= ABS(q)) THEN
      c = e(k) / q
      r = SQRT(c*c+1.)
      e(k) = q * r
      s = 1. / r
      c = c * s
    ELSE
      s = q / e(k)
      r = SQRT(1.+s*s)
      e(k) = e(k) * r
      c = 1. / r
      s = s * c
    END IF
    IF (matz) THEN
      
! *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
      
      km1 = k - 1
      z(km1,m) = -s * z(km1,km1)
      z(km1,km1) = c * z(km1,km1)
      DO  j = kp1, m, 2
        z(j,km1) = s * z(j,m)
        z(j,m) = c * z(j,m)
      END DO
    END IF
    
  END DO
  m = mm1
  mm1 = m - 1
  IF (l0 == m) GO TO 170
  
! *** CHECK FOR CONVERGENCE OR SMALL SUBDIAGONAL ELEMENT
  
  100 DO  i = l0, mm1, 2
    k = mm1 + l0 - i
    l = k + 1
    tmag = ABS(e(l)) + ABS(e(k-1))
    test = tmag + e(k)
    IF (ABS(test - tmag) < eps*tmag) GO TO 120
  END DO
  l = l0
  120 IF (l == m) GO TO 170
  
! *** FORM SHIFT
  
  130 its = its + 1
  IF (its > 30) GO TO 220
  f = e(m-3)
  g = e(m-2)
  c = e(mm1)
  s = e(m)
  p = ((c-f)*(c+f)+(s-g)*(s+g)) / (2.*g*c)
  r = SQRT(p*p+1.)
  q = (g/(p+SIGN(r,p))) - c
  f = e(l)
  lm1 = l - 1
  e(lm1) = ((f-s)*(f+s)+c*q) / f
  
! *** PERFORM ONE IMPLICIT QR ITERATION ON CHOLESKY FACTOR
  
  ls = l0m1
  c = 1.
  s = 1.
  DO  i = l, mm1
    ip1 = i + 1
    im1 = i - 1
    q = s * e(ip1)
    e(ip1) = c * e(ip1)
    IF (ABS(e(im1)) <= ABS(q)) THEN
      c = e(im1) / q
      r = SQRT(c*c+1.)
      e(im1) = q * r
      s = 1. / r
      c = c * s
    ELSE
      s = q / e(im1)
      r = SQRT(1.+s*s)
      e(im1) = e(im1) * r
      c = 1. / r
      s = s * c
    END IF
    f = e(ip1)
    e(ip1) = -s * e(i) + c * f
    e(i) = c * e(i) + s * f
    IF (matz) THEN
      
! *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
      
      DO  j = ls, m0, 2
        f = z(j,ip1)
        z(j,ip1) = -s * z(j,im1) + c * f
        z(j,im1) = c * z(j,im1) + s * f
      END DO
      IF (ls /= l0m1) THEN
        ls = l0m1
      ELSE
        ls = l0
      END IF
    END IF
  END DO
  e(lm1) = 0.
  GO TO 100
  
! *** ITERATION CONVERGED TO ONE ZERO EIGENVALUE
  
  160 e(m) = 0.
  m = mm1
  GO TO 180
  
! *** ITERATION CONVERGED TO EIGENVALUE PAIR
  
  170 e(mm1) = e(m)
  e(m) = -e(m)
  m = m - 2
  
  180 its = 0
  mm1 = m - 1
  IF (m > l0) GO TO 100
  IF (m == l0) GO TO 170
  IF (.NOT.matz) GO TO 30
  IF (skew) GO TO 30
  
! *** COMPUTE EIGENVECTORS FROM ORTHONORMAL COLUMNS OF Z IF NOT SKEW
  
  190 k = m0
  200 IF (e(k) /= 0.) THEN
    km1 = k - 1
    DO  j = l0m1, m0, 2
      z(j,k) = z(j,km1)
      f = z(j+1,k)
      z(j+1,km1) = f
      z(j+1,k) = -f
    END DO
    k = km1
  END IF
  k = k - 1
  IF (k > l0m1) GO TO 200
  IF (ierr /= 0) GO TO 230
  GO TO 30
  
! *** ERROR EXIT
  
  220 ierr = m
  IF (matz) THEN
    IF (.NOT.skew) GO TO 190
  END IF
END IF

230 RETURN
END SUBROUTINE imzd


SUBROUTINE tbakzd(n, a, m, z)

! ****

!  FUNCTION  -  FORMS THE EIGENVECTORS OF A REAL SKEW-SYMMETRIC MATRIX
!                 BY BACK TRANSFORMING THOSE OF THE CORRESPONDING SKEW-
!                 SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRIZD

!  PARAMETERS

!     N        - INPUT INTEGER SPECIFYING THE ORDER OF A

!     A(NA,N)  - INPUT ARRAY CONTAINING INFORMATION ABOUT THE ORTHOGONAL
!                  TRANSFORMATIONS USED IN THE REDUCTION BY  TRIZD  IN
!                  ITS FULL LOWER TRIANGLE

!     M        - INPUT INTEGER SPECIFYING THE NUMBER OF EIGENVECTORS TO
!                  BE BACK TRANSFORMED

!     Z(NZ,M)  - ON INPUT, Z CONTAINS THE REAL AND IMAGINARY (IF COMPLEX)
!                  PARTS OF THE EIGENVECTORS TO BE BACK TRANSFORMED IN ITS
!                  FIRST M COLUMNS
!                ON OUTPUT, Z CONTAINS THE REAL AND IMAGINARY (IF COMPLEX)
!                  PARTS OF THE TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS

! ****

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: a(:,:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: z(:,:)

REAL (dp)  :: h, s
INTEGER    :: i, j, l

IF (m /= 0) THEN
  IF (n /= 1) THEN
    
    DO  i = 2, n
      l = i - 1
      h = a(i,i)
      IF (h /= 0.0_dp) THEN
        
        DO  j = 1, m
          s = DOT_PRODUCT( a(i,1:l), z(1:l,j) )
          
          s = (s/h) / h
          
          z(1:l,j) = z(1:l,j) - s * a(i,1:l)
          
        END DO
      END IF
      
    END DO
  END IF
END IF

RETURN
END SUBROUTINE tbakzd

END MODULE skew_symmetric



PROGRAM driver

! Code converted using TO_F90 by Alan Miller
! Date: 1999-10-29  Time: 18:49:08

! ****

!  FUNCTION  -  THIS IS THE MAIN PROGRAM (DRIVER) FOR ILLUSTRATING THE
!                 USE OF SUBROUTINES TRIZD, IMZD, AND TBAKZD.

!  REFERENCES - EIGENSYSTEM COMPUTATION FOR SKEW-SYMMETRIC MATRICES AND
!                 A CLASS OF SYMMETRIC MATRICES, WARD,R C AND GRAY,L J.
!                 TO APPEAR IN MANUSCRIPT SECTION OF ACM-TOMS.

!               AN ALGORITHM FOR COMPUTING THE EIGENSYSTEM OF SKEW-
!                 SYMMETRIC MATRICES AND A CLASS OF SYMMETRIC MATRICES,
!                 TO APPEAR IN ALGORITHM SECTION OF ACM-TOMS.

!  REQUIRED FUNCTIONS FOR DRIVER AND SUBROUTINES - ABS,SIGN,SQRT,MOD

! ****

USE skew_symmetric

IMPLICIT NONE
REAL (dp)  :: a(6,6), z(6,6), e(6)
REAL (dp)  :: con
INTEGER    :: i, j, n, im1, jp1, ierr
LOGICAL    :: matz, skew

matz = .true.

! *** SET UP SKEW-SYMMETRIC TEST CASE

WRITE (6,5000)
n = 5
skew = .true.

! *** READ AND PRINT TEST MATRIX

WRITE (6,5100)
DO  i = 1, n
  READ (5,5200) a(i,1:n)
  WRITE (6,5200) a(i,1:n)
END DO

! *** COMPUTE EIGENVALUES AND EIGENVECTORS

CALL trizd(n, a, e)
CALL imzd(n, e, matz, skew, z, ierr)
IF (ierr /= 0) WRITE (6,5300) ierr
CALL tbakzd(n, a, n, z)

! *** PRINT EIGENVALUES AND EIGENVECTORS

WRITE (6,5400)
j = 0
20 j = j + 1
WRITE (6,5500)
IF (e(j) /= 0.d0) THEN
  jp1 = j + 1
  WRITE (6,5600) e(j), z(1,j), z(1,jp1)
  WRITE (6,5700) (z(i,j), z(i,jp1), i = 2,n)
  WRITE (6,5500)
  con = -z(1,jp1)
  WRITE (6,5600) e(jp1), z(1,j), con
  DO  i = 2, n
    con = -z(i,jp1)
    WRITE (6,5700) z(i,j), con
  END DO
  j = j + 1
ELSE
  WRITE (6,5800) e(j), z(1,j)
  WRITE (6,5900) z(2:n,j)
END IF
IF (j < n) GO TO 20

! *** SET UP TRIDIAGONAL, SYMMETRIC, ZERO DIAGONAL TEST CASE

WRITE (6,6000)
n = 6
skew = .false.
DO  i = 1, n
  e(i) = 1.
  DO  j = 1, n
    a(i,j) = 0.
  END DO
END DO
DO  i = 2, n
  im1 = i - 1
  a(i,im1) = e(i)
  a(im1,i) = e(i)
END DO

! *** PRINT TEST MATRIX

WRITE (6,6100)
DO  i = 1, n
  WRITE (6,5200) a(i,1:n)
END DO

! *** COMPUTE EIGENVALUES AND EIGENVECTORS

CALL imzd(n, e, matz, skew, z, ierr)
IF (ierr /= 0) WRITE (6,5300) ierr

! *** PRINT EIGENVALUES AND EIGENVECTORS

WRITE (6,6200)
DO  j = 1, n
  WRITE (6,5500)
  WRITE (6,5800) e(j), z(1,j)
  WRITE (6,5900) z(2:n,j)
END DO
STOP

5000 FORMAT ('1', 'EIGENSYSTEM COMPUTATION OF SKEW-SYMMETRIC TEST CASE'//)
5100 FORMAT ('0', '         TEST MATRIX'/)
5200 FORMAT (6F6.0)
5300 FORMAT ('0'/'0', 'IMZD IERR =', i5)
5400 FORMAT ('0'/'0', '  EIGENVALUES', t40, 'EIGENVECTORS')
5500 FORMAT (/)
5600 FORMAT (' ', e15.8, ' * I     ', e15.8, '  +  ', e15.8, ' * I')
5700 FORMAT (t26, e15.8, '  +  ', e15.8, ' * I')
5800 FORMAT (' ', e15.8, t26, e15.8)
5900 FORMAT (t26, e15.8)
6000 FORMAT ('1', 'EIGENSYSTEM COMPUTATION OF TRIDIAGONAL, SYMMETRIC,',  &
             ' ZERO DIAGONAL TEST CASE'//)
6100 FORMAT ('0', t14, 'TEST MATRIX'/)
6200 FORMAT ('0'/'0', '  EIGENVALUES             EIGENVECTORS')
END PROGRAM driver
