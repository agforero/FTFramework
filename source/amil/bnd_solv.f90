MODULE Banded_Solve
! Solve systems of linear banded equations.
! Based upon routines from the Naval Surface Warfare Center Mathematics Library

! If A is an m x n banded matrix, then this module assumes that it is stored
! in banded form.
! e.g. if A has 8 rows and 7 columns but only 4 bands:

!      a11  a12  a13   0    0    0    0
!      a21  a22  a23  a24   0    0    0
!       0   a32  a33  a34  a35   0    0
!       0    0   a43  a44  a45  a46   0
!       0    0    0   a54  a55  a56  a57
!       0    0    0    0   a65  a66  a67
!       0    0    0    0    0   a76  a77
!       0    0    0    0    0    0   a87

! Then the banded form of storage is as:

!       0   a11  a12  a13
!      a21  a22  a23  a24
!      a32  a33  a34  a35
!      a43  a44  a45  a46
!      a54  a55  a56  a57
!      a65  a66  a67   0
!      a76  a77   0    0
!      a87   0    0    0

! The number of lower diagonals in this case is ml = 1.
! The number of upper diagonals is mu = 2.

! N.B. dbslv requires ml extra columns to be added at the right-hand side.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-16  Time: 19:31:30

IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: dbslv

CONTAINS


SUBROUTINE dswap(n, dx, incx, dy, incy)
!
!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: dx(:)
INTEGER, INTENT(IN)        :: incx
REAL (dp), INTENT(IN OUT)  :: dy(:)
INTEGER, INTENT(IN)        :: incy

! Local variables

REAL (dp) :: dtemp
INTEGER   :: i, ix, iy, m, mp1
!
IF (n <= 0) RETURN
IF (incx /= 1 .OR. incy /= 1) THEN
!
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
!
  ix = 1
  iy = 1
  IF (incx < 0) ix = (-n+1) * incx + 1
  IF (incy < 0) iy = (-n+1) * incy + 1
  DO  i = 1, n
    dtemp = dx(ix)
    dx(ix) = dy(iy)
    dy(iy) = dtemp
    ix = ix + incx
    iy = iy + incy
  END DO
  RETURN
END IF
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP
!
m = MOD(n,3)
IF (m /= 0) THEN
  DO  i = 1, m
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
  END DO
  IF (n < 3) RETURN
END IF
mp1 = m + 1
DO  i = mp1, n, 3
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
  dtemp = dx(i+1)
  dx(i+1) = dy(i+1)
  dy(i+1) = dtemp
  dtemp = dx(i+2)
  dx(i+2) = dy(i+2)
  dy(i+2) = dtemp
END DO

RETURN
END SUBROUTINE dswap



SUBROUTINE dbslv(m0, a, ka, n, ml, mu, b, ierr)

! N.B. Argument IWK has been removed.

! ----------------------------------------------------------------------
!  DBSLV EMPLOYS GAUSS ELIMINATION WITH ROW INTERCHANGES TO SOLVE
!  THE NxN BANDED LINEAR SYSTEM AX = B.  THE ARGUMENT M0 SPECIFIES
!  IF DBSLV IS BEING CALLED FOR THE FIRST TIME, OR IF IT IS BEING
!  RECALLED WHERE A IS THE SAME MATRIX BUT B HAS BEEN MODIFIED.
!  ON AN INITIAL CALL TO THE ROUTINE (WHEN M0=0) AN LU DECOMPOSITION
!  OF A IS OBTAINED AND THEN THE EQUATIONS ARE SOLVED.
!  ON SUBSEQUENT CALLS (WHEN M0.NE.0) THE EQUATIONS ARE SOLVED
!  USING THE DECOMPOSITION OBTAINED ON THE INITIAL CALL TO DBSLV.
!
!
!  INPUT ARGUMENTS WHEN M0=0 ---
!
!  A,KA      2-DIMENSIONAL ARRAY OF DIMENSION (KA,M) WHERE
!            KA >= N AND M >= 2*ML + MU + 1.  THE FIRST ML+MU+1
!            COLUMNS CONTAIN THE MATRIX A IN BANDED FORM.
!
!  N         NUMBER OF EQUATIONS AND UNKNOWNS.
!
!  ML        NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!  MU        NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!  B         ARRAY OF N ENTRIES CONTAINING THE RIGHT HAND SIDE DATA.
!
!
!  OUTPUT ARGUMENTS WHEN M0=0 ---
!
!  A         AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!            THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!
!  B         THE SOLUTION OF THE EQUATIONS.
!
!  IWK       ARRAY OF LENGTH N CONTAINING THE PIVOT INDICES.
!
!  IERR      INTEGER SPECIFYING THE STATUS OF THE RESULTS.
!            IERR=0 IF THE SOLUTION OF AX = B IS OBTAINED.
!            OTHERWISE IERR.NE.0.
!
!
!  AFTER AN INITIAL CALL TO DBSLV, THE ROUTINE MAY BE RECALLED WITH M0.NE.0
!  FOR A NEW B.  WHEN M0.NE.0 IT IS ASSUMED THAT A,KA,N,ML,MU,IWK HAVE NOT
!  BEEN MODIFIED.  DBSLV RETRIEVES THE LU DECOMPOSITION WHICH WAS OBTAINED
!  ON THE INITIAL CALL TO DBSLV AND SOLVES THE NEW EQUATIONS AX = B.
!  IN THIS CASE IERR IS NOT REFERENCED.
! ----------------------------------------------------------------------

INTEGER, INTENT(IN)        :: m0
REAL (dp), INTENT(IN OUT)  :: a(:,:)   ! a(ka,*)
INTEGER, INTENT(IN)        :: ka
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ml
INTEGER, INTENT(IN)        :: mu
REAL (dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(OUT)       :: ierr

! Local variables

INTEGER  :: i, iwk(n)

IF (m0 == 0) THEN
!
!     ERROR CHECKING
!
  IF (n <= 0 .OR. n > ka) GO TO 10
  IF (ml < 0 .OR. ml >= n) GO TO 20
  IF (mu < 0 .OR. mu >= n) GO TO 30
!
!     OBTAIN AN LU DECOMPOSITION OF A
!
  CALL dbfa(a, n, ml, mu, iwk, ierr)
  IF (ierr /= 0) RETURN
END IF
!
!     SOLVE THE SYSTEM OF EQUATIONS
!
CALL dbsl(a, n, ml, mu, iwk, b, 0)
RETURN
!
!     ERROR RETURN
!
10 ierr = -1
RETURN

20 ierr = -2
RETURN

30 ierr = -3
RETURN
END SUBROUTINE dbslv



SUBROUTINE dbfa(a, n, ml, mu, ipvt, info)

! N.B. Argument LDA has been removed.

! ----------------------------------------------------------------------
!
!            DBFA FACTORS A REAL BAND MATRIX BY ELIMINATION.
!
!                            --------
!  ON ENTRY
!
!     A       REAL (dp)(LDA, NC)
!             CONTAINS THE MATRIX IN BAND STORAGE.  THE ROWS
!             OF THE ORIGINAL MATRIX ARE STORED IN THE ROWS
!             OF A AND THE DIAGONALS OF THE ORIGINAL MATRIX
!             ARE STORED IN COLUMNS 1 THROUGH ML+MU+1 OF A.
!             NC MUST BE >= 2*ML+MU+1.
!             SEE THE COMMENTS BELOW FOR DETAILS.
!
!     LDA     INTEGER
!             THE LEADING DIMENSION OF THE ARRAY A.  IT IS
!             ASSUMED THAT LDA >= N.
!
!     N       INTEGER
!             THE ORDER OF THE ORIGINAL MATRIX.
!
!     ML      INTEGER
!             NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!             0 <= ML < N .
!
!     MU      INTEGER
!             NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!             0 <= MU < N .
!             MORE EFFICIENT IF ML <= MU .
!
!  ON RETURN
!
!     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE
!             AND THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!             THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!             L IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!             TRIANGULAR MATRICES AND U IS UPPER TRIANGULAR.
!
!     IPVT    INTEGER(N)
!             AN INTEGER VECTOR OF PIVOT INDICES.
!
!     INFO    INTEGER
!             =0  NORMAL VALUE
!             =K  IF  U(K,K) .EQ. 0. THIS IS NOT AN ERROR CONDITION FOR THIS
!                 SUBROUTINE, BUT IT DOES INDICATE THAT DBSL WILL DIVIDE BY
!                 ZERO IF IT IS CALLED.
!
!  BAND STORAGE
!
!        IF A0 IS THE MATRIX THEN THE FOLLOWING CODE WILL STORE
!        A0 IN BAND FORM.
!
!                ML = (BAND WIDTH BELOW THE DIAGONAL)
!                MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                DO 20 I = 1, N
!                   J1 = MAX(1, I-ML)
!                   J2 = MIN(N, I+MU)
!                   DO 10 J = J1, J2
!                      K = J - I + ML + 1
!                      A(I,K) = A0(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!        THIS USES COLUMNS 1 THROUGH  ML + MU + 1  OF A.
!        FURTHERMORE, ML ADDITIONAL COLUMNS ARE NEEDED IN A (STARTING WITH
!        COLUMN  ML+MU+2) FOR ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!        THE TOTAL NUMBER OF COLUMNS NEEDED IN A IS  2*ML+MU+1 .
!
!  EXAMPLE..  IF THE ORIGINAL MATRIX IS
!
!        11 12 13  0  0  0
!        21 22 23 24  0  0
!         0 32 33 34 35  0
!         0  0 43 44 45 46
!         0  0  0 54 55 56
!         0  0  0  0 65 66
!
!   THEN  N = 6, ML = 1, MU = 2, LDA >= 6  AND A SHOULD CONTAIN
!
!           11 12 13  +     , + = USED FOR PIVOTING
!        21 22 23 24  +
!        32 33 34 35  +
!        43 44 45 46  +
!        54 55 56  +  +
!        65 66  +  +  +
!
!  WRITTEN BY E.A.VOORHEES, LOS ALAMOS SCIENTIFIC LABORATORY.
!  MODIFIED BY A.H.MORRIS, NAVAL SURFACE WEAPONS CENTER.
!
!  SUBROUTINES AND FUNCTIONS
!     MIN,DSWAP
! ----------------------------------------------------------------------

REAL (dp), INTENT(IN OUT)  :: a(:,:)   ! a(lda,*)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ml
INTEGER, INTENT(IN)        :: mu
INTEGER, INTENT(OUT)       :: ipvt(:)
INTEGER, INTENT(OUT)       :: info

! Local variables

REAL (dp) :: t
INTEGER   :: i, j, j1, jj, k, l, ll, lm, lm1, lm2, lmk, m, mb, ml1, mp, n1
!
info = 0
IF (ml /= 0) THEN
  m = ml + mu + 1
!
!     SET FILL-IN COLUMNS TO ZERO
!
  DO  j = 1, ml
    jj = m + j
    a(1:n,jj) = 0.d0
  END DO
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
  ml1 = ml + 1
  mb = ml + mu
  n1 = n - 1
!  ldb = lda - 1              ! False dimension
  DO  k = 1, n1
    lm = MIN(n-k,ml)
    lmk = lm + k
    lm1 = lm + 1
    lm2 = ml1 - lm
!
!     SEARCH FOR PIVOT INDEX
!
!     l = -idamax(lm1, a(lmk,lm2), ldb) + lm1 + k
    t = ABS(a(lmk,lm2))
    l = 1
    DO j = 1, lm-1
      IF (ABS(a(lmk-j,lm2+j)) > t) THEN
        t = ABS(a(lmk-j,lm2+j))
        l = j + 1
      END IF
    END DO
    l = lm1 + k - l
    ipvt(k) = l
    mp = MIN(mb, n-k)
!
!     SWAP ROWS IF NECESSARY
!
    ll = ml1 + k - l
    IF (l /= k) CALL dswap(mp+1, a(k,ml1:), 1, a(l,ll:), 1)
!
!     SKIP COLUMN REDUCTION IF PIVOT IS ZERO
!
    IF (a(k,ml1) == 0.d0) THEN
      info = k
    ELSE
!
!     COMPUTE MULTIPLIERS
!
      t = -1.d0 / a(k,ml1)
!      CALL dscal(lm, t, a(lmk,lm2), ldb)
      DO j = 0, lm-1
        a(lmk-j,lm2+j) = t * a(lmk-j,lm2+j)
      END DO
!
!     ROW ELIMINATION WITH COLUMN INDEXING
!
      DO  j = 1, mp
        jj = ml1 + j
        j1 = lm2 + j
!        CALL daxpy(lm, a(k,jj), a(lmk,lm2), ldb, a(lmk,j1), ldb)
        DO i = 0, lm-1
          a(lmk-i,j1+i) = a(lmk-i,j1+i) + a(k,jj) * a(lmk-i,lm2+i)
        END DO
      END DO
    END IF
  END DO
!
  ipvt(n) = n
  IF (a(n,ml1) == 0.d0) info = n
  RETURN
END IF
!
!     CASE WHEN ML = 0
!
DO  k = 1, n
  ipvt(k) = k
  IF (a(k,1) == 0.d0) info = k
END DO

RETURN
END SUBROUTINE dbfa



SUBROUTINE dbsl(a, n, ml, mu, ipvt, b, job)

! N.B. Argument LDA has been removed.

! ----------------------------------------------------------------------
!
!  DBSL SOLVES THE REAL BAND SYSTEM A*X = B OR TRANS(A)*X = B
!  USING THE FACTORS COMPUTED BY DBFA.
!
!                      ----------
!  ON ENTRY
!
!     A       REAL (dp)(LDA, NC)
!             THE OUTPUT FROM DBFA.
!             NC MUST BE >= 2*ML+MU+1 .
!
!     LDA     INTEGER
!             THE LEADING DIMENSION OF THE ARRAY A.
!
!     N       INTEGER
!             THE ORDER OF THE ORIGINAL MATRIX.
!
!     ML      INTEGER
!             NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!     MU      INTEGER
!             NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!     IPVT    INTEGER(N)
!             THE PIVOT VECTOR FROM DBFA.
!
!     B       REAL (dp)(N)
!             THE RIGHT HAND SIDE VECTOR.
!
!     JOB     INTEGER
!             = 0         TO SOLVE  A*X = B .
!             = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                         TRANS(A)  IS THE TRANSPOSE.
!
!  ON RETURN
!
!     B       THE SOLUTION VECTOR  X .
!
!  ERROR CONDITION
!
!     A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A ZERO ON
!     THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY BUT IT IS OFTEN
!     CAUSED BY IMPROPER ARGUMENTS OR IMPROPER SETTING OF LDA.  IT WILL NOT
!     OCCUR IF DBFA AND DBSL ARE CALLED CORRECTLY AND DBFA HAS SET INFO = 0.
!
!  WRITTEN BY E.A. VOORHEES, LOS ALAMOS SCIENTIFIC LABORATORY.
!  MODIFIED BY A.H. MORRIS, NAVAL SURFACE WEAPONS CENTER.
!
!  SUBROUTINES AND FUNCTIONS
!   DAXPY
!
!  FORTRAN  MIN
! ----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: a(:,:)   ! a(lda,*)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ml
INTEGER, INTENT(IN)        :: mu
INTEGER, INTENT(IN)        :: ipvt(:)
REAL (dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(IN)        :: job

! Local variables

REAL (dp) :: t
INTEGER   :: i, ii, k, kb, klm, l, lb, lm, m, ml1, ml2, mlm, nm1
!
m = mu + ml + 1
IF (m /= 1) THEN
  ml1 = ml + 1
  ml2 = ml + 2
  nm1 = n - 1
!  ldb = 1 - lda              ! False dimension
  IF (job == 0) THEN
!
!     JOB = 0 , SOLVE  A * X = B
!       FIRST SOLVE L*Y = B
!
    IF (ml /= 0) THEN
      DO  k = 1, nm1
        lm = MIN(ml,n-k)
        l = ipvt(k)
        t = b(l)
        IF (l /= k) THEN
          b(l) = b(k)
          b(k) = t
        END IF
        klm = k + lm
        mlm = ml1 - lm
!        CALL daxpy(lm, t, a(klm,mlm), ldb, b(k+1), 1)
        DO i = 0, lm-1
          ii = lm-1-i
          b(k+1+i) = b(k+1+i) + t * a(klm-ii,mlm+ii)
        END DO
      END DO
    END IF
!
!       NOW SOLVE  U*X = Y
!
    k = n
    DO  kb = 2, n
      b(k) = b(k) / a(k,ml1)
      lm = MIN(k,m) - 1
      lb = k - lm
      t = -b(k)
!      CALL daxpy(lm, t, a(k-1,ml2), ldb, b(lb), 1)
      DO i = 0, lm-1
        ii = lm-1-i
        b(lb+i) = b(lb+i) + t * a(k-1-ii,ml2+ii)
      END DO
      k = k - 1
    END DO
    b(1) = b(1) / a(1,ml1)
    RETURN
  END IF
!
!     JOB = NONZERO, SOLVE TRANS(A) * X = B
!       FIRST SOLVE  TRANS(U)*Y = B
!
  b(1) = b(1) / a(1,ml1)
  DO  k = 2, n
    lm = MIN(k,m) - 1
    lb = k - lm
!    t = ddot(lm, a(k-1,ml2), ldb, b(lb), 1)
    t = 0.0_dp
    DO i = 0, lm-1
      ii = lm-1-i
      t = t + a(k-1-ii,ml2+ii) * b(lb+i)
    END DO
    b(k) = (b(k)-t) / a(k,ml1)
  END DO
  IF (ml == 0) RETURN
!
!       NOW SOLVE TRANS(L)*X = Y
!
  DO  kb = 1, nm1
    k = n - kb
    lm = MIN(ml,n-k)
    klm = k + lm
    mlm = ml1 - lm
!    b(k) = b(k) + ddot(lm, a(klm,mlm), ldb, b(k+1), 1)
    DO i = 0, lm-1
      ii = lm-1-i
      b(k) = b(k) + a(klm-ii,mlm+ii) * b(k+1+i)
    END DO
    l = ipvt(k)
    IF (l /= k) THEN
      t = b(l)
      b(l) = b(k)
      b(k) = t
    END IF
  END DO
  RETURN
END IF
!
!     CASE WHEN ML = 0 AND MU = 0
!
b(1:n) = b(1:n) / a(1:n,1)

RETURN
END SUBROUTINE dbsl

END MODULE Banded_Solve



PROGRAM t_dbslv
USE Banded_Solve
IMPLICIT NONE

! N.B. Even though A has only 4 bands, the solver needs an extra ML
!      columns to use to store pivots.   Here ML = 1.

REAL (dp)  :: a(8,5), x(8)
INTEGER    :: i, ind, ka = 8, n = 8, ml = 1, mu = 2

! Form the matrix A

a = 0.0_dp
DO i = 1, ka
  IF (i > 1)   a(i,1) = 10*i + i - 1   ! Below the diagonal
               a(i,2) = 10*i + i       ! On the diagonal
  IF (i < n)   a(i,3) = 10*i + i + 1   ! 1st band above the diagonal
  IF (i < n-1) a(i,4) = 10*i + i + 2   ! 2nd band above the diagonal
END DO
x = (/ 12.0_dp, -2.0_dp, 2.0_dp, -2.0_dp, 2.0_dp, -2.0_dp, -77.0_dp, -1.0_dp /)

CALL dbslv(0, a, ka, n, ml, mu, x, ind)
WRITE(*, *) 'IND =', ind
WRITE(*, *) 'Solution:'
WRITE(*, '(8f9.5)') x
WRITE(*, *) 'On exit, matrix A contains:'
DO i = 1, ka
  WRITE(*, '(5f10.5)') a(i,1:5)
END DO

STOP
END PROGRAM t_dbslv
