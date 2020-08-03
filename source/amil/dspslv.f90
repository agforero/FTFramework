MODULE Sparse_Gaussian

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)


CONTAINS

SUBROUTINE dspslv (n, a, ia, ja, b, r, c, MAX, x, itemp, rtemp, ierr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-14  Time: 20:38:15
 
!-----------------------------------------------------------------------
!          SOLUTION OF DOUBLE PRECISION SPARSE EQUATIONS
!-----------------------------------------------------------------------
!  DSPSLV CALLS DNSPIV WHICH USES SPARSE GAUSSIAN ELIMINATION WITH COLUMN
!  INTERCHANGES TO SOLVE THE LINEAR SYSTEM A X = B.  THE ELIMINATION PHASE
!  PERFORMS ROW OPERATIONS ON A AND B TO OBTAIN A UNIT UPPER TRIANGULAR
!  MATRIX U AND A VECTOR Y.  THE SOLUTION PHASE SOLVES U X = Y.

!  INPUT ARGUMENTS---

!  N      INTEGER NUMBER OF EQUATIONS AND UNKNOWNS

!  A      DOUBLE PRECISION ARRAY WITH ONE ENTRY PER NONZERO IN A,
!         CONTAINING THE ACTUAL NONZEROS. (SEE THE MATRIX STORAGE
!         DESCRIPTION BELOW)

!  IA     INTEGER ARRAY OF N+1 ENTRIES CONTAINING ROW POINTERS TO A
!         (SEE MATRIX STORAGE DESCRIPTION BELOW)

!  JA     INTEGER ARRAY WITH ONE ENTRY PER NONZERO IN A, CONTAINING COLUMN
!         NUMBERS OF THE NONZEROS OF A.  (SEE MATRIX STORAGE
!         DESCRIPTION BELOW)

!  B      DOUBLE PRECISION ARRAY OF N ENTRIES CONTAINING THE RIGHT
!         HAND SIDE DATA

!  R      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         ROWS OF A (I.E., THE ELIMINATION ORDER FOR THE EQUATIONS)

!  C      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         COLUMNS OF A.  C IS ALSO AN OUTPUT ARGUMENT

!  MAX    INTEGER NUMBER SPECIFYING MAXIMUM NUMBER OF OFF-DIAGONAL
!         NONZERO ENTRIES OF U WHICH MAY BE STORED

!  ITEMP  INTEGER ARRAY OF 3*N + MAX + 2 ENTRIES, FOR INTERNAL USE

!  RTEMP  DOUBLE PRECISION ARRAY OF N + MAX ENTRIES FOR INTERNAL USE


!  OUTPUT ARGUMENTS---

!  C      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         COLUMNS OF U.  C IS ALSO AN INPUT ARGUMENT

!  X      DOUBLE PRECISION ARRAY OF N ENTRIES CONTAINING THE SOLUTION VECTOR

!  IERR   INTEGER NUMBER WHICH INDICATES ERROR CONDITIONS OR THE ACTUAL
!         NUMBER OF OFF-DIAGONAL ENTRIES IN U (FOR SUCCESSFUL COMPLETION)

!         IERR VALUES ARE---

!         0 < IERR             SUCCESSFUL COMPLETION.  IERR=MAX(1,M)
!                              WHERE M IS THE NUMBER OF OFF-DIAGONAL
!                              NONZERO ENTRIES OF U.

!         IERR = 0             ERROR.  N IS LESS THAN OR EQUAL TO 0

!         -N <= IERR < 0       ERROR.  ROW NUMBER IABS(IERR) OF A IS IS NULL

!         -2*N <= IERR < -N    ERROR.  ROW NUMBER IABS(IERR+N) HAS A
!                              DUPLICATE ENTRY

!         -3*N <= IERR < -2*N  ERROR.  ROW NUMBER IABS(IERR+2*N)
!                              HAS A ZERO PIVOT

!         -4*N <= IERR < -3*N  ERROR.  ROW NUMBER IABS(IERR+3*N)
!                              EXCEEDS STORAGE


!  STORAGE OF SPARSE MATRICES---

!  THE SPARSE MATRIX A IS STORED USING THREE ARRAYS IA, JA, AND A.
!  THE ARRAY A CONTAINS THE NONZEROS OF THE MATRIX ROW-BY-ROW, NOT
!  NECESSARILY IN ORDER OF INCREASING COLUMN NUMBER.  THE ARRAY JA CONTAINS
!  THE COLUMN NUMBERS CORRESPONDING TO THE NONZEROS STORED IN THE ARRAY A
!  (I.E., IF THE NONZERO STORED IN A(K) IS IN COLUMN J, THEN JA(K) = J).
!  THE ARRAY IA CONTAINS POINTERS TO THE ROWS OF NONZEROS/COLUMN INDICES IN
!  THE ARRAY A/JA (I.E., A(IA(I))/JA(IA(I)) IS THE FIRST ENTRY FOR ROW I IN
!  THE ARRAY A/JA).  IA(N+1) IS SET SO THAT IA(N+1) - IA(1) = THE NUMBER OF
!  NONZERO ELEMENTS IN A.
!------------------------

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(*)
INTEGER, INTENT(IN OUT)    :: ia(*)
INTEGER, INTENT(IN OUT)    :: ja(*)
REAL (dp), INTENT(IN OUT)  :: b(n)
INTEGER, INTENT(IN OUT)    :: r(n)
INTEGER, INTENT(IN OUT)    :: c(n)
INTEGER, INTENT(IN OUT)    :: MAX
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(OUT)       :: itemp(*)
REAL (dp), INTENT(IN OUT)  :: rtemp(*)
INTEGER, INTENT(OUT)       :: ierr

INTEGER :: iu, ju, k, l, u, y, p

ierr = 0
IF (n <= 0) RETURN

!  SET INDICES TO DIVIDE TEMPORARY STORAGE FOR DNSPIV

y = 1
u = y + n
p = n + 1
iu = p + n + 1
ju = iu + n + 1

!  COMPUTE THE INVERSE PERMUTATION OF C

DO  k = 1,n
  l = c(k)
  itemp(l) = k
END DO

!  CALL DNSPIV TO PERFORM COMPUTATIONS

CALL dnspiv (n, ia, ja, a, b, MAX, r, c, itemp(1), x, rtemp(y), itemp(p),  &
             itemp(iu), itemp(ju), rtemp(u), ierr)
IF (ierr == 0) ierr = 1
RETURN
END SUBROUTINE dspslv



SUBROUTINE dnspiv(n, ia, ja, a, b, MAX, r, c, ic, x, y, p, iu, ju, u, ierr)
 
!  DNSPIV USES SPARSE GAUSSIAN ELIMINATION WITH COLUMN INTERCHANGES TO SOLVE
!  THE LINEAR SYSTEM A X = B.  THE ELIMINATION PHASE PERFORMS ROW OPERATIONS
!  ON A AND B TO OBTAIN A UNIT UPPER TRIANGULAR MATRIX U AND A VECTOR Y.
!  THE SOLUTION PHASE SOLVES U X = Y.


!  SEE DSPSLV FOR DESCRIPTIONS OF ALL INPUT AND OUTPUT ARGUMENTS
!  OTHER THAN THOSE DESCRIBED BELOW

!  IC  INTEGER ARRAY OF N ENTRIES WHICH IS THE INVERSE OF C
!      (I.E., IC(C(I)) = I). IC IS BOTH AN INPUT AND OUTPUT ARGUMENT.

!  INPUT ARGUMENTS (USED INTERNALLY ONLY)---

!  Y   DOUBLE PRECISION ARRAY OF N ENTRIES USED TO COMPUTE THE UPDATED
!      RIGHT HAND SIDE

!  P   INTEGER ARRAY OF N+1 ENTRIES USED FOR A LINKED LIST.
!      P(N+1) IS THE LIST HEADER, AND THE ENTRY FOLLOWING P(K) IS IN P(P(K)).
!      THUS, P(N+1) IS THE FIRST DATA ITEM, P(P(N+1)) IS THE SECOND, ETC.
!      A POINTER OF N+1 MARKS THE END OF THE LIST

!  IU  INTEGER ARRAY OF N+1 ENTRIES USED FOR ROW POINTERS TO U
!      (SEE MATRIX STORAGE DESCRIPTION BELOW)

!  JU  INTEGER ARRAY OF MAX ENTRIES USED FOR COLUMN NUMBERS OF THE NONZEROS IN
!      THE STRICT UPPER TRIANGLE OF U.  (SEE MATRIX STORAGE DESCRIPTION BELOW)

!  U   DOUBLE PRECISION ARRAY OF MAX ENTRIES USED FOR THE ACTUAL NONZEROS IN
!      THE STRICT UPPER TRIANGLE OF U. (SEE MATRIX STORAGE DESCRIPTION BELOW)


!  STORAGE OF SPARSE MATRICES---

!  THE SPARSE MATRIX A IS STORED USING THREE ARRAYS IA, JA, AND A.
!  THE ARRAY A CONTAINS THE NONZEROS OF THE MATRIX ROW-BY-ROW, NOT
!  NECESSARILY IN ORDER OF INCREASING COLUMN NUMBER.  THE ARRAY JA CONTAINS
!  THE COLUMN NUMBERS CORRESPONDING TO THE NONZEROS STORED IN THE ARRAY A
!  (I.E., IF THE NONZERO STORED IN A(K) IS IN COLUMN J, THEN JA(K) = J).
!  THE ARRAY IA CONTAINS POINTERS TO THE ROWS OF NONZEROS/COLUMN INDICES IN
!  THE ARRAY A/JA (I.E., A(IA(I))/JA(IA(I)) IS THE FIRST ENTRY FOR ROW I IN
!  THE ARRAY A/JA). IA(N+1) IS SET SO THAT IA(N+1) - IA(1) = THE NUMBER OF
!  NONZEROS IN A.  IU, JU, AND U ARE USED IN A SIMILAR WAY TO STORE THE STRICT
!  UPPER TRIANGLE OF U, EXCEPT THAT JU ACTUALLY CONTAINS C(J) INSTEAD OF J


INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ia(*)
INTEGER, INTENT(IN)        :: ja(*)
REAL (dp), INTENT(IN)      :: a(*)
REAL (dp), INTENT(IN)      :: b(n)
INTEGER, INTENT(IN)        :: MAX
INTEGER, INTENT(IN)        :: r(n)
INTEGER, INTENT(IN OUT)    :: c(n)
INTEGER, INTENT(IN OUT)    :: ic(n)
REAL (dp), INTENT(OUT)     :: x(n)
REAL (dp), INTENT(IN OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: p(*)
INTEGER, INTENT(OUT)       :: iu(*)
INTEGER, INTENT(IN OUT)    :: ju(MAX)
REAL (dp), INTENT(IN OUT)  :: u(MAX)
INTEGER, INTENT(OUT)       :: ierr

REAL (dp)  :: dk, lki, xpv, xpvmax, yk

INTEGER    :: ck, i, j, jaj, jmin, jmax, juj, juptr, k, maxc, maxcl, nzcnt,  &
              pk, ppk, pv, v, vi, vj, vk

!  INITIALIZE WORK STORAGE AND POINTERS TO JU

x(1:n) = 0.d0
iu(1) = 1
juptr = 0

!  PERFORM SYMBOLIC AND NUMERIC FACTORIZATION ROW BY ROW
!  VK (VI,VJ) IS THE GRAPH VERTEX FOR ROW K (I,J) OF U

DO  k = 1, n
  
!  INITIALIZE LINKED LIST AND FREE STORAGE FOR THIS ROW
!  THE R(K)-TH ROW OF A BECOMES THE K-TH ROW OF U.
  
  p(n+1) = n + 1
  vk = r(k)
  
!  SET UP ADJACENCY LIST FOR VK, ORDERED IN CURRENT COLUMN ORDER OF U.
!  THE LOOP INDEX GOES DOWNWARD TO EXPLOIT ANY COLUMNS FROM A
!  IN CORRECT RELATIVE ORDER
  
  jmin = ia(vk)
  jmax = ia(vk+1) - 1
  IF (jmin > jmax) GO TO 170
  j = jmax
  20   jaj = ja(j)
  vj = ic(jaj)
  
!  STORE A(K,J) IN WORK VECTOR
  
  x(vj) = a(j)
!  THIS CODE INSERTS VJ INTO ADJACENCY LIST OF VK
  ppk = n + 1
  30   pk = ppk
  ppk = p(pk)
  IF (ppk-vj < 0.0) THEN
    GO TO    30
  ELSE IF (ppk-vj == 0.0) THEN
    GO TO   180
  END IF
  p(vj) = ppk
  p(pk) = vj
  j = j - 1
  IF (j >= jmin) GO TO 20
  
!  THE FOLLOWING CODE COMPUTES THE K-TH ROW OF U
  
  vi = n + 1
  yk = b(vk)
  50   vi = p(vi)
  IF (vi < k) THEN
    
!  VI LT VK -- PROCESS THE L(K,I) ELEMENT AND MERGE THE
!  ADJACENCY OF VI WITH THE ORDERED ADJACENCY OF VK
    
    lki = -x(vi)
    x(vi) = 0.d0
    
!  ADJUST RIGHT HAND SIDE TO REFLECT ELIMINATION
    
    yk = yk + lki * y(vi)
    ppk = vi
    jmin = iu(vi)
    jmax = iu(vi+1) - 1
    IF (jmin > jmax) GO TO 50
    DO  j = jmin, jmax
      juj = ju(j)
      vj = ic(juj)
      
!  IF VJ IS ALREADY IN THE ADJACENCY OF VK, SKIP THE INSERTION
      
      IF (x(vj) == 0.d0) THEN
        
!  INSERT VJ IN ADJACENCY LIST OF VK.
!  RESET PPK TO VI IF WE HAVE PASSED THE CORRECT INSERTION SPOT.
!  (THIS HAPPENS WHEN THE ADJACENCY OF
!  VI IS NOT IN CURRENT COLUMN ORDER DUE TO PIVOTING.)
        
        IF (vj-ppk < 0.0) THEN
          GO TO    60
        ELSE IF (vj-ppk == 0.0) THEN
          GO TO    90
        ELSE
          GO TO    70
        END IF
        60  ppk = vi
        70  pk = ppk
        ppk = p(pk)
        IF (ppk-vj < 0.0) THEN
          GO TO    70
        ELSE IF (ppk-vj == 0.0) THEN
          GO TO    90
        END IF
        p(vj) = ppk
        p(pk) = vj
        ppk = vj
      END IF
      
!  COMPUTE L(K,J) = L(K,J) - L(K,I)*U(I,J) FOR L(K,I) NONZERO
!  COMPUTE U*(K,J) = U*(K,J) - L(K,I)*U(I,J) FOR U(K,J) NONZERO
!  (U*(K,J) = U(K,J)*D(K,K))
      
      90  x(vj) = x(vj) + lki * u(j)
    END DO
    GO TO 50
  END IF
  
!  PIVOT--INTERCHANGE LARGEST ENTRY OF K-TH ROW OF U WITH THE DIAGONAL ENTRY.
  
!  FIND LARGEST ENTRY, COUNTING OFF-DIAGONAL NONZEROS
  
  IF (vi > n) GO TO 190
  xpvmax = ABS(x(vi))
  maxc = vi
  nzcnt = 0
  pv = vi
  110   v = pv
  pv = p(pv)
  IF (pv <= n) THEN
    nzcnt = nzcnt + 1
    xpv = ABS(x(pv))
    IF (xpv <= xpvmax) GO TO 110
    xpvmax = xpv
    maxc = pv
    maxcl = v
    GO TO 110
  END IF
  IF (xpvmax == 0.d0) GO TO 190
  
!  IF VI = K, THEN THERE IS AN ENTRY FOR DIAGONAL
!  WHICH MUST BE DELETED.  OTHERWISE, DELETE THE
!  ENTRY WHICH WILL BECOME THE DIAGONAL ENTRY
  
  IF (vi /= k) THEN
    IF (vi /= maxc) THEN
      p(maxcl) = p(maxc)
      GO TO 120
    END IF
  END IF
  vi = p(vi)
  
!  COMPUTE D(K) = 1/L(K,K) AND PERFORM INTERCHANGE.
  
  120 dk = 1.d0 / x(maxc)
  x(maxc) = x(k)
  i = c(k)
  c(k) = c(maxc)
  c(maxc) = i
  ck = c(k)
  ic(ck) = k
  ic(i) = maxc
  x(k) = 0.d0
  
!  UPDATE RIGHT HAND SIDE.
  
  y(k) = yk * dk
  
!  COMPUTE VALUE FOR IU(K+1) AND CHECK FOR STORAGE OVERFLOW
  
  iu(k+1) = iu(k) + nzcnt
  IF (iu(k+1) > MAX+1) GO TO 200
  
!  MOVE COLUMN INDICES FROM LINKED LIST TO JU.
!  COLUMNS ARE STORED IN CURRENT ORDER WITH ORIGINAL
!  COLUMN NUMBER (C(J)) STORED FOR CURRENT COLUMN J
  
  IF (vi <= n) THEN
    j = vi
    130  juptr = juptr + 1
    ju(juptr) = c(j)
    u(juptr) = x(j) * dk
    x(j) = 0.d0
    j = p(j)
    IF (j <= n) GO TO 130
  END IF
END DO

!  BACKSOLVE U X = Y, AND REORDER X TO CORRESPOND WITH A

k = n
DO  i = 1, n
  yk = y(k)
  jmin = iu(k)
  jmax = iu(k+1) - 1
  IF (jmin <= jmax) THEN
    DO  j = jmin, jmax
      juj = ju(j)
      juj = ic(juj)
      yk = yk - u(j) * y(juj)
    END DO
  END IF
  y(k) = yk
  ck = c(k)
  x(ck) = yk
  k = k - 1
END DO

!  RETURN WITH IERR = NUMBER OF OFF-DIAGONAL NONZEROS IN U

ierr = iu(n+1) - iu(1)
RETURN

!  ERROR RETURNS

!  ROW K OF A IS NULL

170 ierr = -k
RETURN

!  ROW K OF A HAS A DUPLICATE ENTRY

180 ierr = -(n+k)
RETURN

!  ZERO PIVOT IN ROW K

190 ierr = -(2*n+k)
RETURN

!  STORAGE FOR U EXCEEDED ON ROW K

200 ierr = -(3*n+k)
RETURN
END SUBROUTINE dnspiv

END MODULE Sparse_Gaussian
