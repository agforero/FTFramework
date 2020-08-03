SUBROUTINE enum(r, c, n, m, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-02-18  Time: 10:00:55

!    ALGORITHM AS 205 APPL. STATIST. (1984) VOL 33, NO. 3.

!    ENUMERATES ALL R*C CONTINGENCY TABLES WITH GIVEN ROW TOTALS N(I)
!    AND COLUMN TOTALS M(J) AND CALCULATES THE HYPERGEOMETRIC
!    PROBABILITY OF EACH TABLE.

!    FOR TABLES HAVING TWO OR MORE ROW SUMS REPEATED, EQUIVALENT
!    TABLES DIFFERING ONLY BY A ROW PERMUTATION ARE NOT SEPARATELY
!    ENUMERATED.  A REPRESENTATIVE OF EACH EQUIVALENCE CLASS IS ENUMERATED
!    AND THE MULTIPLICITY OF EACH CLASS CALCULATED.

!    FOR EACH TABLE ENUMERATED, SUBROUTINE EVAL IS CALLED TO CARRY OUT
!    CALCULATIONS ON THE TABLE.

IMPLICIT NONE

INTEGER, INTENT(IN OUT)  :: r
INTEGER, INTENT(IN OUT)  :: c
INTEGER, INTENT(IN OUT)  :: n(10)
INTEGER, INTENT(IN OUT)  :: m(10)
INTEGER, INTENT(OUT)     :: ifault

INTERFACE
  SUBROUTINE eval(table, r, c, n, m, prob, mult)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: table(10,10)
    INTEGER, INTENT(IN)  :: r
    INTEGER, INTENT(IN)  :: c
    INTEGER, INTENT(IN)  :: n(10)
    INTEGER, INTENT(IN)  :: m(10)
    REAL, INTENT(IN)     :: prob
    INTEGER, INTENT(IN)  :: mult
  END SUBROUTINE eval
END INTERFACE

INTEGER :: table(10,10), bound(10,10), z(10)
INTEGER :: reps(10), mult(10), repsr, repsc, maxrc, multr, multc
INTEGER :: ntotal, ntot, rm, cm, left, rowbnd, rowsum
INTEGER :: i, j, ii, ij, iip, iim, ip, jj, jjm, jkeep, jm, jnext, jp, k, KEEP
REAL    :: prob(10,10), flogm(201), factlm(201)
REAL    :: prob0
LOGICAL :: rept(10), reptc(10), ieqim

!    LOCAL VARIABLES -
!      TABLE(I,J)   - (I,J)-TH ENTRY OF CURRENT TABLE
!      NTOTAL       - TOTAL OF TABLE ENTRIES
!      BOUND(I,J)   - CURRENT UPPER BOUND ON TABLE(I,J) TO SATISFY
!                     ROW AND COLUMN TOTALS
!      REPT(I)      - LOGICAL = TRUE IF ROW TOTALS N(I), N(I-1) ARE EQUAL
!                             = FALSE OTHERWISE
!      REPS(I)      - NUMBER OF PREVIOUS ROWS EQUAL TO ROW I
!      MULT(I)      - MAXIMUM NUMBER OF EQUIVALENT TABLES
!                     GIVEN FIRST I ROWS
!      Z(J)         - LOWER BOUND ON SUM OF ENTRIES USED BY ALGORITHM C
!      PROB(I,J)    - PARTIAL SUM OF TERMS IN LOG(P)

INTEGER, PARAMETER  :: maxt = 10, nmax = 201, zero = 0.0, one = 1

!      MAXT - MAXIMUM DIMENSION OF TABLE
!      NMAX - MAXIMUM NUMBER OF OBSERVATIONS + 1

!       CHECK INPUT VALUES

ifault = 1
IF (r > MAXT .OR. c > MAXT .OR. r <= 0 .OR. c <= 0) RETURN

ifault = 3
ntotal = 0
DO  i = 1, r
  IF (n(i) <= 0) RETURN
  ntotal = ntotal + n(i)
END DO
ntot = 0
DO  j = 1, c
  IF (m(j) <= 0) RETURN
  ntot = ntot + m(j)
END DO
ifault = 2
IF (ntot /= ntotal) RETURN

ifault = 4
IF (ntotal >= nmax) RETURN

ifault = 0

!       INITIALISE FLOGM(K)=LOG(K-1), FACTLM(K)=LOG(K-1 FACTORIAL)

flogm(1) = zero
factlm(1) = zero
DO  k = 1, ntotal
  flogm(k+1) = LOG(REAL(k))
  factlm(k+1) = factlm(k) + flogm(k+1)
END DO

!       CONSTANTS

rm = r - 1
cm = c - 1

!       SORT ROWS AND COLUMNS INTO ASCENDING ORDER

DO  i = 1, rm
  ip = i + 1
  DO  ii = ip, r
    IF (n(i) > n(ii)) THEN
      KEEP = n(i)
      n(i) = n(ii)
      n(ii) = KEEP
    END IF
  END DO
END DO

DO  j = 1, cm
  jp = j + 1
  DO  jj = jp, c
    IF (m(j) > m(jj)) THEN
      KEEP = m(j)
      m(j) = m(jj)
      m(jj) = KEEP
    END IF
  END DO
END DO

!       CALCULATE MULTIPLICITIES OF ROWS AND COLUMNS

!       REPTC(J) = .TRUE. IF COLUMNS J AND J-1 HAVE THE SAME TOTAL
!       REPT(I)  = .TRUE. IF ROWS I AND I-1 HAVE THE SAME TOTAL

multc = one
repsc = one
reptc(1) = .false.
DO  j = 2, c
  reptc(j) = (m(j) == m(j-1))
  IF (.NOT.reptc(j)) THEN
    repsc = one
  ELSE
    repsc = repsc + one
    multc = multc * repsc
  END IF
END DO

multr = one
repsr = one
rept(1) = .false.
DO  i = 2, r
  rept(i) = (n(i) == n(i-1))
  IF (.NOT.rept(i)) THEN
    repsr = one
  ELSE
    repsr = repsr + one
    multr = multr * repsr
  END IF
END DO

!       IF COLUMN MULTIPLICITY EXCEEDS ROW MULTIPLICITY TRANSPOSE TABLE

IF (multr < multc) THEN
  maxrc = MAX(r,c)
  DO  ij = 1, maxrc
    KEEP = n(ij)
    n(ij) = m(ij)
    m(ij) = KEEP
  END DO
  KEEP = r
  r = c
  c = KEEP
  rm = r - 1
  cm = c - 1
  DO  i = 1, r
    rept(i) = reptc(i)
  END DO
  multr = multc
END IF

!       SET UP INITIAL TABLE

!       MAXIMUM MULTIPLICITY

mult(1) = multr
reps(1) = one

!       CONSTANT TERM IN PROBABILITY

prob0 = -factlm(ntotal+1)
DO  i = 1, r
  ii = n(i)
  prob0 = prob0 + factlm(ii+1)
END DO
DO  j = 1, c
  jj = m(j)
  prob0 = prob0 + factlm(jj+1)
END DO

!       CALCULATE BOUNDS  ON ROW 1

DO  j = 1, c
  bound(1,j) = m(j)
END DO

!       FOR EACH I FIND GREATEST I-TH ROW SATISFYING BOUNDS

DO  i = 1, r
  IF (i /= 1) prob0 = prob(i-1,c)
  left = n(i)
  
!         ELEMENTS OF ROW I
  
  ieqim = rept(i)
  DO  j = 1, cm
    ij = MIN(left,bound(i,j))
    table(i,j) = ij
    IF (j == 1) prob(i,j) = prob0 - factlm(ij+1)
    IF (j /= 1) prob(i,j) = prob(i,j-1) - factlm(ij+1)
    left = left - table(i,j)
    IF (i < r) bound(i+1,j) = bound(i,j) - table(i,j)
    IF (left == 0) GO TO 160
    IF (ieqim) ieqim = table(i,j) == table(i-1,j)
  END DO
  table(i,c) = left
  prob(i,c) = prob(i,cm) - factlm(left+1)
  IF (i < r) bound(i+1,c) = bound(i,c) - left
  GO TO 180

  160 jp = j + 1
  DO  jj = jp, c
    table(i,jj) = 0
    prob(i,jj) = prob(i,jj-1)
    bound(i+1,jj) = bound(i,jj)
  END DO
  180 IF (i /= 1) THEN
    mult(i) = mult(i-1)
    reps(i) = one
    IF (ieqim) THEN
      reps(i) = reps(i-1) + one
      mult(i) = mult(i) / reps(i)
    END IF
  END IF
  
END DO


!       CALL EVAL FOR TABLE 1

CALL eval(table, r, c, n, m, prob(r,c), mult(r))

!       COMMENCE ENUMERATION OF REMAINING TABLES
!       START OF MAIN LOOP

200 i = r
210 i = i - 1

!       IF I = 0 NO MORE TABLES ARE POSSIBLE

IF (i == 0) RETURN

j = cm
left = table(i,c)
rowbnd = bound(i,c)

!       TRY TO DECREASE ELEMENT (I,J)

220 IF (table(i,j) <= 0 .OR. left >= rowbnd) THEN
  
!       ELEMENT (I,J) CANNOT BE DECREASED - TRY (I,J-1)
  
  IF (j == 1) GO TO 210
  left = left + table(i,j)
  rowbnd = rowbnd + bound(i,j)
  j = j - 1
  GO TO 220
END IF

!       DECREASE ELEMENT (I,J)

ij = table(i,j)
prob(i,j) = prob(i,j) + flogm(ij+1)
table(i,j) = table(i,j) - 1
bound(i+1,j) = bound(i+1,j) + 1

!       IF ROW I WAS THE SAME AS ROW I-1 IT IS NO LONGER

IF (reps(i) /= one) THEN
  reps(i) = one
  mult(i) = mult(i-1)
END IF

!       COMPLETE ROW I WITH THE LARGEST POSSIBLE VALUES

ii = i
iip = ii + 1
iim = ii - 1
jnext = j + 1
left = left + 1
GO TO 320

!       FILL UP REMAINING ROWS

230 ii = ii + 1

!       THE LAST ROW IS TREATED SEPARATELY

IF (ii == r) GO TO 370
iip = ii + 1
iim = ii - 1
IF (.NOT.rept(ii)) THEN
  
!       ROW TOTAL N(II) IS NOT A REPEAT - MAKE ROW II AS LARGE AS POSSIBLE
  
  left = n(ii)
  jnext = 1
ELSE
  
!       REPEATED ROW TOTALS
  
!       (I) IF ROW II-1 SATISFIES THE BOUNDS ON ROW II REPEAT IT
  
  DO  j = 1, c
    IF (table(iim,j) > bound(ii,j)) GO TO 250
    ij = table(iim,j)
    table(ii,j) = ij
    bound(iip,j) = bound(ii,j) - table(ii,j)
    IF (j == 1) prob(ii,j) = prob(iim,c) - factlm(ij+1)
    IF (j /= 1) prob(ii,j) = prob(ii,j-1) - factlm(ij+1)
  END DO
  
!       ROW II IS A REPEAT OF ROW II-1
  
  reps(ii) = reps(iim) + one
  mult(ii) = mult(iim) / reps(ii)
  GO TO 230
  
!       ELEMENT J OF ROW II-1 WAS TOO BIG
  
!       CONSTRUCT THE SEQUENCE Z(J) OF LOWER BOUNDS
  
  250 IF (j <= 1) THEN
    
!       IF J=1 THE BOUNDS ARE SATISFIED AUTOMATICALLY
    
    ij = bound(ii,1)
    table(ii,1) = ij
    prob(ii,1) = prob(iim,c) - factlm(ij+1)
    jnext = 2
    left = n(ii) - table(ii,1)
    bound(iip,1) = 0
  ELSE
    z(j) = n(ii)
    jm = j - 1
    IF (j /= c) THEN
      jp = j + 1
      DO  jj = jp, c
        z(j) = z(j) - bound(ii,jj)
      END DO
    END IF
    DO  jjm = 1, jm
      jj = j - jjm
      z(jj) = z(jj+1) - bound(ii,jj+1)
    END DO
    
!       (II) IF THE CUMULATIVE TOTALS OF ROW II-1 ALL EXCEED THE BOUNDS Z(J)
!                      MAKE ELEMENT (II,J) EQUAL TO ITS BOUND
    
    rowsum = 0
    jkeep = 0
    DO  jj = 1, jm
      rowsum = rowsum + table(iim,jj)
      IF (rowsum < z(jj)) GO TO 300
      IF (rowsum > z(jj) .AND. table(iim,jj) > 0) jkeep = jj
    END DO
    table(ii,j) = bound(ii,j)
    bound(iip,j) = 0
    ij = table(ii,j)
    prob(ii,j) = prob(ii,jm) - factlm(ij+1)
    reps(ii) = one
    mult(ii) = mult(iim)
    
!       COMPLETE ROW II WITH THE LARGEST POSSIBLE ELEMENTS
    
    jnext = jp
    left = n(ii)
    DO  jj = 1, j
      left = left - table(ii,jj)
    END DO
    GO TO 320
    
!       (III) THE CUMULATIVE SUMS VIOLATE THE BOUNDS
!       IF NO ELEMENT OF ROW II-1 CAN BE CHANGED TO SATISFY THE BOUNDS
!                 NO SUITABLE ROW II IS POSSIBLE
!       IN THAT CASE GO BACK AND TRY DECREASING ROW II-1
    
    300 IF (jkeep == 0) THEN
      i = ii
      GO TO 210
    END IF
    
!       ELEMENT (II,JKEEP) CAN BE DECREASED
    
    bound(iip,jkeep) = bound(iip,jkeep) + 1
    ij = table(ii,jkeep)
    prob(ii,jkeep) = prob(ii,jkeep) + flogm(ij+1)
    table(ii,jkeep) = table(ii,jkeep) - 1
    
!       COMPLETE THE ROW
    
    jnext = jkeep + 1
    left = n(ii)
    DO  jj = 1, jkeep
      left = left - table(ii,jj)
    END DO
  END IF
END IF

!       ROW II IS COMPLETE UP TO ELEMENT JNEXT-1
!       MAKE THE REMAINING ELEMENTS AS LARGE AS POSSIBLE
!       (THIS SECTION OF CODE IS USED FOR EVERY ROW, REPEATED OR NOT)

320 IF (jnext /= c) THEN
  
  DO  j = jnext, cm
    table(ii,j) = MIN(left,bound(ii,j))
    left = left - table(ii,j)
    bound(iip,j) = bound(ii,j) - table(ii,j)
    ij = table(ii,j)
    IF (j == 1) prob(ii,j) = prob(iim,c) - factlm(ij+1)
    IF (j /= 1) prob(ii,j) = prob(ii,j-1) - factlm(ij+1)
    IF (left == 0) GO TO 340
  END DO
END IF
table(ii,c) = left
prob(ii,c) = prob(ii,cm) - factlm(left+1)
bound(iip,c) = bound(ii,c) - left
GO TO 360

340 jp = j + 1
DO  jj = jp, c
  table(ii,jj) = 0
  prob(ii,jj) = prob(ii,jj-1)
  bound(iip,jj) = bound(ii,jj)
END DO
360 reps(ii) = one
IF (ii > 1) mult(ii) = mult(iim)
GO TO 230

!       THE FINAL ROW

370 IF (.NOT.rept(r)) THEN
  
!       NOT A REPEAT - SET ROW R EQUAL TO ITS BOUNDS
  
  ij = bound(r,1)
  table(r,1) = ij
  prob(r,1) = prob(rm,c) - factlm(ij+1)
  DO  j = 2, c
    ij = bound(r,j)
    table(r,j) = ij
    prob(r,j) = prob(r,j-1) - factlm(ij+1)
  END DO
  mult(r) = mult(rm)
ELSE
  
!       ROW TOTAL R IS A REPEAT - ENSURE THAT IT IS LESS THAN ROW R-1
  
  DO  j = 1, c
    IF (bound(r,j) > table(rm,j)) GO TO 400
    ij = bound(r,j)
    table(r,j) = ij
    IF (j == 1) prob(r,j) = prob(rm,c) - factlm(ij+1)
    IF (j /= 1) prob(r,j) = prob(r,j-1) - factlm(ij+1)
    IF (table(r,j) /= table(rm,j)) GO TO 410
  END DO
  
!       ROW R IS A REPEAT OF ROW R-1
  
  reps(r) = reps(rm) + one
  mult(r) = mult(rm) / reps(r)
  GO TO 430
  
!       IF ROW R WOULD BE BIGGER THAN ROW R-1 GO BACK AND TRY
!          DECREASING ROW R-2
  
  400   i = rm
  GO TO 210
  
!       ROW R IS ALREADY LESS THEN ROW R-1 SO NO MORE CHECKS ARE NEEDED
  
  410   jp = j + 1
  DO  jj = jp, c
    ij = bound(r,jj)
    table(r,jj) = ij
    prob(r,jj) = prob(r,jj-1) - factlm(ij+1)
  END DO
  mult(r) = mult(rm)
END IF

!       THE TABLE IS COMPLETE - CALL SUBROUTINE EVAL

430 CALL eval(table, r, c, n, m, prob(r,c), mult(r))

!       END OF MAIN LOOP

GO TO 200

END SUBROUTINE enum



SUBROUTINE eval(table, r, c, n, m, prob, mult)
IMPLICIT NONE

INTEGER, INTENT(IN)  :: table(10,10)
INTEGER, INTENT(IN)  :: r
INTEGER, INTENT(IN)  :: c
INTEGER, INTENT(IN)  :: n(10)
INTEGER, INTENT(IN)  :: m(10)
REAL, INTENT(IN)     :: prob
INTEGER, INTENT(IN)  :: mult

RETURN

END SUBROUTINE eval
