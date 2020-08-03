MODULE adapt_quad
! Module for adaptive quadrature.
! Adapted from TOMS algorithm 691.
! This version uses the ELF90 subset of Fortran 90.
! Conversion by Alan Miller
! amiller @ bigpond.net.au

CONTAINS

SUBROUTINE qpsrt(limit, last, maxerr, ermax, elist, iord, nrmax)
!     ..................................................................

! 1.   QPSRT
!      ORDERING ROUTINE
!         STANDARD FORTRAN SUBROUTINE
!         REAL VERSION

! 2.   PURPOSE
!         THIS ROUTINE MAINTAINS THE DESCENDING ORDERING IN THE LIST OF THE
!         LOCAL ERROR ESTIMATES RESULTING FROM THE INTERVAL SUBDIVISION
!         PROCESS.  AT EACH CALL TWO ERROR ESTIMATES ARE INSERTED USING THE
!         SEQUENTIAL SEARCH METHOD, TOP-DOWN FOR THE LARGEST ERROR ESTIMATE
!         AND BOTTOM-UP FOR THE SMALLEST ERROR ESTIMATE.

! 3.   CALLING SEQUENCE
!         CALL QPSRT(LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)

!      PARAMETERS (MEANING AT OUTPUT)
!         LIMIT  - INTEGER
!                  MAXIMUM NUMBER OF ERROR ESTIMATES THE LIST CAN CONTAIN

!         LAST   - INTEGER
!                  NUMBER OF ERROR ESTIMATES CURRENTLY IN THE LIST

!         MAXERR - INTEGER
!                  MAXERR POINTS TO THE NRMAX-TH LARGEST ERROR ESTIMATE
!                  CURRENTLY IN THE LIST

!         ERMAX  - REAL
!                  NRMAX-TH LARGEST ERROR ESTIMATE
!                  ERMAX = ELIST(MAXERR)

!         ELIST  - REAL
!                  VECTOR OF DIMENSION LAST CONTAINING THE ERROR ESTIMATES

!         IORD   - INTEGER
!                  VECTOR OF DIMENSION LAST, THE FIRST K ELEMENTS OF
!                  WHICH CONTAIN POINTERS TO THE ERROR ESTIMATES,
!                  SUCH THAT ELIST(IORD(1)), ... , ELIST(IORD(K))
!                  FORM A DECREASING SEQUENCE, WITH K = LAST IF
!                  LAST <= (LIMIT/2+2), AND K = LIMIT+1-LAST OTHERWISE

!         NRMAX  - INTEGER
!                  MAXERR = IORD(NRMAX)

! 4.   NO SUBROUTINES OR FUNCTIONS NEEDED

!     ..................................................................

USE constants_NSWC
IMPLICIT NONE

INTEGER, INTENT(IN)                 :: limit, last
REAL (dp), DIMENSION(:), INTENT(IN) :: elist
INTEGER, INTENT(IN OUT)             :: nrmax
INTEGER, DIMENSION(:), INTENT(OUT)  :: iord
INTEGER, INTENT(OUT)                :: maxerr
REAL (dp), INTENT(OUT)              :: ermax

REAL (dp) :: errmax, errmin
INTEGER   :: i, ibeg, ido, isucc, j, jbnd, jupbn, k

!           CHECK WHETHER THE LIST CONTAINS MORE THAN TWO ERROR ESTIMATES.

!***FIRST EXECUTABLE STATEMENT  QPSRT
IF(last > 2) GO TO 10
iord(1) = 1
iord(2) = 2
GO TO 90

!           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF,
!           DUE TO A DIFFICULT INTEGRAND, SUBDIVISION INCREASED
!           THE ERROR ESTIMATE.   IN THE NORMAL CASE THE INSERT PROCEDURE
!           SHOULD START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.

10 errmax = elist(maxerr)
IF(nrmax == 1) GO TO 30
ido = nrmax-1
DO i = 1, ido
  isucc = iord(nrmax-1)
! ***JUMP OUT OF DO-LOOP
  IF(errmax <= elist(isucc)) EXIT
  iord(nrmax) = isucc
  nrmax = nrmax-1
END DO

!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
!           IN DESCENDING ORDER.  THIS NUMBER DEPENDS ON THE NUMBER OF
!           SUBDIVISIONS STILL ALLOWED.

30 jupbn = last
IF(last > (limit/2+2)) jupbn = limit+3-last
errmin = elist(last)

!           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
!           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).

jbnd = jupbn-1
ibeg = nrmax+1
DO i=ibeg, jbnd
  isucc = iord(i)
! ***JUMP OUT OF DO-LOOP
  IF(errmax >= elist(isucc)) GO TO 60
  iord(i-1) = isucc
END DO
iord(jbnd) = maxerr
iord(jupbn) = last
GO TO 90

!           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.

60 iord(i-1) = maxerr
k = jbnd
DO j=i, jbnd
  isucc = iord(k)
! ***JUMP OUT OF DO-LOOP
  IF(errmin < elist(isucc)) GO TO 80
  iord(k+1) = isucc
  k = k-1
END DO
iord(i) = last
GO TO 90
80 iord(k+1) = last

!           SET MAXERR AND ERMAX.

90 maxerr = iord(nrmax)
ermax = elist(maxerr)
RETURN
END SUBROUTINE qpsrt


SUBROUTINE qelg (n, epstab, result, abserr, res3la, nres, epmach, oflow)
!-----------------------------------------------------------------------

! 1.   PURPOSE
!         THE ROUTINE DETERMINES THE LIMIT OF A GIVEN SEQUENCE OF
!         APPROXIMATIONS, BY MEANS OF THE EPSILON ALGORITHM OF P. WYNN.
!         AN ESTIMATE OF THE ABSOLUTE ERROR IS ALSO GIVEN.
!         THE CONDENSED EPSILON TABLE IS COMPUTED. ONLY THOSE ELEMENTS
!         NEEDED FOR THE COMPUTATION OF THE NEXT DIAGONAL ARE PRESERVED.

! 2.   PARAMETERS
!         N      - INTEGER
!                  EPSTAB(N) CONTAINS THE NEW ELEMENT IN THE
!                  FIRST COLUMN OF THE EPSILON TABLE.

!         EPSTAB - REAL
!                  VECTOR OF DIMENSION 52 CONTAINING THE ELEMENTS OF THE
!                  TWO LOWER DIAGONALS OF THE TRIANGULAR EPSILON TABLE.
!                  THE ELEMENTS ARE NUMBERED STARTING AT THE RIGHT-HAND
!                  CORNER OF THE TRIANGLE.

!         RESULT - REAL
!                  RESULTING APPROXIMATION TO THE INTEGRAL

!         ABSERR - REAL
!                  ESTIMATE OF THE ABSOLUTE ERROR COMPUTED FROM
!                  RESULT AND THE 3 PREVIOUS RESULTS

!         RES3LA - REAL
!                  VECTOR OF DIMENSION 3 CONTAINING THE LAST 3 RESULTS

!         NRES   - INTEGER
!                  NUMBER OF CALLS TO THE ROUTINE
!                  (SHOULD BE ZERO AT FIRST CALL)

!         EPMACH - REAL
!                  THE RELATIVE PRECISION OF THE FLOATING ARITHMETIC
!                  BEING USED.

!         OFLOW  - REAL
!                  THE LARGEST POSITIVE MAGNITUDE.

! 3.   NO SUBROUTINES OR FUNCTIONS USED

!-----------------------------------------------------------------------
USE constants_NSWC
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                 :: n, nres
REAL (dp), INTENT(IN)                   :: epmach, oflow
REAL (dp), INTENT(OUT)                  :: abserr, result
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: epstab, res3la
!---------------------

!      LIST OF MAJOR VARIABLES
!      -----------------------

!      E0     - THE 4 ELEMENTS ON WHICH THE
!      E1       COMPUTATION OF A NEW ELEMENT IN
!      E2       THE EPSILON TABLE IS BASED
!      E3                 E0
!                   E3    E1    NEW
!                         E2
!      NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW DIAGONAL
!      ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!      RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE OF ERROR

!      LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
!      TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
!      DIAGONAL OF THE EPSILON TABLE IS DELETED.

REAL (dp) :: delta1, delta2, delta3, epsinf, error, err1, err2, err3, e0, &
             e1, e1abs, e2, e3, res, ss, tol1, tol2, tol3
INTEGER   :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num

nres = nres + 1
abserr = oflow
result = epstab(n)
IF (n < 3) GO TO 100
limexp = 50
epstab(n + 2) = epstab(n)
newelm = (n - 1)/2
epstab(n) = oflow
num = n
k1 = n
DO i = 1, newelm
  k2 = k1 - 1
  k3 = k1 - 2
  res = epstab(k1 + 2)
  e0 = epstab(k3)
  e1 = epstab(k2)
  e2 = res
  e1abs = ABS(e1)
  delta2 = e2 - e1
  err2 = ABS(delta2)
  tol2 = MAX(ABS(e2),e1abs)*epmach
  delta3 = e1 - e0
  err3 = ABS(delta3)
  tol3 = MAX(e1abs,ABS(e0))*epmach
  IF (err2 > tol2 .OR. err3 > tol3) GO TO 10
  
!      IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
!      ACCURACY, CONVERGENCE IS ASSUMED.
!      RESULT = E2
!      ABSERR = ABS(E1-E0) + ABS(E2-E1)
  
  result = res
  abserr = err2 + err3
! ***JUMP OUT OF DO-LOOP
  GO TO 100
  10 e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1 - e3
  err1 = ABS(delta1)
  tol1 = MAX(e1abs,ABS(e3))*epmach
  
!      IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!      A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
  
  IF (err1 <= tol1 .OR. err2 <= tol2 .OR. err3 <= tol3) GO TO 20
  ss = 1.0D0/delta1 + 1.0D0/delta2 - 1.0D0/delta3
  epsinf = ABS(ss*e1)
  
!      TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND EVENTUALLY
!      OMIT A PART OF THE TABLE ADJUSTING THE VALUE OF N.
  
  IF (epsinf > 0.1D-03) GO TO 30
  20 n = i + i - 1
! ***JUMP OUT OF DO-LOOP
  GO TO 50
  
!      COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST THE VALUE OF RESULT.
  
  30 res = e1 + 1.0D0/ss
  epstab(k1) = res
  k1 = k1 - 2
  error = err2 + ABS(res - e2) + err3
  IF (error > abserr) CYCLE
  abserr = error
  result = res
END DO

!      SHIFT THE TABLE.

50 IF (n == limexp) n = 2*(limexp/2) - 1
ib = 1
IF ((num/2)*2 == num) ib = 2
ie = newelm + 1
DO i = 1, ie
  ib2 = ib + 2
  epstab(ib) = epstab(ib2)
  ib = ib2
END DO
IF (num == n) GO TO 80
indx = num - n + 1
DO i = 1, n
  epstab(i) = epstab(indx)
  indx = indx + 1
END DO
80 IF (nres >= 4) GO TO 90
res3la(nres) = result
abserr = oflow
GO TO 100

!      COMPUTE ERROR ESTIMATE

90 abserr = ABS(result - res3la(3)) + ABS(result - res3la(2)) +  &
            ABS(result - res3la(1))
res3la(1) = res3la(2)
res3la(2) = res3la(3)
res3la(3) = result
100 abserr = MAX(abserr, 5.0D0*epmach*ABS(result))
RETURN
END SUBROUTINE qelg


SUBROUTINE qxgs (f, a, b, epsabs, epsrel, result, abserr, ier, limit, last)

!     THE ROUTINE CALCULATES AN APPROXIMATION RESULT TO A GIVEN
!     DEFINITE INTEGRAL  I = INTEGRAL OF F OVER (A,B),
!     HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
!     ABS(I-RESULT) <= MAX(EPSABS, EPSREL*ABS(I)).

! PARAMETERS
!  ON ENTRY
!     F      - REAL
!              FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!              FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!              DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!     A      - REAL
!              LOWER LIMIT OF INTEGRATION

!     B      - REAL
!              UPPER LIMIT OF INTEGRATION

!     EPSABS - REAL
!              ABSOLUTE ACCURACY REQUESTED

!     EPSREL - REAL
!              RELATIVE ACCURACY REQUESTED

!  ON RETURN
!     RESULT - REAL
!              APPROXIMATION TO THE INTEGRAL

!     ABSERR - REAL
!              ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!              WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

!     IER    - INTEGER
!              IER = 0 NORMAL AND RELIABLE TERMINATION OF THE ROUTINE.
!                      IT IS ASSUMED THAT THE REQUESTED ACCURACY HAS
!                      BEEN ACHIEVED.
!              IER > 0 ABNORMAL TERMINATION OF THE ROUTINE
!                      THE ESTIMATES FOR INTEGRAL AND ERROR ARE
!                      LESS RELIABLE. IT IS ASSUMED THAT THE
!                      REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
!     ERROR MESSAGES
!              IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED HAS BEEN
!                      ACHIEVED.  ONE CAN ALLOW MORE SUB-DIVISIONS BY
!                      INCREASING THE VALUE OF LIMIT (AND TAKING THE ACCORDING
!                      DIMENSION ADJUSTMENTS INTO ACCOUNT.  HOWEVER, IF THIS
!                      YIELDS NO IMPROVEMENT IT IS ADVISED TO ANALYZE THE
!                      INTEGRAND IN ORDER TO DETERMINE THE INTEGRATION
!                      DIFFICULTIES.  IF THE POSITION OF A LOCAL DIFFICULTY
!                      CAN BE DETERMINED (E.G. SINGULARITY, DISCONTINUITY
!                      WITHIN THE INTERVAL) ONE WILL PROBABLY GAIN FROM
!                      SPLITTING UP THE INTERVAL AT THIS POINT AND CALLING THE
!                      INTEGRATOR ON THE SUBRANGES.  IF POSSIBLE, AN
!                      APPROPRIATE SPECIAL-PURPOSE INTEGRATOR SHOULD BE USED,
!                      WHICH IS DESIGNED FOR HANDLING THE TYPE OF DIFFICULTY
!                      INVOLVED.
!                  = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETECTED,
!                      WHICH PREVENTS THE REQUESTED TOLERANCE FROM BEING
!                      ACHIEVED.
!                      THE ERROR MAY BE UNDER-ESTIMATED.
!                  = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR
!                      OCCURS AT SOME POINTS OF THE INTEGRATION INTERVAL.
!                  = 4 THE ALGORITHM DOES NOT CONVERGE.
!                      ROUNDOFF ERROR IS DETECTED IN THE EXTRAPOLATION TABLE.
!                      IT IS PRESUMED THAT THE REQUESTED TOLERANCE CANNOT BE
!                      ACHIEVED, AND THAT THE RETURNED RESULT IS THE BEST
!                      WHICH CAN BE OBTAINED.
!                  = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY CONVERGENT.
!                      IT MUST BE NOTED THAT DIVERGENCE CAN OCCUR WITH ANY
!                      OTHER VALUE OF IER.
!                  = 6 THE INPUT IS INVALID BECAUSE EPSABS OR EPSREL IS
!                      NEGATIVE, LIMIT < 1, LENW < 46*LIMIT, OR LENIW < 3*LIMIT.
!                      RESULT, ABSERR, LAST ARE SET TO ZERO.
!                      EXCEPT WHEN LIMIT OR LENW OR LENIW IS INVALID, IWORK(1),
!                      WORK(LIMIT*2+1) AND WORK(LIMIT*3+1) ARE SET TO ZERO,
!                      WORK(1) IS SET TO A, AND WORK(LIMIT+1) TO B.

!  DIMENSIONING PARAMETERS
!     LIMIT - INTEGER
!             LIMIT DETERMINES THE MAXIMUM NUMBER OF SUBINTERVALS IN THE
!             PARTITION OF THE GIVEN INTEGRATION INTERVAL (A,B), LIMIT >= 1.
!             IF LIMIT < 1, THE ROUTINE WILL END WITH IER = 6.

!     LAST  - INTEGER
!             ON RETURN, LAST EQUALS THE NUMBER OF SUBINTERVALS PRODUCED
!             IN THE SUBDIVISION PROCESS, WHICH DETERMINES THE NUMBER OF
!             SIGNIFICANT ELEMENTS ACTUALLY IN THE WORK ARRAYS.

USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)  :: a, b, epsabs, epsrel
REAL (dp), INTENT(OUT) :: result, abserr
INTEGER, INTENT(IN)    :: limit
INTEGER, INTENT(OUT)   :: ier, last

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

!         CHECK VALIDITY OF LIMIT

ier = 6
last = 0
result = 0.0D0
abserr = 0.0D0
IF (limit < 1) RETURN

!         PREPARE CALL FOR QXGSE.

CALL qxgse(f, a, b, epsabs, epsrel, limit, result, abserr, ier, last)

RETURN
END SUBROUTINE qxgs


SUBROUTINE qxgse(f, a, b, epsabs, epsrel, limit, result, abserr, ier, last)

!       THE ROUTINE CALCULATES AN APPROXIMATION RESULT TO A
!       DEFINITE INTEGRAL I = INTEGRAL OF F OVER (A,B),
!       HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
!       ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

!   PARAMETERS
!    ON ENTRY
!       F      - REAL
!                FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!                FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!                DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!       A      - REAL
!                LOWER LIMIT OF INTEGRATION

!       B      - REAL
!                UPPER LIMIT OF INTEGRATION

!       EPSABS - REAL
!                ABSOLUTE ACCURACY REQUESTED

!       EPSREL - REAL
!                RELATIVE ACCURACY REQUESTED

!       LIMIT  - INTEGER
!                GIVES AN UPPERBOUND ON THE NUMBER OF SUBINTERVALS
!                IN THE PARTITION OF (A,B)

!    ON RETURN
!       RESULT - REAL
!                APPROXIMATION TO THE INTEGRAL

!       ABSERR - REAL
!                ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!                WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

!       IER    - INTEGER
!                IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
!                        ROUTINE. IT IS ASSUMED THAT THE REQUESTED
!                        ACCURACY HAS BEEN ACHIEVED.
!                IER > 0 ABNORMAL TERMINATION OF THE ROUTINE
!                        THE ESTIMATES FOR INTEGRAL AND ERROR ARE
!                        LESS RELIABLE. IT IS ASSUMED THAT THE
!                        REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
!       ERROR MESSAGES
!                    = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
!                        HAS BEEN ACHIEVED. ONE CAN ALLOW MORE SUB-
!                        DIVISIONS BY INCREASING THE VALUE OF LIMIT
!                        (AND TAKING THE ACCORDING DIMENSION
!                        ADJUSTMENTS INTO ACCOUNT). HOWEVER, IF
!                        THIS YIELDS NO IMPROVEMENT IT IS ADVISED
!                        TO ANALYZE THE INTEGRAND IN ORDER TO
!                        DETERMINE THE INTEGRATION DIFFICULTIES. IF
!                        THE POSITION OF A LOCAL DIFFICULTY CAN BE
!                        DETERMINED (E.G. SINGULARITY,
!                        DISCONTINUITY WITHIN THE INTERVAL) ONE
!                        WILL PROBABLY GAIN FROM SPLITTING UP THE
!                        INTERVAL AT THIS POINT AND CALLING THE
!                        INTEGRATOR ON THE SUBRANGES. IF POSSIBLE,
!                        AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
!                        SHOULD BE USED, WHICH IS DESIGNED FOR
!                        HANDLING THE TYPE OF DIFFICULTY INVOLVED.
!                    = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETEC-
!                        TED, WHICH PREVENTS THE REQUESTED
!                        TOLERANCE FROM BEING ACHIEVED.
!                        THE ERROR MAY BE UNDER-ESTIMATED.
!                    = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS AT
!                        SOME POINTS OF THE INTEGRATION INTERVAL.
!                    = 4 THE ALGORITHM DOES NOT CONVERGE.
!                        ROUNDOFF ERROR IS DETECTED IN THE
!                        EXTRAPOLATION TABLE.
!                        IT IS PRESUMED THAT THE REQUESTED TOLERANCE
!                        CANNOT BE ACHIEVED, AND THAT THE RETURNED RESULT
!                        IS THE BEST WHICH CAN BE OBTAINED.
!                    = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY
!                        CONVERGENT.   IT MUST BE NOTED THAT DIVERGENCE
!                        CAN OCCUR WITH ANY OTHER VALUE OF IER.
!                    = 6 THE INPUT IS INVALID BECAUSE EPSABS OR
!                        EPSREL IS NEGATIVE. RESULT, ABSERR,
!                        LAST, RLIST(1), IORD(1), AND ELIST(1)
!                        ARE SET TO ZERO. ALIST(1) AND BLIST(1)
!                        ARE SET TO A AND B RESPECTIVELY.

!       ALIST  - REAL
!                VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                 LAST  ELEMENTS OF WHICH ARE THE LEFT END POINTS
!                OF THE SUBINTERVALS IN THE PARTITION OF THE
!                GIVEN INTEGRATION RANGE (A,B)

!       BLIST  - REAL
!                VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                 LAST  ELEMENTS OF WHICH ARE THE RIGHT END POINTS
!                OF THE SUBINTERVALS IN THE PARTITION OF THE GIVEN
!                INTEGRATION RANGE (A,B)

!       RLIST  - REAL
!                VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST `LAST'
!                ELEMENTS OF WHICH ARE THE INTEGRAL APPROXIMATIONS ON
!                THE SUBINTERVALS

!       ELIST  - REAL
!                VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                 LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
!                ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS

!       IORD   - INTEGER
!                VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST K ELEMENTS
!                OF WHICH ARE POINTERS TO THE ERROR ESTIMATES OVER THE
!                SUBINTERVALS, SUCH THAT ELIST(IORD(1)), ..., ELIST(IORD(K))
!                FORM A DECREASING SEQUENCE, WITH K = LAST IF
!                LAST <= (LIMIT/2+2), AND K = LIMIT+1-LAST OTHERWISE

!       LAST   - INTEGER
!                NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN THE
!                SUBDIVISION PROCESS

!       VALP   - REAL
!       VALN     ARRAYS OF DIMENSION AT LEAST (21,LIMIT) USED TO
!                SAVE THE FUNCTIONAL VALUES

!       LP     - INTEGER
!       LN       VECTORS OF DIMENSION AT LEAST LIMIT, USED TO
!                STORE THE ACTUAL NUMBER OF FUNCTIONAL VALUES
!                SAVED IN THE CORRESPONDING COLUMN OF VALP,VALN

!***ROUTINES CALLED  F, SPMPAR, QELG, QXLQM, QPSRT, QXRRD, QXCPY

USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: a, b, epsabs, epsrel
REAL (dp), INTENT(OUT)  :: result, abserr
INTEGER, INTENT(IN)     :: limit
INTEGER, INTENT(OUT)    :: ier, last

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

REAL (dp) :: abseps, alist(limit), area, area1, area12, area2, a1, a2, b1, &
             b2, blist(limit), correc, defab1, defab2, dres, elist(limit), &
             epmach, erlarg, erlast, errbnd, errmax, error1, error2,       &
             erro12, errsum, ertest, oflow, rerr, resabs, reseps, res3la(3), &
             rlist(limit), rlist2(52), small, t, uflow, valp(21,limit),    &
             valn(21,limit), vp1(21), vp2(21), vn1(21), vn2(21), defabs
INTEGER   :: id, ierro, iord(limit), iroff1, iroff2, iroff3, jupbnd, k,  &
             ksgn, lp(limit), ln(limit), ktmin, maxerr, nres, nrmax, numrl2, &
             lp1, lp2, ln1, ln2
LOGICAL   :: extrap, noext

!            MACHINE DEPENDENT CONSTANTS
!            ---------------------------

!          EPMACH IS THE LARGEST RELATIVE SPACING.
!          UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!          OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

epmach = dpmpar(1)
uflow = dpmpar(2)
oflow = dpmpar(3)

!            TEST ON VALIDITY OF PARAMETERS
!            ------------------------------
last = 0
result = 0.0D0
abserr = 0.0D0
alist(1) = a
blist(1) = b
rlist(1) = 0.0D0
elist(1) = 0.0D0
ier = 6
IF (epsabs < 0.0D0 .OR. epsrel < 0.0D0) GO TO 999
ier = 0
rerr = MAX(epsrel, 50.0D0*epmach)

!      FIRST APPROXIMATION TO THE INTEGRAL
!      -----------------------------------

ierro = 0
lp(1) = 1
ln(1) = 1
valp(1,1) = f((a + b)*0.5D0)
valn(1,1) = valp(1,1)
CALL qxlqm(f, a, b, result, abserr, resabs, defabs, valp(1:,1), valn(1:,1), &
           lp(1), ln(1), 2, epmach, uflow, oflow)

!      TEST ON ACCURACY.

dres = ABS(result)
errbnd = MAX(epsabs,rerr*dres)
last = 1
rlist(1) = result
elist(1) = abserr
iord(1) = 1
IF (abserr <= 100.0D0*epmach*defabs .AND. abserr > errbnd) ier = 2
IF (limit == 1) ier = 1
IF (ier /= 0 .OR. (abserr <= errbnd .AND. abserr /= resabs) .OR.  &
abserr == 0.0D0) GO TO 999

!      INITIALIZATION
!      --------------

rlist2(1) = result
errmax = abserr
maxerr = 1
area = result
errsum = abserr
abserr = oflow
nrmax = 1
nres = 0
numrl2 = 2
ktmin = 0
extrap = .false.
noext = .false.
iroff1 = 0
iroff2 = 0
iroff3 = 0
ksgn = -1
IF (dres >= (1.0D0 - 50.0D0*epmach)*defabs) ksgn = 1
t = 1.0D0 + 100.0D0*epmach

!      MAIN DO-LOOP
!      ------------

DO last = 2, limit
  
!      BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR ESTIMATE.
  
  a1 = alist(maxerr)
  b1 = 0.5D0*(alist(maxerr) + blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  erlast = errmax
  CALL qxrrd(f, valn(1:,maxerr), ln(maxerr), b1, a1, vn1, vp1, ln1, lp1)
  CALL qxlqm(f, a1, b1, area1, error1, resabs, defab1, vp1, vn1, lp1, ln1,  &
             2, epmach, uflow, oflow)
  CALL qxrrd(f, valp(1:,maxerr), lp(maxerr), a2, b2, vp2, vn2, lp2, ln2)
  CALL qxlqm(f, a2, b2, area2, error2, resabs, defab2, vp2, vn2, lp2, ln2,  &
             2, epmach, uflow, oflow)
  
!      IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!      AND ERROR AND TEST FOR ACCURACY.
  
  area12 = area1 + area2
  erro12 = error1 + error2
  errsum = errsum + erro12 - errmax
  area = area + area12 - rlist(maxerr)
  IF (defab1 == error1 .OR. defab2 == error2) GO TO 15
  IF (ABS(rlist(maxerr) - area12) > 0.1D-04*ABS(area12)  &
        .OR. erro12 < 0.99D0*errmax) GO TO 10
  IF (extrap) iroff2 = iroff2 + 1
  IF (.NOT.extrap) iroff1 = iroff1 + 1
  10 IF (last > 10 .AND. erro12 > errmax) iroff3 = iroff3 + 1
  15 rlist(maxerr) = area1
  rlist(last) = area2
  errbnd = MAX(epsabs,rerr*ABS(area))
  
!      TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
  
  IF (iroff1 + iroff2 >= 10 .OR. iroff3 >= 20) ier = 2
  IF (iroff2 >= 5) ierro = 3
  
!      SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS EQUALS LIMIT.
  
  IF (last == limit) ier = 1
  
!      SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!      AT A POINT OF THE INTEGRATION RANGE.
  
  IF (MAX(ABS(a1),ABS(b2)) <= t*(ABS(a2) + 1.d+03*uflow)) ier = 4
  
!      APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
  
  IF (error2 > error1) GO TO 20
  alist(last) = a2
  blist(maxerr) = b1
  blist(last) = b2
  elist(maxerr) = error1
  elist(last) = error2
  CALL qxcpy(valp(1:,maxerr), vp1, lp1)
  lp(maxerr) = lp1
  CALL qxcpy(valn(1:,maxerr), vn1, ln1)
  ln(maxerr) = ln1
  CALL qxcpy(valp(1:,last), vp2, lp2)
  lp(last) = lp2
  CALL qxcpy(valn(1:,last), vn2, ln2)
  ln(last) = ln2
  GO TO 30
  
  20 alist(maxerr) = a2
  alist(last) = a1
  blist(last) = b1
  rlist(maxerr) = area2
  rlist(last) = area1
  elist(maxerr) = error2
  elist(last) = error1
  CALL qxcpy(valp(1:,maxerr), vp2, lp2)
  lp(maxerr) = lp2
  CALL qxcpy(valn(1:,maxerr), vn2, ln2)
  ln(maxerr) = ln2
  CALL qxcpy(valp(1:,last), vp1, lp1)
  lp(last) = lp1
  CALL qxcpy(valn(1:,last), vn1, ln1)
  ln(last) = ln1
  
!      CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
!      IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!      WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
  
  30 CALL qpsrt(limit, last, maxerr, errmax, elist, iord, nrmax)
! ***JUMP OUT OF DO-LOOP
  IF(errsum <= errbnd) GO TO 115
! ***JUMP OUT OF DO-LOOP
  IF (ier /= 0) GO TO 100
  IF (last == 2) GO TO 80
  IF (noext) CYCLE
  erlarg = erlarg - erlast
  IF (ABS(b1 - a1) > small) erlarg = erlarg + erro12
  IF (extrap) GO TO 40
  
!      TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE SMALLEST INTERVAL.
  
  IF (ABS(blist(maxerr) - alist(maxerr)) > small) CYCLE
  extrap = .true.
  nrmax = 2
  
!      THE BOUND 0.3*ERTEST HAS BEEN INTRODUCED TO PERFORM A
!      MORE CAUTIOUS EXTRAPOLATION THAN IN THE ORIGINAL DQAGSE ROUTINE
  
  40 IF (ierro == 3 .OR. erlarg <= 0.3D0*ertest) GO TO 60
  
!      THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!      BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
!      LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
  
  id = nrmax
  jupbnd = last
  IF (last > (2 + limit/2)) jupbnd = limit + 3 - last
  DO k = id, jupbnd
    maxerr = iord(nrmax)
    errmax = elist(maxerr)
! ***JUMP OUT OF DO-LOOP
    IF(ABS(blist(maxerr) - alist(maxerr)) > small) CYCLE
    nrmax = nrmax + 1
  END DO
  
!      PERFORM EXTRAPOLATION.
  
  60 numrl2 = numrl2 + 1
  rlist2(numrl2) = area
  CALL qelg (numrl2, rlist2, reseps, abseps, res3la, nres, epmach, oflow)
  ktmin = ktmin + 1
  IF (ktmin > 5 .AND. abserr < 0.1D-02*errsum) ier = 5
  IF (abseps >= abserr) GO TO 70
  ktmin = 0
  abserr = abseps
  result = reseps
  correc = erlarg
  ertest = MAX(epsabs,rerr*ABS(reseps))
! ***JUMP OUT OF DO-LOOP
  IF (abserr <= ertest) GO TO 100
  
!      PREPARE BISECTION OF THE SMALLEST INTERVAL.
  
  70 IF (numrl2 == 1) noext = .true.
  IF (ier == 5) GO TO 100
  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .false.
  small = small*0.5D0
  erlarg = errsum
  CYCLE
  80 small = ABS(b - a)*0.375D0
  erlarg = errsum
  ertest = errbnd
  rlist2(2) = area
END DO

!      SET FINAL RESULT AND ERROR ESTIMATE.
!      ------------------------------------

100 IF (abserr == oflow) GO TO 115
IF (ier + ierro == 0) GO TO 110
IF (ierro == 3) abserr = abserr + correc
IF (ier == 0) ier = 3
IF (result /= 0.0D0 .AND. area /= 0.0D0) GO TO 105
IF (abserr > errsum) GO TO 115
IF (area == 0.0D0) GO TO 130
GO TO 110
105 IF (abserr/ABS(result) > errsum/ABS(area)) GO TO 115

!      TEST ON DIVERGENCE.

110 IF(ksgn == (-1) .AND. MAX(ABS(result),ABS(area)) <=  &
defabs*0.1D-01) GO TO 130
IF(0.1D-01 > (result/area) .OR. (result/area) > 0.1D+03  &
   .OR. errsum > ABS(area)) ier = 6
GO TO 130

!      COMPUTE GLOBAL INTEGRAL SUM.

115 result = SUM( rlist(1:last) )
abserr = errsum
130 IF (ier > 2) ier = ier - 1
999 RETURN
END SUBROUTINE qxgse


SUBROUTINE qxcpy (a, b, l)

!  TO COPY THE REAL VECTOR B OF LENGTH L   I N T O
!          THE REAL VECTOR A OF LENGTH L

USE constants_NSWC
IMPLICIT NONE

INTEGER, INTENT(IN)   :: l
REAL (dp), DIMENSION(:), INTENT(IN)  :: b
REAL (dp), DIMENSION(:), INTENT(OUT) :: a

a(1:l) = b(1:l)
RETURN
END SUBROUTINE qxcpy


SUBROUTINE qxlqm (f, a, b, result, abserr, resabs, resasc, vr, vs, lr, ls,  &
                  key, epmach, uflow, oflow)

!    TO COMPUTE I = INTEGRAL OF F OVER (A, B), WITH ERROR ESTIMATE
!               J = INTEGRAL OF ABS(F) OVER (A,B)

!   PARAMETERS
!    ON ENTRY
!      F      - REAL
!               FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!               FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!               DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!      A      - REAL
!               LOWER LIMIT OF INTEGRATION

!      B      - REAL
!               UPPER LIMIT OF INTEGRATION

!      VR     - REAL
!               VECTOR OF LENGTH LR CONTAINING THE
!               SAVED  FUNCTIONAL VALUES OF POSITIVE ABSCISSAS

!      VS     - REAL
!               VECTOR OF LENGTH LS CONTAINING THE
!               SAVED  FUNCTIONAL VALUES OF NEGATIVE ABSCISSAS

!      LR     - INTEGER
!      LS       NUMBER OF ELEMENTS IN VR,VS RESPECTIVELY

!    KEY    - INTEGER
!             KEY FOR CHOICE OF LOCAL INTEGRATION RULE
!             RMS FORMULAS ARE USED WITH
!              13 - 19               POINTS IF KEY < 1,
!              13 - 19 - (27)        POINTS IF KEY = 1,
!              13 - 19 - (27) - (41) POINTS IF KEY = 2,
!                   19 -  27  - (41) POINTS IF KEY = 3,
!                         27  -  41  POINTS IF KEY > 3.

!             (RULES) USED IF THE FUNCTION APPEARS REGULAR ENOUGH

!      EPMACH - REAL
!               THE RELATIVE PRECISION OF THE FLOATING
!               ARITHMETIC BEING USED.

!      UFLOW  - REAL
!               THE SMALLEST POSITIVE MAGNITUDE.

!      OFLOW  - REAL
!               THE LARGEST POSITIVE MAGNITUDE.

!    ON RETURN
!      RESULT - REAL
!               APPROXIMATION TO THE INTEGRAL I

!      ABSERR - REAL
!               ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!               WHICH SHOULD NOT EXCEED ABS(I-RESULT)

!      RESABS - REAL
!               APPROXIMATION TO THE INTEGRAL J

!      RESASC - REAL
!               APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A)) OVER (A,B)

!      VR     - REAL
!               VECTOR OF LENGTH LR CONTAINING THE
!               SAVED  FUNCTIONAL VALUES OF POSITIVE ABSCISSAS

!      VS     - REAL
!               VECTOR OF LENGTH LS CONTAINING THE
!               SAVED  FUNCTIONAL VALUES OF NEGATIVE ABSCISSAS

!      LR     - INTEGER
!      LS       NUMBER OF ELEMENTS IN VR,VS RESPECTIVELY

!***ROUTINES CALLED  QXRUL

USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)                   :: a, b, epmach, oflow, uflow
INTEGER, INTENT(IN)                     :: key
REAL (dp), INTENT(OUT)                  :: result, abserr, resabs, resasc
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: vr, vs
INTEGER, INTENT(IN OUT)                 :: lr, ls

REAL (dp) :: t, resg, resk, errold
INTEGER   :: k, k0, k1, k2, key1

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

key1 = MAX(key ,  0)
key1 = MIN(key1,  4)
k0   = MAX(key1-2,0)
k1   = k0 + 1
k2   = MIN(key1+1,3)

CALL qxrul (f, a, b, resg, resabs, resasc, k0, k1, vr, vs, lr, ls)
errold = oflow
t = 10.0D0*epmach
DO k = k1, k2
  CALL qxrul (f, a, b, resk, resabs, resasc, k, k1, vr, vs, lr, ls)
  result = resk
  abserr = ABS(resk - resg)
  IF (resasc /= 0.0D0 .AND. abserr /= 0.0D0)  &
  abserr = resasc*MIN(1.0D0,(200.0D0*abserr/resasc)**1.5D0)
  IF (resabs > uflow/t) abserr = MAX(t*resabs,abserr)
  resg = resk
  IF (abserr > errold*0.16D0) EXIT
  IF (abserr < 1000.0D0*epmach*resabs) CYCLE
  errold = abserr
END DO

RETURN
END SUBROUTINE qxlqm


SUBROUTINE qxrul (f, xl, xu, y, ya, ym, ke, k1, fv1, fv2, l1, l2)

!    TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR ESTIMATE
!    AND CONDITIONALLY COMPUTE
!               J = INTEGRAL OF ABS(F) OVER (A,B)
!               BY USING AN  RMS RULE
!   PARAMETERS
!    ON ENTRY
!      F      - REAL
!               FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!               FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!               DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!      XL     - REAL
!               LOWER LIMIT OF INTEGRATION

!      XU     - REAL
!               UPPER LIMIT OF INTEGRATION

!      KE     - INTEGER
!             KEY FOR CHOICE OF LOCAL INTEGRATION RULE
!             AN RMS RULE IS USED WITH
!                 13      POINTS IF KE  = 2,
!                 19      POINTS IF KE  = 3,
!                 27      POINTS IF KE  = 4,
!                 42      POINTS IF KE  = 5

!      K1     INTEGER
!             VALUE OF KEY FOR WHICH THE ADDITIONAL ESTIMATES
!             YA, YM ARE TO BE COMPUTED

!      FV1    - REAL
!               VECTOR CONTAINING L1
!               SAVED  FUNCTIONAL VALUES OF POSITIVE ABSCISSAS

!      FV2    - REAL
!               VECTOR CONTAINING L2
!               SAVED  FUNCTIONAL VALUES OF NEGATIVE ABSCISSAS

!      L1     - INTEGER
!      L2       NUMBER OF ELEMENTS IN FV1,FV2  RESPECTIVELY

!    ON RETURN
!      Y      - REAL
!               APPROXIMATION TO THE INTEGRAL I
!               RESULT IS COMPUTED BY APPLYING THE REQUESTED RMS RULE

!      YA     - REAL
!               IF KEY = K1  APPROXIMATION TO THE INTEGRAL J
!               ELSE UNCHANGED

!      YM     - REAL
!               IF KEY = K1  APPROXIMATION TO THE INTEGRAL OF
!                              ABS(F-I/(XU-XL)   OVER (XL,XU)
!               ELSE UNCHANGED

!      FV1    - REAL
!               VECTOR CONTAINING L1
!               SAVED  FUNCTIONAL VALUES OF POSITIVE ABSCISSAS

!      FV2    - REAL
!               VECTOR CONTAINING L2
!               SAVED  FUNCTIONAL VALUES OF NEGATIVE ABSCISSAS

!      L1     - INTEGER
!      L2       NUMBER OF ELEMENTS IN FV1,FV2  RESPECTIVELY

!------------------------
USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)                   :: xl, xu
INTEGER, INTENT(IN)                     :: ke, k1
REAL (dp), INTENT(OUT)                  :: y, ya, ym
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: fv1, fv2
INTEGER, INTENT(IN OUT)                 :: l1, l2

REAL (dp) :: ldl, y2, aa, bb, c
INTEGER   :: istart(4) = (/ 0, 7, 17, 31 /), length(4) = (/ 7, 10, 14, 21 /), &
             j, i, is, k, ks

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

!------------------------
REAL (dp) :: xx(21) = (/ 0.D0,                    .25000000000000000000D+00, &
                       .50000000000000000000D+00, .75000000000000000000D+00, &
                       .87500000000000000000D+00, .93750000000000000000D+00, &
                       .10000000000000000000D+01, .37500000000000000000D+00, &
                       .62500000000000000000D+00, .96875000000000000000D+00, &
                       .12500000000000000000D+00, .68750000000000000000D+00, &
                       .81250000000000000000D+00, .98437500000000000000D+00, &
                       .18750000000000000000D+00, .31250000000000000000D+00, &
                       .43750000000000000000D+00, .56250000000000000000D+00, &
                       .84375000000000000000D+00, .90625000000000000000D+00, &
                       .99218750000000000000D+00 /)
REAL (dp) :: ww(52) = (/     &
1.303262173284849021810473057638590518409112513421D-1, 2.390632866847646220320329836544615917290026806242D-1,  &
2.630626354774670227333506083741355715758124943143D-1, 2.186819313830574175167853094864355208948886875898D-1,  &
2.757897646642836865859601197607471574336674206700D-2, 1.055750100538458443365034879086669791305550493830D-1,  &
1.571194260595182254168429283636656908546309467968D-2, 1.298751627936015783241173611320651866834051160074D-1,  &
2.249996826462523640447834514709508786970828213187D-1, 1.680415725925575286319046726692683040162290325505D-1,  &
1.415567675701225879892811622832845252125600939627D-1, 1.006482260551160175038684459742336605269707889822D-1,  &
2.510604860724282479058338820428989444699235030871D-2, 9.402964360009747110031098328922608224934320397592D-3,  &
5.542699233295875168406783695143646338274805359780D-2, 9.986735247403367525720377847755415293097913496236D-2,  &
4.507523056810492466415880450799432587809828791196D-2, 6.300942249647773931746170540321811473310938661469D-2,  &
1.261383225537664703012999637242003647020326905948D-1, 1.273864433581028272878709981850307363453523117880D-1,  &
8.576500414311820514214087864326799153427368592787D-2, 7.102884842310253397447305465997026228407227220665D-2,  &
5.026383572857942403759829860675892897279675661654D-2, 4.683670010609093810432609684738393586390722052124D-3,  &
1.235837891364555000245004813294817451524633100256D-1, 1.148933497158144016800199601785309838604146040215D-1,  &
1.252575774226122633391477702593585307254527198070D-2, 1.239572396231834242194189674243818619042280816640D-1,  &
2.501306413750310579525950767549691151739047969345D-2, 4.915957918146130094258849161350510503556792927578D-2,  &
2.259167374956474713302030584548274729936249753832D-2, 6.362762978782724559269342300509058175967124446839D-2,  &
9.950065827346794643193261975720606296171462239514D-2, 7.048220002718565366098742295389607994441704889441D-2,  &
6.512297339398335645872697307762912795346716454337D-2, 3.998229150313659724790527138690215186863915308702D-2,  &
3.456512257080287509832054272964315588028252136044D-2, 2.212167975884114432760321569298651047876071264944D-3,  &
8.140326425945938045967829319725797511040878579808D-2, 6.583213447600552906273539578430361199084485578379D-2,  &
2.592913726450792546064232192976262988065252032902D-2, 1.187141856692283347609436153545356484256869129472D-1,  &
5.999947605385971985589674757013565610751028128731D-2, 5.500937980198041736910257988346101839062581489820D-2,  &
5.264422421764655969760271538981443718440340270116D-3, 1.533126874056586959338368742803997744815413565014D-2,  &
3.527159369750123100455704702965541866345781113903D-2, 5.000556431653955124212795201196389006184693561679D-2,  &
5.744164831179720106340717579281831675999717767532D-2, 1.598823797283813438301248206397233634639162043386D-2,  &
2.635660410220884993472478832884065450876913559421D-2, 1.196003937945541091670106760660561117114584656319D-2 /)
!------------------------
k = ke + 1
is = istart(k)
ks = length(k)
ldl = xu - xl
bb = ldl*0.5D0
aa = xl + bb

y = 0.0D0
DO i = 1, ks
  c = bb*xx(i)
  IF (i > l1) fv1(i) = f(aa + c)
  IF (i > l2) fv2(i) = f(aa - c)
  j = is + i
  y = y + (fv1(i) + fv2(i))*ww(j)
END DO

y2 = y
y = y*bb
IF (l1 < ks) l1 = ks
IF (l2 < ks) l2 = ks
IF (ke /= k1) RETURN

ya = 0.0D0
DO i = 1, ks
  j = is + i
  ya = ya + (ABS(fv1(i)) + ABS(fv2(i)))*ww(j)
END DO
ya = ya*ABS(bb)

y2 = y2*0.5D0
ym = 0.0D0
DO i = 1, ks
  j = is + i
  ym = ym + (ABS(fv1(i) - y2) + ABS(fv2(i) - y2))*ww(j)
END DO
ym = ym*ABS(bb)
RETURN
END SUBROUTINE qxrul


SUBROUTINE qxrrd (f, z, lz, xl, xu, r, s, lr, ls)

!    TO REORDER THE COMPUTED FUNCTIONAL VALUES BEFORE
!    THE BISECTION OF AN INTERVAL

!   PARAMETERS
!    ON ENTRY
!      F      - REAL
!               FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!               FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!               DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!      XL     - REAL
!               LOWER LIMIT OF INTERVAL

!      XU     - REAL
!               UPPER LIMIT OF INTERVAL

!      Z      - REAL
!               VECTOR CONTAINING LZ SAVED FUNCTIONAL VALUES

!      LZ     - INTEGER
!               NUMBER OF ELEMENTS IN LZ

!    ON RETURN
!      R      - REAL
!      S        VECTORS CONTAINING LR, LS
!               SAVED  FUNCTIONAL VALUES FOR THE NEW INTERVALS

!      LR     - INTEGER
!      LS       NUMBER OF ELEMENTES IN R,S RESPECTIVELY

!***ROUTINES CALLED  F
USE constants_NSWC
IMPLICIT NONE

REAL (dp), DIMENSION(:), INTENT(IN)  :: z
REAL (dp), INTENT(IN)                :: xl, xu
INTEGER, INTENT(IN)                  :: lz
REAL (dp), DIMENSION(:), INTENT(OUT) :: r, s
INTEGER, INTENT(OUT)                 :: lr, ls

REAL (dp) :: dlen, centr

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

dlen = (xu - xl)*0.5D0
centr = xl + dlen
r(1) =  z(3)
r(2) =  z(9)
r(3) =  z(4)
r(4) =  z(5)
r(5) =  z(6)
r(6) =  z(10)
r(7) =  z(7)
s(1) =  z(3)
s(2) =  z(8)
s(3) =  z(2)
s(7) =  z(1)
IF (lz > 11) GO TO 10

r(8) =  f(centr + dlen*0.375D0)
r(9) =  f(centr + dlen*0.625D0)
r(10) = f(centr + dlen*0.96875D0)
lr = 10
IF (lz /= 11) s(4) = f(centr - dlen*0.75D0)
IF (lz == 11) s(4) = z(11)
s(5) =  f(centr - dlen*0.875D0)
s(6) =  f(centr - dlen*0.9375D0)
s(8) =  f(centr - dlen*0.375D0)
s(9) =  f(centr - dlen*0.625D0)
s(10) = f(centr - dlen*0.96875D0)
ls = 10
RETURN

10 r(8) = z(12)
r(9) = z(13)
r(10) = z(14)
lr = 10
s(4) = z(11)
s(5) = f(centr - dlen*0.875D0)
s(6) = f(centr - dlen*0.9375D0)
IF (lz > 14) GO TO 20
s(8)  = f(centr - dlen*0.375D0)
s(9)  = f(centr - dlen*0.625D0)
s(10) = f(centr - dlen*0.96875D0)
ls = 10
RETURN

20 r(11) = z(18)
r(12) = z(19)
r(13) = z(20)
r(14) = z(21)
lr = 14
s(8) = z(16)
s(9) = z(15)
s(10) = f(centr - dlen*0.96875D0)
s(11) = z(17)
ls = 11
RETURN
END SUBROUTINE qxrrd

END MODULE adapt_quad
