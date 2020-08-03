MODULE adapt_quad_infinite

USE constants_NSWC
IMPLICIT NONE


CONTAINS

SUBROUTINE qagi (f, bound, inf, epsabs, epsrel, result, abserr,  &
                 neval, ier, limit, last)
!-----------------------------------------------------------------------

!                   INTEGRATION OVER INFINITE INTERVALS

!-----------------------------------------------------------------------

!  PURPOSE
!     THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
!     INTEGRAL    I = INTEGRAL OF  F  OVER (BOUND,+INFINITY)
!              OR I = INTEGRAL OF  F  OVER (-INFINITY,BOUND)
!              OR I = INTEGRAL OF  F  OVER (-INFINITY,+INFINITY),
!     HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
!     ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

!  PARAMETERS
!   ON ENTRY
!      F      - REAL
!               FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!               FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!               DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!      BOUND  - REAL
!               FINITE BOUND OF INTEGRATION RANGE
!               (HAS NO MEANING IF INTERVAL IS DOUBLY-INFINITE)

!      INF    - INTEGER
!               INDICATING THE KIND OF INTEGRATION RANGE INVOLVED
!               INF = 1 CORRESPONDS TO  (BOUND,+INFINITY),
!               INF = -1            TO  (-INFINITY,BOUND),
!               INF = 2             TO (-INFINITY,+INFINITY).

!      EPSABS - REAL
!               ABSOLUTE ACCURACY REQUESTED

!      EPSREL - REAL
!               RELATIVE ACCURACY REQUESTED

!   ON RETURN
!      RESULT - REAL
!               APPROXIMATION TO THE INTEGRAL

!      ABSERR - REAL
!               ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!               WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

!      NEVAL  - INTEGER
!               NUMBER OF INTEGRAND EVALUATIONS

!      IER    - INTEGER
!               IER = 0 NORMAL AND RELIABLE TERMINATION OF THE ROUTINE.
!                       IT IS ASSUMED THAT THE REQUESTED
!                       ACCURACY HAS BEEN ACHIEVED.
!             - IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE.  THE ESTIMATES
!                       FOR RESULT AND ERROR ARE LESS RELIABLE. IT IS ASSUMED
!                       THAT THE REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.

!               IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED HAS BEEN
!                       ACHIEVED.  ONE CAN ALLOW MORE SUBDIVISIONS BY
!                       INCREASING THE VALUE OF LIMIT (AND TAKING THE
!                       ACCORDING DIMENSION ADJUSTMENTS INTO ACCOUNT).
!                       HOWEVER, IF THIS YIELDS NO IMPROVEMENT IT IS ADVISED
!                       TO ANALYZE THE INTEGRAND IN ORDER TO DETERMINE THE
!                       INTEGRATION DIFFICULTIES.
!                       IF THE POSITION OF A LOCAL DIFFICULTY CAN BE DETERMINED
!                       (E.G. SINGULARITY, DISCONTINUITY WITHIN THE INTERVAL)
!                       ONE WILL PROBABLY GAIN FROM SPLITTING UP THE INTERVAL
!                       AT THIS POINT AND CALLING THE INTEGRATOR ON THE
!                       SUBRANGES.
!                       IF POSSIBLE, AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
!                       SHOULD BE USED, WHICH IS DESIGNED FOR HANDLING THE
!                       TYPE OF DIFFICULTY INVOLVED.
!                   = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETECTED, WHICH
!                       PREVENTS THE REQUESTED TOLERANCE FROM BEING ACHIEVED.
!                       THE ERROR MAY BE UNDER-ESTIMATED.
!                   = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS AT SOME POINTS
!                       OF THE INTEGRATION INTERVAL.
!                   = 4 THE ALGORITHM DOES NOT CONVERGE.
!                       ROUNDOFF ERROR IS DETECTED IN THE EXTRAPOLATION TABLE.
!                       IT IS ASSUMED THAT THE REQUESTED TOLERANCE CANNOT BE
!                       ACHIEVED, AND THAT THE RETURNED RESULT IS THE BEST
!                       WHICH CAN BE OBTAINED.
!                   = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY
!                       CONVERGENT. IT MUST BE NOTED THAT DIVERGENCE CAN OCCUR
!                       WITH ANY OTHER VALUE OF IER.
!                   = 6 THE INPUT IS INVALID BECAUSE EPSABS OR EPSREL IS
!                       NEGATIVE, LIMIT .LT. 1, OR LENW .LT. 4 * LIMIT.
!                       RESULT, ABSERR, NEVAL, LAST ARE SET TO ZERO.

!   DIMENSIONING PARAMETERS
!      LIMIT - INTEGER
!              DIMENSIONING PARAMETER FOR IWORK
!              LIMIT DETERMINES THE MAXIMUM NUMBER OF SUBINTERVALS IN THE
!              PARTITION OF THE GIVEN INTEGRATION INTERVAL (A,B), LIMIT.GE.1.
!              IF LIMIT.LT.1, THE ROUTINE WILL END WITH IER = 6.

!      LENW  - INTEGER
!              DIMENSIONING PARAMETER FOR WORK
!              LENW MUST BE AT LEAST LIMIT*4.
!              IF LENW.LT.LIMIT*4, THE ROUTINE WILL END WITH IER = 6.

!      LAST  - INTEGER
!              ON RETURN, LAST EQUALS THE NUMBER OF SUBINTERVALS PRODUCED IN
!              THE SUBDIVISION PROCESS, WHICH DETERMINES THE NUMBER OF
!              SIGNIFICANT ELEMENTS ACTUALLY IN THE WORK ARRAYS.

!  SUBROUTINES OR FUNCTIONS NEEDED
!        - QAGIE
!        - QK15I
!        - QPSRT
!        - QELG
!        - F (USER PROVIDED FUNCTION)
!        - SPMPAR

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: bound, epsabs, epsrel
INTEGER, INTENT(IN)     :: inf, limit
REAL (dp), INTENT(OUT)  :: result, abserr
INTEGER, INTENT(OUT)    :: neval, ier, last

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fx
  END FUNCTION f
END INTERFACE

!         CHECK VALIDITY OF LIMIT.

ier = 6
neval = 0
last = 0
result = 0.0D0
abserr = 0.0D0
IF (limit < 1) RETURN

!         PREPARE CALL FOR QAGIE.

CALL qagie (f, bound, inf, epsabs, epsrel, limit, result, abserr,  &
            neval, ier, last)
RETURN
END SUBROUTINE qagi



SUBROUTINE qagie (f, bound, inf, epsabs, epsrel, limit, result, abserr,  &
                  neval, ier, last)
!-----------------------------------------------------------------------

!                   INTEGRATION OVER INFINITE INTERVALS

!-----------------------------------------------------------------------

!  PURPOSE
!     THE ROUTINE CALCULATES AN APPROXIMATION RESULT TO A GIVEN INTEGRAL
!                 I = INTEGRAL OF  F  OVER (BOUND,+INFINITY)
!              OR I = INTEGRAL OF  F  OVER (-INFINITY,BOUND)
!              OR I = INTEGRAL OF  F  OVER (-INFINITY,+INFINITY),
!     HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
!     ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

!  PARAMETERS
!   ON ENTRY
!      F      - REAL
!               FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
!               FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
!               DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

!      BOUND  - REAL
!               FINITE BOUND OF INTEGRATION RANGE
!               (HAS NO MEANING IF INTERVAL IS DOUBLY-INFINITE)

!      INF    - INTEGER
!               INDICATING THE KIND OF INTEGRATION RANGE INVOLVED
!               INF = 1 CORRESPONDS TO  (BOUND,+INFINITY),
!               INF = -1            TO  (-INFINITY,BOUND),
!               INF = 2             TO (-INFINITY,+INFINITY).

!      EPSABS - REAL
!               ABSOLUTE ACCURACY REQUESTED

!      EPSREL - REAL
!               RELATIVE ACCURACY REQUESTED

!      LIMIT  - INTEGER
!               GIVES AN UPPER BOUND ON THE NUMBER OF SUBINTERVALS IN THE
!               PARTITION OF (A,B),
!               LIMIT.GE.1

!   ON RETURN
!      RESULT - REAL
!               APPROXIMATION TO THE INTEGRAL

!      ABSERR - REAL
!               ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!               WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

!      NEVAL  - INTEGER
!               NUMBER OF INTEGRAND EVALUATIONS

!      IER    - INTEGER
!               IER = 0 NORMAL AND RELIABLE TERMINATION OF THE ROUTINE.
!                       IT IS ASSUMED THAT THE REQUESTED ACCURACY HAS BEEN
!                       ACHIEVED.
!             - IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE.  THE ESTIMATES
!                       FOR RESULT AND ERROR ARE LESS RELIABLE.  IT IS ASSUMED
!                       THAT THE REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.

!               IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED HAS BEEN
!                       ACHIEVED.  ONE CAN ALLOW MORE SUBDIVISIONS BY
!                       INCREASING THE VALUE OF LIMIT (AND TAKING THE
!                       ACCORDING DIMENSION ADJUSTMENTS INTO ACCOUNT).
!                       HOWEVER, IF THIS YIELDS NO IMPROVEMENT IT IS ADVISED
!                       TO ANALYZE THE INTEGRAND IN ORDER TO DETERMINE THE
!                       INTEGRATION DIFFICULTIES.
!                       IF THE POSITION OF A LOCAL DIFFICULTY CAN BE
!                       DETERMINED (E.G. SINGULARITY, DISCONTINUITY WITHIN THE
!                       INTERVAL) ONE WILL PROBABLY GAIN FROM SPLITTING UP THE
!                       INTERVAL AT THIS POINT AND CALLING THE INTEGRATOR ON
!                       THE SUBRANGES.  IF POSSIBLE, AN APPROPRIATE SPECIAL-
!                       PURPOSE INTEGRATOR SHOULD BE USED, WHICH IS DESIGNED
!                       FOR HANDLING THE TYPE OF DIFFICULTY INVOLVED.
!                   = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETECTED, WHICH
!                       PREVENTS THE REQUESTED TOLERANCE FROM BEING ACHIEVED.
!                       THE ERROR MAY BE UNDER-ESTIMATED.
!                   = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS AT SOME
!                       POINTS OF THE INTEGRATION INTERVAL.
!                   = 4 THE ALGORITHM DOES NOT CONVERGE.
!                       ROUNDOFF ERROR IS DETECTED IN THE EXTRAPOLATION TABLE.
!                       IT IS ASSUMED THAT THE REQUESTED TOLERANCE CANNOT BE
!                       ACHIEVED, AND THAT THE RETURNED RESULT IS THE BEST
!                       WHICH CAN BE OBTAINED.
!                   = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY
!                       CONVERGENT.  IT MUST BE NOTED THAT DIVERGENCE CAN
!                       OCCUR WITH ANY OTHER VALUE OF IER.
!                   = 6 THE INPUT IS INVALID BECAUSE EPSABS OR EPSREL IS
!                       NEGATIVE.
!                       RESULT, ABSERR, NEVAL, LAST, RLIST(1), ELIST(1) AND
!                       IORD(1) ARE SET TO ZERO.
!                       ALIST(1) AND BLIST(1) ARE SET TO 0 AND 1 RESPECTIVELY.

!      ALIST  - REAL
!               VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE LEFT
!               END POINTS OF THE SUBINTERVALS IN THE PARTITION
!               OF THE TRANSFORMED INTEGRATION RANGE (0,1).

!      BLIST  - REAL
!               VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE RIGHT
!               END POINTS OF THE SUBINTERVALS IN THE PARTITION
!               OF THE TRANSFORMED INTEGRATION RANGE (0,1).

!      RLIST  - REAL
!               VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
!               APPROXIMATIONS ON THE SUBINTERVALS

!      ELIST  - REAL
!               VECTOR OF DIMENSION AT LEAST LIMIT,  THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE MODULI
!               OF THE ABSOLUTE ERROR ESTIMATES ON THE
!               SUBINTERVALS

!      IORD   - INTEGER
!               VECTOR OF DIMENSION LIMIT, THE FIRST K ELEMENTS OF WHICH ARE
!               POINTERS TO THE ERROR ESTIMATES OVER THE SUBINTERVALS,
!               SUCH THAT ELIST(IORD(1)), ..., ELIST(IORD(K)) FORM A
!               DECREASING SEQUENCE, WITH K = LAST IF LAST.LE.(LIMIT/2+2),
!               AND K = LIMIT+1-LAST OTHERWISE

!      LAST   - INTEGER
!               NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN THE SUBDIVISION
!               PROCESS

!  SUBROUTINES OR FUNCTIONS NEEDED
!        - QK15I
!        - QPSRT
!        - QELG
!        - F (USER-PROVIDED FUNCTION)
!        - SPMPAR

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: bound, epsabs, epsrel
REAL (dp), INTENT(OUT)  :: result, abserr
INTEGER, INTENT(IN)     :: inf, limit
INTEGER, INTENT(OUT)    :: neval, ier, last

REAL (dp) :: abseps, area, area1, area12, area2, a1, a2, boun, b1, b2,     &
             correc, defabs, defab1, defab2, dres, epmach, erlarg, erlast, &
             errbnd, errmax, error1, error2, erro12, errsum, ertest,       &
             oflow, rerr, resabs, reseps, res3la(3), rlist2(52), small, t, &
             uflow
INTEGER :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin,  &
           maxerr, nres, nrmax, numrl2
LOGICAL :: extrap, noext

REAL (dp), DIMENSION(limit)  :: alist, blist, elist, rlist
INTEGER, DIMENSION(limit)    :: iord

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

!      THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
!      LIMEXP IN SUBROUTINE QELG.


!      LIST OF MAJOR VARIABLES
!      -----------------------

!     ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS CONSIDERED UP
!                 TO NOW
!     BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS CONSIDERED UP
!                 TO NOW
!     RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER (ALIST(I),BLIST(I))
!     RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+2), CONTAINING THE PART
!                 OF THE EPSILON TABLE WHICH IS STILL NEEDED FOR FURTHER
!                 COMPUTATIONS
!     ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!     MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR ESTIMATE
!     ERRMAX    - ELIST(MAXERR)
!     ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
!                 (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!     AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!     ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!     ERRBND    - REQUESTED ACCURACY MAX(EPSABS, EPSREL*ABS(RESULT))
!     *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
!     *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
!     LAST      - INDEX FOR SUBDIVISION
!     NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
!     NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2.  IF AN APPROPRIATE
!                 APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS BEEN OBTAINED,
!                 IT IS PUT IN RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
!                 BY ONE.
!     SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP TO NOW,
!                 MULTIPLIED BY 1.5
!     ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER THAN THE
!                 SMALLEST INTERVAL CONSIDERED UP TO NOW
!     EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS ATTEMPTING TO
!                 PERFORM EXTRAPOLATION. I.E. BEFORE SUBDIVIDING THE SMALLEST
!                 INTERVAL WE TRY TO DECREASE THE VALUE OF ERLARG.
!     NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS NO LONGER
!                 ALLOWED (TRUE-VALUE)

!      MACHINE DEPENDENT CONSTANTS
!      ---------------------------

!     EPMACH IS THE LARGEST RELATIVE SPACING.
!     UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!     OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

epmach = dpmpar(1)
uflow = dpmpar(2)
oflow = dpmpar(3)

!           CHECK EPSABS AND EPSREL
!           -----------------------

neval = 0
last = 0
result = 0.0D0
abserr = 0.0D0
alist(1) = 0.0D0
blist(1) = 1.0D0
rlist(1) = 0.0D0
elist(1) = 0.0D0
iord(1) = 0
ier = 6
IF (epsabs < 0.0D0 .OR. epsrel < 0.0D0) GO TO 999
ier = 0
rerr = MAX(epsrel, 50.0D0*epmach, 0.5D-14)

!     FIRST APPROXIMATION TO THE INTEGRAL
!     -----------------------------------

!     DETERMINE THE INTERVAL TO BE MAPPED ONTO (0,1).
!     IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2, WHERE
!     I1 = INTEGRAL OF F OVER (-INFINITY,0),
!     I2 = INTEGRAL OF F OVER (0,+INFINITY).

boun = bound
IF (inf == 2) boun = 0.0D0
CALL qk15i (f, boun, inf, 0.0D0, 1.0D0, result, abserr, defabs, resabs, &
            epmach, uflow)

!           TEST ON ACCURACY

last = 1
rlist(1) = result
elist(1) = abserr
iord(1) = 1
dres = ABS(result)
errbnd = MAX(epsabs,rerr*dres)
IF (abserr <= 100.0D0*epmach*defabs .AND. abserr > errbnd)  &
ier = 2
IF (limit == 1) ier = 1
IF (ier /= 0 .OR. (abserr <= errbnd .AND. abserr /= resabs)  &
.OR. abserr == 0.0D0) GO TO 130

!           INITIALIZATION
!           --------------

rlist2(1) = result
errmax = abserr
maxerr = 1
area = result
errsum = abserr
abserr = oflow
correc = 0.0D0
nrmax = 1
nres = 0
ktmin = 0
numrl2 = 2
extrap = .false.
noext = .false.
ierro = 0
iroff1 = 0
iroff2 = 0
iroff3 = 0
ksgn = -1
IF (dres >= (1.0D0 - 50.0D0*epmach)*defabs) ksgn = 1
t = 1.0D0 + 100.0D0*epmach

!           MAIN DO-LOOP
!           ------------

DO last = 2,limit
  
!           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
  
  a1 = alist(maxerr)
  b1 = 0.5D0*(alist(maxerr) + blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  erlast = errmax
  CALL qk15i (f, boun, inf, a1, b1, area1, error1,  &
              resabs, defab1, epmach, uflow)
  CALL qk15i (f, boun, inf, a2, b2, area2, error2,  &
              resabs, defab2, epmach, uflow)
  
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR
!           AND TEST FOR ACCURACY.
  
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
  
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
  
  IF (iroff1 + iroff2 >= 10 .OR. iroff3 >= 20) ier = 2
  IF (iroff2 >= 5) ierro = 3
  
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
!           SUBINTERVALS EQUALS LIMIT.
  
  IF (last == limit) ier = 1
  
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT SOME POINTS OF THE INTEGRATION RANGE.
  
  
  IF (MAX(ABS(a1),ABS(b2)) <= t*(ABS(a2) + 0.1D+04*uflow)) ier = 4
  
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
  
  IF (error2 > error1) GO TO 20
  alist(last) = a2
  blist(maxerr) = b1
  blist(last) = b2
  elist(maxerr) = error1
  elist(last) = error2
  GO TO 30
  20 alist(maxerr) = a2
  alist(last) = a1
  blist(last) = b1
  rlist(maxerr) = area2
  rlist(last) = area1
  elist(maxerr) = error2
  elist(last) = error1
  
!           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL WITH
!           NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
  
  30 CALL qpsrt (limit,last,maxerr,errmax,elist,iord,nrmax)
  IF (errsum <= errbnd) GO TO 115
  IF (ier /= 0) GO TO 100
  IF (last == 2) GO TO 80
  IF (noext) CYCLE
  erlarg = erlarg - erlast
  IF (ABS(b1 - a1) > small) erlarg = erlarg + erro12
  IF (extrap) GO TO 40
  
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
  
  IF (ABS(blist(maxerr) - alist(maxerr)) > small) CYCLE
  extrap = .true.
  nrmax = 2
  40 IF (ierro == 3 .OR. erlarg <= ertest) GO TO 60
  
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS
!           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
  
  id = nrmax
  jupbnd = last
  IF (last > (2 + limit/2)) jupbnd = limit + 3 - last
  DO k = id,jupbnd
    maxerr = iord(nrmax)
    errmax = elist(maxerr)
    IF (ABS(blist(maxerr) - alist(maxerr)) > small) CYCLE
    nrmax = nrmax + 1
  END DO
  
!           PERFORM EXTRAPOLATION.
  
  60 numrl2 = numrl2 + 1
  rlist2(numrl2) = area
  CALL qelg (numrl2, rlist2, reseps, abseps, res3la, nres,epmach, oflow)
  ktmin = ktmin + 1
  IF (ktmin > 5 .AND. abserr < 0.1D-02*errsum) ier = 5
  IF (abseps >= abserr) GO TO 70
  ktmin = 0
  abserr = abseps
  result = reseps
  correc = erlarg
  ertest = MAX(epsabs,rerr*ABS(reseps))
  IF (abserr <= ertest) GO TO 100
  
!            PREPARE BISECTION OF THE SMALLEST INTERVAL.
  
  70 IF (numrl2 == 1) noext = .true.
  IF (ier == 5) GO TO 100
  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .false.
  small = small*0.5D0
  erlarg = errsum
  CYCLE
  80 small = 0.375D0
  erlarg = errsum
  ertest = errbnd
  rlist2(2) = area
END DO

!           SET FINAL RESULT AND ERROR ESTIMATE.
!           ------------------------------------

100 IF (abserr == oflow) GO TO 115
IF (ier + ierro == 0) GO TO 110
IF (ierro == 3) abserr = abserr + correc
IF (ier == 0) ier = 3
IF (result /= 0.0D0 .AND. area /= 0.0D0) GO TO 105
IF (abserr > errsum) GO TO 115
IF (area == 0.0D0) GO TO 130
GO TO 110
105 IF (abserr/ABS(result) > errsum/ABS(area)) GO TO 115

!           TEST ON DIVERGENCE

110 IF (ksgn == -1 .AND. MAX(ABS(result),ABS(area)) <=  &
defabs*0.1D-01) GO TO 130
IF (0.1D-01 > (result/area) .OR. (result/area) > 0.1D+03  &
.OR. errsum > ABS(area)) ier = 6
GO TO 130

!           COMPUTE GLOBAL INTEGRAL SUM.

115 result = 0.0D0
DO k = 1,last
  result = result + rlist(k)
END DO
abserr = errsum
130 neval = 30*last - 15
IF (inf == 2) neval = 2*neval
IF (ier > 2) ier = ier - 1
999 RETURN
END SUBROUTINE qagie



SUBROUTINE qk15i (f, boun, inf, a, b, result, abserr, resabs,  &
                  resasc, epmach, uflow)
!-----------------------------------------------------------------------

! 1.  PURPOSE
!        THE ORIGINAL (INFINITE) INTEGRATION RANGE IS MAPPED
!        ONTO THE INTERVAL (0,1) AND (A,B) IS A PART OF (0,1).
!        IT IS THE PURPOSE TO COMPUTE
!        I = INTEGRAL OF TRANSFORMED INTEGRAND OVER (A,B),
!        J = INTEGRAL OF ABS(TRANSFORMED INTEGRAND) OVER (A,B).

! 2.  PARAMETERS
!      ON ENTRY
!        F      - REAL
!                 FUNCTION SUBPROGRAM DEFINING THE INTEGRAND FUNCTION F(X).
!                 THE ACTUAL NAME FOR F NEEDS TO BE DECLARED E X T E R N A L
!                 IN THE CALLING PROGRAM.

!        BOUN   - REAL
!                 FINITE BOUND OF ORIGINAL INTEGRATION
!                 RANGE (SET TO ZERO IF INF = +2)

!        INF    - INTEGER
!                 IF INF = -1, THE ORIGINAL INTERVAL IS
!                             (-INFINITY,BOUND),
!                 IF INF = +1, THE ORIGINAL INTERVAL IS
!                             (BOUND,+INFINITY),
!                 IF INF = +2, THE ORIGINAL INTERVAL IS
!                             (-INFINITY,+INFINITY) AND
!                 THE INTEGRAL IS COMPUTED AS THE SUM OF TWO INTEGRALS,
!                 ONE OVER (-INFINITY,0) AND ONE OVER (0,+INFINITY).

!        A      - REAL
!                 LOWER LIMIT FOR INTEGRATION OVER SUBRANGE OF (0,1)

!        B      - REAL
!                 UPPER LIMIT FOR INTEGRATION OVER SUBRANGE OF (0,1)

!        EPMACH - REAL
!                 THE RELATIVE PRECISION OF THE FLOATING ARITHMETIC BEING USED.

!        UFLOW  - REAL
!                 THE SMALLEST POSITIVE MAGNITUDE.

!      ON RETURN
!        RESULT - REAL
!                 APPROXIMATION TO THE INTEGRAL I
!                 RESULT IS COMPUTED BY APPLYING THE 15-POINT KRONROD RULE
!                 (RESK) OBTAINED BY OPTIMAL ADDITION OF ABSCISSAE TO THE
!                 7-POINT GAUSS RULE(RESG).

!        ABSERR - REAL
!                 ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
!                 WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

!        RESABS - REAL
!                 APPROXIMATION TO THE INTEGRAL J

!        RESASC - REAL
!                 APPROXIMATION TO THE INTEGRAL OF
!                 ABS((TRANSFORMED INTEGRAND)-I/(B-A)) OVER (A,B)

! 3.  SUBROUTINES OR FUNCTIONS NEEDED
!           - F (USER-PROVIDED FUNCTION)

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: boun, a, b, epmach, uflow
INTEGER, INTENT(IN)     :: inf
REAL (dp), INTENT(OUT)  :: result, abserr, resabs, resasc

REAL (dp) :: fv1(7), fv2(7)

INTERFACE
  FUNCTION f(x) RESULT(fx)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION f
END INTERFACE

!     THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
!     (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
!     THEIR CORRESPONDING WEIGHTS ARE GIVEN.

!     XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
!              XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT GAUSS RULE
!              XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!              ADDED TO THE 7-POINT GAUSS RULE

!     WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE

!     WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING TO THE
!              ABSCISSAE XGK(2), XGK(4), ... WG(1), WG(3), ...
!              ARE SET TO ZERO.

REAL (dp) :: absc, absc1, absc2,  centr, dinf, fc, fsum, fval1, fval2,  &
             hlgth, resg, resk, reskh, tabsc1, tabsc2, tol
INTEGER   :: j
REAL (dp), DIMENSION(8) :: xgk = (/  &
                           0.9914553711208126D+00, 0.9491079123427585D+00,  &
                           0.8648644233597691D+00, 0.7415311855993944D+00,  &
                           0.5860872354676911D+00, 0.4058451513773972D+00,  &
                           0.2077849550078985D+00, 0.0000000000000000D+00 /), &
                           wgk = (/  &
                           0.2293532201052922D-01, 0.6309209262997855D-01,  &
                           0.1047900103222502D+00, 0.1406532597155259D+00,  &
                           0.1690047266392679D+00, 0.1903505780647854D+00,  &
                           0.2044329400752989D+00, 0.2094821410847278D+00 /), &
                           wg = (/   &
                           0.0000000000000000D+00, 0.1294849661688697D+00,  &
                           0.0000000000000000D+00, 0.2797053914892767D+00,  &
                           0.0000000000000000D+00, 0.3818300505051189D+00,  &
                           0.0000000000000000D+00, 0.4179591836734694D+00 /)

!     LIST OF MAJOR VARIABLES
!     -----------------------

!     CENTR  - MID POINT OF THE INTERVAL
!     HLGTH  - HALF-LENGTH OF THE INTERVAL
!     ABSC*  - ABSCISSA
!     TABSC* - TRANSFORMED ABSCISSA
!     FVAL*  - FUNCTION VALUE
!     RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
!     RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
!     RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
!              INTEGRAND OVER (A,B), I.E. TO I/(B-A)

dinf = MIN(1,inf)

centr = 0.5D0*(a + b)
hlgth = 0.5D0*(b - a)
tabsc1 = boun + dinf*(1.0D0 - centr)/centr
fval1 = f(tabsc1)
IF (inf == 2) fval1 = fval1 + f(-tabsc1)
fc = (fval1/centr)/centr

!           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL,
!           AND ESTIMATE THE ERROR.

resg = wg(8)*fc
resk = wgk(8)*fc
resabs = ABS(resk)
DO j = 1,7
  absc = hlgth*xgk(j)
  absc1 = centr - absc
  absc2 = centr + absc
  tabsc1 = boun + dinf*(1.0D0 - absc1)/absc1
  tabsc2 = boun + dinf*(1.0D0 - absc2)/absc2
  fval1 = f(tabsc1)
  fval2 = f(tabsc2)
  IF (inf == 2) fval1 = fval1 + f(-tabsc1)
  IF (inf == 2) fval2 = fval2 + f(-tabsc2)
  fval1 = (fval1/absc1)/absc1
  fval2 = (fval2/absc2)/absc2
  fv1(j) = fval1
  fv2(j) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j)*fsum
  resk = resk + wgk(j)*fsum
  resabs = resabs + wgk(j)*(ABS(fval1) + ABS(fval2))
END DO
reskh = resk / 2
resasc = wgk(8)*ABS(fc - reskh)
DO j = 1,7
  resasc = resasc + wgk(j)*(ABS(fv1(j)-reskh) + ABS(fv2(j)-reskh))
END DO
result = resk*hlgth
resasc = resasc*hlgth
resabs = resabs*hlgth
abserr = ABS((resk - resg)*hlgth)
IF (resasc /= 0.0D0 .AND. abserr /= 0.0D0) abserr = resasc*  &
                           MIN(1.0D0, (0.2D+03*abserr/resasc)**1.5D0)
tol = 50.0D0*epmach
IF (resabs > uflow/tol) abserr = MAX(abserr, tol*resabs)

RETURN
END SUBROUTINE qk15i



SUBROUTINE qpsrt(limit, last, maxerr, ermax, elist, iord, nrmax)
!     ..................................................................

! 1.  QPSRT
!     ORDERING ROUTINE
!        STANDARD FORTRAN SUBROUTINE
!        REAL VERSION

! 2.  PURPOSE
!        THIS ROUTINE MAINTAINS THE DESCENDING ORDERING IN THE LIST OF THE
!        LOCAL ERROR ESTIMATES RESULTING FROM THE INTERVAL SUBDIVISION
!        PROCESS.  AT EACH CALL TWO ERROR ESTIMATES ARE INSERTED USING THE
!        SEQUENTIAL SEARCH METHOD, TOP-DOWN FOR THE LARGEST ERROR ESTIMATE
!        AND BOTTOM-UP FOR THE SMALLEST ERROR ESTIMATE.

! 3.  CALLING SEQUENCE
!        CALL QPSRT(LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)

!     PARAMETERS (MEANING AT OUTPUT)
!        LIMIT  - INTEGER
!                 MAXIMUM NUMBER OF ERROR ESTIMATES THE LIST CAN CONTAIN

!        LAST   - INTEGER
!                 NUMBER OF ERROR ESTIMATES CURRENTLY IN THE LIST

!        MAXERR - INTEGER
!                 MAXERR POINTS TO THE NRMAX-TH LARGEST ERROR ESTIMATE
!                 CURRENTLY IN THE LIST

!        ERMAX  - REAL
!                 NRMAX-TH LARGEST ERROR ESTIMATE
!                 ERMAX = ELIST(MAXERR)

!        ELIST  - REAL
!                 VECTOR OF DIMENSION LAST CONTAINING THE ERROR ESTIMATES

!        IORD   - INTEGER
!                 VECTOR OF DIMENSION LAST, THE FIRST K ELEMENTS OF
!                 WHICH CONTAIN POINTERS TO THE ERROR ESTIMATES,
!                 SUCH THAT ELIST(IORD(1)), ... , ELIST(IORD(K))
!                 FORM A DECREASING SEQUENCE, WITH K = LAST IF
!                 LAST <= (LIMIT/2+2), AND K = LIMIT+1-LAST OTHERWISE

!        NRMAX  - INTEGER
!                 MAXERR = IORD(NRMAX)

! 4.  NO SUBROUTINES OR FUNCTIONS NEEDED

!     ..................................................................


INTEGER, INTENT(IN)                  :: limit, last
REAL (dp), DIMENSION(:), INTENT(IN)  :: elist
INTEGER, INTENT(IN OUT)              :: nrmax
INTEGER, DIMENSION(:), INTENT(OUT)   :: iord
INTEGER, INTENT(OUT)                 :: maxerr
REAL (dp), INTENT(OUT)               :: ermax

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

!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO
!           BE MAINTAINED IN DESCENDING ORDER. THIS NUMBER
!           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.

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

! 1.  PURPOSE
!        THE ROUTINE DETERMINES THE LIMIT OF A GIVEN SEQUENCE OF
!        APPROXIMATIONS, BY MEANS OF THE EPSILON ALGORITHM OF P. WYNN.
!        AN ESTIMATE OF THE ABSOLUTE ERROR IS ALSO GIVEN.
!        THE CONDENSED EPSILON TABLE IS COMPUTED. ONLY THOSE ELEMENTS NEEDED
!        FOR THE COMPUTATION OF THE NEXT DIAGONAL ARE PRESERVED.

! 2.  PARAMETERS
!        N      - INTEGER
!                 EPSTAB(N) CONTAINS THE NEW ELEMENT IN THE
!                 FIRST COLUMN OF THE EPSILON TABLE.

!        EPSTAB - REAL
!                 VECTOR OF DIMENSION 52 CONTAINING THE ELEMENTS OF THE TWO
!                 LOWER DIAGONALS OF THE TRIANGULAR EPSILON TABLE.
!                 THE ELEMENTS ARE NUMBERED STARTING AT THE RIGHT-HAND
!                 CORNER OF THE TRIANGLE.

!        RESULT - REAL
!                 RESULTING APPROXIMATION TO THE INTEGRAL

!        ABSERR - REAL
!                 ESTIMATE OF THE ABSOLUTE ERROR COMPUTED FROM
!                 RESULT AND THE 3 PREVIOUS RESULTS

!        RES3LA - REAL
!                 VECTOR OF DIMENSION 3 CONTAINING THE LAST 3 RESULTS

!        NRES   - INTEGER
!                 NUMBER OF CALLS TO THE ROUTINE
!                 (SHOULD BE ZERO AT FIRST CALL)

!        EPMACH - REAL
!                 THE RELATIVE PRECISION OF THE FLOATING ARITHMETIC BEING USED.

!        OFLOW  - REAL
!                 THE LARGEST POSITIVE MAGNITUDE.

! 3.  NO SUBROUTINES OR FUNCTIONS USED

!-----------------------------------------------------------------------

INTEGER, INTENT(IN OUT)                  :: n, nres
REAL (dp), INTENT(IN)                    :: epmach, oflow
REAL (dp), INTENT(OUT)                   :: abserr, result
REAL (dp), DIMENSION(:), INTENT(IN OUT)  :: epstab, res3la
!---------------------

!     LIST OF MAJOR VARIABLES
!     -----------------------

!     E0     - THE 4 ELEMENTS ON WHICH THE
!     E1       COMPUTATION OF A NEW ELEMENT IN
!     E2       THE EPSILON TABLE IS BASED
!     E3                 E0
!                  E3    E1    NEW
!                        E2
!     NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW DIAGONAL
!     ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!     RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE OF ERROR

!     LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON TABLE CAN CONTAIN.
!     IF THIS NUMBER IS REACHED, THE UPPER DIAGONAL OF THE EPSILON TABLE IS
!     DELETED.

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

!           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE ACCURACY,
!           CONVERGENCE IS ASSUMED.
!           RESULT = E2
!           ABSERR = ABS(E1-E0) + ABS(E2-E1)

  result = res
  abserr = err2 + err3
! ***JUMP OUT OF DO-LOOP
  GO TO 100
  10 e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1 - e3
  err1 = ABS(delta1)
  tol1 = MAX(e1abs,ABS(e3))*epmach

!           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N

  IF (err1 <= tol1 .OR. err2 <= tol2 .OR. err3 <= tol3) GO TO 20
  ss = 1.0D0/delta1 + 1.0D0/delta2 - 1.0D0/delta3
  epsinf = ABS(ss*e1)

!           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND EVENTUALLY
!           OMIT A PART OF THE TABLE ADJUSTING THE VALUE OF N.

  IF (epsinf > 0.1D-03) GO TO 30
  20 n = i + i - 1
! ***JUMP OUT OF DO-LOOP
  GO TO 50

!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST THE VALUE OF RESULT.

  30 res = e1 + 1.0D0/ss
  epstab(k1) = res
  k1 = k1 - 2
  error = err2 + ABS(res - e2) + err3
  IF (error > abserr) CYCLE
  abserr = error
  result = res
END DO

!           SHIFT THE TABLE.

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

!           COMPUTE ERROR ESTIMATE

90 abserr = ABS(result - res3la(3)) + ABS(result - res3la(2)) +  &
            ABS(result - res3la(1))
res3la(1) = res3la(2)
res3la(2) = res3la(3)
res3la(3) = result
100 abserr = MAX(abserr,5.0D0*epmach*ABS(result))
RETURN
END SUBROUTINE qelg

END MODULE adapt_quad_infinite
