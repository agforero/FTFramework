MODULE qshep2d

! ALGORITHM 660, COLLECTED ALGORITHMS FROM ACM.

! THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
! VOL. 14, NO. 2, P.149.

! Code converted using TO_F90 by Alan Miller
! Date: 2001-01-18  Time: 10:10:42

IMPLICIT NONE


CONTAINS


SUBROUTINE qshep2(n,x,y,f,nq,nw,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a,ier)
INTEGER, INTENT(IN)   :: n, nq, nw, nr
REAL, INTENT(IN)      :: x(:), y(:), f(:)
INTEGER, INTENT(OUT)  :: lcell(:,:), lnext(:), ier
REAL, INTENT(OUT)     :: xmin, ymin, dx, dy, rmax, rsq(n), a(:,:)

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    01/08/90

!   THIS SUBROUTINE COMPUTES A SET OF PARAMETERS A AND RSQ DEFINING A SMOOTH
! (ONCE CONTINUOUSLY DIFFERENTIABLE) BI-VARIATE FUNCTION Q(X,Y) WHICH
! INTERPOLATES DATA VALUES F AT SCATTERED NODES (X,Y).  THE INTERPOLANT Q MAY
! BE EVALUATED AT AN ARBITRARY POINT BY FUNCTION QS2VAL, AND ITS FIRST
! DERIVATIVES ARE COMPUTED BY SUBROUTINE QS2GRD.

!   THE INTERPOLATION SCHEME IS A MODIFIED QUADRATIC SHEPARD METHOD --

! Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))

! FOR BIVARIATE FUNCTIONS W(K) AND Q(K).  THE NODAL FUNCTIONS ARE GIVEN BY

! Q(K)(X,Y) =  A(1,K)*(X-X(K))**2 + A(2,K)*(X-X(K))*(Y-Y(K))
!            + A(3,K)*(Y-Y(K))**2 + A(4,K)*(X-X(K))
!            + A(5,K)*(Y-Y(K))         + F(K) .

! THUS, Q(K) IS A QUADRATIC FUNCTION WHICH INTERPOLATES THE DATA VALUE AT
! NODE K.  ITS COEFFICIENTS A(,K) ARE OBTAINED BY A WEIGHTED LEAST SQUARES FIT
! TO THE CLOSEST NQ DATA POINTS WITH WEIGHTS SIMILAR TO W(K).
! NOTE THAT THE RADIUS OF INFLUENCE FOR THE LEAST SQUARES FIT IS FIXED FOR
! EACH K, BUT VARIES WITH K.

!   THE WEIGHTS ARE TAKEN TO BE

! W(K)(X,Y) = ( (R(K)-D(K))+ / R(K)*D(K) )**2

! WHERE (R(K)-D(K))+ = 0 IF R(K) <= D(K) AND D(K)(X,Y) IS THE EUCLIDEAN
! DISTANCE BETWEEN (X,Y) AND (X(K),Y(K)).  THE RADIUS OF INFLUENCE R(K) VARIES
! WITH K AND IS CHOSEN SO THAT NW NODES ARE WITHIN THE RADIUS.
! NOTE THAT W(K) IS NOT DEFINED AT NODE (X(K),Y(K)), BUT Q(X,Y) HAS LIMIT F(K)
! AS (X,Y) APPROACHES (X(K),Y(K)).

! ON INPUT --

!        N = NUMBER OF NODES AND ASSOCIATED DATA VALUES.
!            N >= 6.

!        X,Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!              COORDINATES OF THE NODES.

!        F = ARRAY OF LENGTH N CONTAINING THE DATA VALUES
!            IN ONE-TO-ONE CORRESPONDENCE WITH THE NODES.

!        NQ = NUMBER OF DATA POINTS TO BE USED IN THE LEAST SQUARES FIT FOR
!             COEFFICIENTS DEFINING THE NODAL FUNCTIONS Q(K).
!             A HIGHLY RECOMMENDED VALUE IS NQ = 13.
!             5 <= NQ <= MIN(40,N-1).

!        NW = NUMBER OF NODES WITHIN (AND DEFINING) THE RADII OF INFLUENCE
!             R(K) WHICH ENTER INTO THE WEIGHTS W(K).  FOR N SUFFICIENTLY
!             LARGE, A RECOMMENDED VALUE IS NW = 19.
!             1 <= NW <= MIN(40,N-1).

!        NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID DEFINED IN
!             SUBROUTINE STORE2.  A RECTANGLE CONTAINING THE NODES IS
!             PARTITIONED INTO CELLS IN ORDER TO INCREASE SEARCH EFFICIENCY.
!             NR = SQRT(N/3) IS RECOMMENDED.
!             NR >= 1.

! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

!        LCELL = ARRAY OF LENGTH >= NR**2.

!        LNEXT = ARRAY OF LENGTH >= N.

!        RSQ = ARRAY OF LENGTH >= N.

!        A = ARRAY OF LENGTH >= 5N.

! ON OUTPUT --

!        LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED
!                WITH CELLS.  REFER TO STORE2.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!                CES.  REFER TO STORE2.

!        XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL
!                          DIMENSIONS.  REFER TO STORE2.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!               MAXIMUM RADIUS R(K).

!        RSQ = ARRAY CONTAINING THE SQUARES OF THE RADII R(K)
!              WHICH ENTER INTO THE WEIGHTS W(K).

!        A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR
!            QUADRATIC NODAL FUNCTION Q(K) IN COLUMN K.

!   NOTE THAT THE ABOVE OUTPUT PARAMETERS ARE NOT DEFINED
! UNLESS IER = 0.

!        IER = ERROR INDICATOR --
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N, NQ, NW, OR NR IS OUT OF RANGE.
!              IER = 2 IF DUPLICATE NODES WERE ENCOUNTERED.
!              IER = 3 IF ALL NODES ARE COLLINEAR.

! MODULES REQUIRED BY QSHEP2 -- GETNP2, GIVENS, ROTATE,
!                                  SETUP2, STORE2

! INTRINSIC FUNCTIONS CALLED BY QSHEP2 -- ABS, MIN, FLOAT,
!                                            MAX0, MIN0, SQRT

!***********************************************************

INTEGER :: i, ib, ierr, ip1, irm1, irow, j, jp1, k, lmax, lnp, neq,  &
           nn, nnq, nnr, nnw, np, npts(40), nqwmax
REAL :: av, avsq, b(6,6), c, ddx, ddy, dmin, fk, rq, rs, rsmx,  &
        rsold, rws, s, sum, t, xk, xmn, yk, ymn

REAL, PARAMETER  :: dtol = 0.01, rtol = 1.0e-5, sf = 1.0

! LOCAL PARAMETERS --

! AV =               ROOT-MEAN-SQUARE DISTANCE BETWEEN K AND THE NODES IN THE
!                 LEAST SQUARES SYSTEM (UNLESS ADDITIONAL NODES ARE
!                 INTRODUCED FOR STABILITY).
!                 THE FIRST 3 COLUMNS OF THE MATRIX ARE SCALED BY 1/AVSQ,
!                 THE LAST 2 BY 1/AV
! AVSQ =       AV*AV
! B =               TRANSPOSE OF THE AUGMENTED REGRESSION MATRIX
! C =               FIRST COMPONENT OF THE PLANE ROTATION USED TO ZERO THE
!                 LOWER TRIANGLE OF B**T -- COMPUTED BY SUBROUTINE GIVENS
! DDX,DDY =    LOCAL VARIABLES FOR DX AND DY
! DMIN =       MINIMUM OF THE MAGNITUDES OF THE DIAGONAL ELEMENTS OF THE
!                 REGRESSION MATRIX AFTER ZEROS ARE INTRODUCED BELOW THE
!                 DIAGONAL
! DTOL =       TOLERANCE FOR DETECTING AN ILL-CONDITIONED SYSTEM.
!                 THE SYSTEM IS ACCEPTED WHEN DMIN >= DTOL
! FK =               DATA VALUE AT NODE K -- F(K)
! I =               INDEX FOR A, B, AND NPTS
! IB =               DO-LOOP INDEX FOR BACK SOLVE
! IERR =       ERROR FLAG FOR THE CALL TO STORE2
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       ROW INDEX FOR B
! J =               INDEX FOR A AND B
! JP1 =        J+1
! K =               NODAL FUNCTION INDEX AND COLUMN INDEX FOR A
! LMAX =       MAXIMUM NUMBER OF NPTS ELEMENTS (MUST BE CONSISTENT WITH
!                 THE DIMENSION STATEMENT ABOVE)
! LNP =        CURRENT LENGTH OF NPTS
! NEQ =        NUMBER OF EQUATIONS IN THE LEAST SQUARES FIT
! NN,NNQ,NNR = LOCAL COPIES OF N, NQ, AND NR
! NNW =        LOCAL COPY OF NW
! NP =               NPTS ELEMENT
! NPTS =       ARRAY CONTAINING THE INDICES OF A SEQUENCE OF NODES TO BE USED
!                 IN THE LEAST SQUARES FIT OR TO COMPUTE RSQ.  THE NODES ARE
!                 ORDERED BY DISTANCE FROM K AND THE LAST ELEMENT
!                 (USUALLY INDEXED BY LNP) IS USED ONLY TO DETERMINE RQ,
!                 OR RSQ(K) IF NW > NQ
! NQWMAX =     MAX(NQ,NW)
! RQ =               RADIUS OF INFLUENCE WHICH ENTERS INTO THE
!                 WEIGHTS FOR Q(K) (SEE SUBROUTINE SETUP2)
! RS =               SQUARED DISTANCE BETWEEN K AND NPTS(LNP) --
!                 USED TO COMPUTE RQ AND RSQ(K)
! RSMX =       MAXIMUM RSQ ELEMENT ENCOUNTERED
! RSOLD =      SQUARED DISTANCE BETWEEN K AND NPTS(LNP-1) --
!                 USED TO COMPUTE A RELATIVE CHANGE IN RS
!                 BETWEEN SUCCEEDING NPTS ELEMENTS
! RTOL =       TOLERANCE FOR DETECTING A SUFFICIENTLY LARGE
!                 RELATIVE CHANGE IN RS.  IF THE CHANGE IS
!                 NOT GREATER THAN RTOL, THE NODES ARE
!                 TREATED AS BEING THE SAME DISTANCE FROM K
! RWS =        CURRENT VALUE OF RSQ(K)
! S =               SECOND COMPONENT OF THE PLANE GIVENS ROTATION
! SF =               MARQUARDT STABILIZATION FACTOR USED TO DAMP
!                 OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!                 PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!                 IS ILL-CONDITIONED.  AS SF INCREASES, THE
!                 FITTING FUNCTION APPROACHES A LINEAR
! SUM =        SUM OF SQUARED EUCLIDEAN DISTANCES BETWEEN
!                 NODE K AND THE NODES USED IN THE LEAST
!                 SQUARES FIT (UNLESS ADDITIONAL NODES ARE
!                 ADDED FOR STABILITY)
! T =               TEMPORARY VARIABLE FOR ACCUMULATING A SCALAR
!                 PRODUCT IN THE BACK SOLVE
! XK,YK =      COORDINATES OF NODE K -- X(K), Y(K)
! XMN,YMN =    LOCAL VARIABLES FOR XMIN AND YMIN

nn = n
nnq = nq
nnw = nw
nnr = nr
nqwmax = MAX(nnq,nnw)
lmax = MIN(40,nn-1)
IF (5 <= nnq .AND. 1 <= nnw .AND. nqwmax <= lmax .AND. nnr >= 1) THEN
  
! CREATE THE CELL DATA STRUCTURE, AND INITIALIZE RSMX.
  
  CALL store2(nn,x,y,nnr,lcell,lnext,xmn,ymn,ddx,ddy,ierr)
  IF (ierr /= 0) GO TO 150
  rsmx = 0.
  
! OUTER LOOP ON NODE K
  
  DO  k = 1, nn
    xk = x(k)
    yk = y(k)
    fk = f(k)
    
! MARK NODE K TO EXCLUDE IT FROM THE SEARCH FOR NEAREST
!   NEIGHBORS.
    
    lnext(k) = -lnext(k)
    
! INITIALIZE FOR LOOP ON NPTS.
    
    rs = 0.
    sum = 0.
    rws = 0.
    rq = 0.
    lnp = 0
    
! COMPUTE NPTS, LNP, RWS, NEQ, RQ, AND AVSQ.
    
    10 sum = sum + rs
    IF (lnp /= lmax) THEN
      lnp = lnp + 1
      rsold = rs
      CALL getnp2(xk,yk,x,y,nnr,lcell,lnext,xmn,ymn,ddx,ddy,np,rs)
      IF (rs == 0.) GO TO 140
      npts(lnp) = np
      IF ((rs-rsold)/rs < rtol) GO TO 10
      IF (rws == 0. .AND. lnp > nnw) rws = rs
      IF (rq == 0. .AND. lnp > nnq) THEN
        
!   RQ = 0 (NOT YET COMPUTED) AND LNP > NQ.        RQ =
!     SQRT(RS) IS SUFFICIENTLY LARGE TO (STRICTLY) INCLUDE
!     NQ NODES.  THE LEAST SQUARES FIT WILL INCLUDE NEQ =
!     LNP - 1 EQUATIONS FOR 5 <= NQ <= NEQ < LMAX
!     <= N-1.
        
        neq = lnp - 1
        rq = SQRT(rs)
        avsq = sum / REAL(neq)
      END IF
      
!   BOTTOM OF LOOP -- TEST FOR TERMINATION.
      
      IF (lnp > nqwmax) GO TO 20
      GO TO 10
    END IF
    
! ALL LMAX NODES ARE INCLUDED IN NPTS.        RWS AND/OR RQ**2 IS
!   (ARBITRARILY) TAKEN TO BE 10 PERCENT LARGER THAN THE
!   DISTANCE RS TO THE LAST NODE INCLUDED.
    
    IF (rws == 0.) rws = 1.1 * rs
    IF (rq == 0.) THEN
      neq = lmax
      rq = SQRT(1.1*rs)
      avsq = sum / REAL(neq)
    END IF
    
! STORE RSQ(K), UPDATE RSMX IF NECESSARY, AND COMPUTE AV.
    
    20 rsq(k) = rws
    IF (rws > rsmx) rsmx = rws
    av = SQRT(avsq)
    
! SET UP THE AUGMENTED REGRESSION MATRIX (TRANSPOSED) AS THE COLUMNS OF B,
!   AND ZERO OUT THE LOWER TRIANGLE (UPPER TRIANGLE OF B) WITH GIVENS
!   ROTATIONS -- QR DECOMPOSITION WITH ORTHOGONAL MATRIX Q NOT STORED.
    
    i = 0
    30 i = i + 1
    np = npts(i)
    irow = MIN(i,6)
    CALL setup2(xk,yk,fk,x(np),y(np),f(np),av,avsq,rq,b(1,irow))
    IF (i == 1) GO TO 30
    irm1 = irow - 1
    DO  j = 1, irm1
      jp1 = j + 1
      CALL givens(b(j,j),b(j,irow),c,s)
      CALL rotate(6-j,c,s,b(jp1:,j),b(jp1:,irow))
    END DO
    IF (i < neq) GO TO 30
    
! TEST THE SYSTEM FOR ILL-CONDITIONING.
    
    dmin = MIN(ABS(b(1,1)),ABS(b(2,2)),ABS(b(3,3)),ABS(b(4,4)), ABS(b(5,5)))
    IF (dmin*rq < dtol) THEN
      IF (neq /= lmax) THEN
        
! INCREASE RQ AND ADD ANOTHER EQUATION TO THE SYSTEM TO IMPROVE THE
!   CONDITIONING.  THE NUMBER OF NPTS ELEMENTS IS ALSO INCREASED IF NECESSARY.
        
        50 rsold = rs
        neq = neq + 1
        IF (neq /= lmax) THEN
          IF (neq /= lnp) THEN
            
!   NEQ < LNP
            
            np = npts(neq+1)
            rs = (x(np)-xk) ** 2 + (y(np)-yk) ** 2
            IF ((rs-rsold)/rs < rtol) GO TO 50
            rq = SQRT(rs)
            GO TO 30
          END IF
          
!   ADD AN ELEMENT TO NPTS.
          
          lnp = lnp + 1
          CALL getnp2(xk,yk,x,y,nnr,lcell,lnext,xmn,ymn,ddx,ddy, np,rs)
          IF (np == 0) GO TO 140
          npts(lnp) = np
          IF ((rs-rsold)/rs < rtol) GO TO 50
          rq = SQRT(rs)
          GO TO 30
        END IF
        
        rq = SQRT(1.1*rs)
        GO TO 30
      END IF
      
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS -- ADD MULTIPLES OF THE
!   FIRST THREE UNIT VECTORS TO THE FIRST THREE EQUATIONS.
      
      DO  i = 1, 3
        b(i,6) = sf
        ip1 = i + 1
        DO  j = ip1, 6
          b(j,6) = 0.
        END DO
        DO  j = i, 5
          jp1 = j + 1
          CALL givens(b(j,j),b(j,6),c,s)
          CALL rotate(6-j,c,s,b(jp1:,j),b(jp1:,6))
        END DO
      END DO
      
! TEST THE STABILIZED SYSTEM FOR ILL-CONDITIONING.
      
      dmin = MIN(ABS(b(1,1)),ABS(b(2,2)),ABS(b(3,3)),  &
             ABS(b(4,4)),ABS(b(5,5)))
      IF (dmin*rq < dtol) GO TO 150
    END IF
    
! SOLVE THE 5 BY 5 TRIANGULAR SYSTEM FOR THE COEFFICIENTS
    
    DO  ib = 1, 5
      i = 6 - ib
      t = 0.
      IF (i /= 5) THEN
        ip1 = i + 1
        DO  j = ip1, 5
          t = t + b(j,i) * a(j,k)
        END DO
      END IF
      a(i,k) = (b(6,i)-t) / b(i,i)
    END DO
    
! SCALE THE COEFFICIENTS TO ADJUST FOR THE COLUMN SCALING.
    
    DO  i = 1, 3
      a(i,k) = a(i,k) / avsq
    END DO
    a(4,k) = a(4,k) / av
    a(5,k) = a(5,k) / av
    
! UNMARK K AND THE ELEMENTS OF NPTS.
    
    lnext(k) = -lnext(k)
    DO  i = 1, lnp
      np = npts(i)
      lnext(np) = -lnext(np)
    END DO
  END DO
  
! NO ERRORS ENCOUNTERED.
  
  xmin = xmn
  ymin = ymn
  dx = ddx
  dy = ddy
  rmax = SQRT(rsmx)
  ier = 0
  RETURN
END IF

! N, NQ, NW, OR NR IS OUT OF RANGE.

ier = 1
RETURN

! DUPLICATE NODES WERE ENCOUNTERED BY GETNP2.

140 ier = 2
RETURN

! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES.

150 xmin = xmn
ymin = ymn
dx = ddx
dy = ddy
ier = 3

RETURN
END SUBROUTINE qshep2



FUNCTION qs2val(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a) &
         RESULT(fn_val)
INTEGER, INTENT(IN)  :: n, nr, lcell(:,:), lnext(:)
REAL, INTENT(IN)     :: px, py, x(:), y(:), f(:), xmin, ymin, dx, dy,  &
                        rmax, rsq(:), a(:,:)
REAL                 :: fn_val

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    10/28/87

!   THIS FUNCTION RETURNS THE VALUE Q(PX,PY) WHERE Q IS THE WEIGHTED SUM
! OF QUADRATIC NODAL FUNCTIONS DEFINED IN SUBROUTINE QSHEP2.
! QS2GRD MAY BE CALLED TO COMPUTE A GRADIENT OF Q ALONG WITH THE VALUE,
! AND/OR TO TEST FOR ERRORS.

! ON INPUT --

!        PX,PY = CARTESIAN COORDINATES OF THE POINT P AT
!                WHICH Q IS TO BE EVALUATED.

!        N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!            N >= 6.

!        X,Y,F = ARRAYS OF LENGTH N CONTAINING THE NODES AND
!                DATA VALUES INTERPOLATED BY Q.

!        NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!             REFER TO STORE2.  NR >= 1.

!        LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
!                REFER TO STORE2.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES.
!                REFER TO STORE2.

!        XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.
!                          DX AND DY MUST BE POSITIVE.  REFER TO STORE2.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ -- MAXIMUM RADIUS.

!        RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!              WHICH ENTER INTO THE WEIGHTS DEFINING Q.

!        A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!            NODAL FUNCTIONS DEFINING Q.

!   INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.   THE PARAMETERS OTHER
! THAN PX AND PY SHOULD BE INPUT UNALTERED FROM THEIR VALUES ON OUTPUT FROM
! QSHEP2.  THIS FUNCTION SHOULD NOT BE CALLED IF A NONZERO ERROR FLAG WAS
! RETURNED BY QSHEP2.

! ON OUTPUT --

!        QS2VAL = FUNCTION VALUE Q(PX,PY) UNLESS N, NR, DX, DY, OR RMAX IS
!                 INVALID, IN WHICH CASE NO VALUE IS RETURNED.

! MODULES REQUIRED BY QS2VAL -- NONE

! INTRINSIC FUNCTIONS CALLED BY QS2VAL -- IFIX, SQRT

!***********************************************************

! Local variables
REAL     :: delx, dely, ds, dxsq, dysq, rd, rds, rs, sw, swq, w, xp, yp
INTEGER  :: i, imax, imin, j, jmax, jmin, k, kp

xp = px
yp = py
IF (n < 6 .OR. nr < 1 .OR. dx <= 0. .OR. dy <= 0. .OR. rmax < 0.) RETURN

! SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING THE RANGE OF THE
!   SEARCH FOR NODES WHOSE RADII INCLUDE P.
!   THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTERSECTED BY (OR CONTAINED
!   IN) A CIRCLE OF RADIUS RMAX CENTERED AT P.

imin = (xp-xmin-rmax)/dx + 1
imax = (xp-xmin+rmax)/dx + 1
IF (imin < 1) imin = 1
IF (imax > nr) imax = nr
jmin = (yp-ymin-rmax)/dy + 1
jmax = (yp-ymin+rmax)/dy + 1
IF (jmin < 1) jmin = 1
IF (jmax > nr) jmax = nr

! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE OF RADIUS RMAX.

IF (imin <= imax .AND. jmin <= jmax) THEN
  
! ACCUMULATE WEIGHT VALUES IN SW AND WEIGHTED NODAL FUNCTION VALUES IN SWQ.
!   THE WEIGHTS ARE W(K) = ((R-D)+/(R*D))**2 FOR R**2 = RSQ(K) AND
!   D = DISTANCE BETWEEN P AND NODE K.
  
  sw = 0.
  swq = 0.
  
! OUTER LOOP ON CELLS (I,J).
  
  DO  j = jmin, jmax
    DO  i = imin, imax
      k = lcell(i,j)
      IF (k /= 0) THEN
        
! INNER LOOP ON NODES K.
        
        10 delx = xp - x(k)
        dely = yp - y(k)
        dxsq = delx * delx
        dysq = dely * dely
        ds = dxsq + dysq
        rs = rsq(k)
        IF (ds < rs) THEN
          IF (ds == 0.) GO TO 40
          rds = rs * ds
          rd = SQRT(rds)
          w = (rs+ds-rd-rd) / rds
          sw = sw + w
          swq = swq + w * (a(1,k)*dxsq + a(2,k)*delx*dely + a(3,k)*dysq +  &
                a(4,k)*delx + a(5,k)*dely + f(k))
        END IF
        
! BOTTOM OF LOOP ON NODES IN CELL (I,J).
        
        kp = k
        k = lnext(kp)
        IF (k /= kp) GO TO 10
      END IF
    END DO
  END DO
  
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.
  
  IF (sw == 0.) GO TO 50
  fn_val = swq / sw
  RETURN
  
! (PX,PY) = (X(K),Y(K))
  
  40 fn_val = f(k)
  RETURN
END IF

! ALL WEIGHTS ARE 0 AT P.

50 fn_val = 0.

RETURN
END FUNCTION qs2val



SUBROUTINE qs2grd(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,  &
                  rmax,rsq,a,q,qx,qy,ier)
INTEGER, INTENT(IN)   :: n, nr, lcell(:,:), lnext(:)
REAL, INTENT(IN)      :: px, py, x(:), y(:), f(:), xmin, ymin, dx, dy,  &
                         rmax, rsq(:), a(:,:)
REAL, INTENT(OUT)     :: q, qx, qy
INTEGER, INTENT(OUT)  :: ier

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    10/28/87

!   THIS SUBROUTINE COMPUTES THE VALUE AND GRADIENT AT (PX,PY) OF THE
! INTERPOLATORY FUNCTION Q DEFINED IN SUBROUTINE QSHEP2.
! Q(X,Y) IS A WEIGHTED SUM OF QUADRATIC NODAL FUNCTIONS.

! ON INPUT --

!        PX,PY = CARTESIAN COORDINATES OF THE POINT AT WHICH
!                Q AND ITS PARTIALS ARE TO BE EVALUATED.

!        N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!            N >= 6.

!        X,Y,F = ARRAYS OF LENGTH N CONTAINING THE NODES AND
!                DATA VALUES INTERPOLATED BY Q.

!        NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!             REFER TO STORE2.  NR >= 1.

!        LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
!                REFER TO STORE2.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES.
!                REFER TO STORE2.

!        XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.
!                          DX AND DY MUST BE POSITIVE.  REFER TO STORE2.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ -- MAXIMUM RADIUS.

!        RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!              WHICH ENTER INTO THE WEIGHTS DEFINING Q.

!        A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!            NODAL FUNCTIONS DEFINING Q.

!   INPUT PARAMETERS ARE NOT ALTERED BY THIS SUBROUTINE.
! THE PARAMETERS OTHER THAN PX AND PY SHOULD BE INPUT UNALTERED FROM THEIR
! VALUES ON OUTPUT FROM QSHEP2.  THIS SUBROUTINE SHOULD NOT BE CALLED IF A
! NONZERO ERROR FLAG WAS RETURNED BY QSHEP2.

! ON OUTPUT --

!        Q = VALUE OF Q AT (PX,PY) UNLESS IER .EQ. 1, IN
!            WHICH CASE NO VALUES ARE RETURNED.

!        QX,QY = FIRST PARTIAL DERIVATIVES OF Q AT (PX,PY) UNLESS IER .EQ. 1.

!        IER = ERROR INDICATOR
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N, NR, DX, DY OR RMAX IS INVALID.
!              IER = 2 IF NO ERRORS WERE ENCOUNTERED BUT (PX,PY) IS NOT WITHIN
!                      THE RADIUS R(K) FOR ANY NODE K (AND THUS Q=QX=QY=0).

! MODULES REQUIRED BY QS2GRD -- NONE

! INTRINSIC FUNCTIONS CALLED BY QS2GRD -- IFIX, SQRT

!***********************************************************

! Local variables
REAL     :: delx, dely, ds, dxsq, dysq, qk, qkx, qky, rd, rds, rs,  &
            sw, swq, swqx, swqy, sws, swx, swy, t, w, wx, wy, xp, yp
INTEGER  :: i, imax, imin, j, jmax, jmin, k, kp

xp = px
yp = py
IF (n >= 6 .AND. nr >= 1 .AND. dx > 0. .AND. dy > 0. .AND. rmax >= 0.) THEN
  
! SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING P.
!   THE RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE THE CELLS WHICH MUST
!   BE SEARCHED ARE THOSE INTERSECTED BY (OR CONTAINED IN) A CIRCLE OF RADIUS
!   RMAX CENTERED AT P.
  
  imin = (xp-xmin-rmax)/dx + 1
  imax = (xp-xmin+rmax)/dx + 1
  IF (imin < 1) imin = 1
  IF (imax > nr) imax = nr
  jmin = (yp-ymin-rmax)/dy + 1
  jmax = (yp-ymin+rmax)/dy + 1
  IF (jmin < 1) jmin = 1
  IF (jmax > nr) jmax = nr
  
! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE OF RADIUS RMAX.
  
  IF (imin > imax .OR. jmin > jmax) GO TO 50
  
! Q = SWQ/SW = SUM(W(K)*Q(K))/SUM(W(K)) WHERE THE SUM IS FROM K = 1 TO N,
!   Q(K) IS THE QUADRATIC NODAL FUNCTION, AND W(K) = ((R-D)+/(R*D))**2 FOR
!   RADIUS R(K) AND DISTANCE D(K).
!   THUS
  
!         QX = (SWQX*SW - SWQ*SWX)/SW**2  AND
!         QY = (SWQY*SW - SWQ*SWY)/SW**2
  
!   WHERE SWQX AND SWX ARE PARTIAL DERIVATIVES WITH RESPECT TO X OF SWQ
!   AND SW, RESPECTIVELY.  SWQY AND SWY ARE DEFINED SIMILARLY.
  
  sw = 0.
  swx = 0.
  swy = 0.
  swq = 0.
  swqx = 0.
  swqy = 0.
  
! OUTER LOOP ON CELLS (I,J).
  
  DO  j = jmin, jmax
    DO  i = imin, imax
      k = lcell(i,j)
      IF (k /= 0) THEN
        
! INNER LOOP ON NODES K.
        
        10 delx = xp - x(k)
        dely = yp - y(k)
        dxsq = delx * delx
        dysq = dely * dely
        ds = dxsq + dysq
        rs = rsq(k)
        IF (ds < rs) THEN
          IF (ds == 0.) GO TO 40
          rds = rs * ds
          rd = SQRT(rds)
          w = (rs+ds-rd-rd) / rds
          t = 2. * (rd-rs) / (ds*rds)
          wx = delx * t
          wy = dely * t
          qkx = 2. * a(1,k) * delx + a(2,k) * dely
          qky = a(2,k) * delx + 2. * a(3,k) * dely
          qk = (qkx*delx+qky*dely) / 2.
          qkx = qkx + a(4,k)
          qky = qky + a(5,k)
          qk = qk + a(4,k) * delx + a(5,k) * dely + f(k)
          sw = sw + w
          swx = swx + wx
          swy = swy + wy
          swq = swq + w * qk
          swqx = swqx + wx * qk + w * qkx
          swqy = swqy + wy * qk + w * qky
        END IF
        
! BOTTOM OF LOOP ON NODES IN CELL (I,J).
        
        kp = k
        k = lnext(kp)
        IF (k /= kp) GO TO 10
      END IF
    END DO
  END DO
  
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.
  
  IF (sw == 0.) GO TO 50
  q = swq / sw
  sws = sw * sw
  qx = (swqx*sw-swq*swx) / sws
  qy = (swqy*sw-swq*swy) / sws
  ier = 0
  RETURN
  
! (PX,PY) = (X(K),Y(K))
  
  40 q = f(k)
  qx = a(4,k)
  qy = a(5,k)
  ier = 0
  RETURN
END IF

! INVALID INPUT PARAMETER.

ier = 1
RETURN

! NO CELLS CONTAIN A POINT WITHIN RMAX OF P, OR
!   SW = 0 AND THUS DS >= RSQ(K) FOR ALL K.

50 q = 0.
qx = 0.
qy = 0.
ier = 2

RETURN
END SUBROUTINE qs2grd



SUBROUTINE getnp2(px,py,x,y,nr,lcell,lnext,xmin,ymin,dx,dy,np,dsq)
INTEGER, INTENT(IN)      :: nr, lcell(:,:)
INTEGER, INTENT(IN OUT)  :: lnext(:)
REAL, INTENT(IN)         :: px, py, x(1), y(1), xmin, ymin, dx, dy
INTEGER, INTENT(OUT)     :: np
REAL, INTENT(OUT)        :: dsq

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   GIVEN A SET OF N NODES AND THE DATA STRUCTURE DEFINED IN SUBROUTINE STORE2,
! THIS SUBROUTINE USES THE CELL METHOD TO FIND THE CLOSEST UNMARKED NODE NP TO
! A SPECIFIED POINT P.
! NP IS THEN MARKED BY SETTING LNEXT(NP) TO -LNEXT(NP).  (A NODE IS MARKED IF
! AND ONLY IF THE CORRESPONDING LNEXT ELEMENT IS NEGATIVE.   THE ABSOLUTE
! VALUES OF LNEXT ELEMENTS, HOWEVER, MUST BE PRESERVED.)   THUS, THE CLOSEST
! M NODES TO P MAY BE DETERMINED BY A SEQUENCE OF M CALLS TO THIS ROUTINE.
! NOTE THAT IF THE NEAREST NEIGHBOR TO NODE K IS TO BE DETERMINED (PX = X(K)
! AND PY = Y(K)), THEN K SHOULD BE MARKED BEFORE THE CALL TO THIS ROUTINE.

!   THE SEARCH IS BEGUN IN THE CELL CONTAINING (OR CLOSEST TO) P AND PROCEEDS
! OUTWARD IN RECTANGULAR LAYERS UNTIL ALL CELLS WHICH CONTAIN POINTS WITHIN
! DISTANCE R OF P HAVE BEEN SEARCHED, WHERE R IS THE DISTANCE FROM P TO THE
! FIRST UNMARKED NODE ENCOUNTERED (INFINITE IF NO UNMARKED NODES ARE PRESENT).

! ON INPUT --

!        PX,PY = CARTESIAN COORDINATES OF THE POINT P WHOSE
!                NEAREST UNMARKED NEIGHBOR IS TO BE FOUND.

!        X,Y = ARRAYS OF LENGTH N, FOR N >= 2, CONTAINING
!              THE CARTESIAN COORDINATES OF THE NODES.

!        NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!             NR >= 1.

!        LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES
!                (OR THEIR NEGATIVES).

!        XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.
!                          DX AND DY MUST BE POSITIVE.

!   INPUT PARAMETERS OTHER THAN LNEXT ARE NOT ALTERED BY THIS ROUTINE. WITH
! THE EXCEPTION OF (PX,PY) AND THE SIGNS OF LNEXT ELEMENTS, THESE PARAMETERS
! SHOULD BE UNALTERED FROM THEIR VALUES ON OUTPUT FROM SUBROUTINE STORE2.

! ON OUTPUT --

!        NP = INDEX (FOR X AND Y) OF THE NEAREST UNMARKED NODE TO P,
!             OR 0 IF ALL NODES ARE MARKED OR NR < 1 OR DX <= 0 OR
!             DY <= 0.
!             LNEXT(NP) < 0 IF NP .NE. 0.

!        DSQ = SQUARED EUCLIDEAN DISTANCE BETWEEN P AND NODE
!              NP, OR 0 IF NP = 0.

! MODULES REQUIRED BY GETNP2 -- NONE

! INTRINSIC FUNCTIONS CALLED BY GETNP2 -- ABS, IFIX, SQRT

!***********************************************************

! Local variables
LOGICAL  :: first
REAL     :: delx, dely, r, rsmin, rsq, xp, yp
INTEGER  :: i, i0, i1, i2, imax, imin, j, j0, j1, j2, jmax, jmin, l, lmin, ln

xp = px
yp = py

! TEST FOR INVALID INPUT PARAMETERS.

IF (nr >= 1 .AND. dx > 0. .AND. dy > 0.) THEN
  
! INITIALIZE PARAMETERS --
  
!   FIRST = TRUE IFF THE FIRST UNMARKED NODE HAS YET TO BE ENCOUNTERED,
!   IMIN,IMAX,JMIN,JMAX = CELL INDICES DEFINING THE RANGE OF THE SEARCH,
!   DELX,DELY = PX-XMIN AND PY-YMIN,
!   I0,J0 = CELL CONTAINING OR CLOSEST TO P,
!   I1,I2,J1,J2 = CELL INDICES OF THE LAYER WHOSE INTERSECTION WITH THE RANGE
!                  DEFINED BY IMIN,..., JMAX IS CURRENTLY BEING SEARCHED.
  
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  delx = xp - xmin
  dely = yp - ymin
  i0 = delx/dx + 1
  IF (i0 < 1) i0 = 1
  IF (i0 > nr) i0 = nr
  j0 = dely/dy + 1
  IF (j0 < 1) j0 = 1
  IF (j0 > nr) j0 = nr
  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
  
! OUTER LOOP ON LAYERS, INNER LOOP ON LAYER CELLS, EXCLUDING
!   THOSE OUTSIDE THE RANGE (IMIN,IMAX) X (JMIN,JMAX).

  10 i = 0
  loop40:  DO  j = j1, j2
    IF (j > jmax) GO TO 50
    IF (j >= jmin) THEN
      DO  i = i1, i2
        IF (i > imax) CYCLE loop40
        IF (i >= imin) THEN
          IF (j == j1 .OR. j == j2 .OR. i == i1 .OR. i == i2) THEN
            
! SEARCH CELL (I,J) FOR UNMARKED NODES L.
            
            l = lcell(i,j)
            IF (l /= 0) THEN
              
!   LOOP ON NODES IN CELL (I,J).
              
              20 ln = lnext(l)
              IF (ln >= 0) THEN
                
!   NODE L IS NOT MARKED.
                
                rsq = (x(l)-xp) ** 2 + (y(l)-yp) ** 2
                IF (first) THEN
                  
!   NODE L IS THE FIRST UNMARKED NEIGHBOR OF P ENCOUNTERED.
!     INITIALIZE LMIN TO THE CURRENT CANDIDATE FOR NP, AND RSMIN TO THE
!     SQUARED DISTANCE FROM P TO LMIN.  IMIN, IMAX, JMIN, AND JMAX ARE UPDATED
!     TO DEFINE THE SMALLEST RECTANGLE CONTAINING A CIRCLE OF RADIUS R =
!     SQRT(RSMIN) CENTERED AT P, AND CONTAINED IN (1,NR) X (1,NR) (EXCEPT
!     THAT, IF P IS OUTSIDE THE RECTANGLE DEFINED BY THE NODES, IT IS POSSIBLE
!     THAT IMIN > NR, IMAX < 1, JMIN > NR, OR JMAX < 1).
!     FIRST IS RESET TO FALSE.
                  
                  lmin = l
                  rsmin = rsq
                  r = SQRT(rsmin)
                  imin = (delx-r)/dx + 1
                  IF (imin < 1) imin = 1
                  imax = (delx+r)/dx + 1
                  IF (imax > nr) imax = nr
                  jmin = (dely-r)/dy + 1
                  IF (jmin < 1) jmin = 1
                  jmax = (dely+r)/dy + 1
                  IF (jmax > nr) jmax = nr
                  first = .false.
                ELSE
                  
!   TEST FOR NODE L CLOSER THAN LMIN TO P.
                  
                  IF (rsq < rsmin) THEN
                    
!   UPDATE LMIN AND RSMIN.
                    
                    lmin = l
                    rsmin = rsq
                  END IF
                END IF
              END IF
              
!   TEST FOR TERMINATION OF LOOP ON NODES IN CELL (I,J).
              
              IF (ABS(ln) /= l) THEN
                l = ABS(ln)
                GO TO 20
              END IF
            END IF
          END IF
        END IF
      END DO
    END IF
  END DO loop40
  
! TEST FOR TERMINATION OF LOOP ON CELL LAYERS.
  
  50 IF (i1 > imin .OR. i2 < imax .OR. j1 > jmin .OR. j2 < jmax) THEN
    i1 = i1 - 1
    i2 = i2 + 1
    j1 = j1 - 1
    j2 = j2 + 1
    GO TO 10
  END IF
  
! UNLESS NO UNMARKED NODES WERE ENCOUNTERED, LMIN IS THE
!   CLOSEST UNMARKED NODE TO P.
  
  IF (.NOT.first) THEN
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
    RETURN
  END IF
END IF

! ERROR -- NR, DX, OR DY IS INVALID OR ALL NODES ARE MARKED.

np = 0
dsq = 0.

RETURN
END SUBROUTINE getnp2



SUBROUTINE givens(a,b,c,s)
REAL, INTENT(IN OUT)  :: a, b
REAL, INTENT(OUT)     :: c, s

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --
!     ( C  S)
! G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND ENTRY OF THE
!     (-S  C)
! 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS IS NORMALLY FOLLOWED BY A CALL
! TO ROTATE WHICH APPLIES THE TRANSFORMATION TO A 2 BY N MATRIX.
! THIS ROUTINE WAS TAKEN FROM LINPACK.

! ON INPUT --

!        A,B = COMPONENTS OF THE 2-VECTOR TO BE ROTATED.

! ON OUTPUT --

!        A = VALUE OVERWRITTEN BY R = +/-SQRT(A*A + B*B)

!        B = VALUE OVERWRITTEN BY A VALUE Z WHICH ALLOWS C
!            AND S TO BE RECOVERED AS FOLLOWS --
!              C = SQRT(1-Z*Z), S=Z     IF ABS(Z) <= 1.
!              C = 1/Z, S = SQRT(1-C*C) IF ABS(Z) > 1.

!        C = +/-(A/R)

!        S = +/-(B/R)

! MODULES REQUIRED BY GIVENS -- NONE

! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT

!***********************************************************

REAL  :: aa, bb, r, u, v

! LOCAL PARAMETERS --

! AA,BB = LOCAL COPIES OF A AND B
! R =          C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   VARIABLES USED TO SCALE A AND B FOR COMPUTING R

aa = a
bb = b
IF (ABS(aa) > ABS(bb)) THEN
  
! ABS(A) > ABS(B)
  
  u = aa + aa
  v = bb / u
  r = SQRT(.25+v*v) * u
  c = aa / r
  s = v * (c+c)
  
! NOTE THAT R HAS THE SIGN OF A, C > 0, AND S HAS SIGN(A)*SIGN(B).
  
  b = s
  a = r
  RETURN
END IF

! ABS(A) <= ABS(B)

IF (bb /= 0.) THEN
  u = bb + bb
  v = aa / u
  
! STORE R IN A.
  
  a = SQRT(.25+v*v) * u
  s = bb / a
  c = v * (s+s)
  
! NOTE THAT R HAS THE SIGN OF B, S > 0, AND C HAS SIGN(A)*SIGN(B).
  
  b = 1.
  IF (c /= 0.) b = 1. / c
  RETURN
END IF

! A = B = 0.

c = 1.
s = 0.

RETURN
END SUBROUTINE givens



SUBROUTINE rotate(n,c,s,x,y)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: c, s
REAL, INTENT(IN OUT)  :: x(:), y(:)

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!                                             ( C  S)
!   THIS ROUTINE APPLIES THE GIVENS ROTATION (           ) TO THE
!                                             (-S  C)
!                (X(1) ... X(N))
! 2 BY N MATRIX (               ).
!                (Y(1) ... Y(N))

! ON INPUT --

!        N = NUMBER OF COLUMNS TO BE ROTATED.

!        C,S = ELEMENTS OF THE GIVENS ROTATION.   THESE MAY BE DETERMINED
!              BY SUBROUTINE GIVENS.

!        X,Y = ARRAYS OF LENGTH >= N CONTAINING THE VECTORS TO BE ROTATED.

! PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS ROUTINE.

! ON OUTPUT --

!        X,Y = ROTATED VECTORS.

! MODULES REQUIRED BY ROTATE -- NONE

!***********************************************************

INTEGER  :: i
REAL     :: xi, yi

! LOCAL PARAMETERS --

! I =          DO-LOOP INDEX
! XI,YI = X(I), Y(I)

IF (n <= 0 .OR. (c == 1. .AND. s == 0.)) RETURN
DO  i = 1, n
  xi = x(i)
  yi = y(i)
  x(i) = c * xi + s * yi
  y(i) = -s * xi + c * yi
END DO

RETURN
END SUBROUTINE rotate



SUBROUTINE setup2(xk,yk,fk,xi,yi,fi,s1,s2,r,row)
REAL, INTENT(IN)   :: xk, yk, fk, xi, yi, fi, s1, s2, r
REAL, INTENT(OUT)  :: row(6)

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   THIS ROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED REGRESSION MATRIX FOR A
! WEIGHTED LEAST-SQUARES FIT OF A QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA
! VALUES F, WHERE Q(XK,YK) = FK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS) ARE
! SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LINEAR TERMS) ARE SCALED
! BY 1/S1.  THE WEIGHT IS (R-D)/(R*D) IF R > D AND 0 IF R <= D, WHERE D
! IS THE DISTANCE BETWEEN NODES I AND K.

! ON INPUT --

!        XK,YK,FK = COORDINATES AND DATA VALUE AT NODE K -- INTERPOLATED BY Q.

!        XI,YI,FI = COORDINATES AND DATA VALUE AT NODE I.

!        S1,S2 = RECIPROCALS OF THE SCALE FACTORS.

!        R = RADIUS OF INFLUENCE ABOUT NODE K DEFINING THE WEIGHT.

!        ROW = ARRAY OF LENGTH 6.

! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

! ON OUTPUT --

!        ROW = VECTOR CONTAINING A ROW OF THE AUGMENTED REGRESSION MATRIX.

! MODULES REQUIRED BY SETUP2 -- NONE

! INTRINSIC FUNCTION CALLED BY SETUP2 -- SQRT

!***********************************************************

REAL  :: d, dx, dy, dxsq, dysq, w, w1, w2

! LOCAL PARAMETERS -

! I =         DO-LOOP INDEX
! DX =         XI - XK
! DY =         YI - YK
! DXSQ = DX*DX
! DYSQ = DY*DY
! D =         DISTANCE BETWEEN NODES K AND I
! W =         WEIGHT ASSOCIATED WITH THE ROW
! W1 =         W/S1
! W2 =         W/S2

dx = xi - xk
dy = yi - yk
dxsq = dx * dx
dysq = dy * dy
d = SQRT(dxsq+dysq)
IF (d > 0. .AND. d < r) THEN
  w = (r-d) / r / d
  w1 = w / s1
  w2 = w / s2
  row(1) = dxsq * w2
  row(2) = dx * dy * w2
  row(3) = dysq * w2
  row(4) = dx * w1
  row(5) = dy * w1
  row(6) = (fi-fk) * w
  RETURN
END IF

! NODES K AND I COINCIDE OR NODE I IS OUTSIDE OF THE RADIUS OF INFLUENCE.
!   SET ROW TO THE ZERO VECTOR.

row(1:6) = 0.

RETURN
END SUBROUTINE setup2



SUBROUTINE store2(n,x,y,nr,lcell,lnext,xmin,ymin,dx,dy,ier)
INTEGER, INTENT(IN)      :: n, nr
INTEGER, INTENT(IN OUT)  :: lcell(:,:), lnext(:)
REAL, INTENT(IN)         :: x(:), y(:)
REAL, INTENT(OUT)        :: xmin, ymin, dx, dy
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   GIVEN A SET OF N ARBITRARILY DISTRIBUTED NODES IN THE PLANE, THIS
! SUBROUTINE CREATES A DATA STRUCTURE FOR A CELL-BASED METHOD OF SOLVING
! CLOSEST-POINT PROBLEMS.   THE SMALLEST RECTANGLE CONTAINING THE NODES IS
! PARTITIONED INTO AN NR BY NR UNIFORM GRID OF CELLS, AND NODES ARE
! ASSOCIATED WITH CELLS.   IN PARTICULAR, THE DATA STRUCTURE STORES THE
! INDICES OF THE NODES CONTAINED IN EACH CELL.
! FOR A UNIFORM RANDOM DISTRIBUTION OF NODES, THE NEAREST NODE TO AN
! ARBITRARY POINT CAN BE DETERMINED IN CONSTANT EXPECTED TIME.

! ON INPUT --

!        N = NUMBER OF NODES.  N >= 2.

!        X,Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!              COORDINATES OF THE NODES.

!        NR = NUMBER OF ROWS AND COLUMNS IN THE GRID.  THE CELL DENSITY
!             (AVERAGE NUMBER OF NODES PER CELL) IS D = N/(NR**2).
!             A RECOMMENDED VALUE, BASED ON EMPIRICAL EVIDENCE,
!             IS D = 3 -- NR = SQRT(N/3).  NR >= 1.

! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

!        LCELL = ARRAY OF LENGTH >= NR**2.

!        LNEXT = ARRAY OF LENGTH >= N.

! ON OUTPUT --

!        LCELL = NR BY NR CELL ARRAY SUCH THAT LCELL(I,J) CONTAINS THE INDEX
!                (FOR X AND Y) OF THE FIRST NODE (NODE WITH SMALLEST INDEX) IN
!                CELL (I,J), OR LCELL(I,J) = 0 IF NO NODES ARE CONTAINED IN
!                THE CELL.  THE UPPER RIGHT CORNER OF CELL (I,J) HAS
!                COORDINATES (XMIN+I*DX,YMIN+J*DY).
!                LCELL IS NOT DEFINED IF IER .NE. 0.

!        LNEXT = ARRAY OF NEXT-NODE INDICES SUCH THAT LNEXT(K) CONTAINS THE
!                INDEX OF THE NEXT NODE IN THE CELL WHICH CONTAINS NODE K, OR
!                LNEXT(K) = K IF K IS THE LAST NODE IN THE CELL FOR
!                K = 1,...,N.  (THE NODES CONTAINED IN A CELL ARE ORDERED BY
!                THEIR INDICES.)
!                IF, FOR EXAMPLE, CELL (I,J) CONTAINS NODES 2, 3, AND 5 (AND
!                NO OTHERS), THEN LCELL(I,J) = 2, LNEXT(2) = 3, LNEXT(3) = 5,
!                AND LNEXT(5) = 5.  LNEXT IS NOT DEFINED IF IER .NE. 0.

!        XMIN,YMIN = CARTESIAN COORDINATES OF THE LOWER LEFT CORNER OF THE
!                    RECTANGLE DEFINED BY THE NODES (SMALLEST NODAL
!                    COORDINATES) UNLESS IER = 1.  THE UPPER RIGHT CORNER IS
!                    (XMAX,YMAX) FOR XMAX = XMIN + NR*DX AND
!                    YMAX = YMIN + NR*DY.

!        DX,DY = DIMENSIONS OF THE CELLS UNLESS IER = 1.  DX = (XMAX-XMIN)/NR
!                AND DY = (YMAX-YMIN)/NR WHERE XMIN, XMAX, YMIN, AND YMAX ARE
!                THE EXTREMA OF X AND Y.

!        IER = ERROR INDICATOR --
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N < 2 OR NR < 1.
!              IER = 2 IF DX = 0 OR DY = 0.

! MODULES REQUIRED BY STORE2 -- NONE

! INTRINSIC FUNCTIONS CALLED BY STORE2 -- FLOAT, IFIX

!***********************************************************

! Local variables
INTEGER  :: i, j, k, kb, l, nn, nnr, np1
REAL     :: delx, dely, xmn, xmx, ymn, ymx

nn = n
nnr = nr
IF (nn >= 2 .AND. nnr >= 1) THEN
  
! COMPUTE THE DIMENSIONS OF THE RECTANGLE CONTAINING THE NODES.
  
  xmn = x(1)
  xmx = xmn
  ymn = y(1)
  ymx = ymn
  DO  k = 2, nn
    IF (x(k) < xmn) xmn = x(k)
    IF (x(k) > xmx) xmx = x(k)
    IF (y(k) < ymn) ymn = y(k)
    IF (y(k) > ymx) ymx = y(k)
  END DO
  xmin = xmn
  ymin = ymn
  
! COMPUTE CELL DIMENSIONS AND TEST FOR ZERO AREA.
  
  delx = (xmx-xmn) / REAL(nnr)
  dely = (ymx-ymn) / REAL(nnr)
  dx = delx
  dy = dely
  IF (delx == 0. .OR. dely == 0.) GO TO 50
  
! INITIALIZE LCELL.
  
  DO  j = 1, nnr
    DO  i = 1, nnr
      lcell(i,j) = 0
    END DO
  END DO
  
! LOOP ON NODES, STORING INDICES IN LCELL AND LNEXT.
  
  np1 = nn + 1
  DO  k = 1, nn
    kb = np1 - k
    i = (x(kb)-xmn)/delx + 1
    IF (i > nnr) i = nnr
    j = (y(kb)-ymn)/dely + 1
    IF (j > nnr) j = nnr
    l = lcell(i,j)
    lnext(kb) = l
    IF (l == 0) lnext(kb) = kb
    lcell(i,j) = kb
  END DO
  
! NO ERRORS ENCOUNTERED
  
  ier = 0
  RETURN
END IF

! INVALID INPUT PARAMETER

ier = 1
RETURN

! DX = 0 OR DY = 0

50 ier = 2
RETURN
END SUBROUTINE store2

END MODULE qshep2d



PROGRAM qs2test
!                           QS2TEST

!   THIS PROGRAM TESTS THE SCATTERED DATA INTERPOLATION PACKAGE QSHEP2D
! BY PRINTING THE MAXIMUM ERRORS ASSOCIATED WITH INTERPOLATED VALUES AND
! GRADIENTS ON A 10 BY 10 UNIFORM GRID IN THE UNIT SQUARE.  THE DATA SET
! CONSISTS OF 36 NODES WITH DATA VALUES TAKEN FROM A QUADRATIC FUNCTION FOR
! WHICH THE METHOD IS EXACT.  THE RATIO OF MAXIMUM INTERPOLATION ERROR
! RELATIVE TO THE MACHINE PRECISION IS ALSO PRINTED.  THIS SHOULD BE O(1).
! THE INTERPOLATED VALUES FROM QS2VAL AND QS2GRD ARE COMPARED FOR AGREEMENT.

USE qshep2d
IMPLICIT NONE
INTEGER :: i, ier, j, k, lcell(3,3), lnext(36)
REAL    :: a(5,36), dx, dy, eps, eq, eqx, eqy, f(36), p(10), px, py, q, q1, &
           qx, qy, rmax, rq, rsq(36), x(36), xmin, y(36), yk, ymin

! QSHEP2 PARAMETERS AND LOGICAL UNIT FOR OUTPUT

INTEGER, PARAMETER  :: lout = 6, n = 36, nq = 13, nr = 3, nw = 19

! GENERATE A 6 BY 6 GRID OF NODES IN THE UNIT SQUARE WITH THE NATURAL ORDERING.

k = 0
DO  j = 1, 6
  yk = (6-j) / 5.
  DO  i = 1, 6
    k = k + 1
    x(k) = (i-1) / 5.
    y(k) = yk
  END DO
END DO

! COMPUTE THE DATA VALUES.

DO  k = 1, n
  f(k) = fq(x(k),y(k))
END DO

! COMPUTE PARAMETERS DEFINING THE INTERPOLANT Q.

CALL qshep2(n,x,y,f,nq,nw,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a,ier)
IF (ier == 0) THEN

! GENERATE A 10 BY 10 UNIFORM GRID OF INTERPOLATION POINTS (P(I),P(J)) IN THE
!   UNIT SQUARE.  THE FOUR CORNERS COINCIDE WITH NODES.

  DO  i = 1, 10
    p(i) = (i-1) / 9.
  END DO

! COMPUTE THE MACHINE PRECISION EPS.

  eps = EPSILON(1.0)

! COMPUTE INTERPOLATION ERRORS AND TEST FOR AGREEMENT IN THE
!   Q VALUES RETURNED BY QS2VAL AND QS2GRD.

  EQ = 0.
  eqx = 0.
  eqy = 0.
  DO  j = 1, 10
    py = p(j)
    DO  i = 1, 10
      px = p(i)
      q1 = qs2val(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a)
      CALL qs2grd(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,  &
                  rmax,rsq,a,q,qx,qy,ier)
      IF (ier /= 0) GO TO 80
      IF (ABS(q1-q) > 3.*ABS(q)*eps) GO TO 90
      EQ = MAX(EQ,ABS(fq(px,py)-q))
      eqx = MAX(eqx,ABS(fx(px,py)-qx))
      eqy = MAX(eqy,ABS(fy(px,py)-qy))
    END DO
  END DO

! PRINT ERRORS AND THE RATIO EQ/EPS.

  rq = EQ / eps
  WRITE (lout,5000)
  WRITE (lout,5100) EQ, rq
  WRITE (lout,5200) eqx
  WRITE (lout,5300) eqy
  STOP
END IF

! ERROR IN QSHEP2

WRITE (lout,5400) ier
STOP

! ERROR IN QS2GRD

80 WRITE (lout,5500) ier
STOP

! VALUES RETURNED BY QS2VAL AND QS2GRD DIFFER BY A RELATIVE
!   AMOUNT GREATER THAN 3*EPS.

90 WRITE (lout,5600) q1, q
STOP

5000 FORMAT(//' ','MAXIMUM ABSOLUTE ERRORS IN THE INTERPOLANT Q AND PARTIAL'/ &
            ' DERIVATIVES QX AND QY RELATIVE TO MACHINE PRECISION EPS'//  &
            t12, 'FUNCTION   MAX ERROR   MAX ERROR/EPS'/)
5100 FORMAT (t15, 'Q       ', e9.3, '       ', f4.2)
5200 FORMAT (t15, 'QX      ', e9.3)
5300 FORMAT (t15, 'QY      ', e9.3)
5400 FORMAT (///' *** ERROR IN QSHEP2 -- IER =', i2, ' ***')
5500 FORMAT (///' *** ERROR IN QS2GRD -- IER =', i2, ' ***')
5600 FORMAT (///' *** ERROR -- INTERPOLATED VALUES ',  &
             'Q1 (QS2VAL) AND Q2 (QS2GRD) DIFFER --'//  &
             t7, 'Q1 = ', e21.14, '     Q2 = ', e21.14)

CONTAINS

! QUADRATIC TEST FUNCTION AND PARTIAL DERIVATIVES

! fq(xx,yy) = ((xx + 2.*yy)/3.) ** 2
! fx(xx,yy) = 2. * (xx + 2.*yy) / 9.
! fy(xx,yy) = 4. * (xx + 2.*yy) / 9.


FUNCTION fq(xx, yy) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy
REAL              :: fn_val

fn_val = ((xx + 2.*yy)/3.) ** 2
RETURN
END FUNCTION fq



FUNCTION fx(xx, yy) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy
REAL              :: fn_val

fn_val = 2. * (xx + 2.*yy) / 9.
RETURN
END FUNCTION fx



FUNCTION fy(xx, yy) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy
REAL              :: fn_val

fn_val = 4. * (xx + 2.*yy) / 9.
RETURN
END FUNCTION fy

END PROGRAM qs2test
