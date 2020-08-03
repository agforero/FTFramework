MODULE qshep3d

!      ALGORITHM 661, COLLECTED ALGORITHMS FROM ACM.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-01-18  Time: 17:33:38
 
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 14, NO. 2, P.151.

IMPLICIT NONE


CONTAINS


SUBROUTINE qshep3(n,x,y,z,f,nq,nw,nr,lcell,lnext,xyzmin,xyzdel,  &
                  rmax,rsq,a,ier)
INTEGER, INTENT(IN)   :: n, nq, nw, nr
REAL, INTENT(IN)      :: x(:), y(:), z(:), f(:)
REAL, INTENT(OUT)     :: xyzmin(3), xyzdel(3), rmax, rsq(:), a(:,:)
INTEGER, INTENT(OUT)  :: lcell(:,:,:), lnext(:), ier

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    01/08/90

!   THIS SUBROUTINE COMPUTES A SET OF PARAMETERS A AND RSQ DEFINING A SMOOTH
! (ONCE CONTINUOUSLY DIFFERENTIABLE) TRI-VARIATE FUNCTION Q(X,Y,Z) WHICH
! INTERPOLATES DATA VALUES F AT SCATTERED NODES (X,Y,Z).  THE INTERPOLANT Q
! MAY BE EVALUATED AT AN ARBITRARY POINT BY FUNCTION QS3VAL, AND ITS FIRST
! DERIVATIVES ARE COMPUTED BY SUBROUTINE QS3GRD.

!   THE INTERPOLATION SCHEME IS A MODIFIED QUADRATIC SHEPARD METHOD --

! Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))

! FOR TRIVARIATE FUNCTIONS W(K) AND Q(K).  THE NODAL FUNCTIONS ARE GIVEN BY

! Q(K)(X,Y,Z) = A(1,K)*DX**2 + A(2,K)*DX*DY + A(3,K)*DY**2 +
!                A(4,K)*DX*DZ + A(5,K)*DY*DZ + A(6,K)*DZ**2 +
!                A(7,K)*DX + A(8,K)*DY + A(9,K)*DZ + F(K)

! WHERE DX = (X-X(K)), DY = (Y-Y(K)), AND DZ = (Z-Z(K)).
! THUS, Q(K) IS A QUADRATIC FUNCTION WHICH INTERPOLATES THE DATA VALUE AT
! NODE K.  ITS COEFFICIENTS A(,K) ARE OBTAINED BY A WEIGHTED LEAST SQUARES
! FIT TO THE CLOSEST NQ DATA POINTS WITH WEIGHTS SIMILAR TO W(K).
! NOTE THAT THE RADIUS OF INFLUENCE FOR THE LEAST SQUARES FIT IS FIXED FOR
! EACH K, BUT VARIES WITH K.

!   THE WEIGHTS ARE TAKEN TO BE

! W(K)(X,Y,Z) = ( (R(K)-D(K))+ / R(K)*D(K) )**2

! WHERE (R(K)-D(K))+ = 0 IF R(K) .LE. D(K), AND D(K)(X,Y,Z) IS THE EUCLIDEAN
! DISTANCE BETWEEN (X,Y,Z) AND NODE K.  THE RADIUS OF INFLUENCE R(K) VARIES
! WITH K AND IS CHOSEN SO THAT NW NODES ARE WITHIN THE RADIUS.
! NOTE THAT W(K) IS NOT DEFINED AT NODE (X(K),Y(K),Z(K)), BUT Q(X,Y,Z) HAS
! LIMIT F(K) AS (X,Y,Z) APPROACHES (X(K),Y(K),Z(K)).

! ON INPUT --

!        N = NUMBER OF NODES AND ASSOCIATED DATA VALUES.
!            N .GE. 10.

!        X,Y,Z = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!                COORDINATES OF THE NODES.

!        F = ARRAY OF LENGTH N CONTAINING THE DATA VALUES
!            IN ONE-TO-ONE CORRESPONDENCE WITH THE NODES.

!        NQ = NUMBER OF DATA POINTS TO BE USED IN THE LEAST SQUARES FIT FOR
!             COEFFICIENTS DEFINING THE NODAL FUNCTIONS Q(K).
!             A RECOMMENDED VALUE IS NQ = 17.
!             9 .LE. NQ .LE. MIN(40,N-1).

!        NW = NUMBER OF NODES WITHIN (AND DEFINING) THE RADII OF INFLUENCE
!             R(K) WHICH ENTER INTO THE WEIGHTS W(K).
!             FOR N SUFFICIENTLY LARGE, A RECOMMENDED VALUE IS NW = 32.
!             1 .LE. NW .LE. MIN(40,N-1).

!        NR = NUMBER OF ROWS, COLUMNS, AND PLANES IN THE CELL GRID DEFINED
!             IN SUBROUTINE STORE3.  A BOX CONTAINING THE NODES IS PARTITIONED
!             INTO CELLS IN ORDER TO INCREASE SEARCH EFFICIENCY.
!             NR = (N/3)**(1/3) IS RECOMMENDED.  NR .GE. 1.

! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

!        LCELL = ARRAY OF LENGTH .GE. NR**3.

!        LNEXT = ARRAY OF LENGTH .GE. N.

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH .GE. 3.

!        RSQ = ARRAY OF LENGTH .GE. N.

!        A = ARRAY OF LENGTH .GE. 9N.

! ON OUTPUT --

!        LCELL = NR BY NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
!                REFER TO STORE3.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES.
!                REFER TO STORE3.

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH 3 CONTAINING MINIMUM NODAL
!                        COORDINATES AND CELL DIMENSIONS, RESPECTIVELY.
!                        REFER TO STORE3.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!               MAXIMUM RADIUS R(K).

!        RSQ = ARRAY CONTAINING THE SQUARES OF THE RADII R(K)
!              WHICH ENTER INTO THE WEIGHTS W(K).

!        A = 9 BY N ARRAY CONTAINING THE COEFFICIENTS FOR
!            QUADRATIC NODAL FUNCTION Q(K) IN COLUMN K.

!   NOTE THAT THE ABOVE OUTPUT PARAMETERS ARE NOT DEFINED UNLESS IER = 0.

!        IER = ERROR INDICATOR --
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N, NQ, NW, OR NR IS OUT OF RANGE.
!              IER = 2 IF DUPLICATE NODES WERE ENCOUNTERED.
!              IER = 3 IF ALL NODES ARE COPLANAR.

! MODULES REQUIRED BY QSHEP3 -- GETNP3, GIVENS, ROTATE,
!                                  SETUP3, STORE3

! INTRINSIC FUNCTIONS CALLED BY QSHEP3 -- ABS, MIN, FLOAT, MAX, MIN, SQRT

!***********************************************************

INTEGER :: i, ib, ierr, ip1, irm1, irow, j, jp1, k, lmax, lnp, neq,  &
           nn, nnq, nnr, nnw, np, npts(40), nqwmax
REAL :: av, avsq, b(10,10), c, dmin, fk, rq, rs, rsmx, rsold,  &
        rws, s, sum, t, xk, yk, zk, xyzdl(3), xyzmn(3)

REAL, PARAMETER  :: dtol = 0.01, rtol = 1.e-5, sf = 1.0

! LOCAL PARAMETERS --

! AV =         ROOT-MEAN-SQUARE DISTANCE BETWEEN K AND THE NODES IN THE
!                 LEAST SQUARES SYSTEM (UNLESS ADDITIONAL NODES ARE INTRODUCED
!                 FOR STABILITY).
!                 THE FIRST 6 COLUMNS OF THE MATRIX ARE SCALED BY 1/AVSQ,
!                 THE LAST 3 BY 1/AV
! AVSQ =       AV*AV
! B =          TRANSPOSE OF THE AUGMENTED REGRESSION MATRIX
! C =          FIRST COMPONENT OF THE PLANE ROTATION USED TO ZERO THE LOWER
!                 TRIANGLE OF B**T -- COMPUTED BY SUBROUTINE GIVENS
! DMIN =       MINIMUM OF THE MAGNITUDES OF THE DIAGONAL ELEMENTS OF THE
!                 REGRESSION MATRIX AFTER ZEROS ARE INTRODUCED BELOW THE
!                 DIAGONAL
! DTOL =       TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!                 SYSTEM.  THE SYSTEM IS ACCEPTED WHEN DMIN .GE. DTOL
! FK =         DATA VALUE AT NODE K -- F(K)
! I =          INDEX FOR A, B, NPTS, XYZMIN, XYZMN, XYZDEL,
!                 AND XYZDL
! IB =         DO-LOOP INDEX FOR BACK SOLVE
! IERR =       ERROR FLAG FOR THE CALL TO STORE3
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       ROW INDEX FOR B
! J =          INDEX FOR A AND B
! JP1 =        J+1
! K =          NODAL FUNCTION INDEX AND COLUMN INDEX FOR A
! LMAX =       MAXIMUM NUMBER OF NPTS ELEMENTS (MUST BE CONSISTENT WITH
!                 THE DIMENSION STATEMENT ABOVE)
! LNP =        CURRENT LENGTH OF NPTS
! NEQ =        NUMBER OF EQUATIONS IN THE LEAST SQUARES FIT
! NN,NNQ,NNR = LOCAL COPIES OF N, NQ, AND NR
! NNW =        LOCAL COPY OF NW
! NP =         NPTS ELEMENT
! NPTS =       ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!                 NODES TO BE USED IN THE LEAST SQUARES FIT
!                 OR TO COMPUTE RSQ.  THE NODES ARE ORDERED
!                 BY DISTANCE FROM K AND THE LAST ELEMENT
!                 (USUALLY INDEXED BY LNP) IS USED ONLY TO
!                 DETERMINE RQ, OR RSQ(K) IF NW .GT. NQ
! NQWMAX =     MAX(NQ,NW)
! RQ =         RADIUS OF INFLUENCE WHICH ENTERS INTO THE
!                 WEIGHTS FOR Q(K) (SEE SUBROUTINE SETUP3)
! RS =         SQUARED DISTANCE BETWEEN K AND NPTS(LNP) --
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
! S =          SECOND COMPONENT OF THE PLANE GIVENS ROTATION
! SF =         MARQUARDT STABILIZATION FACTOR USED TO DAMP
!                 OUT THE FIRST 6 SOLUTION COMPONENTS (SECOND
!                 PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!                 IS ILL-CONDITIONED.  AS SF INCREASES, THE
!                 FITTING FUNCTION APPROACHES A LINEAR
! SUM =        SUM OF SQUARED EUCLIDEAN DISTANCES BETWEEN
!                 NODE K AND THE NODES USED IN THE LEAST
!                 SQUARES FIT (UNLESS ADDITIONAL NODES ARE
!                 ADDED FOR STABILITY)
! T =               TEMPORARY VARIABLE FOR ACCUMULATING A SCALAR
!                 PRODUCT IN THE BACK SOLVE
! XK,YK,ZK =   COORDINATES OF NODE K -- X(K), Y(K), Z(K)
! XYZDL =      LOCAL VARIABLES FOR XYZDEL
! XYZMN =      LOCAL VARIABLES FOR XYZMIN

nn = n
nnq = nq
nnw = nw
nnr = nr
nqwmax = MAX(nnq,nnw)
lmax = MIN(40,nn-1)
IF (9 <= nnq .AND. 1 <= nnw .AND. nqwmax <= lmax .AND. nnr >= 1) THEN
  
! CREATE THE CELL DATA STRUCTURE, AND INITIALIZE RSMX.
  
  CALL store3(nn,x,y,z,nnr,lcell,lnext,xyzmn,xyzdl,ierr)
  IF (ierr /= 0) GO TO 160
  rsmx = 0.
  
! OUTER LOOP ON NODE K
  
  DO  k = 1, nn
    xk = x(k)
    yk = y(k)
    zk = z(k)
    fk = f(k)
    
! MARK NODE K TO EXCLUDE IT FROM THE SEARCH FOR NEAREST NEIGHBORS.
    
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
      CALL getnp3(xk,yk,zk,x,y,z,nnr,lcell,lnext,xyzmn,xyzdl,np,rs )
      IF (rs == 0.) GO TO 150
      npts(lnp) = np
      IF ((rs-rsold)/rs < rtol) GO TO 10
      IF (rws == 0. .AND. lnp > nnw) rws = rs
      IF (rq == 0. .AND. lnp > nnq) THEN
        
!   RQ = 0 (NOT YET COMPUTED) AND LNP .GT. NQ.        RQ = SQRT(RS) IS
!     SUFFICIENTLY LARGE TO (STRICTLY) INCLUDE NQ NODES.  THE LEAST SQUARES
!     FIT WILL INCLUDE NEQ = LNP-1 EQUATIONS FOR 9 .LE. NQ .LE. NEQ .LT. LMAX
!     .LE. N-1.
        
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
    
! SET UP THE AUGMENTED REGRESSION MATRIX (TRANSPOSED) AS THE
!   COLUMNS OF B, AND ZERO OUT THE LOWER TRIANGLE (UPPER
!   TRIANGLE OF B) WITH GIVENS ROTATIONS -- QR DECOMPOSITION
!   WITH ORTHOGONAL MATRIX Q NOT STORED.
    
    i = 0
    30 i = i + 1
    np = npts(i)
    irow = MIN(i,10)
    CALL setup3(xk,yk,zk,fk,x(np),y(np),z(np),f(np),av,avsq,rq, b(1,irow))
    IF (i == 1) GO TO 30
    irm1 = irow - 1
    DO  j = 1, irm1
      jp1 = j + 1
      CALL givens(b(j,j),b(j,irow),c,s)
      CALL rotate(10-j,c,s,b(jp1:,j),b(jp1:,irow))
    END DO
    IF (i < neq) GO TO 30
    
! TEST THE SYSTEM FOR ILL-CONDITIONING.
    
    dmin = MIN(ABS(b(1,1)), ABS(b(2,2)), ABS(b(3,3)), ABS(b(4,4)),  &
               ABS(b(5,5)), ABS(b(6,6)), ABS(b(7,7)), ABS(b(8,8)),  &
               ABS(b(9,9)))
    IF (dmin*rq < dtol) THEN
      IF (neq /= lmax) THEN
        
! INCREASE RQ AND ADD ANOTHER EQUATION TO THE SYSTEM TO IMPROVE THE
!   CONDITIONING.  THE NUMBER OF NPTS ELEMENTS IS ALSO INCREASED IF NECESSARY.
        
        50 rsold = rs
        neq = neq + 1
        IF (neq /= lmax) THEN
          IF (neq /= lnp) THEN
            
!   NEQ .LT. LNP
            
            np = npts(neq+1)
            rs = (x(np)-xk) ** 2 + (y(np)-yk) ** 2 + (z(np)-zk) ** 2
            IF ((rs-rsold)/rs < rtol) GO TO 50
            rq = SQRT(rs)
            GO TO 30
          END IF
          
!   ADD AN ELEMENT TO NPTS.
          
          lnp = lnp + 1
          CALL getnp3(xk,yk,zk,x,y,z,nnr,lcell,lnext,xyzmn,xyzdl, np,rs)
          IF (np == 0) GO TO 150
          npts(lnp) = np
          IF ((rs-rsold)/rs < rtol) GO TO 50
          rq = SQRT(rs)
          GO TO 30
        END IF
        
        rq = SQRT(1.1*rs)
        GO TO 30
      END IF
      
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS -- ADD
!   MULTIPLES OF THE FIRST SIX UNIT VECTORS TO THE FIRST
!   SIX EQUATIONS.
      
      DO  i = 1, 6
        b(i,10) = sf
        ip1 = i + 1
        DO  j = ip1, 10
          b(j,10) = 0.
        END DO
        DO  j = i, 9
          jp1 = j + 1
          CALL givens(b(j,j),b(j,10),c,s)
          CALL rotate(10-j,c,s,b(jp1:,j),b(jp1:,10))
        END DO
      END DO
      
! TEST THE STABILIZED SYSTEM FOR ILL-CONDITIONING.
      
      dmin = MIN(ABS(b(1,1)), ABS(b(2,2)), ABS(b(3,3)), ABS(b(4,4)), &
                 ABS(b(5,5)), ABS(b(6,6)), ABS(b(7,7)), ABS(b(8,8)), &
                 ABS(b(9,9)))
      IF (dmin*rq < dtol) GO TO 160
    END IF
    
! SOLVE THE 9 BY 9 TRIANGULAR SYSTEM FOR THE COEFFICIENTS
    
    DO  ib = 1, 9
      i = 10 - ib
      t = 0.
      IF (i /= 9) THEN
        ip1 = i + 1
        DO  j = ip1, 9
          t = t + b(j,i) * a(j,k)
        END DO
      END IF
      a(i,k) = (b(10,i)-t) / b(i,i)
    END DO
    
! SCALE THE COEFFICIENTS TO ADJUST FOR THE COLUMN SCALING.
    
    DO  i = 1, 6
      a(i,k) = a(i,k) / avsq
    END DO
    a(7,k) = a(7,k) / av
    a(8,k) = a(8,k) / av
    a(9,k) = a(9,k) / av
    
! UNMARK K AND THE ELEMENTS OF NPTS.
    
    lnext(k) = -lnext(k)
    DO  i = 1, lnp
      np = npts(i)
      lnext(np) = -lnext(np)
    END DO
  END DO
  
! NO ERRORS ENCOUNTERED.
  
  DO  i = 1, 3
    xyzmin(i) = xyzmn(i)
    xyzdel(i) = xyzdl(i)
  END DO
  rmax = SQRT(rsmx)
  ier = 0
  RETURN
END IF

! N, NQ, NW, OR NR IS OUT OF RANGE.

ier = 1
RETURN

! DUPLICATE NODES WERE ENCOUNTERED BY GETNP3.

150 ier = 2
RETURN

! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES.

160 DO  i = 1, 3
  xyzmin(i) = xyzmn(i)
  xyzdel(i) = xyzdl(i)
END DO
ier = 3
RETURN
END SUBROUTINE qshep3



FUNCTION qs3val(px,py,pz,n,x,y,z,f,nr,lcell,lnext,xyzmin,xyzdel,rmax,rsq,a) &
         RESULT(fn_val)
INTEGER, INTENT(IN)  :: n, nr, lcell(:,:,:), lnext(:)
REAL, INTENT(IN)     :: px, py, pz, x(:), y(:), z(:), f(:), xyzmin(3),  &
                        xyzdel(3), rmax, rsq(:), a(9,n)
REAL                 :: fn_val

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    10/28/87

!   THIS FUNCTION RETURNS THE VALUE Q(PX,PY,PZ) WHERE Q IS
! THE WEIGHTED SUM OF QUADRATIC NODAL FUNCTIONS DEFINED IN
! SUBROUTINE QSHEP3.  QS3GRD MAY BE CALLED TO COMPUTE A
! GRADIENT OF Q ALONG WITH THE VALUE, AND/OR TO TEST FOR
! ERRORS.

! ON INPUT --

!        PX,PY,PZ = CARTESIAN COORDINATES OF THE POINT P AT
!                   WHICH Q IS TO BE EVALUATED.

!        N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!            N .GE. 10.

!        X,Y,Z,F = ARRAYS OF LENGTH N CONTAINING THE NODES
!                  AND DATA VALUES INTERPOLATED BY Q.

!        NR = NUMBER OF ROWS, COLUMNS, AND PLANES IN THE CELL
!             GRID.  REFER TO STORE3.  NR .GE. 1.

!        LCELL = NR BY NR BY NR ARRAY OF NODAL INDICES ASSO-
!                CIATED WITH CELLS.  REFER TO STORE3.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!                CES.  REFER TO STORE3.

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH 3 CONTAINING MINI-
!                        MUM NODAL COORDINATES AND CELL DIM-
!                        ENSIONS, RESPECTIVELY.        XYZDEL
!                        ELEMENTS MUST BE POSITIVE.  REFER TO
!                        STORE3.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!               MAXIMUM RADIUS.

!        RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!              WHICH ENTER INTO THE WEIGHTS DEFINING Q.

!        A = 9 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!            NODAL FUNCTIONS DEFINING Q.

!   INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.        THE
! PARAMETERS OTHER THAN PX, PY AND PZ SHOULD BE INPUT UNAL-
! TERED FROM THEIR VALUES ON OUTPUT FROM QSHEP3.  THIS FUNC-
! TION SHOULD NOT BE CALLED IF A NONZERO ERROR FLAG WAS
! RETURNED BY QSHEP3.

! ON OUTPUT --

!        QS3VAL = FUNCTION VALUE Q(PX,PY,PZ) UNLESS N, NR,
!                 XYZDEL, OR RMAX IS INVALID, IN WHICH CASE
!                 NO VALUE IS RETURNED.

! MODULES REQUIRED BY QS3VAL -- NONE

! INTRINSIC FUNCTIONS CALLED BY QS3VAL -- IFIX, SQRT

!***********************************************************

! Local variables
REAL     :: delx, dely, delz, ds, dx, dxsq, dy, dysq, dz, dzsq,  &
            rd, rds, rs, sw, swq, w, xmin, xp, ymin, yp, zmin, zp
INTEGER  :: i, imax, imin, j, jmax, jmin, k, kmax, kmin, l, lp

xp = px
yp = py
zp = pz
xmin = xyzmin(1)
ymin = xyzmin(2)
zmin = xyzmin(3)
dx = xyzdel(1)
dy = xyzdel(2)
dz = xyzdel(3)
IF (n < 10 .OR. nr < 1 .OR. dx <= 0. .OR. dy <= 0. .OR. dz <= 0. .OR. rmax  &
     < 0.) RETURN

! SET IMIN, IMAX, JMIN, JMAX, KMIN, AND KMAX TO CELL INDICES DEFINING THE
!   RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE P.
!   THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTERSECTED BY (OR CONTAINED
!   IN) A SPHERE OF RADIUS RMAX CENTERED AT P.

imin = (xp-xmin-rmax)/dx + 1
imax = (xp-xmin+rmax)/dx + 1
IF (imin < 1) imin = 1
IF (imax > nr) imax = nr
jmin = (yp-ymin-rmax)/dy + 1
jmax = (yp-ymin+rmax)/dy + 1
IF (jmin < 1) jmin = 1
IF (jmax > nr) jmax = nr
kmin = (zp-zmin-rmax)/dz + 1
kmax = (zp-zmin+rmax)/dz + 1
IF (kmin < 1) kmin = 1
IF (kmax > nr) kmax = nr

! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE SPHERE OF RADIUS RMAX.

IF (imin <= imax .AND. jmin <= jmax .AND. kmin <= kmax) THEN
  
! ACCUMULATE WEIGHT VALUES IN SW AND WEIGHTED NODAL FUNCTION
!   VALUES IN SWQ.  THE WEIGHTS ARE W(L) = ((R-D)+/(R*D))**2
!   FOR R**2 = RSQ(L) AND D = DISTANCE BETWEEN P AND NODE L.
  
  sw = 0.
  swq = 0.
  
! OUTER LOOP ON CELLS (I,J,K).
  
  DO  k = kmin, kmax
    DO  j = jmin, jmax
      DO  i = imin, imax
        l = lcell(i,j,k)
        IF (l /= 0) THEN
          
! INNER LOOP ON NODES L.
          
          10 delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)
          dxsq = delx * delx
          dysq = dely * dely
          dzsq = delz * delz
          ds = dxsq + dysq + dzsq
          rs = rsq(l)
          IF (ds < rs) THEN
            IF (ds == 0.) GO TO 50
            rds = rs * ds
            rd = SQRT(rds)
            w = (rs+ds-rd-rd) / rds
            sw = sw + w
            swq = swq + w * (a(1,l)*dxsq+a(2,l)*delx*dely+a(3,l)*  &
                dysq+a(4,l)*delx*delz+a(5,l)*dely*delz+a(6,l)*dzsq  &
                +a(7,l)*delx+a(8,l)*dely+a(9,l)*delz+f(l))
          END IF
          
! BOTTOM OF LOOP ON NODES IN CELL (I,J,K).
          
          lp = l
          l = lnext(lp)
          IF (l /= lp) GO TO 10
        END IF
      END DO
    END DO
  END DO
  
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(L) FOR ANY NODE L.
  
  IF (sw == 0.) GO TO 60
  fn_val = swq / sw
  RETURN
  
! (PX,PY,PZ) = (X(L),Y(L),Z(L))
  
  50 fn_val = f(l)
  RETURN
END IF

! ALL WEIGHTS ARE 0 AT P.

60 fn_val = 0.
RETURN
END FUNCTION qs3val



SUBROUTINE qs3grd(px,py,pz,n,x,y,z,f,nr,lcell,lnext,xyzmin,xyzdel,  &
                  rmax,rsq,a,q,qx,qy,qz,ier)
INTEGER, INTENT(IN)   :: n, nr, lcell(:,:,:), lnext(:)
REAL, INTENT(IN)      :: px, py, pz, x(:), y(:), z(:), f(:), xyzmin(3), &
                         xyzdel(3), rmax, rsq(:), a(9,n)
INTEGER, INTENT(OUT)  :: ier
REAL, INTENT(OUT)     :: q, qx, qy, qz

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767
!                                                    10/28/87

!   THIS SUBROUTINE COMPUTES THE VALUE AND GRADIENT AT (PX,PY,PZ) OF THE
! INTERPOLATORY FUNCTION Q DEFINED IN SUBROUTINE QSHEP3.  Q(X,Y,Z) IS A
! WEIGHTED SUM OF QUADRATIC NODAL FUNCTIONS.

! ON INPUT --

!        PX,PY,PZ = CARTESIAN COORDINATES OF THE POINT P AT WHICH Q AND
!                   ITS PARTIALS ARE TO BE EVALUATED.

!        N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!            N .GE. 10.

!        X,Y,Z,F = ARRAYS OF LENGTH N CONTAINING THE NODES
!                  AND DATA VALUES INTERPOLATED BY Q.

!        NR = NUMBER OF ROWS, COLUMNS AND PLANES IN THE CELL GRID.
!             REFER TO STORE3.  NR .GE. 1.

!        LCELL = NR BY NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
!                REFER TO STORE3.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES.
!                REFER TO STORE3.

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH 3 CONTAINING MINIMUM NODAL
!                        COORDINATES AND CELL DIMENSIONS, RESPECTIVELY.
!                        XYZDEL ELEMENTS MUST BE POSITIVE.  REFER TO STORE3.

!        RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ -- MAXIMUM RADIUS.

!        RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!              WHICH ENTER INTO THE WEIGHTS DEFINING Q.

!        A = 9 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!            NODAL FUNCTIONS DEFINING Q.

!   INPUT PARAMETERS ARE NOT ALTERED BY THIS SUBROUTINE.
! THE PARAMETERS OTHER THAN PX, PY, AND PZ SHOULD BE INPUT UNALTERED FROM
! THEIR VALUES ON OUTPUT FROM QSHEP3.  THIS SUBROUTINE SHOULD NOT BE CALLED
! IF A NONZERO ERROR FLAG WAS RETURNED BY QSHEP3.

! ON OUTPUT --

!        Q = VALUE OF Q AT (PX,PY,PZ) UNLESS IER .EQ. 1, IN
!            WHICH CASE NO VALUES ARE RETURNED.

!        QX,QY,QZ = FIRST PARTIAL DERIVATIVES OF Q AT
!                   (PX,PY,PZ) UNLESS IER .EQ. 1.

!        IER = ERROR INDICATOR
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N, NR, XYZDEL, OR RMAX IS INVALID.
!              IER = 2 IF NO ERRORS WERE ENCOUNTERED BUT
!                      (PX,PY,PZ) IS NOT WITHIN THE RADIUS
!                      R(K) FOR ANY NODE K (AND THUS Q = QX =
!                      QY = QZ = 0).

! MODULES REQUIRED BY QS3GRD -- NONE

! INTRINSIC FUNCTIONS CALLED BY QS3GRD -- IFIX, SQRT

!***********************************************************

! Local variables
REAL     :: delx, dely, delz, ds, dx, dxsq, dy, dysq, dz, dzsq,  &
            ql, qlx, qly, qlz, rd, rds, rs, sw, swq, swqx, swqy, swqz, &
            sws, swx, swy, swz, t, w, wx, wy, wz, xmin, xp, ymin, yp,  &
            zmin, zp
INTEGER  :: i, imax, imin, j, jmax, jmin, k, kmax, kmin, l, lp

xp = px
yp = py
zp = pz
xmin = xyzmin(1)
ymin = xyzmin(2)
zmin = xyzmin(3)
dx = xyzdel(1)
dy = xyzdel(2)
dz = xyzdel(3)
IF (n >= 10 .AND. nr >= 1 .AND. dx > 0. .AND. dy > 0. .AND. dz > 0.  &
       .AND. rmax >= 0.) THEN
  
! SET IMIN, IMAX, JMIN, JMAX, KMIN, AND KMAX TO CELL INDICES DEFINING THE
!   RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE P.
!   THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTERSECTED BY (OR CONTAINED
!   IN) A SPHERE OF RADIUS RMAX CENTERED AT P.
  
  imin = (xp-xmin-rmax)/dx + 1
  imax = (xp-xmin+rmax)/dx + 1
  IF (imin < 1) imin = 1
  IF (imax > nr) imax = nr
  jmin = (yp-ymin-rmax)/dy + 1
  jmax = (yp-ymin+rmax)/dy + 1
  IF (jmin < 1) jmin = 1
  IF (jmax > nr) jmax = nr
  kmin = (zp-zmin-rmax)/dz + 1
  kmax = (zp-zmin+rmax)/dz + 1
  IF (kmin < 1) kmin = 1
  IF (kmax > nr) kmax = nr
  
! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE SPHERE OF RADIUS RMAX.
  
  IF (imin > imax .OR. jmin > jmax .OR. kmin > kmax) GO TO 60
  
! Q = SWQ/SW = SUM(W(L)*Q(L))/SUM(W(L)) WHERE THE SUM IS
!   FROM L = 1 TO N, Q(L) IS THE QUADRATIC NODAL FUNCTION,
!   AND W(L) = ((R-D)+/(R*D))**2 FOR RADIUS R(L) AND DIST-
!   ANCE D(L).        THUS
  
!         QX = (SWQX*SW - SWQ*SWX)/SW**2
!         QY = (SWQY*SW - SWQ*SWY)/SW**2
!         QZ = (SWQZ*SW - SWQ*SWZ)/SW**2
  
!   WHERE SWQX AND SWX ARE PARTIAL DERIVATIVES WITH RESPECT
!   TO X OF SWQ AND SW, RESPECTIVELY.  SWQY, SWY, SWQZ, AND
!   SWZ ARE DEFINED SIMILARLY.
  
  sw = 0.
  swx = 0.
  swy = 0.
  swz = 0.
  swq = 0.
  swqx = 0.
  swqy = 0.
  swqz = 0.
  
! OUTER LOOP ON CELLS (I,J,K).
  
  DO  k = kmin, kmax
    DO  j = jmin, jmax
      DO  i = imin, imax
        l = lcell(i,j,k)
        IF (l /= 0) THEN
          
! INNER LOOP ON NODES L.
          
          10 delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)
          dxsq = delx * delx
          dysq = dely * dely
          dzsq = delz * delz
          ds = dxsq + dysq + dzsq
          rs = rsq(l)
          IF (ds < rs) THEN
            IF (ds == 0.) GO TO 50
            rds = rs * ds
            rd = SQRT(rds)
            w = (rs+ds-rd-rd) / rds
            t = 2. * (rd-rs) / (ds*rds)
            wx = delx * t
            wy = dely * t
            wz = delz * t
            qlx = 2. * a(1,l) * delx + a(2,l) * dely + a(4,l) * delz
            qly = a(2,l) * delx + 2. * a(3,l) * dely + a(5,l) * delz
            qlz = a(4,l) * delx + a(5,l) * dely + 2. * a(6,l) * delz
            ql = (qlx*delx+qly*dely+qlz*delz) / 2. + a(7,l) * delx  &
                + a(8,l) * dely + a(9,l) * delz + f(l)
            qlx = qlx + a(7,l)
            qly = qly + a(8,l)
            qlz = qlz + a(9,l)
            sw = sw + w
            swx = swx + wx
            swy = swy + wy
            swz = swz + wz
            swq = swq + w * ql
            swqx = swqx + wx * ql + w * qlx
            swqy = swqy + wy * ql + w * qly
            swqz = swqz + wz * ql + w * qlz
          END IF
          
! BOTTOM OF LOOP ON NODES IN CELL (I,J,K).
          
          lp = l
          l = lnext(lp)
          IF (l /= lp) GO TO 10
        END IF
      END DO
    END DO
  END DO
  
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(L) FOR ANY NODE L.
  
  IF (sw == 0.) GO TO 60
  q = swq / sw
  sws = sw * sw
  qx = (swqx*sw-swq*swx) / sws
  qy = (swqy*sw-swq*swy) / sws
  qz = (swqz*sw-swq*swz) / sws
  ier = 0
  RETURN
  
! (PX,PY,PZ) = (X(L),Y(L),Z(L))
  
  50 q = f(l)
  qx = a(7,l)
  qy = a(8,l)
  qz = a(9,l)
  ier = 0
  RETURN
END IF

! INVALID INPUT PARAMETER.

ier = 1
RETURN

! NO CELLS CONTAIN A POINT WITHIN RMAX OF P, OR
!   SW = 0 AND THUS DS .GE. RSQ(L) FOR ALL L.

60 q = 0.
qx = 0.
qy = 0.
qz = 0.
ier = 2
RETURN
END SUBROUTINE qs3grd



SUBROUTINE getnp3(px,py,pz,x,y,z,nr,lcell,lnext,xyzmin,xyzdel,np,dsq)
INTEGER, INTENT(IN)      :: nr, lcell(:,:,:)
INTEGER, INTENT(IN OUT)  :: lnext(:)
REAL, INTENT(IN)         :: px, py, pz, x(:), y(:), z(:), xyzmin(3), xyzdel(3)
REAL, INTENT(OUT)        :: dsq
INTEGER, INTENT(OUT)     :: np

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   GIVEN A SET OF N NODES AND THE DATA STRUCTURE DEFINED IN SUBROUTINE
! STORE3, THIS SUBROUTINE USES THE CELL METHOD TO FIND THE CLOSEST UNMARKED
! NODE NP TO A SPECIFIED POINT P.  NP IS THEN MARKED BY SETTING LNEXT(NP) TO
! -LNEXT(NP).  (A NODE IS MARKED IF AND ONLY IF THE CORRESPONDING LNEXT
! ELEMENT IS NEGATIVE.  THE ABSOLUTE VALUES OF LNEXT ELEMENTS, HOWEVER, MUST
! BE PRESERVED.)   THUS, THE CLOSEST M NODES TO P MAY BE DETERMINED BY A
! SEQUENCE OF M CALLS TO THIS ROUTINE.  NOTE THAT IF THE NEAREST NEIGHBOR TO
! NODE K IS TO BE DETERMINED (PX = X(K), PY = Y(K), AND PZ = Z(K)), THEN
! K SHOULD BE MARKED BEFORE THE CALL TO THIS ROUTINE.

!   THE SEARCH IS BEGUN IN THE CELL CONTAINING (OR CLOSEST TO) P AND PROCEEDS
! OUTWARD IN BOX-SHAPED LAYERS UNTIL ALL CELLS WHICH CONTAIN POINTS WITHIN
! DISTANCE R OF P HAVE BEEN SEARCHED, WHERE R IS THE DISTANCE FROM P TO THE
! FIRST UNMARKED NODE ENCOUNTERED (INFINITE IF NO UNMARKED NODES ARE PRESENT).

! ON INPUT --

!        PX,PY,PZ = CARTESIAN COORDINATES OF THE POINT P WHOSE NEAREST
!                   UNMARKED NEIGHBOR IS TO BE FOUND.

!        X,Y,Z = ARRAYS OF LENGTH N, FOR N .GE. 2, CONTAINING
!                THE CARTESIAN COORDINATES OF THE NODES.

!        NR = NUMBER OF ROWS, COLUMNS, AND PLANES IN THE CELL GRID.
!             NR .GE. 1.

!        LCELL = NR BY NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.

!        LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDICES
!                (OR THEIR NEGATIVES).

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH 3 CONTAINING MINIMUM NODAL
!                        COORDINATES AND CELL DIMENSIONS, RESPECTIVELY.
!                        XYZDEL ELEMENTS MUST BE POSITIVE.

!   INPUT PARAMETERS OTHER THAN LNEXT ARE NOT ALTERED BY THIS ROUTINE.
! WITH THE EXCEPTION OF (PX,PY,PZ) AND THE SIGNS OF LNEXT ELEMENTS, THESE
! PARAMETERS SHOULD BE UNALTERED FROM THEIR VALUES ON OUTPUT FROM SUBROUTINE
! STORE3.

! ON OUTPUT --

!        NP = INDEX (FOR X, Y, AND Z) OF THE NEAREST UNMARKED
!             NODE TO P, OR 0 IF ALL NODES ARE MARKED OR NR
!             .LT. 1 OR AN ELEMENT OF XYZDEL IS NOT POSITIVE.
!             LNEXT(NP) .LT. 0 IF NP .NE. 0.

!        DSQ = SQUARED EUCLIDEAN DISTANCE BETWEEN P AND NODE
!              NP, OR 0 IF NP = 0.

! MODULES REQUIRED BY GETNP3 -- NONE

! INTRINSIC FUNCTIONS CALLED BY GETNP3 -- ABS, IFIX, SQRT

!***********************************************************

LOGICAL  :: first
REAL     :: delx, dely, delz, dx, dy, dz, r, rsmin, rsq, xp, yp, zp
INTEGER  :: i, i0, i1, i2, imax, imin, j, j0, j1, j2, jmax, jmin,  &
            k, k0, k1, k2, kmax, kmin, l, lmin, ln

xp = px
yp = py
zp = pz
dx = xyzdel(1)
dy = xyzdel(2)
dz = xyzdel(3)

! TEST FOR INVALID INPUT PARAMETERS.

IF (nr >= 1 .AND. dx > 0. .AND. dy > 0. .AND. dz > 0.) THEN
  
! INITIALIZE PARAMETERS --
  
!   FIRST = TRUE IFF THE FIRST UNMARKED NODE HAS YET TO BE
!            ENCOUNTERED,
!   IMIN,...,KMAX = CELL INDICES DEFINING THE RANGE OF THE
!                    SEARCH,
!   DELX,DELY,DELZ = PX-XYZMIN(1), PY-XYZMIN(2), AND
!                     PZ-XYZMIN(3),
!   I0,J0,K0 = CELL CONTAINING OR CLOSEST TO P,
!   I1,...,K2 = CELL INDICES OF THE LAYER WHOSE INTERSECTION
!                WITH THE RANGE DEFINED BY IMIN,...,KMAX IS
!                CURRENTLY BEING SEARCHED.
  
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  kmin = 1
  kmax = nr
  delx = xp - xyzmin(1)
  dely = yp - xyzmin(2)
  delz = zp - xyzmin(3)
  i0 = delx/dx + 1
  IF (i0 < 1) i0 = 1
  IF (i0 > nr) i0 = nr
  j0 = dely/dy + 1
  IF (j0 < 1) j0 = 1
  IF (j0 > nr) j0 = nr
  k0 = delz/dz + 1
  IF (k0 < 1) k0 = 1
  IF (k0 > nr) k0 = nr
  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
  k1 = k0
  k2 = k0
  
! OUTER LOOP ON LAYERS, INNER LOOP ON LAYER CELLS, EXCLUDING
!   THOSE OUTSIDE THE RANGE (IMIN,IMAX) X (JMIN,JMAX) X (KMIN,KMAX).
  10 i = 0
  loop50:  DO  k = k1, k2
    IF (k > kmax) GO TO 60
    IF (k >= kmin) THEN
      loop40:  DO  j = j1, j2
        IF (j > jmax) CYCLE loop50
        IF (j >= jmin) THEN
          DO  i = i1, i2
            IF (i > imax) CYCLE loop40
            IF (i >= imin) THEN
              IF (k == k1 .OR. k == k2 .OR. j == j1 .OR. j == j2 .OR. i  &
                     == i1 .OR. i == i2) THEN
                
! SEARCH CELL (I,J,K) FOR UNMARKED NODES L.
                
                l = lcell(i,j,k)
                IF (l /= 0) THEN
                  
!   LOOP ON NODES IN CELL (I,J,K).
                  
                  20 ln = lnext(l)
                  IF (ln >= 0) THEN
                    
!   NODE L IS NOT MARKED.
                    
                    rsq = (x(l)-xp) ** 2 + (y(l)-yp) ** 2 + ( z(l)-zp) ** 2
                    IF (first) THEN
                      
!   NODE L IS THE FIRST UNMARKED NEIGHBOR OF P ENCOUNTERED.
!     INITIALIZE LMIN TO THE CURRENT CANDIDATE FOR NP, AND
!     RSMIN TO THE SQUARED DISTANCE FROM P TO LMIN.  IMIN,
!     IMAX, JMIN, JMAX, KMIN, AND KMAX ARE UPDATED TO DEFINE
!     THE SMALLEST RECTANGLE CONTAINING A SPHERE OF RADIUS
!     R = SQRT(RSMIN) CENTERED AT P, AND CONTAINED IN (1,NR)
!     X (1,NR) X (1,NR) (EXCEPT THAT, IF P IS OUTSIDE THE
!     BOX DEFINED BY THE NODES, IT IS POSSIBLE THAT IMIN
!     .GT. NR OR IMAX .LT. 1, ETC.).  FIRST IS RESET TO
!     FALSE.
                      
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
                      kmin = (delz-r)/dz + 1
                      IF (kmin < 1) kmin = 1
                      kmax = (delz+r)/dz + 1
                      IF (kmax > nr) kmax = nr
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
                  
!   TEST FOR TERMINATION OF LOOP ON NODES IN CELL (I,J,K).
                  
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
    END IF
  END DO loop50
  
! TEST FOR TERMINATION OF LOOP ON CELL LAYERS.
  
  60 IF (i1 > imin .OR. i2 < imax .OR. j1 > jmin .OR. j2 < jmax .OR. k1  &
         > kmin .OR. k2 < kmax) THEN
    i1 = i1 - 1
    i2 = i2 + 1
    j1 = j1 - 1
    j2 = j2 + 1
    k1 = k1 - 1
    k2 = k2 + 1
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

! ERROR -- NR OR XYZDEL IS INVALID OR ALL NODES ARE MARKED.

np = 0
dsq = 0.
RETURN
END SUBROUTINE getnp3



SUBROUTINE givens(a,b,c,s)
REAL, INTENT(IN OUT)  :: a, b
REAL, INTENT(OUT)     :: c, s

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --
!     ( C  S)
! G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND
!     (-S  C)
! ENTRY OF THE 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS
! IS NORMALLY FOLLOWED BY A CALL TO ROTATE WHICH APPLIES
! THE TRANSFORMATION TO A 2 BY N MATRIX.  THIS ROUTINE WAS
! TAKEN FROM LINPACK.

! ON INPUT --

!        A,B = COMPONENTS OF THE 2-VECTOR TO BE ROTATED.

! ON OUTPUT --

!        A = VALUE OVERWRITTEN BY R = +/-SQRT(A*A + B*B)

!        B = VALUE OVERWRITTEN BY A VALUE Z WHICH ALLOWS C
!            AND S TO BE RECOVERED AS FOLLOWS --
!              C = SQRT(1-Z*Z), S=Z     IF ABS(Z) .LE. 1.
!              C = 1/Z, S = SQRT(1-C*C) IF ABS(Z) .GT. 1.

!        C = +/-(A/R)

!        S = +/-(B/R)

! MODULES REQUIRED BY GIVENS -- NONE

! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT

!***********************************************************

REAL :: aa, bb, r, u, v

! LOCAL PARAMETERS --

! AA,BB = LOCAL COPIES OF A AND B
! R =          C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   VARIABLES USED TO SCALE A AND B FOR COMPUTING R

aa = a
bb = b
IF (ABS(aa) > ABS(bb)) THEN
  
! ABS(A) .GT. ABS(B)
  
  u = aa + aa
  v = bb / u
  r = SQRT(.25+v*v) * u
  c = aa / r
  s = v * (c+c)
  
! NOTE THAT R HAS THE SIGN OF A, C .GT. 0, AND S HAS
!   SIGN(A)*SIGN(B).
  
  b = s
  a = r
  RETURN
END IF

! ABS(A) .LE. ABS(B)

IF (bb /= 0.) THEN
  u = bb + bb
  v = aa / u
  
! STORE R IN A.
  
  a = SQRT(.25+v*v) * u
  s = bb / a
  c = v * (s+s)
  
! NOTE THAT R HAS THE SIGN OF B, S .GT. 0, AND C HAS
!   SIGN(A)*SIGN(B).
  
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
! 2 BY N MATRIX (              ).
!                (Y(1) ... Y(N))

! ON INPUT --

!        N = NUMBER OF COLUMNS TO BE ROTATED.

!        C,S = ELEMENTS OF THE GIVENS ROTATION.    THESE MAY BE
!              DETERMINED BY SUBROUTINE GIVENS.

!        X,Y = ARRAYS OF LENGTH .GE. N CONTAINING THE VECTORS TO BE ROTATED.

! PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS ROUTINE.

! ON OUTPUT --

!        X,Y = ROTATED VECTORS.

! MODULES REQUIRED BY ROTATE -- NONE

!***********************************************************

INTEGER :: i
REAL :: xi, yi

! LOCAL PARAMETERS --

! I =  DO-LOOP INDEX
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



SUBROUTINE setup3(xk,yk,zk,fk,xi,yi,zi,fi,s1,s2,r,row)
REAL, INTENT(IN)      :: xk, yk, zk, fk, xi, yi, zi, fi, s1, s2, r
REAL, INTENT(IN OUT)  :: row(10)

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   THIS ROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED REGRESSION MATRIX FOR
! A WEIGHTED LEAST-SQUARES FIT OF A QUADRATIC FUNCTION Q(X,Y,Z) TO A SET OF
! DATA VALUES F, WHERE Q(XK,YK,ZK) = FK.  THE FIRST 6 COLUMNS (QUADRATIC
! TERMS) ARE SCALED BY 1/S2, AND COLUMNS 7, 8, AND 9 (LINEAR TERMS) ARE
! SCALED BY 1/S1.  THE WEIGHT IS (R-D)/(R*D) IF R .GT. D, AND 0 IF R .LE. D,
! WHERE D IS THE DISTANCE BETWEEN NODES I AND K.

! ON INPUT --

!        XK,YK,ZK,FK = COORDINATES AND DATA VALUE AT NODE K
!                      (INTERPOLATED BY Q).

!        XI,YI,ZI,FI = COORDINATES AND DATA VALUE AT NODE I.

!        S1,S2 = RECIPROCALS OF THE SCALE FACTORS.

!        R = RADIUS OF INFLUENCE ABOUT NODE K DEFINING THE WEIGHT.

!        ROW = ARRAY OF LENGTH 10.

! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

! ON OUTPUT --

!        ROW = ARRAY CONTAINING A ROW OF THE AUGMENTED REGRESSION MATRIX.

! MODULES REQUIRED BY SETUP3 -- NONE

! INTRINSIC FUNCTION CALLED BY SETUP3 -- SQRT

!***********************************************************

INTEGER :: i
REAL :: dx, dy, dz, dxsq, dysq, dzsq, d, w, w1, w2

! LOCAL PARAMETERS -

! I =         DO-LOOP INDEX
! DX =         XI - XK
! DY =         YI - YK
! DZ =         ZI - ZK
! DXSQ = DX*DX
! DYSQ = DY*DY
! DZSQ = DZ*DZ
! D =         DISTANCE BETWEEN NODES K AND I
! W =         WEIGHT ASSOCIATED WITH THE ROW
! W1 =         W/S1
! W2 =         W/S2

dx = xi - xk
dy = yi - yk
dz = zi - zk
dxsq = dx * dx
dysq = dy * dy
dzsq = dz * dz
d = SQRT(dxsq+dysq+dzsq)
IF (d > 0. .AND. d < r) THEN
  w = (r-d) / r / d
  w1 = w / s1
  w2 = w / s2
  row(1) = dxsq * w2
  row(2) = dx * dy * w2
  row(3) = dysq * w2
  row(4) = dx * dz * w2
  row(5) = dy * dz * w2
  row(6) = dzsq * w2
  row(7) = dx * w1
  row(8) = dy * w1
  row(9) = dz * w1
  row(10) = (fi-fk) * w
  RETURN
END IF

! NODES K AND I COINCIDE OR NODE I IS OUTSIDE OF THE RADIUS
!   OF INFLUENCE.  SET ROW TO THE ZERO VECTOR.

DO  i = 1, 10
  row(i) = 0.
END DO
RETURN
END SUBROUTINE setup3



SUBROUTINE store3(n,x,y,z,nr,lcell,lnext,xyzmin,xyzdel,ier)
INTEGER, INTENT(IN)   :: n, nr
REAL, INTENT(IN)      :: x(:), y(:), z(:)
INTEGER, INTENT(OUT)  :: lcell(:,:,:), lnext(:), ier
REAL, INTENT(OUT)     :: xyzmin(3), xyzdel(3)

!***********************************************************

!                                                ROBERT RENKA
!                                        UNIV. OF NORTH TEXAS
!                                              (817) 565-2767

!   GIVEN A SET OF N ARBITRARILY DISTRIBUTED NODES IN THREE-SPACE, THIS
! SUBROUTINE CREATES A DATA STRUCTURE FOR A CELL-BASED METHOD OF SOLVING
! CLOSEST-POINT PROBLEMS.   THE SMALLEST BOX CONTAINING THE NODES IS
! PARTITIONED INTO AN NR BY NR BY NR UNIFORM GRID OF CELLS, AND NODES ARE
! ASSOCIATED WITH CELLS.    IN PARTICULAR, THE DATA STRUCTURE STORES THE
! INDICES OF THE NODES CONTAINED IN EACH CELL.
! FOR A UNIFORM RANDOM DISTRIBUTION OF NODES, THE NEAREST NODE TO AN
! ARBITRARY POINT CAN BE DETERMINED IN CONSTANT EXPECTED TIME.

! ON INPUT --

!        N = NUMBER OF NODES.  N .GE. 2.

!        X,Y,Z = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!                COORDINATES OF THE NODES.

!        NR = NUMBER OF ROWS, COLUMNS, AND PLANES IN THE GRID.
!             THE CELL DENSITY (AVERAGE NUMBER OF NODES PER CELL) IS
!             D = N/(NR**3).    A RECOMMENDED VALUE, BASED ON EMPIRICAL
!             EVIDENCE, IS D = 3
!             -- NR = (N/3)**(1/3).  NR .GE. 1.

! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

!        LCELL = ARRAY OF LENGTH .GE. NR**3.

!        LNEXT = ARRAY OF LENGTH .GE. N.

!        XYZMIN,XYZDEL = ARRAYS OF LENGTH .GE. 3.

! ON OUTPUT --

!        LCELL = NR BY NR BY NR CELL ARRAY SUCH THAT
!                LCELL(I,J,K) CONTAINS THE INDEX (FOR X, Y,
!                AND Z) OF THE FIRST NODE (NODE WITH SMALLEST
!                INDEX) IN CELL (I,J,K), OR LCELL(I,J,K) = 0
!                IF NO NODES ARE CONTAINED IN THE CELL.        THE
!                CORNER OF CELL (I,J,K) WHICH IS FARTHEST
!                FROM THE BOX CORNER DEFINED BY XYZMIN HAS
!                COORDINATES (XMIN+I*DX,YMIN+J*DY,ZMIN+K*DZ),
!                WHERE (XMIN,YMIN,ZMIN) ARE THE ELEMENTS OF
!                XYZMIN.  LCELL IS NOT DEFINED IF IER .NE. 0.

!        LNEXT = ARRAY OF NEXT-NODE INDICES SUCH THAT
!                LNEXT(L) CONTAINS THE INDEX OF THE NEXT NODE
!                IN THE CELL WHICH CONTAINS NODE L, OR
!                LNEXT(L) = L IF L IS THE LAST NODE IN THE
!                CELL FOR L = 1,...,N.  (THE NODES CONTAINED
!                IN A CELL ARE ORDERED BY THEIR INDICES.)
!                IF, FOR EXAMPLE, CELL (I,J,K) CONTAINS NODES
!                2, 3, AND 5 (AND NO OTHERS), THEN
!                LCELL(I,J,K) = 2, LNEXT(2) = 3, LNEXT(3) =
!                5, AND LNEXT(5) = 5.  LNEXT IS NOT DEFINED
!                IF IER .NE. 0.

!        XYZMIN = ARRAY OF LENGTH 3 CONTAINING THE MINIMUM
!                 NODAL COORDINATES XMIN, YMIN, AND ZMIN (IN
!                 THAT ORDER) UNLESS IER = 1.  THE OPPOSITE
!                 CORNER OF THE BOX DEFINED BY THE NODES IS
!                 (XMIN+NR*DX,YMIN+NR*DY,ZMIN+NR*DZ).

!        XYZDEL = ARRAY OF LENGTH 3 CONTAINING THE DIMENSIONS
!                 OF THE CELLS UNLESS IER = 1.  XYZDEL(1) =
!                 (XMAX-XMIN)/NR, XYZDEL(2) = (YMAX-YMIN)/NR,
!                 AND XYZDEL(3) = (ZMAX-ZMIN)/NR, WHERE XMIN,
!                 XMAX, YMIN, YMAX, ZMIN, AND ZMAX ARE THE
!                 EXTREMA OF X, Y, AND Z.

!        IER = ERROR INDICATOR --
!              IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!              IER = 1 IF N .LT. 2 OR NR .LT. 1.
!              IER = 2 IF A COMPONENT OF XYZDEL IS NOT
!                      POSITIVE.

! MODULES REQUIRED BY STORE3 -- NONE

! INTRINSIC FUNCTIONS CALLED BY STORE3 -- FLOAT, IFIX

!***********************************************************

! Local variables
INTEGER  :: i, j, k, l, lb, ll, nn, nnr, np1
REAL     :: delx, dely, delz, xmn, xmx, ymn, ymx, zmn, zmx

nn = n
nnr = nr
IF (nn >= 2 .AND. nnr >= 1) THEN
  
! COMPUTE THE DIMENSIONS OF THE RECTANGLE CONTAINING THE
!   NODES.
  
  xmn = x(1)
  xmx = xmn
  ymn = y(1)
  ymx = ymn
  zmn = z(1)
  zmx = zmn
  DO  l = 2, nn
    IF (x(l) < xmn) xmn = x(l)
    IF (x(l) > xmx) xmx = x(l)
    IF (y(l) < ymn) ymn = y(l)
    IF (y(l) > ymx) ymx = y(l)
    IF (z(l) < zmn) zmn = z(l)
    IF (z(l) > zmx) zmx = z(l)
  END DO
  xyzmin(1) = xmn
  xyzmin(2) = ymn
  xyzmin(3) = zmn
  
! COMPUTE CELL DIMENSIONS AND TEST FOR ZERO AREA.
  
  delx = (xmx-xmn) / REAL(nnr)
  dely = (ymx-ymn) / REAL(nnr)
  delz = (zmx-zmn) / REAL(nnr)
  xyzdel(1) = delx
  xyzdel(2) = dely
  xyzdel(3) = delz
  IF (delx == 0. .OR. dely == 0. .OR. delz == 0.) GO TO 60
  
! INITIALIZE LCELL.
  
  lcell(1:nnr,1:nnr,1:nnr) = 0
  
! LOOP ON NODES, STORING INDICES IN LCELL AND LNEXT.
  
  np1 = nn + 1
  DO  ll = 1, nn
    lb = np1 - ll
    i = (x(lb)-xmn)/delx + 1
    IF (i > nnr) i = nnr
    j = (y(lb)-ymn)/dely + 1
    IF (j > nnr) j = nnr
    k = (z(lb)-zmn)/delz + 1
    IF (k > nnr) k = nnr
    l = lcell(i,j,k)
    lnext(lb) = l
    IF (l == 0) lnext(lb) = lb
    lcell(i,j,k) = lb
  END DO
  
! NO ERRORS ENCOUNTERED
  
  ier = 0
  RETURN
END IF

! INVALID INPUT PARAMETER

ier = 1
RETURN

! NONPOSITIVE XYZDEL COMPONENT

60 ier = 2
RETURN
END SUBROUTINE store3

END MODULE qshep3d



PROGRAM qs3test

!                           QS3TEST

!   THIS PROGRAM TESTS THE SCATTERED DATA INTERPOLATION
! PACKAGE QSHEP3D BY PRINTING THE MAXIMUM ERRORS ASSOCIATED
! WITH INTERPOLATED VALUES AND GRADIENTS ON A 5 BY 5 BY 5
! UNIFORM GRID IN THE UNIT CUBE.  THE DATA SET CONSISTS
! OF 64 NODES WITH DATA VALUES TAKEN FROM A QUADRATIC FUNC-
! TION FOR WHICH THE METHOD IS EXACT.  THE RATIO OF MAXIMUM
! INTERPOLATION ERROR RELATIVE TO THE MACHINE PRECISION IS
! ALSO PRINTED.  THIS SHOULD BE O(1).  THE INTERPOLATED
! VALUES FROM QS3VAL AND QS3GRD ARE COMPARED FOR AGREEMENT.

USE qshep3d
IMPLICIT NONE
INTEGER  :: i, ier, j, k, l, lcell(3,3,3), lnext(64)
REAL     :: a(9,64), eps, eq, eqx, eqy, eqz, f(64), px, py, pz,  &
            q, q1, qx, qy, qz, rmax, rq, rsq(64), x(64), xyzmin(3),  &
            xyzdel(3), y(64), yl, z(64), zl
REAL     :: p(5)

! QSHEP3 PARAMETERS AND LOGICAL UNIT FOR OUTPUT

INTEGER, PARAMETER  :: lout = 6, n = 64, nq = 17, nr = 3, nw = 32

! GENERATE A 4 BY 4 BY 4 GRID OF NODES IN THE UNIT CUBE.

l = 0
DO  k = 1, 4
  zl = (k-1) / 3.
  DO  j = 1, 4
    yl = (j-1) / 3.
    DO  i = 1, 4
      l = l + 1
      x(l) = (i-1) / 3.
      y(l) = yl
      z(l) = zl
    END DO
  END DO
END DO

! COMPUTE THE DATA VALUES.

DO  l = 1, n
  f(l) = fq(x(l),y(l),z(l))
END DO

! COMPUTE PARAMETERS DEFINING THE INTERPOLANT Q.

CALL qshep3(n,x,y,z,f,nq,nw,nr,lcell,lnext,xyzmin,xyzdel,rmax,rsq,a,ier)
IF (ier == 0) THEN

! GENERATE A 5 BY 5 BY 5 UNIFORM GRID OF INTERPOLATION POINTS (P(I),P(J),P(K))
!   IN THE UNIT CUBE.  THE EIGHT CORNERS COINCIDE WITH NODES.

  DO  i = 1, 5
    p(i) = (i-1) / 4.
  END DO

! COMPUTE THE MACHINE PRECISION EPS.

  eps = EPSILON(1.0)

! COMPUTE INTERPOLATION ERRORS AND TEST FOR AGREEMENT IN THE
!   Q VALUES RETURNED BY QS3VAL AND QS3GRD.

  EQ = 0.
  eqx = 0.
  eqy = 0.
  eqz = 0.
  DO  k = 1, 5
    pz = p(k)
    DO  j = 1, 5
      py = p(j)
      DO  i = 1, 5
        px = p(i)
        q1 = qs3val(px,py,pz,n,x,y,z,f,nr,lcell,lnext,xyzmin,  &
                    xyzdel,rmax,rsq,a)
        CALL qs3grd(px,py,pz,n,x,y,z,f,nr,lcell,lnext,xyzmin,  &
                    xyzdel,rmax,rsq,a,q,qx,qy,qz,ier)
        IF (ier /= 0) GO TO 100
        IF (ABS(q1-q) > 3.*ABS(q)*eps) GO TO 110
        EQ = MAX(EQ,ABS(fq(px,py,pz)-q))
        eqx = MAX(eqx,ABS(fx(px,py,pz)-qx))
        eqy = MAX(eqy,ABS(fy(px,py,pz)-qy))
        eqz = MAX(eqz,ABS(fz(px,py,pz)-qz))
      END DO
    END DO
  END DO

! PRINT ERRORS AND THE RATIO EQ/EPS.

  rq = EQ / eps
  WRITE (lout,5000)
  WRITE (lout,5100) EQ, rq
  WRITE (lout,5200) eqx
  WRITE (lout,5300) eqy
  WRITE (lout,5400) eqz
  STOP
END IF

! ERROR IN QSHEP3

WRITE (lout,5500) ier
STOP

! ERROR IN QS3GRD

100 WRITE (lout,5600) ier
STOP

! VALUES RETURNED BY QS3VAL AND QS3GRD DIFFER BY A RELATIVE
!   AMOUNT GREATER THAN 3*EPS.

110 WRITE (lout,5700) q1, q
STOP

5000 FORMAT (///' MAXIMUM ABSOLUTE ERRORS IN THE ',  &
             'INTERPOLANT Q AND PARTIAL'/' ',  &
             'DERIVATIVES (QX,QY,QZ) RELATIVE TO MACHINE PRECISION EPS'//  &
             t12, 'FUNCTION   MAX ERROR   MAX ERROR/EPS'/)
5100 FORMAT (t15, 'Q       ', e9.3, '       ', f4.2)
5200 FORMAT (t15, 'QX      ', e9.3)
5300 FORMAT (t15, 'QY      ', e9.3)
5400 FORMAT (t15, 'QZ      ', e9.3)
5500 FORMAT (//' *** ERROR IN QSHEP3 -- IER =',i2,' ***')
5600 FORMAT (//' *** ERROR IN QS3GRD -- IER =',i2,' ***')
5700 FORMAT (//' *** ERROR -- INTERPOLATED VALUES ',  &
             'Q1 (QS3VAL) AND Q2 (QS3GRD) DIFFER --'//  &
             t7, 'Q1 = ', e21.14, '     Q2 = ', e21.14)

CONTAINS

! QUADRATIC TEST FUNCTION AND PARTIAL DERIVATIVES

! fq(xx,yy,zz) = ((xx+2.*yy+3.*zz)/6.) ** 2
! fx(xx,yy,zz) = (xx+2.*yy+3.*zz) / 18.
! fy(xx,yy,zz) = (xx+2.*yy+3.*zz) / 9.
! fz(xx,yy,zz) = (xx+2.*yy+3.*zz) / 6.

FUNCTION fq(xx, yy, zz) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy, zz
REAL              :: fn_val

fn_val = ((xx+2.*yy+3.*zz)/6.) ** 2
RETURN
END FUNCTION fq



FUNCTION fx(xx, yy, zz) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy, zz
REAL              :: fn_val

fn_val = (xx+2.*yy+3.*zz) / 18.
RETURN
END FUNCTION fx



FUNCTION fy(xx, yy, zz) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy, zz
REAL              :: fn_val

fn_val = (xx+2.*yy+3.*zz) / 9.
RETURN
END FUNCTION fy



FUNCTION fz(xx, yy, zz) RESULT(fn_val)
REAL, INTENT(IN)  :: xx, yy, zz
REAL              :: fn_val

fn_val = (xx+2.*yy+3.*zz) / 6.
RETURN
END FUNCTION fz

END PROGRAM qs3test
