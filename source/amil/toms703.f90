!      ALGORITHM 703 , COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 2, JUNE, 1992, PP. 156-158.

! This ELF90-compatible version by Alan Miller
! amiller@bigpond.net.au
! http://users.bigpond.net.au/amiller
! Latest revision - 11 August 2003

MODULE common703
IMPLICIT NONE

! Contains the COMMON declarations from the Fortran 77 version,
! and the precision declaration.

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

!   COMMON BLOCKS OF INTEREST TO THE USER.

!   /MEBDF2/

!   HUSED     LAST STEPSIZE SUCCESSFULLY USED BY THE INTEGRATOR
!   NQUSED    ORDER LAST SUCCESSFULLY USED
!   NSTEP     NUMBER OF STEPS TAKEN SO FAR
!   NFE       NUMBER OF FUNCTION EVALUATIONS USED SO FAR
!   NJE       NUMBER OF JACOBIAN EVALUATIONS USED SO FAR
!   NDEC      NUMBER OF LU DECOMPOSITIONS USED SO FAR
!   NBSOL     NUMBER OF 'BACKSOLVES' USED SO FAR
!   NPSET     NUMBER OF TIMES A NEW COEFFICIENT MATRIX HAS BEEN FORMED SO FAR
!   NCOSET    NUMBER OF TIMES THE ORDER OF THE METHOD USED HAS BEEN
!             CHANGED SO FAR
!   MAXORD    THE MAXIMUM ORDER USED SO FAR IN THE INTEGRATION.

!   PUT THE FOLLOWING STATEMENT IN THE MAIN PROGRAM:
!        USE common703

! BLOCK DATA
!     MACHINE DEPENDENT CONSTANTS ARE SET HERE
!     UROUND  THE SMALLEST +VE NUMBER SUCH THAT 1.0 + UROUND > 1.0
!     EPSJAC  THIS IS USED IN FORMING A NUMERICAL JACOBIAN AND IS
!                EPSJAC = SQRT( UROUND )

! COMMON /machne/uround,epsjac

REAL (dp), PARAMETER :: uround = EPSILON(1.0_dp)
REAL (dp), SAVE      :: epsjac

! COMMON /fails/nfail1,nfail2,nt
! COMMON /mebdf1/t,h,hmin,hmax,epscpy,ncopy,mfcopy,kflag,jstart,toutcp
! COMMON /mebdf2/hused,nqused,nstep,nfe,nje,ndec,nbsol,npset,ncoset,maxord
! COMMON /mebdf3/iweval,newpar,nrenew,jsinup,jsnold,idoub,jchang,l,nq
! COMMON /mebdf4/m1,m2,meqc1,meqc2,mq1tmp,mq2tmp,isamp
! COMMON /mebdf5/qi,rc,rh,rmax,told
! COMMON /mebdf6/crate1,crate2,tcrat1,tcrat2,avnew2,avold2,avnewj,avoldj
! COMMON /mebdfl/cfail,jnewim,sample
! COMMON /stpsze/hstpsz,top

INTEGER, SAVE   :: idoub,isamp,iweval,jchang,jsinup,jsnold,jstart,kflag,l,  &
                   m1,m2,maxord,meqc1,meqc2,mfcopy,mq1tmp,mq2tmp,nbsol,     &
                   ncopy,ncoset,ndec,newpar,nfail1,nfail2,nfe,nje,npset,nq, &
                   nqused,nrenew,nstep,nt
LOGICAL, SAVE   :: cfail,jnewim,sample
REAL (dp), SAVE :: avnew2,avold2,avnewj,avoldj,epscpy,crate1,crate2,h,hmin, &
                   hmax,hstpsz(2,14),hused,qi,rc,rh,rmax,t,tcrat1,tcrat2,   &
                   told,top,toutcp

END MODULE common703



MODULE toms703
USE common703
IMPLICIT NONE

INTERFACE
  SUBROUTINE pderv(t,y,pw)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:,:)
    REAL (dp), INTENT(OUT)  :: pw(:,:)     ! pw(4,4)
  END SUBROUTINE pderv
END INTERFACE


CONTAINS


SUBROUTINE mebdf(n, t0, h0, y0, tout, tend, eps, mf, INDEX, lout, fcn)

IMPLICIT NONE
INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: t0
REAL (dp), INTENT(IN OUT) :: h0
REAL (dp), INTENT(IN OUT) :: y0(:)
REAL (dp), INTENT(IN)     :: tout
REAL (dp), INTENT(IN)     :: tend
REAL (dp), INTENT(IN)     :: eps
INTEGER, INTENT(IN)       :: mf
INTEGER, INTENT(IN OUT)   :: INDEX
INTEGER, INTENT(IN)       :: lout

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

!     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     THIS IS THE JULY 1990 VERSION OF OVDRIV, A PACKAGE FOR THE SOLUTION OF
!     THE INITIAL VALUE PROBLEM FOR SYSTEMS OF ORDINARY DIFFERENTIAL EQUATIONS
!     DY/DT = F(Y,T),    Y=(Y(1),Y(2),Y(3), . . . ,Y(N))
!     SUBROUTINE OVDRIV IS A DRIVER ROUTINE FOR THIS PACKAGE

!                    REFERENCES

!     1.  J.R. CASH, ON THE INTEGRATION OF STIFF SYSTEMS OF O.D.E.S
!         USING EXTENDED BACKWARD DIFFERENTIATION FORMULAE,
!         NUMER. MATH. 34, 235-246 (1980)
!     2.  THE INTEGRATION OF STIFF INITIAL VALUE PROBLEMS IN O.D.E.S
!         USING MODIFIED EXTENDED BACKWARD DIFFERENTIATION FORMULAE,
!         COMP. AND MATHS. WITH APPLICS., 9, 645-657, (1983).
!     3.  J.R. CASH, STABLE RECURSIONS WITH APPLICATIONS TO THE
!         NUMERICAL SOLUTION OF STIFF SYSTEMS, ACADEMIC PRESS.(1979)
!     4.  A.C. HINDMARSH, GEAR.. ORDINARY DIFFERENTIAL EQUATION
!         SYSTEM SOLVER, UCID-30001 REV. 3, LAWRENCE LIVERMORE
!         LABORATORY, P.O. BOX 808, LIVERMORE, CA 94550, DEC. 1974.
!     5.  J.R. CASH AND S. CONSIDINE, AN MEBDF CODE FOR STIFF INITIAL
!         VALUE PROBLEMS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE, 1992.

!     ----------------------------------------------------------------
!     OVDRIV IS TO BE CALLED ONCE FOR EACH OUTPUT VALUE OF T, AND
!     IN TURN MAKES REPEATED CALLS TO THE CORE INTEGRATOR STIFF.

!     THE INPUT PARAMETERS ARE ..
!     N     =  THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
!     T0    =  THE INITIAL VALUE OF T, THE INDEPENDENT VARIABLE
!              (USED ONLY ON THE FIRST CALL)
!     H0    =  THE NEXT STEP SIZE IN T (USED FOR INPUT ONLY ON THE FIRST CALL)
!     Y0    =  A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF Y
!              (USED FOR INPUT ONLY ON FIRST CALL)
!     TEND  =  END OF THE RANGE OF INTEGRATION.
!     TOUT  =  THE VALUE OF T AT WHICH OUTPUT IS DESIRED NEXT.
!              INTEGRATION WILL NORMALLY GO SLIGHTLY BEYOND TOUT AND THE
!              PACKAGE WILL INTERPOLATE TO T = TOUT MUST HAVE TEND >= TOUT IF
!              H0 IS POSITIVE AND TOUT >= TEND OTHERWISE.
!     EPS   =  THE RELATIVE ERROR BOUND (USED ONLY ON THE FIRST CALL,
!              UNLESS INDEX = -1).   SINGLE STEP ERROR ESTIMATES DIVIDED BY
!              YMAX(I) WILL BE KEPT LESS THAN EPS IN ROOT-MEAN-SQUARE NORM
!              (I.E. EUCLIDEAN NORM DIVIDED BY SQRT(N) ).
!              THE VECTOR YMAX OF WEIGHTS IS COMPUTED IN OVDRIV.
!              INITIALLY YMAX(I) IS ABS(Y(I)) WITH A DEFAULT VALUE OF 1 IF
!              Y(I)=0 INITIALLY.
!              THEREAFTER, YMAX(I) IS THE LARGEST VALUE OF ABS(Y(I))
!              SEEN SO FAR, OR THE INITIAL VALUE YMAX(I) IF THAT IS LARGER
!     MF    =  THE METHOD FLAG.  AT PRESENT ONLY MF=21 OR 22 IS ALLOWED
!              WHICH IS EXTENDED BACKWARD DIFFERENTIATION FORMULAE USING THE
!              CHORD METHOD WITH ANALYTIC OR NUMERICAL JACOBIAN.
!              THE USER NEEDS TO SPECIFY SUBROUTINE PDERV IF MF=21 (SEE BELOW)
!     INDEX =  INTEGER USED ON INPUT TO INDICATE TYPE OF CALL,
!              1   THIS IS THE FIRST CALL FOR THE PROBLEM
!                  INDEX HAS TO BE SET = 0 AFTER FIRST CALL.
!              0   THIS IS NOT THE FIRST CALL FOR THIS PROBLEM
!                  AND INTEGRATION IS TO CONTINUE UP TO T=TOUT OR BEYOND.
!             -1   THIS IS NOT THE FIRST CALL FOR THE PROBLEM,
!                  AND THE USER HAS RESET N, EPS, AND/OR MF.
!              2   SAME AS 0 EXCEPT THAT TOUT HAS TO BE HIT
!                  EXACTLY (NO INTERPOLATION IS DONE).
!                  ASSUMES TOUT >= THE CURRENT T.
!              3   ONE STEP MODE. SAME AS 0 EXCEPT CONTROL
!                  RETURNS TO CALLING PROGRAM AFTER ONE STEP.
!              SINCE THE NORMAL OUTPUT VALUE OF INDEX IS 0,
!              IT NEED NOT BE RESET FOR NORMAL CONTINUATION.
!     LOUT  =  LOGICAL OUTPUT DEVICE.

!     AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURED AND A NORMAL
!     CONTINUATION IS DESIRED, SIMPLY RESET TOUT AND CALL AGAIN.
!     ALL OTHER PARAMETERS WILL BE READY FOR THE NEXT CALL.
!     A CHANGE OF PARAMETERS WITH INDEX = -1 CAN BE MADE AFTER
!     EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN.

!     THE OUTPUT PARAMETERS ARE..
!     T0    =  THE VALUE OF T WHICH RELATES TO THE CURRENT SOLUTION POINT Y0()
!     H0    =  THE STEPSIZE H USED LAST, WHETHER SUCCESSFULLY OR NOT.
!     Y0    =  THE COMPUTED VALUES OF Y AT T = TOUT
!     TOUT  =  UNCHANGED FROM ITS INPUT VALUE.
!     INDEX =  INTEGER USED ON OUTPUT TO INDICATE RESULTS, WITH
!              THE FOLLOWING VALUES AND MEANINGS..

!      3   ONE STEP MODE IS BEING USED AND A SUCCESSFUL STEP HAS BEEN
!          COMPLETED.  THE USER SHOULD CALL MEBDF AGAIN WITH INDEX = 3.
!      1   ONE SUCCESSFUL STEP HAS BEEN COMPLETED AND THE ONE STEP MODE
!          IS NOT BEING USED.  THE USER MUST NOW SET INDEX = 0 AND CALL
!          MEBDF AGAIN.
!      0   INTEGRATION WAS COMPLETED TO TOUT OR BEYOND.

!     -1   THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE ERROR TEST
!          EVEN AFTER REDUCING H BY A FACTOR OF 1.E10 FROM ITS INITIAL VALUE.

!     -2   AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS HALTED EITHER BY
!          REPEATED ERROR TEST FAILURES OR BY A TEST ON EPS.
!          TOO MUCH ACCURACY HAS BEEN REQUESTED.

!     -3   THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE CORRECTOR
!          CONVERGENCE EVEN AFTER REDUCING H BY A FACTOR OF 1.E10 FROM ITS
!          INITIAL VALUE.

!     -4   IMMEDIATE HALT BECAUSE OF ILLEGAL VALUES OF INPUT PARAMETERS.
!          SEE PRINTED MESSAGE.

!     -5   INDEX WAS -1 ON INPUT, BUT THE DESIRED CHANGES OF PARAMETERS WERE
!          NOT IMPLEMENTED BECAUSE TOUT WAS NOT BEYOND T.
!          INTERPOLATION AT T = TOUT WAS PERFORMED AS ON A NORMAL RETURN.
!          TO TRY AGAIN, SIMPLY CALL AGAIN WITH INDEX = -1 AND A NEW TOUT.

!     IN ADDITION TO OVDRIV, THE FOLLOWING ROUTINES ARE PROVIDED IN THE PACKAGE..

!     INTERP( - )   INTERPOLATES TO GET THE OUTPUT VALUES
!                   AT T=TOUT FROM THE DATA IN THE Y ARRAY.
!     STIFF( - )    IS THE CORE INTEGRATOR ROUTINE.  IT PERFORMS A
!                   SINGLE STEP AND ASSOCIATED ERROR CONTROL.
!     COSET( - )    SETS COEFFICIENTS FOR BACKWARD DIFFERENTIATION
!                   SCHEMES FOR USE IN THE CORE INTEGRATOR.
!     PSET( - )     COMPUTES AND PROCESSES THE JACOBIAN MATRIX J = DF/DY
!     DEC( - )      PERFORMS AN LU DECOMPOSITION ON A MATRIX.
!     SOL( - )      SOLVES LINEAR SYSTEMS A*X = B AFTER DEC
!                   HAS BEEN CALLED FOR THE MATRIX A

!     THE FOLLOWING ROUTINES ARE TO BE SUPPLIED BY THE USER..

!     FCN(T,Y,YDOT)  COMPUTES THE FUNCTION YDOT = F(Y,T), THE
!                         RIGHT HAND SIDE OF THE O.D.E.
!                         HERE Y AND YDOT ARE VECTORS OF LENGTH N
!     PDERV(T,Y,PD)  COMPUTES THE N*N JACOBIAN MATRIX OF PARTIAL DERIVATIVES
!                         AND STORES IT IN PD AS AN N BY N ARRAY.
!                         PD(I,J) IS TO BE SET TO THE PARTIAL DERIVATIVE OF
!                         YDOT(I) WITH RESPECT TO Y(J).
!                         PDERV IS CALLED ONLY IF MITER = 1.
!                         OTHERWISE A DUMMY ROUTINE CAN BE SUBSTITUTED.

!     THE DIMENSIONS OF PW BELOW MUST BE AT LEAST N x N.  THE DIMENSIONS OF
!     YMAX, ERROR, SAVE1, SAVE2, IPIV AND THE FIRST DIMENSION OF Y SHOULD ALL
!     BE AT LEAST N.

!  **************************************************************************

epsjac = SQRT(uround)

IF (INDEX == 1) THEN
  IF (n <= 0) THEN
    WRITE (lout,9020) n
    INDEX = -4
    
  END IF
  
  IF (INDEX < 0) RETURN
END IF

CALL ovdriv(n, t0, h0, y0, tout, tend, eps, mf, INDEX, lout, fcn)

RETURN

9020 FORMAT (// ' ***** ERROR ***** INTEGRATION HALTED IN MEBDF'//  &
                '   >>> ILLEGAL VALUE FOR NUMBER OF EQUATIONS <<< ',  &
                '                     WITH N = ',i6)

END SUBROUTINE mebdf


SUBROUTINE ovdriv(n, t0, h0, y0, tout, tend, eps, mf, INDEX, lout, fcn)

IMPLICIT NONE
INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: t0
REAL (dp), INTENT(IN OUT) :: h0
REAL (dp), INTENT(IN OUT) :: y0(:)
REAL (dp), INTENT(IN)     :: tout
REAL (dp), INTENT(IN)     :: tend
REAL (dp), INTENT(IN)     :: eps
INTEGER, INTENT(IN)       :: mf
INTEGER, INTENT(IN OUT)   :: INDEX
INTEGER, INTENT(IN)       :: lout

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

!     START OF PROGRAM PROPER

!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i, kgo, nhcut
REAL (dp) :: ayi, d
!     ..
! Local arrays
REAL (dp) :: y(n,15), yhold(n,15), pwcopy(n*n)
! If n > 50, increase the 50's in the two lines below
REAL (dp), SAVE  :: ymax(50), pw(50,50)
INTEGER, SAVE    :: ipiv(50)
!     ..
IF (INDEX == 0) THEN
!        I.E. NORMAL CONTINUATION OF INTEGRATION
  t0=t
  hmax = ABS(tend-t0)*10.0_dp
  IF ((t-tout)*h >= 0.0_dp) THEN
!           HAVE OVERSHOT THE ENDPOINT, SO INTERPOLATE
    CALL interp(n, jstart, h, t, y, tout, y0)
    INDEX = kflag
    t0 = tout
    h0 = h
    RETURN
    
  END IF
  
ELSE IF (INDEX == 2) THEN
!        I.E. CONTINUING INTEGRATION BUT WISH TO HIT TOUT
  t0 = t
  hmax = ABS(tend-t0)*10.0_dp
  IF (((t+h)-tout)*h > 0.0_dp) THEN
!           WE HAVE ALREADY OVERSHOT THE END OR WE WILL
!           DO SO ON THE NEXT STEP
    IF (((t-tout)*h >= 0.0_dp) .OR. (ABS(t-tout) <=  &
          100.0_dp*uround*hmax)) THEN
!              HAVE OVERSHOT, SO INTERPOLATE
      CALL interp(n, jstart, h, t, y, tout, y0)
      t0 = tout
      h0 = h
      INDEX = kflag
      RETURN
      
    ELSE
!              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
!              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
      h = (tout-t)* (1.0_dp-4.0_dp*uround)
      jstart = -1
    END IF
    
  END IF
  
ELSE IF (INDEX == -1) THEN
!        NOT FIRST CALL BUT PARAMETERS RESET
  t0 = t
  IF ((t-tout)*h >= 0.0_dp) THEN
!           HAVE OVERSHOT TOUT
    WRITE (lout,9080) t,tout,h
    CALL interp(n, jstart, h, t, y, tout, y0)
    h0 = h
    t0 = tout
    INDEX = -5
    RETURN
    
  ELSE
    jstart = -1
    ncopy = n
    epscpy = eps
    mfcopy = mf
  END IF
  
ELSE IF (INDEX == 3) THEN
  t0 = t
  IF ((t-tout)*h >= 0.0_dp) THEN
!           HAVE OVERSHOT,SO INTERPOLATE
    CALL interp(n, jstart, h, t, y, tout, y0)
    INDEX = kflag
    t0 = tout
    h0 = h
    RETURN
    
  END IF
  
ELSE
!        INDEX SHOULD BE 1 AND THIS IS THE FIRST CALL FOR THIS PROBLEM
!        CHECK THE ARGUMENTS THAT WERE PASSED FOR CORRECTNESS
  IF (INDEX /= 1) THEN
!           VALUE OF INDEX NOT ALLOWED
    WRITE (lout,9070) INDEX
    INDEX = -4
  END IF
  
  IF (eps <= 0.0_dp) THEN
!           ILLEGAL VALUE FOR RELATIVE ERROR TOLERENCE
    WRITE (lout,9040)
    INDEX = -4
  END IF
  
  IF (n <= 0) THEN
!           ILLEGAL VALUE FOR THE NUMBER OF EQUATIONS
    WRITE (lout,9050)
    INDEX = -4
  END IF
  
  IF ((t0-tout)*h0 >= 0.0_dp) THEN
!           PARAMETERS FOR INTEGRATION ARE ILLEGAL
    WRITE (lout,9060)
    INDEX = -4
  END IF
  
  IF((tend-tout)*h0 < -100.0_dp*uround) THEN
!      TOUT IS BEYOND TEND.
    WRITE(lout,9110) tout,tend
    INDEX = -4
  END IF
  IF ((mf /= 21) .AND. (mf /= 22)) THEN
!           ILLEGAL VALUE FOR METHOD FLAG
    WRITE (lout,9090) mf
    INDEX = -4
  END IF
  
  IF (INDEX /= 1) THEN
    RETURN
    
  ELSE
!           THE INITIAL PARAMETERS ARE O.K. SO INITIALISE EVERYTHING
!           ELSE NECESSARY FOR THE INTEGRATION.
!           IF VALUES OF YMAX OTHER THAN THOSE SET BELOW ARE DESIRED,
!           THEY SHOULD BE SET HERE. ALL YMAX(I) MUST BE POSITIVE. IF
!           VALUES FOR HMIN OR HMAX, THE BOUNDS ON ABS(H), OTHER THAN
!           THOSE BELOW ARE DESIRED, THEY SHOULD BE SET BELOW
    DO i = 1,n
      ymax(i) = ABS(y0(i))
      IF (ymax(i) == 0.0_dp) ymax(i) = 1.0_dp
      y(i,1) = y0(i)
    END DO
    ncopy = n
    t = t0
    h = h0
    hmin = ABS(h0)
    hmax = ABS(t0-tend)*10.0_dp
    epscpy = eps
    mfcopy = mf
    jstart = 0
    nhcut = 0
  END IF
  
END IF
!     <<<<<<<<<<<<<<<<<
!     <  TAKE A STEP  >
!     <<<<<<<<<<<<<<<<<
20 IF ((t+h) == t) THEN
  WRITE (lout,9000)
  nt = nt + 1
END IF

CALL stiff(eps, h, hmax, hmin, jstart, kflag, mf, t, tend, y, n, ymax,  &
           pw, pwcopy, yhold, ipiv, lout, fcn)

kgo = 1 - kflag
IF (kgo == 1) THEN
!        NORMAL RETURN FROM STIFF
  GO TO 30
  
ELSE IF (kgo == 2) THEN
!        COULD NOT ACHIEVE ERROR WITH HMIN
!        SO CHOP HMIN IF WE HAVEN'T DONE SO 10 TIMES
  GO TO 60
  
ELSE IF (kgo == 3) THEN
!        ERROR SMALLER THAN CAN BE HANDLED FOR THIS PROBLEM
  WRITE (lout, 9010) t, h
  GO TO 70
  
ELSE IF (kgo == 4) THEN
!        COULD NOT ACHIEVE CONVERGENCE WITH HMIN
  WRITE (lout, 9030) t
  GO TO 60
  
END IF

! ---------------------------------------------------------------------
!     NORMAL RETURN FROM THE INTEGRATOR.

!     THE WEIGHTS YMAX(I) ARE UPDATED.  IF DIFFERENT VALUES ARE DESIRED,
!     THEY SHOULD BE SET HERE.  A TEST IS MADE FOR EPS BEING TOO SMALL
!     FOR THE MACHINE PRECISION.

!     ANY OTHER TESTS OR CALCULATIONS THAT ARE REQUIRED AFTER EVERY
!     STEP SHOULD BE INSERTED HERE.

!     IF INDEX = 3, Y0 IS SET TO THE CURRENT Y VALUES ON RETURN.
!     IF INDEX = 2, H IS CONTROLLED TO HIT TOUT (WITHIN ROUNDOFF
!     ERROR), AND THEN THE CURRENT Y VALUES ARE PUT IN Y0 ON RETURN.
!     FOR ANY OTHER VALUE OF INDEX, CONTROL RETURNS TO THE INTEGRATOR
!     UNLESS TOUT HAS BEEN REACHED.  THEN INTERPOLATED VALUES OF Y ARE
!     COMPUTED AND STORED IN Y0 ON RETURN.
!     IF INTERPOLATION IS NOT DESIRED, THE CALL TO INTERP SHOULD BE
!     REMOVED AND CONTROL TRANSFERRED TO STATEMENT 500 INSTEAD OF 520.
! --------------------------------------------------------------------
30 d = 0.0_dp
DO i = 1,n
  ayi = ABS(y(i,1))
  ymax(i) = MAX(ymax(i),ayi)
  d = d + (ayi/ymax(i))**2
END DO
d = d* (uround/eps)**2
IF (d > DBLE(n)) THEN
!        EPS SMALLER THAN MACHINE PRECISION
  WRITE (lout,9020) t
  kflag = -2
  GO TO 70
  
END IF
IF (INDEX == 3 .OR. INDEX == 1) GO TO 70

IF (ABS(t-tout) <= ABS(15.0_dp*uround*tout)) THEN
!        EFFECTIVELY WE HAVE HIT TOUT
  INDEX = kflag
  t0 = tout
  DO i = 1,n
    y0(i) = y(i,1)
  END DO
  h0 = h
  RETURN
  
END IF

IF (INDEX == 2) THEN
!        CONTINUING INTEGRATION BUT MUST HIT TOUT EXACTLY
  IF (((t+h)-tout)*h > 0.0_dp) THEN
!           WE HAVE ALREADY OVERSHOT THE END OR WE WILL DO
!           SO ON THE NEXT STEP
    IF (((t-tout)*h >= 0.0_dp) .OR. (ABS(t-tout) <=  &
          100.0_dp*uround*hmax)) THEN
!              HAVE OVERSHOT, SO INTERPOLATE
      CALL interp(n, jstart, h, t, y, tout, y0)
      t0 = tout
      h0 = h
      INDEX = kflag
      RETURN
      
    ELSE
!              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
!              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
      h = (tout-t)* (1.0_dp-4.0_dp*uround)
      jstart = -1
    END IF
    
  END IF
  
ELSE IF ((t-tout)*h >= 0.0_dp) THEN
!        HAVE OVERSHOT, SO INTERPOLATE
  CALL interp(n, jstart, h, t, y, tout, y0)
  INDEX = kflag
  h0 = h
  t0 = tout
  RETURN
  
END IF

GO TO 20

! -------------------------------------------------------------------
!    ON AN ERROR RETURN FROM THE INTEGRATOR, AN IMMEDIATE RETURN OCCURS
!    IF KFLAG = -2, AND RECOVERY ATTEMPTS ARE MADE OTHERWISE.
!    H AND HMIN ARE REDUCED BY A FACTOR OF .1 UP TO 10 TIMES BEFORE GIVING UP.
! --------------------------------------------------------------------
60 IF (nhcut == 10) THEN
!        HAVE REDUCED H TEN TIMES
  WRITE (lout,9100)
  GO TO 70
  
END IF

nhcut = nhcut + 1
hmin = 0.1_dp*hmin
h = 0.1_dp*h
IF (nstep > 0) THEN
  jstart = -1
  
ELSE
!        THE INITIAL STEPSIZE WAS TOO LARGE. RESET CERTAIN VARIABLES
!        AND TRY AGAIN.
  nje = 0
  nfe = 1
  cfail = .true.
  newpar = 0
  mq1tmp = 0
  mq2tmp = 0
  meqc1 = 0
  meqc2 = 0
  tcrat1 = 0.0_dp
  tcrat2 = 0.0_dp
  crate1 = 1.0_dp
  crate2 = 1.0_dp
  nt = 0
  nstep = 0
  nbsol = 0
  npset = 0
  ncoset = 0
  ndec = 0
  jstart = -1
END IF

GO TO 20

70 IF(ABS(t-tout) > 1000.0_dp*uround) THEN
  y0(1:n) = y(1:n,1)
  t0 = t
  
ELSE
!        HAVE EITHER PASSED TOUT OR WE ARE EXTREMELY CLOSE TO IT
!        SO INTERPOLATE.
  CALL interp(n, jstart, h, t, y, tout, y0)
  t0 = tout
  INDEX = kflag
END IF

h0 = h
IF(kflag /= 0) INDEX = kflag
RETURN
! -------------------------- END OF SUBROUTINE OVDRIV -----------------
9000 FORMAT (' WARNING..  T + H = T ON NEXT STEP.')
9010 FORMAT (//' KFLAG = -2 FROM INTEGRATOR AT T = ',e16.8,'  H =',  &
    e16.8/ '  THE REQUESTED ERROR IS SMALLER THAN CAN BE HANDLED'//)
9020 FORMAT (//' INTEGRATION HALTED BY MEBDF AT T = ',e16.8/  &
    '  EPS TOO SMALL TO BE ATTAINED FOR THE MACHINE PRECISION', /)
9030 FORMAT (//' KFLAG = -3 FROM INTEGRATOR AT T = ',e16.8/  &
    '  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED'/)
9040 FORMAT (//' ILLEGAL INPUT.. EPS .LE. 0.'//)
9050 FORMAT (//' ILLEGAL INPUT.. N .LE. 0'//)
9060 FORMAT (//' ILLEGAL INPUT.. (T0-TOUT)*H >= 0.'//)
9070 FORMAT (//' ILLEGAL INPUT.. INDEX =',i5//)
9080 FORMAT (//' INDEX = -1 ON INPUT WITH (T-TOUT)*H >= 0.'/  &
    ' T =',e16.8,'   TOUT =',e16.8,'   H =',e16.8/  &
    ' INTERPOLATION WAS DONE AS ON NORMAL RETURN.'/  &
    ' DESIRED PARAMETER CHANGES WERE NOT MADE.')
9090 FORMAT (//' ILLEGAL INPUT.. METHOD FLAG, MF, = ',i6/  &
    '         ALLOWED VALUES ARE 21 OR 22',/)
9100 FORMAT (//' PROBLEM APPEARS UNSOLVABLE WITH GIVEN INPUT'/  &
    '         HMIN REDUCED BY A FACTOR OF 1.0E10'//)
9110 FORMAT(//' TOUT IS BEYOND TEND' /  &
    'TOUT =' ,e16.8, '    TEND = ', e16.8/  &
    'RESET TOUT,TEND SO THAT (TEND-TOUT)*H0 >=0' )
END SUBROUTINE ovdriv


SUBROUTINE interp(n, jstart, h, t, y, tout, y0)

IMPLICIT NONE
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: jstart
REAL (dp), INTENT(IN)   :: h
REAL (dp), INTENT(IN)   :: t
REAL (dp), INTENT(IN)   :: y(:,:)     ! y(n,15)
REAL (dp), INTENT(IN)   :: tout
REAL (dp), INTENT(OUT)  :: y0(n)

!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i, j, l
REAL (dp) :: s, s1
!     ..
DO i = 1,n
  y0(i) = y(i,1)
END DO
l = jstart + 2
s = (tout-t)/h
s1 = 1.0_dp
DO j = 2,l
  s1 = s1*(s + (j-2))/DBLE(j-1)
  DO i = 1,n
    y0(i) = y0(i) + s1*y(i,j)
  END DO
END DO
RETURN
! -------------- END OF SUBROUTINE INTERP ---------------------------
END SUBROUTINE interp


SUBROUTINE coset(nq, el, elst, tq, maxder)

IMPLICIT NONE
INTEGER, INTENT(IN)        :: nq
REAL (dp), INTENT(IN OUT)  :: el(:)
REAL (dp), INTENT(IN OUT)  :: elst(:)
REAL (dp), INTENT(OUT)     :: tq(:)
INTEGER, INTENT(OUT)       :: maxder

! -------------------------------------------------------------------------
!   COSET IS CALLED BY THE INTEGRATOR AND SETS THE COEFFICIENTS USED BY THE
!   CONVENTIONAL BACKWARD DIFFERENTIATION SCHEME.  THE VECTOR EL OF LENGTH
!   NQ+1 DETERMINES THE BASIC METHOD.  THE VECTOR ELST OF DIMENSION NQ+2
!   DETERMINES THE MODIFIED EXTENDED BACKWARD DIFFERENTIATION FORMULAE
!   COEFFICIENTS.  THE VECTOR TQ OF LENGTH 4 IS INVOLVED IN ADJUSTING THE
!   STEPSIZE IN RELATION TO THE TRUNCATION ERROR.  ITS VALUES ARE GIVEN BY
!   THE PERTST ARRAY.  THE VECTORS EL AND TQ BOTH DEPEND ON METH AND NQ.
!   COSET ALSO SETS MAXDER, THE MAXIMUM ORDER OF THE METHOD AVAILABLE.  THE
!   COEFFICIENTS IN PERTST NEED TO BE GIVEN TO ONLY ABOUT ONE PERCENT
!   ACCURACY.  THE ORDER IN WHICH THE GROUPS APPEAR BELOW IS COEFFICIENTS
!   FOR ORDER NQ-1, COEFFICIENTS FOR ORDER NQ, COEFFICIENTS FOR ORDER NQ+1.
! ------------------------------------------------------------------------
!     ..
!     .. LOCAL ARRAYS ..
REAL (dp), SAVE :: pertst(8,3) = RESHAPE( (/ 1.0_dp, 2.0_dp, 4.5_dp,    &
                   7.333_dp, 10.42_dp, 13.7_dp, 17.15_dp, 20.74_dp,     &
                   2.0_dp, 4.5_dp, 7.333_dp, 10.42_dp, 13.7_dp, 17.15_dp, &
                   20.74_dp, 24.46_dp, 4.5_dp, 7.333_dp, 10.42_dp,      &
                   13.7_dp, 17.15_dp, 20.74_dp, 24.46_dp, 1.0_dp /),    &
                   (/ 8, 3 /) )
!     ..
! -------------------------------------------------------------------
!     THE FOLLOWING COEFFICIENTS SHOULD BE DEFINED TO MACHINE ACCURACY.
!     THEIR DERIVATION IS GIVEN IN REFERENCE 2.
! -------------------------------------------------------------------
IF (nq > maxord) maxord = nq
maxder = 6
ncoset = ncoset + 1

SELECT CASE (nq)
  CASE (1)
    el(1) = 1.0_dp
    elst(1) = 1.5_dp
    elst(3) = -0.5_dp
  CASE (2)
    el(1) = 6.6666666666667D-01
    el(3) = 3.3333333333333D-01
    elst(1) = 9.5652173913043D-01
    elst(3) = 2.1739130434782D-01
    elst(4) = -1.7391304347826D-01
  CASE (3)
    el(1) = 5.4545454545455D-01
    el(3) = 4.5454545454545D-01
    el(4) = 1.8181818181818D-01
    elst(1) = 7.6142131979695D-01
    elst(3) = 3.2994923857868D-01
    elst(4) = 8.6294416243654D-02
    elst(5) = -9.1370558375634D-02
  CASE (4)
    el(1) = 0.48_dp
    el(3) = 0.52_dp
    el(4) = 0.28_dp
    el(5) = 0.12_dp
    elst(1) = 6.5733706517393D-01
    elst(3) = 4.0023990403838D-01
    elst(4) = 1.5793682526989D-01
    elst(5) = 4.4382247101159D-02
    elst(6) = -5.7576969212315D-02
  CASE (5)
    el(1) = 4.3795620437956D-01
    el(3) = 5.62043795620436D-01
    el(4) = 3.43065693430656D-01
    el(5) = 1.97080291970802D-01
    el(6) = 8.75912408759123D-02
    elst(1) = 5.9119243917152D-01
    elst(3) = 4.4902473356122D-01
    elst(4) = 2.1375427307460D-01
    elst(5) = 9.0421610027481503D-02
    elst(6) = 2.6409276761177D-02
    elst(7) = -4.0217172732757D-02
  CASE (6)
    el(1) = 4.08163265306120D-01
    el(3) = 5.91836734693874D-01
    el(4) = 3.87755102040813D-01
    el(5) = 2.51700680272107D-01
    el(6) = 1.49659863945577D-01
    el(7) = 6.80272108843534D-02
    elst(1) = 5.4475876041119D-01
    elst(3) = 4.8525549636077D-01
    elst(4) = 2.5789750131312D-01
    elst(5) = 1.3133738525800D-01
    elst(6) = 5.7677396763462D-02
    elst(7) = 1.7258197643881D-02
    elst(8) = -3.0014256771967D-02
  CASE (7)
    el(1) = 3.85674931129476D-01
    el(3) = 6.14325068870521D-01
    el(4) = 4.21487603305783D-01
    el(5) = 2.9292929292929D-01
    el(6) = 1.96510560146923D-01
    el(7) = 1.19375573921028D-01
    el(8) = 5.50964187327820D-02
    elst(1) = 5.0999746293734D-01
    elst(3) = 5.1345839935281D-01
    elst(4) = 2.9364346131937D-01
    elst(5) = 1.6664672120553D-01
    elst(6) = 8.8013735242353D-02
    elst(7) = 3.9571794884069D-02
    elst(8) = 1.2039080338722D-02
    elst(9) = -2.3455862290154D-02
END SELECT

tq(1:3) = pertst(nq,1:3)
tq(4) = 0.5_dp*tq(2) / nq
RETURN
! --------------------- END OF SUBROUTINE COSET ---------------------
END SUBROUTINE coset


SUBROUTINE pset(y, n, h, t, uround, epsjac, con, miter, ier, nrenew, ymax,  &
                save2, pw, pwcopy, wrkspc, ipiv, fcn)

IMPLICIT NONE
REAL (dp), INTENT(IN OUT)  :: y(:,:)     ! y(n,15)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: h
REAL (dp), INTENT(IN OUT)  :: t
REAL (dp), INTENT(IN)      :: uround
REAL (dp), INTENT(IN)      :: epsjac
REAL (dp), INTENT(IN)      :: con
INTEGER, INTENT(IN OUT)    :: miter
INTEGER, INTENT(OUT)       :: ier
INTEGER, INTENT(IN OUT)    :: nrenew
REAL (dp), INTENT(IN)      :: ymax(:)
REAL (dp), INTENT(IN)      :: save2(:)
REAL (dp), INTENT(OUT)     :: pw(:,:)
REAL (dp), INTENT(IN OUT)  :: pwcopy(:)
REAL (dp), INTENT(IN OUT)  :: wrkspc(:)
INTEGER, INTENT(IN OUT)    :: ipiv(:)

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

! -------------------------------------------------------------------
!     PSET IS CALLED BY STIFF TO COMPUTE AND PROCESS THE COEFFICIENT
!     MATRIX I - H*EL(1)*J  WHERE J IS AN APPROXIMATION TO
!     THE RELEVANT JACOBIAN.  THIS MATRIX IS THEN SUBJECTED TO LU
!     DECOMPOSITION IN PREPARATION FOR LATER SOLUTION OF LINEAR SYSTEMS
!     OF ALGEBRAIC EQUATIONS WITH LU AS THE COEFFICIENT MATRIX.  THE
!     MATRIX J IS FOUND BY THE USER-SUPPLIED ROUTINE PDERV IF MITER=1
!     OR BY FINITE DIFFERENCING IF MITER = 2.
!     IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION WITH
!     PSET USES THE FOLLOWING ..
!     EPSJAC = DSQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.
! *******************************************************************
!     THE ARGUMENT NRENEW IS USED TO SIGNAL WHETHER OR NOT
!     WE REQUIRE A NEW JACOBIAN TO BE CALCULATED.

!        IF NRENEW > 0 THEN WE REQUIRE A NEW J TO BE COMPUTED
!                  = 0 THEN USE A COPY OF THE LAST J COMPUTED
! *******************************************************************
!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i, j, j1, jjkk
REAL (dp) :: d, r, r0, tempry, yj

npset = npset + 1
IF (nrenew == 0) THEN
  pw(1:n,1:n) = RESHAPE( pwcopy(1:n*n), (/ n, n /) )
  pw = pw * con
  GO TO 70
END IF

IF (miter == 2) GO TO 30
CALL pderv(t,y,pw)
pwcopy(1:n*n) = RESHAPE( pw(1:n,1:n), (/ n*n /) )
pw = con * pw
nje = nje + 1
GO TO 70

30 nje = nje + 1
d = 0.0_dp
DO i = 1,n
  d = d + save2(i)**2
END DO
r0 = ABS(h)*SQRT(d)*1000.0_dp*uround
j1 = 0
DO j = 1,n
  yj = y(j,1)
  r = epsjac*ymax(j)
  r = MAX(r,r0)
  y(j,1) = y(j,1) + r
  d = con/r
  CALL fcn(t, y(:,1), wrkspc)
  DO i = 1,n
    jjkk = i + j1
    tempry = (wrkspc(i) - save2(i))
    pwcopy(jjkk) = tempry/r
    pw(i,j) = tempry*d
  END DO
  y(j,1) = yj
  j1 = j1 + n
END DO
nfe = nfe + n
70 DO i = 1,n
  pw(i,i) = pw(i,i) + 1.0_dp
END DO
CALL dec(n, pw, ipiv, ier)
!      NDEC = NDEC + 1
RETURN
! ---------------------- END OF SUBROUTINE PSET ---------------------
END SUBROUTINE pset


SUBROUTINE dec(n, a, ip, ier)

IMPLICIT NONE
INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: a(:,:)    ! a(ndim,n)
INTEGER, INTENT(OUT)      :: ip(:)
INTEGER, INTENT(OUT)      :: ier

! -------------------------------------------------------------------
!     MATRIX TRIANGULARISATION BY GAUSSIAN ELIMINATION
!     INPUT..
!     N    = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF ARRAY A.
!     A    = MATRIX TO BE TRIANGULARISED.
!     OUTPUT..
!     A(I,J),  I.LE.J = UPPER TRIANGULAR FACTOR, U.
!     A(I,J),  I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I-L.
!     IP(K), K.LT.N = INDEX OF KTH PIVOT ROW.
!     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
!     IER = 0 IF MATRIX IS NON-SINGULAR, OR K IF FOUND TO BE SINGULAR
!                  AT STAGE K.
!     USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!     DETERM(A) = IP(N)*A(1,1)*A(2,2)* . . . *A(N,N).
!     IF IP(N) = 0, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.

!     REFERENCE.
!     C.B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, C.A.C.M
!     15 (1972), P.274.
!     ------------------------------------------------------------------

!     ..
!     .. LOCAL SCALARS ..
INTEGER :: i, j, k, kp1, m, nm1
!     ..
ier = 0
ip(n) = 1
IF (n == 1) GO TO 70
nm1 = n - 1
DO k = 1,nm1
  kp1 = k + 1
  m = k
  DO i = kp1,n
    IF (ABS(a(i,k)) > ABS(a(m,k))) m = i
  END DO
  ip(k) = m
  t = a(m,k)
  IF (m == k) GO TO 20
  ip(n) = -ip(n)
  a(m,k) = a(k,k)
  a(k,k) = t
  20 IF (t == 0.0_dp) GO TO 80
  t = 1.0_dp/t
  DO i = kp1,n
    a(i,k) = -a(i,k)*t
  END DO
  DO j = kp1,n
    t = a(m,j)
    a(m,j) = a(k,j)
    a(k,j) = t
    IF (t == 0.0_dp) CYCLE
    DO i = kp1,n
      a(i,j) = a(i,j) + a(i,k)*t
    END DO
  END DO
END DO
70 k = n
IF (a(n,n) == 0.0_dp) GO TO 80
RETURN

80 ier = k
ip(n) = 0
RETURN
!     --------------------- END OF SUBROUTINE DEC ----------------------
END SUBROUTINE dec


SUBROUTINE sol(n, a, b, ip)

IMPLICIT NONE
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: a(:,:)     ! a(ndim,n)
REAL (dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(IN)        :: ip(:)

!     ..
!     .. LOCAL SCALARS ..
INTEGER :: i, k, kb, km1, kp1, m, nm1
!     ..
!     ------------------------------------------------------------------
!     SOLUTION OF LINEAR SYSTEM, A*X = B.
!     INPUT ..
!     N = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF MATRIX A.
!     A = TRIANGULARISED MATRIX OBTAINED FROM DEC.
!     B = RIGHT HAND SIDE VECTOR.
!     IP = PIVOT VECTOR OBTAINED FROM DEC.
!     DO NOT USE IF DEC HAS SET IER .NE. 0
!     OUTPUT..
!     B = SOLUTION VECTOR, X.
!     ------------------------------------------------------------------
IF (n == 1) GO TO 50
nm1 = n - 1
DO k = 1,nm1
  kp1 = k + 1
  m = ip(k)
  t = b(m)
  b(m) = b(k)
  b(k) = t
  DO i = kp1,n
    b(i) = b(i) + a(i,k)*t
  END DO
END DO
DO kb = 1,nm1
  km1 = n - kb
  k = km1 + 1
  b(k) = b(k)/a(k,k)
  t = -b(k)
  DO i = 1,km1
    b(i) = b(i) + a(i,k)*t
  END DO
END DO
50 b(1) = b(1)/a(1,1)
RETURN
!     ------------------------- END OF SUBROUTINE SOL ------------------
END SUBROUTINE sol


SUBROUTINE errors(n, tq, eps, edn, e, eup, bnd)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: tq(:)
REAL (dp), INTENT(IN)   :: eps
REAL (dp), INTENT(OUT)  :: edn
REAL (dp), INTENT(OUT)  :: e
REAL (dp), INTENT(OUT)  :: eup
REAL (dp), INTENT(OUT)  :: bnd

!     ***************************************************

!     THIS ROUTINE CALCULATES ERRORS USED IN TESTS IN STIFF .

!     ***************************************************
!     ..
!     .. LOCAL SCALARS ..
!     ..
REAL (dp) :: sqhol

sqhol = n*eps*eps
edn = tq(1)*tq(1)*sqhol

!     ** ERROR ASSOCIATED WITH LOWER ORDER METHOD

e = tq(2)*tq(2)*sqhol

!     ** ERROR ASSOCIATED WITH PRESENT ORDER

eup = tq(3)*tq(3)*sqhol

!     ** ERROR ASSOCIATED WITH HIGHER ORDER METHOD

bnd = tq(4)*tq(4)*sqhol*0.3_dp
RETURN

END SUBROUTINE errors


SUBROUTINE prdict(t, h, y, l, n, yprime, fcn)

IMPLICIT NONE
REAL (dp), INTENT(IN OUT) :: t
REAL (dp), INTENT(IN)     :: h
REAL (dp), INTENT(OUT)    :: y(:,:)    ! y(n,15)
INTEGER, INTENT(IN)       :: l
INTEGER, INTENT(IN)       :: n
REAL (dp), INTENT(IN OUT) :: yprime(:)

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

! **********************************************************************
!     PREDICTS A VALUE FOR Y AT (T+H) GIVEN THE HISTORY ARRAY AT Y(T)
!     THEN EVALUATES THE DERIVATIVE AT THIS POINT, THE RESULT OF THIS
!     EVALUATION BEING STORED IN YPRIME()
! **********************************************************************

!     ..
!     .. LOCAL SCALARS ..
INTEGER :: i,j2
!     ..

DO j2 = 2,l
  DO i = 1,n
    y(i,1) = y(i,1) + y(i,j2)
  END DO
END DO
t = t + h
CALL fcn(t, y(:,1), yprime)
nfe = nfe + 1
RETURN

END SUBROUTINE prdict


SUBROUTINE itrat2(y, n, t, hbeta, errbnd, arh, crate, tcrate, m, worked,  &
                  ymax, error, save1, save2, pw, ipiv, lmb, fcn)

IMPLICIT NONE
REAL (dp), INTENT(IN)      :: y(:,:)     ! y(n,15)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: t
REAL (dp), INTENT(IN)      :: hbeta
REAL (dp), INTENT(IN OUT)  :: errbnd
REAL (dp), INTENT(IN)      :: arh(:)
REAL (dp), INTENT(IN OUT)  :: crate
REAL (dp), INTENT(IN OUT)  :: tcrate
INTEGER, INTENT(IN OUT)    :: m
LOGICAL, INTENT(OUT)       :: worked
REAL (dp), INTENT(IN)      :: ymax(:)
REAL (dp), INTENT(OUT)     :: error(:)
REAL (dp), INTENT(IN OUT)  :: save1(:)
REAL (dp), INTENT(IN OUT)  :: save2(:)
REAL (dp), INTENT(IN OUT)  :: pw(:,:)
INTEGER, INTENT(IN OUT)    :: ipiv(:)
INTEGER, INTENT(IN)        :: lmb

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i
REAL (dp) :: d, d1
!     ..
!     .. DATA STATEMENTS ..

REAL (dp), PARAMETER :: zero = 0.0_dp
!     ..
IF(lmb == 1) GO TO 25

DO i = 1,n
  save1(i) = -save1(i) + hbeta*save2(i) + arh(i)
END DO
CALL sol(n, pw, save1, ipiv)
nbsol = nbsol + 1
d = zero
DO i = 1,n
  error(i) = error(i) + save1(i)
  d = d + (save1(i)/ymax(i))**2
  save1(i) = y(i,1) + error(i)
END DO
tcrate = tcrate + crate
d1 = d
m = 1
CALL fcn(t,save1,save2)
nfe = nfe + 1

25 worked = .true.
30 DO i = 1,n
  save1(i) = -save1(i) + hbeta*save2(i) + arh(i)
END DO

!     IF WE ARE HERE THEN PARTIALS ARE O.K.

CALL sol(n, pw, save1, ipiv)
nbsol = nbsol + 1

!     WE NOW CALCULATE A WEIGHTED RMS TYPE NORM

d = zero
DO i = 1,n
  error(i) = error(i) + save1(i)
  d = d + (save1(i)/ymax(i))**2
  save1(i) = y(i,1) + error(i)
END DO
! -------------------------------------------------------------------
!     TEST FOR CONVERGENCE.  IF M.GT.0 , AN ESTIMATE OF THE CONVERGENCE
!     RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
! -------------------------------------------------------------------
IF (m /= 0) THEN
  IF (d1 /= zero) crate = MAX(0.9_dp*crate, d/d1)
END IF

tcrate = tcrate + crate
IF ((d*MIN(1.0_dp, 2.0_dp*crate)) < errbnd) RETURN
IF (m /= 0) THEN
  IF (d > d1) THEN
    worked = .false.
    RETURN
  END IF
END IF

d1 = d
IF (m == 3) THEN
  worked = .false.
  RETURN
END IF

m = m + 1
CALL fcn(t, save1, save2)
nfe = nfe + 1
GO TO 30

END SUBROUTINE itrat2


SUBROUTINE stiff(eps, h, hmax, hmin, jstart, kflag, mf, t, tend, y, n, ymax, &
                 pw, pwcopy, yhold, ipiv, lout, fcn)

IMPLICIT NONE
REAL (dp), INTENT(IN)      :: eps
REAL (dp), INTENT(IN OUT)  :: h
REAL (dp), INTENT(IN)      :: hmax
REAL (dp), INTENT(IN)      :: hmin
INTEGER, INTENT(IN OUT)    :: jstart
INTEGER, INTENT(OUT)       :: kflag
INTEGER, INTENT(IN)        :: mf
REAL (dp), INTENT(IN OUT)  :: t
REAL (dp), INTENT(IN)      :: tend
REAL (dp), INTENT(OUT)     :: y(:,:)        ! y(n,15)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: ymax(:)
REAL (dp), INTENT(IN OUT)  :: pw(:,:)
REAL (dp), INTENT(IN OUT)  :: pwcopy(:)
REAL (dp), INTENT(IN OUT)  :: yhold(:,:)    ! yhold(n,15)
INTEGER, INTENT(IN OUT)    :: ipiv(:)
INTEGER, INTENT(IN)        :: lout

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

!     ------------------------------------------------------------------
!     STIFF PERFORMS ONE STEP OF THE INTEGRATION OF AN INITIAL VALUE
!     PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS.
!     COMMUNICATION WITH STIFF IS DONE WITH THE FOLLOWING VARIABLES..
!     Y      AN N BY LMAX+4 ARRAY CONTAINING THE DEPENDENT VARIABLES
!              AND THEIR BACKWARD DIFFERENCES.  LMAX-1 =MAXDER IS THE
!              MAXIMUM ORDER AVAILABLE.  SEE SUBROUTINE COSET.
!              Y(I,J+1) CONTAINS THE JTH BACKWARD DIFFERENCE OF Y(I)
!     T      THE INDEPENDENT VARIABLE. T IS UPDATED ON EACH STEP TAKEN.
!     H      THE STEPSIZE TO BE ATTEMPTED ON THE NEXT STEP.
!              H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING
!              THE PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE BUT
!              ITS SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.
!     HMIN   THE MINIMUM AND MAXIMUM ABSOLUTE VALUE OF THE STEPSIZE
!     HMAX   TO BE USED FOR THE STEP.  THESE MAY BE CHANGED AT ANY
!              TIME BUT WILL NOT TAKE EFFECT UNTIL THE NEXT H CHANGE.
!     EPS    THE RELATIVE ERROR BOUND.  SEE DESCRIPTION IN MEBDF
!     UROUND THE UNIT ROUNDOFF OF THE MACHINE.
!     N      THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
!     MF     THE METHOD FLAG.  MUST BE SET TO 21 OR 22 AT PRESENT
!     KFLAG  A COMPLETION CODE WITH THE FOLLOWING MEANINGS..
!                  0  THE STEP WAS SUCCESSFUL
!                 -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED
!                       WITH ABS(H) = HMIN.
!                 -2  THE REQUESTED ERROR IS SMALLER THAN CAN
!                       BE HANDLED FOR THIS PROBLEM.
!                 -3  CORRECTOR CONVERGENCE COULD NOT BE
!                       ACHIEVED FOR ABS(H)=HMIN.
!        ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF T AND THE Y ARRAY ARE
!        AS AT THE BEGINNING OF THE LAST STEP ATTEMPTED, AND H IS THE LAST
!        STEP SIZE ATTEMPTED.
!     JSTART  AN INTEGER USED ON INPUT AND OUTPUT.
!       ON INPUT IT HAS THE FOLLOWING VALUES AND MEANINGS..
!              0  PERFORM THE FIRST STEP.
!          .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST
!          .LT.0  TAKE THE NEXT STEP WITH A NEW VALUE OF H, EPS OR N
!       ON EXIT, JSTART IS NQUSED, THE ORDER OF THE METHOD LAST USED.
!     YMAX     AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL
!              ERRORS IN Y ARE COMPARED
!     PW       A BLOCK OF LOCATIONS USED FOR PARTIAL DERIVATIVES
!     IPIV     AN INTEGER ARRAY OF LENGTH N USED FOR PIVOT INFORMATION.

!     JNEWIM   IS TO INDICATE IF PRESENT ITERATION MATRIX
!                WAS FORMED USING A NEW J OR OLD J.
!     JSNOLD   KEEPS TRACK OF NUMBER OF STEPS TAKEN WITH
!                PRESENT ITERATION MATRIX (BE IT FORMED BY
!                A NEW J OR NOT).
!     AVNEWJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
!                MATRIX WAS FORMED BY A NEW J.
!     AVOLDJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
!                MATRIX WAS FORMED BY AN OLD J.
!     NRENEW   FLAG THAT IS USED IN COMMUNICATION WITH SUBROUTINE PSET.
!                IF  NRENEW > 0  THEN FORM A NEW JACOBIAN BEFORE COMPUTING THE
!                                COEFFICIENT MATRIX FOR THE NEWTON ITERATION
!                           = 0  FORM THE COEFFICIENT MATRIX USING A
!                                COPY OF AN OLD JACOBIAN
!     NEWPAR   FLAG USED IN THIS SUBROUTINE TO INDICATE IF A JACOBIAN
!              HAS BEEN EVALUATED FOR THE CURRENT STEP
! **********************************************************************

!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i, ier, iiter, iiter2, ijus, itst, j, j1, j2, kfail, ll,  &
             m3step, maxder, meth, miter, newq
LOGICAL   :: finish, worked
REAL (dp) :: atol, d, ddown, dstep2, dup, ebdf1, ebdf2, efail, enq1, enq2,  &
             enq3, prbdf1, prbdf2, prfail, pr1, pr2, pr3, qq, red, rhbdf1,  &
             rhbdf2, rhlmit, trange
!     ..
!     .. LOCAL ARRAYS ..
REAL (dp), SAVE  :: el(10), elst(10), tq(4)
REAL (dp)        :: error(n), save1(n), save2(n), ynhold(n,2)
! If N > 50, increase the 50 in the line below
REAL (dp), SAVE  :: arh(50)
!     ..
!     .. SAVE STATEMENT ..
REAL (dp), SAVE :: epsold, hold, edn, eup, bnd, e
INTEGER, SAVE   :: mfold, lmax
!     ..
!     .. DATA STATEMENTS ..

REAL (dp), SAVE      :: oldlo = 1.0_dp
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp
!     ..

IF (jstart == 0) arh = 0.0_dp
el(2) = one
elst(2) = one
told = t
kflag = 0
IF (jstart > 0) GO TO 60
IF (jstart /= 0) GO TO 30
!     ------------------------------------------------------------------
!     ON THE FIRST CALL, THE ORDER IS SET TO 1 AND THE INITIAL YDOT
!     IS CALCULATED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE
!     INCREASED IN A SINGLE STEP.  RMAX IS SET EQUAL TO 1.D4 INITIALLY
!     TO COMPENSATE FOR THE SMALL INITIAL H, BUT THEN IS NORMALLY = 10.
!     IF A FAILURE OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST),
!     RMAX IS SET AT 2. FOR THE NEXT INCREASE.
!     ------------------------------------------------------------------
CALL fcn(t, y(:,1), save1)
y(1:n,2) = h*save1(1:n)
meth = 2
miter = mf - 10*meth
nq = 1
l = 2
idoub = 3
kfail = 0
rmax = 10000.0_dp
rc = zero
crate1 = one
crate2 = one
jsnold = 0
jnewim = .true.
tcrat1 = zero
tcrat2 = zero
top=0
atol=eps/10.0_dp
DO i=1,12
  hstpsz(1,i)=1.0_dp
  hstpsz(2,i)=atol
END DO
hold = h
mfold = mf
nstep = 0
nfe = 1
nje = 0
ndec = 0
npset = 0
ncoset = 0
maxord = 1
nfail1 = 0
nfail2 = 0
cfail = .true.
avnewj = zero
avoldj = zero
avnew2 = zero
avold2 = zero
sample = .false.
isamp = 0
!     **************************************************
!     CFAIL=.TRUE. ENSURES THAT WE CALCULATE A NEW
!     J ON THE FIRST CALL.
!     **************************************************
meqc1 = 0
meqc2 = 0
mq1tmp = 0
mq2tmp = 0
nbsol = 0
!     -----------------------------------------------------------------
!     IF THE CALLER HAS CHANGED N OR EPS, THE CONSTANTS E, EDN, EUP
!     AND BND MUST BE RESET.  E IS A COMPARISON FOR ERRORS AT THE
!     CURRENT ORDER NQ.  EUP IS TO TEST FOR INCREASING THE ORDER,
!     EDN FOR DECREASING THE ORDER.  BND IS USED TO TEST FOR CONVERGENCE
!     OF THE CORRECTOR ITERATES.   IF THE CALLER HAS CHANGED H, Y MUST
!     BE RE-SCALED.  IF H IS CHANGED, IDOUB IS SET TO L+1 TO PREVENT
!     FURTHER CHANGES IN H FOR THAT MANY STEPS.
!     -----------------------------------------------------------------
CALL coset(nq, el, elst, tq, maxder)
lmax = maxder + 1
rc = rc*el(1)/oldlo
oldlo = el(1)
iweval = miter
nrenew = 1
newpar = 0
!     *****************************************************
!     NRENEW AND NEWPAR ARE TO INSTRUCT ROUTINE THAT
!     WE WISH A NEW J TO BE CALCULATED FOR THIS STEP.
!     *****************************************************
CALL errors(n, tq, eps, edn, e, eup, bnd)
epsold = eps
arh(1:n) = el(2)*y(1:n,1)

yhold(1:n,1:l) = y(1:n,1:l)
qi = h*el(1)
qq = one/qi
CALL prdict(t, h, y, l, n, save2, fcn)
GO TO 110
!     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DIFFERENT PARAMETERS ON THIS CALL        <
!     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
30 y(1:n,1:l) = yhold(1:n,1:l)
IF (mf /= mfold) THEN
  meth = mf/10
  miter = mf - 10*meth
  mfold = mf
  iweval = miter
END IF

IF (eps /= epsold) THEN
  CALL errors(n, tq, eps, edn, e, eup, bnd)
  epsold = eps
END IF

IF (h /= hold) THEN
  rh = h/hold
  h = hold
  GO TO 50
  
ELSE
  GO TO 60
  
END IF

!     *********************************************
!     RE-SCALE Y AFTER A CHANGE OF STEPSIZE   *
!     *********************************************
40 rh = MAX(rh,hmin/ABS(h))
50 rh = MIN(rh,hmax/ABS(h),rmax)
CALL rscale(n, l, rh, y)
rmax = 10.0_dp
jchang = 1
h = h*rh
rc = rc*rh
IF (jsnold > 40) THEN
  cfail = .true.
  newpar = 0
  rc = zero
! **********************************************************************
!        CFAIL=TRUE AND NEWPAR=0 SHOULD FORCE A NEW J TO BE EVALUATED
!        AFTER 42 STEPS WITH AN OLD J, IF WE HAVE HAD A FAILURE OF ANY
!        KIND ON THE FIRST, SECOND OR THIRD STAGE OF THE CURRENT STEP
! **********************************************************************
END IF

idoub = l + 1
yhold(1:n,1:l) = y(1:n,1:l)

60 IF (ABS(rc-one) > 0.3_dp) iweval = miter
!     ------------------------------------------------------------------
!     THIS SECTION COMPUTES THE PREDICTED VALUES OF Y
!     AND THE RHS, ARH, FOR USE IN THE NEWTON ITERATION SCHEME.
!     RC IS THE RATIO OF THE NEW TO OLD VALUES OF THE COEFFICIENT H*EL(1).
!     WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IWEVAL IS
!     SET TO MITER TO FORCE THE PARTIALS TO BE UPDATED.
!     ------------------------------------------------------------------
qi = h*el(1)
qq = one/qi
arh(1:n) = el(2)*y(1:n,1)
DO j1 = 2,nq
  DO i = 1,n
    arh(i) = arh(i) + el(j1+1)*y(i,j1)
  END DO
END DO
IF (jchang == 1) THEN
!        IF WE HAVE CHANGED STEPIZE THEN PREDICT A VALUE FOR Y(T+H)
!        AND EVALUATE THE DERIVATIVE THERE (STORED IN SAVE2())
  CALL prdict(t, h, y, l, n, save2, fcn)
  
ELSE
!        ELSE USE THE VALUES COMPUTED FOR THE SECOND BDF FROM THE LAST STEP.
!        Y( ,LMAX+3) HOLDS THE VALUE FOR THE DERIVATIVE AT (T+H)
!        AND Y( ,LMAX+4) HOLDS THE APPROXIMATION TO Y AT THIS POINT.
  DO i = 1,n
    save2(i) = y(i,lmax+3)
    y(i,1) = y(i,lmax+4)
  END DO
  t = t + h
END IF

110 IF (iweval <= 0) GO TO 120
! -------------------------------------------------------------------
!     IF INDICATED, THE MATRIX P = I - H*EL(1)*J IS RE-EVALUATED BEFORE
!     STARTING THE CORRECTOR ITERATION.  IWEVAL IS SET = 0 TO INDICATE
!     THAT THIS HAS BEEN DONE.  P IS COMPUTED AND PROCESSED IN PSET.
!     THE PROCESSED MATRIX IS STORED IN PW
! -------------------------------------------------------------------
iweval = 0
rc = one
iiter = meqc1 - mq1tmp
iiter2 = meqc2 - mq2tmp
IF (jnewim) THEN
  IF (jsnold >= 3) THEN
    avnewj = tcrat1 / iiter
    avnew2 = tcrat2 / iiter2
    
  ELSE
    avnewj = one
    avnew2 = one
  END IF
  
ELSE
  
!          MATRIX P WAS FORMED WITH A COPY OF J
  
  IF (jsnold >= 3) THEN
    avoldj = tcrat1 / iiter
    avold2 = tcrat2 / iiter2
    IF (avoldj < avnewj) THEN
      avnewj = avoldj
      
    ELSE IF (ABS(avoldj-avnewj) > 0.3_dp .OR.  &
            (avoldj > 0.85_dp .AND. avoldj /= one)) THEN
      
!              SINCE IN CERTAIN INSTANCES AVOLDJ WILL
!              BE 1.0 AND THERE WILL BE NO NEED TO UPDATE J.
      
      cfail = .true.
      crate1 = one
      crate2 = one
    END IF
    
  ELSE
    cfail = .true.
    crate1 = one
    crate2 = one
    
!           *********************************************************
!           IF WE HAVE REACHED HERE THINGS MUST HAVE GONE WRONG
!           *********************************************************
    
  END IF
  
END IF

tcrat1 = zero
tcrat2 = zero
IF (cfail) THEN
  IF (newpar == 1) THEN
    nrenew = 0
    jnewim = .true.
    
  ELSE
    nrenew = 1
    newpar = 1
    jsinup = -1
    jnewim = .true.
  END IF
  
ELSE
  nrenew = 0
  jnewim = .false.
END IF

cfail = .false.
jsnold = 0
mq1tmp = meqc1
mq2tmp = meqc2
CALL pset(y, n, h, t, uround, epsjac, -qi, miter, ier, nrenew, ymax, save2,  &
          pw, pwcopy, error, ipiv, fcn)
!     NOTE THAT ERROR() IS JUST BEING USED AS A WORKSPACE BY PSET
IF (ier /= 0) THEN
!     IF IER>0 THEN WE HAVE  A  SINGULARITY IN THE COEFFICIENT MATRIX
  ijus=1
  red=0.5_dp
  GO TO 450
  
END IF


120 DO i = 1,n
  save1(i) = y(i,1)
  error(i) = zero
END DO
m1 = 0
! **********************************************************************
!     UP TO 4 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS
!     MADE ON THE R.M.S. NORM OF EACH CORRECTION ,USING BND, WHICH
!     DEPENDS ON EPS.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
!     VECTOR  ERROR(I).  THE Y ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.
!     THE UPDATED Y VECTOR IS STORED TEMPORARILY IN SAVE1.
! **********************************************************************
IF (.NOT.sample) THEN
  CALL itrat2(y, n, t, qi, bnd, arh, crate1, tcrat1, m1, worked, ymax,  &
              error, save1, save2, pw, ipiv, 1, fcn)
  itst = 2
  
ELSE
  CALL itrat2(y, n, t, qi, bnd, arh, crate1, tcrat1, m1, worked, ymax,  &
              error, save1, save2, pw, ipiv, 0, fcn)
  itst = 3
END IF

meqc1 = meqc1 + m1 + 1

!       NOW TEST TO SEE IF IT WAS SUCCESSFUL OR NOT

IF (.NOT.worked) THEN
  nfail1 = nfail1 + 1
! **********************************************************************
!        THE CORRECTOR ITERATION FAILED TO CONVERGE IN 4 TRIES. IF
!        PARTIALS ARE NOT UP TO DATE, THEY ARE RE-EVALUATED FOR THE
!        NEXT TRY. OTHERWISE THE Y ARRAY IS RETRACTED TO ITS VALUES
!        BEFORE PREDICTION AND H IS REDUCED IF POSSIBLE. IF NOT A
!        NON-CONVERGENCE EXIT IS TAKEN
! **********************************************************************
  IF (iweval == -1) THEN
!           HAVE BEEN USING OLD PARTIALS, UPDATE THEM AND TRY AGAIN
    iweval = miter
    cfail = .true.
    CALL fcn(t, y(:,1), save2)
    nfe = nfe + 1
    GO TO 110
    
  END IF
  
  ijus = 0
  red = 0.5_dp
  GO TO 450
  
END IF

iweval = -1
hused = h
nqused = nq
DO i = 1,n
  save2(i) = (save1(i) - arh(i))*qq
  y(i,1) = save1(i)
END DO

!     UPDATE THE DIFFERENCES AT N+1

DO j = 2,l
  DO i = 1,n
    y(i,j) = y(i,j-1) - yhold(i,j-1)
  END DO
END DO

!     COMPUTE ERROR IN THE SOLUTION

d = zero
DO i = 1,n
  d = d + ((y(i,l) - yhold(i,l))/ymax(i))**2
END DO

!     STORING Y FROM FIRST STEP FOR USE IN THIRD STEP

DO i = 1,n
  ynhold(i,1) = y(i,1)
  ynhold(i,2) = save2(i)
END DO
IF (d > e) GO TO 330
IF ((m1+m2) >= itst) THEN
  m2 = 0
  ebdf1 = 0.5_dp / l
  prbdf1 = ((d/ (e*0.7_dp))**ebdf1)*1.5_dp + 1.6D-6
  rhbdf1 = one/prbdf1
  CALL hchose(rhbdf1, h)
  rhlmit = 1 - (5+itst-m1-m2)*0.5D-1
  IF (rhbdf1 < rhlmit) THEN
    ijus = 1
    red = rhbdf1
    GO TO 450
  END IF
  
END IF

kfail = 0
! ----------------------------------------------------------------------
DO i = 1,n
  arh(i) = el(2)*y(i,1)
END DO
DO j1 = 2,nq
  DO i = 1,n
    arh(i) = arh(i) + el(j1+1)*y(i,j1)
  END DO
END DO
CALL prdict(t,h,y,l,n,save2,fcn)
DO i = 1,n
  save1(i) = y(i,1)
  error(i) = zero
END DO
m2 = 0

!     FOR NOW WILL ASSUME THAT WE DO NOT WISH TO SAMPLE
!     AT THE N+2 STEP POINT

CALL itrat2(y, n, t, qi, bnd, arh, crate2, tcrat2, m2, worked, ymax, error,  &
            save1, save2, pw, ipiv, 1, fcn)
meqc2 = meqc2 + m2 + 1

!       NOW CHECK TO SEE IF IT WAS SUCCESSFUL OR NOT

IF (.NOT.worked) THEN
  nfail2 = nfail2 + 1
  ijus = 0
  red = 0.5_dp
  GO TO 450
  
END IF

!        IF WE ARE DOWN TO HERE THEN THINGS MUST HAVE CONVERGED

DO i = 1,n
  y(i,lmax+3) = (save1(i)-arh(i))*qq
  y(i,lmax+4) = save1(i)
END DO
IF ((m2 >= 2) .AND. (crate2 < 0.35_dp)) THEN
  DO j = 2,nq
    DO i = 1,n
      save1(i) = save1(i) - y(i,j)
    END DO
  END DO
  dstep2 = zero
  DO i = 1,n
    dstep2 = dstep2 + ((save1(i)-ynhold(i,1)-y(i,l))/ymax(i))**2
  END DO
  ebdf2 = 0.5_dp / l
  prbdf2 = ((dstep2/ (e*0.7_dp))**ebdf2)*1.5_dp + 1.6D-6
  rhbdf2 = one/prbdf2
  CALL hchose(rhbdf2, h)
  IF (rhbdf2 < 0.9_dp) THEN
    ijus = 1
    red = rhbdf2
    GO TO 450
  END IF
  
END IF


!     WE ARE NOW COMPUTING THE THIRD STAGE

ll = l + 1
t = told + h
DO i = 1,n
  arh(i) = h* (elst(nq+2)*y(i,lmax+3) + (elst(1)-el(1))*ynhold(i,2))
  DO j1 = 1,nq
    arh(i) = arh(i) + elst(j1+1)*yhold(i,j1)
  END DO
END DO
DO i = 1,n
  save2(i) = ynhold(i,2)
  y(i,1) = ynhold(i,1)
END DO
m3step = 0
300 DO i = 1,n
  save1(i) = -y(i,1) + qi*save2(i) + arh(i)
END DO
CALL sol(n, pw, save1, ipiv)
nbsol = nbsol + 1
d = zero
DO i = 1,n
  d = d + (save1(i)/ymax(i))**2
  y(i,1) = y(i,1) + save1(i)
END DO
IF ((d*MIN(one, 2.0_dp*crate1)) <= bnd) GO TO 360
IF (m3step == 4) THEN
  WRITE (lout,9000)
  ijus=1
  red=0.5_dp
  GO TO 450
END IF

m3step = m3step + 1
CALL fcn(t, y(:,1), save2)
nfe = nfe + 1
GO TO 300

330 kfail = kfail - 1
! **************************************************************************
!   THE ERROR TEST FAILED. KFAIL KEEPS TRACK OF MULTIPLE FAILURES.  RESTORE T
!   AND THE Y ARRAY TO THEIR PREVIOUS VALUES AND PREPARE TO TRY THE STEP
!   AGAIN.  COMPUTE THE OPTIMAL STEP SIZE FOR THIS ORDER AND ONE ORDER LOWER.
! **************************************************************************
t = told
hold = h
efail = 0.5_dp / l
y(1:n,1:l) = yhold(1:n,1:l)
rmax = 2.0_dp
IF (ABS(h) <= hmin*1.00001_dp) THEN
  
!        REQUESTED ERROR NOT POSSIBLE WITH GIVEN HMIN
  
  kflag = -1
  hold = h
  RETURN
  
END IF

IF (kfail <= -3) GO TO 340

!     PREDICTING A NEW H AFTER INSUFFICIENT ACCURACY

prfail = ((d/e)**efail)*1.5_dp + 1.6D-6
newq = nq
rh = one / (prfail*DBLE(-kfail))
CALL hchose(rh, h)
GO TO 40
! **********************************************************************
!     CONTROL REACHES THIS STAGE IF 3 OR MORE FAILURES HAVE OCCURED.
!     IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE Y
!     ARRAY HAVE ERRORS OF THE WRONG ORDER. HENCE THE FIRST DERIVATIVE
!     IS RE-COMPUTED, AND THE ORDER IS SET TO 1. THEN H IS REDUCED BY A
!     FACTOR OF 10, AND THE STEP IS RETRIED. AFTER A TOTAL OF 7
!     FAILURES AN EXIT IS TAKEN WITH KFLAG=-2.
! **********************************************************************
340 IF (kfail == -7) THEN
!        ERROR SMALLER THAN CAN BE HANDLED FOR PROBLEM
  kflag = -2
  hold = h
  RETURN
  
END IF
!     *********************************
!     START FROM ORDER 1 AGAIN    *
!     *********************************
jchang = 1
rh = MAX(hmin/ABS(h),0.1_dp)
CALL hchose(rh, h)
h = h*rh
CALL fcn(t, yhold(:,1), save1)
nfe = nfe + 1
DO i = 1,n
  y(i,1) = yhold(i,1)
  y(i,2) = h*save1(i)
  yhold(i,2) = y(i,2)
END DO
iweval = miter
cfail = .true.
!     SINCE WE HAVE HAD PROBLEMS PROCEED WITH THIS ORDER
!     FOR 10 STEPS (IF WE CAN)
idoub = 10
IF (nq == 1) GO TO 60
nq = 1
l = 2
!     RESET ORDER, RECALCULATE ERROR BOUNDS
CALL coset(nq,el,elst,tq,maxder)
lmax = maxder + 1
rc = rc*el(1)/oldlo
oldlo = el(1)
CALL errors(n, tq, eps, edn, e, eup, bnd)
!     NOW JUMP TO NORMAL CONTINUATION POINT
GO TO 60

! **********************************************************************
!     THE ITERATION FOR THE CORRECTED SOLUTION HAS CONVERGED.
!     UPDATE THE Y ARRAY.
! **********************************************************************
360 DO j2 = 2,ll
  DO i = 1,n
    y(i,j2) = y(i,j2-1) - yhold(i,j2-1)
  END DO
END DO
! ----------------------------------------------------------------------
!     IF THE COUNTER IDOUB EQUALS 2 , STORE Y(I,LMAX+5) WHICH IS USED IN
!     ASSESSING THE POSSIBILITY OF INCREASING THE ORDER. IF IDOUB = 0
!     CONTROL PASSES TO 480 WHERE AN ATTEMPT TO CHANGE THE STEPSIZE AND
!     ORDER IS MADE.
! ----------------------------------------------------------------------
IF (idoub == 2) THEN
  DO i = 1,n
    y(i,lmax+5) = y(i,ll)
  END DO
END IF

idoub = idoub - 1
trange=(tend-told-h)*h
IF(trange < 0.0_dp) GO TO 440
jchang = 0
IF (idoub == 0) THEN
  sample = .false.
  isamp = isamp + 1
  IF (isamp == 4) THEN
    sample = .true.
    isamp = 0
  END IF
  
! **********************************************************************
!        NOW COMPUTE THE FACTORS PR1, PR2 AND PR3, BY WHICH
!        H COULD BE DIVIDED AT ORDER NQ-1, ORDER NQ AND ORDER NQ+1
!        RESPECTIVELY. THE SMALLEST OF THESE IS DETERMINED AND THE NEW
!        ORDER CHOSEN ACCORDINGLY. IF THE ORDER IS TO BE INCREASED WE
!        MUST COMPUTE ONE MORE BACKWARD DIFFERENCE.
! **********************************************************************
  pr3 = 1.d+20
  IF (l /= lmax) THEN
    dup = zero
    DO i = 1,n
      dup = dup + ((y(i,ll) - y(i,lmax+5)) / ymax(i))**2
    END DO
    enq3 = 0.5_dp / (l+1)
    pr3 = ((dup/eup)**enq3)*1.7_dp + 1.8D-6
  END IF
  
  enq2 = 0.5_dp / l
  d = zero
  DO i = 1,n
    d = d + ((y(i,l) - yhold(i,l)) / ymax(i))**2
  END DO
  pr2 = ((d/e)**enq2)*1.5_dp + 1.6D-6
  pr1 = 1.d+20
  IF (nq > 1) THEN
    ddown = zero
    DO i = 1,n
      ddown = ddown + (y(i,l)/ymax(i))**2
    END DO
    enq1 = 0.5_dp / nq
    pr1 = ((ddown/edn)**enq1)*1.6_dp + 1.7D-6
  END IF
  
  IF (pr2 <= pr3) THEN
    IF (pr2 > pr1) THEN
      newq = nq - 1
      rh = 1.0_dp/pr1
      
    ELSE
      newq = nq
      rh = 1.0_dp/pr2
    END IF
    
  ELSE IF (pr3 < pr1) THEN
    newq = l
    rh = 1.0_dp/pr3
    
  ELSE
    newq = nq - 1
    rh = 1.0_dp/pr1
  END IF
  
  rh = MIN(rh,rmax)
  CALL hchose(rh, h)
  IF ((jsinup <= 20).AND.(kflag == 0).AND.(rh < 1.1_dp)) THEN
!           WE HAVE RUN INTO PROBLEMS
    idoub = 10
    nq = nqused
    l = nq + 1
    GO TO 440
    
  END IF
! **********************************************************************
!        IF THERE IS A CHANGE IN ORDER, RESET NQ, L AND THE
!        COEFFICIENTS. IN ANY CASE H IS RESET ACCORDINGLY TO RH AND THE
!        Y ARRAY IS RE-SCALED
! **********************************************************************
  IF (newq /= nq) THEN
    IF (newq > nq) THEN
!              ADD AN EXTRA TERM TO THE HISTORY ARRAY
      DO i = 1,n
        y(i,ll) = y(i,l) - yhold(i,l)
      END DO
    END IF
    
    nq = newq
    l = nq + 1
!           RESET ORDER,RECALCULATE ERROR BOUNDS
    CALL coset(nq,el,elst,tq,maxder)
    lmax = maxder + 1
    rc = rc*el(1)/oldlo
    oldlo = el(1)
    CALL errors(n, tq, eps, edn, e, eup, bnd)
  END IF
  
!        NOW RESCALE THE HISTORY ARRAY FOR THE NEW STEPSIZE
  rh = MAX(rh,hmin/ABS(h))
  rh = MIN(rh,hmax/ABS(h),rmax)
  CALL rscale(n,l,rh,y)
  rmax = 10.0_dp
  jchang = 1
  h = h*rh
  rc = rc*rh
  IF (jsnold > 40) THEN
    rc = zero
!           IF WE HAVE BEEN USING THE SAME COEFFICIENT MATRIX FOR MORE
!           THAN 40 STEPS FORCE A NEW ONE TO BE FORMED ON THE NEXT STEP.
!           CODE WILL ATTEMPT TO USE A COPY OF THE LAST JACOBIAN FORMED
  END IF
  
  idoub = l + 1
  
END IF

! ----------------------------------------------------------------------
!     STORE THE Y ARRAY IN THE MATRIX YHOLD.  STORE IN THE Y ARRAY THE
!     INFORMATION NECESSARY TO PERFORM AN INTERPOLATION TO FIND THE
!     SOLUTION AT THE SPECIFIED OUTPUT POINT IF APPROPRIATE.
! ----------------------------------------------------------------------
440 yhold(1:n,1:l) = y(1:n,1:l)
nstep = nstep + 1
jsinup = jsinup + 1
jsnold = jsnold + 1
jstart = nqused
t = told + hused
hold = h
kfail = 0
newpar = 0
cfail = .false.
RETURN

450 finish = .false.
t = told
rmax = 2.0_dp
DO j1=1,l
  DO i=1,n
    y(i,j1) = yhold(i,j1)
  END DO
END DO
IF(ABS(h) <= hmin*1.00001_dp) THEN
  
!   CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED
  
  IF(nstep == 0) THEN
    kflag = -1
  ELSE
    kflag = -3
  END IF
  
!    TO SUPPRESS ERROR MESSAGES AT START AS H MAY
!    HAVE BEEN TOO LARGE ON THE FIRST STEP.
  
  hold = h
  finish = .true.
END IF
rh = red

!     TRY AGAIN WITH UPDATED PARTIALS

IF(ijus == 0) CALL hchose(rh, h)
IF(.NOT.finish) THEN
  GO TO 40
END IF

RETURN
! ------------------- END OF SUBROUTINE STIFF --------------------------
9000 FORMAT ('  CORRECTOR HAS NOT CONVERGED')
END SUBROUTINE stiff


SUBROUTINE rscale(n, l, rh, y)

IMPLICIT NONE
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: l
REAL (dp), INTENT(IN)      :: rh
REAL (dp), INTENT(IN OUT)  :: y(:,:)      ! y(n,15)

!     ..
!     .. LOCAL SCALARS ..
INTEGER   :: i, j
REAL (dp) :: ta, tb, tc, td, te, tf, zz
!     ..
!     .. LOCAL ARRAYS ..
REAL (dp) :: di(8,8)

!    SUBROUTINE IS FOR RESCALING THE HISTORY ARRAY AFTER A CHANGE IN STEPSIZE

!    N      ORDER OF THE PROBLEM
!    L      NUMBER OF TERMS IN THE HISTORY ARRAY TO BE RESCALED
!    RH     RATIO OF THE STEPSIZE CHANGE (I.E. RH = HNEW/HOLD)
!    Y()    THE HISTORY ARRAY

!     ..

di(2,2) = rh
IF (l > 2) THEN
  ta = rh*rh
  di(2,3) = rh*(1.0_dp - rh)/2.0_dp
  di(3,3) = ta
  IF (l > 3) THEN
    tb = ta*rh
    di(2,4) = rh*((rh - 3.0_dp)*rh + 2.0_dp)/6.0_dp
    di(3,4) = ta*(1.0_dp - rh)
    di(4,4) = tb
    IF (l > 4) THEN
      tc = tb*rh
      di(2,5) = - (((rh - 6.0_dp)*rh + 11.0_dp)*rh - 6.0_dp)*rh/ 24.0_dp
      di(3,5) = ta*((7.0_dp*rh - 18.0_dp)*rh + 11.0_dp)/12.0_dp
      di(4,5) = 1.5_dp*tb* (1.0_dp - rh)
      di(5,5) = tc
      IF (l > 5) THEN
        td = tc*rh
        di(2,6) = ((((rh - 10.0_dp)*rh + 35.0_dp)*rh - 50.0_dp)  &
                  *rh + 24.0_dp)*rh / 120.0_dp
        di(3,6) = - (((3.0_dp*rh - 14.0_dp)*rh + 21.0_dp)*rh - 10.0_dp)*ta/12.0_dp
        di(4,6) = ((5.0_dp*rh - 12.0_dp)*rh + 7.0_dp)*tb/4.0_dp
        di(5,6) = 2.0_dp*tc*(1.0_dp - rh)
        di(6,6) = td
        IF (l > 6) THEN
          te = td*rh
          di(2,7) = -rh*(rh - 1.0_dp)*(rh - 2.0_dp)*  &
                    (rh - 3.0_dp)*(rh - 4.0_dp)*(rh - 5.0_dp)/ 720.0_dp
          di(3,7) = ta* ((((62.0_dp*rh - 450.0_dp)*rh +  &
                    1190.0_dp)*rh - 1350.0_dp)*rh + 548.0_dp) /720.0_dp
          di(4,7) = tb* (((-18.0_dp*rh + 75.0_dp)*rh  &
                    - 102.0_dp)*rh + 45.0_dp)/ 24.0_dp
          di(5,7) = tc* ((13.0_dp*rh - 30.0_dp)*rh + 17.0_dp) /6.0_dp
          di(6,7) = 2.5_dp*td*(1.0_dp - rh)
          di(7,7) = te
          IF (l > 7) THEN
            tf = te*rh
            di(2,8) = rh*(rh - 1.0_dp)*(rh - 2.0_dp)*(rh - 3.0_dp)*  &
                      (rh - 4.0_dp)*(rh - 5.0_dp)*(rh - 6.0_dp)/ 5040.0_dp
            di(3,8) = ta*(((((-126.0_dp*rh + 1302.0_dp)*rh - 5250.0_dp)*  &
                      rh + 10290.0_dp)*rh - 9744.0_dp )*rh + 3528.0_dp)/5040.0_dp
            di(4,8) = tb*((((43.0_dp*rh - 270.0_dp)*rh +  &
                      625.0_dp)*rh - 630.0_dp)*rh + 232.0_dp) /120.0_dp
            di(5,8) = tc*(((-10.0_dp*rh + 39.0_dp)*rh -  &
                      50.0_dp)*rh + 21.0_dp)/ 6.0_dp
            di(6,8) = td* ((20.0_dp*rh - 45.0_dp)*rh + 25.0_dp )/ 6.0_dp
            di(7,8) = 3.0_dp*te*(1.0_dp - rh)
            di(8,8) = tf
          END IF
          
        END IF
        
      END IF
      
    END IF
    
  END IF
  
END IF

DO i = 1,n
  DO j = 2,l
    zz = DOT_PRODUCT( di(j,j:l), y(i,j:l) )
    y(i,j) = zz
  END DO
END DO
RETURN

END SUBROUTINE rscale


SUBROUTINE hchose(rh, h)

IMPLICIT NONE
REAL (dp), INTENT(IN OUT) :: rh
REAL (dp), INTENT(IN)     :: h

! Local variables
INTEGER :: i, i2

!     FIRST MOVE ALL ELEMENTS DOWN ONE PLACE

IF (h /= hstpsz(2,1)) THEN
  DO i=12,2,-1
    i2 = i - 1
    hstpsz(1,i) = hstpsz(1,i2)
    hstpsz(2,i) = hstpsz(2,i2)
  END DO
  
!          NOW INSERT VALUE OF H USED BEFORE THIS CALL
  
  hstpsz(1,2) = h/hstpsz(2,1)
  hstpsz(2,1) = h
END IF

!     NOW DECIDE ON THE NEW CHANGE

IF (rh > 1.0_dp) THEN
ELSE IF (hstpsz(1,2) <= 1.0_dp) THEN
ELSE IF ((rh*h) <= hstpsz(2,2)) THEN
ELSE
  rh = hstpsz(2,2)/h
END IF
hstpsz(1,1) = rh
RETURN

END SUBROUTINE hchose

END MODULE toms703
