FUNCTION BVALU ( T, A, N, K, IDERIV, X, INBV, WORK )

!*****************************************************************************80
!
!! BVALU evaluates a B-spline or its derivatives.
!
!  Discussion:
!
!    BVALU evaluates the B-representation (T,A,N,K) of a B-spline
!    at X for the function value on IDERIV=0 or any of its
!    derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!    (right derivatives) are returned except at the right end
!    point X=T(N+1) where left limiting values are computed.  The
!    spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
!    a fatal error message when X is outside of this interval.
!
!    To compute left derivatives or left limiting values at a
!    knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Daniel Amos
!
!  Reference:
!
!    Carl de Boor, 
!    Package for calculating with B-splines,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 3, June 1977, pages 441-472.
!
!  Parameters:
!
!         Input      T,A,X are real
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV = 0 returns the B-spline value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time BVALU is called.
!
!         Output     WORK,BVALU are real
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          BVALU  - value of the IDERIV-th derivative at X
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  real ( kind = 4 ) bvalu
  integer ( kind = 4 ) I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ
  integer ( kind = 4 ) IP1, IP1MJ, J, JJ, J1, J2, KMIDER, KMJ, KM1, KPK, MFLAG
  real ( kind = 4 ) A(*), FKMJ, T(N+K), WORK(3*K), X

  BVALU = 0.0E0

  IF(K.LT.1) then
    CALL XERMSG ('SLATEC', 'BVALU', 'K DOES NOT SATISFY K.GE.1', 2,1)
    RETURN
  end if

  IF(N.LT.K) then
    CALL XERMSG ('SLATEC', 'BVALU', 'N DOES NOT SATISFY N.GE.K', 2,1)
    RETURN
  end if

  IF(IDERIV.LT.0 .OR. IDERIV.GE.K) then
    CALL XERMSG ('SLATEC', 'BVALU', &
      'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
    RETURN
  end if

  KMIDER = K - IDERIV
!
!  FIND I IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!  (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
!
  KM1 = K - 1
  CALL INTRV(T, N+1, X, INBV, I, MFLAG)

  IF (X.LT.T(K)) then
    CALL XERMSG ('SLATEC', 'BVALU', &
      'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
    RETURN
  end if

  IF (MFLAG.EQ.0) GO TO 20

  IF (X.GT.T(I)) then
    CALL XERMSG ('SLATEC', 'BVALU', &
      'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
    RETURN
  end if

   10 continue

  IF (I.EQ.K) GO TO 140

  I = I - 1

  IF (X.EQ.T(I)) GO TO 10
!
!  DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
!  WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
!
   20 continue

      IMK = I - K
      DO J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
      end do

      DO J = 1, IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO JJ = 1, KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
        end do
      end do
!
!  COMPUTE VALUE AT X IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!  GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
!
   60 continue

      IF (IDERIV.EQ.KM1) GO TO 100

      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
      end do
      IDERP1 = IDERIV + 1
      DO J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ) &
                    *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
        end do
      end do

  100 continue

      BVALU = WORK(1)

      RETURN

  140 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU', &
        'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)

  RETURN
END
FUNCTION CHFCM ( D1, D2, DELTA )

!*****************************************************************************80
!
!! CHFCM checks a single cubic for monotonicity.
!
!  Discussion:
!
!    Called by PCHCM to determine the monotonicity properties of the
!    cubic with boundary derivative values D1,D2 and chord slope DELTA.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     D1,D2:IN  are the derivative values at the ends of an interval.
!
!     DELTA:IN  is the data slope over that interval.
!
!     ISMON : indicates the monotonicity of the cubic segment:
!             ISMON = -3  if function is probably decreasing;
!             ISMON = -1  if function is strictly decreasing;
!             ISMON =  0  if function is constant;
!             ISMON =  1  if function is strictly increasing;
!             ISMON =  2  if function is non-monotonic;
!             ISMON =  3  if function is probably increasing.
!           If ABS(ISMON)=3, the derivative values are too close to the
!           boundary of the monotonicity region to declare monotonicity
!           in the presence of roundoff error.
!
  implicit none

  integer ( kind = 4 ) chfcm
  real ( kind = 4 ) D1, D2, DELTA
  integer ( kind = 4 )  ISMON, ITRUE
  real ( kind = 4 ) A, B, EPS, ONE, PHI
  real ( kind = 4 ) r1mach

  SAVE ONE

  DATA  ONE /1.0/
!
!  MACHINE-DEPENDENT PARAMETER, SHOULD BE ABOUT 10*UROUND.
!
  EPS = 10.0E+00 * R1MACH ( 4 )
!
!  MAKE THE CHECK.
!
  IF ( DELTA == 0.0E+00 )  THEN
!
!  CASE OF CONSTANT DATA.
!
     IF ((D1 == 0.0E+00 ) .AND. (D2 == 0.0E+00 ))  THEN
        ISMON = 0
     ELSE
        ISMON = 2
     end if
  ELSE
!
!  DATA IS NOT CONSTANT, PICK UP SIGN.
!
     ITRUE = SIGN (ONE, DELTA)
     A = D1/DELTA
     B = D2/DELTA
     IF ((A.LT. 0.0E+00 ) .OR. (B.LT. 0.0E+00 ))  THEN
        ISMON = 2
     ELSE IF ((A.LE. 3.0E+00 -EPS) .AND. (B.LE. 3.0E+00 -EPS))  THEN
!
!  INSIDE SQUARE (0,3)X(0,3) IMPLIES OK.
!
        ISMON = ITRUE
     ELSE IF ((A.GT. 4.0E+00 +EPS) .AND. (B.GT. 4.0E+00 +EPS))  THEN
!
!  OUTSIDE SQUARE (0,4)X(0,4) IMPLIES NONMONOTONIC.
!
        ISMON = 2
     ELSE
!
!  MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
!
        A = A - 2.0E+00
        B = B - 2.0E+00
        PHI = ((A*A + B*B) + A*B) - 3.0E+00 
        IF (PHI .LT. -EPS)  THEN
           ISMON = ITRUE
        ELSE IF (PHI .GT. EPS)  THEN
           ISMON = 2
        ELSE
!
!  TOO CLOSE TO BOUNDARY TO TELL, IN THE PRESENCE OF ROUND-OFF ERRORS.
!
           ISMON = 3*ITRUE
        end if
     end if
  end if
!
!  RETURN VALUE.
!
  CHFCM = ISMON

  RETURN
END
SUBROUTINE CHFDV ( X1, X2, F1, F2, D1, D2, NE, XE, FE, DE, NEXT, IERR )

!*****************************************************************************80
!
!! CHFDV evaluates a cubic polynomial and derivative in Hermite form.
!
!  Discussion:
!
!    Evaluate a cubic polynomial given in Hermite form and its
!    first derivative at an array of points.  
!
!    Evaluates the cubic polynomial determined by function values
!    F1, F2 and derivatives D1, D2 on interval (X1,X2), together with
!    its first derivative, at the points  XE(J), J=1(1)NE.
!
!    While designed for use by PCHFD, it may be useful directly as an evaluator
!    for a piecewise cubic Hermite function in applications,
!    such as graphing, where the interval is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  NE, NEXT(2), IERR
!        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!        CALL  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!  Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer ( kind = 4 ) array indicating number 
!     of extrapolation points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1 == X2 .
!                (Output arrays have not been changed in either case.)
!
  implicit none

  integer ( kind = 4 )  NE, NEXT(2), IERR
  real ( kind = 4 ) X1, X2, F1, F2, D1, D2, XE(*), FE(*), DE(*)
  integer ( kind = 4 )  I
  real ( kind = 4 ) C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA
!
!  CHECK ARGUMENTS.
!
  IF ( NE .LT. 1 ) then
    IERR = -1
    CALL XERMSG ('SLATEC', 'CHFDV', &
      'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
    RETURN
  end if

  H = X2 - X1

  IF ( H == 0.0E+00 ) then
    IERR = -2
    CALL XERMSG ('SLATEC', 'CHFDV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)
    return
  end if
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN( 0.0E+00, H)
  XMA = MAX( 0.0E+00, H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!
!  DELTA IS NO LONGER NEEDED.
!
  C2 = -(DEL1+DEL1 + DEL2)
  C2T2 = C2 + C2
  C3 = (DEL1 + DEL2)/H
!
!  H, DEL1 AND DEL2 ARE NO LONGER NEEDED.
!
  C3T3 = C3+C3+C3
!
!  EVALUATION LOOP.
!
  DO I = 1, NE

     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
     DE(I) = D1 + X*(C2T2 + X*C3T3)
!
!  COUNT EXTRAPOLATION POINTS.
!
     IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
     IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!
!  NOTE REDUNDANCY.  IF EITHER CONDITION IS TRUE, OTHER IS FALSE.
!
  end do

  RETURN
END
SUBROUTINE CHFEV ( X1, X2, F1, F2, D1, D2, NE, XE, FE, NEXT, IERR )

!*****************************************************************************80
!
!! CHFEV evaluates a cubic polynomial in Hermite form.
!
!  Discussion:
!
!    Evaluate a cubic polynomial given in Hermite form at an
!    array of points.  
!
!    Evaluates the cubic polynomial determined by function values
!    F1, F2 and derivatives D1, D2 on interval (X1,X2) at the points
!    XE(J), J=1(1)NE.
!
!    While designed for use by PCHFE, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  NE, NEXT(2), IERR
!        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!        CALL  CHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
!
!  Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the function is to be
!           evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     NEXT -- (output) integer ( kind = 4 ) array indicating number 
!     of extrapolation points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1 == X2 .
!                (The FE-array has not been changed in either case.)
!
  implicit none

  integer ( kind = 4 )  NE, NEXT(2), IERR
  real ( kind = 4 ) X1, X2, F1, F2, D1, D2, XE(*), FE(*)
  integer ( kind = 4 )  I
  real ( kind = 4 ) C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA
!
!  CHECK ARGUMENTS.
!
  IF (NE .LT. 1) then
    IERR = -1
    CALL XERMSG ('SLATEC', 'CHFEV', &
      'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
    RETURN
  end if

  H = X2 - X1

  IF (H == 0.0E+00 ) then
    IERR = -2
    CALL XERMSG ('SLATEC', 'CHFEV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)
    RETURN
  end if
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN( 0.0E+00, H)
  XMA = MAX( 0.0E+00, H)
!
!  COMPUTE CUBIC COEFFICIENTS EXPANDED ABOUT X1.
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!
!  DELTA IS NO LONGER NEEDED.
!
  C2 = -(DEL1+DEL1 + DEL2)
  C3 = (DEL1 + DEL2)/H
!
!  H, DEL1 AND DEL2 ARE NO LONGER NEEDED.
!
!  EVALUATION LOOP.
!
  DO I = 1, NE
     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
!
!  COUNT EXTRAPOLATION POINTS.
!
     IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
     IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!
!  NOTE REDUNDANCY.  IF EITHER CONDITION IS TRUE, OTHER IS FALSE.
!
  end do

  RETURN
END
FUNCTION CHFIE ( X1, X2, F1, F2, D1, D2, A, B )

!*****************************************************************************80
!
!! CHFIE evaluates the integral of a single cubic for PCHIA.
!
!  Discussion:
!
!    Called by PCHIA to evaluate the integral of a single cubic (in
!    Hermite form) over an arbitrary interval (A,B).
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        REAL  X1, X2, F1, F2, D1, D2, A, B
!        REAL  VALUE, CHFIE
!
!        VALUE = CHFIE (X1, X2, F1, F2, D1, D2, A, B)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     X1,X2 -- (input) endpoints if interval of definition of cubic.
!
!     F1,F2 -- (input) function values at the ends of the interval.
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!
!     A,B -- (input) endpoints of interval of integration.
!
  implicit none

  real chfie
  real ( kind = 4 ) X1, X2, F1, F2, D1, D2, A, B
  real ( kind = 4 ) DTERM, FTERM, H, HALF, PHIA1, PHIA2, PHIB1, PHIB2
  real ( kind = 4 ) PSIA1, PSIA2, PSIB1, PSIB2, SIX, TA1, TA2, TB1, TB2
  real ( kind = 4 ) UA1, UA2, UB1, UB2

  SAVE HALF,SIX
!
!  INITIALIZE.
!
  DATA  HALF /0.5/, SIX /6./

  IF (X1 == X2)  THEN

     CHFIE = 0.0E+00

  ELSE

     H = X2 - X1
     TA1 = (A - X1) / H
     TA2 = (X2 - A) / H
     TB1 = (B - X1) / H
     TB2 = (X2 - B) / H

     UA1 = TA1**3
     PHIA1 = UA1 * ( 2.0E+00 - TA1)
     PSIA1 = UA1 * ( 3.0E+00 *TA1 - 4.0E+00 )
     UA2 = TA2**3
     PHIA2 =  UA2 * ( 2.0E+00 - TA2)
     PSIA2 = -UA2 * ( 3.0E+00 *TA2 - 4.0E+00 )

     UB1 = TB1**3
     PHIB1 = UB1 * ( 2.0E+00 - TB1)
     PSIB1 = UB1 * ( 3.0E+00 *TB1 - 4.0E+00 )
     UB2 = TB2**3
     PHIB2 =  UB2 * ( 2.0E+00 - TB2)
     PSIB2 = -UB2 * ( 3.0E+00 *TB2 - 4.0E+00 )

     FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
     DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)

     CHFIE = (HALF*H) * (FTERM + DTERM)

  end if

  RETURN
END
FUNCTION COMP ( IERACT, IEREXP, LOUT, KPRINT )

!*****************************************************************************80
!
!! COMP compares actual and expected values of error flag.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IERACT, the actual error flag value.
!
!    Input, integer ( kind = 4 ) IEREXP, the expected error flag value.
!
!    Input, integer ( kind = 4 ) LOUT, the output device.
!
!    Input, integer ( kind = 4 ) KPRINT, the printing level.
!
!    Output, logical COMP, is true if the flags match.
!
  implicit none

  integer ( kind = 4 ) IERACT, IEREXP, LOUT, KPRINT
  logical comp

  IF (IERACT .EQ. IEREXP)  THEN
    COMP = .TRUE.
    IF (KPRINT .GE. 3) then
      WRITE ( LOUT, '(a)' ) '     OK.'
    end if
  ELSE
    COMP = .FALSE.
    IF ( KPRINT .GE. 3 ) then
      WRITE ( LOUT, '(a,i5)' ) '  COMPARE FAILED -- IERR = ', ieract
    end if
  ENDIF

  RETURN
END
function d1mach ( i )

!*****************************************************************************80
!
!! D1MACH returns real ( kind = 8 ) real machine constants.
!
!  Discussion:
!
!    Assuming that the internal representation of a real ( kind = 8 ) real
!    number is in base B, with T the number of base-B digits in the mantissa,
!    and EMIN the smallest possible exponent and EMAX the largest possible 
!    exponent, then
!
!      D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      D1MACH(3) = B^(-T), the smallest relative spacing.
!      D1MACH(4) = B^(1-T), the largest relative spacing.
!      D1MACH(5) = log10(B).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
!
  implicit none

  real ( kind = 8 ) d1mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    d1mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    d1mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    d1mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    d1mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    d1mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  end if

  return
end
FUNCTION DBVALU ( T, A, N, K, IDERIV, X, INBV, WORK )

!*****************************************************************************80
!
!! DBVALU evaluates a B-spline or its derivatives.
!
!  Discussion:
!
!    DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
!    at X for the function value on IDERIV=0 or any of its
!    derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!    (right derivatives) are returned except at the right end
!    point X=T(N+1) where left limiting values are computed.  The
!    spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
!    a fatal error message when X is outside of this interval.
!
!    To compute left derivatives or left limiting values at a
!    knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Daniel Amos
!
!  Reference:
!
!    Carl de Boor, 
!    Package for calculating with B-splines,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 3, June 1977, pages 441-472.
!
!  Parameters:
!
!         Input      T,A,X are double precision
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV = 0 returns the B-spline value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time DBVALU is called.
!
!         Output     WORK,DBVALU are double precision
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          DBVALU  - value of the IDERIV-th derivative at X
!
  implicit none

  real ( kind = 8 ) dbvalu
  integer ( kind = 4 ) I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ
  integer ( kind = 4 ) IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER
  integer ( kind = 4 ) KMJ, KM1, KPK, MFLAG, N
  real ( kind = 8 ) A(*), FKMJ, T(*), WORK(*), X

  DBVALU = 0.0D0

  IF(K.LT.1) then
    CALL XERMSG ('SLATEC', 'DBVALU', 'K DOES NOT SATISFY K.GE.1', 2,1)
    RETURN
  end if

  IF ( N .LT. K ) then
    CALL XERMSG ('SLATEC', 'DBVALU', 'N DOES NOT SATISFY N.GE.K', 2, 1)
    RETURN
  end if

  IF(IDERIV.LT.0 .OR. IDERIV.GE.K) then
    CALL XERMSG ('SLATEC', 'DBVALU', &
      'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
    RETURN
  end if

  KMIDER = K - IDERIV
!
!  FIND I IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!  (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
!
  KM1 = K - 1

  CALL DINTRV ( T, N+1, X, INBV, I, MFLAG )

  IF (X.LT.T(K)) then
   CALL XERMSG ('SLATEC', 'DBVALU', &
     'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
    RETURN
  end if

  IF (MFLAG.EQ.0) GO TO 20

  IF (X.GT.T(I)) GO TO 130

10 continue

  IF (I.EQ.K) GO TO 140

  I = I - 1

  IF (X.EQ.T(I)) GO TO 10
!
!  DIFFERENCE THE COEFFICIENTS IDERIV TIMES
!  WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
!
20 continue

  IMK = I - K

  DO J=1,K
    IMKPJ = IMK + J
    WORK(J) = A(IMKPJ)
  end do

  DO J=1,IDERIV
    KMJ = K - J
    FKMJ = KMJ
    DO JJ = 1, KMJ
      IHI = I + JJ
      IHMKMJ = IHI - KMJ
      WORK(JJ) = ( WORK(JJ+1) - WORK(JJ) ) / ( T(IHI) - T(IHMKMJ) ) * FKMJ
    end do
  end do
!
!  COMPUTE VALUE AT X IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!  GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
!
  IF ( IDERIV /= KM1 ) then

    IP1 = I + 1
    KPK = K + K
    J1 = K + 1
    J2 = KPK + 1

    DO J=1,KMIDER
      IPJ = I + J
      WORK(J1) = T(IPJ) - X
      IP1MJ = IP1 - J
      WORK(J2) = X - T(IP1MJ)
      J1 = J1 + 1
      J2 = J2 + 1
    end do

    DO J=IDERIV + 1, KM1
      KMJ = K - J
      ILO = KMJ
      DO JJ=1,KMJ
        WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ) &
                  *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
        ILO = ILO - 1
      end do
    end do

  end if

  DBVALU = WORK(1)

  RETURN

  130 CONTINUE
  CALL XERMSG ('SLATEC', 'DBVALU', &
        'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
  RETURN

  140 CONTINUE
  CALL XERMSG ('SLATEC', 'DBVALU', &
        'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
  RETURN
END
FUNCTION DCHFCM ( D1, D2, DELTA )

!*****************************************************************************80
!
!! DCHFCM checks a single cubic for monotonicity.
!
!  Discussion:
!
!    Called by DPCHCM to determine the monotonicity properties of the
!    cubic with boundary derivative values D1, D2 and chord slope DELTA.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     D1,D2:IN  are the derivative values at the ends of an interval.
!
!     DELTA:IN  is the data slope over that interval.
!
!     ISMON : indicates the monotonicity of the cubic segment:
!             ISMON = -3  if function is probably decreasing;
!             ISMON = -1  if function is strictly decreasing;
!             ISMON =  0  if function is constant;
!             ISMON =  1  if function is strictly increasing;
!             ISMON =  2  if function is non-monotonic;
!             ISMON =  3  if function is probably increasing.
!           If ABS(ISMON)=3, the derivative values are too close to the
!           boundary of the monotonicity region to declare monotonicity
!           in the presence of roundoff error.
!
  implicit none

  integer ( kind = 4 ) dchfcm
  real ( kind = 8 )  D1, D2, DELTA, D1MACH
  integer ( kind = 4 ) ISMON, ITRUE
  real ( kind = 8 )  A, B, EPS, ONE, PHI

  SAVE ONE
!
!  INITIALIZE.
!
  DATA ONE/1.D0/
!
!  MACHINE-DEPENDENT PARAMETER.  SHOULD BE ABOUT 10*UROUND.
!
  EPS = 10.0D+00 * D1MACH(4)
!
!  MAKE THE CHECK.
!
  IF (DELTA == 0.0D+00 )  THEN
!
!  CASE OF CONSTANT DATA.
!
     IF ((D1 ==  0.0D+00 ) .AND. (D2 ==  0.0D+00 ))  THEN
        ISMON = 0
     ELSE
        ISMON = 2
     end if
  ELSE
!
!  DATA IS NOT CONSTANT.  PICK UP SIGN.
!
     ITRUE = DSIGN (ONE, DELTA)
     A = D1/DELTA
     B = D2/DELTA
     IF ((A.LT. 0.0D+00 ) .OR. (B.LT. 0.0D+00 ))  THEN
        ISMON = 2
     ELSE IF ((A.LE. 3.0D+00 -EPS) .AND. (B.LE. 3.0D+00 -EPS))  THEN
!
!  INSIDE SQUARE (0,3)X(0,3) IMPLIES OK.
!
        ISMON = ITRUE
     ELSE IF ((A.GT. 4.0D+00 +EPS) .AND. (B.GT. 4.0D+00 +EPS))  THEN
!
!  OUTSIDE SQUARE (0,4)X(0,4) IMPLIES NONMONOTONIC.
!
        ISMON = 2
     ELSE
!
!  MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
!
        A = A - 2.0D+00
        B = B - 2.0D+00
        PHI = ((A*A + B*B) + A*B) - 3.0D+00 
        IF (PHI .LT. -EPS)  THEN
           ISMON = ITRUE
        ELSE IF (PHI .GT. EPS)  THEN
           ISMON = 2
        ELSE
!
!  TOO CLOSE TO BOUNDARY TO TELL, IN THE PRESENCE OF ROUND-OFF ERRORS.
!
           ISMON = 3*ITRUE
        end if
     end if
  end if

  DCHFCM = ISMON

  RETURN
END
SUBROUTINE DCHFDV ( X1, X2, F1, F2, D1, D2, NE, XE, FE, DE, NEXT, IERR )

!*****************************************************************************80
!
!! DCHFDV evaluates a cubic Hermite polynomial and derivative at many points.
!
!  Discussion:
!
!    Evaluate a cubic polynomial given in Hermite form and its
!    first derivative at an array of points.  
!
!    Evaluates the cubic polynomial determined by function values
!    F1, F2 and derivatives D1, D2 on interval (X1,X2), together with
!    its first derivative, at the points  XE(J), J=1(1)NE.
!
!    While designed for use by DPCHFD, it may be useful directly as an 
!    evaluator for a piecewise cubic Hermite function in applications,
!    such as graphing, where the interval is known in advance.
!
!    If only function values are required, use DCHFEV instead.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  NE, NEXT(2), IERR
!        real ( kind = 8 )  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE),
!                          DE(NE)
!
!        CALL  DCHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!  Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer ( kind = 4 ) array indicating number of 
!     extrapolation points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1 == X2 .
!                (Output arrays have not been changed in either case.)
!
  implicit none

  integer ( kind = 4 )  NE, NEXT(2), IERR
  real ( kind = 8 )  X1, X2, F1, F2, D1, D2, XE(*), FE(*), DE(*)
  integer ( kind = 4 )  I
  real ( kind = 8 )  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA
!
!  CHECK ARGUMENTS.
!
  IF (NE .LT. 1)  GO TO 5001
  H = X2 - X1
  IF (H == 0.0D+00 )  GO TO 5002
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN ( 0.0D+00, H)
  XMA = MAX ( 0.0D+00, H)
!
!  COMPUTE CUBIC COEFFICIENTS EXPANDED ABOUT X1.
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!
!  DELTA IS NO LONGER NEEDED.
!
  C2 = -(DEL1+DEL1 + DEL2)
  C2T2 = C2 + C2
  C3 = (DEL1 + DEL2)/H
!
!  H, DEL1 AND DEL2 ARE NO LONGER NEEDED.
!
  C3T3 = C3+C3+C3
!
!  EVALUATION LOOP.
!
  DO I = 1, NE
     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
     DE(I) = D1 + X*(C2T2 + X*C3T3)
!
!  COUNT EXTRAPOLATION POINTS.
!
     IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
     IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!
!  NOTE REDUNDANCY.  IF EITHER CONDITION IS TRUE, OTHER IS FALSE.
!
  end do

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  NE.LT.1 RETURN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DCHFDV', &
     'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  X1 == X2 RETURN.
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'DCHFDV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)

  RETURN
END
SUBROUTINE DCHFEV ( X1, X2, F1, F2, D1, D2, NE, XE, FE, NEXT, IERR )

!*****************************************************************************80
!
!! DCHFEV evaluates a cubic Hermite polynomial at many points.
!
!  Discussion:
!
!    Evaluate a cubic polynomial given in Hermite form at an
!    array of points.  
!
!    Evaluates the cubic polynomial determined by function values
!    F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
!    XE(J), J=1(1)NE.
!
!    While designed for use by DPCHFE, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  NE, NEXT(2), IERR
!        real ( kind = 8 )  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!        CALL  DCHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
!
!  Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     NEXT -- (output) integer ( kind = 4 ) array indicating number of 
!     extrapolation points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1 == X2 .
!                (The FE-array has not been changed in either case.)
!
  implicit none

  integer ( kind = 4 )  NE, NEXT(2), IERR
  real ( kind = 8 )  X1, X2, F1, F2, D1, D2, XE(*), FE(*)
  integer ( kind = 4 )  I
  real ( kind = 8 )  C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA
!
!  CHECK ARGUMENTS.
!
  IF (NE .LT. 1)  GO TO 5001
  H = X2 - X1
  IF (H == 0.0D+00 )  GO TO 5002
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN ( 0.0D+00, H)
  XMA = MAX ( 0.0D+00, H)
!
!  COMPUTE CUBIC COEFFICIENTS EXPANDED ABOUT X1.
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!
!  DELTA IS NO LONGER NEEDED.
!
  C2 = -(DEL1+DEL1 + DEL2)
  C3 = (DEL1 + DEL2)/H
!
!  H, DEL1 AND DEL2 ARE NO LONGER NEEDED.
!
!  EVALUATION LOOP.
!
  DO I = 1, NE
     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
!
!  COUNT EXTRAPOLATION POINTS.
!
     IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
     IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!
!  NOTE REDUNDANCY.  IF EITHER CONDITION IS TRUE, OTHER IS FALSE.
!
  end do

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  NE.LT.1 RETURN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DCHFEV', &
   'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  X1 == X2 RETURN.
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'DCHFEV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)

  RETURN
END
FUNCTION DCHFIE ( X1, X2, F1, F2, D1, D2, A, B )

!*****************************************************************************80
!
!! DCHFIE evaluates the integral of a single cubic for DPCHIA.
!
!  Discussion:
!
!    Called by DPCHIA to evaluate the integral of a single cubic (in
!    Hermite form) over an arbitrary interval (A,B).
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        real ( kind = 8 )  X1, X2, F1, F2, D1, D2, A, B
!        real ( kind = 8 )  VALUE, DCHFIE
!
!        VALUE = DCHFIE (X1, X2, F1, F2, D1, D2, A, B)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     X1,X2 -- (input) endpoints if interval of definition of cubic.
!
!     F1,F2 -- (input) function values at the ends of the interval.
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!
!     A,B -- (input) endpoints of interval of integration.
!
  implicit none

  real ( kind = 8 ) dchfie
  real ( kind = 8 )  X1, X2, F1, F2, D1, D2, A, B
  real ( kind = 8 )  DTERM, FTERM, H, HALF, PHIA1, PHIA2
  real ( kind = 8 )  PHIB1, PHIB2, PSIA1, PSIA2, PSIB1, PSIB2, SIX, TA1, TA2
  real ( kind = 8 )  TB1, TB2, UA1, UA2, UB1, UB2

  SAVE HALF, SIX
!
!  INITIALIZE.
!
  DATA  HALF/.5D0/, SIX/6.D0/

  IF (X1 == X2)  THEN

     DCHFIE = 0.0D+00

  ELSE

     H = X2 - X1
     TA1 = (A - X1) / H
     TA2 = (X2 - A) / H
     TB1 = (B - X1) / H
     TB2 = (X2 - B) / H

     UA1 = TA1**3
     PHIA1 = UA1 * ( 2.0D+00 - TA1)
     PSIA1 = UA1 * ( 3.0D+00 *TA1 - 4.0D+00 )
     UA2 = TA2**3
     PHIA2 =  UA2 * ( 2.0D+00 - TA2)
     PSIA2 = -UA2 * ( 3.0D+00 *TA2 - 4.0D+00 )

     UB1 = TB1**3
     PHIB1 = UB1 * ( 2.0D+00 - TB1)
     PSIB1 = UB1 * ( 3.0D+00 *TB1 - 4.0D+00 )
     UB2 = TB2**3
     PHIB2 =  UB2 * ( 2.0D+00 - TB2)
     PSIB2 = -UB2 * ( 3.0D+00 *TB2 - 4.0D+00 )

     FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
     DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)

     DCHFIE = (HALF*H) * (FTERM + DTERM)

  end if

  RETURN
END
SUBROUTINE DEVCHK ( LOUT, KPRINT, NPTS, XEV, FEV, DEV, FEV2, FAIL )

!*****************************************************************************80
!
!! DEVCHK tests evaluation accuracy of DCHFDV and DCHFEV for DPCHQ1.
!
!  Discussion:
!
!    USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
!    DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
!    1. CHECKS THAT DCHFDV AND DCHFEV BOTH REPRODUCE ENDPOINT VALUES.
!    2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
!    AND:
!    A. CHECKS ACCURACY OF DCHFDV FUNCTION AND DERIVATIVE VALUES
!       AGAINST EXACT VALUES.
!    B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
!    C. CHECKS THAT FUNCTION VALUES FROM DCHFEV AGREE WITH THOSE
!       FROM DCHFDV.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    LOUT
!
!    KPRINT
!
!    NPTS
!
!    XEV
!
!    FEV
!
!    DEV
!
!    FEV2
!
!    FAIL
!
  implicit none

  integer ( kind = 4 )  LOUT, KPRINT, NPTS
  real ( kind = 8 ) XEV(*), FEV(*), DEV(*), FEV2(*)
  LOGICAL  FAIL

  integer ( kind = 4 ) I, IERR, IINT, NEXT(2), NEXT2(2), NINT
  real ( kind = 8 ) &
       AED, AED2, AEDMAX, AEDMIN, AEF, AEF2, AEFMAX, AEFMIN, &
       CHECK(2), CHECKF(2), CHECKD(2), D1, D2, DERMAX, DTRUE, DX, &
       EPS1, EPS2, F1, F2, FACT, FERMAX, FLOORD, FLOORF, FOUR, &
       FTRUE, LEFT(3), MACHEP, &
       ONE, RED, RED2, REDMAX, REDMIN, REF, REF2, REFMAX, &
       REFMIN, RIGHT(3), SMALL, TEN, TOL1, TOL2, &
       X1, X2, XADMAX, XADMIN, XAFMAX, XAFMIN, XRDMAX, &
       XRDMIN, XRFMAX, XRFMIN, ZERO
  LOGICAL  FAILOC, FAILNX

  real ( kind = 8 )  D1MACH

  real ( kind = 4 )  RAND
  EXTERNAL  RAND
!
!  DEFINE RELATIVE ERROR WITH FLOOR.
!
  real ( kind = 8 ) RERR, ERR, VALUE, FLOOR
  RERR(ERR,VALUE,FLOOR) = ERR / MAX(ABS(VALUE), FLOOR)
!
!  INITIALIZE.
!
  DATA  ZERO /0.D0/, ONE /1.D0/, FOUR /4.D0/, TEN /10.D0/
  DATA  SMALL  /1.0D-10/
  DATA  NINT /3/
  DATA   LEFT /-1.5D0, 2.0D-10, 1.0D0 /
  DATA  RIGHT / 2.5D0, 3.0D-10, 1.0D+8/

  MACHEP = D1MACH ( 4 )
  EPS1 = FOUR * MACHEP
  EPS2 = TEN * MACHEP

  FAIL = .FALSE.

  IF (KPRINT .GE. 2) then
    WRITE (LOUT, 3000)
  end if
!
!  CYCLE OVER INTERVALS.
!
  DO IINT = 1, NINT

    X1 =  LEFT(IINT)
    X2 = RIGHT(IINT)

    FACT = MAX ( SQRT ( X2 - X1 ), ONE )
    TOL1 = EPS1 * FACT
    TOL2 = EPS2 * FACT
!
!  COMPUTE AND PRINT ENDPOINT VALUES.
!
    CALL DFDTRU (X1, F1, D1)
    CALL DFDTRU (X2, F2, D2)

    IF (KPRINT .GE. 3)  THEN
      IF (IINT .EQ. 1)  WRITE (LOUT, 2000)
      WRITE (LOUT, '(/)')
      WRITE (LOUT, 2001)  'X1', X1, 'X2', X2
      WRITE (LOUT, 2001)  'F1', F1, 'F2', F2
      WRITE (LOUT, 2001)  'D1', D1, 'D2', D2
   END IF

    IF (KPRINT .GE. 2)  WRITE (LOUT, 3001)  X1, X2
!
!  COMPUTE FLOORS FOR RELATIVE ERRORS.
!
    FLOORF = MAX( MIN(ABS(F1),ABS(F2)), SMALL)
    FLOORD = MAX( MIN(ABS(D1),ABS(D2)), SMALL)
!
!  CHECK REPRODUCTION OF ENDPOINT VALUES.
!
    XEV(1) = X1
    XEV(2) = X2

    CALL DCHFDV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECKF, CHECKD, &
        NEXT, IERR)

    AEF  = CHECKF(1)-F1
    REF  = RERR(AEF , F1, FLOORF)
    AEF2 = CHECKF(2)-F2
    REF2 = RERR(AEF2, F2, FLOORF)
    AED  = CHECKD(1)-D1
    RED  = RERR(AED , D1, FLOORD)
    AED2 = CHECKD(2)-D2
    RED2 = RERR(AED2, D2, FLOORD)

    FAILOC = MAX(ABS(REF),ABS(REF2),ABS(RED),ABS(RED2)) .GT. TOL1
    FAIL = FAIL .OR. FAILOC

    IF (KPRINT .GE. 3)  THEN
      WRITE (LOUT, 2002)  NEXT, AEF, AEF2, AED, AED2
      WRITE (LOUT, 2003)  REF, REF2, RED, RED2
    ENDIF

    IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3002)
!
!  DCHFEV SHOULD AGREE EXACTLY WITH DCHFDV.
!
    CALL DCHFEV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECK, NEXT, IERR)
  
    FAILOC = (CHECK(1).NE.CHECKF(1)) .OR. (CHECK(2).NE.CHECKF(2))
    FAIL = FAIL .OR. FAILOC

    IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3003)
!
!  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
!  THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
!  TO LEFT AND 6 TO RIGHT OF [X1,X2].
!
    DX = (X2-X1)/(NPTS-10)
    DO I = 1, NPTS
      XEV(I) = (X1 + (I-5)*DX) + DX*RAND(0.0E0)
    end do

    CALL DCHFDV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV, DEV, NEXT, IERR)

    IF (IERR .NE. 0)  THEN
      FAILOC = .TRUE.
      IF (KPRINT .GE. 2)  WRITE (LOUT, 4003)  IERR
    ELSE
!
!  CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
!
      DO I = 1, NPTS

        CALL DFDTRU (XEV(I), FTRUE, DTRUE)
        AEF = FEV(I) - FTRUE
        REF = RERR(AEF, FTRUE, FLOORF)
        AED = DEV(I) - DTRUE
        RED = RERR(AED, DTRUE, FLOORD)

        IF (I .EQ. 1)  THEN
!
!  INITIALIZE.
!
          AEFMIN = AEF
          AEFMAX = AEF
          AEDMIN = AED
          AEDMAX = AED
          REFMIN = REF
          REFMAX = REF
          REDMIN = RED
          REDMAX = RED
          XAFMIN = XEV(1)
          XAFMAX = XEV(1)
          XADMIN = XEV(1)
          XADMAX = XEV(1)
          XRFMIN = XEV(1)
          XRFMAX = XEV(1)
          XRDMIN = XEV(1)
          XRDMAX = XEV(1)

        ELSE
!
!  SELECT.
!
          IF (AEF .LT. AEFMIN)  THEN
            AEFMIN = AEF
            XAFMIN = XEV(I)
          ELSE IF (AEF .GT. AEFMAX)  THEN
            AEFMAX = AEF
            XAFMAX = XEV(I)
          ENDIF

          IF (AED .LT. AEDMIN)  THEN
            AEDMIN = AED
            XADMIN = XEV(I)
          ELSE IF (AED .GT. AEDMAX)  THEN
            AEDMAX = AED
            XADMAX = XEV(I)
          ENDIF

          IF (REF .LT. REFMIN)  THEN
            REFMIN = REF
            XRFMIN = XEV(I)
          ELSE IF (REF .GT. REFMAX)  THEN
            REFMAX = REF
            XRFMAX = XEV(I)
          ENDIF

          IF (RED .LT. REDMIN)  THEN
            REDMIN = RED
            XRDMIN = XEV(I)
          ELSE IF (RED .GT. REDMAX)  THEN
            REDMAX = RED
            XRDMAX = XEV(I)
          ENDIF

        ENDIF

      end do

      FERMAX = MAX (ABS(REFMAX), ABS(REFMIN))
      DERMAX = MAX (ABS(REDMAX), ABS(REDMIN))

      FAILNX = (NEXT(1) + NEXT(2)) .NE. 10
      FAILOC = FAILNX .OR. (MAX(FERMAX, DERMAX) .GT. TOL2)

    ENDIF

    FAIL = FAIL .OR. FAILOC
!
!  PRINT SUMMARY.
!
    IF (KPRINT .GE. 3)  THEN
      WRITE (LOUT, 2004)  NPTS-10, NEXT
      WRITE (LOUT, 2005)  'MIN', AEFMIN, REFMIN, AEDMIN, REDMIN
      WRITE (LOUT, 2006) XAFMIN, XRFMIN, XADMIN, XRDMIN
      WRITE (LOUT, 2005)  'MAX', AEFMAX, REFMAX, AEDMAX, REDMAX
      WRITE (LOUT, 2006) XAFMAX, XRFMAX, XADMAX, XRDMAX
    ENDIF

    IF (KPRINT .GE. 2)  THEN
      IF (FAILOC) THEN
        IF (FERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'F', FERMAX, TOL2
        IF (DERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'D', DERMAX, TOL2
        IF (FAILNX)  WRITE (LOUT, 4006)  NEXT
      ELSE
        WRITE (LOUT, 5006)
      ENDIF
    ENDIF
!
!  CHECK THAT DCHFEV AGREES WITH DCHFDV.
!
    CALL DCHFEV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV2, NEXT2, IERR)

    IF (IERR .NE. 0)  THEN
      FAILOC = .TRUE.
      IF (KPRINT .GE. 2)  WRITE (LOUT, 3007)  IERR
    ELSE
      AEFMAX = ABS(FEV2(1) - FEV(1))
      XAFMAX = XEV(1)
      DO I = 2, NPTS
        AEF = ABS(FEV2(I) - FEV(I))
        IF (AEF .GT. AEFMAX)  THEN
          AEFMAX = AEF
          XAFMAX = XEV(I)
        ENDIF
      end do
      FAILNX = (NEXT2(1).NE.NEXT(1)) .OR. (NEXT2(2).NE.NEXT(2))
      FAILOC = FAILNX .OR. (AEFMAX.NE.ZERO)
      IF (KPRINT .GE. 2)  THEN
        IF (FAILOC)  THEN
          WRITE (LOUT, 3008)
          IF (AEFMAX.NE.ZERO)  WRITE (LOUT, 3009)  AEFMAX, XAFMAX
          IF (FAILNX)  WRITE (LOUT, 4009)  NEXT2, NEXT
        ELSE
          WRITE (LOUT, 5009)
        ENDIF
      ENDIF
    ENDIF

    FAIL = FAIL .OR. FAILOC
!
!  GO BACK FOR ANOTHER INTERVAL.
!
  end do

  RETURN
!
!  FORMATS.
!
 2000 FORMAT (/10X,'DCHFDV ACCURACY TEST')
 2001 FORMAT (10X,A2,' =',1P,D18.10,5X,A2,' =',D18.10)
 2002 FORMAT (/' ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')' &
     // 1P,4X,'F1:',D13.5,4X,'F2:',D13.5, &
       4X,'D1:',D13.5,4X,'D2:',D13.5) 
 2003 FORMAT (1P,4(7X,D13.5))
 2004 FORMAT (/' ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:', &
         15X,'(NEXT =',2I3,')' &
     //30X,'FUNCTION',17X,'DERIVATIVE' &
      /15X,2(11X,'ABS',9X,'REL') )
 2005 FORMAT (/5X,A3,'IMUM ERROR:  ',1P,2D12.4,2X,2D12.4)
 2006 FORMAT ( 5X,'LOCATED AT X =  ',1P,2D12.4,2X,2D12.4)
 3000 FORMAT (//10X,'DEVCHK RESULTS'/10X,'--------------')
 3001 FORMAT (/10X,'INTERVAL = (',1P,D12.5,',',D12.5,' ):' )
 3002 FORMAT (/' ***** DCHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
 3003 FORMAT (/' ***** DCHFEV DOES NOT AGREE WITH DCHFDV AT ENDPOINTS.')
 3006 FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',1P,D12.5,',' &
     /    17X,'EXCEEDS TOLERANCE =',D12.5)
 3007 FORMAT (/' ***** ERROR ***** DCHFEV RETURNED IERR =',I5)
 3008 FORMAT (/' ***** DCHFEV DID NOT AGREE WITH DCHFDV:')
 3009 FORMAT (7X,'MAXIMUM DIFFERENCE ',1P,D12.5, &
         '; OCCURRED AT X =',D12.5)
 4003 FORMAT (/' ***** ERROR ***** DCHFDV RETURNED IERR =',I5)
 4006 FORMAT (/' ***** REPORTED NEXT =',2I5,'   RATHER THAN    4    6')
 4009 FORMAT (7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
 5006 FORMAT (/' DCHFDV RESULTS OK.')
 5009 FORMAT (/' DCHFEV AGREES WITH DCHFDV.')
END
SUBROUTINE DEVERK (LOUT, KPRINT, FAIL)

!*****************************************************************************80
!
!! DEVERK tests error returns from DPCHIP evaluators for DPCHQ1.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOUT, the output unit.
!
!    KPRINT
!
!    FAIL
!
  implicit none

  integer ( kind = 4 ) LOUT, KPRINT
  LOGICAL  FAIL

  integer ( kind = 4 ) I, IERR, KONTRL, N, NERR, NEXT(2)
  real ( kind = 8 ) D(10), DUM, F(10), TEMP, X(10)
  LOGICAL  COMP, SKIP

  PARAMETER (N = 10)

  NERR = 0

  CALL XGETF (KONTRL)
  IF (KPRINT .LE. 2) THEN
     CALL XSETF (0)
  ELSE
     CALL XSETF (1)
  ENDIF

  IF (KPRINT .GE. 3)  WRITE (LOUT, 2000)
  IF (KPRINT .GE. 2)  WRITE (LOUT, 5000)
!
!  FIRST, TEST DCHFEV AND DCHFDV.
!
  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  CALL DCHFEV (0.D0, 1.D0, 3.D0, 7.D0, 3.D0, 6.D0, 0, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
  CALL DCHFEV (1.D0, 1.D0, 3.D0, 7.D0, 3.D0, 6.D0, 1, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  CALL DCHFDV (0.D0, 1.D0, 3.D0, 7.D0, 3.D0, 6.D0, 0, DUM, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
  CALL DCHFDV (1.D0, 1.D0, 3.D0, 7.D0, 3.D0, 6.D0, 1, DUM, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1
!
!  SET UP PCH DEFINITION.
!
  DO I = 1, N
    X(I) = I
    F(I) = I + 2
    D(I) = 1.D0
  end do
!
!  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
!
  TEMP = X(4)
  X(4) = X(7)
  X(7) = TEMP
!
!  NOW, TEST DPCHFE AND DPCHFD.
!
  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  SKIP = .FALSE.
  CALL DPCHFE (1, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
  SKIP = .FALSE.
  CALL DPCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
  SKIP = .TRUE.
  CALL DPCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  SKIP = .FALSE.
  CALL DPCHFD (1, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
  SKIP = .FALSE.
  CALL DPCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
  SKIP = .TRUE.
  CALL DPCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
  IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1
!
!  SUMMARIZE RESULTS.
!
  IF (KPRINT .GT. 2)  CALL XERDMP
  IF (NERR .EQ. 0)  THEN
    FAIL = .FALSE.
    IF (KPRINT .GE. 2)  WRITE (LOUT, 5002)
  ELSE
    FAIL = .TRUE.
    IF (KPRINT .GE. 2)  WRITE (LOUT, 5003)  NERR
  ENDIF

  CALL XSETF (KONTRL)

  RETURN
 2000 FORMAT ('1'//10X,'TEST ERROR RETURNS')
 5000 FORMAT (//10X,'DEVERK RESULTS'/10X,'--------------')
 5001 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
 5002 FORMAT (/' ALL ERROR RETURNS OK.')
 5003 FORMAT (//' ***** TROUBLE IN DEVERK *****' &
     //5X,I5,' TESTS FAILED TO GIVE EXPECTED RESULTS.')
END
SUBROUTINE DEVPCK (LOUT, KPRINT, X, Y, F, FX, FY, XE, YE, FE, DE, FE2, FAIL)

!*****************************************************************************80
!
!! DEVPCK tests usage of increment argument in DPCHFD and DPCHFE for DPCHQ1.
!  
!  Discussion:
!
!    EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
!    ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
!
!    INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
!    SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
!
!    ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
!    TEST ROUTINES.
!
!    RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
!   
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOUT, the output unit.
!
  implicit none

  integer ( kind = 4 ) LOUT, KPRINT
  LOGICAL  FAIL
  real ( kind = 8 ) X(10), Y(10), F(10,10), FX(10,10), FY(10,10)
  real ( kind = 8 ) XE(51), YE(51), FE(51), DE(51), FE2(51)

  integer ( kind = 4 )I, IER2, IERR, INC, J, K, NE, NERR, NMAX, NX, NY
  LOGICAL  FAILD, FAILE, FAILOC, SKIP
  real ( kind = 8 ) DERMAX, DERR, DTRUE, DX, FDIFF, FDIFMX, FERMAX, FERR
  real ( kind = 8 ) FTRUE, MACHEP, TOL, PDERMX, PDIFMX, PFERMX, ZERO
  real ( kind = 8 ) D1MACH
!
!  DEFINE TEST FUNCTION AND DERIVATIVES.
!
  real ( kind = 8 ) AX, AY, FCN, DFDX, DFDY
  FCN (AX,AY) =  AX*(AY*AY)*(AX*AX + 1.D0)
  DFDX(AX,AY) = (AY*AY)*(3.D0*AX*AX + 1.D0)
  DFDY(AX,AY) =   2.D0*AX*AY*(AX*AX + 1.D0)

  DATA  NMAX /10/,  NX /4/,  NY /6/
  DATA  NE /51/
  DATA  ZERO /0.D0/

  MACHEP = D1MACH(4)
!
!  Following tolerance is looser than S.P. version to avoid
!  spurious failures on some systems.
!
  TOL = 25.D0*MACHEP

  FAIL = .FALSE.
!
!  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
!  X =  0.25(0.25)1.   ;
!  Y = -0.75(0.5 )1.75 .
!
  DO I = 1, NX-1
    X(I) = 0.25D0*I
  end do
  X(NX) = 1.D0

  DO J = 1, NY
    Y(J) = 0.5D0*J - 1.25D0
    DO I = 1, NX
      F(I,J) = FCN (X(I), Y(J))
      FX(I,J) = DFDX(X(I), Y(J))
      FY(I,J) = DFDY(X(I), Y(J))
    end do
  end do
!
!  SET UP EVALUATION POINTS:
!  XE =  0.(0.02)1. ;
!  YE = -2.(0.08)2. .
!
  DX = 1.D0/(NE-1)
  DO K = 1, NE-1
    XE(K) = DX*(K-1)
    YE(K) = 4.D0*XE(K) - 2.D0
  end do
  XE(NE) = 1.D0
  YE(NE) = 2.D0

  IF (KPRINT .GE. 3)  WRITE (LOUT, 1000)
  IF (KPRINT .GE. 2)  WRITE (LOUT, 1001)
!
!  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING).
!
  NERR = 0
  INC = 1
  SKIP = .FALSE.

  DO J = 1, NY

     CALL DPCHFD (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE, DE, IERR)

     IF (KPRINT .GE. 3) then
       WRITE (LOUT, 2000)  INC, 'J', J, 'Y', Y(J), IERR
     end if

     IF (IERR .LT. 0)  GO TO 15

     IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'X'
!
!  DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
!
     CALL DPCHFE (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE2, IER2)

     DO K = 1, NE

       FTRUE =  FCN(XE(K), Y(J))
       FERR = FE(K) - FTRUE
       DTRUE = DFDX(XE(K), Y(J))
       DERR = DE(K) - DTRUE
       IF (KPRINT .GT. 3) then
         WRITE (LOUT, 2002)  XE(K), FTRUE, FE(K), FERR, DTRUE, DE(K), DERR
       end if

       IF (K .EQ. 1)  THEN
!
!  INITIALIZE.
!
       FERMAX = ABS(FERR)
       PFERMX = XE(1)
       DERMAX = ABS(DERR)
       PDERMX = XE(1)
       FDIFMX = ABS(FE2(1) - FE(1))
       PDIFMX = XE(1)
    ELSE
!
!  SELECT.
!
       FERR = ABS(FERR)
       IF (FERR .GT. FERMAX)  THEN
      FERMAX = FERR
      PFERMX = XE(K)
       ENDIF
       DERR = ABS(DERR)
       IF (DERR .GT. DERMAX)  THEN
      DERMAX = DERR
      PDERMX = XE(K)
       ENDIF
       FDIFF = ABS(FE2(K) - FE(K))
       IF (FDIFF .GT. FDIFMX)  THEN
      FDIFMX = FDIFF
      PDIFMX = XE(K)
       ENDIF
    ENDIF

     end do

     FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
     FAILE = FDIFMX .NE. ZERO
     FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.13) .OR. (IER2.NE.IERR)

     IF (FAILOC .AND. (KPRINT.GE.2)) then
       WRITE (LOUT, 2003)  'J', J, 'Y', Y(J)
     end if

     IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) ) then
    WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
     end if

     IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL

     IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) ) then
       WRITE (LOUT, 2005)  FDIFMX, PDIFMX
     end if

     IF ((IERR.NE.13) .AND. (KPRINT.GE.2)) then
       WRITE (LOUT, 2006)  'D', IERR, 13
     end if

     IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2)) then
       WRITE (LOUT, 2006)  'E', IER2, IERR
     end if

     GO TO 19

   15    CONTINUE

     FAILOC = .TRUE.
     IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR

   19    CONTINUE

     IF (FAILOC)  NERR = NERR + 1
     FAIL = FAIL .OR. FAILOC

  end do

  IF (KPRINT .GE. 2)  THEN
     IF (NERR .GT. 0)  THEN
    WRITE (LOUT, 3001)  NERR, 'J'
     ELSE
    WRITE (LOUT, 4000)  'J'
     ENDIF
  ENDIF
!
!  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING).
!
  NERR = 0
  INC = NMAX
  SKIP = .FALSE.

  DO I = 1, NX

    CALL DPCHFD (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE, DE, IERR)

    IF (KPRINT .GE. 3) then
      WRITE (LOUT, 2000)  INC, 'I', I, 'X', X(I), IERR
    end if

    IF (IERR .LT. 0)  GO TO 35
    IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'Y'
!
!  DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
!
    CALL DPCHFE (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE2, IER2)

    DO K = 1, NE

      FTRUE =  FCN(X(I), YE(K))
      FERR = FE(K) - FTRUE
      DTRUE = DFDY(X(I), YE(K))
      DERR = DE(K) - DTRUE
      IF (KPRINT .GT. 3) then
        WRITE (LOUT, 2002)  YE(K), FTRUE, FE(K), FERR, DTRUE, DE(K), DERR
      end if

      IF (K .EQ. 1)  THEN
!
!  INITIALIZE.
!
        FERMAX = ABS(FERR)
        PFERMX = YE(1)
        DERMAX = ABS(DERR)
        PDERMX = YE(1)
        FDIFMX = ABS(FE2(1) - FE(1))
        PDIFMX = YE(1)
      ELSE
!
!  SELECT.
!
        FERR = ABS(FERR)
        IF (FERR .GT. FERMAX)  THEN
          FERMAX = FERR
          PFERMX = YE(K)
        ENDIF
        DERR = ABS(DERR)
        IF (DERR .GT. DERMAX)  THEN
          DERMAX = DERR
          PDERMX = YE(K)
        ENDIF
        FDIFF = ABS(FE2(K) - FE(K))
        IF (FDIFF .GT. FDIFMX)  THEN
          FDIFMX = FDIFF
          PDIFMX = YE(K)
        ENDIF
      ENDIF

    end do

    FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
    FAILE = FDIFMX .NE. ZERO
    FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.20) .OR. (IER2.NE.IERR)

    IF (FAILOC .AND. (KPRINT.GE.2)) then
      WRITE (LOUT, 2003)  'I', I, 'X', X(I)
    end if

    IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) ) then
      WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
    end if

    IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL

    IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) ) then
      WRITE (LOUT, 2005)  FDIFMX, PDIFMX
    end if

    IF ((IERR.NE.20) .AND. (KPRINT.GE.2)) then
      WRITE (LOUT, 2006)  'D', IERR, 20
    end if

    IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2)) then
      WRITE (LOUT, 2006)  'E', IER2, IERR
    end if

    GO TO 39

   35    CONTINUE

     FAILOC = .TRUE.
     IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR

   39    CONTINUE
     IF (FAILOC)  NERR = NERR + 1
     FAIL = FAIL .OR. FAILOC

  end do

  IF (KPRINT .GE. 2)  THEN
    IF (NERR .GT. 0)  THEN
      WRITE (LOUT, 3001)  NERR, 'I'
    ELSE
      WRITE (LOUT, 4000)  'I'
    ENDIF
  ENDIF

  RETURN

 1000 FORMAT ('1'//10X,'TEST DPCHFE AND DPCHFD')
 1001 FORMAT (//10X,'DEVPCK RESULTS'/10X,'--------------')
 2000 FORMAT (//20X,'DPCHFD INCREMENT TEST -- INCFD = ',I2 &
      /15X,'ON ',A1,'-LINE ',I2,',  ',A1,' =',F8.4, &
         '  --  IERR =',I3)
 2001 FORMAT ( /3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF', &
         13X,'D',8X,'DE',9X,'DIFF')
 2002 FORMAT (F7.2,2(2X,2F10.5,1P,D15.5,0P))
 2003 FORMAT (/' ***** DPCHFD AND/OR DPCHFE FAILED ON ',A1,'-LINE ',I1, &
              ',  ',A1,' =',F8.4)
 2004 FORMAT (/19X,'  MAXIMUM ERROR IN FUNCTION =',1P, &
                1P,D13.5,0P,' (AT',F6.2,'),' &
     /33X,    'IN DERIVATIVE =',1P,D13.5,0P,' (AT',F6.2,').' )
 2005 FORMAT ( '  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =', &
                1P,D13.5,0P,' (AT',F6.2,').' )
 2006 FORMAT (/'  DPCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
 2014 FORMAT ('  *** BOTH SHOULD BE .LE. TOL =',1P,D12.5,' ***')
 3000 FORMAT (//' ***** ERROR ***** DPCHFD RETURNED IERR =',I5//)
 3001 FORMAT (//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON',I2, &
             1X, A1,'-LINES.'//)
 4000 FORMAT (/' DPCHFD AND DPCHFE OK ON ',A1,'-LINES.')
END
SUBROUTINE DFDTRU ( X, F, D )

!*****************************************************************************80
!
!! DFDTRU computes exact function values for DEVCHK.
!
!  Discussion:
!
!    F(X) = X * ( X + 1 ) * ( X - 2 ).
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F, D, the values of the function and
!    derivative at X.
!
  implicit none

  real ( kind = 8 ) D
  real ( kind = 8 ) F
  real ( kind = 8 ) FACT1
  real ( kind = 8 ) FACT2
  real ( kind = 8 ) X
  real ( kind = 8 ) XX

  XX = X
  FACT1 = XX + 1.0D+00
  FACT2 = XX - 2.0D+00
  F = XX * FACT1 * FACT2
  D = FACT1 * FACT2 + XX * ( FACT1 + FACT2 )

  RETURN
END
SUBROUTINE DINTRV ( XT, LXT, X, ILO, ILEFT, MFLAG )

!*****************************************************************************80
!
!! DINTRV finds the interval containing a point.
!
!  Discussion:
!
!    Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
!    such that XT(ILEFT) .LE. X where XT(*) is a subdivision of the X interval.
!
!    DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!    LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!    the X interval.  Precisely,
!
!                      X .LT. XT(1)                1         -1
!       if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!         XT(LXT) .LE. X                         LXT        1,
!
!    That is, when multiplicities are present in the break point
!    to the left of X, the largest index is taken for ILEFT.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Daniel Amos
!
!  Reference:
!
!    Carl de Boor, 
!    Package for calculating with B-splines,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 3, June 1977, pages 441-472.
!
!  Parameters:
!
!       Input      XT,X are double precision
!        XT      - XT is a knot or break point vector of length LXT
!        LXT     - length of the XT vector
!        X       - argument
!        ILO     - an initialization parameter which must be set
!                  to 1 the first time the spline array XT is
!                  processed by DINTRV.
!
!       Output
!        ILO     - ILO contains information for efficient process-
!                  ing after the initial call and ILO must not be
!                  changed by the user.  Distinct splines require
!                  distinct ILO parameters.
!        ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!        MFLAG   - signals when X lies out of bounds
!
  implicit none

  integer ( kind = 4 ) IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
  real ( kind = 8 ) X, XT(*)

  IHI = ILO + 1
  IF (IHI.LT.LXT) GO TO 10
  IF (X.GE.XT(LXT)) GO TO 110
  IF (LXT.LE.1) GO TO 90
  ILO = LXT - 1
  IHI = LXT

   10 continue

  IF (X.GE.XT(IHI)) GO TO 40
  IF (X.GE.XT(ILO)) GO TO 100
!
!  NOW X .LT. XT(IHI) . FIND LOWER BOUND
!
  ISTEP = 1
   20 continue
  IHI = ILO
  ILO = IHI - ISTEP
  IF (ILO.LE.1) GO TO 30
  IF (X.GE.XT(ILO)) GO TO 70
  ISTEP = ISTEP*2
  GO TO 20
   30 continue
  ILO = 1
  IF (X.LT.XT(1)) GO TO 90
  GO TO 70
!
!  NOW X .GE. XT(ILO) . FIND UPPER BOUND
!
   40 continue
  ISTEP = 1
   50 continue
  ILO = IHI
  IHI = ILO + ISTEP
  IF (IHI.GE.LXT) GO TO 60
  IF (X.LT.XT(IHI)) GO TO 70
  ISTEP = ISTEP*2
  GO TO 50
   60 continue
  IF (X.GE.XT(LXT)) GO TO 110
  IHI = LXT
!
!  NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
!
   70 continue
  MIDDLE = (ILO+IHI)/2
  IF (MIDDLE.EQ.ILO) GO TO 100
!
!  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
!
  IF (X.LT.XT(MIDDLE)) GO TO 80
  ILO = MIDDLE
  GO TO 70
   80 continue
  IHI = MIDDLE
  GO TO 70
!
!  SET OUTPUT AND RETURN
!
   90 continue
  MFLAG = -1
  ILEFT = 1
  RETURN
  100 continue
  MFLAG = 0
  ILEFT = ILO
  RETURN
  110 continue
  MFLAG = 1
  ILEFT = LXT

  RETURN
END
SUBROUTINE DPCHBS ( N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF, NDIM, KORD, &
  IERR )

!*****************************************************************************80
!
!! DPCHBS: Piecewise cubic Hermite to B-Spline converter.
!
!  Discussion:
!
!    DPCHBS computes the B-spline representation of the PCH function
!    determined by N,X,F,D.  To be compatible with the rest of PCHIP,
!    DPCHBS includes INCFD, the increment between successive values of
!    the F- and D-arrays.
!
!    The output is the B-representation for the function:  NKNOTS, T,
!    BCOEF, NDIM, KORD.
!
!    Restrictions/assumptions:
!     1. N.GE.2 .  (not checked)
!     2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
!     3. INCFD.GT.0 .  (not checked)
!     4. KNOTYP.LE.2 .  (error return if not)
!    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
!    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
!
!    * Indicates this applies only if KNOTYP.LT.0 .
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch,
!    Representations for Parametric Cubic Splines,
!    Computer Aided Geometric Design,
!    Volume 6, 1989, pages 79-82.
!
! *Usage:
!
!        integer ( kind = 4 )  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
!        PARAMETER  (INCFD = ...)
!        real ( kind = 8 )  X(nmax), F(INCFD,nmax), D(INCFD,nmax),
!       *      T(2*nmax+4), BCOEF(2*nmax)
!
!        CALL DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
!       *             NDIM, KORD, IERR)
!
!  Parameters:
!
!     N:IN  is the number of data points, N.ge.2 .  (not checked)
!
!     X:IN  is the real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.   (not checked)
!           nmax, the dimension of X, must be .ge.N.
!
!     F:IN  is the real array of dependent variable values.
!           F(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of F, must be .ge.N.
!
!     D:IN  is the real array of derivative values at the data points.
!           D(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of D, must be .ge.N.
!
!     INCFD:IN  is the increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           It may have the value 1 for one-dimensional applications,
!           in which case F and D may be singly-subscripted arrays.
!
!     KNOTYP:IN  is a flag to control the knot sequence.
!           The knot sequence T is normally computed from X by putting
!           a double knot at each X and setting the end knot pairs
!           according to the value of KNOTYP:
!              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
!              KNOTYP = 1:  Replicate lengths of extreme subintervals:
!                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
!                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
!              KNOTYP = 2:  Periodic placement of boundary knots:
!                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
!                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
!              Here M=NDIM=2*N.
!           If the input value of KNOTYP is negative, however, it is
!           assumed that NKNOTS and T were set in a previous call.
!           This option is provided for improved efficiency when used
!           in a parametric setting.
!
!     NKNOTS:INOUT  is the number of knots.
!           If KNOTYP.GE.0, then NKNOTS will be set to NDIM+4.
!           If KNOTYP.LT.0, then NKNOTS is an input variable, and an
!              error return will be taken if it is not equal to NDIM+4.
!
!     T:INOUT  is the array of 2*N+4 knots for the B-representation.
!           If KNOTYP.GE.0, T will be returned by DPCHBS with the
!              interior double knots equal to the X-values and the
!              boundary knots set as indicated above.
!           If KNOTYP.LT.0, it is assumed that T was set by a
!              previous call to DPCHBS.  (This routine does **not**
!              verify that T forms a legitimate knot sequence.)
!
!     BCOEF:OUT  is the array of 2*N B-spline coefficients.
!
!     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
!
!     KORD:OUT  is the order of the B-spline.  (Set to 4.)
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -4  if KNOTYP.GT.2 .
!              IERR = -5  if KNOTYP.LT.0 and NKNOTS.NE.(2*N+4).
!
  implicit none

  integer ( kind = 4 )  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*), T(*), BCOEF(*)
  integer ( kind = 4 )  K, KK
  real ( kind = 8 )  DOV3, HNEW, HOLD
  CHARACTER*8  LIBNAM, SUBNAM
!
!  Initialize.
!
  NDIM = 2*N
  KORD = 4
  IERR = 0
  LIBNAM = 'SLATEC'
  SUBNAM = 'DPCHBS'
!
!  Check argument validity.  Set up knot sequence if OK.
!
  IF ( KNOTYP.GT.2 )  THEN
     IERR = -1
     CALL XERMSG (LIBNAM, SUBNAM, 'KNOTYP GREATER THAN 2', IERR, 1)
     RETURN
  end if
  IF ( KNOTYP.LT.0 )  THEN
     IF ( NKNOTS.NE.NDIM+4 )  THEN
        IERR = -2
        CALL XERMSG (LIBNAM, SUBNAM, &
          'KNOTYP.LT.0 AND NKNOTS.NE.(2*N+4)', IERR, 1)
        RETURN
     end if
  ELSE
!
!  Set up knot sequence.
!
     NKNOTS = NDIM + 4
     CALL DPCHKT (N, X, KNOTYP, T)
  end if
!
!  Compute B-spline coefficients.
!
  HNEW = T(3) - T(1)
  DO K = 1, N
     KK = 2*K
     HOLD = HNEW
!
!  The following requires mixed mode arithmetic.
!
     DOV3 = D(1,K)/3
     BCOEF(KK-1) = F(1,K) - HOLD*DOV3
!
!  The following assumes T(2*K+1) = X(K).
!
     HNEW = T(KK+3) - T(KK+1)
     BCOEF(KK) = F(1,K) + HNEW*DOV3
  end do

  RETURN
END
SUBROUTINE DPCHCE ( IC, VC, N, X, H, SLOPE, D, INCFD, IERR )

!*****************************************************************************80
!
!! DPCHCE sets boundary conditions for DPCHIC.
!
!  Discussion:
!
!    Called by DPCHIC to set end derivatives as requested by the user.
!    It must be called after interior derivative values have been set.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, IERR
!        real ( kind = 8 )  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!           ( see prologue to DPCHIC for details. )
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!                    VC(2) need be set only if IC(2) = 2 or 3 .
!
!     N -- (input) number of data points.  (assumes N.GE.2)
!
!     X -- (input) real*8 array of independent variable values.  (the
!           elements of X are assumed to be strictly increasing.)
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real*8 array of derivative values at the data points.
!           The value corresponding to X(I) must be stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!          (output) the value of D at X(1) and/or X(N) is changed, if
!           necessary, to produce the requested boundary conditions.
!           no other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!
  implicit none

  integer ( kind = 4 )  IC(2), N, INCFD, IERR
  real ( kind = 8 )  VC(2), X(*), H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  IBEG, IEND, IERF, INDEX, J, K
  real ( kind = 8 )  HALF, STEMP(3), XTEMP(4)
  SAVE HALF
  real ( kind = 8 )  DPCHDF, DPCHST
!
!  INITIALIZE.
!
  DATA HALF/.5D0/

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
  IF ( ABS(IBEG).GT.N )  IBEG = 0
  IF ( ABS(IEND).GT.N )  IEND = 0
!
!  TREAT BEGINNING BOUNDARY CONDITION.
!
  IF (IBEG == 0)  GO TO 2000
  K = ABS(IBEG)
  IF (K == 1)  THEN
!
!  BOUNDARY VALUE PROVIDED.
!
     D(1,1) = VC(1)
  ELSE IF (K == 2)  THEN
!
!  BOUNDARY SECOND DERIVATIVE PROVIDED.
!
     D(1,1) = HALF*( ( 3.0D+00 *SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
  ELSE IF (K .LT. 5)  THEN
!
!  USE K-POINT DERIVATIVE FORMULA.
!  PICK UP FIRST K POINTS, IN REVERSE ORDER.
!
     DO J = 1, K
        INDEX = K-J+1
!
!  INDEX RUNS FROM K DOWN TO 1.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
     end do

     D(1,1) = DPCHDF (K, XTEMP, STEMP, IERF)

     IF (IERF .NE. 0)  GO TO 5001
  ELSE
!
!  USE 'NOT A KNOT' CONDITION.
!
     D(1,1) = ( 3.0D+00 *(H(1)*SLOPE(2) + H(2)*SLOPE(1)) &
              - 2.0D+00 * (H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
  end if
!
  IF (IBEG .GT. 0)  GO TO 2000
!
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!
  IF (SLOPE(1) == 0.0D+00 )  THEN
     IF (D(1,1) .NE. 0.0D+00 )  THEN
        D(1,1) = 0.0D+00
        IERR = IERR + 1
     end if
  ELSE IF ( DPCHST(D(1,1),SLOPE(1)) .LT. 0.0D+00 )  THEN
     D(1,1) = 0.0D+00
     IERR = IERR + 1
  ELSE IF ( ABS(D(1,1)) .GT. 3.0D+00 *ABS(SLOPE(1)) )  THEN
     D(1,1) = 3.0D+00 *SLOPE(1)
     IERR = IERR + 1
  end if
!
!  TREAT END BOUNDARY CONDITION.
!
 2000 CONTINUE

  IF ( IEND == 0 ) then
    return
  end if

  K = ABS(IEND)
  IF (K == 1)  THEN
!
!  BOUNDARY VALUE PROVIDED.
!
     D(1,N) = VC(2)
  ELSE IF (K == 2)  THEN
!
!  BOUNDARY SECOND DERIVATIVE PROVIDED.
!
     D(1,N) = HALF*( ( 3.0D+00 *SLOPE(N-1) - D(1,N-1)) + &
                                            HALF*VC(2)*H(N-1) )
  ELSE IF (K .LT. 5)  THEN
!
!  USE K-POINT DERIVATIVE FORMULA.  PICK UP LAST K POINTS.
!
     DO J = 1, K
        INDEX = N-K+J
!
!  INDEX RUNS FROM N+1-K UP TO N.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
     end do

     D(1,N) = DPCHDF (K, XTEMP, STEMP, IERF)

     IF (IERF .NE. 0)  GO TO 5001
  ELSE
!
!  USE 'NOT A KNOT' CONDITION.
!
     D(1,N) = ( 3.0D+00 *(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1)) &
              - 2.0D+00 * (H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) ) &
                                                          / H(N-2)
  end if

  IF ( IEND .GT. 0 ) then
    return
  end if
!
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!
  IF (SLOPE(N-1) == 0.0D+00 )  THEN
     IF (D(1,N) .NE. 0.0D+00 )  THEN
        D(1,N) = 0.0D+00
        IERR = IERR + 2
     end if
  ELSE IF ( DPCHST(D(1,N),SLOPE(N-1)) .LT. 0.0D+00 )  THEN
     D(1,N) = 0.0D+00
     IERR = IERR + 2
  ELSE IF ( ABS(D(1,N)) .GT. 3.0D+00 *ABS(SLOPE(N-1)) )  THEN
     D(1,N) = 3.0D+00 *SLOPE(N-1)
     IERR = IERR + 2
  end if

  RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!
!  ERROR RETURN FROM DPCHDF.  THIS CASE SHOULD NEVER OCCUR.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHCE', 'ERROR RETURN FROM DPCHDF', IERR, 1)

  RETURN
END
SUBROUTINE DPCHCI ( N, H, SLOPE, D, INCFD )

!*****************************************************************************80
!
!! DPCHCI sets interior derivatives for DPCHIC.
!
!  Discussion:
!
!    Called by DPCHIC to set derivatives needed to determine a monotone
!    piecewise cubic Hermite interpolant to the data.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by DPCHIM.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N
!        real ( kind = 8 )  H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
!
!  Parameters:
!
!     N -- (input) number of data points.
!           If N=2, simply does linear interpolation.
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (output) real*8 array of derivative values at data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
  implicit none

  integer ( kind = 4 ) incfd

  real ( kind = 8 )  DPCHST
  integer ( kind = 4 )  N
  real ( kind = 8 )  H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  I, NLESS1
  real ( kind = 8 )  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM
  real ( kind = 8 ) HSUMT3, W1, W2

  NLESS1 = N - 1
  DEL1 = SLOPE(1)
!
!  SPECIAL CASE N=2, USE LINEAR INTERPOLATION.
!
  IF ( NLESS1 .le. 1 )  then
    D(1,1) = DEL1
    D(1,N) = DEL1
    return
  end if
!
!  NORMAL CASE  ( 3 <= N ).
!
  DEL2 = SLOPE(2)
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  HSUM = H(1) + H(2)
  W1 = (H(1) + HSUM)/HSUM
  W2 = -H(1)/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2

  IF ( DPCHST(D(1,1),DEL1) .LE. 0.0D+00 )  THEN
    D(1,1) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
    DMAX = 3.0D+00 *DEL1
    IF (ABS(D(1,1)) .GT. ABS(DMAX)) then
      D(1,1) = DMAX
    end if
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO I = 2, NLESS1

    IF ( 2 < I ) then
      HSUM = H(I-1) + H(I)
      DEL1 = DEL2
      DEL2 = SLOPE(I)
    end if
!
!  SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!  USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
    IF ( DPCHST(DEL1,DEL2) .gt. 0.0D+00 ) then
      HSUMT3 = HSUM+HSUM+HSUM
      W1 = (HSUM + H(I-1))/HSUMT3
      W2 = (HSUM + H(I)  )/HSUMT3
      DMAX = MAX( ABS(DEL1), ABS(DEL2) )
      DMIN = MIN( ABS(DEL1), ABS(DEL2) )
      DRAT1 = DEL1/DMAX
      DRAT2 = DEL2/DMAX
      D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
    else
      D(1,I) = 0.0D+00
    end if

  end do
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  W1 = -H(N-1)/HSUM
  W2 = (H(N-1) + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2

  IF ( DPCHST(D(1,N),DEL2) .LE. 0.0D+00 )  THEN
    D(1,N) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
    DMAX = 3.0D+00 *DEL2
    IF (ABS(D(1,N)) .GT. ABS(DMAX)) then
      D(1,N) = DMAX
    end if
  end if

  RETURN
END
SUBROUTINE DPCHCM ( N, X, F, D, INCFD, SKIP, ISMON, IERR )

!*****************************************************************************80
!
!! DPCHCM checks a cubic Hermite function for monotonicity.
!
!  Discussion:
!
!    Checks the piecewise cubic Hermite function defined by N, X, F, D
!    for monotonicity.
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
! *Usage:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, ISMON(N), IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        CALL  DPCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
!
!  Parameters:
!
!     N:IN  is the number of data points.  (Error return if N.LT.2 .)
!
!     X:IN  is a real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F:IN  is a real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D:IN  is a real*8 array of derivative values.  D(1+(I-1)*INCFD) is
!           is the value corresponding to X(I).
!
!     INCFD:IN  is the increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP:INOUT  is a logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed.
!           SKIP will be set to .TRUE. on normal return.
!
!     ISMON:OUT  is an integer ( kind = 4 ) array indicating on which intervals
!     the PCH function defined by  N, X, F, D  is monotonic.
!           For data interval [X(I),X(I+1)],
!             ISMON(I) = -3  if function is probably decreasing;
!             ISMON(I) = -1  if function is strictly decreasing;
!             ISMON(I) =  0  if function is constant;
!             ISMON(I) =  1  if function is strictly increasing;
!             ISMON(I) =  2  if function is non-monotonic;
!             ISMON(I) =  3  if function is probably increasing.
!                If ABS(ISMON)=3, this means that the D-values are near
!                the boundary of the monotonicity region.  A small
!                increase produces non-monotonicity; decrease, strict
!                monotonicity.
!           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
!              the entire function is monotonic on [X(1),X(N)].
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!          (The ISMON-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 ) N

  integer ( kind = 4 ) IERR
  integer ( kind = 4 ) INCFD
  integer ( kind = 4 ) ISMON(N)
  real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N)
  LOGICAL  SKIP
  integer ( kind = 4 ) I, NSEG
  real ( kind = 8 )  DELTA
  integer ( kind = 4 ) DCHFCM
!
!  CHECK ARGUMENTS.
!
  IF ( .not. skip ) then

    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do

  end if

  SKIP = .TRUE.
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  NSEG = N - 1

  DO I = 1, NSEG

     DELTA = (F(1,I+1)-F(1,I))/(X(I+1)-X(I))

     ISMON(I) = DCHFCM (D(1,I), D(1,I+1), DELTA)

     IF ( I == 1 ) THEN
        ISMON(N) = ISMON(1)
     ELSE
!
!  Need to figure out cumulative monotonicity from following
!  "multiplication table":
!
!                    +        I S M O N (I)
!                     +  -3  -1   0   1   3   2
!                      +------------------------+
!               I   -3 I -3  -3  -3   2   2   2 I
!               S   -1 I -3  -1  -1   2   2   2 I
!               M    0 I -3  -1   0   1   3   2 I
!               O    1 I  2   2   1   1   3   2 I
!               N    3 I  2   2   3   3   3   2 I
!              (N)   2 I  2   2   2   2   2   2 I
!                      +------------------------+
!
!  Note that the 2 row and column are out of order so as not
!  to obscure the symmetry in the rest of the table.
!
!  No change needed if equal or constant on this interval or
!  already declared nonmonotonic.
!
        IF ( (ISMON(I).NE.ISMON(N)) .AND. (ISMON(I).NE.0) &
                                    .AND. (ISMON(N).NE.2) )  THEN
           IF ( (ISMON(I) == 2) .OR. (ISMON(N) == 0) )  THEN
              ISMON(N) =  ISMON(I)
           ELSE IF (ISMON(I)*ISMON(N) .LT. 0)  THEN
!
!  This interval has opposite sense from curve so far.
!
              ISMON(N) = 2
           ELSE
!
!  At this point, both are nonzero with same sign, and
!  we have already eliminated case both +-1.
!
              ISMON(N) = ISIGN (3, ISMON(N))
           end if
        end if
     end if

  end do

  IERR = 0
  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHCM', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  INCFD.LT.1 RETURN.
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHCM', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN

 5003 CONTINUE
!
!  X-ARRAY NOT STRICTLY INCREASING.
!
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHCM', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)

  RETURN
END
SUBROUTINE DPCHCS ( SWITCH, N, H, SLOPE, D, INCFD, IERR )

!*****************************************************************************80
!
!! DPCHCS adjusts derivative values for DPCHIC.
!
!  Discussion:
!
!    Called by DPCHIC to adjust the values of D in the vicinity of a
!    switch in direction of monotonicity, to produce a more "visually
!    pleasing" curve than that given by DPCHIM.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        real ( kind = 8 )  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!  Parameters:
!
!     SWITCH -- (input) indicates the amount of control desired over
!           local excursions from data.
!
!     N -- (input) number of data points.  (assumes N.GT.2 .)
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real*8 array of derivative values at the data points,
!           as determined by DPCHCI.
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually
!           pleasing" curve.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.  should be zero.
!           If negative, trouble in DPCHSW.  (should never happen.)
!
  implicit none

  integer incfd

  integer ( kind = 4 )  N, IERR
  real ( kind = 8 )  SWITCH, H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  I, INDX, K, NLESS1
  real ( kind = 8 )  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE
  real ( kind = 8 ) SLMAX, WTAVE(2)

  SAVE ONE, FUDGE
  real ( kind = 8 )  DPCHST
!
!  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!
  real ( kind = 8 )  DPCHSD, S1, S2, H1, H2
  DPCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!
!  INITIALIZE.
!
  DATA ONE/1.D0/
  DATA FUDGE /4.D0/

  IERR = 0
  NLESS1 = N - 1
!
!  LOOP OVER SEGMENTS.
!
  DO I = 2, NLESS1

     IF ( DPCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900

  100    CONTINUE
!
!  SLOPE SWITCHES MONOTONICITY AT I-TH POINT.
!
!  DO NOT CHANGE D IF 'UP-DOWN-UP'.
!
        IF (I .GT. 2)  THEN
           IF ( DPCHST(SLOPE(I-2),SLOPE(I)) .GT. 0.0D+00 )  GO TO 900

        end if
        IF (I .LT. NLESS1)  THEN
           IF ( DPCHST(SLOPE(I+1),SLOPE(I-1)) .GT. 0.0D+00 )  GO TO 900
        end if
!
!  COMPUTE PROVISIONAL VALUE FOR D(1,I).
!
        DEXT = DPCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!
!  DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!
        IF ( DPCHST(DEXT, SLOPE(I-1)) )  200, 900, 250

  200       CONTINUE
!
!  DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!  EXTREMUM IS IN (X(I-1),X(I)).
!
           K = I-1
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
!
           WTAVE(2) = DEXT

           IF (K .GT. 1) then
              WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
           end if

           GO TO 400

  250       CONTINUE
!
!  DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS.  EXTREMUM IS IN (X(I),X(I+1)).
!
           K = I
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
!
           WTAVE(1) = DEXT

           IF (K .LT. NLESS1) then
              WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
           end if

           GO TO 400

  300    CONTINUE
!
!  AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO.
!  CHECK FOR FLAT-TOPPED PEAK.
!
        IF (I == NLESS1)  GO TO 900
        IF ( DPCHST(SLOPE(I-1), SLOPE(I+1)) .GE. 0.0D+00 )  GO TO 900
!
!  WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
!
        K = I
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
!
        WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
        WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))

  400    CONTINUE
!
!  AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!  ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!  WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K), IF K.GT.1
!  WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1), IF K.LT.N-1
!
     SLMAX = ABS(SLOPE(K))

     IF (K .GT. 1) then
       SLMAX = MAX( SLMAX, ABS(SLOPE(K-1)) )
     end if

     IF (K.LT.NLESS1) SLMAX = MAX( SLMAX, ABS(SLOPE(K+1)) )

     IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
     DEL(2) = SLOPE(K) / SLMAX
     IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX

     IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
!
!  NORMAL CASE.  EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
!
        FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
        D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
        FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
        D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
     ELSE
!
!  SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!  K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
!
        FACT = FUDGE* ABS(DEL(2))
        D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
!
!  NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!  I-K+1 = 2 IF K=I-1(=1).
!
     end if
!
!  ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!
     IF (SWITCH .LE. 0.0D+00 )  GO TO 900

     DFLOC = H(K)*ABS(SLOPE(K))
     IF (K .GT. 1)    DFLOC = MAX( DFLOC, H(K-1)*ABS(SLOPE(K-1)) )
     IF (K.LT.NLESS1) DFLOC = MAX( DFLOC, H(K+1)*ABS(SLOPE(K+1)) )
     DFMX = SWITCH*DFLOC
     INDX = I-K+1
!
!  INDX = 1 IF K=I, 2 IF K=I-1.
!
     CALL DPCHSW(DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)

     IF (IERR .NE. 0) then
       RETURN
     end if

  900    CONTINUE

  end do

  RETURN
END
FUNCTION DPCHDF ( K, X, S, IERR )

!*****************************************************************************80
!
!! DPCHDF computes divided differences for DPCHCE and DPCHSP.
!
!  Discussion:
!
!    Uses a divided difference formulation to compute a K-point approx-
!    imation to the derivative at X(K) based on the data in X and S.
!
!    Called by DPCHCE and DPCHSP to compute 3- and 4-point boundary
!    derivative approximations.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K.LT.2 .
!        DPCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
  implicit none

  integer ( kind = 4 )  K

  real ( kind = 8 ) dpchdf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(k)
  real ( kind = 8 ) value
  real ( kind = 8 )  X(K)

  ierr = 0
!
!  CHECK FOR LEGAL VALUE OF K.
!
  IF ( K .LT. 3 ) then
    IERR = -1
    CALL XERMSG ('SLATEC', 'DPCHDF', 'K LESS THAN THREE', IERR, 1)
    DPCHDF = 0.0D+00
    return
  end if
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
  DO J = 2, K-1
    DO I = 1, K-J
      S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    end do
  end do
!
!  EVALUATE DERIVATIVE AT X(K).
!
  VALUE = S(1)
  DO I = 2, K-1
    VALUE = S(I) + VALUE*(X(K)-X(I))
  end do

  DPCHDF = VALUE

  RETURN
END
SUBROUTINE DPCHEV ( N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR )

!*****************************************************************************80
!
!! DPCHEV: value and derivative of spline or cubic Hermite at many points.
!
!  Discussion:
!
!    Evaluates the function and first derivative of the cubic Hermite
!    or spline function defined by  N, X, F, D, at the array of points XVAL.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence: CALL  DPCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
!
!     integer ( kind = 4 )  N, NVAL, IERR
!     real ( kind = 8 )  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real ( kind = 8 ) array of independent variable 
!           values.  The elements of X must be strictly increasing:
!             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
!
!     F -- (input) real ( kind = 8 ) array of function values.  F(I) is
!           the value corresponding to X(I).
!
!     D -- (input) real ( kind = 8 ) array of derivative values.  
!          D(I) is the value corresponding to X(I).
!
!  NVAL -- (input) number of points at which the functions are to be
!           evaluated. ( Error return if NVAL.LT.1 )
!
!  XVAL -- (input) real ( kind = 8 ) array of points at which the 
!          functions are to be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XVAL are increasing relative to X;
!              that is,   XVAL(J) .GE. X(I)
!              implies    XVAL(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XVAL are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  FVAL -- (output) real ( kind = 8 ) array of values of the cubic 
!          Hermite function defined by  N, X, F, D  at the points  XVAL.
!
!  DVAL -- (output) real ( kind = 8 ) array of values of the 
!          first derivative of the same function at the points  XVAL.
!
!  IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NVAL.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine DCHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
  implicit none

  integer ( kind = 4 )  N, NVAL, IERR
  real ( kind = 8 )  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)

  integer ( kind = 4 ) INCFD
  LOGICAL SKIP

  skip = .true.
  incfd = 1

  CALL DPCHFD ( N, X, F, D, INCFD, SKIP, NVAL, XVAL, FVAL, DVAL, IERR )

  RETURN
END
SUBROUTINE DPCHEZ ( N, X, F, D, SPLINE, WK, LWK, IERR )

!*****************************************************************************80
!
!! DPCHEZ sets up a spline or cubic Hermite interpolant.
!
!  Discussion:
!
!    Sets derivatives for spline (two continuous derivatives) or
!    Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data. This is an easy to use driver for the routines
!    by F. N. Fritsch. Various boundary
!    conditions are set to default values by DPCHEZ. Many other choices
!    are available in the subroutines PCHIC, DPCHIM and DPCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence:   CALL  DPCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
!
!     integer ( kind = 4 )  N, IERR,  LWK
!     real ( kind = 8 )  X(N), F(N), D(N), WK(*)
!     LOGICAL SPLINE
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real ( kind = 8 ) array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real ( kind = 8 ) array of dependent variable values to be 
!    interpolated.  F(I) is value corresponding to X(I).
!
!     D -- (output) real ( kind = 8 ) array of derivative values at 
!     the data points.
!
!     SPLINE -- (input) logical variable to specify if the interpolant
!           is to be a spline with two continuous derivaties
!           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
!           continuous derivative (set SPLINE=.FALSE.).
!        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
!           default "not-a-knot" boundary condition, with a continuous
!           third derivative at X(2) and X(N-1).
!              If SPLINE=.FALSE. the interpolating Hermite cubic will be
!           monotone if the input data is monotone. Boundary conditions
!           computed from the derivative of a local quadratic unless thi
!           alters monotonicity.
!
!     WK -- (scratch) real ( kind = 8 ) work array, which must be declared by 
!     the calling program to be at least 2*N if SPLINE is .TRUE. and not used
!           otherwise.
!
!     LWK -- (input) length of work array WK. (Error return if
!           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means tha
!                 IERR switches in the direction of monotonicity were de
!                 When SPLINE=.FALSE.,  DPCHEZ guarantees that if the inp
!                 data is monotone, the interpolant will be too. This wa
!                 is to alert you to the fact that the input data was no
!                 monotone.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  integer ( kind = 4 ) IERR
  real ( kind = 8 )  X(N), F(N), D(N), WK(LWK)
  LOGICAL SPLINE
  integer ( kind = 4 ) IC(2), INCFD
  real ( kind = 8 )  VC(2)

  ic(1) = 0
  ic(2) = 0
  incfd = 1

  IF ( SPLINE ) THEN
    CALL  DPCHSP ( IC, VC, N, X, F, D, INCFD, WK, LWK, IERR )
  ELSE
    CALL  DPCHIM ( N, X, F, D, INCFD, IERR )
  end if

  RETURN
END
SUBROUTINE DPCHFD ( N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR )

!*****************************************************************************80
!
!! DPCHFD evaluates a piecewise cubic Hermite and derivative at many points.
!
!  Discussion:
!
!    Evaluate a piecewise cubic Hermite function and its first
!    derivative at an array of points.  May be used by itself
!    for Hermite interpolation, or as an evaluator for DPCHIM
!    or DPCHIC.
!
!    Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!    gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!    If only function values are required, use DPCHFE, instead.
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, NE, IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
!                          DE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine DCHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, NE, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*), DE(*)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then

    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do

    skip = .true.

  end if
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ( NE.LT.1 )  GO TO 5004
  IERR = 0
!
!  LOOP OVER INTERVALS.
!  INTERVAL INDEX IS  IL = IR-1.
!  INTERVAL IS X(IL).LE.X.LT.X(IR).
!
  JFIRST = 1
  IR = 2

   10 CONTINUE
!
!  SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
     IF (JFIRST .GT. NE) then
       return
     end if
!
!  LOCATE ALL POINTS IN INTERVAL.
!
     DO J = JFIRST, NE
        IF (XE(J) .GE. X(IR))  GO TO 30
     end do

     J = NE + 1
     GO TO 40
!
!  HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     IF (IR == N)  J = NE + 1

   40    CONTINUE
     NJ = J - JFIRST
!
!  SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
     IF (NJ == 0)  GO TO 50
!
!  EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
    CALL DCHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
               NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)

     IF (IERC .LT. 0)  GO TO 5005

     IF (NEXT(2) == 0)  GO TO 42
!
!  IF (NEXT(2) .GT. 0)  THEN IN THE CURRENT SET OF XE-POINTS, 
!  THERE ARE NEXT(2) TO THE RIGHT OF X(IR).
!
        IF (IR .LT. N)  GO TO 41
!
!  IF ( IR == N ) THEN THESE ARE ACTUALLY EXTRAPOLATION POINTS.
!
           IERR = IERR + NEXT(2)
           GO TO 42
   41       CONTINUE
!
!  ELSE WE SHOULD NEVER HAVE GOTTEN HERE.
!
           GO TO 5005

   42    CONTINUE

     IF (NEXT(1) == 0)  GO TO 49
!
!  IF (NEXT(1) .GT. 0) THEN IN THE CURRENT SET OF XE-POINTS, 
!  THERE ARE NEXT(1) TO THE LEFT OF X(IR-1).
!
        IF (IR .GT. 2)  GO TO 43
!
!  IF (IR == 2) THEN THESE ARE ACTUALLY EXTRAPOLATION POINTS.
!
           IERR = IERR + NEXT(1)
           GO TO 49
   43       CONTINUE
!
!  ELSE XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!  EVALUATION INTERVAL.
!
!  FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
!
           DO I = JFIRST, J-1
              IF (XE(I) .LT. X(IR-1))  GO TO 45
           end do
!
!  CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR IN DCHFDV.
!
           GO TO 5005

   45          CONTINUE
!
!  RESET J.  (THIS WILL BE THE NEW JFIRST.)
!
           J = I
!
!  NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
!
           DO I = 1, IR-1
              IF (XE(J) .LT. X(I)) GO TO 47
           end do
!
!  CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!
!  AT THIS POINT, EITHER  XE(J) .LT. X(1)
!  OR      X(I-1) .LE. XE(J) .LT. X(I) .
!  RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!  CYCLING.
!
           IR = MAX(1, I-1)
!           end if
!        end if
   49    CONTINUE

     JFIRST = J

   50 CONTINUE
  IR = IR + 1
  IF (IR .LE. N)  GO TO 10

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHFD', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHFD', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHFD', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN

 5004 CONTINUE
!     NE.LT.1 RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'DPCHFD', &
    'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM DCHFDV.  THIS CASE SHOULD NEVER OCCUR.
  IERR = -5
  CALL XERMSG ('SLATEC', 'DPCHFD', &
    'ERROR RETURN FROM DCHFDV -- FATAL', IERR, 2)
  RETURN
END
SUBROUTINE DPCHFE ( N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR )

!*****************************************************************************80
!
!! DPCHFE evaluates a piecewise cubic Hermite function at many points.
!
!  Discussion:
!
!    DPCHFE evaluates a piecewise cubic Hermite function at an array of
!    points.  It may be used by itself for Hermite interpolation,
!    or as an evaluator for DPCHIM or DPCHIC.
!
!    Evaluates the cubic Hermite function defined by  N, X, F, D  at
!    the points  XE(J), J=1(1)NE.
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, NE, IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!             (The FE-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 ) incfd

  integer ( kind = 4 )  N, NE, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then

    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do

  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ( NE.LT.1 )  GO TO 5004
  IERR = 0
!
!  LOOP OVER INTERVALS.
!  INTERVAL INDEX IS  IL = IR-1  .
!  INTERVAL IS X(IL).LE.X.LT.X(IR) .
!
  JFIRST = 1
  IR = 2
   10 CONTINUE
!
!  SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
     IF (JFIRST .GT. NE) then
       return
     end if
!
!  LOCATE ALL POINTS IN INTERVAL.
!
     DO J = JFIRST, NE
        IF (XE(J) .GE. X(IR))  GO TO 30
     end do

     J = NE + 1
     GO TO 40
!
!  HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     IF (IR == N)  J = NE + 1

   40    CONTINUE
     NJ = J - JFIRST
!
!  SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
     IF (NJ == 0)  GO TO 50
!
!  EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
    CALL DCHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
      NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)

     IF (IERC .LT. 0)  GO TO 5005

     IF (NEXT(2) == 0)  GO TO 42
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
        IF (IR .LT. N)  GO TO 41
!           IF (IR == N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(2)
           GO TO 42
   41       CONTINUE
!
!  ELSE WE SHOULD NEVER HAVE GOTTEN HERE.
!
           GO TO 5005
!           end if
!        end if
   42    CONTINUE
!
     IF (NEXT(1) == 0)  GO TO 49
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
        IF (IR .GT. 2)  GO TO 43
!           IF (IR == 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(1)
           GO TO 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!  FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
!
           DO I = JFIRST, J-1
              IF (XE(I) .LT. X(IR-1))  GO TO 45
           end do
!
!  CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR IN DCHFEV.
!
           GO TO 5005

   45          CONTINUE
!
!  RESET J.  THIS WILL BE THE NEW JFIRST.
!
           J = I
!
!  NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
!
           DO I = 1, IR-1
              IF (XE(J) .LT. X(I)) GO TO 47
           end do
!
!  CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!
!  AT THIS POINT, EITHER  XE(J) .LT. X(1)
!  OR      X(I-1) .LE. XE(J) .LT. X(I) .
!  RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE CYCLING.
!
           IR = MAX(1, I-1)
!           end if
!        end if
   49    CONTINUE

     JFIRST = J

   50 CONTINUE
  IR = IR + 1
  IF (IR .LE. N)  GO TO 10

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHFE', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  INCFD.LT.1 RETURN.
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHFE', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHFE', 'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'DPCHFE', &
    'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM DCHFEV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -5
  CALL XERMSG ('SLATEC', 'DPCHFE', &
    'ERROR RETURN FROM DCHFEV -- FATAL', IERR, 2)
  RETURN
END
FUNCTION DPCHIA ( N, X, F, D, INCFD, SKIP, A, B, IERR )

!*****************************************************************************80
!
!! DPCHIA evaluates the definite integral of a piecewise cubic Hermite function.
!
!  Discussion:
!
!    Evaluate the definite integral of a piecewise cubic
!    Hermite function over an arbitrary interval.
!
!    Evaluates the definite integral of the cubic Hermite function
!    defined by  N, X, F, D  over the interval [A, B].
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N), A, B
!        real ( kind = 8 )  VALUE, DPCHIA
!        LOGICAL  SKIP
!
!        VALUE = DPCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR.GE.0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -4  in case of an error return from DPCHID (which
!                         should never occur).
!
  implicit none

  real ( kind = 8 ) dpchia
  integer ( kind = 4 )  N, INCFD, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*), A, B
  LOGICAL  SKIP

  integer ( kind = 4 )  I, IA, IB, IERD, IL, IR
  real ( kind = 8 )  VALUE, XA, XB
  real ( kind = 8 )  DCHFIE, DPCHID

  VALUE = 0.0D+00
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0

  IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) ) then
    IERR = IERR + 1
  end if

  IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) ) then
    IERR = IERR + 2
  end if
!
!  COMPUTE INTEGRAL VALUE.
!
  IF (A .NE. B)  THEN
     XA = MIN (A, B)
     XB = MAX (A, B)
     IF (XB .LE. X(2))  THEN
!
!  INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!
        VALUE = DCHFIE (X(1),X(2), F(1,1),F(1,2), &
                                   D(1,1),D(1,2), A, B)

     ELSE IF (XA .GE. X(N-1))  THEN
!
!  INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!
        VALUE = DCHFIE(X(N-1),X(N), F(1,N-1),F(1,N), &
                                  D(1,N-1),D(1,N), A, B)
     ELSE
!
!  NORMAL CASE: XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
!  LOCATE IA AND IB SUCH THAT
!  X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
!
        IA = 1
        DO I = 1, N-1
           IF (XA .GT. X(I))  IA = I + 1
        end do
!
!  IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
!  IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
!
        IB = N
        DO I = N, IA, -1
           IF (XB .LT. X(I))  IB = I - 1
        end do
!
!  IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
!  IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
!
!  COMPUTE THE INTEGRAL.
!
        IF (IB .LT. IA)  THEN
!
!  THIS MEANS IB = IA-1 AND (A,B) IS A SUBSET OF (X(IB),X(IA)).
!
           VALUE = DCHFIE (X(IB),X(IA), F(1,IB),F(1,IA), &
                                      D(1,IB),D(1,IA), A, B)

        ELSE
!
!  FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
!  (Case (IB == IA) is taken care of by initialization
!  of VALUE to ZERO.)
!
           IF (IB .GT. IA)  THEN
              VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)

              IF (IERD .LT. 0)  GO TO 5004
           end if
!
!  THEN ADD ON INTEGRAL OVER (XA,X(IA)).
!
           IF (XA .LT. X(IA))  THEN
              IL = MAX(1, IA-1)
              IR = IL + 1
              VALUE = VALUE + DCHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                       D(1,IL),D(1,IR), XA, X(IA))
           end if
!
!  THEN ADD ON INTEGRAL OVER (X(IB),XB).
!
           IF (XB .GT. X(IB))  THEN
              IR = MIN (IB+1, N)
              IL = IR - 1
              VALUE = VALUE + DCHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                        D(1,IL),D(1,IR), X(IB), XB)
           end if
!
!  FINALLY, ADJUST SIGN IF NECESSARY.
!
           IF (A .GT. B)  VALUE = -VALUE
        end if
     end if
  end if

  DPCHIA = VALUE

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHIA', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHIA', 'INCREMENT LESS THAN ONE', IERR, 1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHIA', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return

 5004 CONTINUE
!
!  TROUBLE IN DPCHID.  (SHOULD NEVER OCCUR.)
!
  IERR = -4
  CALL XERMSG ('SLATEC', 'DPCHIA', 'TROUBLE IN DPCHID', IERR, 1)

  return
END
SUBROUTINE DPCHIC ( IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR )

!*****************************************************************************80
!
!! DPCHIC sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Sets derivatives needed to determine a piecewise monotone 
!    piecewise cubic interpolant to the data given in X and F satisfying the
!    boundary conditions specified by IC and VC.
!
!    User control is available over boundary conditions and/or
!    treatment of points where monotonicity switches direction,
!    using argument SWITCH.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by DPCHFE or DPCHFD.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, NWK, IERR
!        real ( kind = 8 )  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!                          WK(NWK)
!
!        CALL DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  for the default boundary condition (the same as
!                     used by DPCHIM).
!           If IBEG.NE.0, then its sign indicates whether the boundary
!                     derivative is to be adjusted, if necessary, to be
!                     compatible with monotonicity:
!              IBEG.GT.0  if no adjustment is to be performed.
!              IBEG.LT.0  if the derivative is to be adjusted for
!                     monotonicity.
!
!           Allowable values for the magnitude of IBEG are:
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!           IBEG = 5  to set D(1) so that the second derivative is con-
!              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!              This option is somewhat analogous to the "not a knot"
!              boundary condition provided by DPCHSP.
!
!          NOTES (IBEG):
!           1. An error return is taken if ABS(IBEG).GT.5 .
!           2. Only in case  IBEG.LE.0  is it guaranteed that the
!              interpolant will be monotonic in the first interval.
!              If the returned value of D(1) lies between zero and
!              3*SLOPE(1), the interpolant will be monotonic.  This
!              is **NOT** checked if IBEG.GT.0 .
!           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
!              tonicity, a warning error is returned.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES (IEND):
!           1. An error return is taken if ABS(IEND).GT.5 .
!           2. Only in case  IEND.LE.0  is it guaranteed that the
!              interpolant will be monotonic in the last interval.
!              If the returned value of D(1+(N-1)*INCFD) lies between
!              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!              This is **NOT** checked if IEND.GT.0 .
!           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
!              achieve monotonicity, a warning error is returned.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     SWITCH -- (input) indicates desired treatment of points where
!           direction of monotonicity switches:
!           Set SWITCH to zero if interpolant is required to be mono-
!           tonic in each interval, regardless of monotonicity of data.
!             NOTES:
!              1. This will cause D to be set to zero at all switch
!                 points, thus forcing extrema there.
!              2. The result of using this option with the default boun-
!                 dary conditions will be identical to using DPCHIM, but
!                 will generally cost more compute time.
!                 This option is provided only to facilitate comparison
!                 of different switch and/or boundary conditions.
!           Set SWITCH nonzero to use a formula based on the 3-point
!              difference formula in the vicinity of switch points.
!           If SWITCH is positive, the interpolant on each interval
!              containing an extremum is controlled to not deviate from
!              the data by more than SWITCH*DFLOC, where DFLOC is the
!              maximum of the change of F on this interval and its two
!              immediate neighbors.
!           If SWITCH is negative, no such control is to be imposed.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine a monotone cubic
!           Hermite function on each subinterval on which the data
!           are monotonic, except possibly adjacent to switches in
!           monotonicity. The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real*8 array of working storage.  The user may
!           wish to know that the returned values are:
!              WK(I)     = H(I)     = X(I+1) - X(I) ;
!              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!           for  I = 1(1)N-1.
!
!     NWK -- (input) length of work array.
!           (Error return if  NWK.LT.2*(N-1) .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if ABS(IBEG).GT.5 .
!              IERR = -5  if ABS(IEND).GT.5 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK.LT.2*(N-1) .
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) nwk

  integer ( kind = 4 )  IC(2), N, IERR
  real ( kind = 8 )  VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*), WK(NWK)
  integer ( kind = 4 )  I, IBEG, IEND, NLESS1
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO I = 2, N
     IF ( X(I).LE.X(I-1) )  GO TO 5003
  end do

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  IF (ABS(IBEG) .GT. 5)  IERR = IERR - 1
  IF (ABS(IEND) .GT. 5)  IERR = IERR - 2
  IF (IERR .LT. 0)  GO TO 5004
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  NLESS1 = N - 1
  IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
!
!  SET UP H AND SLOPE ARRAYS.
!
  DO I = 1, NLESS1
     WK(I) = X(I+1) - X(I)
     WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
  end do
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF (NLESS1 .GT. 1)  GO TO 1000
  D(1,1) = WK(2)
  D(1,N) = WK(2)
  GO TO 3000
!
!  NORMAL CASE  (N .GE. 3) .
!
 1000 CONTINUE
!
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!
  CALL DPCHCI (N, WK(1), WK(N), D, INCFD)
!
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!
  IF (SWITCH == 0.0D+00 )  GO TO 3000
  CALL DPCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR)
  IF (IERR .NE. 0)  GO TO 5008
!
!  SET END CONDITIONS.
!
 3000 CONTINUE

  IF ( (IBEG == 0) .AND. (IEND == 0) ) then
    return
  end if

  CALL DPCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR)
  IF (IERR .LT. 0)  GO TO 5009

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  N.LT.2 RETURN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHIC', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHIC', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHIC', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  CALL XERMSG ('SLATEC', 'DPCHIC', 'IC OUT OF RANGE', IERR, 1)
  RETURN
!
 5007 CONTINUE
!
!     NWK .LT. 2*(N-1)  RETURN.
!
  IERR = -7
  CALL XERMSG ('SLATEC', 'DPCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
  RETURN
!
 5008 CONTINUE
!     ERROR RETURN FROM DPCHCS.
  IERR = -8
  CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCS', IERR, 1)
  RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM DPCHCE.  THIS CASE SHOULD NEVER OCCUR.
  IERR = -9
  CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCE', IERR, 1)
  RETURN
END
FUNCTION DPCHID ( N, X, F, D, INCFD, SKIP, IA, IB, IERR )

!*****************************************************************************80
!
!! DPCHID: integral of a piecewise cubic Hermite function between data points.
!
!  Discussion:
!
!    Evaluate the definite integral of a piecewise cubic
!    Hermite function over an interval whose endpoints are data
!    points.
!
!    Evaluates the definite integral of the cubic Hermite function
!    defined by  N, X, F, D  over the interval [X(IA), X(IB)].
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IA, IB, IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
!
!     IA,IB -- (input) indices in X-array for the limits of integration.
!           both must be in the range [1,N].  (Error return if not.)
!           No restrictions on their relative values.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IA or IB is out of range.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  real ( kind = 8 ) dpchid
  integer ( kind = 4 )  N, INCFD, IA, IB, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IUP, LOW
  real ( kind = 8 )  H, HALF, SIX, SUM, VALUE
  SAVE HALF, SIX
!
!  INITIALIZE.
!
  DATA HALF/.5D0/, SIX/6.D0/

  VALUE = 0.0D+00
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
  IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
  IERR = 0
!
!  COMPUTE INTEGRAL VALUE.
!
  IF (IA .NE. IB)  THEN
     LOW = MIN(IA, IB)
     IUP = MAX(IA, IB) - 1
     SUM = 0.0D+00
     DO I = LOW, IUP
        H = X(I+1) - X(I)
        SUM = SUM + H*( (F(1,I) + F(1,I+1)) + &
                       (D(1,I) - D(1,I+1))*(H/SIX) )
     end do
     VALUE = HALF * SUM
     IF (IA .GT. IB)  VALUE = -VALUE
  end if

  DPCHID = VALUE

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  N.LT.2 RETURN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHID', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHID', 'INCREMENT LESS THAN ONE', IERR, 1)
  return

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHID', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return

 5004 CONTINUE
!     IA OR IB OUT OF RANGE RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'DPCHID', 'IA OR IB OUT OF RANGE', IERR, 1)

  return
END
SUBROUTINE DPCHIM ( N, X, F, D, INCFD, IERR )

!*****************************************************************************80
!
!! DPCHIM sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.  Boundary values
!    are provided which are compatible with monotonicity.  The
!    interpolant will have an extremum at each point where 
!    monotonicity switches direction.  See DPCHIC if user control
!    is desired over boundary or switch conditions.
!
!    Sets derivatives needed to determine a monotone piecewise cubic
!    Hermite interpolant to the data given in X and F.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  (See DPCHIC if user control of boundary con-
!    ditions is desired.)
!
!    If the data are only piecewise monotonic, the interpolant will
!    have an extremum at each point where monotonicity switches direc-
!    tion.  (See DPCHIC if user control is desired in such cases.)
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by DPCHFE or DPCHFD.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        real ( kind = 8 )  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).  DPCHIM is designed for monotonic data, but it will
!           work for any F-array.  It will force extrema at points where
!           monotonicity switches direction.  If some other treatment of
!           switch points is desired, DPCHIC should be used instead.
! 
!     D -- (output) real*8 array of derivative values at the data
!           points.  If the data are monotonic, these values will
!           determine a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, IERR
  real ( kind = 8 )  X(*), F(INCFD,*), D(INCFD,*)
  integer ( kind = 4 )  I, NLESS1
  real ( kind = 8 )  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE
  real ( kind = 8 ) H1, H2, HSUM, HSUMT3, W1, W2
  real ( kind = 8 )  DPCHST
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO I = 2, N
     IF ( X(I).LE.X(I-1) )  GO TO 5003
  end do
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0
  NLESS1 = N - 1
  H1 = X(2) - X(1)
  DEL1 = (F(1,2) - F(1,1))/H1
  DSAVE = DEL1
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF ( NLESS1 .le. 1) then
    D(1,1) = DEL1
    D(1,N) = DEL1
    return
  end if
!
!  NORMAL CASE  (N .GE. 3).
!
  H2 = X(3) - X(2)
  DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  HSUM = H1 + H2
  W1 = (H1 + HSUM)/HSUM
  W2 = -H1/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2
  IF ( DPCHST(D(1,1),DEL1) .LE. 0.0D+00 )  THEN
     D(1,1) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0D+00 *DEL1
     IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO I = 2, NLESS1

     IF ( 2 < I ) then
       H1 = H2
       H2 = X(I+1) - X(I)
       HSUM = H1 + H2
       DEL1 = DEL2
       DEL2 = (F(1,I+1) - F(1,I))/H2
     end if
!
!  SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = 0.0D+00
     IF ( DPCHST(DEL1,DEL2) )  42, 41, 45
!
!  COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
     IF (DEL2 == 0.0D+00 )  GO TO 50
     IF ( DPCHST(DSAVE,DEL2) .LT. 0.0D+00 )  IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50

   42    CONTINUE
     IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50
!
!  USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H1)/HSUMT3
     W2 = (HSUM + H2)/HSUMT3
     DMAX = MAX( ABS(DEL1), ABS(DEL2) )
     DMIN = MIN( ABS(DEL1), ABS(DEL2) )
     DRAT1 = DEL1/DMAX
     DRAT2 = DEL2/DMAX
     D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)

   50    CONTINUE

  end do
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  W1 = -H2/HSUM
  W2 = (H2 + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2
  IF ( DPCHST(D(1,N),DEL2) .LE. 0.0D+00 )  THEN
     D(1,N) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0D+00 *DEL2
     IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
  end if

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHIM', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHIM', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHIM', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
END
SUBROUTINE DPCHKT ( N, X, KNOTYP, T )

!*****************************************************************************80
!
!! DPCHKT computes the B-spline knot sequence for DPCHBS.
!
!  Discussion:
!
!    Set a knot sequence for the B-spline representation of a PCH
!    function with breakpoints X.  All knots will be at least double.
!
!    Endknots are set as:
!    (1) quadruple knots at endpoints if KNOTYP=0;
!    (2) extrapolate the length of end interval if KNOTYP=1;
!    (3) periodic if KNOTYP=2.
!
!    Restrictions/assumptions:
!    1. N.GE.2 .  (not checked)
!    2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
!    3. 0.LE.KNOTYP.LE.2 .  (Acts like KNOTYP=0 for any other value.)
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!  Input arguments:  N, X, KNOTYP.
!  Output arguments:  T.
!
  implicit none

  integer ( kind = 4 )  N, KNOTYP
  real ( kind = 8 )  X(*), T(*)
  integer ( kind = 4 )  J, K, NDIM
  real ( kind = 8 )  HBEG, HEND

  NDIM = 2 * N
!
!  Set interior knots.
!
  J = 1
  DO K = 1, N
     J = J + 2
     T(J) = X(K)
     T(J+1) = T(J)
  end do
!
!  Assertion:  
!  At this point T(3),...,T(NDIM+2) have been set and J=NDIM+1.
!
!  Set end knots according to KNOTYP.
!
  HBEG = X(2) - X(1)
  HEND = X(N) - X(N-1)
  IF (KNOTYP == 1 )  THEN
!
!  Extrapolate.
!
     T(2) = X(1) - HBEG
     T(NDIM+3) = X(N) + HEND
  ELSE IF ( KNOTYP == 2 )  THEN
!
!  Periodic.
!
     T(2) = X(1) - HEND
     T(NDIM+3) = X(N) + HBEG
  ELSE
!
!  Quadruple end knots.
!
     T(2) = X(1)
     T(NDIM+3) = X(N)
  end if
  T(1) = T(2)
  T(NDIM+4) = T(NDIM+3)

  RETURN
END
SUBROUTINE DPCHQ1 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! DPCHQ1: Test the PCHIP evaluators DCHFDV, DCHFEV, DPCHFD, DPCHFE.
!
!  Discussion:
!
!    This routine carries out three tests of the PCH evaluators:
!    DEVCHK tests the single-cubic evaluators.
!    DEVPCK tests the full PCH evaluators.
!    DEVERK exercises the error returns in all evaluators.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!     SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!     IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  integer ( kind = 4 ) I1, I2, I3, I4, I5, I6, I7, I8, I9, IFAIL, NPTS
  real ( kind = 8 ) WORK (4000)
  LOGICAL  FAIL

  IF (KPRINT .GE. 2) then
    WRITE (LUN, 1000) KPRINT
  end if
!
!  TEST DCHFDV AND DCHFEV.
!
  IFAIL = 0
  NPTS = 1000
  I1 = 1  + NPTS
  I2 = I1 + NPTS
  I3 = I2 + NPTS

  CALL DEVCHK ( LUN, KPRINT, NPTS, WORK(1), WORK(I1), WORK(I2), &
    WORK(I3), FAIL )

  IF (FAIL) then
    IFAIL = IFAIL + 1
  end if
!
!  TEST DPCHFD AND DPCHFE.
!
  I1 = 1  +  10
  I2 = I1 +  10
  I3 = I2 + 100
  I4 = I3 + 100
  I5 = I4 + 100
  I6 = I5 +  51
  I7 = I6 +  51
  I8 = I7 +  51
  I9 = I8 +  51

  CALL DEVPCK (LUN, KPRINT, WORK(1), WORK(I1), WORK(I2), WORK(I3), &
    WORK(I4), WORK(I5), WORK(I6), WORK(I7), WORK(I8), WORK(I9), FAIL)

  IF (FAIL) then
    IFAIL = IFAIL + 2
  end if
!
!  TEST ERROR RETURNS.
!
  CALL DEVERK (LUN, KPRINT, FAIL)

  IF (FAIL) then
    IFAIL = IFAIL + 4
  end if
!
!  PRINT SUMMARY AND TERMINATE.
!  At this point, IFAIL has the following value:
!  IFAIL = 0  IF ALL TESTS PASSED.
!  IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
!  IFAIL=1  IF SINGLE CUBIC  TEST FAILED. (SEE DEVCHK OUTPUT.)
!  IFAIL=2  IF DPCHFD/DPCHFE TEST FAILED. (SEE DEVPCK OUTPUT.)
!  IFAIL=4  IF ERROR RETURN  TEST FAILED. (SEE DEVERK OUTPUT.)
!
  IF ((KPRINT.GE.2).AND.(IFAIL.NE.0)) then
    WRITE (LUN, 3001)  IFAIL
  end if

  IF (IFAIL.EQ.0)  THEN
    IPASS = 1
    IF (KPRINT.GE.2) WRITE(LUN,99998)
  ELSE
    IPASS = 0
    IF (KPRINT.GE.1) WRITE(LUN,99999)
  ENDIF

  RETURN
 1000 FORMAT ('1'/' ------------ DPCHIP QUICK CHECK OUTPUT', &
         ' ------------' //20X,'( KPRINT =',I2,' )')
 3001 FORMAT (/' *** TROUBLE ***',I5,' EVALUATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL EVALUATION TESTS', &
         ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME EVALUATION TESTS', &
         ' ************')
END
SUBROUTINE DPCHQ2 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! DPCHQ2 tests the PCHIP integrators DPCHIA and DPCHID.
!
!  Discussion:
!
!    This routine constructs data from a cubic, integrates it with DPCHIA
!    and compares the results with the correct answer.
!    Since DPCHIA calls DPCHID, this tests both integrators.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    LUN   :IN  is the unit number to which output is to be written.
!
!    KPRINT:IN  controls the amount of output, as specified in the
!        SLATEC Guidelines.
!
!    IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!        IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS

  integer ( kind = 4 ) I, IEREXP(17), IERR, IFAIL, N, NPAIRS
  real ( kind = 8 ) A(17), B(17), CALC, D(7), ERRMAX, ERROR, F(7), MACHEP
  real ( kind = 8 ) ONE, THRQTR, TOL, TRUE, TWO, X(7)
  LOGICAL  FAIL, SKIP

  real ( kind = 8 ) DPCHIA, D1MACH

  real ( kind = 8 ) AX, FCN, DERIV, ANTDER
  FCN(AX) = 3.0D+00 *AX*AX*(AX-TWO)
  DERIV(AX) = 3.0D+00 *AX*(TWO*(AX-TWO) + AX)
  ANTDER(AX) = AX**3 * (THRQTR*AX - TWO)

  DATA  THRQTR /0.75D0/,  ONE /1.D0/,  TWO /2.D0/
  DATA  N /7/
  DATA  X /-4.D0, -2.D0, -0.9D0, 0.D0, 0.9D0, 2.D0, 4.D0/
  DATA  NPAIRS /17/
  DATA  A /-3.0D0, 3.0D0,-0.5D0,-0.5D0,-0.5D0,-4.0D0,-4.0D0, 3.0D0, &
       -5.0D0,-5.0D0,-6.0D0, 6.0D0,-1.5D0,-1.5D0,-3.0D0, 3.0D0, 0.5D0/
  DATA  B / 3.0D0,-3.0D0, 1.0D0, 2.0D0, 5.0D0,-0.5D0, 4.0D0, 5.0D0, &
       -3.0D0, 5.0D0,-5.0D0, 5.0D0,-0.5D0,-1.0D0,-2.5D0, 3.5D0, 0.5D0/
  DATA  IEREXP /0,0,0,0,2,0,0,2,1,3,3,3,0,0,0,0,0/
!
!  SET PASS/FAIL TOLERANCE.
!
  MACHEP = D1MACH(4)
  TOL = 100.D0*MACHEP
!
!  SET UP PCH FUNCTION DEFINITION.
!
  DO I = 1, N
    F(I) =   FCN(X(I))
    D(I) = DERIV(X(I))
  end do

  IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
  IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
  IF (KPRINT .GE. 3)  WRITE (LUN, 1002)  (X(I), F(I), D(I), I=1,N)
!
!  LOOP OVER (A,B)-PAIRS.
!
  IF (KPRINT .GE. 3)  WRITE (LUN, 2000)

  IFAIL = 0

  SKIP = .FALSE.

  DO I = 1, NPAIRS

    CALC = DPCHIA (N, X, F, D, 1, SKIP, A(I), B(I), IERR)

    IF (IERR .GE. 0)  THEN

      FAIL = IERR .NE. IEREXP(I)
      TRUE = ANTDER(B(I)) - ANTDER(A(I))
      ERROR = CALC - TRUE

      IF (KPRINT .GE. 3)  THEN
        IF (FAIL)  THEN
          WRITE (LUN, 2001) A(I), B(I), IERR, TRUE, CALC, ERROR, IEREXP(I)
        ELSE
         WRITE (LUN, 2002) A(I), B(I), IERR, TRUE, CALC, ERROR
        ENDIF
      ENDIF

      ERROR = ABS(ERROR) / MAX(ONE, ABS(TRUE))

      IF (FAIL .OR. (ERROR.GT.TOL))  IFAIL = IFAIL + 1

      IF (I .EQ. 1)  THEN
        ERRMAX = ERROR
      ELSE
        ERRMAX = MAX(ERRMAX, ERROR)
      ENDIF

    ELSE

      IF (KPRINT .GE. 3)  WRITE (LUN, 2002)  A(I), B(I), IERR
      IFAIL = IFAIL + 1

    ENDIF

  end do
!
!  PRINT SUMMARY.
!
  IF (KPRINT .GE. 2)  THEN
    WRITE (LUN, 2003)  ERRMAX, TOL
    IF (IFAIL .NE. 0)  WRITE (LUN, 3001)  IFAIL
  ENDIF

  IF (IFAIL.EQ.0)  THEN
    IPASS = 1
    IF (KPRINT.GE.2) WRITE(LUN,99998)
  ELSE
    IPASS = 0
    IF (KPRINT.GE.1) WRITE(LUN,99999)
  ENDIF

  RETURN
 1000 FORMAT ('1'//10X,'TEST DPCHIP INTEGRATORS')
 1001 FORMAT (//10X,'DPCHQ2 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' //11X,'X',9X,'F',9X,'D' /(5X,3F10.3) )
 2000 FORMAT (// 5X,'TEST RESULTS:' &
         //'    A     B    ERR     TRUE',16X,'CALC',15X,'ERROR')
 2001 FORMAT (2F6.1,I5,1P,2D20.10,D15.5,'  (',I1,') *****' )
 2002 FORMAT (2F6.1,I5,1P,2D20.10,D15.5)
 2003 FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,D15.5, &
                ',   TOLERANCE:',1P,D15.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS', &
         ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME INTEGRATION TESTS', &
         ' ************')
END
SUBROUTINE DPCHQ3 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! DPCHQ3 tests the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
!
!  Discussion:
!
!    This routine interpolates a constructed data set with all three
!    DPCHIP interpolators and compares the results with those obtained
!    on a Cray X/MP.  Two different values of the DPCHIC parameter SWITCH
!    are used.
!
!    1. The Cray results are given only to nine significant figures,
!       so don't expect them to match to more.
!    2. The results will depend to some extent on the accuracy of
!       the EXP function.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!        SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!        IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  LOGICAL  COMP
  real ( kind = 8 ) D1MACH
  integer ( kind = 4 ) I, IC(2), IERR, IFAIL, N, NBAD, NBADZ, NWK
  PARAMETER  (N = 9,  NWK = 2*N)
  real ( kind = 8 ) D(N), DC(N), DC5, DC6, DM(N), DS(N), ERR, F(N)
  real ( kind = 8 ) MONE, TOL, TOLD, TOLZ, VC(2), X(N), WK(NWK), ZERO
  PARAMETER  (ZERO = 0.0D0,  MONE = -1.0D0)
  CHARACTER*6  RESULT

  DATA  IC /0, 0/
  DATA  X /-2.2D0,-1.2D0,-1.0D0,-0.5D0,-0.01D0, 0.5D0, 1.0D0, &
       2.0D0, 2.2D0/
!
!  Results generated on Cray X/MP (9 sign. figs.)
!
  DATA  DM / 0.    , 3.80027352D-01, 7.17253009D-01, &
         5.82014161D-01, 0.    ,-5.68208031D-01, &
        -5.13501618D-01,-7.77910977D-02,-2.45611117D-03/
  DATA  DC5,DC6 / 1.76950158D-02,-5.69579814D-01/
  DATA  DS /-5.16830792D-02, 5.71455855D-01, 7.40530225D-01, &
         7.63864934D-01, 1.92614386D-02,-7.65324380D-01, &
        -7.28209035D-01,-7.98445427D-02,-2.85983446D-02/
  IFAIL = 0
!
!  Set tolerances.
!
  TOL  = 10*D1MACH(4)
  TOLD = MAX( 1.0D-7, 10*TOL )
  TOLZ = ZERO

  IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
  IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
!
!  Set up data.
!
  DO I = 1, N
     F(I) = EXP(-X(I)**2)
  end do

  IF (KPRINT .GE. 3)  THEN
    WRITE (LUN, 1002)
    DO I = 1, 4
      WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
    end do
    WRITE (LUN, 1011)  X(5), F(5), DM(5), DC5, DS(5)
    WRITE (LUN, 1011)  X(6), F(6), DM(6), DC6, DS(6)
    DO I = 7, N
      WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
    end do
  ENDIF
!
!  Test DPCHIM.
!
  IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IM'

  CALL DPCHIM (N, X, F, D, 1, IERR)
!
!  Expect IERR=1 (one monotonicity switch).
!
  IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 1
  IF ( .NOT.COMP (IERR, 1, LUN, KPRINT) )  THEN
     IFAIL = IFAIL + 1
  ELSE
     IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
     NBAD = 0
     NBADZ = 0

     DO I = 1, N
    RESULT = '  OK'
!
!  D-values should agree with stored values.
!  (Zero values should agree exactly.)
!
    IF ( DM(I).EQ.ZERO )  THEN
       ERR = ABS( D(I) )
       IF ( ERR.GT.TOLZ )  THEN
      NBADZ = NBADZ + 1
      RESULT = '**BADZ'
       ENDIF
    ELSE
       ERR = ABS( (D(I)-DM(I))/DM(I) )
       IF ( ERR.GT.TOLD )  THEN
      NBAD = NBAD + 1
      RESULT = '**BAD'
       ENDIF
    ENDIF
    IF (KPRINT.GE.3) then
      WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
    end if
     end do

     IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
    IFAIL = IFAIL + 1
    IF ((NBADZ.NE.0).AND.(KPRINT.GE.2)) then
       WRITE (LUN, 2004)  NBAD
    end if
    IF ((NBAD.NE.0).AND.(KPRINT.GE.2)) then
       WRITE (LUN, 2005)  NBAD, 'IM', TOLD
    end if
     ELSE
    IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IM'
     ENDIF

  ENDIF
!
!  Test DPCHIC -- options set to reproduce DPCHIM.
!
  IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'

  CALL DPCHIC (IC, VC, ZERO, N, X, F, DC, 1, WK, NWK, IERR)
!
!  Expect IERR=0 .
!
  IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
  IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
     IFAIL = IFAIL + 1
  ELSE
     IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
     NBAD = 0
     DO I = 1, N
    RESULT = '  OK'
!
!  D-values should agree exactly with those computed by DPCHIM.
!  (To be generous, will only test to machine precision.)
!
    ERR = ABS( D(I)-DC(I) )
    IF ( ERR.GT.TOL )  THEN
       NBAD = NBAD + 1
       RESULT = '**BAD'
    ENDIF
    IF (KPRINT.GE.3) then
       WRITE (LUN, 2003)  I, X(I), DC(I), ERR, RESULT
    end if
     end do

     IF ( NBAD.NE.0 )  THEN
    IFAIL = IFAIL + 1
    IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'IC', TOL
     ELSE
    IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
     ENDIF
  ENDIF
!
!  Test DPCHIC -- default nonzero switch derivatives.
!
  IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'

  CALL DPCHIC (IC, VC, MONE, N, X, F, D, 1, WK, NWK, IERR)
!
!    Expect IERR=0 .
!
  IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
  IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
     IFAIL = IFAIL + 1
  ELSE
     IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
     NBAD = 0
     NBADZ = 0

     DO I = 1, N
    RESULT = '  OK'
!
!  D-values should agree exactly with those computed in
!  previous call, except at points 5 and 6.
!
    IF ( (I.LT.5).OR.(I.GT.6) )  THEN
       ERR = ABS( D(I)-DC(I) )
       IF ( ERR.GT.TOLZ )  THEN
      NBADZ = NBADZ + 1
      RESULT = '**BADA'
       ENDIF
    ELSE
       IF ( I.EQ.5 )  THEN
      ERR = ABS( (D(I)-DC5)/DC5 )
       ELSE
      ERR = ABS( (D(I)-DC6)/DC6 )
       ENDIF
       IF ( ERR.GT.TOLD )  THEN
      NBAD = NBAD + 1
      RESULT = '**BAD'
       ENDIF
    ENDIF
    IF (KPRINT.GE.3) then
       WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
    end if
     end do

     IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
    IFAIL = IFAIL + 1
    IF ((NBADZ.NE.0).AND.(KPRINT.GE.2)) then
       WRITE (LUN, 2007)  NBAD
    end if
    IF ((NBAD.NE.0).AND.(KPRINT.GE.2)) then
       WRITE (LUN, 2005)  NBAD, 'IC', TOLD
    end if
     ELSE
    IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
     ENDIF
  ENDIF
!
!  Test DPCHSP.
!
  IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'SP'

  CALL DPCHSP (IC, VC, N, X, F, D, 1, WK, NWK, IERR)
!
!  Expect IERR=0 .
!
  IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
  IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
     IFAIL = IFAIL + 1
  ELSE
     IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
     NBAD = 0

     DO I = 1, N
    RESULT = '  OK'
!
!  D-values should agree with stored values.
!
    ERR = ABS( (D(I)-DS(I))/DS(I) )
    IF ( ERR.GT.TOLD )  THEN
       NBAD = NBAD + 1
       RESULT = '**BAD'
    ENDIF
    IF (KPRINT.GE.3) then
      WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
    end if
     end do

     IF ( NBAD.NE.0 )  THEN
    IFAIL = IFAIL + 1
    IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'SP', TOLD
     ELSE
    IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'SP'
     ENDIF
  ENDIF
!
!  PRINT SUMMARY AND TERMINATE.
!
  IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

  IF (IFAIL.EQ.0)  THEN
     IPASS = 1
     IF (KPRINT.GE.2) WRITE(LUN,99998)
  ELSE
     IPASS = 0
     IF (KPRINT.GE.1) WRITE(LUN,99999)
  ENDIF

  RETURN

 1000 FORMAT ('1'//10X,'TEST DPCHIP INTERPOLATORS')
 1001 FORMAT (//10X,'DPCHQ3 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' &
      /39X,'---------- EXPECTED D-VALUES ----------' &
      /12X,'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
 1010 FORMAT (5X,F10.2,1P,D15.5,4X,D15.5,15X,D15.5)
 1011 FORMAT (5X,F10.2,1P,D15.5,4X,3D15.5)
 2000 FORMAT (/5X,'DPCH',A2,' TEST:')
 2001 FORMAT (15X,'EXPECT  IERR =',I5)
 2002 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
 2003 FORMAT (5X,I5,F10.2,1P,2D15.5,2X,A)
 2004 FORMAT (/'    **',I5,' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
 2005 FORMAT (/'    **',I5,' DPCH',A2,' RESULTS FAILED TOLERANCE TEST.', &
          '  TOL =',1P,D10.3)
 2006 FORMAT (/5X,'  ALL DPCH',A2,' RESULTS OK.')
 2007 FORMAT (/'    **',I5,' DPCHIC RESULTS FAILED TO AGREE WITH', &
           ' PREVIOUS CALL.')
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS', &
     ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS', &
     ' ************')
END
SUBROUTINE DPCHQ4 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! DPCHQ4 tests the PCHIP monotonicity checker DPCHCM.
!
!  Discussion:
!
!    This routine tests a constructed data set with three different
!    INCFD settings and compares with the expected results.  It then
!    runs a special test to check for bug in overall monotonicity found
!    in DPCHMC.  Finally, it reverses the data and repeats all tests.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  integer ( kind = 4 ) MAXN, MAXN2, MAXN3, NB
  PARAMETER  (MAXN = 16,  MAXN2 = 8,  MAXN3 = 6,  NB = 7)
  integer ( kind = 4 ) I, IERR, IFAIL, INCFD, ISMEX1(MAXN), ISMEX2(MAXN2)
  integer ( kind = 4 ) ISMEX3(MAXN3), ISMEXB(NB), ISMON(MAXN), K, N, NS(3)
  real ( kind = 8 ) D(MAXN), DB(NB), F(MAXN), FB(NB), X(MAXN)
  LOGICAL  SKIP
!
!  DEFINE EXPECTED RESULTS.
!
  DATA  ISMEX1 / 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 2/
  DATA  ISMEX2 / 1, 2, 2, 1, 2, 2, 1, 2/
  DATA  ISMEX3 / 1, 1, 1, 1, 1, 1/
  DATA  ISMEXB / 1, 3, 1, -1, -3, -1, 2/
!
!  DEFINE TEST DATA.
!
  DATA  NS /16, 8, 6/

      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
!
!  Define X, F, D.
!
      DO I = 1, MAXN
         X(I) = I
         D(I) = 0.D0
      end do

      DO I = 2, MAXN, 3
         D(I) = 2.D0
      end do

      DO I = 1, 3
         F(I) = X(I)
         F(I+ 3) = F(I  ) + 1.D0
         F(I+ 6) = F(I+3) + 1.D0
         F(I+ 9) = F(I+6) + 1.D0
         F(I+12) = F(I+9) + 1.D0
      end do

      F(16) = 6.D0
!
!  Define FB, DB.
!
      FB(1) = 0.D0
      FB(2) = 2.D0
      FB(3) = 3.D0
      FB(4) = 5.D0
      DB(1) = 1.D0
      DB(2) = 3.D0
      DB(3) = 3.D0
      DB(4) = 0.D0
      DO I = 1, 3
         FB(NB-I+1) =  FB(I)
         DB(NB-I+1) = -DB(I)
      end do
!
!  INITIALIZE.
!
      IFAIL = 0

      IF (KPRINT .GE. 3)  THEN
         WRITE (LUN, 1002)
         DO I = 1, NB
            WRITE (LUN, 1010)  I, X(I), F(I), D(I), FB(I), DB(I)
         end do

         DO I = NB+1, MAXN
            WRITE (LUN, 1010)  I, X(I), F(I), D(I)
        end do
      ENDIF
!
!  TRANSFER POINT FOR SECOND SET OF TESTS.
!
   25 CONTINUE
!
!  Loop over a series of values of INCFD.
!
      DO INCFD = 1, 3

         N = NS(INCFD)
         SKIP = .FALSE.

         CALL DPCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)

         IF (KPRINT.GE.3) then
            WRITE (LUN, 2000)  INCFD, IERR, (ISMON(I), I=1,N)
         end if

         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN,2001)
         ELSE
            DO I = 1, N
               IF (INCFD.EQ.1)  THEN
                  IF ( ISMON(I).NE.ISMEX1(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                        WRITE (LUN, 2002)  (ISMEX1(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ELSE IF (INCFD.EQ.2) THEN
                  IF ( ISMON(I).NE.ISMEX2(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                        WRITE (LUN, 2002)  (ISMEX2(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ELSE
                  IF ( ISMON(I).NE.ISMEX3(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                        WRITE (LUN, 2002)  (ISMEX3(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ENDIF
            end do
         ENDIF

30      continue

      end do
!
!  Test for -1,3,1 bug.
!
      SKIP = .FALSE.

      CALL DPCHCM (NB, X, FB, DB, 1, SKIP, ISMON, IERR)

      IF (KPRINT.GE.3) then
         WRITE (LUN, 2030)  IERR, (ISMON(I), I=1,NB)
      end if

      IF ( IERR.NE.0 )  THEN
         IFAIL = IFAIL + 1
         IF (KPRINT.GE.3)  WRITE (LUN,2001)
      ELSE
         DO I = 1, NB
            IF ( ISMON(I).NE.ISMEXB(I) )  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.3) then
                  WRITE (LUN, 2002)  (ISMEXB(K),K=1,NB)
               end if
               GO TO 35
            ENDIF
         end do
      ENDIF
   35 CONTINUE

      IF (F(1).LT.0.)  GO TO 90
!
!  Change sign and do again.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2050)

      DO I = 1, MAXN
         F(I) = -F(I)
         D(I) = -D(I)
         IF ( ISMEX1(I).NE.2 )  ISMEX1(I) = -ISMEX1(I)
      end do

      DO I = 1, MAXN2
         IF ( ISMEX2(I).NE.2 )  ISMEX2(I) = -ISMEX2(I)
      end do

      DO I = 1, MAXN3
         IF ( ISMEX3(I).NE.2 )  ISMEX3(I) = -ISMEX3(I)
      end do

      DO I = 1, NB
         FB(I) = -FB(I)
         DB(I) = -DB(I)
         IF ( ISMEXB(I).NE.2 )  ISMEXB(I) = -ISMEXB(I)
      end do
      GO TO 25
!
!  PRINT SUMMARY AND TERMINATE.
!
   90 CONTINUE

      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST DPCHIP MONOTONICITY CHECKER')
 1001 FORMAT (//10X,'DPCHQ4 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' &
             // 9X,'I',4X,'X',5X,'F',5X,'D',5X,'FB',4X,'DB')
 1010 FORMAT (5X,I5,5F6.1)
 2000 FORMAT (/4X,'INCFD =',I2,':  IERR =',I3/15X,'ISMON =',16I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (' *** Failed -- expect:',16I3)
 2030 FORMAT (/4X,' Bug test:  IERR =',I3/15X,'ISMON =',7I3)
 2050 FORMAT (/4X,'Changing sign of data.....')
 3001 FORMAT (/' *** TROUBLE ***',I5,' MONOTONICITY TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL MONOTONICITY TESTS', &
             ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME MONOTONICITY TESTS', &
             ' ************')
END
SUBROUTINE DPCHQ5 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! DPCHQ5 tests the PCH to B-spline conversion routine DPCHBS.
!
!  Discussion:
!
!    This routine tests a constructed data set with four different
!    KNOTYP settings.  It computes the function and derivatives of the
!    resulting B-representation via DBVALU and compares with PCH data.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  real ( kind = 8 ) DBVALU, D1MACH
  EXTERNAL  DBVALU, DPCHBS, D1MACH
  integer ( kind = 4 ) I, IERR, IFAIL, INBV, J, KNOTYP, K, N, NDIM, NKNOTS
  PARAMETER  (N = 9)
  real ( kind = 8 ) BCOEF(2*N), D(N), DCALC, DERR, DERMAX, F(N)
  real ( kind = 8 ) FCALC, FERR, FERMAX, T(2*N+4), TERR, TERMAX, TOL, TOLZ
  real ( kind = 8 ) TSAVE(2*N+4), WORK(16*N), X(N), ZERO
  PARAMETER  (ZERO = 0.0D0)
  LOGICAL  FAIL
  real ( kind = 8 ) ANS, ERR, RELERR
  RELERR (ERR, ANS) = ABS(ERR) / MAX(1.0D-5,ABS(ANS))
!
!  Define test data.
!
      DATA  X /-2.2D0,   -1.2D0,   -1.0D0,   -0.5D0,   -0.01D0, &
                0.5D0,    1.0D0,    2.0D0,    2.2D0/
      DATA  F / 0.0079D0, 0.2369D0, 0.3679D0, 0.7788D0, 0.9999D0, &
                0.7788D0, 0.3679D0, 0.1083D0, 0.0079D0/
      DATA  D / 0.0000D0, 0.3800D0, 0.7173D0, 0.5820D0, 0.0177D0, &
               -0.5696D0,-0.5135D0,-0.0778D0,-0.0025D0/

      IFAIL = 0
      TOL = 100*D1MACH(4)
      TOLZ = ZERO

      IF (KPRINT.GE.3)  WRITE (LUN, 1000)
      IF (KPRINT.GE.2)  WRITE (LUN, 1001)
!
!  Loop over a series of values of KNOTYP.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 1010)

      DO KNOTYP = 2, -1, -1

         CALL DPCHBS (N, X, F, D, 1, KNOTYP, NKNOTS, T, BCOEF, NDIM, K, IERR)

         IF (KPRINT.GE.3) then
           WRITE (LUN, 2000) KNOTYP, NKNOTS, NDIM, K, IERR
         end if

         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN, 2001)
         ELSE
!
!  Compare evaluated results with inputs to DPCHBS.
!
            INBV = 1
            FERMAX = ZERO
            DERMAX = ZERO

            IF (KPRINT.GE.3)  THEN
               WRITE (LUN, 2002)
               WRITE (LUN, 2003)  T(1), T(2)
               J = 1
            ENDIF

            DO I = 1, N
               FCALC = DBVALU (T, BCOEF, NDIM, K, 0, X(I), INBV, WORK)
               FERR = F(I) - FCALC
               FERMAX = MAX(FERMAX, RELERR(FERR,F(I)) )
               DCALC = DBVALU (T, BCOEF, NDIM, K, 1, X(I), INBV, WORK)
               DERR = D(I) - DCALC
               DERMAX = MAX(DERMAX, RELERR(DERR,D(I)) )
               IF (KPRINT.GE.3)  THEN
                  J = J + 2
                  WRITE (LUN, 2004)  X(I), T(J), T(J+1), F(I), FERR, D(I), DERR
               ENDIF
            end do

            IF (KPRINT.GE.3)  THEN
               J = J + 2
               WRITE (LUN, 2003)  T(J), T(J+1)
            ENDIF

            FAIL = (FERMAX.GT.TOL).OR.(DERMAX.GT.TOL)
            IF (FAIL)  IFAIL = IFAIL + 1

            IF ((KPRINT.GE.3).OR.(KPRINT.GE.2).AND.FAIL) then
              WRITE (LUN, 2005)  FERMAX, DERMAX, TOL
            end if

         ENDIF
!
!  Special check for KNOTYP=-1.
!
         IF (KNOTYP.EQ.0)  THEN
!
!  Save knot vector for next test.
!
            DO I = 1, NKNOTS
               TSAVE(I) = T(I)
            end do

         ELSE IF (KNOTYP.EQ.-1)  THEN
!
!  Check that knot vector is unchanged.
!
            TERMAX = ZERO
            DO I = 1, NKNOTS
               TERR = ABS(T(I) - TSAVE(I))
               TERMAX = MAX(TERMAX, TERR)
            end do

            IF (TERMAX.GT.TOLZ)  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.2)  WRITE (LUN, 2007)  TERMAX, TOLZ
            ENDIF

         ENDIF

      end do
!
!  PRINT SUMMARY AND TERMINATE.
!
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCH TO B-SPLINE CONVERTER')
 1001 FORMAT (//10X,'DPCHQ5 RESULTS'/10X,'--------------')
 1010 FORMAT (/4X,'(Results should be the same for all KNOTYP values.)')
 2000 FORMAT (/4X,'KNOTYP =',I2,':  NKNOTS =',I3,',  NDIM =',I3, &
                              ',  K =',I2,',  IERR =',I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (/15X,'X',9X,'KNOTS',10X,'F',7X,'FERR',8X,'D',7X,'DERR')
 2003 FORMAT (18X,2F8.2)
 2004 FORMAT (10X,3F8.2,F10.4,1P,D10.2,0P,F10.4,1P,D10.2)
 2005 FORMAT (/5X,'Maximum relative errors:' &
            /15X,'F-error =',1P,D13.5,5X,'D-error =',D13.5 &
             /5X,'Both should be less than  TOL =',D13.5)
 2007 FORMAT (/' *** T-ARRAY MAXIMUM CHANGE =',1P,D13.5, &
                ';  SHOULD NOT EXCEED TOLZ =',D13.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' CONVERSION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL CONVERSION TESTS', &
             ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME CONVERSION TESTS', &
             ' ************')
END
FUNCTION DPCHQA ( N, X, F, D, A, B, IERR )

!*****************************************************************************80
!
!! DPCHQA: definite integral of spline or piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Evaluates the definite integral of the cubic Hermite or spline
!    function defined by  N, X, F, D  over the interval [A, B].  This
!    is an easy to use driver for the routine DPCHIA by F.N. Fritsch.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence:
!
!           VALUE = DPCHQA (N, X, F, D, A, B, IERR)
!
!     integer ( kind = 4 )  N, IERR
!     real ( kind = 8 )  X(N), F(N), D(N), A, B
!
!  Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real ( kind = 8 ) array of independent variable
!           values.  The elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real ( kind = 8 ) array of function values.
!           F(I) is the value corresponding to X(I).
!
!     D -- (input) real ( kind = 8 ) array of derivative values.  D(I) is
!           the value corresponding to X(I).
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (Value has not been computed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  real ( kind = 8 ) dpchqa
  integer ( kind = 4 )  N, IERR
  real ( kind = 8 )  X(N), F(N), D(N), A, B
  integer ( kind = 4 )  INCFD
  real ( kind = 8 )  DPCHIA
  LOGICAL SKIP

  incfd = 1
  skip = .true.

  DPCHQA  =  DPCHIA ( N, X, F, D, INCFD, SKIP, A, B, IERR )

  RETURN
END
SUBROUTINE DPCHSP ( IC, VC, N, X, F, D, INCFD, WK, NWK, IERR )

!*****************************************************************************80
!
!! DPCHSP: derivatives for Hermite representation of cubic spline interpolant.
!
!  Discussion:
!
!    Set derivatives needed to determine the Hermite represen-
!    tation of the cubic spline interpolant to given data, with
!    specified boundary conditions.
!
!    Computes the Hermite representation of the cubic spline inter-
!    polant to the data given in X and F satisfying the boundary
!    conditions specified by IC and VC.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by DPCHFE or DPCHFD.
!
!    This is a modified version of C. de Boor's cubic spline
!    routine CUBSPL.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, NWK, IERR
!        real ( kind = 8 )  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor's cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine the cubic spline
!           interpolant with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real*8 array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK.LT.2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
!              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
  implicit none

  integer ( kind = 4 )  IC(2), N, INCFD, NWK, IERR
  real ( kind = 8 )  VC(2), X(*), F(INCFD,*), D(INCFD,*), WK(2,*)
  integer ( kind = 4 )  IBEG, IEND, INDEX, J, NM1
  real ( kind = 8 )  G, HALF, ONE, STEMP(3), XTEMP(4)

  SAVE HALF, ONE
  real ( kind = 8 )  DPCHDF

  DATA  HALF/.5D0/, ONE/1.D0/
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO J = 2, N
     IF ( X(J).LE.X(J-1) )  GO TO 5003
  end do

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
  IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
  IF ( IERR.LT.0 )  GO TO 5004
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  IF ( NWK .LT. 2*N )  GO TO 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
!
  DO J=2,N
     WK(1,J) = X(J) - X(J-1)
     WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
  end do
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
  IF ( IBEG.GT.N )  IBEG = 0
  IF ( IEND.GT.N )  IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
  IF ( (IBEG == 1).OR.(IBEG == 2) )  THEN
     D(1,1) = VC(1)
  ELSE IF (IBEG .GT. 2)  THEN
!
!  PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
!
     DO  J = 1, IBEG
        INDEX = IBEG-J+1
!
!  INDEX RUNS FROM IBEG DOWN TO 1.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
     end do

     D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)

     IF (IERR .NE. 0)  GO TO 5009
     IBEG = 1
  end if

  IF ( (IEND == 1).OR.(IEND == 2) )  THEN
     D(1,N) = VC(2)
  ELSE IF (IEND .GT. 2)  THEN
!
!  PICK UP LAST IEND POINTS.
!
     DO J = 1, IEND
        INDEX = N-IEND+J
!
!  INDEX RUNS FROM N+1-IEND UP TO N.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
     end do

     D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)

     IF (IERR .NE. 0)  GO TO 5009
     IEND = 1
  end if
!
!  BEGIN CODING FROM CUBSPL
!
!  A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!  WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!  WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
  IF (IBEG == 0)  THEN
     IF (N == 2)  THEN
!
!  NO CONDITION AT LEFT END AND N = 2.
!
        WK(2,1) = ONE
        WK(1,1) = ONE
        D(1,1) = 2.0D+00 * WK(2,2)
     ELSE
!
!  NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
!
        WK(2,1) = WK(1,3)
        WK(1,1) = WK(1,2) + WK(1,3)
        D(1,1) =((WK(1,2) + 2.0D+00 * WK(1,1))*WK(2,2)*WK(1,3) &
                          + WK(1,2)**2*WK(2,3)) / WK(1,1)
     end if
  ELSE IF (IBEG == 1)  THEN
!
!  SLOPE PRESCRIBED AT LEFT END.
!
     WK(2,1) = ONE
     WK(1,1) = 0.0D+00
  ELSE
!
!  SECOND DERIVATIVE PRESCRIBED AT LEFT END.
!
     WK(2,1) = 2.0D+00
     WK(1,1) = ONE
     D(1,1) = 3.0D+00 *WK(2,2) - HALF*WK(1,2)*D(1,1)
  end if
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
  NM1 = N-1
  IF (NM1 .GT. 1)  THEN
     DO J=2,NM1
        IF (WK(2,J-1) == 0.0D+00 )  GO TO 5008
        G = -WK(1,J+1)/WK(2,J-1)
        D(1,J) = G*D(1,J-1) &
                    + 3.0D+00 *(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
        WK(2,J) = G*WK(1,J-1) + 2.0D+00 * (WK(1,J) + WK(1,J+1))
     end do
  end if
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!  (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!  IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!  SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!  AT THIS POINT.
!
  IF (IEND == 1)  GO TO 30

  IF (IEND == 0)  THEN
     IF (N == 2 .AND. IBEG == 0)  THEN
!
!  NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
!
        D(1,2) = WK(2,2)
        GO TO 30
     ELSE IF ((N == 2) .OR. (N == 3 .AND. IBEG == 0))  THEN
!
!  EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!  NOT-A-KNOT AT LEFT END POINT).
!
        D(1,N) = 2.0D+00 * WK(2,N)
        WK(2,N) = ONE
        IF (WK(2,N-1) == 0.0D+00 )  GO TO 5008
        G = -ONE/WK(2,N-1)
     ELSE
!
!  NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
!  KNOT AT LEFT END POINT.
!
        G = WK(1,N-1) + WK(1,N)
!
!  DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
!
        D(1,N) = ((WK(1,N)+2.0D+00 * G)*WK(2,N)*WK(1,N-1) &
                    + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
        IF (WK(2,N-1) == 0.0D+00 )  GO TO 5008
        G = -G/WK(2,N-1)
        WK(2,N) = WK(1,N-1)
     end if
  ELSE
!
!  SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
!
     D(1,N) = 3.0D+00 *WK(2,N) + HALF*WK(1,N)*D(1,N)
     WK(2,N) = 2.0D+00
     IF (WK(2,N-1) == 0.0D+00 )  GO TO 5008
     G = -ONE/WK(2,N-1)
  end if
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
  WK(2,N) = G*WK(1,N-1) + WK(2,N)
  IF (WK(2,N) == 0.0D+00 )   GO TO 5008
  D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
   30 CONTINUE
  DO J=NM1,1,-1
     IF (WK(2,J) == 0.0D+00 )  GO TO 5008
     D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
  end do
!
!  END  CODING FROM CUBSPL.
!
  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHSP', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHSP', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'DPCHSP', &
     'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  CALL XERMSG ('SLATEC', 'DPCHSP', 'IC OUT OF RANGE', IERR, 1)
  RETURN
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
  IERR = -7
  CALL XERMSG ('SLATEC', 'DPCHSP', 'WORK ARRAY TOO SMALL', IERR, 1)
  RETURN
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES
!   ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3).
  IERR = -8
  CALL XERMSG ('SLATEC', 'DPCHSP', 'SINGULAR LINEAR SYSTEM', IERR, 1)
  RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM DPCHDF.  THIS CASE SHOULD NEVER OCCUR.
  IERR = -9
  CALL XERMSG ('SLATEC', 'DPCHSP', 'ERROR RETURN FROM DPCHDF', &
   IERR, 1)

  RETURN
END
FUNCTION DPCHST (ARG1, ARG2)

!*****************************************************************************80
!
!! DPCHST: carry out a sign test.
!
!  Discussion:
!
!    The object is to do this without multiplying ARG1*ARG2, to avoid
!    possible over/underflow problems.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real ( kind = 8 )  ARG1, ARG2
  real ( kind = 8 ) dpchst

  DPCHST = SIGN ( 1.0D+00, ARG1 ) * SIGN ( 1.0D+00, ARG2 )

  IF ((ARG1 == 0.0D+00) .OR. (ARG2 == 0.0D+00)) then
    DPCHST = 0.0D+00
  end if

  RETURN
END
SUBROUTINE DPCHSW ( DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR )

!*****************************************************************************80
!
!! DPCHSW limits excursion from data for DPCHCS.
!
!  Discussion:
!
!    Called by DPCHCS to adjust D1 and D2 if necessary to insure that
!    the extremum on this interval is not further than DFMAX from the
!    extreme data value.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  IEXTRM, IERR
!        real ( kind = 8 )  DFMAX, D1, D2, H, SLOPE
!
!        CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!  Parameters:
!
!     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!           the cubic determined by derivative values D1,D2.  (assumes
!           DFMAX.GT.0.)
!
!     IEXTRM -- (input) index of the extreme data value.  (assumes
!           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!           (Assumes D1*D2 .LE. 0.)
!          (output) may be modified if necessary to meet the restriction
!           imposed by DFMAX.
!
!     H -- (input) interval length.  (Assumes  H.GT.0.)
!
!     SLOPE -- (input) data slope on the interval.
!
!     IERR -- (output) error flag.  should be zero.
!           If IERR=-1, assumption on D1 and D2 is not satisfied.
!           If IERR=-2, quadratic equation locating extremum has
!                       negative discriminant (should never occur).
!
  implicit none

  integer ( kind = 4 )  IEXTRM, IERR
  real ( kind = 8 )  DFMAX, D1, D2, H, SLOPE
  real ( kind = 8 )  CP, FACT, HPHI, LAMBDA, NU, ONE, PHI, RADCAL
  real ( kind = 8 )  RHO, SIGMA, SMALL, THAT, THIRD

  SAVE ONE, FACT
  SAVE THIRD
  real ( kind = 8 )  D1MACH

  DATA ONE /1.D0/
  data FACT /100.D0/
!        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
  DATA  THIRD /0.33333D0/
!
!  NOTATION AND GENERAL REMARKS.
!
!  RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!  LAMBDA IS THE RATIO OF D2 TO D1.
!  THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!  PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!  WHERE  THAT = (XHAT - X1)/H .
!  THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!  SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!
!  SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!
  SMALL = FACT * D1MACH(4)
!
!  MAIN CALCULATION.
!
  IF (D1 == 0.0D+00 )  THEN
!
!  SPECIAL CASE -- D1 == ZERO .
!
!  IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
!
     IF (D2 == 0.0D+00 )  GO TO 5001

     RHO = SLOPE/D2
!
!  EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
!
     IF (RHO .GE. THIRD) then
       return
     end if

     THAT = (2.0D+00 * ( 3.0D+00 *RHO-ONE)) / ( 3.0D+00 *( 2.0D+00 * RHO-ONE))
     PHI = THAT**2 * (( 3.0D+00 *RHO-ONE)/ 3.0D+00 )
!
!  CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
!
     IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!  TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
!
     HPHI = H * ABS(PHI)
     IF (HPHI*ABS(D2) .GT. DFMAX)  THEN
!
!  AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
!
        D2 = SIGN (DFMAX/HPHI, D2)
     end if
  ELSE

     RHO = SLOPE/D1
     LAMBDA = -D2/D1
     IF (D2 == 0.0D+00 )  THEN
!
!  SPECIAL CASE:  D2 == ZERO .
!  EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
!
        IF (RHO .GE. THIRD) then
          ierr = 0
          return
        end if

        CP = 2.0D+00 - 3.0D+00 *RHO
        NU = ONE - 2.0D+00 * RHO
        THAT = ONE / ( 3.0D+00 *NU)
     ELSE
        IF (LAMBDA .LE. 0.0D+00 )  GO TO 5001
!
!  NORMAL CASE:  D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!
        NU = ONE - LAMBDA - 2.0D+00 * RHO
        SIGMA = ONE - RHO
        CP = NU + SIGMA
        IF (ABS(NU) .GT. SMALL)  THEN
           RADCAL = (NU - ( 2.0D+00 * RHO+ONE))*NU + SIGMA**2
           IF (RADCAL .LT. 0.0D+00 )  GO TO 5002
           THAT = (CP - SQRT(RADCAL)) / ( 3.0D+00 *NU)
        ELSE
           THAT = ONE/( 2.0D+00 * SIGMA)
        end if
     end if
     PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!
!  CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
!
     IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!  TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
!
     HPHI = H * ABS(PHI)
     IF (HPHI*ABS(D1) .GT. DFMAX)  THEN
!
!  AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
!
        D1 = SIGN (DFMAX/HPHI, D1)
        D2 = -LAMBDA*D1
     end if
  end if

  IERR = 0

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'DPCHSW', 'D1 AND/OR D2 INVALID', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'DPCHSW', 'NEGATIVE RADICAL', IERR, 1)

  RETURN
END
SUBROUTINE EVCHCK (LOUT, KPRINT, NPTS, XEV, FEV, DEV, FEV2, FAIL)

!*****************************************************************************80
!
!! EVCHK tests evaluation accuracy of CHFDV and DHFEV for PCHQK1.
!
!  Discussion:
!
!    USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
!    DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
!    1. CHECKS THAT CHFDV AND CHFEV BOTH REPRODUCE ENDPOINT VALUES.
!    2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
!    AND:
!    A. CHECKS ACCURACY OF CHFDV FUNCTION AND DERIVATIVE VALUES
!       AGAINST EXACT VALUES.
!    B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
!    C. CHECKS THAT FUNCTION VALUES FROM CHFEV AGREE WITH THOSE
!       FROM CHFDV.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    LOUT
!
!    KPRINT
!
!    NPTS
!
!    XEV
!
!    FEV
!
!    DEV
!
!    FEV2
!
!    FAIL
!
  implicit none

  integer ( kind = 4 ) LOUT, KPRINT, NPTS
  real ( kind = 4 ) XEV(*), FEV(*), DEV(*), FEV2(*)
  LOGICAL  FAIL
  integer ( kind = 4 ) I, IERR, IINT, NEXT(2), NEXT2(2), NINT
  real ( kind = 4 ) AED, AED2, AEDMAX, AEDMIN, AEF, AEF2, AEFMAX, AEFMIN
  real ( kind = 4 ) CHECK(2), CHECKF(2), CHECKD(2), D1, D2, DERMAX, DTRUE, DX
  real ( kind = 4 ) EPS1, EPS2, F1, F2, FACT, FERMAX, FLOORD, FLOORF, FOUR
  real ( kind = 4 ) FTRUE, LEFT(3), MACHEP
  real ( kind = 4 ) ONE, RED, RED2, REDMAX, REDMIN, REF, REF2, REFMAX
  real ( kind = 4 ) REFMIN, RIGHT(3), SMALL, TEN, TOL1, TOL2
  real ( kind = 4 ) X1, X2, XADMAX, XADMIN, XAFMAX, XAFMIN, XRDMAX
  real ( kind = 4 ) XRDMIN, XRFMAX, XRFMIN, ZERO
  LOGICAL  FAILOC, FAILNX

  real ( kind = 4 ) R1MACH
  real ( kind = 4 ) RAND
  EXTERNAL  RAND
!
!  DEFINE RELATIVE ERROR WITH FLOOR.
!
      REAL  RERR, ERR, VALUE, FLOOR
      RERR(ERR,VALUE,FLOOR) = ERR / MAX(ABS(VALUE), FLOOR)

      DATA  ZERO /0.E0/, ONE /1.E0/, FOUR /4.E0/, TEN /10.E0/
      DATA  SMALL  /1.0E-10/
      DATA  NINT /3/
      DATA   LEFT /-1.5E0, 2.0E-10, 1.0E0 /
      DATA  RIGHT / 2.5E0, 3.0E-10, 1.0E+8/

      MACHEP = R1MACH(4)
      EPS1 = FOUR*MACHEP
      EPS2 = TEN*MACHEP

      FAIL = .FALSE.

      IF (KPRINT .GE. 2)  WRITE (LOUT, 3000)
!
!  CYCLE OVER INTERVALS.
!
      DO IINT = 1, NINT

      X1 =  LEFT(IINT)
      X2 = RIGHT(IINT)

      FACT = MAX(SQRT(X2-X1), ONE)
      TOL1 = EPS1 * FACT
      TOL2 = EPS2 * FACT
!
!  COMPUTE AND PRINT ENDPOINT VALUES.
!
      CALL FDTRUE (X1, F1, D1)
      CALL FDTRUE (X2, F2, D2)

      IF (KPRINT .GE. 3)  THEN
         IF (IINT .EQ. 1)  WRITE (LOUT, 2000)
         WRITE (LOUT, '(/)')
         WRITE (LOUT, 2001)  'X1', X1, 'X2', X2
         WRITE (LOUT, 2001)  'F1', F1, 'F2', F2
         WRITE (LOUT, 2001)  'D1', D1, 'D2', D2
      ENDIF

      IF (KPRINT .GE. 2)  WRITE (LOUT, 3001)  X1, X2
!
!  COMPUTE FLOORS FOR RELATIVE ERRORS.
!
      FLOORF = MAX( MIN(ABS(F1),ABS(F2)), SMALL)
      FLOORD = MAX( MIN(ABS(D1),ABS(D2)), SMALL)
!
!  CHECK REPRODUCTION OF ENDPOINT VALUES.
!
      XEV(1) = X1
      XEV(2) = X2

      CALL CHFDV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECKF, CHECKD, NEXT, IERR)

      AEF  = CHECKF(1)-F1
      REF  = RERR(AEF , F1, FLOORF)
      AEF2 = CHECKF(2)-F2
      REF2 = RERR(AEF2, F2, FLOORF)
      AED  = CHECKD(1)-D1
      RED  = RERR(AED , D1, FLOORD)
      AED2 = CHECKD(2)-D2
      RED2 = RERR(AED2, D2, FLOORD)

      FAILOC = MAX(ABS(REF),ABS(REF2),ABS(RED),ABS(RED2)) .GT. TOL1
      FAIL = FAIL .OR. FAILOC

      IF (KPRINT .GE. 3)  THEN
         WRITE (LOUT, 2002)  NEXT, AEF, AEF2, AED, AED2
         WRITE (LOUT, 2003)  REF, REF2, RED, RED2
      ENDIF

      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3002)
!
!  CHFEV SHOULD AGREE EXACTLY WITH CHFDV.
!
      CALL CHFEV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECK, NEXT, IERR)

      FAILOC = (CHECK(1).NE.CHECKF(1)) .OR. (CHECK(2).NE.CHECKF(2))
      FAIL = FAIL .OR. FAILOC
      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LOUT, 3003)
!
!  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
!  THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
!  TO LEFT AND 6 TO RIGHT OF [X1,X2].
!
      DX = (X2-X1)/(NPTS-10)
      DO I = 1, NPTS
         XEV(I) = (X1 + (I-5)*DX) + DX*RAND(ZERO)
      end do

      CALL CHFDV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV, DEV, NEXT, IERR)

      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 4003)  IERR
      ELSE
!
!  CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
!
      DO I = 1, NPTS

         CALL FDTRUE (XEV(I), FTRUE, DTRUE)
         AEF = FEV(I) - FTRUE
         REF = RERR(AEF, FTRUE, FLOORF)
         AED = DEV(I) - DTRUE
         RED = RERR(AED, DTRUE, FLOORD)

         IF (I .EQ. 1)  THEN
!
!  INITIALIZE.
!
            AEFMIN = AEF
            AEFMAX = AEF
            AEDMIN = AED
            AEDMAX = AED
            REFMIN = REF
            REFMAX = REF
            REDMIN = RED
            REDMAX = RED
            XAFMIN = XEV(1)
            XAFMAX = XEV(1)
            XADMIN = XEV(1)
            XADMAX = XEV(1)
            XRFMIN = XEV(1)
            XRFMAX = XEV(1)
            XRDMIN = XEV(1)
            XRDMAX = XEV(1)
         ELSE
!
!  SELECT.
!
            IF (AEF .LT. AEFMIN)  THEN
               AEFMIN = AEF
               XAFMIN = XEV(I)
            ELSE IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
            IF (AED .LT. AEDMIN)  THEN
               AEDMIN = AED
               XADMIN = XEV(I)
            ELSE IF (AED .GT. AEDMAX)  THEN
               AEDMAX = AED
               XADMAX = XEV(I)
            ENDIF
            IF (REF .LT. REFMIN)  THEN
               REFMIN = REF
               XRFMIN = XEV(I)
            ELSE IF (REF .GT. REFMAX)  THEN
               REFMAX = REF
               XRFMAX = XEV(I)
            ENDIF
            IF (RED .LT. REDMIN)  THEN
               REDMIN = RED
               XRDMIN = XEV(I)
            ELSE IF (RED .GT. REDMAX)  THEN
               REDMAX = RED
               XRDMAX = XEV(I)
            ENDIF
         ENDIF

         end do

         FERMAX = MAX (ABS(REFMAX), ABS(REFMIN))
         DERMAX = MAX (ABS(REDMAX), ABS(REDMIN))

         FAILNX = (NEXT(1) + NEXT(2)) .NE. 10
         FAILOC = FAILNX .OR. (MAX(FERMAX, DERMAX) .GT. TOL2)
      ENDIF
      FAIL = FAIL .OR. FAILOC
!
!  PRINT SUMMARY.
!
      IF (KPRINT .GE. 3)  THEN
         WRITE (LOUT, 2004)  NPTS-10, NEXT
         WRITE (LOUT, 2005)  'MIN', AEFMIN, REFMIN, AEDMIN, REDMIN
         WRITE (LOUT, 2006) XAFMIN, XRFMIN, XADMIN, XRDMIN
         WRITE (LOUT, 2005)  'MAX', AEFMAX, REFMAX, AEDMAX, REDMAX
         WRITE (LOUT, 2006) XAFMAX, XRFMAX, XADMAX, XRDMAX
      ENDIF

      IF (KPRINT .GE. 2)  THEN
         IF (FAILOC) THEN
            IF (FERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'F', FERMAX, TOL2
            IF (DERMAX .GT. TOL2)  WRITE (LOUT, 3006) 'D', DERMAX, TOL2
            IF (FAILNX)  WRITE (LOUT, 4006)  NEXT
         ELSE
            WRITE (LOUT, 5006)
         ENDIF
      ENDIF
!
!  CHECK THAT CHFEV AGREES WITH CHFDV.
!
      CALL CHFEV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV2, NEXT2, IERR)

      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3007)  IERR
      ELSE
         AEFMAX = ABS(FEV2(1) - FEV(1))
         XAFMAX = XEV(1)
         DO I = 2, NPTS
            AEF = ABS(FEV2(I) - FEV(I))
            IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
         end do
         FAILNX = (NEXT2(1).NE.NEXT(1)) .OR. (NEXT2(2).NE.NEXT(2))
         FAILOC = FAILNX .OR. (AEFMAX.NE.ZERO)
         IF (KPRINT .GE. 2)  THEN
            IF (FAILOC)  THEN
               WRITE (LOUT, 3008)
               IF (AEFMAX.NE.ZERO)  WRITE (LOUT, 3009)  AEFMAX, XAFMAX
               IF (FAILNX)  WRITE (LOUT, 4009)  NEXT2, NEXT
            ELSE
               WRITE (LOUT, 5009)
            ENDIF
         ENDIF
      ENDIF

      FAIL = FAIL .OR. FAILOC
!
!  GO BACK FOR ANOTHER INTERVAL.
!
      end do

      RETURN

 2000 FORMAT (/10X,'CHFDV ACCURACY TEST')
 2001 FORMAT (10X,A2,' =',1P,E18.10,5X,A2,' =',E18.10)
 2002 FORMAT (/' ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')' &
             // 1P,4X,'F1:',E13.5,4X,'F2:',E13.5, &
                  4X,'D1:',E13.5,4X,'D2:',E13.5)
 2003 FORMAT (1P,4(7X,E13.5))
 2004 FORMAT (/' ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:', &
                     15X,'(NEXT =',2I3,')' &
             //30X,'FUNCTION',17X,'DERIVATIVE' &
              /15X,2(11X,'ABS',9X,'REL') )
 2005 FORMAT (/5X,A3,'IMUM ERROR:  ',1P,2E12.4,2X,2E12.4)
 2006 FORMAT ( 5X,'LOCATED AT X =  ',1P,2E12.4,2X,2E12.4)
 3000 FORMAT (//10X,'EVCHCK RESULTS'/10X,'--------------')
 3001 FORMAT (/10X,'INTERVAL = (',1P,E12.5,',',E12.5,' ):' )
 3002 FORMAT (/' ***** CHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
 3003 FORMAT (/' ***** CHFEV DOES NOT AGREE WITH CHFDV AT ENDPOINTS.')
 3006 FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',1P,E12.5,',' &
             /        17X,'EXCEEDS TOLERANCE =',E12.5)
 3007 FORMAT (/' ***** ERROR ***** CHFEV RETURNED IERR =',I5)
 3008 FORMAT (/' ***** CHFEV DID NOT AGREE WITH CHFDV:')
 3009 FORMAT (7X,'MAXIMUM DIFFERENCE ',1P,E12.5, &
                    '; OCCURRED AT X =',E12.5)
 4003 FORMAT (/' ***** ERROR ***** CHFDV RETURNED IERR =',I5)
 4006 FORMAT (/' ***** REPORTED NEXT =',2I5,'   RATHER THAN    4    6')
 4009 FORMAT (7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
 5006 FORMAT (/' CHFDV RESULTS OK.')
 5009 FORMAT (/' CHFEV AGREES WITH CHFDV.')

END
SUBROUTINE EVERCK (LOUT, KPRINT, FAIL)

!*****************************************************************************80
!
!! EVERCK tests error returns from PCHIP evaluators for PCHQK1.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    LOUT
!
!    KPRINT
!
!    FAIL
!
  implicit none

  integer ( kind = 4 ) LOUT, KPRINT
  LOGICAL  FAIL
  integer ( kind = 4 ) I, IERR, KONTRL, N, NERR, NEXT(2)
  real ( kind = 4 ) D(10), DUM, F(10), TEMP, X(10)
  LOGICAL  COMP, SKIP
  PARAMETER (N = 10)

  NERR = 0

  CALL XGETF (KONTRL)
  IF (KPRINT .LE. 2) THEN
    CALL XSETF (0)
  ELSE
    CALL XSETF (1)
  ENDIF

  IF (KPRINT .GE. 3)  WRITE (LOUT, 2000)
  IF (KPRINT .GE. 2)  WRITE (LOUT, 5000)
!
!  FIRST, TEST CHFEV AND CHFDV.
!
  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  CALL CHFEV (0.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 0, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
  CALL CHFEV (1.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 1, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
  CALL CHFDV (0.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 0, DUM, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

  IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -2
  CALL CHFDV (1.E0, 1.E0, 3.E0, 7.E0, 3.E0, 6.E0, 1, DUM, DUM, DUM, NEXT, IERR)
  IF (.NOT. COMP (IERR, -2, LOUT, KPRINT) )  NERR = NERR + 1
!
!  SET UP PCH DEFINITION.
!
      DO I = 1, N
         X(I) = I
         F(I) = I + 2
         D(I) = 1.E0
      end do
!
!  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
!
      TEMP = X(4)
      X(4) = X(7)
      X(7) = TEMP
!
!  NOW, TEST PCHFE AND PCHFD.
!
      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      SKIP = .FALSE.
      CALL PCHFE (1, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
      SKIP = .FALSE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1

      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
      SKIP = .TRUE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1

      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -1
      SKIP = .FALSE.
      CALL PCHFD (1, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LOUT, KPRINT) )  NERR = NERR + 1

      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -3
      SKIP = .FALSE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LOUT, KPRINT) )  NERR = NERR + 1

      IF (KPRINT .GE. 3)  WRITE (LOUT, 5001)  -4
      SKIP = .TRUE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LOUT, KPRINT) )  NERR = NERR + 1
!
!  SUMMARIZE RESULTS.
!
      IF (KPRINT .GT. 2)  CALL XERDMP
      IF (NERR .EQ. 0)  THEN
         FAIL = .FALSE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 5002)
      ELSE
         FAIL = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 5003)  NERR
      ENDIF

  CALL XSETF (KONTRL)
  RETURN

 2000 FORMAT ('1'//10X,'TEST ERROR RETURNS')
 5000 FORMAT (//10X,'EVERCK RESULTS'/10X,'--------------')
 5001 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
 5002 FORMAT (/' ALL ERROR RETURNS OK.')
 5003 FORMAT (//' ***** TROUBLE IN EVERCK *****' &
             //5X,I5,' TESTS FAILED TO GIVE EXPECTED RESULTS.')

END
SUBROUTINE EVPCCK (LOUT, KPRINT, X, Y, F, FX, FY, XE, YE, FE, DE, FE2, FAIL)

!*****************************************************************************80
!
!! EVPCCK tests usage of increment argument in PCHFD and PCHFE for PCHQK1.
!
!  Discussion:
!
!    EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
!    ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
!
!    INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
!    SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
!
!    ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
!    TEST ROUTINES.
!
!    RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
!   
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
  implicit none

  integer ( kind = 4 ) LOUT, KPRINT
  LOGICAL  FAIL
  real ( kind = 4 ) X(10), Y(10), F(10,10), FX(10,10), FY(10,10)
  real ( kind = 4 ) XE(51), YE(51), FE(51), DE(51), FE2(51)
  integer ( kind = 4 ) I, IER2, IERR, INC, J, K, NE, NERR, NMAX, NX, NY
  LOGICAL  FAILD, FAILE, FAILOC, SKIP
  real ( kind = 4 ) DERMAX, DERR, DTRUE, DX, FDIFF, FDIFMX, FERMAX, FERR
  real ( kind = 4 ) FTRUE, MACHEP, TOL, PDERMX, PDIFMX, PFERMX, ZERO
  real ( kind = 4 ) R1MACH
!
!  DEFINE TEST FUNCTION AND DERIVATIVES.
!
  real ( kind = 4 ) AX, AY, FCN, DFDX, DFDY
  FCN(AX,AY)  =  AX*(AY*AY)*(AX*AX + 1.E0)
  DFDX(AX,AY) = (AY*AY)*(3.E0*AX*AX + 1.E0)
  DFDY(AX,AY) =   2.E0*AX*AY*(AX*AX + 1.E0)

  DATA  NMAX /10/,  NX /4/,  NY /6/
  DATA  NE /51/
  DATA  ZERO /0.E0/

  MACHEP = R1MACH(4)
  TOL = 10.E0*MACHEP

      FAIL = .FALSE.
!
!  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
!  X =  0.25(0.25)1.   ;
!  Y = -0.75(0.5 )1.75 .
!
      DO I = 1, NX-1
         X(I) = 0.25E0*I
      end do
      X(NX) = 1.E0

      DO J = 1, NY
         Y(J) = 0.5E0*J - 1.25E0
         DO I = 1, NX
             F(I,J) = FCN (X(I), Y(J))
            FX(I,J) = DFDX(X(I), Y(J))
            FY(I,J) = DFDY(X(I), Y(J))
         end do
      end do
!
!  SET UP EVALUATION POINTS:
!  XE =  0.(0.02)1. ;
!  YE = -2.(0.08)2. .
!
      DX = 1.E0/(NE-1)
      DO K = 1, NE-1
         XE(K) = DX*(K-1)
         YE(K) = 4.E0*XE(K) - 2.E0
      end do
      XE(NE) = 1.E0
      YE(NE) = 2.E0

      IF (KPRINT .GE. 3)  WRITE (LOUT, 1000)
      IF (KPRINT .GE. 2)  WRITE (LOUT, 1001)
!
!  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING)
!
      NERR = 0
      INC = 1
      SKIP = .FALSE.

      DO J = 1, NY

         CALL PCHFD (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE, DE, IERR)

         IF (KPRINT .GE. 3) then
            WRITE (LOUT, 2000)  INC, 'J', J, 'Y', Y(J), IERR
         end if

         IF (IERR .LT. 0)  GO TO 15
         IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'X'
!
!  PCHFE SHOULD AGREE EXACTLY WITH PCHFD.
!
         CALL PCHFE (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE2, IER2)

         DO K = 1, NE

            FTRUE =  FCN(XE(K), Y(J))
            FERR = FE(K) - FTRUE
            DTRUE = DFDX(XE(K), Y(J))
            DERR = DE(K) - DTRUE
            IF (KPRINT .GT. 3) then
               WRITE (LOUT, 2002)  XE(K), FTRUE, FE(K), FERR, &
                                          DTRUE, DE(K), DERR
            end if

            IF (K .EQ. 1)  THEN
!
!  INITIALIZE.
!
               FERMAX = ABS(FERR)
               PFERMX = XE(1)
               DERMAX = ABS(DERR)
               PDERMX = XE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = XE(1)
            ELSE
!
!  SELECT.
!
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = XE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = XE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = XE(K)
               ENDIF
            ENDIF

         end do

         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.13) .OR. (IER2.NE.IERR)

         IF (FAILOC .AND. (KPRINT.GE.2)) then
            WRITE (LOUT, 2003)  'J', J, 'Y', Y(J)
         end if

         IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) ) then
            WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
         end if

         IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL

         IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) ) then
           WRITE (LOUT, 2005)  FDIFMX, PDIFMX
         end if

         IF ((IERR.NE.13) .AND. (KPRINT.GE.2)) then
           WRITE (LOUT, 2006)  'D', IERR, 13
         end if

         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2)) then
           WRITE (LOUT, 2006)  'E', IER2, IERR
         end if

         GO TO 19

   15    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR

   19    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC

      end do

      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LOUT, 3001)  NERR, 'J'
         ELSE
            WRITE (LOUT, 4000)  'J'
         ENDIF
      ENDIF
!
!  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING).
!
      NERR = 0
      INC = NMAX
      SKIP = .FALSE.

      DO I = 1, NX

         CALL PCHFD (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE, DE, IERR)

         IF (KPRINT .GE. 3) then
             WRITE (LOUT, 2000)  INC, 'I', I, 'X', X(I), IERR
         end if

         IF (IERR .LT. 0)  GO TO 35
         IF (KPRINT .GT. 3)  WRITE (LOUT, 2001)  'Y'
!
!  PCHFE SHOULD AGREE EXACTLY WITH PCHFD.
!
         CALL PCHFE (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE2, IER2)

         DO K = 1, NE

            FTRUE =  FCN(X(I), YE(K))
            FERR = FE(K) - FTRUE
            DTRUE = DFDY(X(I), YE(K))
            DERR = DE(K) - DTRUE
            IF (KPRINT .GT. 3) then
               WRITE (LOUT, 2002)  YE(K), FTRUE, FE(K), FERR, &
                                          DTRUE, DE(K), DERR
            end if

            IF (K .EQ. 1)  THEN
!
!  INITIALIZE.
!
               FERMAX = ABS(FERR)
               PFERMX = YE(1)
               DERMAX = ABS(DERR)
               PDERMX = YE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = YE(1)
            ELSE
!
!  SELECT.
!
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = YE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = YE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = YE(K)
               ENDIF
            ENDIF

         end do

         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.20) .OR. (IER2.NE.IERR)

         IF (FAILOC .AND. (KPRINT.GE.2)) then
           WRITE (LOUT, 2003)  'I', I, 'X', X(I)
         end if

         IF ((KPRINT.GE.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) ) then
           WRITE (LOUT, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
         end if

         IF (FAILD .AND. (KPRINT.GE.2))  WRITE (LOUT, 2014)  TOL

         IF ((KPRINT.GE.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) ) then
           WRITE (LOUT, 2005)  FDIFMX, PDIFMX
         end if

         IF ((IERR.NE.20) .AND. (KPRINT.GE.2)) then
           WRITE (LOUT, 2006)  'D', IERR, 20
         end if

         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2)) then
           WRITE (LOUT, 2006)  'E', IER2, IERR
         end if

         GO TO 39

   35    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LOUT, 3000) IERR

   39    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC

      end do

      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LOUT, 3001)  NERR, 'I'
         ELSE
            WRITE (LOUT, 4000)  'I'
         ENDIF
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCHFE AND PCHFD')
 1001 FORMAT (//10X,'EVPCCK RESULTS'/10X,'--------------')
 2000 FORMAT (//20X,'PCHFD INCREMENT TEST -- INCFD = ',I2 &
             /15X,'ON ',A1,'-LINE ',I2,',  ',A1,' =',F8.4, &
                '  --  IERR =',I3)
 2001 FORMAT ( /3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF', &
                         13X,'D',8X,'DE',9X,'DIFF')
 2002 FORMAT (F7.2,2(2X,2F10.5,1P,E15.5,0P))
 2003 FORMAT (/' ***** PCHFD AND/OR PCHFE FAILED ON ',A1,'-LINE ',I1, &
                                  ',  ',A1,' =',F8.4)
 2004 FORMAT (/17X,'  MAXIMUM ERROR IN FUNCTION =',1P, &
                                        1P,E13.5,0P,' (AT',F6.2,'),' &
             /31X,    'IN DERIVATIVE =',1P,E13.5,0P,' (AT',F6.2,').' )
 2005 FORMAT ( '  MAXIMUM DIFFERENCE BETWEEN PCHFE AND PCHFD =', &
                                        1P,E13.5,0P,' (AT',F6.2,').' )
 2006 FORMAT (/'  PCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
 2014 FORMAT ('  *** BOTH SHOULD BE .LE. TOL =',1P,E12.5,' ***')
 3000 FORMAT (//' ***** ERROR ***** PCHFD RETURNED IERR =',I5//)
 3001 FORMAT (//' ***** ERROR ***** PCHFD AND/OR PCHFE FAILED ON',I2, &
                                     1X,A1,'-LINES.'//)
 4000 FORMAT (/' PCHFD AND PCHFE OK ON ',A1,'-LINES.')
END
SUBROUTINE FDTRUE (X, F, D)

!*****************************************************************************80
!
!! FDTRUE computes exact function values for EVCHCK.
!
!  Discussion:
!
!    F(X) = X * ( X + 1 ) * ( X - 2 ).
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the evaluation point.
!
!    Output, real ( kind = 4 ) F, D, the values of the function and
!    derivative at X.
!
  implicit none

  real ( kind = 4 ) X, F, D
  real ( kind = 8 ) FACT1, FACT2, XX

  XX = real ( X, kind = 8 )
  FACT1 = XX + 1.0D+00
  FACT2 = XX - 2.0D+00
  F = real ( XX * FACT1 * FACT2, kind = 4 )
  D = real ( FACT1*FACT2 + XX*(FACT1 + FACT2), kind = 4 )

  RETURN
END
SUBROUTINE FDUMP ( )

!*****************************************************************************80
!
!! FDUMP can be used to print a symbolic dump.
!
!  Discussion:
!
!    FDUMP is intended to be replaced by a locally written
!    version which produces a symbolic dump.  Failing this,
!    it should be replaced by a version which prints the
!    subprogram nesting list.  Note that this dump must be
!    printed on each of up to five files, as indicated by the
!    XGETUA routine.  See XSETUA and XGETUA for details.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
  implicit none

  RETURN
END
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical ( kind = 4 ) lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i1mach ( i )

!*****************************************************************************80
!
!! I1MACH returns integer machine constants.
!
!  Discussion:
!
!    Input/output unit numbers.
!
!      I1MACH(1) = the standard input unit.
!      I1MACH(2) = the standard output unit.
!      I1MACH(3) = the standard punch unit.
!      I1MACH(4) = the standard error message unit.
!
!    Words.
!
!      I1MACH(5) = the number of bits per integer storage unit.
!      I1MACH(6) = the number of characters per integer storage unit.
!
!    Integers.
!
!    Assume integer ( kind = 4 )s are represented in the S digit base A form:
!
!      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
!
!    where 0 <= X(1:S-1) < A.
!
!      I1MACH(7) = A, the base.
!      I1MACH(8) = S, the number of base A digits.
!      I1MACH(9) = A^S-1, the largest integer.
!
!    Floating point numbers
!
!    Assume floating point numbers are represented in the T digit 
!    base B form:
!
!      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
!
!    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
!
!      I1MACH(10) = B, the base.
!
!    Single precision
!
!      I1MACH(11) = T, the number of base B digits.
!      I1MACH(12) = EMIN, the smallest exponent E.
!      I1MACH(13) = EMAX, the largest exponent E.
!
!    real ( kind = 8 )
!
!      I1MACH(14) = T, the number of base B digits.
!      I1MACH(15) = EMIN, the smallest exponent E.
!      I1MACH(16) = EMAX, the largest exponent E.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 16.
!
!    Output, integer ( kind = 4 ) I1MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
    stop
  else if ( i == 1 ) then
    i1mach = 5
  else if ( i == 2 ) then
    i1mach = 6
  else if ( i == 3 ) then
    i1mach = 7
  else if ( i == 4 ) then
    i1mach = 6
  else if ( i == 5 ) then
    i1mach = 32
  else if ( i == 6 ) then
    i1mach = 4
  else if ( i == 7 ) then
    i1mach = 2
  else if ( i == 8 ) then
    i1mach = 31
  else if ( i == 9 ) then
    i1mach = 2147483647
  else if ( i == 10 ) then
    i1mach = 2
  else if ( i == 11 ) then
    i1mach = 24
  else if ( i == 12 ) then
    i1mach = -125
  else if ( i == 13 ) then
    i1mach = 128
  else if ( i == 14 ) then
    i1mach = 53
  else if ( i == 15 ) then
    i1mach = -1021
  else if ( i == 16 ) then
    i1mach = 1024
  else if ( 16 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
    stop
  end if

  return
end
SUBROUTINE INTRV (XT, LXT, X, ILO, ILEFT, MFLAG)

!*****************************************************************************80
!
!! INTRV seeks the interval containing or nearest to a given point.
!
!  Discussion:
!
!    Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
!    such that XT(ILEFT) .LE. X where XT(*) is a subdivision of the X interval.
!
!    INTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!    LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!    the X interval.  Precisely,
!
!                      X .LT. XT(1)                1         -1
!       if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!         XT(LXT) .LE. X                         LXT        1,
!
!    That is, when multiplicities are present in the break point
!    to the left of X, the largest index is taken for ILEFT.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Daniel Amos
!
!  Reference:
!
!    Carl de Boor, 
!    Package for calculating with B-splines,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 3, June 1977, pages 441-472.
!
!  Parameters:
!
!       Input      XT,X are real
!        XT      - XT is a knot or break point vector of length LXT
!        LXT     - length of the XT vector
!        X       - argument
!        ILO     - an initialization parameter which must be set
!                  to 1 the first time the spline array XT is
!                  processed by DINTRV.
!
!       Output
!        ILO     - ILO contains information for efficient process-
!                  ing after the initial call and ILO must not be
!                  changed by the user.  Distinct splines require
!                  distinct ILO parameters.
!        ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!        MFLAG   - signals when X lies out of bounds
!
  implicit none

  integer ( kind = 4 ) IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
  real ( kind = 4 ) X, XT(*)

      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT

   10 continue

      IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
!
!  NOW X .LT. XT(IHI) . FIND LOWER BOUND.
!
      ISTEP = 1
   20 continue
      IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 continue
      ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
!
!  NOW X .GE. XT(ILO) . FIND UPPER BOUND.
!
   40 continue
      ISTEP = 1
   50 continue
      ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 continue
      IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
!
!  NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL.
!
   70 continue
      MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
!
!  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1.
!
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
!
!  SET OUTPUT AND RETURN.
!
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT

  RETURN
END
FUNCTION J4SAVE ( IWHICH, IVALUE, ISET )

!*****************************************************************************80
!
!! J4SAVE saves or recalls global variables needed by the error handler.
!
!  Discussion:
!
!    J4SAVE saves and recalls several global variables needed
!    by the library error handling routines.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
  implicit none

  integer ( kind = 4 ) IPARAM(9)
  logical iset
  integer ( kind = 4 ) ivalue
  integer ( kind = 4 ) iwhich
  integer ( kind = 4 ) j4save

  SAVE IPARAM

  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/

  J4SAVE = IPARAM(IWHICH)

  IF ( ISET ) then
    IPARAM(IWHICH) = IVALUE
  end if

  RETURN
END
SUBROUTINE PCHBS ( N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF, NDIM, KORD, &
  IERR )

!*****************************************************************************80
!
!! PCHBS: Piecewise cubic Hermite to B-Spline converter.
!
!  Discussion:
!
!    PCHBS computes the B-spline representation of the PCH function
!    determined by N,X,F,D.  To be compatible with the rest of PCHIP,
!    PCHBS includes INCFD, the increment between successive values of
!    the F- and D-arrays.
!
!    The output is the B-representation for the function:  NKNOTS, T,
!    BCOEF, NDIM, KORD.
!
!    Since it is assumed that the input PCH function has been
!    computed by one of the other routines in the package PCHIP,
!    input arguments N, X, INCFD are **not** checked for validity.
!
!    Restrictions/assumptions include:
!     1. N.GE.2 .  (not checked)
!     2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
!     3. INCFD.GT.0 .  (not checked)
!     4. KNOTYP.LE.2 .  (error return if not)
!    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
!    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
!
!    * Indicates this applies only if KNOTYP.LT.0 .
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch,
!    Representations for Parametric Cubic Splines,
!    Computer Aided Geometric Design,
!    Volume 6, 1989, pages 79-82.
!
! *Usage:
!
!        integer ( kind = 4 )  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
!        PARAMETER  (INCFD = ...)
!        REAL  X(nmax), F(INCFD,nmax), D(INCFD,nmax), T(2*nmax+4),
!       *      BCOEF(2*nmax)
!
!        CALL  PCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
!       *             NDIM, KORD, IERR)
!
!  Parameters:
!
!     N:IN  is the number of data points, N.ge.2 .  (not checked)
!
!     X:IN  is the real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.   (not checked)
!           nmax, the dimension of X, must be .ge.N.
!
!     F:IN  is the real array of dependent variable values.
!           F(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of F, must be .ge.N.
!
!     D:IN  is the real array of derivative values at the data points.
!           D(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of D, must be .ge.N.
!
!     INCFD:IN  is the increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           It may have the value 1 for one-dimensional applications,
!           in which case F and D may be singly-subscripted arrays.
!
!     KNOTYP:IN  is a flag to control the knot sequence.
!           The knot sequence T is normally computed from X by putting
!           a double knot at each X and setting the end knot pairs
!           according to the value of KNOTYP:
!              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
!              KNOTYP = 1:  Replicate lengths of extreme subintervals:
!                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
!                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
!              KNOTYP = 2:  Periodic placement of boundary knots:
!                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
!                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
!              Here M=NDIM=2*N.
!           If the input value of KNOTYP is negative, however, it is
!           assumed that NKNOTS and T were set in a previous call.
!           This option is provided for improved efficiency when used
!           in a parametric setting.
!
!     NKNOTS:INOUT  is the number of knots.
!           If KNOTYP.GE.0, then NKNOTS will be set to NDIM+4.
!           If KNOTYP.LT.0, then NKNOTS is an input variable, and an
!              error return will be taken if it is not equal to NDIM+4.
!
!     T:INOUT  is the array of 2*N+4 knots for the B-representation.
!           If KNOTYP.GE.0, T will be returned by PCHBS with the
!              interior double knots equal to the X-values and the
!              boundary knots set as indicated above.
!           If KNOTYP.LT.0, it is assumed that T was set by a
!              previous call to PCHBS.  (This routine does **not**
!              verify that T forms a legitimate knot sequence.)
!
!     BCOEF:OUT  is the array of 2*N B-spline coefficients.
!
!     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
!
!     KORD:OUT  is the order of the B-spline.  (Set to 4.)
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -4  if KNOTYP.GT.2 .
!              IERR = -5  if KNOTYP.LT.0 and NKNOTS.NE.(2*N+4).
!
  implicit none

  integer ( kind = 4 )  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*), T(*), BCOEF(*)
  integer ( kind = 4 )  K, KK
  real ( kind = 4 ) DOV3, HNEW, HOLD
  CHARACTER*8  LIBNAM, SUBNAM
!
!  Initialize.
!
  NDIM = 2*N
  KORD = 4
  IERR = 0
  LIBNAM = 'SLATEC'
  SUBNAM = 'PCHBS'
!
!  Check argument validity.  Set up knot sequence if OK.
!
  IF ( KNOTYP.GT.2 )  THEN
    IERR = -1
    CALL XERMSG (LIBNAM, SUBNAM, 'KNOTYP GREATER THAN 2', IERR, 1)
    RETURN
  end if

  IF ( KNOTYP.LT.0 )  THEN

    IF ( NKNOTS.NE.NDIM+4 )  THEN
      IERR = -2
      CALL XERMSG (LIBNAM, SUBNAM, &
        'KNOTYP.LT.0 AND NKNOTS.NE.(2*N+4)', IERR, 1)
      RETURN
    end if

  ELSE
!
!  Set up knot sequence.
!
   NKNOTS = NDIM + 4
   CALL PCHKT (N, X, KNOTYP, T)

  end if
!
!  Compute B-spline coefficients.
!
  HNEW = T(3) - T(1)
  DO K = 1, N
     KK = 2*K
     HOLD = HNEW
!
!  The following requires mixed mode arithmetic.
!
     DOV3 = D(1,K)/3
     BCOEF(KK-1) = F(1,K) - HOLD*DOV3
!
!  The following assumes T(2*K+1) = X(K).
!
     HNEW = T(KK+3) - T(KK+1)
     BCOEF(KK) = F(1,K) + HNEW*DOV3

  end do

  RETURN
END
SUBROUTINE PCHCE ( IC, VC, N, X, H, SLOPE, D, INCFD, IERR )

!*****************************************************************************80
!
!! PCHCE sets boundary conditions for PCHIC.
!
!  Discussion:
!
!    Called by PCHIC to set end derivatives as requested by the user.
!    It must be called after interior derivative values have been set.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, IERR
!        REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!           ( see prologue to PCHIC for details. )
!
!     VC -- (input) real array of length 2 specifying desired boundary
!           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!                    VC(2) need be set only if IC(2) = 2 or 3 .
!
!     N -- (input) number of data points.  (assumes N.GE.2)
!
!     X -- (input) real array of independent variable values.  (the
!           elements of X are assumed to be strictly increasing.)
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real array of derivative values at the data points.
!           The value corresponding to X(I) must be stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!          (output) the value of D at X(1) and/or X(N) is changed, if
!           necessary, to produce the requested boundary conditions.
!           no other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!
  implicit none

  integer ( kind = 4 )  IC(2), N, INCFD, IERR
  real ( kind = 4 ) VC(2), X(*), H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  IBEG, IEND, IERF, INDEX, J, K
  real ( kind = 4 ) HALF, STEMP(3), XTEMP(4)
  SAVE HALF
  real ( kind = 4 ) PCHDF, PCHST
!
!  INITIALIZE.
!
  DATA  HALF /0.5/

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
  IF ( ABS(IBEG).GT.N )  IBEG = 0
  IF ( ABS(IEND).GT.N )  IEND = 0
!
!  TREAT BEGINNING BOUNDARY CONDITION.
!
  IF (IBEG == 0)  GO TO 2000
  K = ABS(IBEG)
  IF (K == 1)  THEN
!
!  BOUNDARY VALUE PROVIDED.
!
     D(1,1) = VC(1)
  ELSE IF (K == 2)  THEN
!
!  BOUNDARY SECOND DERIVATIVE PROVIDED.
!
     D(1,1) = HALF*( ( 3.0E+00 *SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
  ELSE IF (K .LT. 5)  THEN
!
!  USE K-POINT DERIVATIVE FORMULA.
!  PICK UP FIRST K POINTS, IN REVERSE ORDER.
!
     DO J = 1, K
        INDEX = K-J+1
!
!  INDEX RUNS FROM K DOWN TO 1.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
    end do

     D(1,1) = PCHDF (K, XTEMP, STEMP, IERF)

     IF (IERF .NE. 0)  GO TO 5001
  ELSE
!
!  USE 'NOT A KNOT' CONDITION.
!
     D(1,1) = ( 3.0E+00 *(H(1)*SLOPE(2) + H(2)*SLOPE(1)) &
               - 2.0E+00 *(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
  end if
!
  IF (IBEG .GT. 0)  GO TO 2000
!
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!
  IF (SLOPE(1) == 0.0E+00 )  THEN
     IF (D(1,1) .NE. 0.0E+00 )  THEN
        D(1,1) = 0.0E+00
        IERR = IERR + 1
     end if
  ELSE IF ( PCHST(D(1,1),SLOPE(1)) .LT. 0.0E+00 )  THEN
     D(1,1) = 0.0E+00
     IERR = IERR + 1
  ELSE IF ( ABS(D(1,1)) .GT. 3.0E+00 *ABS(SLOPE(1)) )  THEN
     D(1,1) = 3.0E+00 *SLOPE(1)
     IERR = IERR + 1
  end if
!
!  TREAT END BOUNDARY CONDITION.
!
 2000 CONTINUE

  IF (IEND == 0) then
    return
  end if

  K = ABS(IEND)
  IF (K == 1)  THEN
!
!  BOUNDARY VALUE PROVIDED.
!
     D(1,N) = VC(2)
  ELSE IF (K == 2)  THEN
!
!  BOUNDARY SECOND DERIVATIVE PROVIDED.
!
     D(1,N) = HALF*( ( 3.0E+00 *SLOPE(N-1) - D(1,N-1)) + &
                                            HALF*VC(2)*H(N-1) )
  ELSE IF (K .LT. 5)  THEN
!
!  USE K-POINT DERIVATIVE FORMULA.  PICK UP LAST K POINTS.
!
     DO J = 1, K
        INDEX = N-K+J
!
!  INDEX RUNS FROM N+1-K UP TO N.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
     end do

     D(1,N) = PCHDF (K, XTEMP, STEMP, IERF)

     IF (IERF .NE. 0)  GO TO 5001
  ELSE
!
!  USE 'NOT A KNOT' CONDITION.
!
     D(1,N) = ( 3.0E+00 *(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1)) &
               - 2.0E+00 * (H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) ) &
                                                           / H(N-2)
  end if

  IF (IEND .GT. 0) then
    return
  end if
!
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!
  IF (SLOPE(N-1) == 0.0E+00)  THEN
     IF (D(1,N) .NE. 0.0E+00)  THEN
        D(1,N) = 0.0E+00
        IERR = IERR + 2
     end if
  ELSE IF ( PCHST(D(1,N),SLOPE(N-1)) .LT. 0.0E+00)  THEN
     D(1,N) = 0.0E+00
     IERR = IERR + 2
  ELSE IF ( ABS(D(1,N)) .GT. 3.0E+00 *ABS(SLOPE(N-1)) )  THEN
     D(1,N) = 3.0E+00 *SLOPE(N-1)
     IERR = IERR + 2
  end if

  RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!     ERROR RETURN FROM PCHDF.
!  THIS CASE SHOULD NEVER OCCUR
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHCE', 'ERROR RETURN FROM PCHDF', IERR, 1)

  RETURN
END
SUBROUTINE PCHCI ( N, H, SLOPE, D, INCFD )

!*****************************************************************************80
!
!! PCHCI sets interior derivatives for PCHIC.
!
!  Discussion:
!
!    Called by PCHIC to set derivatives needed to determine a monotone
!    piecewise cubic Hermite interpolant to the data.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by PCHIM.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N
!        REAL  H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  PCHCI (N, H, SLOPE, D, INCFD)
!
!  Parameters:
!
!     N -- (input) number of data points.
!           If N=2, simply does linear interpolation.
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
  implicit none

  integer ( kind = 4 )  N, INCFD
  real ( kind = 4 ) H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  I, NLESS1
  real ( kind = 4 ) DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM, HSUMT3
  real ( kind = 4 ) W1, W2
  real ( kind = 4 ) PCHST

  NLESS1 = N - 1
  DEL1 = SLOPE(1)
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF ( NLESS1 .le. 1) then
    D(1,1) = DEL1
    D(1,N) = DEL1
    return
  end if
!
!  NORMAL CASE  (N .GE. 3).
!
  DEL2 = SLOPE(2)
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  HSUM = H(1) + H(2)
  W1 = (H(1) + HSUM)/HSUM
  W2 = -H(1)/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2
  IF ( PCHST(D(1,1),DEL1) .LE. 0.0E+00)  THEN
     D(1,1) = 0.0E+00
  ELSE IF ( PCHST(DEL1,DEL2) .LT. 0.0E+00)  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0E+00 *DEL1
     IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO I = 2, NLESS1

     IF ( 2 < I ) then
       HSUM = H(I-1) + H(I)
       DEL1 = DEL2
       DEL2 = SLOPE(I)
     end if
!
!  SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = 0.0E+00
     IF ( PCHST(DEL1,DEL2) .LE. 0.0E+00)  GO TO 50
!
!  USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H(I-1))/HSUMT3
     W2 = (HSUM + H(I)  )/HSUMT3
     DMAX = MAX( ABS(DEL1), ABS(DEL2) )
     DMIN = MIN( ABS(DEL1), ABS(DEL2) )
     DRAT1 = DEL1/DMAX
     DRAT2 = DEL2/DMAX
     D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)

   50    CONTINUE

  end do
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  W1 = -H(N-1)/HSUM
  W2 = (H(N-1) + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2
  IF ( PCHST(D(1,N),DEL2) .LE. 0.0E+00)  THEN
     D(1,N) = 0.0E+00
  ELSE IF ( PCHST(DEL1,DEL2) .LT. 0.0E+00)  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0E+00 *DEL2
     IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX

  end if

  RETURN
END
SUBROUTINE PCHCM ( N, X, F, D, INCFD, SKIP, ISMON, IERR )

!*****************************************************************************80
!
!! PCHCM checks a cubic Hermite function for monotonicity.
!
!  Discussion:
!
!    Checks the piecewise cubic Hermite function defined by  N,X,F,D
!    for monotonicity.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!***DESCRIPTION
!
! *Usage:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, ISMON(N), IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        CALL  PCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
!
! *Arguments:
!
!     N:IN  is the number of data points.  (Error return if N.LT.2 .)
!
!     X:IN  is a real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F:IN  is a real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D:IN  is a real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD:IN  is the increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP:INOUT  is a logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed.
!           SKIP will be set to .TRUE. on normal return.
!
!     ISMON:OUT  is an integer ( kind = 4 ) array indicating on which 
!     intervals the PCH function defined by  N, X, F, D  is monotonic.
!           For data interval [X(I),X(I+1)],
!             ISMON(I) = -3  if function is probably decreasing;
!             ISMON(I) = -1  if function is strictly decreasing;
!             ISMON(I) =  0  if function is constant;
!             ISMON(I) =  1  if function is strictly increasing;
!             ISMON(I) =  2  if function is non-monotonic;
!             ISMON(I) =  3  if function is probably increasing.
!                If ABS(ISMON)=3, this means that the D-values are near
!                the boundary of the monotonicity region.  A small
!                increase produces non-monotonicity; decrease, strict
!                monotonicity.
!           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
!              the entire function is monotonic on [X(1),X(N)].
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!          (The ISMON-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, ISMON(N), IERR
  real ( kind = 4 ) X(N), F(INCFD,N), D(INCFD,N)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, NSEG
  real ( kind = 4 ) DELTA
  integer ( kind = 4 )  CHFCM
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  SKIP = .TRUE.
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  NSEG = N - 1

  DO I = 1, NSEG

     DELTA = (F(1,I+1)-F(1,I))/(X(I+1)-X(I))

     ISMON(I) = CHFCM (D(1,I), D(1,I+1), DELTA)

     IF (I == 1)  THEN
        ISMON(N) = ISMON(1)
     ELSE
!
!  Need to figure out cumulative monotonicity from following
!  "multiplication table":
!
!                    +        I S M O N (I)
!                     +  -3  -1   0   1   3   2
!                      +------------------------+
!               I   -3 I -3  -3  -3   2   2   2 I
!               S   -1 I -3  -1  -1   2   2   2 I
!               M    0 I -3  -1   0   1   3   2 I
!               O    1 I  2   2   1   1   3   2 I
!               N    3 I  2   2   3   3   3   2 I
!              (N)   2 I  2   2   2   2   2   2 I
!                      +------------------------+
!  Note that the 2 row and column are out of order so as not
!  to obscure the symmetry in the rest of the table.
!
!  No change needed if equal or constant on this interval or
!  already declared nonmonotonic.
!
        IF ( (ISMON(I).NE.ISMON(N)) .AND. (ISMON(I).NE.0) &
                                   .AND. (ISMON(N).NE.2) )  THEN
           IF ( (ISMON(I) == 2) .OR. (ISMON(N) == 0) )  THEN
              ISMON(N) =  ISMON(I)
           ELSE IF (ISMON(I)*ISMON(N) .LT. 0)  THEN
!
!  This interval has opposite sense from curve so far.
!
              ISMON(N) = 2
           ELSE
!
!  At this point, both are nonzero with same sign, and
!  we have already eliminated case both +-1.
!
              ISMON(N) = ISIGN (3, ISMON(N))
           end if
        end if
     end if

  end do

  IERR = 0
  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHCM', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHCM', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHCM', 'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)

  RETURN
END
SUBROUTINE PCHCS ( SWITCH, N, H, SLOPE, D, INCFD, IERR )

!*****************************************************************************80
!
!! PCHCS adjusts derivative values for PCHIC.
!
!  Discussion:
!
!    Called by PCHIC to adjust the values of D in the vicinity of a
!    switch in direction of monotonicity, to produce a more "visually
!    pleasing" curve than that given by PCHIM.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!  Parameters:
!
!     SWITCH -- (input) indicates the amount of control desired over
!           local excursions from data.
!
!     N -- (input) number of data points.  (assumes N.GT.2 .)
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real array of derivative values at the data points,
!           as determined by PCHCI.
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually
!           pleasing" curve.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.  should be zero.
!           If negative, trouble in PCHSW.  (should never happen.)
!
  implicit none

  integer ( kind = 4 )  N, INCFD, IERR
  real ( kind = 4 ) SWITCH, H(*), SLOPE(*), D(INCFD,*)
  integer ( kind = 4 )  I, INDX, K, NLESS1
  real ( kind = 4 ) DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE, SLMAX
  real WTAVE(2)

  SAVE ONE, FUDGE
  real ( kind = 4 ) PCHST
!
!  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!
  REAL  PCHSD, S1, S2, H1, H2
  PCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!
!  INITIALIZE.
!
  DATA ONE /1./
  DATA  FUDGE /4./

  IERR = 0
  NLESS1 = N - 1
!
!  LOOP OVER SEGMENTS.
!
  DO I = 2, NLESS1

     IF ( PCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900

  100    CONTINUE
!
!  SLOPE SWITCHES MONOTONICITY AT I-TH POINT.
!
!  DO NOT CHANGE D IF 'UP-DOWN-UP'.
!
        IF (I .GT. 2)  THEN
           IF ( PCHST(SLOPE(I-2),SLOPE(I)) .GT. 0.0E+00)  GO TO 900
        end if
        IF (I .LT. NLESS1)  THEN
           IF ( PCHST(SLOPE(I+1),SLOPE(I-1)) .GT. 0.0E+00)  GO TO 900
        end if
!
!  COMPUTE PROVISIONAL VALUE FOR D(1,I).
!
        DEXT = PCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!
!  DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!
        IF ( PCHST(DEXT, SLOPE(I-1)) )  200, 900, 250

  200       CONTINUE
!
!  DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS.  EXTREMUM IS IN (X(I-1),X(I)).
!
           K = I-1
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
!
           WTAVE(2) = DEXT
           IF (K .GT. 1) then
              WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
           end if
           GO TO 400

  250       CONTINUE
!
!  DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS.  EXTREMUM IS IN (X(I),X(I+1)).
!
           K = I
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
!
           WTAVE(1) = DEXT
           IF (K .LT. NLESS1) then
              WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
           end if
           GO TO 400

  300    CONTINUE
!
!  AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO.  
!  CHECK FOR FLAT-TOPPED PEAK.
!
        IF (I == NLESS1)  GO TO 900
        IF ( PCHST(SLOPE(I-1), SLOPE(I+1)) .GE. 0.0E+00)  GO TO 900
!
!  WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
!
        K = I
!
!  SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
!
        WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
        WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))

  400    CONTINUE
!
!  AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!  ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!  WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K), IF K.GT.1
!  WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1), IF K.LT.N-1
!
     SLMAX = ABS(SLOPE(K))
     IF (K .GT. 1)    SLMAX = MAX( SLMAX, ABS(SLOPE(K-1)) )
     IF (K.LT.NLESS1) SLMAX = MAX( SLMAX, ABS(SLOPE(K+1)) )

     IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
     DEL(2) = SLOPE(K) / SLMAX
     IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX

     IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
!
!  NORMAL CASE: EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
!
        FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
        D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
        FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
        D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
     ELSE
!
!  SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!  K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
!
        FACT = FUDGE* ABS(DEL(2))
        D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
!
!  NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!  I-K+1 = 2 IF K=I-1(=1).
!
     end if
!
!  ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!
     IF (SWITCH .LE. 0.0E+00)  GO TO 900

     DFLOC = H(K)*ABS(SLOPE(K))
     IF (K .GT. 1)    DFLOC = MAX( DFLOC, H(K-1)*ABS(SLOPE(K-1)) )
     IF (K.LT.NLESS1) DFLOC = MAX( DFLOC, H(K+1)*ABS(SLOPE(K+1)) )
     DFMX = SWITCH*DFLOC
     INDX = I-K+1
!
!  INDX = 1 IF K=I, 2 IF K=I-1.
!
     CALL PCHSW (DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)

     IF (IERR .NE. 0)  RETURN

  900    CONTINUE

  end do

  RETURN
END
FUNCTION PCHDF ( K, X, S, IERR )

!*****************************************************************************80
!
!! PCHDF computes divided differences for PCHCE and PCHSP.
!
!  Discussion:
!
!    Uses a divided difference formulation to compute a K-point approx-
!    imation to the derivative at X(K) based on the data in X and S.
!
!    Called by PCHCE and PCHSP to compute 3- and 4-point boundary
!    derivative approximations.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K.LT.2 .
!        PCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
  implicit none

  integer ( kind = 4 )  K, IERR
  real ( kind = 4 ) pchdf
  real ( kind = 4 ) X(K), S(K)
  integer ( kind = 4 )  I, J
  real ( kind = 4 ) VALUE
!
!  CHECK FOR LEGAL VALUE OF K.
!
  IF (K .LT. 3)  GO TO 5001
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
  DO J = 2, K-1
     DO I = 1, K-J
        S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
     end do
  end do
!
!  EVALUATE DERIVATIVE AT X(K).
!
  VALUE = S(1)
  DO I = 2, K-1
     VALUE = S(I) + VALUE*(X(K)-X(I))
  end do

  IERR = 0
  PCHDF = VALUE

  RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!
!  K.LT.3 RETURN.
!
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHDF', 'K LESS THAN THREE', IERR, 1)
  PCHDF = 0.0E+00

  RETURN
END
SUBROUTINE PCHEV ( N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR )

!*****************************************************************************80
!
!! PCHEV: value and derivative of spline or cubic Hermite at many points.
!
!  Discussion:
!
!    Evaluates the function and first derivative of the cubic Hermite
!    or spline function defined by  N, X, F, D, at the array of points
!    XVAL.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence: CALL PCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
!
!     integer ( kind = 4 )  N, NVAL, IERR
!     real  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable 
!           values.  The elements of X must be strictly increasing:
!             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
!
!     F -- (input) real array of function values.  F(I) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  
!          D(I) is the value corresponding to X(I).
!
!  NVAL -- (input) number of points at which the functions are to be
!           evaluated. ( Error return if NVAL.LT.1 )
!
!  XVAL -- (input) real array of points at which the 
!          functions are to be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XVAL are increasing relative to X;
!              that is,   XVAL(J) .GE. X(I)
!              implies    XVAL(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XVAL are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  FVAL -- (output) real array of values of the cubic 
!          Hermite function defined by  N, X, F, D  at the points  XVAL.
!
!  DVAL -- (output) real array of values of the 
!          first derivative of the same function at the points  XVAL.
!
!  IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NVAL.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
  implicit none

  integer ( kind = 4 )  N, NVAL, IERR
  real ( kind = 4 ) X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
  integer ( kind = 4 ) INCFD
  LOGICAL SKIP

  incfd = 1
  skip = .true.

  CALL PCHFD ( N, X, F, D, INCFD, SKIP, NVAL, XVAL, FVAL, DVAL, IERR )

  RETURN
END
SUBROUTINE PCHEZ ( N, X, F, D, SPLINE, WK, LWK, IERR )

!*****************************************************************************80
!
!! PCHEZ sets up a spline or cubic Hermite interpolant.
!
!  Discussion:
!
!    Sets derivatives for spline (two continuous derivatives) or
!    Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data.  Various boundary
!    conditions are set to default values by PCHEZ. Many other choices
!    are available in the subroutines PCHIC, PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence:   CALL  PCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
!
!     integer ( kind = 4 )  N, IERR,  LWK
!     real  X(N), F(N), D(N), WK(*)
!     LOGICAL SPLINE
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(I) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!
!     SPLINE -- (input) logical variable to specify if the interpolant
!           is to be a spline with two continuous derivaties
!           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
!           continuous derivative (set SPLINE=.FALSE.).
!        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
!           default "not-a-knot" boundary condition, with a continuous
!           third derivative at X(2) and X(N-1).
!              If SPLINE=.FALSE. the interpolating Hermite cubic will be
!           monotone if the input data is monotone. Boundary conditions
!           computed from the derivative of a local quadratic unless thi
!           alters monotonicity.
!
!     WK -- (scratch) real work array, which must be declared by the cal
!           program to be at least 2*N if SPLINE is .TRUE. and not used
!           otherwise.
!
!     LWK -- (input) length of work array WK. (Error return if
!           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means tha
!                 IERR switches in the direction of monotonicity were de
!                 When SPLINE=.FALSE.,  PCHEZ guarantees that if the inp
!                 data is monotone, the interpolant will be too. This wa
!                 is to alert you to the fact that the input data was no
!                 monotone.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, LWK, IERR
  real ( kind = 4 ) X(N), F(N), D(N), WK(LWK)
  LOGICAL SPLINE
  integer ( kind = 4 ) IC(2), INCFD
  real ( kind = 4 ) VC(2)

  ic(1) = 0
  ic(2) = 0
  incfd = 1

  IF ( SPLINE ) THEN
    CALL PCHSP ( IC, VC, N, X, F, D, INCFD, WK, LWK, IERR )
  ELSE
    CALL PCHIM ( N, X, F, D, INCFD, IERR )
  end if

  RETURN
END
SUBROUTINE PCHFD ( N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR )

!*****************************************************************************80
!
!! PCHFD evaluates a piecewise cubic Hermite and derivative at many points.
!
!  Discussion:
!
!    Evaluate a piecewise cubic Hermite function and its first
!    derivative at an array of points.  May be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.
!
!    Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!    gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!    If only function values are required, use PCHFE, instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, NE, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
!        LOGICAL  SKIP
!
!        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, NE, IERR
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*), DE(*)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ( NE.LT.1 )  GO TO 5004
  IERR = 0
!
!  LOOP OVER INTERVALS.
!  INTERVAL INDEX IS  IL = IR-1  .
!  INTERVAL IS X(IL).LE.X.LT.X(IR) .
!
  JFIRST = 1
  IR = 2
   10 CONTINUE
!
!  SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
     IF (JFIRST .GT. NE) then
       return
     end if
!
!  LOCATE ALL POINTS IN INTERVAL.
!
     DO J = JFIRST, NE
        IF (XE(J) .GE. X(IR))  GO TO 30
     end do

     J = NE + 1
     GO TO 40
!
!  HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     IF (IR == N)  J = NE + 1

   40    CONTINUE
     NJ = J - JFIRST
!
!  SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
     IF (NJ == 0)  GO TO 50
!
!  EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
    CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
              NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
     IF (IERC .LT. 0)  GO TO 5005

     IF (NEXT(2) == 0)  GO TO 42
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
        IF (IR .LT. N)  GO TO 41
!           IF (IR == N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(2)
           GO TO 42
   41       CONTINUE
!
!  ELSE WE SHOULD NEVER HAVE GOTTEN HERE.
!
           GO TO 5005
!           end if
!        end if
   42    CONTINUE

     IF (NEXT(1) == 0)  GO TO 49
!
!  IF (NEXT(1) .GT. 0)  THEN
!  IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!  LEFT OF X(IR-1).
!
        IF (IR .GT. 2)  GO TO 43
!
! IF (IR == 2) THEN THESE ARE ACTUALLY EXTRAPOLATION POINTS.
!
           IERR = IERR + NEXT(1)
           GO TO 49
   43       CONTINUE
!
!  ELSE XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!  EVALUATION INTERVAL.
!
!  FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
!
           DO I = JFIRST, J-1
              IF (XE(I) .LT. X(IR-1))  GO TO 45
           end do
!
!  CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR IN CHFDV.
!
           GO TO 5005

   45          CONTINUE
!
!  RESET J.  (THIS WILL BE THE NEW JFIRST.)
!
           J = I
!
!  NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
!
           DO I = 1, IR-1
              IF (XE(J) .LT. X(I)) GO TO 47
           end do
!
!  CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!
!  AT THIS POINT, EITHER  XE(J) .LT. X(1)
!  OR      X(I-1) .LE. XE(J) .LT. X(I) .
!  RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE CYCLING.
!
           IR = MAX(1, I-1)
!           end if
!        end if
   49    CONTINUE

     JFIRST = J

   50 CONTINUE
  IR = IR + 1
  IF (IR .LE. N)  GO TO 10

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHFD', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHFD', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHFD', 'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'PCHFD', &
    'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN
!
 5005 CONTINUE
!
!     ERROR RETURN FROM CHFDV.  THIS CASE SHOULD NEVER OCCUR.
!
  IERR = -5
  CALL XERMSG ('SLATEC', 'PCHFD', 'ERROR RETURN FROM CHFDV -- FATAL', IERR, 2)

  RETURN
END
SUBROUTINE PCHFE ( N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR )

!*****************************************************************************80
!
!! PCHFE evaluates a piecewise cubic Hermite function at many points.
!
!  Discussion:
!
!    Evaluates the cubic Hermite function defined by  N, X, F, D  at
!    the points XE(J), J=1(1)NE.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, NE, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
!        LOGICAL  SKIP
!
!        CALL  PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the function is to be
!           evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!             (The FE-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, NE, IERR
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*)
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ( NE.LT.1 )  GO TO 5004
  IERR = 0
!
!  LOOP OVER INTERVALS.
!  INTERVAL INDEX IS  IL = IR-1  .
!  INTERVAL IS X(IL).LE.X.LT.X(IR) .
!
  JFIRST = 1
  IR = 2
   10 CONTINUE
!
!  SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
     IF (JFIRST .GT. NE) then
       return
     end if
!
!  LOCATE ALL POINTS IN INTERVAL.
!
     DO J = JFIRST, NE
        IF (XE(J) .GE. X(IR))  GO TO 30
     end do

     J = NE + 1
     GO TO 40
!
!  HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     IF (IR == N)  J = NE + 1

   40    CONTINUE
     NJ = J - JFIRST
!
!  SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
     IF (NJ == 0)  GO TO 50
!
!  EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
    CALL CHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
               NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)
     IF (IERC .LT. 0)  GO TO 5005

     IF (NEXT(2) == 0)  GO TO 42
!
!  IF (NEXT(2) .GT. 0)  THEN IN THE CURRENT SET OF XE-POINTS, THERE ARE 
!  NEXT(2) TO THE RIGHT OF X(IR).
!
        IF (IR .LT. N)  GO TO 41
!
! IF (IR == N)  THEN THESE ARE ACTUALLY EXTRAPOLATION POINTS.
!
           IERR = IERR + NEXT(2)
           GO TO 42
   41       CONTINUE
!
!  ELSE WE SHOULD NEVER HAVE GOTTEN HERE.
!
           GO TO 5005
!           end if
!        end if
   42    CONTINUE

     IF (NEXT(1) == 0)  GO TO 49
!
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
        IF (IR .GT. 2)  GO TO 43
!           IF (IR == 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(1)
           GO TO 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!  FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
!
           DO I = JFIRST, J-1
              IF (XE(I) .LT. X(IR-1))  GO TO 45
           end do
!
!  CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR IN CHFEV.
!
           GO TO 5005

   45          CONTINUE
!
!  RESET J.  (THIS WILL BE THE NEW JFIRST.)
!
           J = I
!
!  NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
!
           DO I = 1, IR-1
              IF (XE(J) .LT. X(I)) GO TO 47
           end do
!
!  CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!
!  AT THIS POINT, EITHER  XE(J) .LT. X(1)
!  OR      X(I-1) .LE. XE(J) .LT. X(I) .
!  RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!  CYCLING.
!
           IR = MAX(1, I-1)
!           end if
!        end if
   49    CONTINUE

     JFIRST = J

   50 CONTINUE
  IR = IR + 1
  IF (IR .LE. N)  GO TO 10

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHFE', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHFE', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHFE', 'X-ARRAY NOT STRICTLY INCREASING' &
    , IERR, 1)
  RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'PCHFE', &
    'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM CHFEV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -5
  CALL XERMSG ('SLATEC', 'PCHFE', &
    'ERROR RETURN FROM CHFEV -- FATAL', IERR, 2)
  RETURN
END
FUNCTION PCHIA ( N, X, F, D, INCFD, SKIP, A, B, IERR )

!*****************************************************************************80
!
!! PCHIA evaluates the definite integral of a piecewise cubic Hermite function.
!
!  Discussion:
!
!    Evaluates the definite integral of the cubic Hermite function
!    defined by  N, X, F, D  over the interval [A, B].
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
!        REAL  VALUE, PCHIA
!        LOGICAL  SKIP
!
!        VALUE = PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on return with IERR.GE.0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -4  in case of an error return from PCHID (which
!                         should never occur).
!
  implicit none

  integer ( kind = 4 )  N, INCFD, IERR
  real ( kind = 4 ) pchia
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*), A, B
  LOGICAL  SKIP
  integer ( kind = 4 )  I, IA, IB, IERD, IL, IR
  real ( kind = 4 ) VALUE, XA, XB
  real ( kind = 4 ) CHFIE, PCHID

  VALUE = 0.0E+00
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0
  IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
  IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
!
!  COMPUTE INTEGRAL VALUE.
!
  IF (A .NE. B)  THEN
     XA = MIN (A, B)
     XB = MAX (A, B)
     IF (XB .LE. X(2))  THEN
!
!  INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!
        VALUE = CHFIE (X(1),X(2), F(1,1),F(1,2), &
                                 D(1,1),D(1,2), A, B)

     ELSE IF (XA .GE. X(N-1))  THEN
!
!  INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!
        VALUE = CHFIE(X(N-1),X(N), F(1,N-1),F(1,N), &
                                  D(1,N-1),D(1,N), A, B)

     ELSE
!
!  NORMAL CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
!  LOCATE IA AND IB SUCH THAT
!  X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
!
        IA = 1
        DO I = 1, N-1
           IF (XA .GT. X(I))  IA = I + 1
        end do
!
!  IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
!  IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
!
        IB = N
        DO I = N, IA, -1
           IF (XB .LT. X(I))  IB = I - 1
        end do
!
!  IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
!  IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
!
!  COMPUTE THE INTEGRAL.
!
        IF (IB .LT. IA)  THEN
!
!  THIS MEANS IB = IA-1 AND (A,B) IS A SUBSET OF (X(IB),X(IA)).
!
           VALUE = CHFIE (X(IB),X(IA), F(1,IB),F(1,IA), &
                                      D(1,IB),D(1,IA), A, B)
        ELSE
!
!  FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
!  (Case (IB == IA) is taken care of by initialization
!  of VALUE to ZERO.)
!
           IF (IB .GT. IA)  THEN
              VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
              IF (IERD .LT. 0)  GO TO 5004
           end if
!
!  THEN ADD ON INTEGRAL OVER (XA,X(IA)).
!
           IF (XA .LT. X(IA))  THEN
              IL = MAX(1, IA-1)
              IR = IL + 1
              VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                       D(1,IL),D(1,IR), XA, X(IA))
           end if
!
!  THEN ADD ON INTEGRAL OVER (X(IB),XB).
!
           IF (XB .GT. X(IB))  THEN
              IR = MIN (IB+1, N)
              IL = IR - 1
              VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                       D(1,IL),D(1,IR), X(IB), XB)
           end if
!
!  FINALLY, ADJUST SIGN IF NECESSARY.
!
           IF (A .GT. B)  VALUE = -VALUE
        end if
     end if
  end if

  PCHIA = VALUE

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHIA', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHIA', 'INCREMENT LESS THAN ONE', IERR, 1)
  return

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHIA', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return
!
 5004 CONTINUE
!     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
  IERR = -4
  CALL XERMSG ('SLATEC', 'PCHIA', 'TROUBLE IN PCHID', IERR, 1)

  return
END
SUBROUTINE PCHIC ( IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR )

!*****************************************************************************80
!
!! PCHIC sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Sets derivatives needed to determine a piecewise monotone piece-
!    wise cubic interpolant to the data given in X and F satisfying the
!    boundary conditions specified by IC and VC.
!
!    The treatment of points where monotonicity switches direction is
!    controlled by argument SWITCH.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, NWK, IERR
!        REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  for the default boundary condition (the same as
!                     used by PCHIM).
!           If IBEG.NE.0, then its sign indicates whether the boundary
!                     derivative is to be adjusted, if necessary, to be
!                     compatible with monotonicity:
!              IBEG.GT.0  if no adjustment is to be performed.
!              IBEG.LT.0  if the derivative is to be adjusted for
!                     monotonicity.
!
!           Allowable values for the magnitude of IBEG are:
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!           IBEG = 5  to set D(1) so that the second derivative is con-
!              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!              This option is somewhat analogous to the "not a knot"
!              boundary condition provided by PCHSP.
!
!          NOTES (IBEG):
!           1. An error return is taken if ABS(IBEG).GT.5 .
!           2. Only in case  IBEG.LE.0  is it guaranteed that the
!              interpolant will be monotonic in the first interval.
!              If the returned value of D(1) lies between zero and
!              3*SLOPE(1), the interpolant will be monotonic.  This
!              is **NOT** checked if IBEG.GT.0 .
!           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
!              tonicity, a warning error is returned.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES (IEND):
!           1. An error return is taken if ABS(IEND).GT.5 .
!           2. Only in case  IEND.LE.0  is it guaranteed that the
!              interpolant will be monotonic in the last interval.
!              If the returned value of D(1+(N-1)*INCFD) lies between
!              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!              This is **NOT** checked if IEND.GT.0 .
!           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
!              achieve monotonicity, a warning error is returned.
!
!     VC -- (input) real array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     SWITCH -- (input) indicates desired treatment of points where
!           direction of monotonicity switches:
!           Set SWITCH to zero if interpolant is required to be mono-
!           tonic in each interval, regardless of monotonicity of data.
!             NOTES:
!              1. This will cause D to be set to zero at all switch
!                 points, thus forcing extrema there.
!              2. The result of using this option with the default boun-
!                 dary conditions will be identical to using PCHIM, but
!                 will generally cost more compute time.
!                 This option is provided only to facilitate comparison
!                 of different switch and/or boundary conditions.
!           Set SWITCH nonzero to use a formula based on the 3-point
!              difference formula in the vicinity of switch points.
!           If SWITCH is positive, the interpolant on each interval
!              containing an extremum is controlled to not deviate from
!              the data by more than SWITCH*DFLOC, where DFLOC is the
!              maximum of the change of F on this interval and its two
!              immediate neighbors.
!           If SWITCH is negative, no such control is to be imposed.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!           These values will determine a monotone cubic Hermite func-
!           tion on each subinterval on which the data are monotonic,
!           except possibly adjacent to switches in monotonicity.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real array of working storage.  The user may wish
!           to know that the returned values are:
!              WK(I)     = H(I)     = X(I+1) - X(I) ;
!              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!           for  I = 1(1)N-1.
!
!     NWK -- (input) length of work array.
!           (Error return if  NWK.LT.2*(N-1) .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if ABS(IBEG).GT.5 .
!              IERR = -5  if ABS(IEND).GT.5 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK.LT.2*(N-1) .
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer incfd

  integer ( kind = 4 )  IC(2), N, NWK, IERR
  real ( kind = 4 ) VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*), WK(NWK)
  integer ( kind = 4 )  I, IBEG, IEND, NLESS1
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO I = 2, N
     IF ( X(I).LE.X(I-1) )  GO TO 5003
  end do

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  IF (ABS(IBEG) .GT. 5)  IERR = IERR - 1
  IF (ABS(IEND) .GT. 5)  IERR = IERR - 2
  IF (IERR .LT. 0)  GO TO 5004
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  NLESS1 = N - 1
  IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
!
!  SET UP H AND SLOPE ARRAYS.
!
  DO I = 1, NLESS1
     WK(I) = X(I+1) - X(I)
     WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
  end do
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF ( NLESS1 .le. 1)  then
    D(1,1) = WK(2)
    D(1,N) = WK(2)
    GO TO 3000
  end if
!
!  NORMAL CASE  (N .GE. 3) .
!
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!
  CALL PCHCI (N, WK(1), WK(N), D, INCFD)
!
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!
  IF (SWITCH == 0.0E+00)  GO TO 3000

  CALL PCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR)

  IF (IERR .NE. 0)  GO TO 5008
!
!  SET END CONDITIONS.
!
 3000 CONTINUE
  IF ( (IBEG == 0) .AND. (IEND == 0) ) then
    return
  end if

  CALL PCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR)

  IF (IERR .LT. 0)  GO TO 5009

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHIC', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!
!     INCFD.LT.1 RETURN.
!
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHIC', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHIC', 'X-ARRAY NOT STRICTLY INCREASING' &
    , IERR, 1)
  RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  CALL XERMSG ('SLATEC', 'PCHIC', 'IC OUT OF RANGE', IERR, 1)
  RETURN
!
 5007 CONTINUE
!     NWK .LT. 2*(N-1)  RETURN.
  IERR = -7
  CALL XERMSG ('SLATEC', 'PCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
  RETURN
!
 5008 CONTINUE
!     ERROR RETURN FROM PCHCS.
  IERR = -8
  CALL XERMSG ('SLATEC', 'PCHIC', 'ERROR RETURN FROM PCHCS', IERR, 1)
  RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM PCHCE.
!  THIS CASE SHOULD NEVER OCCUR ***
  IERR = -9
  CALL XERMSG ('SLATEC', 'PCHIC', 'ERROR RETURN FROM PCHCE', IERR, 1)
  RETURN
END
FUNCTION PCHID ( N, X, F, D, INCFD, SKIP, IA, IB, IERR )

!*****************************************************************************80
!
!! PCHID: integral of a piecewise cubic Hermite function between data points.
!
!  Discussion:
!
!    Evaluates the definite integral of the cubic Hermite function
!    defined by  N, X, F, D  over the interval [X(IA), X(IB)].
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F- and D-arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IA, IB, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
!
!  Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
!
!     IA,IB -- (input) indices in X-array for the limits of integration.
!           both must be in the range [1,N].  (Error return if not.)
!           No restrictions on their relative values.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IA or IB is out of range.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer incfd

  integer ( kind = 4 )  N, IA, IB, IERR
  real ( kind = 4 ) pchid
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*)
  LOGICAL  SKIP

  integer ( kind = 4 )  I, IUP, LOW
  REAL  H, HALF, SIX, SUM, VALUE
  SAVE HALF, SIX

  DATA HALF /0.5/,  SIX /6./

  VALUE = 0.0E+00
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then
    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do
  end if

  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
  IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
  IERR = 0
!
!  COMPUTE INTEGRAL VALUE.
!
  IF (IA .NE. IB)  THEN
     LOW = MIN(IA, IB)
     IUP = MAX(IA, IB) - 1
     SUM = 0.0E+00
     DO I = LOW, IUP
        H = X(I+1) - X(I)
        SUM = SUM + H*( (F(1,I) + F(1,I+1)) + &
                       (D(1,I) - D(1,I+1))*(H/SIX) )
     end do
     VALUE = HALF * SUM
     IF (IA .GT. IB)  VALUE = -VALUE

  end if

  PCHID = VALUE

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHID', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHID', 'INCREMENT LESS THAN ONE', IERR, 1)
  return

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHID', &
    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return
!
 5004 CONTINUE
!     IA OR IB OUT OF RANGE RETURN.
  IERR = -4
  CALL XERMSG ('SLATEC', 'PCHID', 'IA OR IB OUT OF RANGE', IERR, 1)

  return
END
SUBROUTINE PCHIM ( N, X, F, D, INCFD, IERR )

!*****************************************************************************80
!
!! PCHIM sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.
!
!    Sets derivatives needed to determine a monotone piecewise cubic
!    Hermite interpolant to the data given in X and F.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  (See PCHIC if user control of boundary con-
!    ditions is desired.)
!
!    If the data are only piecewise monotonic, the interpolant will
!    have an extremum at each point where monotonicity switches direc-
!    tion.  (See PCHIC if user control is desired in such cases.)
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  PCHIM (N, X, F, D, INCFD, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!           PCHIM is designed for monotonic data, but it will work for
!           any F-array.  It will force extrema at points where mono-
!           tonicity switches direction.  If some other treatment of
!           switch points is desired, PCHIC should be used instead.
!
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, INCFD, IERR
  real ( kind = 4 ) X(*), F(INCFD,*), D(INCFD,*)
  integer ( kind = 4 )  I, NLESS1
  real ( kind = 4 ) DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE
  real ( kind = 4 ) H1, H2, HSUM, HSUMT3, W1, W2
  real ( kind = 4 ) PCHST
!
!  VALIDITY-CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO I = 2, N
     IF ( X(I).LE.X(I-1) )  GO TO 5003
  end do
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0
  NLESS1 = N - 1
  H1 = X(2) - X(1)
  DEL1 = (F(1,2) - F(1,1))/H1
  DSAVE = DEL1
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF (NLESS1 .le. 1) then
    D(1,1) = DEL1
    D(1,N) = DEL1
    return
  end if
!
!  NORMAL CASE  (N .GE. 3).
!
  H2 = X(3) - X(2)
  DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  HSUM = H1 + H2
  W1 = (H1 + HSUM)/HSUM
  W2 = -H1/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2
  IF ( PCHST(D(1,1),DEL1) .LE. 0.0E+00)  THEN
     D(1,1) = 0.0E+00
  ELSE IF ( PCHST(DEL1,DEL2) .LT. 0.0E+00)  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0E+00 *DEL1
     IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO I = 2, NLESS1

     IF ( 2 < I ) then
       H1 = H2
       H2 = X(I+1) - X(I)
       HSUM = H1 + H2
       DEL1 = DEL2
       DEL2 = (F(1,I+1) - F(1,I))/H2
     end if
!
!  SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = 0.0E+00
     IF ( PCHST(DEL1,DEL2) )  42, 41, 45
!
!  COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
     IF (DEL2 == 0.0E+00)  GO TO 50
     IF ( PCHST(DSAVE,DEL2) .LT. 0.0E+00)  IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50

   42    CONTINUE
     IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50
!
!  USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H1)/HSUMT3
     W2 = (HSUM + H2)/HSUMT3
     DMAX = MAX( ABS(DEL1), ABS(DEL2) )
     DMIN = MIN( ABS(DEL1), ABS(DEL2) )
     DRAT1 = DEL1/DMAX
     DRAT2 = DEL2/DMAX
     D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)

   50    CONTINUE

  end do
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  W1 = -H2/HSUM
  W2 = (H2 + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2
  IF ( PCHST(D(1,N),DEL2) .LE. 0.0E+00)  THEN
     D(1,N) = 0.0E+00
  ELSE IF ( PCHST(DEL1,DEL2) .LT. 0.0E+00)  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0E+00 *DEL2
     IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
  end if

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHIM', &
    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHIM', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHIM', 'X-ARRAY NOT STRICTLY INCREASING' &
    , IERR, 1)

  RETURN
END
SUBROUTINE PCHKT ( N, X, KNOTYP, T )

!*****************************************************************************80
!
!! PCHKT computes the B-spline knot sequence for PCHBS.
!
!  Discussion:
!
!    Set a knot sequence for the B-spline representation of a PCH
!    function with breakpoints X.  All knots will be at least double.
!    Endknots are set as:
!    (1) quadruple knots at endpoints if KNOTYP=0;
!    (2) extrapolate the length of end interval if KNOTYP=1;
!    (3) periodic if KNOTYP=2.
!
!    Restrictions/assumptions:
!    1. N.GE.2 .  (not checked)
!    2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
!    3. 0.LE.KNOTYP.LE.2 .  (Acts like KNOTYP=0 for any other value.)
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!  Input arguments:  N, X, KNOTYP.
!  Output arguments:  T.
!
  implicit none

  integer ( kind = 4 )  N, KNOTYP
  real ( kind = 4 ) X(*), T(*)
  integer ( kind = 4 )  J, K, NDIM
  real ( kind = 4 ) HBEG, HEND
!
!  Initialize.
!
  NDIM = 2*N
!
!  Set interior knots.
!
  J = 1
  DO K = 1, N
     J = J + 2
     T(J) = X(K)
     T(J+1) = T(J)
  end do
!
!  Assertion:
!  At this point T(3),...,T(NDIM+2) have been set and J=NDIM+1.
!
!  Set end knots according to KNOTYP.
!
  HBEG = X(2) - X(1)
  HEND = X(N) - X(N-1)
  IF (KNOTYP == 1 )  THEN
!
!  Extrapolate.
!
     T(2) = X(1) - HBEG
     T(NDIM+3) = X(N) + HEND
  ELSE IF ( KNOTYP == 2 )  THEN
!
!  Periodic.
!
     T(2) = X(1) - HEND
     T(NDIM+3) = X(N) + HBEG
  ELSE
!
!  Quadruple end knots.
!
     T(2) = X(1)
     T(NDIM+3) = X(N)
  end if
  T(1) = T(2)
  T(NDIM+4) = T(NDIM+3)

  RETURN
END
FUNCTION PCHQA ( N, X, F, D, A, B, IERR )

!*****************************************************************************80
!
!! PCHQA: definite integral of spline or piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Evaluates the definite integral of a piecewise cubic Hermite
!    or spline function over an arbitrary interval, easy to use.
!
!    Evaluates the definite integral of the cubic Hermite or spline
!    function defined by  N, X, F, D  over the interval [A, B].
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    David Kahaner
!
!  Reference:
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Calling sequence:
!
!           VALUE = PCHQA (N, X, F, D, A, B, IERR)
!
!     integer ( kind = 4 )  N, IERR
!     real  X(N), F(N), D(N), A, B
!
!  Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable
!           values.  The elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.
!           F(I) is the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(I) is
!           the value corresponding to X(I).
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (Value has not been computed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer ( kind = 4 )  N, IERR
  real ( kind = 4 ) pchqa
  real ( kind = 4 ) X(N), F(N), D(N), A, B

  integer ( kind = 4 )  INCFD
  real ( kind = 4 ) PCHIA
  LOGICAL SKIP

  incfd = 1
  skip = .true.

  PCHQA = PCHIA ( N, X, F, D, INCFD, SKIP, A, B, IERR )

  RETURN
END
SUBROUTINE PCHQK1 ( LUN, KPRINT, IPASS )

!*****************************************************************************80
!
!! PCHQK1 tests the PCHIP evaluators CHFDV, CHFEV, PCHFD, PCHFE.
!
!  Discussion:
!
!    This routine carries out three tests of the PCH evaluators:
!    DEVCHK tests the single-cubic evaluators.
!    DEVPCK tests the full PCH evaluators.
!    DEVERK exercises the error returns in all evaluators.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!     SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!     IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  integer ( kind = 4 ) I1, I2, I3, I4, I5, I6, I7, I8, I9, IFAIL, NPTS
  real ( kind = 4 ) WORK (4000)
  LOGICAL  FAIL

      IF (KPRINT .GE. 2)  WRITE (LUN, 1000) KPRINT
!
!  TEST CHFDV AND CHFEV.
!
      IFAIL = 0
      NPTS = 1000
      I1 = 1  + NPTS
      I2 = I1 + NPTS
      I3 = I2 + NPTS
      CALL EVCHCK (LUN, KPRINT, NPTS, WORK(1), WORK(I1), WORK(I2), &
                                                WORK(I3), FAIL)
      IF (FAIL)  IFAIL = IFAIL + 1
!
!  TEST PCHFD AND PCHFE.
!
      I1 = 1  +  10
      I2 = I1 +  10
      I3 = I2 + 100
      I4 = I3 + 100
      I5 = I4 + 100
      I6 = I5 +  51
      I7 = I6 +  51
      I8 = I7 +  51
      I9 = I8 +  51
      CALL EVPCCK (LUN, KPRINT, WORK(1), WORK(I1), WORK(I2), WORK(I3), &
                   WORK(I4), WORK(I5), WORK(I6), WORK(I7), WORK(I8), &
                   WORK(I9), FAIL)
      IF (FAIL)  IFAIL = IFAIL + 2
!
!  TEST ERROR RETURNS.
!
      CALL EVERCK (LUN, KPRINT, FAIL)
      IF (FAIL)  IFAIL = IFAIL + 4
!
!  PRINT SUMMARY AND TERMINATE.
!  At this point, IFAIL has the following value:
!  IFAIL = 0  IF ALL TESTS PASSED.
!  IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
!  IFAIL=1  IF SINGLE CUBIC TEST FAILED. (SEE EVCHCK OUTPUT.)
!  IFAIL=2  IF PCHFD/PCHFE  TEST FAILED. (SEE EVPCCK OUTPUT.)
!  IFAIL=4  IF ERROR RETURN TEST FAILED. (SEE EVERCK OUTPUT.)
!
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN
 1000 FORMAT ('1'/' ------------  PCHIP QUICK CHECK OUTPUT', &
              ' ------------' //20X,'( KPRINT =',I2,' )')
 3001 FORMAT (/' *** TROUBLE ***',I5,' EVALUATION TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL EVALUATION TESTS', &
              ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME EVALUATION TESTS', &
              ' ************')
END
SUBROUTINE PCHQK2 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! PCHQK2 tests the PCHIP integrators PCHIA and PCHID.
!
!  Discussion:
!
!    This routine constructs data from a cubic, integrates it with PCHIA
!    and compares the results with the correct answer.
!    Since PCHIA calls PCHID, this tests both integrators.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!    LUN   :IN  is the unit number to which output is to be written.
!
!    KPRINT:IN  controls the amount of output, as specified in the
!        SLATEC Guidelines.
!
!    IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!        IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  integer ( kind = 4 ) I, IEREXP(17), IERR, IFAIL, N, NPAIRS
  real ( kind = 4 ) A(17), B(17), CALC, D(7), ERRMAX, ERROR, F(7), MACHEP
  real ( kind = 4 ) ONE, THRQTR, TOL, TRUE, TWO, X(7)
  LOGICAL  FAIL, SKIP
  real ( kind = 4 ) PCHIA, R1MACH
!
!  DEFINE TEST FUNCTIONS.
!
      REAL  AX, FCN, DERIV, ANTDER
      FCN(AX) = 3.0E+00 *AX*AX*(AX-TWO)
      DERIV(AX) = 3.0E+00 *AX*(TWO*(AX-TWO) + AX)
      ANTDER(AX) = AX**3 * (THRQTR*AX - TWO)
!
!  INITIALIZE.
!
      DATA  THRQTR /0.75E0/,  ONE /1.E0/,  TWO /2.E0/
      DATA  N /7/
      DATA  X /-4.E0, -2.E0, -0.9E0, 0.E0, 0.9E0, 2.E0, 4.E0/
      DATA  NPAIRS /17/
      DATA  A /-3.0E0, 3.0E0,-0.5E0,-0.5E0,-0.5E0,-4.0E0,-4.0E0, 3.0E0, &
        -5.0E0,-5.0E0,-6.0E0, 6.0E0,-1.5E0,-1.5E0,-3.0E0, 3.0E0, 0.5E0/
      DATA  B / 3.0E0,-3.0E0, 1.0E0, 2.0E0, 5.0E0,-0.5E0, 4.0E0, 5.0E0, &
        -3.0E0, 5.0E0,-5.0E0, 5.0E0,-0.5E0,-1.0E0,-2.5E0, 3.5E0, 0.5E0/
      DATA  IEREXP /0,0,0,0,2,0,0,2,1,3,3,3,0,0,0,0,0/
!
!  SET PASS/FAIL TOLERANCE.
!
      MACHEP = R1MACH(4)
      TOL = 100.E0*MACHEP
!
!  SET UP PCH FUNCTION DEFINITION.
!
      DO I = 1, N
         F(I) =   FCN(X(I))
         D(I) = DERIV(X(I))
      end do

      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
      IF (KPRINT .GE. 3)  WRITE (LUN, 1002)  (X(I), F(I), D(I), I=1,N)
!
!  LOOP OVER (A,B)-PAIRS.
!
      IF (KPRINT .GE. 3)  WRITE (LUN, 2000)

      IFAIL = 0

      SKIP = .FALSE.

      DO I = 1, NPAIRS

         CALC = PCHIA (N, X, F, D, 1, SKIP, A(I), B(I), IERR)

         IF (IERR .GE. 0)  THEN
            FAIL = IERR .NE. IEREXP(I)
            TRUE = ANTDER(B(I)) - ANTDER(A(I))
            ERROR = CALC - TRUE
            IF (KPRINT .GE. 3)  THEN
               IF (FAIL)  THEN
                 WRITE (LUN, 2001) A(I), B(I), IERR, TRUE, CALC, ERROR, &
                                                IEREXP(I)
               ELSE
                 WRITE (LUN, 2002) A(I), B(I), IERR, TRUE, CALC, ERROR
               ENDIF
            ENDIF

            ERROR = ABS(ERROR) / MAX(ONE, ABS(TRUE))
            IF (FAIL .OR. (ERROR.GT.TOL))  IFAIL = IFAIL + 1
            IF (I .EQ. 1)  THEN
               ERRMAX = ERROR
            ELSE
               ERRMAX = MAX(ERRMAX, ERROR)
            ENDIF
         ELSE
            IF (KPRINT .GE. 3)  WRITE (LUN, 2002)  A(I), B(I), IERR
            IFAIL = IFAIL + 1
         ENDIF

      end do
!
!  PRINT SUMMARY.
!
      IF (KPRINT .GE. 2)  THEN
         WRITE (LUN, 2003)  ERRMAX, TOL
         IF (IFAIL .NE. 0)  WRITE (LUN, 3001)  IFAIL
      ENDIF
!
!  TERMINATE.
!
      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCHIP INTEGRATORS')
 1001 FORMAT (//10X,'PCHQK2 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' //11X,'X',9X,'F',9X,'D' /(5X,3F10.3) )
 2000 FORMAT (// 5X,'TEST RESULTS:' &
              //'    A     B    ERR     TRUE',16X,'CALC',15X,'ERROR')
 2001 FORMAT (2F6.1,I5,1P,2E20.10,E15.5,'  (',I1,') *****' )
 2002 FORMAT (2F6.1,I5,1P,2E20.10,E15.5)
 2003 FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,E15.5, &
                             ',   TOLERANCE:',1P,E15.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL INTEGRATION TESTS', &
              ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME INTEGRATION TESTS', &
              ' ************')
END
SUBROUTINE PCHQK3 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! PCHQK3 tests the PCHIP interpolators PCHIC, PCHIM, PCHSP.
!
!  Discussion:
!
!    This routine interpolates a constructed data set with all three
!    PCHIP interpolators and compares the results with those obtained
!    on a Cray X/MP.  Two different values of the PCHIC parameter SWITCH
!    are used.
!
!    1. The Cray results are given only to nine significant figures,
!       so don't expect them to match to more.
!    2. The results will depend to some extent on the accuracy of
!       the EXP function.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!        SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!        IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  LOGICAL  COMP
  real ( kind = 4 ) R1MACH
  integer ( kind = 4 ) I, IC(2), IERR, IFAIL, N, NBAD, NBADZ, NWK
  PARAMETER  (N = 9,  NWK = 2*N)
  real ( kind = 4 ) D(N), DC(N), DC5, DC6, DM(N), DS(N), ERR, F(N), MONE, TOL
  real ( kind = 4 ) TOLD, TOLZ, VC(2), X(N), WK(NWK), ZERO
  PARAMETER  (ZERO = 0.0E0,  MONE = -1.0E0)
  CHARACTER*6  RESULT
!
!  Initialize.
!
      DATA  IC /0, 0/
      DATA  X /-2.2E0,-1.2E0,-1.0E0,-0.5E0,-0.01E0, 0.5E0, 1.0E0, &
                2.0E0, 2.2E0/
!
!  Results generated on Cray X/MP (9 sign. figs.)
!
      DATA  DM / 0.            , 3.80027352E-01, 7.17253009E-01, &
                 5.82014161E-01, 0.            ,-5.68208031E-01, &
                -5.13501618E-01,-7.77910977E-02,-2.45611117E-03/
      DATA  DC5,DC6 / 1.76950158E-02,-5.69579814E-01/
      DATA  DS /-5.16830792E-02, 5.71455855E-01, 7.40530225E-01, &
                 7.63864934E-01, 1.92614386E-02,-7.65324380E-01, &
                -7.28209035E-01,-7.98445427E-02,-2.85983446E-02/
      IFAIL = 0
!
!  Set tolerances.
!
      TOL  = 10*R1MACH(4)
      TOLD = MAX( 1.0E-7, 10*TOL )
      TOLZ = ZERO

      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
!
!  Set up data.
!
      DO I = 1, N
         F(I) = EXP(-X(I)**2)
      end do

      IF (KPRINT .GE. 3)  THEN
         WRITE (LUN, 1002)
         DO I = 1, 4
            WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
         end do
         WRITE (LUN, 1011)  X(5), F(5), DM(5), DC5, DS(5)
         WRITE (LUN, 1011)  X(6), F(6), DM(6), DC6, DS(6)
         DO I = 7, N
            WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
         end do
      ENDIF
!
!  Test PCHIM.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IM'

      CALL PCHIM (N, X, F, D, 1, IERR)
!
!  Expect IERR=1 (one monotonicity switch).
!
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 1
      IF ( .NOT.COMP (IERR, 1, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         NBADZ = 0
         DO I = 1, N
            RESULT = '  OK'
!
!  D-values should agree with stored values.
!  (Zero values should agree exactly.)
!
            IF ( DM(I).EQ.ZERO )  THEN
               ERR = ABS( D(I) )
               IF ( ERR.GT.TOLZ )  THEN
                  NBADZ = NBADZ + 1
                  RESULT = '**BADZ'
               ENDIF
            ELSE
               ERR = ABS( (D(I)-DM(I))/DM(I) )
               IF ( ERR.GT.TOLD )  THEN
                  NBAD = NBAD + 1
                  RESULT = '**BAD'
               ENDIF
            ENDIF
            IF (KPRINT.GE.3) then
               WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
            end if

         end do

         IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
            IFAIL = IFAIL + 1
            IF ((NBADZ.NE.0).AND.(KPRINT.GE.2)) then
               WRITE (LUN, 2004)  NBAD
            end if
            IF ((NBAD.NE.0).AND.(KPRINT.GE.2)) then
               WRITE (LUN, 2005)  NBAD, 'IM', TOLD
            end if
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IM'
         ENDIF
      ENDIF
!
!  Test PCHIC -- options set to reproduce PCHIM.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'

      CALL PCHIC (IC, VC, ZERO, N, X, F, DC, 1, WK, NWK, IERR)
!
!  Expect IERR=0 .
!
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE

         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         DO I = 1, N
            RESULT = '  OK'
!
!  D-values should agree exactly with those computed by PCHIM.
!  (To be generous, will only test to machine precision.)
!
            ERR = ABS( D(I)-DC(I) )
            IF ( ERR.GT.TOL )  THEN
               NBAD = NBAD + 1
               RESULT = '**BAD'
            ENDIF
            IF (KPRINT.GE.3) then
              WRITE (LUN, 2003)  I, X(I), DC(I), ERR, RESULT
            end if
         end do

         IF ( NBAD.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'IC', TOL
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
         ENDIF
      ENDIF
!
!  Test PCHIC -- default nonzero switch derivatives.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'

      CALL PCHIC (IC, VC, MONE, N, X, F, D, 1, WK, NWK, IERR)
!
!  Expect IERR=0 .
!
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         NBADZ = 0
         DO I = 1, N
            RESULT = '  OK'
!
!  D-values should agree exactly with those computed in
!  previous call, except at points 5 and 6.
!
            IF ( (I.LT.5).OR.(I.GT.6) )  THEN
               ERR = ABS( D(I)-DC(I) )
               IF ( ERR.GT.TOLZ )  THEN
                  NBADZ = NBADZ + 1
                  RESULT = '**BADA'
               ENDIF
            ELSE
               IF ( I.EQ.5 )  THEN
                  ERR = ABS( (D(I)-DC5)/DC5 )
               ELSE
                  ERR = ABS( (D(I)-DC6)/DC6 )
               ENDIF
               IF ( ERR.GT.TOLD )  THEN
                  NBAD = NBAD + 1
                  RESULT = '**BAD'
               ENDIF
            ENDIF
            IF (KPRINT.GE.3) then
              WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
            end if
         end do

         IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
            IFAIL = IFAIL + 1
            IF ((NBADZ.NE.0).AND.(KPRINT.GE.2)) then
              WRITE (LUN, 2007)  NBAD
            end if
            IF ((NBAD.NE.0).AND.(KPRINT.GE.2)) then
              WRITE (LUN, 2005)  NBAD, 'IC', TOLD
            end if
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
         ENDIF
      ENDIF
!
!  Test PCHSP.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'SP'

      CALL PCHSP (IC, VC, N, X, F, D, 1, WK, NWK, IERR)
!
!  Expect IERR=0 .
!
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         DO I = 1, N
            RESULT = '  OK'
!
!  D-values should agree with stored values.
!
            ERR = ABS( (D(I)-DS(I))/DS(I) )
            IF ( ERR.GT.TOLD )  THEN
               NBAD = NBAD + 1
               RESULT = '**BAD'
            ENDIF
            IF (KPRINT.GE.3) then
              WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
            end if
         end do

         IF ( NBAD.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'SP', TOLD
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'SP'
         ENDIF
      ENDIF
!
!  PRINT SUMMARY AND TERMINATE.
!
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCHIP INTERPOLATORS')
 1001 FORMAT (//10X,'PCHQK3 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' &
              /39X,'---------- EXPECTED D-VALUES ----------' &
              /12X,'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
 1010 FORMAT (5X,F10.2,1P,E15.5,4X,E15.5,15X,E15.5)
 1011 FORMAT (5X,F10.2,1P,E15.5,4X,3E15.5)
 2000 FORMAT (/5X,'PCH',A2,' TEST:')
 2001 FORMAT (15X,'EXPECT  IERR =',I5)
 2002 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
 2003 FORMAT (5X,I5,F10.2,1P,2E15.5,2X,A)
 2004 FORMAT (/'    **',I5,'  PCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
 2005 FORMAT (/'    **',I5,'  PCH',A2,' RESULTS FAILED TOLERANCE TEST.', &
                          '  TOL =',1P,E10.3)
 2006 FORMAT (/5X,'  ALL  PCH',A2,' RESULTS OK.')
 2007 FORMAT (/'    **',I5,'  PCHIC RESULTS FAILED TO AGREE WITH', &
                           ' PREVIOUS CALL.')
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL INTERPOLATION TESTS', &
             ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME INTERPOLATION TESTS', &
             ' ************')
END
SUBROUTINE PCHQK4 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! PCHQK4 tests the PCHIP monotonicity checker PCHCM.
!
!  Discussion:
!
!    This routine tests a constructed data set with three different
!    INCFD settings and compares with the expected results.  It then
!    runs a special test to check for bug in overall monotonicity found
!    in PCHMC.  Finally, it reverses the data and repeats all tests.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS

  integer ( kind = 4 ) MAXN, MAXN2, MAXN3, NB
  PARAMETER  (MAXN = 16,  MAXN2 = 8,  MAXN3 = 6,  NB = 7)
  integer ( kind = 4 ) I, IERR, IFAIL, INCFD, ISMEX1(MAXN), ISMEX2(MAXN2)
  integer ( kind = 4 ) ISMEX3(MAXN3), ISMEXB(NB), ISMON(MAXN), K, N, NS(3)
  real ( kind = 4 ) D(MAXN), DB(NB), F(MAXN), FB(NB), X(MAXN)
  LOGICAL  SKIP
!
!  DEFINE EXPECTED RESULTS.
!
      DATA  ISMEX1 / 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 2/
      DATA  ISMEX2 / 1, 2, 2, 1, 2, 2, 1, 2/
      DATA  ISMEX3 / 1, 1, 1, 1, 1, 1/
      DATA  ISMEXB / 1, 3, 1, -1, -3, -1, 2/
!
!  DEFINE TEST DATA.
!
      DATA  NS /16, 8, 6/

      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
!
!  Define X, F, D.
!
      DO I = 1, MAXN
         X(I) = I
         D(I) = 0.E0
      end do

      DO I = 2, MAXN, 3
         D(I) = 2.E0
      end do

      DO I = 1, 3
         F(I) = X(I)
         F(I+ 3) = F(I  ) + 1.E0
         F(I+ 6) = F(I+3) + 1.E0
         F(I+ 9) = F(I+6) + 1.E0
         F(I+12) = F(I+9) + 1.E0
      end do

      F(16) = 6.E0
!
!  Define FB, DB.
!
      FB(1) = 0.E0
      FB(2) = 2.E0
      FB(3) = 3.E0
      FB(4) = 5.E0
      DB(1) = 1.E0
      DB(2) = 3.E0
      DB(3) = 3.E0
      DB(4) = 0.E0
      DO I = 1, 3
         FB(NB-I+1) =  FB(I)
         DB(NB-I+1) = -DB(I)
      end do
!
!  INITIALIZE.
!
      IFAIL = 0

      IF (KPRINT .GE. 3)  THEN
         WRITE (LUN, 1002)
         DO I = 1, NB
            WRITE (LUN, 1010)  I, X(I), F(I), D(I), FB(I), DB(I)
         end do
         DO I = NB+1, MAXN
            WRITE (LUN, 1010)  I, X(I), F(I), D(I)
         end do
      ENDIF
!
!  TRANSFER POINT FOR SECOND SET OF TESTS.
!
   25 CONTINUE
!
!  Loop over a series of values of INCFD.
!
      DO INCFD = 1, 3

         N = NS(INCFD)
         SKIP = .FALSE.
         CALL PCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
         IF (KPRINT.GE.3) then
           WRITE (LUN, 2000)  INCFD, IERR, (ISMON(I), I=1,N)
         end if
         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN,2001)
         ELSE
            DO I = 1, N
               IF (INCFD.EQ.1)  THEN
                  IF ( ISMON(I).NE.ISMEX1(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                       WRITE (LUN, 2002)  (ISMEX1(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ELSE IF (INCFD.EQ.2) THEN
                  IF ( ISMON(I).NE.ISMEX2(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                        WRITE (LUN, 2002)  (ISMEX2(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ELSE
                  IF ( ISMON(I).NE.ISMEX3(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3) then
                       WRITE (LUN, 2002)  (ISMEX3(K),K=1,N)
                     end if
                     GO TO 30
                  ENDIF
               ENDIF
            end do
         ENDIF

   30 CONTINUE

      end do
!
!  Test for -1,3,1 bug.
!
      SKIP = .FALSE.

      CALL PCHCM (NB, X, FB, DB, 1, SKIP, ISMON, IERR)

      IF (KPRINT.GE.3) then
        WRITE (LUN, 2030)  IERR, (ISMON(I), I=1,NB)
      end if
      IF ( IERR.NE.0 )  THEN
         IFAIL = IFAIL + 1
         IF (KPRINT.GE.3)  WRITE (LUN,2001)
      ELSE
         DO I = 1, NB
            IF ( ISMON(I).NE.ISMEXB(I) )  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.3) then
                 WRITE (LUN, 2002)  (ISMEXB(K),K=1,NB)
               end if
               GO TO 35
            ENDIF
         end do
      ENDIF
   35 CONTINUE

      IF (F(1).LT.0.)  GO TO 90
!
!  Change sign and do again.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 2050)
      DO I = 1, MAXN
         F(I) = -F(I)
         D(I) = -D(I)
         IF ( ISMEX1(I).NE.2 )  ISMEX1(I) = -ISMEX1(I)
      end do
      DO I = 1, MAXN2
         IF ( ISMEX2(I).NE.2 )  ISMEX2(I) = -ISMEX2(I)
      end do
      DO I = 1, MAXN3
         IF ( ISMEX3(I).NE.2 )  ISMEX3(I) = -ISMEX3(I)
      end do
      DO I = 1, NB
         FB(I) = -FB(I)
         DB(I) = -DB(I)
         IF ( ISMEXB(I).NE.2 )  ISMEXB(I) = -ISMEXB(I)
      end do
      GO TO 25
!
!  PRINT SUMMARY AND TERMINATE.
!
   90 CONTINUE
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCHIP MONOTONICITY CHECKER')
 1001 FORMAT (//10X,'PCHQK4 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:' &
             // 9X,'I',4X,'X',5X,'F',5X,'D',5X,'FB',4X,'DB')
 1010 FORMAT (5X,I5,5F6.1)
 2000 FORMAT (/4X,'INCFD =',I2,':  IERR =',I3/15X,'ISMON =',16I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (' *** Failed -- expect:',16I3)
 2030 FORMAT (/4X,' Bug test:  IERR =',I3/15X,'ISMON =',7I3)
 2050 FORMAT (/4X,'Changing sign of data.....')
 3001 FORMAT (/' *** TROUBLE ***',I5,' MONOTONICITY TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL MONOTONICITY TESTS', &
             ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME MONOTONICITY TESTS', &
             ' ************')
END
SUBROUTINE PCHQK5 (LUN, KPRINT, IPASS)

!*****************************************************************************80
!
!! PCHQK5 tests the PCH to B-spline conversion routine PCHBS.
!
!  Discussion:
!
!    This routine tests a constructed data set with four different
!    KNOTYP settings.  It computes the function and derivatives of the
!    resulting B-representation via BVALU and compares with PCH data.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
  implicit none

  integer ( kind = 4 ) LUN, KPRINT, IPASS
  real ( kind = 4 ) BVALU, R1MACH
  EXTERNAL  BVALU, PCHBS, R1MACH
  integer ( kind = 4 ) I, IERR, IFAIL, INBV, J, KNOTYP, K, N, NDIM, NKNOTS
  PARAMETER  (N = 9)
  real ( kind = 4 ) BCOEF(2*N), D(N), DCALC, DERR, DERMAX, F(N), FCALC, FERR
  real ( kind = 4 ) FERMAX, T(2*N+4), TERR, TERMAX, TOL, TOLZ, TSAVE(2*N+4)
  real ( kind = 4 ) WORK(16*N), X(N), ZERO
  PARAMETER  (ZERO = 0.0E0)
  LOGICAL  FAIL
!
!  Define relative error function.
!
  real ( kind = 4 ) ANS, ERR, RELERR
  RELERR (ERR, ANS) = ABS(ERR) / MAX(1.0E-5,ABS(ANS))
!
!  Define test data.
!
      DATA  X /-2.2E0,   -1.2E0,   -1.0E0,   -0.5E0,   -0.01E0, &
                0.5E0,    1.0E0,    2.0E0,    2.2E0/
      DATA  F / 0.0079E0, 0.2369E0, 0.3679E0, 0.7788E0, 0.9999E0, &
                0.7788E0, 0.3679E0, 0.1083E0, 0.0079E0/
      DATA  D / 0.0000E0, 0.3800E0, 0.7173E0, 0.5820E0, 0.0177E0, &
               -0.5696E0,-0.5135E0,-0.0778E0,-0.0025E0/
!
!  Initialize.
!
      IFAIL = 0
      TOL = 100*R1MACH(4)
      TOLZ = ZERO

      IF (KPRINT.GE.3)  WRITE (LUN, 1000)
      IF (KPRINT.GE.2)  WRITE (LUN, 1001)
!
!  Loop over a series of values of KNOTYP.
!
      IF (KPRINT.GE.3)  WRITE (LUN, 1010)

      DO KNOTYP = 2, -1, -1

         CALL PCHBS (N, X, F, D, 1, KNOTYP, NKNOTS, T, BCOEF, NDIM, K, IERR)

         IF (KPRINT.GE.3) then
            WRITE (LUN, 2000) KNOTYP, NKNOTS, NDIM, K, IERR
         end if

         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN, 2001)
         ELSE
!
!  Compare evaluated results with inputs to PCHBS.
!
            INBV = 1
            FERMAX = ZERO
            DERMAX = ZERO
            IF (KPRINT.GE.3)  THEN
               WRITE (LUN, 2002)
               WRITE (LUN, 2003)  T(1), T(2)
               J = 1
            ENDIF

            DO I = 1, N
               FCALC = BVALU (T, BCOEF, NDIM, K, 0, X(I), INBV, WORK)
               FERR = F(I) - FCALC
               FERMAX = MAX(FERMAX, RELERR(FERR,F(I)) )
               DCALC = BVALU (T, BCOEF, NDIM, K, 1, X(I), INBV, WORK)
               DERR = D(I) - DCALC
               DERMAX = MAX(DERMAX, RELERR(DERR,D(I)) )
               IF (KPRINT.GE.3)  THEN
                  J = J + 2
                  WRITE (LUN, 2004)  X(I), T(J), T(J+1), F(I), FERR, &
                                                         D(I), DERR
               ENDIF
            end do

            IF (KPRINT.GE.3)  THEN
               J = J + 2
               WRITE (LUN, 2003)  T(J), T(J+1)
            ENDIF
            FAIL = (FERMAX.GT.TOL).OR.(DERMAX.GT.TOL)
            IF (FAIL)  IFAIL = IFAIL + 1
            IF ((KPRINT.GE.3).OR.(KPRINT.GE.2).AND.FAIL) then
               WRITE (LUN, 2005)  FERMAX, DERMAX, TOL
            end if
         ENDIF
!
!  Special check for KNOTYP=-1.
!
         IF (KNOTYP.EQ.0)  THEN
!
!  Save knot vector for next test.
!
            DO I = 1, NKNOTS
               TSAVE(I) = T(I)
            end do

         ELSE IF (KNOTYP.EQ.-1)  THEN
!
!  Check that knot vector is unchanged.
!
            TERMAX = ZERO
            DO I = 1, NKNOTS
               TERR = ABS(T(I) - TSAVE(I))
               TERMAX = MAX(TERMAX, TERR)
            end do

            IF (TERMAX.GT.TOLZ)  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.2)  WRITE (LUN, 2007)  TERMAX, TOLZ
            ENDIF
         ENDIF

      end do
!
!  PRINT SUMMARY AND TERMINATE.
!
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL

      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF

      RETURN

 1000 FORMAT ('1'//10X,'TEST PCH TO B-SPLINE CONVERTER')
 1001 FORMAT (//10X,'PCHQK5 RESULTS'/10X,'--------------')
 1010 FORMAT (/4X,'(Results should be the same for all KNOTYP values.)')
 2000 FORMAT (/4X,'KNOTYP =',I2,':  NKNOTS =',I3,',  NDIM =',I3, &
                              ',  K =',I2,',  IERR =',I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (/15X,'X',9X,'KNOTS',10X,'F',7X,'FERR',8X,'D',7X,'DERR')
 2003 FORMAT (18X,2F8.2)
 2004 FORMAT (10X,3F8.2,F10.4,1P,E10.2,0P,F10.4,1P,E10.2)
 2005 FORMAT (/5X,'Maximum relative errors:' &
             /15X,'F-error =',1P,E13.5,5X,'D-error =',E13.5 &
              /5X,'Both should be less than  TOL =',E13.5)
 2007 FORMAT (/' *** T-ARRAY MAXIMUM CHANGE =',1P,E13.5, &
                 ';  SHOULD NOT EXCEED TOLZ =',E13.5)
 3001 FORMAT (/' *** TROUBLE ***',I5,' CONVERSION TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL CONVERSION TESTS', &
              ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME CONVERSION TESTS', &
              ' ************')
END
SUBROUTINE PCHSP ( IC, VC, N, X, F, D, INCFD, WK, NWK, IERR )

!*****************************************************************************80
!
!! PCHSP: derivatives for Hermite representation of cubic spline interpolant.
!
!  Discussion:
!
!    Computes the Hermite representation of the cubic spline inter-
!    polant to the data given in X and F satisfying the boundary
!    conditions specified by IC and VC.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of C. de Boor's cubic spline
!    routine CUBSPL.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer ( kind = 4 )  IC(2), N, NWK, IERR
!        REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!  Parameters:
!
!     IC -- (input) integer ( kind = 4 ) array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor's cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!           These values will determine the cubic spline interpolant
!           with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK.LT.2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
!              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
  implicit none

  integer ( kind = 4 )  IC(2), N, INCFD, NWK, IERR
  real ( kind = 4 ) VC(2), X(*), F(INCFD,*), D(INCFD,*), WK(2,*)
  integer ( kind = 4 )  IBEG, IEND, INDEX, J, NM1
  real ( kind = 4 ) G, HALF, ONE, STEMP(3), XTEMP(4)
  SAVE HALF, ONE
  real ( kind = 4 ) PCHDF

  DATA  HALF /0.5/,  ONE /1./
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO J = 2, N
     IF ( X(J).LE.X(J-1) )  GO TO 5003
  end do

  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
  IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
  IF ( IERR.LT.0 )  GO TO 5004
!
!  FUNCTION DEFINITION IS OK.  GO ON.
!
  IF ( NWK .LT. 2*N )  GO TO 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
!
  DO J=2,N
     WK(1,J) = X(J) - X(J-1)
     WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
  end do
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
  IF ( IBEG.GT.N )  IBEG = 0
  IF ( IEND.GT.N )  IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
  IF ( (IBEG == 1).OR.(IBEG == 2) )  THEN
     D(1,1) = VC(1)
  ELSE IF (IBEG .GT. 2)  THEN
!
!  PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
!
     DO J = 1, IBEG
        INDEX = IBEG-J+1
!
!  INDEX RUNS FROM IBEG DOWN TO 1.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
     end do

     D(1,1) = PCHDF (IBEG, XTEMP, STEMP, IERR)

     IF (IERR .NE. 0)  GO TO 5009
     IBEG = 1
  end if

  IF ( (IEND == 1).OR.(IEND == 2) )  THEN
     D(1,N) = VC(2)
  ELSE IF (IEND .GT. 2)  THEN
!
!  PICK UP LAST IEND POINTS.
!
     DO J = 1, IEND
        INDEX = N-IEND+J
!
!  INDEX RUNS FROM N+1-IEND UP TO N.
!
        XTEMP(J) = X(INDEX)
        IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
     end do

     D(1,N) = PCHDF (IEND, XTEMP, STEMP, IERR)

     IF (IERR .NE. 0)  GO TO 5009
     IEND = 1
  end if
!
!  BEGIN CODING FROM CUBSPL
!
!  A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS
!  ELIMINATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!  WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!  WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
  IF (IBEG == 0)  THEN
     IF (N == 2)  THEN
!
!  NO CONDITION AT LEFT END AND N = 2.
!
        WK(2,1) = ONE
        WK(1,1) = ONE
        D(1,1) = 2.0E+00 * WK(2,2)
     ELSE
!
!  NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
!
        WK(2,1) = WK(1,3)
        WK(1,1) = WK(1,2) + WK(1,3)
        D(1,1) =((WK(1,2) + 2.0E+00 * WK(1,1))*WK(2,2)*WK(1,3) &
                          + WK(1,2)**2*WK(2,3)) / WK(1,1)
     end if
  ELSE IF (IBEG == 1)  THEN
!
!  SLOPE PRESCRIBED AT LEFT END.
!
     WK(2,1) = ONE
     WK(1,1) = 0.0E+00
  ELSE
!
!  SECOND DERIVATIVE PRESCRIBED AT LEFT END.
!
     WK(2,1) = 2.0E+00
     WK(1,1) = ONE
     D(1,1) = 3.0E+00 *WK(2,2) - HALF*WK(1,2)*D(1,1)
  end if
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
  NM1 = N-1
  IF (NM1 .GT. 1)  THEN
     DO J=2,NM1
        IF (WK(2,J-1) == 0.0E+00)  GO TO 5008
        G = -WK(1,J+1)/WK(2,J-1)
        D(1,J) = G*D(1,J-1) &
                   + 3.0E+00 *(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
        WK(2,J) = G*WK(1,J-1) + 2.0E+00 * (WK(1,J) + WK(1,J+1))
     end do
  end if
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!  (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!  IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!  SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!  AT THIS POINT.
!
  IF (IEND == 1)  GO TO 30

  IF (IEND == 0)  THEN
     IF (N == 2 .AND. IBEG == 0)  THEN
!
!  NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
!
        D(1,2) = WK(2,2)
        GO TO 30
     ELSE IF ((N == 2) .OR. (N == 3 .AND. IBEG == 0))  THEN
!
!  EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!  NOT-A-KNOT AT LEFT END POINT).
!
        D(1,N) = 2.0E+00 * WK(2,N)
        WK(2,N) = ONE
        IF (WK(2,N-1) == 0.0E+00)  GO TO 5008
        G = -ONE/WK(2,N-1)
     ELSE
!
!  NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
!  KNOT AT LEFT END POINT.
!
        G = WK(1,N-1) + WK(1,N)
!
!  DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
!
        D(1,N) = ((WK(1,N) + 2.0E+00 * G)*WK(2,N)*WK(1,N-1) &
                    + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
        IF (WK(2,N-1) == 0.0E+00)  GO TO 5008
        G = -G/WK(2,N-1)
        WK(2,N) = WK(1,N-1)
     end if
  ELSE
!
!  SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
!
     D(1,N) = 3.0E+00 *WK(2,N) + HALF*WK(1,N)*D(1,N)
     WK(2,N) = 2.0E+00
     IF (WK(2,N-1) == 0.0E+00)  GO TO 5008
     G = -ONE/WK(2,N-1)
  end if
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
  WK(2,N) = G*WK(1,N-1) + WK(2,N)
  IF (WK(2,N) == 0.0E+00)   GO TO 5008
  D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
   30 CONTINUE

  DO J=NM1,1,-1
     IF (WK(2,J) == 0.0E+00)  GO TO 5008
     D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
  end do
!
!  END CODING FROM CUBSPL.
!

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHSP', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHSP', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  CALL XERMSG ('SLATEC', 'PCHSP', 'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  CALL XERMSG ('SLATEC', 'PCHSP', 'IC OUT OF RANGE', IERR, 1)
  RETURN
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
  IERR = -7
  CALL XERMSG ('SLATEC', 'PCHSP', 'WORK ARRAY TOO SMALL', IERR, 1)
  RETURN
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES
!   ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3).
  IERR = -8
  CALL XERMSG ('SLATEC', 'PCHSP', 'SINGULAR LINEAR SYSTEM', IERR, 1)
  RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM PCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -9
  CALL XERMSG ('SLATEC', 'PCHSP', 'ERROR RETURN FROM PCHDF', IERR, 1)

  RETURN
END
FUNCTION PCHST (ARG1, ARG2)

!*****************************************************************************80
!
!! PCHST: carry out a sign test.
!
!  Discussion:
!
!    The object is to do this without multiplying ARG1*ARG2, to avoid
!    possible over/underflow problems.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Parameters:
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real ( kind = 4 ) ARG1, ARG2
  real ( kind = 4 ) pchst

  PCHST = SIGN ( 1.0E+00, ARG1 ) * SIGN ( 1.0E+00, ARG2 )

  IF ((ARG1 == 0.0E+00) .OR. (ARG2 == 0.0E+00)) then
    PCHST = 0.0E+00
  end if

  RETURN
END
SUBROUTINE PCHSW ( DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR )

!*****************************************************************************80
!
!! PCHSW limits excursion from data for PCHCS.
!
!  Discussion:
!
!    Called by  PCHCS  to adjust D1 and D2 if necessary to insure that
!    the extremum on this interval is not further than DFMAX from the
!    extreme data value.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer ( kind = 4 )  IEXTRM, IERR
!        REAL  DFMAX, D1, D2, H, SLOPE
!
!        CALL  PCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!  Parameters:
!
!     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!           the cubic determined by derivative values D1,D2.  (assumes
!           DFMAX.GT.0.)
!
!     IEXTRM -- (input) index of the extreme data value.  (assumes
!           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!           (Assumes D1*D2 .LE. 0.)
!          (output) may be modified if necessary to meet the restriction
!           imposed by DFMAX.
!
!     H -- (input) interval length.  (Assumes  H.GT.0.)
!
!     SLOPE -- (input) data slope on the interval.
!
!     IERR -- (output) error flag.  should be zero.
!           If IERR=-1, assumption on D1 and D2 is not satisfied.
!           If IERR=-2, quadratic equation locating extremum has
!                       negative discriminant (should never occur).
!
  implicit none

  integer ( kind = 4 )  IEXTRM, IERR
  real ( kind = 4 ) DFMAX, D1, D2, H, SLOPE
  real ( kind = 4 ) CP, FACT, HPHI, LAMBDA, NU, ONE, PHI, RADCAL, RHO, SIGMA
  real ( kind = 4 ) SMALL, THAT, THIRD

  SAVE ONE, FACT
  SAVE THIRD
  real ( kind = 4 ) R1MACH

  DATA ONE /1./, FACT /100./
!        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
  DATA  THIRD /0.33333/
!
!  NOTATION AND GENERAL REMARKS.
!
!  RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!  LAMBDA IS THE RATIO OF D2 TO D1.
!  THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!  PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!  WHERE  THAT = (XHAT - X1)/H .
!  THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!  SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!
!  SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!
  SMALL = FACT * R1MACH(4)
!
!  DO MAIN CALCULATION.
!
  IF (D1 == 0.0E+00)  THEN
!
!  SPECIAL CASE: D1 == ZERO .
!
!  IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
!
     IF (D2 == 0.0E+00)  GO TO 5001

     RHO = SLOPE/D2
!
!  EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
!
     IF (RHO .GE. THIRD) then
       return
     end if

     THAT = ( 2.0E+00 * ( 3.0E+00 *RHO-ONE)) / ( 3.0E+00 *( 2.0E+00 * RHO-ONE))
     PHI = THAT**2 * (( 3.0E+00 *RHO-ONE)/ 3.0E+00 )
!
!  CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
!
     IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!  TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
!
     HPHI = H * ABS(PHI)
     IF (HPHI*ABS(D2) .GT. DFMAX)  THEN
!
!  AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
!
        D2 = SIGN (DFMAX/HPHI, D2)
     end if
  ELSE

     RHO = SLOPE/D1
     LAMBDA = -D2/D1
     IF (D2 == 0.0E+00)  THEN
!
!  SPECIAL CASE: D2 == ZERO .
!
!  EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
!
        IF (RHO .GE. THIRD) then
          ierr = 0
          return
        end if

        CP = 2.0E+00 - 3.0E+00 *RHO
        NU = ONE - 2.0E+00 * RHO
        THAT = ONE / ( 3.0E+00 *NU)
     ELSE
        IF (LAMBDA .LE. 0.0E+00)  GO TO 5001
!
!  NORMAL CASE.  D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!
        NU = ONE - LAMBDA - 2.0E+00 * RHO
        SIGMA = ONE - RHO
        CP = NU + SIGMA
        IF (ABS(NU) .GT. SMALL)  THEN
           RADCAL = (NU - ( 2.0E+00 * RHO+ONE))*NU + SIGMA**2
           IF (RADCAL .LT. 0.0E+00)  GO TO 5002
           THAT = (CP - SQRT(RADCAL)) / ( 3.0E+00 *NU)
        ELSE
           THAT = ONE/( 2.0E+00 * SIGMA)
        end if
     end if
     PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!
!  CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
!
     IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!  TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
!
     HPHI = H * ABS(PHI)
     IF (HPHI*ABS(D1) .GT. DFMAX)  THEN
!
!  AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
!
        D1 = SIGN (DFMAX/HPHI, D1)
        D2 = -LAMBDA*D1
     end if
  end if

  IERR = 0

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
  IERR = -1
  CALL XERMSG ('SLATEC', 'PCHSW', 'D1 AND/OR D2 INVALID', IERR, 1)
  RETURN

 5002 CONTINUE
!     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
  IERR = -2
  CALL XERMSG ('SLATEC', 'PCHSW', 'NEGATIVE RADICAL', IERR, 1)

  RETURN
END
function r1mach ( i )

!*****************************************************************************80
!
!! R1MACH returns single precision real machine constants.
!
!  Discussion:
!
!    Assume that single precision real numbers are stored with a mantissa 
!    of T digits in base B, with an exponent whose value must lie 
!    between EMIN and EMAX.  Then for values of I between 1 and 5, 
!    R1MACH will return the following values:
!
!      R1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      R1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      R1MACH(3) = B^(-T), the smallest relative spacing.
!      R1MACH(4) = B^(1-T), the largest relative spacing.
!      R1MACH(5) = log10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 4 ) R1MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r1mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  else if ( i == 1 ) then
    r1mach = 1.1754944E-38
  else if ( i == 2 ) then
    r1mach = 3.4028235E+38
  else if ( i == 3 ) then
    r1mach = 5.9604645E-08
  else if ( i == 4 ) then
    r1mach = 1.1920929E-07
  else if ( i == 5 ) then
    r1mach = 0.3010300E+00
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  end if
 
  return
end
FUNCTION RAND ( R )

!*****************************************************************************80
!
!! RAND generates a uniformly distributed random number.
!
!  Discussion:
!
!    This pseudo-random number generator is portable among a wide
!    variety of computers.  RAND(R) undoubtedly is not as good as many
!    readily available installation dependent versions, and so this
!    routine is not recommended for widespread usage.  Its redeeming
!    feature is that the exact same random numbers (to within final round-
!    off error) can be generated from machine to machine.  Thus, programs
!    that make use of random numbers can be easily transported to and
!    checked in a new environment.
!
!    The random numbers are generated by the linear congruential
!    method described, e.g., by Knuth in Seminumerical Methods (p.9),
!    Addison-Wesley, 1969.  Given the I-th number of a pseudo-random
!    sequence, the I+1 -st number is generated from
!      X(I+1) = (A*X(I) + C) MOD M,
!    where here M = 2**22 = 4194304, C = 1731 and several suitable values
!    of the multiplier A are discussed below.  Both the multiplier A and
!    random number X are represented in double precision as two 11-bit
!    words.  The constants are chosen so that the period is the maximum
!    possible, 4194304.
!
!    In order that the same numbers be generated from machine to
!    machine, it is necessary that 23-bit integers be reducible modulo
!    2**11 exactly, that 23-bit integers be added exactly, and that 11-bit
!    integers be multiplied exactly.  Furthermore, if the restart option
!    is used (where R is between 0 and 1), then the product R*2**22 =
!    R*4194304 must be correct to the nearest integer.
!
!    The first four random numbers should be .0004127026,
!    .6750836372, .1614754200, and .9086198807.  The tenth random number
!    is .5527787209, and the hundredth is .3600893021 .  The thousandth
!    number should be .2176990509 .
!
!    In order to generate several effectively independent sequences
!    with the same generator, it is necessary to know the random number
!    for several widely spaced calls.  The I-th random number times 2**22,
!    where I=K*P/8 and P is the period of the sequence (P = 2**22), is
!    still of the form L*P/8.  In particular we find the I-th random
!    number multiplied by 2**22 is given by
!    I   =  0  1*P/8  2*P/8  3*P/8  4*P/8  5*P/8  6*P/8  7*P/8  8*P/8
!    RAND=  0  5*P/8  2*P/8  7*P/8  4*P/8  1*P/8  6*P/8  3*P/8  0
!    Thus the 4*P/8 = 2097152 random number is 2097152/2**22.
!
!    Several multipliers have been subjected to the spectral test
!    (see Knuth, p. 82).  Four suitable multipliers roughly in order of
!    goodness according to the spectral test are
!    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5
!    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5
!    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5
!    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1
!
!    In the table below LOG10(NU(I)) gives roughly the number of
!    random decimal digits in the random numbers considered I at a time.
!    C is the primary measure of goodness.  In both cases bigger is better.
!
!                   LOG10 NU(I)              C(I)
!       A       I=2  I=3  I=4  I=5    I=2  I=3  I=4  I=5
!
!    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
!    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
!    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
!    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
!    Best
!    Possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Wayne Fullerton
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) R, used for a seed value.
!    If R=0., the next random number of the sequence is generated.
!    If R .LT. 0., the last generated number will be returned for
!    possible use in a restart procedure.
!    If R .GT. 0., the sequence of random numbers will start with
!    the seed R mod 1.  This seed is also returned as the value of
!    RAND provided the arithmetic is done exactly.
!
!    Output, real ( kind = 8 ) RAND, a pseudo-random number between 0. and 1.
!
  implicit none

  integer ( kind = 4 ) ia0
  integer ( kind = 4 ) ia1
  integer ( kind = 4 ) ia1ma0
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ix0
  integer ( kind = 4 ) ix1
  integer ( kind = 4 ) iy0
  integer ( kind = 4 ) iy1
  real ( kind = 4 ) r
  real ( kind = 4 ) rand

  SAVE IA1, IA0, IA1MA0, IC, IX1, IX0

  DATA IA1, IA0, IA1MA0 /1536, 1029, 507/
  DATA IC /1731/
  DATA IX1, IX0 /0, 0/

  IF ( R.LT.0.0E+00 ) GO TO 10

  IF (R.GT.0.0E+00 ) GO TO 20
!
!  A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
!  + IA0*IX0) + IA0*IX0
!
  IY0 = IA0*IX0
  IY1 = IA1*IX1 + IA1MA0*(IX0-IX1) + IY0
  IY0 = IY0 + IC
  IX0 = MOD (IY0, 2048)
  IY1 = IY1 + (IY0-IX0)/2048
  IX1 = MOD (IY1, 2048)

 10 continue
   
  RAND = IX1*2048 + IX0
  RAND = RAND / 4194304.0E+00
  RETURN

 20 continue

  IX1 = MOD(R,1.)*4194304. + 0.5E+00
  IX0 = MOD (IX1, 2048)
  IX1 = (IX1-IX0)/2048
  GO TO 10

END
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
SUBROUTINE XERCNT ( LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL )

!*****************************************************************************80
!
!! XERCNT allows user control over handling of errors.
!
!  Discussion:
!
!    Allows user control over handling of individual errors.
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCNT.
!    If the user has provided his own version of XERCNT, he
!    can then override the value of KONTROL used in processing
!    this message by redefining its value.
!    KONTRL may be set to any value from -2 to 2.
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  CHARACTER ( len = * ) LIBRAR
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  character ( len = *  ) SUBROU

  RETURN
END
SUBROUTINE XERDMP ( )

!*****************************************************************************80
!
!! XERDMP prints the error tables, then clears them.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
  implicit none

  integer ( kind = 4 ) kount

  CALL XERSVE (' ',' ',' ',0,0,0,KOUNT)

  RETURN
END
SUBROUTINE XERHLT ( MESSG )

!*****************************************************************************80
!
!! XERHLT halts program execution after a fatal error.
!
!  Discussion:
!
!    XERHLT aborts the execution of the program.
!    The error message causing the abort is given in the calling
!    sequence, in case one needs it for printing on a dayfile,
!    for example.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!    MESSG is as in XERMSG.
!
  implicit none

  CHARACTER ( len = * ) MESSG

  STOP 1
END
SUBROUTINE XERMSG ( LIBRAR, SUBROU, MESSG, NERR, LEVEL )

!*****************************************************************************80
!
!! XERMSG processes error messages.
!
!  Discussion:
!
!    XERMSG processes a diagnostic message in a manner determined by the
!    value of LEVEL and the current value of the library error control
!    flag, KONTRL.  See subroutine XSETF for details.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Kirby Fong
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer ( kind = 4 ) value that is chosen by the library 
!             routine's author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer ( kind = 4 ).
!
!    LEVEL    An integer ( kind = 4 ) value in the range 0 to 2 that indicates
!             the level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) kdummy
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) lerr
  integer ( kind = 4 ) level
  character * ( 20 ) lfirst
  integer ( kind = 4 ) lkntrl
  CHARACTER * ( * ) LIBRAR
  integer ( kind = 4 ) llevel
  integer ( kind = 4 ) ltemp
  integer ( kind = 4 ) maxmes
  character * ( * ) messg
  integer ( kind = 4 ) mkntrl
  integer ( kind = 4 ) nerr
  character * ( * ) SUBROU
  character * ( 72 ) temp
  CHARACTER * ( 8 ) XLIBR
  character * ( 8 ) XSUBR

  LKNTRL = J4SAVE (2, 0, .FALSE.)
  MAXMES = J4SAVE (4, 0, .FALSE.)
!
!  LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!  MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!  SHOULD BE PRINTED.
!
!  WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!  CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!  AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
  IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR == 0 .OR. &
    LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
     CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
        'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
        'JOB ABORT DUE TO FATAL ERROR.', 72)
     CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
     CALL XERHLT (' ***XERMSG -- INVALID INPUT')
     RETURN
  end if
!
!  RECORD THE MESSAGE.
!
  I = J4SAVE (1, NERR, .TRUE.)
  CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!  HANDLE PRINT-ONCE WARNING MESSAGES.
!
  IF (LEVEL == -1 .AND. KOUNT.GT.1) RETURN
!
!  ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
  XLIBR  = LIBRAR
  XSUBR  = SUBROU
  LFIRST = MESSG
  LERR   = NERR
  LLEVEL = LEVEL
  CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)

  LKNTRL = MAX(-2, MIN(2,LKNTRL))
  MKNTRL = ABS(LKNTRL)
!
!  SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!  ZERO AND THE ERROR IS NOT FATAL.
!
  IF (LEVEL.LT.2 .AND. LKNTRL == 0) GO TO 30
  IF (LEVEL == 0 .AND. KOUNT.GT.MAXMES) GO TO 30
  IF (LEVEL == 1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL == 1) GO TO 30
  IF (LEVEL == 2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
!
!  ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!  MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!  AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!  IS NOT ZERO.
!
  IF (LKNTRL .NE. 0) THEN
     TEMP(1:21) = 'MESSAGE FROM ROUTINE '
     I = MIN(LEN(SUBROU), 16)
     TEMP(22:21+I) = SUBROU(1:I)
     TEMP(22+I:33+I) = ' IN LIBRARY '
     LTEMP = 33 + I
     I = MIN(LEN(LIBRAR), 16)
     TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
     TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
     LTEMP = LTEMP + I + 1
     CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!  IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!  PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!  FROM EACH OF THE FOLLOWING THREE OPTIONS.
!  1.  LEVEL OF THE MESSAGE
!         'INFORMATIVE MESSAGE'
!         'POTENTIALLY RECOVERABLE ERROR'
!         'FATAL ERROR'
!  2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!         'PROG CONTINUES'
!         'PROG ABORTED'
!  3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!      MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!      WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!         'TRACEBACK REQUESTED'
!         'TRACEBACK NOT REQUESTED'
!  NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!  EXCEED 74 CHARACTERS.
!  WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
  IF (LKNTRL .GT. 0) THEN
!
!  THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
     IF (LEVEL .LE. 0) THEN
        TEMP(1:20) = 'INFORMATIVE MESSAGE,'
        LTEMP = 20
     ELSEIF (LEVEL == 1) THEN
        TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
        LTEMP = 30
     ELSE
        TEMP(1:12) = 'FATAL ERROR,'
        LTEMP = 12
     end if
!
!  THEN WHETHER THE PROGRAM WILL CONTINUE.
!
     IF ((MKNTRL == 2 .AND. LEVEL.GE.1) .OR. &
         (MKNTRL == 1 .AND. LEVEL == 2)) THEN
        TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
        LTEMP = LTEMP + 14
     ELSE
        TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
        LTEMP = LTEMP + 16
     end if
!
!  FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
     IF (LKNTRL .GT. 0) THEN
        TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
        LTEMP = LTEMP + 20
     ELSE
        TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
        LTEMP = LTEMP + 24
     end if
     CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!  NOW SEND OUT THE MESSAGE.
!
  CALL XERPRN (' *  ', -1, MESSG, 72)
!
!  IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!  TRACEBACK.
!
  IF (LKNTRL .GT. 0) THEN
     WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
     DO I=16,22
        IF (TEMP(I:I) .NE. ' ') GO TO 20
     end do

   20    continue

     CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
     CALL FDUMP
  end if
!
!  IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
  IF (LKNTRL .NE. 0) THEN
     CALL XERPRN (' *  ', -1, ' ', 72)
     CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
     CALL XERPRN ('    ',  0, ' ', 72)
  end if
!
!  IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!  CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 continue

   IF (LEVEL.LE.0 .OR. (LEVEL == 1 .AND. MKNTRL.LE.1)) RETURN
!
!  THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!  FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!  SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
  IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
     IF (LEVEL == 1) THEN
        CALL XERPRN (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
     ELSE
        CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
     end if
     CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
     CALL XERHLT (' ')
  ELSE
     CALL XERHLT (MESSG)
  end if

  RETURN
END
SUBROUTINE XERPRN ( PREFIX, NPREF, MESSG, NWRAP )

!*****************************************************************************80
!
!! XERPRN prints error messages from XERMSG.
!
!  Discussion:
!
!    This routine sends one or more lines to each of the (up to five)
!    logical units to which error messages are to be sent.  This routine
!    is called several times by XERMSG, sometimes with a single line to
!    print and sometimes with a (potentially very long) message that may
!    wrap around into multiple lines.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Kirby Fong
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type integer ( kind = 4 ).  This argument is the 
!         number of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type integer ( kind = 4 ).  This gives the maximum 
!         size piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
  implicit none

  CHARACTER*148 CBUFF
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) idelta
  integer ( kind = 4 ) IU(5)
  integer ( kind = 4 ) lenmsg
  integer ( kind = 4 ) lpiece
  integer ( kind = 4 ) lpref
  integer ( kind = 4 ) lwrap
  character * ( * ) messg
  integer ( kind = 4 ) n
  CHARACTER*2 NEWLIN
  PARAMETER (NEWLIN = '$$')
  integer ( kind = 4 ) nextc
  integer ( kind = 4 ) npref
  integer ( kind = 4 ) NUNIT
  integer ( kind = 4 ) nwrap
  CHARACTER*(*) PREFIX

  CALL XGETUA(IU,NUNIT)
!
!  A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!  ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!  ERROR MESSAGE UNIT.
!
  N = I1MACH(4)
  DO I=1,NUNIT
     IF (IU(I) == 0) IU(I) = N
  end do
!
!  LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!  BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!  THE REST OF THIS ROUTINE.
!
  IF ( NPREF .LT. 0 ) THEN
     LPREF = LEN(PREFIX)
  ELSE
     LPREF = NPREF
  end if

  LPREF = MIN(16, LPREF)
  IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
!
!  LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!  TIME FROM MESSG TO PRINT ON ONE LINE.
!
  LWRAP = MAX(16, MIN(132, NWRAP))
!
!  SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
  LENMSG = LEN(MESSG)
  N = LENMSG
  DO I=1,N
     IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
     LENMSG = LENMSG - 1
  end do
   30 CONTINUE
!
!  IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
  IF (LENMSG == 0) THEN
     CBUFF(LPREF+1:LPREF+1) = ' '
     DO I=1,NUNIT
        WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
     end do
     RETURN
  end if
!
!  SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!  STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!  WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!  WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!  WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!  INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!  OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!  OF THE SECOND ARGUMENT.
!
!  THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!  FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!  OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!  POSITION NEXTC.
!
!  LPIECE == 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                REMAINDER OF THE CHARACTER STRING.  LPIECE
!                SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                WHICHEVER IS LESS.
!
!  LPIECE == 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                SHOULD BE INCREMENTED BY 2.
!
!  LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!  ELSE          THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
!                RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                PROPERLY HANDLES THE END CASE WHERE LPIECE  == 
!                LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                AT THE END OF A LINE.
!
  NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
  IF (LPIECE == 0) THEN
!
!  THERE WAS NO NEW LINE SENTINEL FOUND.
!
     IDELTA = 0
     LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
     IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
        DO I=LPIECE+1,2,-1
           IF (MESSG(NEXTC+I-1:NEXTC+I-1) == ' ') THEN
              LPIECE = I-1
              IDELTA = 1
              GOTO 54
           end if
        end do
     end if
   54    continue

     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSEIF (LPIECE == 1) THEN
!
!  WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!  DON'T PRINT A BLANK LINE.
!
     NEXTC = NEXTC + 2
     GO TO 50
  ELSEIF (LPIECE .GT. LWRAP+1) THEN
!
!  LPIECE SHOULD BE SET DOWN TO LWRAP.
!
     IDELTA = 0
     LPIECE = LWRAP
     DO I=LPIECE+1,2,-1
        IF (MESSG(NEXTC+I-1:NEXTC+I-1) == ' ') THEN
           LPIECE = I-1
           IDELTA = 1
           GOTO 58
        end if
     end do
   58    continue
     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSE
!
!  IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
!  WE SHOULD DECREMENT LPIECE BY ONE.
!
     LPIECE = LPIECE - 1
     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC  = NEXTC + LPIECE + 2
  end if
!
!  PRINT
!
  DO I=1,NUNIT
     WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
  end do

  IF (NEXTC .LE. LENMSG) GO TO 50

  RETURN
END
SUBROUTINE XERSVE ( LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT )

!*****************************************************************************80
!
!! XERSVE records that an error has occurred.
!
!  Discussion:
!
!    Record that this error occurred and possibly dump and clear the
!    tables.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
! *Usage:
!
!        integer ( kind = 4 )  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
!  Parameters:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
  implicit none

  integer ( kind = 4 ) lentab
  PARAMETER (LENTAB=10)

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kflag
  integer ( kind = 4 ) kount(lentab)
  integer ( kind = 4 ) kountx
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ) levtab(lentab)
  CHARACTER*8 LIB
  CHARACTER*8  LIBTAB(LENTAB)
  CHARACTER*(*) LIBRAR
  integer ( kind = 4 ) LUN(5)
  CHARACTER*20 MES
  CHARACTER*(*) MESSG
  CHARACTER*20 MESTAB(LENTAB)
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nertab(lentab)
  integer ( kind = 4 ) nmsg
  integer ( kind = 4 ) nunit
  CHARACTER*8 SUB
  CHARACTER*(*) SUBROU
  CHARACTER*8 SUBTAB(LENTAB)

  SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG

  DATA KOUNTX/0/, NMSG/0/

  IF (KFLAG.LE.0) THEN
!
!  Dump the table.
!
     IF (NMSG == 0) RETURN
!
!  Print to each unit.
!
     CALL XGETUA (LUN, NUNIT)

     DO KUNIT = 1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT == 0) IUNIT = I1MACH(4)
!
!  Print the table header.
!
        WRITE (IUNIT,9000)
!
!  Print body of table.
!
        DO I = 1,NMSG
           WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
             NERTAB(I),LEVTAB(I),KOUNT(I)
        end do
!
!  Print number of other errors.
!
        IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
        WRITE (IUNIT,9030)

     end do
!
!  Clear the error tables.
!
     IF (KFLAG == 0) THEN
        NMSG = 0
        KOUNTX = 0
     end if
  ELSE
!
!  PROCESS A MESSAGE.
!  SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!  OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
     LIB = LIBRAR
     SUB = SUBROU
     MES = MESSG
     DO I = 1,NMSG
        IF (LIB == LIBTAB(I) .AND. SUB == SUBTAB(I) .AND. &
           MES == MESTAB(I) .AND. NERR == NERTAB(I) .AND. &
           LEVEL == LEVTAB(I)) THEN
              KOUNT(I) = KOUNT(I) + 1
              ICOUNT = KOUNT(I)
              RETURN
        end if
     end do

     IF (NMSG.LT.LENTAB) THEN
!
!  Empty slot found for new message.
!
        NMSG = NMSG + 1
        LIBTAB(I) = LIB
        SUBTAB(I) = SUB
        MESTAB(I) = MES
        NERTAB(I) = NERR
        LEVTAB(I) = LEVEL
        KOUNT (I) = 1
        ICOUNT    = 1
     ELSE
!
!  Table is full.
!
        KOUNTX = KOUNTX+1
        ICOUNT = 0
     end if
  end if
  RETURN
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
     ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
     '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
END
SUBROUTINE XGETF ( KONTRL )

!*****************************************************************************80
!
!! XGETF returns the current value of the error control flag.
!
!  Discussion:
!
!    XGETF returns the current value of the error control flag
!    in KONTRL.  See subroutine XSETF for flag value meanings.
!    KONTRL is an output parameter only.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KONTRL, the current value of the
!    error control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) kontrl

  KONTRL = J4SAVE(2,0,.FALSE.)

  RETURN
END
SUBROUTINE XGETUA ( IUNITA, N )

!*****************************************************************************80
!
!! XGETUA returns error unit numbers.
!
!  Discussion:
!
!    XGETUA may be called to determine the unit number or numbers
!    to which error messages are being sent.
!    These unit numbers may have been set by a call to XSETUN,
!    or a call to XSETUA, or may be a default value.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) iunita(5)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) n

  N = J4SAVE(5,0,.FALSE.)

  DO I=1,N
    INDEX = I+4
    IF (I == 1) then
      INDEX = 3
    end if
    IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
  end do

  RETURN
END
SUBROUTINE XSETF ( KONTRL )

!*****************************************************************************80
!
!! XSETF sets the error control flag.
!
!  Discussion:
!
!    XSETF sets the error control flag value to KONTRL.
!    (KONTRL is an input parameter only.)
!    The following table shows how each message is treated,
!    depending on the values of KONTRL and LEVEL.  (See XERMSG
!    for description of LEVEL.)
!
!    If KONTRL is zero or negative, no information other than the
!    message itself (including numeric values, if any) will be
!    printed.  If KONTRL is positive, introductory messages,
!    trace-backs, etc., will be printed in addition to the message.
!
!          ABS(KONTRL)
!    LEVEL        0              1              2
!    value
!      2        fatal          fatal          fatal
!
!      1     not printed      printed         fatal
!
!      0     not printed      printed        printed
!
!     -1     not printed      printed        printed
!                              only           only
!                              once           once
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KONTRL, the error control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kontrl
  CHARACTER ( len = 8 ) XERN1

  IF ( ABS ( KONTRL ) .GT. 2) THEN
    WRITE (XERN1, '(I8)') KONTRL
    CALL XERMSG ('SLATEC', 'XSETF', &
      'INVALID ARGUMENT = ' // XERN1, 1, 2)
    RETURN
  ENDIF

  JUNK = J4SAVE(2,KONTRL,.TRUE.)

  RETURN
END
