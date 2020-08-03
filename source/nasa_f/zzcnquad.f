C$Procedure  ZZCNQUAD ( Solve quadratic equation for cone intercept )
 
      SUBROUTINE ZZCNQUAD ( A, B, C, UB, N, R1, R2 ) 
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Solve a quadratic equation using an upper bound for the absolute
C     value of the roots. Only real roots are computed. This routine
C     addresses the case of a small, non-zero leading coefficient.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     EQUATION
C     MATH
C     QUADRATIC
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      DOUBLE PRECISION      UB
      INTEGER               N
      DOUBLE PRECISION      R1
      DOUBLE PRECISION      R2
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     A,         
C     B,
C     C          I   Coefficients of quadratic equation.
C     UB         I   Upper bound on magnitude of roots.
C     N          O   Number of real roots within bound.
C     R1         O   Root of smaller absolute value.
C     R2         O   Root of larger absolute value.
C     BIG        P   Limit on magnitude of coefficients and UB.
C 
C$ Detailed_Input
C
C     A,
C     B,
C     C          are the real coefficients of a quadratic equation.
C                The equation is
C
C                      2
C                   A X   +  B X  +  C  =  0
C     
C
C     UB         is an upper bound on the absolute value of real
C                roots to be found by this routine. Roots having
C                absolute value larger than UB are not returned.
C
C$ Detailed_Output
C
C
C     N          is the number of real roots found that satisfy the
C                magnitude bound constraint.
C
C                Degenerate cases: if A = B = 0, then
C
C                   If C  = 0, N is set to -1; R1 and R2 are set to
C                   zero.
C
C                   If C != 0, N is set to -2; R1 and R2 are set to
C                   zero.
C
C     R1,
C     R2        are roots returned in increasing order of absolute
C               value. If there is one real root of multiplicity 2, 
C               both R1 and R2 are set to that root, provided that
C               the absolute value of the root does not exceed UB.
C
C               If the absolute value of a root exceeds UB the 
C               corresponding output argument is set to zero. If the
C               roots are complex, both R1 and R2 are set to zero.
C
C$ Parameters
C
C     BIG            is a limit on the absolute value of the input
C                    coefficients and on the magnitude limit for the
C                    roots. BIG is set to
C
C                       SQRT( DPMAX() ) / 100
C
C$ Exceptions
C
C     1)  If any of the input coefficients have absolute value larger
C         than the parameter BIG, the error is diagnosed by a routine
C         in the call tree of this routine.
C
C     2)  If UB is non-positive or larger than the parameter BIG, the
C         error is diagnosed by a routine in the call tree of this
C         routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports INCNSG.
C
C     This routine makes use of ZZBQUAD; it calls ZZBQUAD for
C     cases that can be handled accurately by that routine.
C     Other cases are handled in-line.
C
C$ Examples
C 
C     See usage in INCNSG.
C
C$ Restrictions
C
C     1)  This is a private routine; it should not be called by
C         non-SPICELIB code.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 23-JAN-2017 (NJB)
C     
C        Previous version 22-JUN-2016 (NJB)
C
C-&

C$ Index_Entries
C
C     solve quadratic equation for cone intercept
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      TOUCHD

      LOGICAL               RETURN

C
C     Local parameters
C
      DOUBLE PRECISION      SMALL
      PARAMETER           ( SMALL = 1.D-8 )
      
C
C     Local variables
C
      DOUBLE PRECISION      COEFFS ( 3 )
      DOUBLE PRECISION      INV1
      DOUBLE PRECISION      INV2
      DOUBLE PRECISION      INVUB
      DOUBLE PRECISION      MAXMAG

      INTEGER               I
      INTEGER               MAXIX
      INTEGER               NX

      LOGICAL               FIRST

C
C     Saved values
C     
      SAVE                  FIRST
      SAVE                  INVUB

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /
      DATA                  INVUB / -1.D0  /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZCNQUAD' )

C
C     On the first pass, set the upper bound for the reciprocal
C     solution.
C
      IF ( FIRST ) THEN

         INVUB = SQRT( DPMAX() ) / 200.D0

         FIRST = .FALSE.

      END IF

C
C     Handle the degenerate cases first.
C
      IF ( ( A .EQ. 0.D0 ) .AND. ( B .EQ. 0.D0 )  ) THEN

         R1 = 0.D0
         R2 = 0.D0

         IF ( C .EQ. 0.D0 ) THEN
            
            N = -1
         ELSE
            N = -2
         END IF

         CALL CHKOUT ( 'ZZCNQUAD' )
         RETURN

      END IF

C
C     Scale the input coefficients.
C     
      MAXMAG = MAX( ABS(A), ABS(B), ABS(C) )

      COEFFS(1) = TOUCHD(  A / MAXMAG )
      COEFFS(2) = TOUCHD(  B / MAXMAG )
      COEFFS(3) = TOUCHD(  C / MAXMAG )

C
C     Identify the coefficient of largest magnitude.
C
      MAXIX = 1

      DO I = 2, 3

         IF ( ABS(COEFFS(I)) .GT. ABS(COEFFS(MAXIX)) ) THEN
C
C           Record the index of the maximum magnitude.
C
            MAXIX  = I

         END IF

      END DO

C
C     Make sure the value of maximum magnitude is +/- 1.
C
      COEFFS(MAXIX) = SIGN( 1.D0, COEFFS(MAXIX) )

C
C     Find roots in a manner suited to the coefficients we have.
C
      IF (       (  ABS( COEFFS(1) ) .GE.  SMALL  ) 
     .     .OR.  (       COEFFS(1)   .EQ.  0.D0   )  ) THEN
C
C        This is a numerically well-behaved case. Delegate the
C        job to ZZBQUAD.
C        
         CALL ZZBQUAD ( COEFFS(1), COEFFS(2), COEFFS(3), UB,
     .                  N,         NX,        R1,        R2  )
         


      ELSE IF (  ABS( COEFFS(3) ) .GE.  SMALL  ) THEN

C
C        The zero-order coefficient has magnitude >= SMALL.
C
C        The original equation
C
C              2
C           a x  + b x + c = 0
C
C        can be replaced by
C
C              2
C           c y  + b y + a = 0
C
C        where 
C
C           y = 1/x
C
C        Here 
C
C          |c| >= SMALL
C          |c| <= 1
C 
C          |a|  < SMALL
C
C
C        Because the quadratic coefficient is bounded away from zero,
C        the roots of the reciprocal equation are not in danger of
C        overflowing. So we can safely solve for 1/x. We might have
C        complex roots; these are rejected.
C  
C        The roots of the transformed equation don't have a maximum
C        magnitude restriction imposed by UB. We set the upper bound
C        to a value that ZZBQUAD will allow.
C
         CALL ZZBQUAD ( COEFFS(3), COEFFS(2), COEFFS(1), INVUB, 
     .                  N,         NX,        INV1,      INV2  )

         IF ( N .EQ. 1 ) THEN
C
C           We have one real root. Make sure we can invert it.
C           
            IF (  ABS( INV1 * UB )  .GE.  1.D0  ) THEN
C
C
C              |1/INV1| <= UB
C
C
               R1 = 1.D0 / INV1

            ELSE
C
C              There are no real roots having magnitude within the
C              bound.
C
               N  = 0

            END IF

C
C           There is no second root.
C
            R2 = 0.D0

         ELSE IF ( N .EQ. 2 ) THEN
C
C           We have two real roots. The one of larger magnitude is
C           the second one. The reciprocal of this root will be
C           the smaller root of the original equation, as long
C           as the reciprocal is within bounds.
C           
            IF (  ABS( INV2 * UB )  .GE.  1.D0  ) THEN
C
C
C              |1/INV2| <= UB
C
C
               R1 = 1.D0 / INV2

C
C              Proceed to the first root of the transformed equation.
C
               IF (  ABS( INV1 * UB )  .GE.  1.D0  ) THEN
C
C
C                 |1/INV1| <= UB
C
C
                  R2 = 1.D0 / INV1

               ELSE
C
C                 Only the second root qualifies for inversion.
C
                  N  = 1
                  
                  R2 = 0.D0

               END IF


            ELSE
C
C              The reciprocal of the larger root is too big; the
C              reciprocal of the smaller root will be even larger.
C              There are no real roots having magnitude within the
C              bound.
C               
               N  = 0

               R1 = 0.D0
               R2 = 0.D0

            END IF


         ELSE
C
C           We have no viable roots of the transformed equation, so
C           we have no viable roots of the original one.
C
            N  = 0

            R1 = 0.D0
            R2 = 0.D0

         END IF



      ELSE

C
C        The linear coefficient B has the greatest magnitude, which
C        is 1. The quadratic coefficient A is "small":  0 < |A| < 1.D-8.
C        The zero-order coefficient is "small" as well.
C
C        It will be convenient to make B equal to 1; do this now.
C
         IF ( B .LT. 0.D0 ) THEN

            COEFFS(1) = -COEFFS(1)
            COEFFS(2) = -COEFFS(2)
            COEFFS(3) = -COEFFS(3)

         END IF

C
C        In this case we use a low-order Taylor expansion about
C        x = 0 for the square root term of the formula for the roots:
C
C                                  inf
C                                  __
C                           1/2    \    (k)     k
C           T(x) = ( 1 + x )    =  /_  f   (0) x / (k!)
C
C                                  k=0
C
C
C                              2      3         4
C                =  1 + x/2 - x /8 + x /16 + O(x )
C
C                
C        Apply this formula to that for the solution having the 
C        positive square root term. Here let `x' be
C
C
C                   2
C           -4ac / b
C
C        which equals
C
C           -4ac
C
C        since we've set b = 1.
C
C
C        Then the root is
C
C                              
C                  -1 + sqrt( 1 - 4ac )
C           x  =   --------------------
C            1             2a
C
C
C                                      2 2          3
C                  -1 + ( 1 - 2ac - 16a c /8 + O((ac)) )
C              =   -------------------------------------
C                                   2a
C
C        Discarding the high-order terms in a, we have
C
C           
C           x  ~=  ( -1 + 1 - 2ac ) / 2a  =  -c
C            1
C
C        Similarly, we have
C
C
C           x  ~=  ( -1 - 1 + 2ac ) / 2a  =  ( ac - 1 )/a  = c - 1/a
C            2
C
C
C        Based on the conditions that got us here, we know 
C
C           |c| < 1
C
C           |c - 1/a| ~= |1/a| > 1.e8
C        

         N  = 0
         R1 = 0.D0
         R2 = 0.D0

         IF ( ABS( COEFFS(3) ) .LE. UB ) THEN

            R1 = -COEFFS(3)

            N  = 1

            IF (        ABS( COEFFS(1)*COEFFS(3) - 1.D0 ) 
     .            .LT.  ABS( COEFFS(1)*UB               )  ) THEN

               R2 = COEFFS(3) - 1.D0/COEFFS(1)

               N  = 2

            END IF

         END IF

      END IF

      CALL CHKOUT ( 'ZZCNQUAD' )
      RETURN
      END
