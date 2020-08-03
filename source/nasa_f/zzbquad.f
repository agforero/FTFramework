C$Procedure  ZZBQUAD ( Solve quadratic equation with bounds )
 
      SUBROUTINE ZZBQUAD ( A, B, C, UB, N, NX, R1, R2 )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Solve a quadratic equation using an upper bound for the absolute
C     value of the roots. Only real roots are computed.
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
      INTEGER               NX
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
C     NX         O   Number of real roots exceeding bound.
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
C                   If C  = 0, N is set to -1; NX is set to 0; R1 and
C                   R2 are set to zero.
C
C                   If C != 0, N is set to -2; NX is set to 0; R1 and
C                   R2 are set to zero.
C
C
C     NX        is the number of real roots having absolute values that
C               exceed the bound. 
C
C               If the roots are complex, both N and NX are set to zero.
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
C         than the parameter BIG, the error SPICE(VALUEOUTOFRANGE) is
C         signaled.
C
C     2)  If UB is non-positive or larger than the parameter BIG, the
C         error SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     See ZZCNQUAD for a more robust implementation that makes use
C     of this one.
C
C$ Examples
C 
C     See usage in ZZCNQUAD
C
C$ Restrictions
C
C     1)  This routine may suffer from loss of precision for coefficient
C         sets having small, non-zero leading coefficients. See ZZCNQUAD
C         for a more robust implementation.
C
C     2)  This is a private routine; it should not be called by
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
C-    SPICELIB Version 1.0.0, 21-JUN-2016 (NJB)
C
C        Made minor updates to header comments.
C
C        24-SEP-2014 (NJB)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     solve quadratic equation with bounds on roots
C
C-&
 
C$ Revisions
C
C     None.
C
C-& 

C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      TOUCHD

      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      BIG
      DOUBLE PRECISION      DENOM
      DOUBLE PRECISION      DSCRIM
      DOUBLE PRECISION      NUM1
      DOUBLE PRECISION      NUM2
      DOUBLE PRECISION      SQDISC

      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  BIG
      SAVE                  FIRST
      
C
C     Initial values
C
      DATA                  BIG   / 0.D0   /
      DATA                  FIRST / .TRUE. /

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF ( FIRST ) THEN
         
         BIG   = SQRT( DPMAX() ) / 100
         FIRST = .FALSE.

      END IF

C
C     Set invalid counts to start out. Initialize R1 and R2.
C
      N   = -3
      NX  = -3
      R1  = 0.D0
      R2  = 0.D0

C
C     Reject all large magnitude coefficients.
C
      IF (     ( ABS(A) .GT. BIG ) 
     .    .OR. ( ABS(B) .GT. BIG )
     .    .OR. ( ABS(C) .GT. BIG ) ) THEN

         CALL CHKIN  ( 'ZZBQUAD'                                      )
         CALL SETMSG ( 'Coefficients must have magnitude less than '
     .   //            'or equal to #, but were A = #; B = #; C = #.' )
         CALL ERRDP  ( '#',  BIG                                      )
         CALL ERRDP  ( '#',  A                                        )
         CALL ERRDP  ( '#',  B                                        )
         CALL ERRDP  ( '#',  C                                        )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                       )
         CALL CHKOUT ( 'ZZBQUAD'                                      )
         RETURN

      END IF

C
C     Reject large magnitude upper bounds as well.
C
      IF ( ABS(UB) .GT. BIG ) THEN

         CALL CHKIN  ( 'ZZBQUAD'                                      )
         CALL SETMSG ( 'Upper bounds must have magnitude less than '
     .   //            'or equal to #, but was #.'                    )
         CALL ERRDP  ( '#',  BIG                                      )
         CALL ERRDP  ( '#',  UB                                       )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                       )
         CALL CHKOUT ( 'ZZBQUAD'                                      )
         RETURN

      END IF

C
C     The upper bound must be positive.
C
      IF ( UB .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZBQUAD'                                 )
         CALL SETMSG ( 'Upper bound must be positive but was #.' )
         CALL ERRDP  ( '#',  UB                                  )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                  )
         CALL CHKOUT ( 'ZZBQUAD'                                 )
         RETURN

      END IF

C
C     Handle the degenerate cases first.
C
      IF ( A .EQ. 0.D0 ) THEN

         IF ( B .EQ. 0.D0 ) THEN
C
C           The equation is of the form
C
C              C = 0
C
            IF ( C .EQ. 0.D0 ) THEN
C
C              The equation is satisfied for all real numbers.
C               
               N  = -1
               NX =  0

            ELSE
C
C              There are no solutions.
C
               N  = -2
               NX =  0
            
            END IF

         ELSE
C
C           The equation is first-order:
C
C              B*X + C = 0
C
C           In this branch, B is non-zero.
C
            IF ( ABS(C) .LE. ABS(UB * B) ) THEN

               N   =   1
               NX  =   0

               R1  =  -C / B
               R2  =   R1      

            ELSE
C              
C              The magnitude of the solution is too large.
C              
               N   =  0
               NX  =  1

            END IF

         END IF

      
      ELSE
C
C        The leading coefficient of the equation is non-zero.
C
C        We can safely compute the discriminant now, due the
C        check we've already performed.
C
         DSCRIM = TOUCHD(  ( B*B ) - ( 4*A*C )  )

         IF ( DSCRIM .LT. 0.D0 ) THEN
C
C           We have complex roots, so we're done.
C
            N  = 0
            NX = 0

         ELSE IF ( DSCRIM .EQ. 0.D0 ) THEN
C
C           We have a single real root of multiplicity 2.
C
C           Compare the magnitude of the root to the upper bound.
C
            NUM1   =  -B
            DENOM  =   2 * A

            IF ( ABS(NUM1) .GE. ABS(DENOM * UB) ) THEN
C
C              The root is too large; we won't compute it.
C
               N   =  0
               NX  =  1

            ELSE
C
C              Set both roots to the same value. In this branch,
C              A is non-zero.
C
               N   =  1
               NX  =  0

               R1  =  ( NUM1 / A ) / 2
               R2  =  R1

            END IF


         ELSE
C
C           We have two nominally distinct real roots. Whether
C           they're distinct double precision numbers depends
C           on the relative magnitudes of A and DSCRIM.
C
            DENOM  = 2 * A
            SQDISC = SQRT( DSCRIM )

            IF ( B .GT. 0.D0 ) THEN

               NUM2 =  -B - SQDISC
               NUM1 =  -B + SQDISC
            ELSE
               NUM2 =  -B + SQDISC
               NUM1 =  -B - SQDISC
            END IF

C
C           See whether the root of larger magnitude is computable.
C
            IF ( ABS(NUM2) .LE. ABS(UB * DENOM) ) THEN
C
C              The root is computable. 

               N  = 2
               NX = 0
C
C              In this branch, A is non-zero.
C
               R2 = ( NUM2 / A ) / 2

               IF ( ABS(R2) .GT. 0.D0 ) THEN
C
C                 Compute R1 using R2 and C; this avoids loss
C                 of precision that may occur when NUM1 is computed.
C                
C                 We know R1 has smaller magnitude than R2 and R2
C                 is computable, and we know A is non-zero, so R1
C                 can be computed without a divide-by-zero error,
C                 and it is computable as long as no intermediate
C                 results overflow. The bounds on A and R2 ensure
C                 that A*R2 is computable.
C
                  R1 =  C / ( A * R2 )

               ELSE
C
C                 The root of larger magnitude has magnitude 0. This
C                 doesn't leave many possible values for the root of
C                 smaller magnitude.
C
                  R1 = 0.D0

               END IF 


            ELSE
C
C              The root of larger magnitude is not computable.
C              Check the root of smaller magnitude.
C
               IF ( ABS(NUM1) .LE. ABS(UB * DENOM) ) THEN
C
C                 The root is computable.
C
                  N  = 1
                  NX = 1

                  R1 = ( NUM1 / A ) / 2

               ELSE
C
C                 Neither root is computable.
C
                  N  = 0
                  NX = 2
 
               END IF
            
            END IF

         END IF

      END IF 
      
      RETURN
      END
