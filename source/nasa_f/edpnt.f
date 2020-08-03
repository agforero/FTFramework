C$Procedure EDPNT ( Ellipsoid point  )
 
      SUBROUTINE EDPNT ( P, A, B, C, EP )

C$ Abstract
C
C     Scale a point so that it lies on the surface of a specified
C     triaxial ellipsoid that is centered at the origin and aligned
C     with the Cartesian coordinate axes.
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
C     ELLIPSOID
C     GEOMETRY
C     MATH
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      P  ( 3 )
      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      DOUBLE PRECISION      EP ( 3 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   A point in three-dimensional space.
C     A          I   Semi-axis length in the X direction.
C     B          I   Semi-axis length in the Y direction.
C     C          I   Semi-axis length in the Z direction.
C     EP         O   Point on ellipsoid.
C
C$ Detailed_Input
C
C     P          is a non-zero point in three-dimensional space.
C
C     A,
C     B,
C     C          are, respectively, the semi-axis lengths of a triaxial
C                ellipsoid in the X, Y, and Z directions. The axes of
C                the ellipsoid are aligned with the axes of the
C                Cartesian coordinate system.
C
C$ Detailed_Output
C
C     EP         is the result of scaling the input point P so that
C                it lies on the surface of the triaxial ellipsoid
C                defined by the input semi-axis lengths.
C                
C$ Parameters
C
C     None.
C
C$ Exceptions
C 
C     1)  If any of the target ellipsoid's semi-axis lengths is
C         non-positive, the error SPICE(INVALIDAXES) is signaled.
C
C     2)  If P is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C
C     3)  If the level surface parameter of the input point
C         underflows, the error SPICE(POINTTOOSMALL) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine efficiently computes the ellipsoid surface point
C     corresponding to a specified ray emanating from the origin.
C     Practical examples of this computation occur in the SPICELIB
C     routines LATSRF and SRFREC.
C     
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C
C     1)  Find the surface intercept point on an ellipsoid having radii
C
C            ( 3, 2, 1 )
C
C         of the ray emanating from the origin and having direction
C         vector
C
C            ( 1, 1, 1 )
C
C
C     Example code begins here.
C
C
C           PROGRAM EX1
C           IMPLICIT NONE
C
C           CHARACTER*(*)         FMT1
C           PARAMETER           ( FMT1 = '(A,3F17.14)' )
C
C           DOUBLE PRECISION      A
C           DOUBLE PRECISION      B
C           DOUBLE PRECISION      C
C           DOUBLE PRECISION      V      ( 3 )
C           DOUBLE PRECISION      EP     ( 3 )
C           DOUBLE PRECISION      LEVEL
C
C           A = 3.D0
C           B = 2.D0
C           C = 1.D0
C
C           CALL VPACK ( 1.D0, 1.D0, 1.D0, V )
C
C           CALL EDPNT ( V, A, B, C, EP )
C
C           WRITE (*,FMT1) 'EP    = ', EP
C     C
C     C     Verify that EP is on the ellipsoid.
C     C
C           LEVEL = (EP(1)/A)**2 + (EP(2)/B)**2 + (EP(3)/C)**2
C
C           WRITE (*,FMT1) 'LEVEL = ', LEVEL
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output was:
C
C        EP    =  0.85714285714286 0.85714285714286 0.85714285714286
C        LEVEL =  1.00000000000000
C       
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.0.0, 19-APR-2016 (NJB) (EDW)
C
C-&
 
C$ Index_Entries
C
C     scale point to lie on ellipsoid 
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      TOUCHD

      LOGICAL               FAILED
      LOGICAL               VZERO

C
C     Local variables
C 
      DOUBLE PRECISION      LEVEL
      DOUBLE PRECISION      SQ

C
C     Use discovery check-in.
C
      IF ( FAILED() ) THEN
         RETURN
      END IF

      IF (      ( A .LE. 0.D0 )
     .     .OR. ( B .LE. 0.D0 )
     .     .OR. ( C .LE. 0.D0 ) ) THEN

         CALL CHKIN  ( 'EDPNT'                            )
         CALL SETMSG ( 'Ellipsoid radii must be strictly '
     .   //            'positive but are (#, #, #).'      )
         CALL ERRDP  ( '#', A                             )
         CALL ERRDP  ( '#', B                             )
         CALL ERRDP  ( '#', C                             )
         CALL SIGERR ( 'SPICE(INVALIDRADII)'              )
         CALL CHKOUT ( 'EDPNT'                            )
         RETURN

      END IF

C
C     The input point must be non-zero, or we can't scale it
C     to the ellipsoid.
C     
      IF ( VZERO(P) ) THEN

         CALL CHKIN  ( 'EDPNT'                            )
         CALL SETMSG ( 'Input point was the zero vector. '
     .   //            'A non-zero vector is required.'   )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                )
         CALL CHKOUT ( 'EDPNT'                            )
         RETURN

      END IF

C
C     Find the level surface parameter of the input point with respect
C     to the scaled ellipsoid.
C
      LEVEL = TOUCHD( ( P(1)/A )**2 + ( P(2)/B )**2 + ( P(3)/C )**2  )

      IF ( LEVEL .LE. 0.D0 ) THEN
C
C        We expect that LEVEL will be non-negative, but it could
C        be zero. We check for negative values as a precaution.
C
         CALL CHKIN  ( 'EDPNT'                                    )
         CALL SETMSG ( 'Input point''s level surface parameter '
     .   //            'was non-positive. The point is too close '
     .   //            'to the origin to be scaled to the '
     .   //            'ellipsoid. The point was (#, #, #).'      )
         CALL ERRDP  ( '#', P(1)                                  )
         CALL ERRDP  ( '#', P(2)                                  )
         CALL ERRDP  ( '#', P(3)                                  )
         CALL SIGERR ( 'SPICE(POINTTOOSMALL)'                     )
         CALL CHKOUT ( 'EDPNT'                                    )
         RETURN

      END IF

C
C     Scale the point to one for which the level surface parameter is 1.
C
      SQ    = SQRT(LEVEL) 

      EP(1) = P(1) / SQ
      EP(2) = P(2) / SQ
      EP(3) = P(3) / SQ
     
      RETURN
      END 
