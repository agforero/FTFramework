C$Procedure ZZELNAXX ( ellipse normal axis intercepts )
 
      SUBROUTINE ZZELNAXX ( A, B, LAT, XXPT, YXPT )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Given semi-axis lengths of an ellipse and a planetodetic
C     latitude, find the X-axis and Y-axis intercepts of the normal
C     line corresponding to the point at that latitude in the right
C     (positive X) half-plane.
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
C     DSK
C     ELLIPSE
C     GEOMETRY
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      XXPT
      DOUBLE PRECISION      YXPT

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     A          I   Ellipse's semi-axis length in the X direction.
C     B          I   Ellipse's semi-axis length in the Y direction.
C     LAT        I   Planetodetic latitude in radians.
C     XXPT       O   Normal line intercept on the X-axis.
C     YXPT       O   Normal line intercept on the Y-axis.
C
C$ Detailed_Input
C
C     A,
C     B          are, respectively, an ellipse's semi-axis lengths
C                in the X and Y directions. The ellipse lies in
C                two-dimensional Euclidean space.
C
C                Any of the relationships
C
C                   A < B,  A = B,  A > B
C
C                are allowed.
C
C
C     LAT        is a planetodetic latitude corresponding to 
C                some point P on the ellipse defined by A and B,
C                where 
C
C                   P(1) >= 0
C
C                Units are radians.
C
C$ Detailed_Output
C
C     XXPT       is the X-intercept of a line passing through the
C                ellipse and normal to the ellipse at the point having
C                the given latitude.
C
C                If LAT = 0 radians, XXPT is defined to be 
C               
C                        2
C                   A - B / A
C
C
C     YXPT       is the Y-intercept of a line passing through the
C                ellipse and normal to the ellipse at the point having
C                the given latitude.
C
C                If LAT = pi/2 radians, YXPT is defined to be 
C
C                         2
C                    B - A / B
C
C                If LAT = -pi/2 radians, YXPT is defined to be 
C
C                         2
C                   -B + A / B
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C 
C     1)  If the X semi-axis length A is non-positive, the
C         error SPICE(NONPOSITIVEAXIS) will be signaled.
C
C     2)  If the Y semi-axis length B is non-positive, the
C         error SPICE(NONPOSITIVEAXIS) will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C
C     The outputs of this routine can be used to determine the
C     apex of a cone of constant planetodetic latitude.
C
C
C     The axis intercepts computed by this routine are defined
C     as described below.
C
C     Let M be the normal vector's slope. Let DELTAX be
C     X - XXPT. Then, given that the normal direction is
C     parallel to 
C
C             2     2
C        ( X/A,  Y/B  )
C
C     we have, for non-zero DELTAX, 
C
C        Y = M * DELTAX
C
C               2   2
C          = ( A / B )  * ( Y / X ) * ( X - XXPT )
C
C     So
C             2    2
C        X ( B  / A ) = X - XXPT  
C
C     and
C                         2   2
C        XXPT = X ( 1 - (B / A ) )
C
C
C     The X intercept is a linear function of the X-coordinate
C     of the point on the ellipsoid. Define the intercept for
C     
C        X = 0
C
C     as the limit as X -> 0 of the expression for XXPT above. 
C     The expression is continuous and defined at X = 0, so 
C     we can just evaluate the expression at X = 0 as for other
C     values of X. 
C
C
C     Using the definition of M above, we have for non-zero X:
C
C        Y - YXPT = M * X
C
C     or
C
C        YXPT = Y - M*X
C
C                      2   2
C             = Y - ( A / B ) (Y/X) * X
C
C                           2   2
C             = Y (  1 - ( A / B )  )
C
C     As above, we define YXPT at X = 0 to be the limit as X -> 0
C     of the expression above, which is 
C
C                    2   2
C        B (  1 - ( A / B )  ) 
C
C
C$ Examples
C
C     See usage in ZZRYXPDT.
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 09-MAR-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find x and y intercepts of ellipse normal line
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      EPT    ( 3 )
      DOUBLE PRECISION      NORMAL ( 3 )


      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF (  ( A .LE. 0.D0 ) .OR. ( B .LE. 0.D0 )  ) THEN

         CALL CHKIN  ( 'ZZELNAXX'                             )
         CALL SETMSG ( 'Semi-axis lengths were A = #; B = #. '
     .   //            'Both must be positive.'               )
         CALL ERRDP  ( '#',  A                                )
         CALL ERRDP  ( '#',  B                                )
         CALL SIGERR ( 'SPICE(NONPOSITIVEAXIS)'               )
         CALL CHKOUT ( 'ZZELNAXX'                             )
         RETURN

      END IF

C
C     Find the point lying on the positive X portion of the ellipsoid
C     and having the input planetodetic latitude.
C
C     To start, create a normal vector pointing in the direction
C     indicated by the latitude. We'll work in three dimensions in
C     order to take advantage of existing code. The third coordinates
C     of all participating vectors will be zero.
C
      NORMAL(1) = COS( LAT )
      NORMAL(2) = SIN( LAT )
      NORMAL(3) = 0.D0

      CALL EDNMPT ( A, B, B, NORMAL, EPT )
 
C
C     Compute the X-axis and Y-axis intercepts of the line
C     passing through EPT and parallel to a normal vector
C     at EPT. Refer to the Particulars above for details.
C
      XXPT = (  1.D0 - ( (B/A)**2 )  ) * EPT(1) 

      YXPT = (  1.D0 - ( (A/B)**2 )  ) * EPT(2) 

      RETURN
      END
