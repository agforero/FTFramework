C$Procedure ZZPDPLTC (Planetodetic coordinates, point latitude check)
 
      LOGICAL FUNCTION ZZPDPLTC ( RE, F, P, LAT )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Indicate whether a given point on a planetodetic latitude cone
C     has the correct latitude sign.
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
C     COORDINATES
C     GEOMETRY
C     PLANETODETIC
C     LATITUDE
C     MATH
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      RE
      DOUBLE PRECISION      F
      DOUBLE PRECISION      P   ( 3 )
      DOUBLE PRECISION      LAT

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     RE         I   Equatorial radius.
C     F          I   Flattening coefficient.
C     P          I   Three-dimensional point.
C     LAT        I   Planetodetic latitude.
C
C     The function returns .TRUE. if the sign of the planetodetic
C     latitude of the input point matches that of the input 
C     planetodetic latitude.
C     
C$ Detailed_Input
C
C     RE,
C     F          are, respectively, the equatorial radius
C                and flattening coefficient of a biaxial 
C                spheroid. 
C
C                The polar radius RP of the spheroid is
C                
C                   RP = RE * ( 1 - F )
C
C                RP may be less than, equal to, or greater than RE.
C
C
C     P          is a point (equivalently, a vector) in
C                three-dimensional space. P is expressed in Cartesian
C                coordinates. P must lie on the planetodetic
C                latitude cone corresponding to LAT. 
C
C                The units of P must be consistent with those of RE.
C
C
C     LAT        is a planetodetic latitude value defining a cone.
C                Units are radius.
C
C
C$ Detailed_Output
C
C     The function returns .TRUE. if any of the following
C     are true:
C
C        - The input spheroid is prolate or spherical.
C
C        - The input latitude is zero.
C
C        - The input point is determined to have planetodetic
C          latitude of the same sign as the input latitude.
C
C    Otherwise the function returns .FALSE.
C
C    See Particulars below for details.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If either the equatorial radius or flattening coefficient is
C         invalid, the error SPICE(VALUEOUTOFRANGE) will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports computation of the intersection between
C     a ray and a planetodetic volume element. The routine is used
C     to determine whether the intersection of the ray and latitude
C     cone is actually on the element's boundary.
C
C     This function serves as a "macro" that executes a logical test
C     that would be awkward to perform in-line.
C
C     For a given reference spheroid, all points having a given 
C     planetodetic latitude lie on a cone. The axis of the cone
C     coincides with the Z axis of the coordinate system.
C     
C     The possibility that a point on a given latitude can have a Z
C     coordinate with the wrong sign exists for oblate reference
C     spheroids. For these shapes, the vertex of a positive latitude
C     cone has a negative Z component, and the vertex of a negative
C     latitude cone has a positive Z component.
C
C     The purpose of the function is to determine, for a point that
C     lies on a given planetodetic latitude cone, whether that point is
C     on the correct side of the X-Y plane: that is, the side
C     corresponding to the sign of the input latitude value.
C
C     This check is not as simple as checking the sign of the Z
C     component of the input point. For values of LAT that are non-zero
C     but have very small magnitude, points that are nominally on the
C     corresponding latitude cone may have Z components of the wrong
C     sign due to round-off errors that occurred in the process of
C     computing those points. For such cases, the input point is
C     checked by comparing its distance from the Z axis to the X
C     intercept of a line normal to the reference spheroid at a point
C     having planetodetic latitude LAT and longitude 0. This check can
C     be performed accurately even when the Z component of the input
C     point consists of noise.
C
C     This routine does not test the input point to determine whether
C     it lies on the latitude cone defined by the input argument LAT.
C     The caller must ensure that this is the case.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     This routine does not test the input point to determine whether
C     it lies on the latitude cone defined by the input argument LAT.
C     The caller must ensure that this is the case.
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
C-    SPICELIB Version 1.0.0, 13-JAN-2017 (NJB)
C
C-&
 
C$ Index_Entries
C
C     is point on planetodetic latitude boundary
C
C-&




C
C     SPICELIB functions
C      
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      DOUBLE PRECISION      LIMIT
      PARAMETER           ( LIMIT = 1.D-2 )

C
C     Local variables
C
      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      R
      DOUBLE PRECISION      R2
      DOUBLE PRECISION      XXPT
      DOUBLE PRECISION      YXPT

C
C     Give the function a default value.
C
      ZZPDPLTC = .FALSE.

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZPDPLTC' )


C
C     The equatorial radius must be greater than zero.
C
      IF ( RE .LE. 0.0D0 ) THEN

          CALL SETMSG ( 'Equatorial radius was *.' )
          CALL ERRDP  ( '*', RE                    )
          CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'   )
          CALL CHKOUT ( 'ZZPDPLTC'                 )
          RETURN

      END IF

C
C     If the flattening coefficient is greater than one, the polar
C     radius computed below is negative. If it's equal to one, the
C     polar radius is zero. Either case is a problem, so signal an
C     error and check out.
C
      IF ( F .GE. 1 ) THEN

          CALL SETMSG ( 'Flattening coefficient was *.'  )
          CALL ERRDP  ( '*', F                           )
          CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
          CALL CHKOUT ( 'ZZPDPLTC'                       )
          RETURN

      END IF

C
C     The input point is assumed to be on the cone
C     corresponding to the input latitude. 
C
C     If the reference spheroid is prolate or spherical, there's
C     nothing to do: the point is automatically on the correct side of
C     the X-Y plane.
C
      IF ( F .LE. 0.D0 ) THEN

         ZZPDPLTC = .TRUE.

      ELSE
C
C        This is the oblate case.
C
C
C        If the point is on the "correct" side of the X-Y plane---
C        that is, its Z component as the same sign as LAT, the
C        point is considered to have the correct latitude.
C
C        If the point is on the X-Y plane, or if LAT is zero, the point
C        is considered to have the indicated latitude. We condense
C        these cases by requiring that
C
C              LAT * P(3) >= 0
C
C           rather than
C
C              LAT * P(3) > 0
C
C
         IF (  ( P(3) * LAT )  .GE.  0.D0  ) THEN

            ZZPDPLTC = .TRUE.


         ELSE IF ( ABS(LAT) .GE. LIMIT ) THEN
C
C           Ideally, the input point is considered to have the given
C           latitude if the point is on the side of the X-Y plane
C           corresponding to the sign of the input latitude. The
C           problem with this criterion is that it can't be applied
C           correctly when LAT has very small magnitude.
C
C           If the magnitude of LAT is above the limit, it's ok to 
C           use the sign of P(3) to determine whether the point 
C           has the given latitude.
C
C           The point has the indicated latitude if both LAT and P(3)
C           have the same sign. In this case, we know they have the
C           opposite sign.
C
            ZZPDPLTC = .FALSE.

         ELSE
C
C           At this point we know LAT is non-zero, so the cone
C           corresponding to LAT has its vertex on the opposite side of
C           the X-Y plane from any point having latitude LAT. So it's
C           possible for a point to be on the cone but not have the
C           correct latitude.
C
C           We're in the special case where the point's Z component
C           has the opposite sign as LAT, and the magnitude of LAT
C           is below the limit. We don't automatically reject the
C           point in this case: we'll accept it if it is far enough
C           from the Z axis to be outside the portion of the latitude
C           cone on the wrong side of the X-Y plane.
C
            A = RE
            B = A * ( 1.D0 - F )

C   
C           Compute the intercepts of a normal vector of a point
C           at latitude LAT, longitude 0, with the X and Y axes.
C
            CALL ZZELNAXX ( A, B, LAT, XXPT, YXPT )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZPDPLTC' )
               RETURN
            END IF

C
C           We check the point's distance from the Z axis. This can be
C           done accurately even when the Z component of P consists of
C           noise.
C
C           The point is considered to have the correct latitude when
C           it is farther from the Z axis than the intercept on the X
C           axis of a normal line passing through a point having
C           latitude LAT and longitude 0. Ideally, a point that is on
C           the latitude cone and that satisfies this criterion must be
C           on the correct side of the X-Y plane.
C           

            R2       = ( P(1) * P(1) ) +  ( P(2) * P(2) )

            R        = SQRT(  MAX( R2, 0.D0 )  )

            ZZPDPLTC = R .GE. XXPT

         END IF

      END IF      


      CALL CHKOUT ( 'ZZPDPLTC' )
      RETURN
      END
