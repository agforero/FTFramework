C$Procedure ZZELLBDS ( Create bounding ellipsoids )
 
      SUBROUTINE ZZELLBDS ( A, B, HMAX, HMIN, AMAX, BMAX, AMIN, BMIN )
      IMPLICIT NONE
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Given an oblate spheroid and upper and lower height bounds
C     relative to that spheroid, determine radii of inner and outer
C     bounding spheroids.
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
C     DEM
C     DSK
C     ELLIPSOID
C     GEOMETRY
C     TOPOGRAPHY
C
C$ Declarations

      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      HMAX
      DOUBLE PRECISION      HMIN
      DOUBLE PRECISION      AMAX
      DOUBLE PRECISION      BMAX
      DOUBLE PRECISION      AMIN
      DOUBLE PRECISION      BMIN

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     A          I   Input spheroid's semi-major axis length.
C     B          I   Input spheroid's semi-minor axis length.
C     HMAX       I   Maximum height relative to input spheroid.
C     HMIN       I   Minimum height relative to input spheroid.
C     AMAX       O   Outer spheroid's semi-major axis length.
C     BMAX       O   Outer spheroid's semi-minor axis length.
C     AMIN       O   Inner spheroid's semi-major axis length.
C     BMIN       O   Inner spheroid's semi-minor axis length.
C
C$ Detailed_Input
C
C     A          is the semi-major axis length of a reference
C                spheroid. A may have any units; B, HMAX, and
C                HMIN must have the same units.
C
C     B          is the semi-minor axis length of a reference
C                spheroid. B must not exceed A.
C
C     HMAX       is an upper bound for a set of heights relative
C                the reference spheroid defined by A and B. HMAX
C                is a signed quantity.
C
C     HMIN       is a lower bound for a set of heights relative
C                the reference spheroid defined by A and B. HMIN
C                is a signed quantity. HMIN must not exceed HMAX.
C                HMIN must be greater than -B.
C
C$ Detailed_Output
C
C     AMAX,
C     BMAX       are, respectively, semi-major and semi-minor axis
C                lengths for an outer bounding spheroid. The set of
C                points at height HMAX relative to the input spheroid
C                defined by A and B lies on or below the spheroid
C                defined by AMAX and BMAX.
C
C                When HMAX is non-negative
C
C                   AMAX  =  A + HMAX
C                   BMAX  =  B + HMAX*(A/B)
C
C                When HMAX is negative
C                  
C                   AMAX  =  A + HMAX*(B/A)
C                   BMAX  =  B + HMAX
C     AMIN,
C     BMIN       are, respectively, semi-major and semi-minor axis
C                lengths for an inner bounding spheroid. The set of
C                points at height HMIN relative to the input spheroid
C                defined by A and B lies on or above the spheroid
C                defined by AMIN and BMIN.
C
C                When HMIN is non-positive
C
C                   AMIN  =  A + HMIN
C                   BMIN  =  B + HMIN*(A/B)
C
C                When HMIN is positive
C                  
C                   AMIN  =  A + HMIN*(B/A)
C                   BMIN  =  B + HMIN
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the semi-minor axis length B is non-positive, the
C         error SPICE(NONPOSITIVERADIUS) will be signaled.
C
C     2)  If the semi-major axis length A is less than the 
C         semi-minor axis length B, the error SPICE(RADIIOUTOFORDER)
C         will be signaled.
C
C     3)  If HMIN is less than or equal to -B, the error
C         SPICE(LOWERBOUNDTOOLOW) will be signaled.
C
C     4)  If the lower height bound HMIN is greater than the 
C         upper height bound HMAX, the error SPICE(BOUNDSOUTOFORDER)
C         will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Ellipsoidal bounding surfaces may be used for volume elements in
C     planetodetic coordinates. The combination of these surfaces,
C     together with cones of constant latitude and planes of constant
C     longitude, may be used to create a simple set of bounding
C     surfaces for such an element. These surfaces can be used to
C     rapidly determine that a given ray does not intersect the volume
C     element. 
C
C     For a surface represented by a digital elevation model (DEM),
C     efficient solution of the ray-surface intercept problem is
C     enhanced by determination of a line segment---also called a
C     "chord"---outside of which no solution is possible. One way of
C     generating a chord is to determine the intersection of the ray
C     with inner or outer bounding surfaces. An inner bounding surface
C     is one which has, for every latitude and longitude covered by the
C     DEM, altitude less than or equal to that of the DEM. An outer
C     bounding surface has, for every latitude and longitude covered by
C     the DEM, altitude greater than or equal to the DEM.
C
C     In order for bounding surfaces to be useful, it must be possible
C     to rapidly compute intersections of rays with these surfaces. 
C     The bounding surfaces must also lie "close" to the DEM surface;
C     otherwise chords formed by intersections of rays with the bounding
C     surfaces may be longer than necessary.
C
C     Spheroids are a natural choice for bounding surfaces, since DEMs
C     are typically referenced to a spheroidal surface, and because
C     ray-ellipsoid intersections can be rapidly computed by closed-form
C     algebraic methods.
C
C     It might seem that, given a spheroid with semi-major axis length
C     A and semi-minor axis length B, and given a surface having maximum
C     and minimum spheroid-relative heights HMAX and HMIN respectively,
C     that the level surfaces having heights HMAX and HMIN relative to
C     the spheroid would be candidates for the inner and outer bounding
C     spheroids. The outer spheroid would have radii A+HMAX and
C     B+HMAX, while the inner spheroid would have radii A+HMIN and
C     B+HMIN. However, it can be shown (by numerical comparison, for
C     example), that these spheroids do not bound the level surfaces
C     at heights HMAX and HMIN, except in certain special cases such as
C     that of a spherical reference surface.
C
C     This routine generates semi-axis lengths of inner and outer
C     bounding spheroids that are valid for all eccentricities of the
C     reference spheroid (as long as A >= B) and, for reference
C     spheroids of low eccentricity, lie close to the level surfaces at
C     heights HMIN and HMAX. A discussion of the method follows.
C
C
C     Explanation
C     ===========
C
C     Since we're working with spheroids, we can reduce the problem
C     to a two-dimensional one. We'll compute inner and outer bounding
C     ellipses for an ellipse with positive semi-axis lengths A and B.
C     Revolving the bounding ellipses about the Z axis creates bounding
C     ellipsoids.
C
C        Consider an ellipse E with vertical semi-minor axis length
C        B and horizontal semi-major axis A. For a point (x,y) on E,
C        let N be the normal vector
C
C                  x      y
C           N = ( ---- , ---- )
C                   2      2
C                  A      B
C
C        Let LAMBDA be a constant, and Let E' be the curve
C
C           { (x,y) + LAMBDA*N }
C
C        Then E' is an ellipse (since it's produced by a linear
C        transformation of E) with semi-axis lengths
C
C           (A + LAMBDA/A), (B + LAMBDA/B)
C
C
C        For any point (x,y) on E, the height of E' above that point is
C        
C           LAMBDA*||N||
C
C        The square of this height HSQ is 
C
C
C                       2
C           HSQ = LAMBDA  * <N, N>
C
C               
C                               2      2
C                       2      x      y
C               = LAMBDA  * ( ---- + ---- )
C                               4      4
C                              A      B
C
C                               2     2       2  2
C                       2      x     B ( 1 - x /A )
C               = LAMBDA  * ( ---- + -------------- )
C                               4           4
C                              A           B
C
C                               2                   
C                       2      x     1        1        1
C               = LAMBDA  * ( ---- (----  -  ----) +  ---- )
C                               2     2        2        2
C                              A     A        B        B
C
C
C        If A = B, this expression is constant.
C
C        If A > B, then the term
C
C           1/A**2 - 1/B**2         
C
C        is negative, so HSQ is a decreasing function of x for 
C
C           0 <= x <= A
C
C        This implies that for x in the above range, E' is closest
C        to E when x = A. We'll use this fact to generate bounding
C        ellipsoids with the following properties:
C
C           If HMAX >= 0    E' will have height HMAX at x = A and
C                           height >= HMAX if 0 <= x < A.
C           
C           If HMAX <  0    E' will have height HMAX at x = 0 and
C                           height >= HMAX if 0 < x <= A.
C
C           If HMIN <= 0    E' will have height HMIN at x = A and
C                           height <= HMIN if 0 <= x < A.
C
C           If HMIN >  0    E' will have height HMIN at x = 0 and
C                           height <= HMIN if 0 < x <= A.
C
C
C
C     Application to prolate spheroids
C     ================================
C
C     For a spheroid having semi-axes A, B, C, for which 
C
C        A = B < C
C
C     This routine can be applied to the semi-axis lengths
C     A' and B', where
C         
C        A'  =  C
C        B'  =  A
C
C     The outer bounding surface is a prolate spheroid with semi-axis
C     lengths
C
C        BMAX, BMAX, AMAX
C
C     The inner bounding surface is a prolate spheroid with semi-axis
C     lengths
C
C        BMIN, BMIN, AMIN
C        
C
C$ Examples
C
C     See usage in ZZRYXPDT.
C
C$ Restrictions
C
C     1)  The bounding spheroids generated by this routine may not be
C         suitable for reference spheroids with high eccentricity.
C
C     2)  Callers of this routine should not rely explicitly on the
C         formulas shown in the Detailed Output section of this
C         routine's header comments. The formulas could be upgraded
C         in a future version of this routine.       
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
C-    SPICELIB Version 1.0.0, 13-JAN-2017 (NJB)
C
C        Added error check for oblate case where lower bounding
C        ellipsoid has negative polar radius.
C
C        Updated $Particulars; added mention of use for volume elements
C        in planetodetic coordinate systems. Remarks about application
C        to prolate spheroids were added. Changed contents of $Files
C        section to "None." Made miscellaneous small changes to
C        comments.
C
C     Original version 27-NOV-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     create bounding ellipsoids for a dem
C     bounding ellipsoids for a digital elevation model
C     create bounding ellipsoids for a dsk type 4 segment
C
C-&

 
C
C
C     Use discovery check-in.
C
      IF ( B .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZELLBDS'                                )
         CALL SETMSG ( 'This routine requires B > 0, but B = #.' )
         CALL ERRDP  ( '#',  B                                   )
         CALL SIGERR ( 'SPICE(NONPOSITIVERADIUS)'                )
         CALL CHKOUT ( 'ZZELLBDS'                                )
         RETURN

      END IF

      IF ( B .GT. A ) THEN

         CALL CHKIN  ( 'ZZELLBDS'                                 )
         CALL SETMSG ( 'This routine requires A >= B, but A = #; '
     .   //            'B = #.'                                   )
         CALL ERRDP  ( '#',  A                                    )
         CALL ERRDP  ( '#',  B                                    )
         CALL SIGERR ( 'SPICE(RADIIOUTOFORDER)'                   )
         CALL CHKOUT ( 'ZZELLBDS'                                 )
         RETURN

      END IF

      IF ( ( B + HMIN ) .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZELLBDS'                                 )
         CALL SETMSG ( 'This routine requires B + HMIN > 0, but '
     .   //            'B = #; HMIN = #, B+HMIN = #.'             )
         CALL ERRDP  ( '#',  B                                    )
         CALL ERRDP  ( '#',  HMIN                                 )
         CALL ERRDP  ( '#',  B + HMIN                             )
         CALL SIGERR ( 'SPICE(LOWERBOUNDTOOLOW)'                  )
         CALL CHKOUT ( 'ZZELLBDS'                                 )
         RETURN

      END IF

      IF ( HMIN .LT. 0.D0 ) THEN

         IF (  ( B + (A/B)*HMIN ) .LE. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZELLBDS'                                 )
            CALL SETMSG ( 'For oblate spheroids and HMIN < 0, '
     .      //            'This routine requires B + (A/B)HMIN > 0, '
     .      //            'but A = #, B = #; HMIN = #, '
     .      //            'B+(A/B)HMIN = #.'                         )
            CALL ERRDP  ( '#',  A                                    )
            CALL ERRDP  ( '#',  B                                    )
            CALL ERRDP  ( '#',  HMIN                                 )
            CALL ERRDP  ( '#',  B + (A/B)*HMIN                       )
            CALL SIGERR ( 'SPICE(LOWERBOUNDTOOLOW)'                  )
            CALL CHKOUT ( 'ZZELLBDS'                                 )
            RETURN

         END IF

      END IF


      IF ( HMIN .GT. HMAX ) THEN

         CALL CHKIN  ( 'ZZELLBDS'                                 )
         CALL SETMSG ( 'This routine requires HMAX >= HMIN, but '
     .   //            'HMIN = #; HMAX = #.'                      )
         CALL ERRDP  ( '#',  HMIN                                 )
         CALL ERRDP  ( '#',  HMAX                                 )
         CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)'                  )
         CALL CHKOUT ( 'ZZELLBDS'                                 )
         RETURN

      END IF

C
C     In the following comments, N, E, E', and LAMBDA are
C     defined as in the Particulars section above.
C
C           
C     Generate radii of the outer bounding ellipsoid.
C
      IF ( HMAX .GE. 0.D0 ) THEN
C
C        Pick radii of E' so that E' matches
C
C           E + HMAX * N / ||N||
C
C        that is, E' has height HMAX above E, at x=A. 
C
C        For smaller x, the height of E' above E will
C        will be greater than or equal to HMAX.
C
C        Set LAMBDA = A * HMAX.
C
C        Then the radii of E' are
C
C                             |
C           x + LAMBDA*||N||  |
C                             |x=A,y=0
C
C        and
C                             |
C           y + LAMBDA*||N||  |
C                             |x=0,y=B
C
C        so the radii of E', AMAX and BMAX, are:
C
C           AMAX =  A + LAMBDA*A/A**2  =  A + HMAX
C           BMAX =  B + LAMBDA*B/B**2  =  B + HMAX*(A/B)
C
C
         AMAX = A  +  HMAX
         BMAX = B  +  HMAX * (A/B)

      ELSE
C
C        HMAX < 0.
C
C        In this case the outer bounding ellipse should match E+HMAX
C        at x = 0. The ellipse will be closer to E for x > 0.
C
C        Set LAMBDA = B * HMAX. Then 
C
C           AMAX =  A + LAMBDA*A/A**2  =  A + HMAX * (B/A)
C           BMAX =  B + LAMBDA*B/B**2  =  B + HMAX

         AMAX = A  +  HMAX * (B/A)
         BMAX = B  +  HMAX

      END IF
      
C
C     Find radii of the inner bounding ellipsoid.
C
      IF ( HMIN .LE. 0 ) THEN
C
C        This case is similar to that of the outer bounding
C        ellipsoid for HMAX >= 0. We can create an ellipse
C        that has height HMIN at x = A and that is further 
C        from E for x < A.        
C
C        Set LAMBDA = A * HMAX. Then 
C
C           AMAX =  A + LAMBDA*A/A**2  =  A + HMAX
C           BMAX =  B + LAMBDA*B/B**2  =  B + HMAX*(A/B)

         AMIN = A  +  HMIN
         BMIN = B  +  HMIN * (A/B)

      ELSE 
C
C        HMIN > 0.
C
C        In this case the inner bounding ellipse should match E+HMIN
C        at x = 0. The ellipse will be closer to E for x > 0.
C
C        Set LAMBDA = B * HMAX. Then 
C
C           AMIN =  A + LAMBDA*A/A**2  =  A + HMIN * (B/A)
C           BMIN =  B + LAMBDA*B/B**2  =  B + HMIN
C
         AMIN = A  +  HMIN * (B/A)
         BMIN = B  +  HMIN

      END IF

      END

