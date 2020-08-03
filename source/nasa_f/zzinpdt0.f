C$Procedure ZZINPDT0 ( DSK, in planetodetic element, w/o margin? )

      SUBROUTINE ZZINPDT0 ( P, LON, BOUNDS, CORPAR, EXCLUD, INSIDE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of planetodetic coordinates for
C     inclusion in a specified planetodetic volume element. Ellipsoidal
C     surfaces are used in place of surfaces of constant altitude. The
C     test is performed without using the "greedy" margins for the
C     element. The built-in angular rounding tolerance is still used
C     for the longitude comparison.
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
C     DSK
C
C$ Keywords
C     
C     DSK
C     GEOMETRY
C     INTERSECTION
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsktol.inc'

      DOUBLE PRECISION      P      ( 3 )
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      CORPAR ( * )
      INTEGER               EXCLUD
      LOGICAL               INSIDE


      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

      INTEGER               LONIDX
      PARAMETER           ( LONIDX = 1 )

      INTEGER               LATIDX
      PARAMETER           ( LATIDX = 2 )

      INTEGER               ALTIDX
      PARAMETER           ( ALTIDX = 3 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   Cartesian coordinates of point.
C     LON        I   Planetocentric longitude of point.
C     BOUNDS     I   Coordinate bounds of element.
C     CORPAR     I   Coordinate parameters.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C     ANGMRG     P   Angular rounding margin.
C     LONALI     P   Longitude alias margin.
C     NONE       P   Meaning: exclude nothing.
C     LONIDX     P   Longitude index. Also used for exclusion.
C     LATIDX     P   Latitude index. Also used for exclusion.
C     ALTIDX     P   Altitude index. Also used for exclusion.
C
C$ Detailed_Input
C
C     P          is a point expressed in Cartesian coordinates. The
C                reference frame is that of the segment of interest.
C                The point is to be checked to determine whether it is
C                inside the volume element specified by BOUNDS.
C
C
C     LON        is the planetocentric longitude of the point P.
C
C
C     BOUNDS     is an 2x3 array containing the bounds of a volume
C                element expressed in planetodetic coordinates. BOUNDS
C                defines the volume element used in the comparison. In
C                the element
C
C                   BOUNDS(I,J) 
C
C                J is the coordinate index. J is one of
C
C                   { LONIDX, LATIDX, ALTIDX }
C
C                I is the bound index.
C
C                   I = 1   ->   lower bound
C                   I = 2   ->   upper bound
C
C                If necessary, the longitude bounds are normalized
C                for the comparison. The normalized bounds differ
C                by no more than 2*pi and are in increasing order.
C                See ZZNRMLON for details.
C
C                The aliasing margin LONALI and the angular rounding
C                margin ANGMRG are used for longitude comparisons.
C                However, the greedy margin is not used. See the
C                Parameters section below for further information.
C                
C
C     CORPAR     is an array of coordinate parameters; normally these
C                are obtained from a DSK descriptor. The first element
C                is the equatorial radius of the reference spheroid.
C                The second element is the flattening coefficient.
C                Any other elements are ignored by this routine.
C
C
C     EXCLUD     is either a coordinate index or the parameter NONE.
C
C                If EXCLUD is set to one of
C
C                   { LONIDX, LATIDX, ALTIDX }
C
C                then the indicated coordinate is excluded from
C                comparison with the corresponding volume element
C                boundaries.
C
C                If EXCLUD is set to NONE, all coordinates are
C                compared.
C
C                Exclusion of coordinates is used in cases where a
C                point is known to be on a level surface of a given
C                coordinate. For example, if a point is on the surface
C                of altitude equal to the upper altitude bound,
C                altitude need not be used in the comparison and in
C                fact can't be meaningfully compared, due to round-off
C                errors.
C
C$ Detailed_Output
C
C     INSIDE     is a logical flag that is set to .TRUE. if and
C                only if the input coordinates represent a 
C                point inside or on the surface of the volume
C                element, according to the comparisons that are
C                performed.
C
C                The value of INSIDE is not affected by the value
C                of any excluded coordinate.
C
C$ Parameters
C
C     LONALI     is a longitude shift margin. If the input longitude
C                coordinate LON is less than the longitude lower bound
C                by more than this amount, the value LON + 2*pi is used
C                in the longitude comparison. If the input longitude
C                exceeds the upper longitude bound by this amount, the
C                value LON - 2*pi is used in the comparison.
C
C     ANGMRG     is an angular rounding margin. If the input longitude
C                is outside of the longitude bounds by less than this
C                amount, the longitude is considered to be within the
C                bounds. The margin is used analogously for latitude.
C
C                This margin is distinct from the "greedy" segment
C                selection margin, which is used to expand a volume
C                element for the purpose of determining whether a given
C                point lies in that element.
C
C                See dsktol.inc for details.
C
C     NONE       when used as a value of the input argument EXCLUD,
C                indicates that no coordinates should be excluded from
C                comparison.
C
C     LONIDX     is the index of longitude in the second dimension of
C                BOUNDS. When used as a value of the input argument
C                EXCLUD, indicates that longitude should be excluded
C                from comparison.
C
C     LATIDX     is the index of latitude in the second dimension of
C                BOUNDS. When used as a value of the input argument
C                EXCLUD, indicates that latitude should be excluded
C                from comparison.
C
C
C     ALTIDX     is the index of altitude in the second dimension of
C                BOUNDS. When used as a value of the input argument
C                EXCLUD, indicates that altitude should be excluded
C                from comparison.
C
C$ Exceptions
C
C     Error free.
C
C$ Files
C
C     None.
C     
C$ Particulars
C
C     This routine uses ellipsoidal bounding surfaces, in place of 
C     surfaces of constant altitude, to represent the "outer" and
C     "inner" segment boundaries. These surfaces must be consistent
C     with those used by the ray-planetodetic volume element intercept
C     routine ZZRYTPDT.
C
C     The radii of these ellipsoidal surfaces are those produced by
C     the routine ZZELLBDS.
C     
C     Note that a margin is used for longitude shifting, although
C     none is used to expand the volume element.
C
C$ Examples
C
C     See usage in ZZINPDT.
C
C$ Restrictions
C      
C     1)  Validity of the ellipsoid shape parameters must be ensured by
C         the caller.
C
C     2)  This is a private routine. It is meant to be used only by the
C         DSK subsystem.
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
C-    SPICELIB Version 1.0.0, 16-AUG-2016 (NJB) 
C     
C        Now uses ANGMRG for latitude comparisons.
C
C        25-MAY-2016 (NJB) 
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     test point against planetodetic element without margin
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      TWOPI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C     
      INTEGER               LOWER
      PARAMETER           ( LOWER = 1 )
      
      INTEGER               UPPER
      PARAMETER           ( UPPER = 2 )

C
C     Numeric relation codes returned by ZZPDCMPL:
C     
      INTEGER               LT
      PARAMETER           ( LT = -1 )
      
      INTEGER               EQ
      PARAMETER           ( EQ =  0 )

      INTEGER               GT
      PARAMETER           ( GT =  1 )


C
C     Local variables
C
      DOUBLE PRECISION      EMAX
      DOUBLE PRECISION      EMIN
      DOUBLE PRECISION      F
      DOUBLE PRECISION      LEVEL
      DOUBLE PRECISION      LOCLON
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXALT
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINALT
      DOUBLE PRECISION      PI2
      DOUBLE PRECISION      PMAX
      DOUBLE PRECISION      PMIN
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RP

      INTEGER               MAXREL
      INTEGER               MINREL

      LOGICAL               FIRST

      
      
C
C     Saved variables
C      
      SAVE                  FIRST
      SAVE                  PI2

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /
      DATA                  PI2   / -1.D0  /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZINPDT0' )

      IF ( FIRST ) THEN

         PI2   = TWOPI()
         FIRST = .FALSE.

      END IF

C
C     Unpack the shape parameters. Set the polar axis length for
C     later use.
C
      RE = CORPAR(1)
      F  = CORPAR(2)
      RP = RE * ( 1.D0 - F )
      
C
C     Assume the point is outside to start. This allows us
C     to skip setting INSIDE when we find a boundary test
C     failure.
C
      INSIDE = .FALSE.
C 
C     Compare coordinates of the input point and the segment
C     bounds. Deal with altitude last, since it involves the
C     most work. We may be able to exit before performing the 
C     altitude tests.
C
      IF ( EXCLUD .NE. LATIDX ) THEN
C
C        Compare the point's latitude to the segment's latitude bounds.
C
         MINLAT = MAX( BOUNDS(LOWER, LATIDX) - ANGMRG, -HALFPI() )
         MAXLAT = MIN( BOUNDS(UPPER, LATIDX) + ANGMRG,  HALFPI() )

         CALL ZZPDCMPL ( RE, F, P, MINLAT, MINREL )
         CALL ZZPDCMPL ( RE, F, P, MAXLAT, MAXREL )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZINPDT0' )
            RETURN
         END IF
 

         IF (      ( MINREL .EQ. LT ) 
     .        .OR. ( MAXREL .EQ. GT )  ) THEN
C
C           The point's latitude is outside of the segment's range.
C            
            CALL CHKOUT ( 'ZZINPDT0' )
            RETURN

         END IF

      END IF

C
C     Move the longitude of the input point into the interval 
C        
C        [ MINLON, MAXLON ]
C
C     if necessary and if possible.
C
      IF ( EXCLUD .NE. LONIDX ) THEN
C
C        Put the local longitude bounds in order, if necessary.
C
         CALL ZZNRMLON ( BOUNDS(LOWER,LONIDX),  BOUNDS(UPPER,LONIDX), 
     .                   ANGMRG,       MINLON,  MAXLON               )
C
C        Compare the point's longitude to the segment's longitude
C        bounds.
C      
         LOCLON = LON

         IF ( LON .LT. ( MINLON - LONALI ) ) THEN
C
C           If the point's longitude is less than the segment's
C           longitude by more than a small margin, shift the longitude
C           right by 2*pi.

            LOCLON = LON + PI2         

         ELSE IF ( LON .GT. ( MAXLON + LONALI ) ) THEN
C
C           If the point's longitude is greater than the segment's
C           longitude by more than a small margin, shift the longitude
C           left by 2*pi.

            LOCLON = LON - PI2         

         END IF

         
         IF (        ( LOCLON .LT. (MINLON-ANGMRG) )
     .          .OR. ( LOCLON .GT. (MAXLON+ANGMRG) )  ) THEN
C
C           The point's longitude, adjusted if necessary for
C           comparison, is outside of the segment's range.
C
            CALL CHKOUT ( 'ZZINPDT0' )
            RETURN

         END IF


      END IF

      IF ( EXCLUD .NE. ALTIDX ) THEN
C
C        Extract altitude bounds from the segment descriptor.
C
         MINALT = BOUNDS(LOWER, ALTIDX)
         MAXALT = BOUNDS(UPPER, ALTIDX)
C
C        Find the semi-axes of the bounding spheroids.
C        
         IF ( F .GE. 0.D0 ) THEN
C
C           This is the oblate case. Get the semi-axis lengths
C           for the inner and outer bounding spheroids.
C
            CALL ZZELLBDS ( RE,   RP,   MAXALT, MINALT, 
     .                      EMAX, PMAX, EMIN,   PMIN   )

         ELSE
C
C           This is the prolate case. RP is the longer axis. Get the
C           semi-axis lengths for the inner and outer bounding
C           spheroids.
C 
C           In this call, we'll store the radii associated with
C           the longer axis in PMAX and PMIN.
C
            CALL ZZELLBDS ( RP,   RE,   MAXALT, MINALT, 
     .                      PMAX, EMAX, PMIN,   EMIN   )

         END IF

C
C        Compute the input point's level surface parameters
C        for the inner and outer bounding ellipsoids. Do these
C        computations one at a time, so the second one can be
C        skipped if the first one shows the point is outside
C        the outer ellipsoid.
C        
         LEVEL = ( P(1)/EMAX )**2 + ( P(2)/EMAX )**2 + ( P(3)/PMAX )**2 

         IF ( LEVEL .GT. 1.D0 ) THEN
C
C           The point is outside of the outer bounding ellipsoid.
C
            CALL CHKOUT ( 'ZZINPDT0' )
            RETURN

         END IF


         LEVEL = ( P(1)/EMIN )**2 + ( P(2)/EMIN )**2 + ( P(3)/PMIN )**2 

         IF ( LEVEL .LT. 1.D0 ) THEN
C
C           The point is inside the inner bounding ellipsoid, which
C           implies it is outside of the segment's boundaries.
C
            CALL CHKOUT ( 'ZZINPDT0' )
            RETURN

         END IF

      END IF

C
C     Getting to this point means the input point is inside
C     the segment. Being on the boundary counts as inside.
C  
      INSIDE = .TRUE.

    
      CALL CHKOUT ( 'ZZINPDT0' )  
      RETURN
      END 
