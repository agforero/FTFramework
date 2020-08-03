C$Procedure ZZINLAT0 ( DSK, in latitudinal element, w/o margin? )

      SUBROUTINE ZZINLAT0 ( R, LON, LAT, BOUNDS, EXCLUD, INSIDE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of latitudinal coordinates for
C     inclusion in a specified latitudinal volume element. The test is
C     performed without using the "greedy" margins for the element.
C     The built-in angular rounding tolerance is still used for the
C     longitude comparison.
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

      DOUBLE PRECISION      R
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      INTEGER               EXCLUD
      LOGICAL               INSIDE

      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

      INTEGER               LONIDX
      PARAMETER           ( LONIDX = 1 )

      INTEGER               LATIDX
      PARAMETER           ( LATIDX = 2 )

      INTEGER               RADIDX
      PARAMETER           ( RADIDX = 3 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     R          I   Radius of point.   
C     LON        I   Planetocentric longitude of point.
C     LAT        I   Planetocentric latitude of point.
C     BOUNDS     I   Coordinate bounds of element.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C     ANGMRG     P   Angular rounding margin.
C     LONALI     P   Longitude shift margin.
C     NONE       P   Meaning: exclude nothing.
C     LONIDX     P   Longitude index. Also used for exclusion.
C     LATIDX     P   Latitude index. Also used for exclusion.
C     RADIDX     P   Radius index. Also used for exclusion.
C
C$ Detailed_Input
C
C     R,
C     LON,
C     LAT        are, respectively, the radius, longitude, and
C                latitude of a point. The point is to be checked
C                to determine whether it is inside the volume
C                element specified by BOUNDS.
C
C     BOUNDS     is an 2x3 array containing the bounds of a volume
C                element expressed in latitudinal coordinates. BOUNDS
C                defines the volume element used in the comparison. In
C                the element
C
C                   BOUNDS(I,J) 
C
C                J is the coordinate index. J is one of
C
C                   { LONIDX, LATIDX, RADIDX }
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
C     EXCLUD     is either a coordinate index or the parameter NONE.
C
C                If EXCLUD is set to one of
C
C                   { LONIDX, LATIDX, RADIDX }
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
C                coordinate. For example, if a point is on the sphere
C                of radius equal to the upper radius bound, radius need
C                not be used in the comparison and in fact can't be
C                meaningfully compared, due to round-off errors.
C
C$ Detailed_Output
C
C     INSIDE     is a logical flag that is set to .TRUE. if and only if
C                the input coordinates represent a point inside or on
C                the surface of the volume element, according to the
C                comparisons that are performed.
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
C                See dsktol.inc for details.
C
C     ANGMRG     is an angular rounding margin. If the input longitude
C                is outside of the longitude bounds by less than this
C                amount, the longitude is considered to be within the
C                bounds.
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
C     RADIDX     is the index of radius in the second dimension of
C                BOUNDS. When used as a value of the input argument
C                EXCLUD, indicates that radius should be excluded
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
C     Note that a margin is used for longitude shifting, although
C     none is used to expand the volume element.
C
C$ Examples
C
C     See usage in ZZINLAT.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
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
C-    SPICELIB Version 1.0.0, 07-AUG-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     test point against latitudinal element without margin
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      TWOPI

C
C     Local parameters
C     
      INTEGER               LOWER
      PARAMETER           ( LOWER = 1 )
      
      INTEGER               UPPER
      PARAMETER           ( UPPER = 2 )

C
C     Local variables
C
      DOUBLE PRECISION      LOCLON
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      PI2

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



      IF ( FIRST ) THEN

         PI2   = TWOPI()
         FIRST = .FALSE.

      END IF

C
C     Assume the point is outside to start. This allows us
C     to skip setting INSIDE when we find a boundary test
C     failure.
C
      INSIDE = .FALSE.
C 
C     Compare coordinates of the input point and the segment
C     bounds.
C
C     Special case: if the input point is at the origin, and the
C     volume element contains the origin, the point is inside.
C     
      IF ( R .EQ. 0.D0 ) THEN

         IF ( BOUNDS(LOWER, RADIDX) .EQ. 0.D0 ) THEN

            INSIDE = .TRUE.

            RETURN

         END IF

      END IF

         

      IF ( EXCLUD .NE. RADIDX ) THEN
C
C        Compare the point's radius to the segment's radius bounds.
C        
         MINR = BOUNDS(LOWER, RADIDX)
         MAXR = BOUNDS(UPPER, RADIDX)

         IF ( ( R .LT. MINR ) .OR. ( R .GT. MAXR ) ) THEN
C
C           The point's radius is outside of the segment's range.
C
            RETURN

         END IF

      END IF


      IF ( EXCLUD .NE. LATIDX ) THEN
C
C        Compare the point's latitude to the segment's latitude bounds.
C
         MINLAT = BOUNDS(LOWER, LATIDX)
         MAXLAT = BOUNDS(UPPER, LATIDX)

         IF ( ( LAT .LT. MINLAT ) .OR. ( LAT .GT. MAXLAT ) ) THEN
C
C           The point's latitude is outside of the segment's range.
C
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
         CALL ZZNRMLON ( BOUNDS(1,LONIDX),  BOUNDS(2,LONIDX), 
     .                   ANGMRG,   MINLON,  MAXLON           )
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
            RETURN

         END IF


      END IF
C
C     Getting to this point means the input point is inside
C     the segment. Being on the boundary counts as inside.
C  
      INSIDE = .TRUE.
      
      RETURN
      END 
