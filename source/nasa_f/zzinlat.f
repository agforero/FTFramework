C$Procedure ZZINLAT ( DSK, in latitudinal element? )

      SUBROUTINE ZZINLAT ( P, BOUNDS, MARGIN, EXCLUD, INSIDE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of latitudinal coordinates for
C     inclusion in a specified latitudinal volume element. The test is
C     performed using margins for the element.
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
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      MARGIN
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

      DOUBLE PRECISION      LATMRG
      PARAMETER           ( LATMRG = 1.D-8 )


C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   Input point.
C     BOUNDS     I   Coordinate bounds of element.
C     MARGIN     I   Margin used for inclusion testing.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C     LATMRG     P   Margin for latitude near the poles.
C     LONALI     P   Longitude shift margin.
C     ANGMRG     P   Angular rounding margin.
C     NONE       P   Meaning: exclude nothing.
C     LONIDX     P   Longitude index. Also used for exclusion.
C     LATIDX     P   Latitude index. Also used for exclusion.
C     RADIDX     P   Radius index. Also used for exclusion.
C
C$ Detailed_Input
C
C     P          is a point expressed in Cartesian coordinates. The
C                point is to be checked to determine whether it is
C                inside the volume element specified by BOUNDS.
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
C                The aliasing margin LONALI, the greedy margin, and the
C                angular rounding margin ANGMRG are used for longitude
C                comparisons. ANGMRG is used if its magnitude exceeds
C                that of the greedy margin. See the Parameters section
C                below for further information.
C                
C
C     MARGIN     is a fraction used to expand the volume element for
C                inclusion testing. 
C
C                   - The outer radius is scaled by (1+MARGIN).
C                     The inner radius is scaled by (1-MARGIN).
C
C                   - Latitude bounds are extended by MARGIN radians,
C                     within the interval -pi/2 : pi/2.
C
C                   - For any input point having latitude "not close"
C                     to pi/2 or -pi/2, the element's longitude
C                     interval is extended by
C
C                        MARGIN / cos(LAT)
C
C                     where LAT is the latitude of the input point.
C         
C                     Here "close to" means "within LATMRG radians of."
C
C                   - For any input point having latitude close to 
C                     pi/2 or -pi/2, comparison against the longitude
C                     bounds of the volume element is omitted if
C                     MARGIN > 0.
C   
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
C
C     LATMRG     Margin for determining whether the input latitude is
C                "close" to the poles. Units are radians. This value
C                should be kept in sync with the default "greedy"
C                margin.
C
C     LONALI     is a longitude shift margin. If the input longitude 
C                coordinate LON is less than the longitude lower bound
C                by more than this amount, the value LON + 2*pi is
C                used in the longitude comparison, provided that the
C                comparison is not excluded.
C
C                See dsktol.inc for details.
C
C     ANGMRG     is an angular rounding margin. If the input longitude
C                is outside of the longitude bounds by less than this
C                amount, the longitude is considered to be within the
C                bounds. ANGMRG is used when the magnitude of the
C                "greedy" margin is less than ANGMRG. ANGMRG is not
C                used for computations involving coordinates other than
C                longitude.
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
C     1)  If MARGIN is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     2)  If EXCLUD is outside of the range 0:3, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C     
C$ Particulars
C
C     None.
C
C$ Examples
C
C     See usage in ZZRYTLAT.
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
C-    SPICELIB Version 1.0.0, 03-JUN-2016 (NJB) 
C
C        Original version 03-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     test point against latitudinal element using margin
C
C-&




C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      TWOPI

      LOGICAL               RETURN
      
C
C     Local parameters 
C

C
C     Element boundary indices:
C     
      INTEGER               WEST
      PARAMETER           ( WEST   = 1 )

      INTEGER               EAST
      PARAMETER           ( EAST   = 2 )

      INTEGER               SOUTH
      PARAMETER           ( SOUTH  = 1 )

      INTEGER               NORTH
      PARAMETER           ( NORTH  = 2 )

      INTEGER               INNER
      PARAMETER           ( INNER  = 1 )

      INTEGER               OUTER
      PARAMETER           ( OUTER  = 2 )


C
C     Local variables
C
      DOUBLE PRECISION      AMAXR
      DOUBLE PRECISION      AMINR
      DOUBLE PRECISION      AMAXLT
      DOUBLE PRECISION      AMINLT
      DOUBLE PRECISION      AMAXLO
      DOUBLE PRECISION      AMINLO
      DOUBLE PRECISION      DLON
      DOUBLE PRECISION      HPI
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      LONMRG
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      PI2
      DOUBLE PRECISION      R
      DOUBLE PRECISION      SMAX
      DOUBLE PRECISION      SMIN

      LOGICAL               FIRST

      
      SAVE                  FIRST
      SAVE                  PI2
      SAVE                  HPI

      DATA                  FIRST / .TRUE. /
      DATA                  PI2   / -1.D0  /
      DATA                  HPI   / -1.D0  /

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF ( FIRST ) THEN

         PI2   = TWOPI()
         HPI   = HALFPI()

         FIRST = .FALSE.

      END IF

C
C     Get the latitudinal coordinates of the input point.
C
      CALL RECLAT ( P, R, LON, LAT )


C
C     Handle the simpler zero-margin case separately.
C
      IF ( MARGIN .EQ. 0.D0 ) THEN

         CALL ZZINLAT0 ( R, LON, LAT, BOUNDS, EXCLUD, INSIDE )

         RETURN

      ELSE IF ( MARGIN .LT. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZINLAT'                                )
         CALL SETMSG ( 'Margin must be non-negative but was #.' )
         CALL ERRDP  ( '#', MARGIN                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINLAT'                                )
         RETURN

      END IF


      IF (  ( EXCLUD .LT. 0 ) .OR. ( EXCLUD .GT. 3 )  ) THEN

         CALL CHKIN  ( 'ZZINLAT'                                    )
         CALL SETMSG ( 'EXCLUD must be in the range 0:3 but was #.' )
         CALL ERRINT ( '#', EXCLUD                                  )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                     )
         CALL CHKOUT ( 'ZZINLAT'                                    )
         RETURN

      END IF

C
C     Special case: if the input point is within distance MARGIN
C     from the origin, and the minimum radius of the volume element
C     is less than or equal to MARGIN, the point is inside.
C     
      IF ( R .LE. MARGIN ) THEN

         IF ( BOUNDS(INNER, RADIDX) .LE. MARGIN ) THEN

            INSIDE = .TRUE.

            RETURN

         END IF

      END IF


C
C     Assume the point is outside to start. This allows us
C     to skip setting INSIDE when we find a boundary test
C     failure.
C
      INSIDE = .FALSE.

C
C     Get local copies of the coordinate bounds. Don't normalize the
C     longitude bounds until we know we need them.
C 
      MINR   = BOUNDS( INNER, RADIDX ) 
      MAXR   = BOUNDS( OUTER, RADIDX ) 

      MINLAT = BOUNDS( SOUTH, LATIDX )
      MAXLAT = BOUNDS( NORTH, LATIDX )

C
C     Compare coordinates to adjusted latitude and radius
C     boundaries.
C
      IF ( EXCLUD .NE. RADIDX ) THEN
C
C        Create adjusted radius bounds.
C 
         SMAX  = 1.D0 + MARGIN
         SMIN  = 1.D0 - MARGIN

         AMINR = MAX( 0.D0,  MINR*SMIN )
         AMAXR = MAXR * SMAX

         IF ( ( R .LT. AMINR ) .OR. ( R .GT. AMAXR ) ) THEN
             RETURN
         END IF

      END IF


      IF ( EXCLUD .NE. LATIDX ) THEN
C
C        Create adjusted latitude bounds.
C
         AMINLT = MAX( -HPI, MINLAT - MARGIN )
         AMAXLT = MIN(  HPI, MAXLAT + MARGIN )

         IF ( ( LAT .LT. AMINLT ) .OR. ( LAT .GT. AMAXLT ) ) THEN
            RETURN
         END IF

      END IF

C
C     At this point, the input radius and latitude are within the
C     adjusted bounds, if their tests haven't been excluded by
C     the caller.
C
C     Perform longitude tests, unless they're excluded by the 
C     caller.
C
      IF ( EXCLUD .NE. LONIDX ) THEN   

         CALL ZZNRMLON ( BOUNDS(WEST,LONIDX), BOUNDS(EAST,LONIDX),
     .                   ANGMRG,              MINLON,  
     .                   MAXLON                                   )
C
C        Set the margin to be used for longitude interval 
C        inclusion tests.
C
         LONMRG = MAX( ABS(ANGMRG), ABS(MARGIN) )
         
C
C        We have a special case for segments that include the poles. If
C        the input point is close enough to a pole contained in the
C        segment, we consider the point to be included in the segment,
C        regardless of the point's longitude. All other points get the
C        normal longitude test.
C
         IF (       ( LAT .LE. (  HPI - LATMRG ) )
     .        .AND. ( LAT .GE. ( -HPI + LATMRG ) )  ) THEN
C
C           This is the usual case: the latitude of the input point
C           is bounded away from the poles.
C
C           Check the point's longitude against the segment's
C           longitude bounds.
C
C           We'll scale the longitude margin to compensate for the
C           latitude of the input point. Note that the division
C           below is safe; presuming a reasonable value of MARGIN;
C           we know that 
C
C              DLON << 1     
C         
            DLON   = LONMRG / MAX( ABS(COS(LAT)), LATMRG )

            AMINLO = MINLON - DLON
            AMAXLO = MAXLON + DLON

C
C           Now move the input point's longitude into range, if
C           necessary.
C     
            IF ( LON .LT. AMINLO ) THEN

               IF (  LON  .LT. ( AMINLO - LONALI )  ) THEN
C
C                 See whether an aliased version of LON is a match.
C
                  LON = LON + PI2

               ELSE
C
C                 Consider LON to be a match with the lower bound.
C
                  LON = AMINLO

               END IF

            ELSE IF ( LON .GT. AMAXLO ) THEN

               IF (  LON  .GT. ( AMAXLO + LONALI )  ) THEN
C
C                 See whether an aliased version of LON is a match.
C
                  LON = LON - PI2

               ELSE
C
C                 Consider LON to be a match with the upper bound.
C
                  LON = AMAXLO

               END IF

            END IF

C
C           Compare the adjusted longitude of the input point to the
C           adjusted longitude bounds.
C
            IF ( ( LON .LT. AMINLO ) .OR. ( LON .GT. AMAXLO ) ) THEN

               RETURN

            END IF              


         ELSE
C
C           The latitude of the input point is close to one of the
C           poles.
C
C           This is a no-op case.
C           
C           The input point has already passed whichever of the radius
C           and latitude tests that were not excluded.
C
C           If the element has a non-degenerate latitude boundary
C           having the same sign as the latitude of the input point,
C           and if latitude is excluded because the input point is
C           already nominally on that boundary, then passing the radius
C           check implies that the point is close to the element.
C
C           If the element has a degenerate latitude boundary having
C           the same sign as the latitude of the input point---namely,
C           the element contains the pole, and latitude is excluded,
C           then then passing the radius check implies that the point
C           is close to the portion of the Z-axis contained in the
C           element.
C
C           If the radius check has been excluded because the point is
C           already nominally on one of the element's radius boundaries,
C           the passing the latitude test implies the point is close
C           to the element.
C
C           In all cases, as long as EXCLUD has been set appropriately,
C           the point is close to the element. We consider the point to
C           be in the expanded element.
C            
         END IF

      END IF
  
C
C     All tests that were commanded have been passed. The input
C     point is considered to be contained in the expanded volume
C     element.
C
      INSIDE = .TRUE.

      END 
