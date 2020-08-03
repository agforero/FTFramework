C$Procedure ZZINPDT ( DSK, in planetodetic element? )

      SUBROUTINE ZZINPDT ( P, BOUNDS, CORPAR, MARGIN, EXCLUD, INSIDE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of planetodetic coordinates for
C     inclusion in a specified planetodetic volume element. Ellipsoidal
C     surfaces are used in place of surfaces of constant altitude. The
C     test is performed using margins for the element.
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
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      MARGIN
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
 

      DOUBLE PRECISION      LATMRG
      PARAMETER           ( LATMRG = 1.D-8 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   Input point.
C     BOUNDS     I   Coordinate bounds of element.
C     CORPAR     I   Coordinate system parameters.
C     MARGIN     I   Margin used for inclusion testing.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C     LONALI     P   Longitude shift margin.
C     NONE       P   Meaning: exclude nothing.
C     LONIDX     P   Longitude index. Also used for exclusion.
C     LATIDX     P   Latitude index. Also used for exclusion.
C     ALTIDX     P   Altitude index. Also used for exclusion.
C     LATMRG     P   High/low latitude margin.
C
C$ Detailed_Input
C
C     P          is a point expressed in Cartesian coordinates. The
C                point is to be checked to determine whether it is
C                inside the volume element specified by BOUNDS.
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
C                If the longitude upper bound is not greater than the
C                longitude lower bound, a value greater than the upper
C                bound by 2*pi is used for the comparison.
C                
C
C     CORPAR     is an array of coordinate parameters; normally these
C                are obtained from a DSK descriptor. The first element
C                is the equatorial radius of the reference spheroid.
C                The second element is the flattening coefficient.
C                Any other elements are ignored by this routine.
C
C
C     MARGIN     is a fraction used to expand the volume element for
C                inclusion testing. 
C
C                   - The upper altitude bound is increased by 
C
C                        MARGIN*ABS( BOUNDS(2,3) )
C
C                     The lower altitude bound is decreased by 
C
C                        MARGIN*ABS( BOUNDS(1,3) )
C
C                     Note that these bounds are used as inputs to
C                     compute the radii of upper and lower bounding
C                     ellipsoids. The input point's position is tested
C                     for inclusion within or exclusion from these
C                     ellipsoids.
C
C
C                   - Latitude bounds are extended by MARGIN radians,
C                     within the interval -pi/2 : pi/2.
C
C                   - For any input point having latitude "not close"
C                     to pi/2 or -pi/2, the element's longitude
C                     interval is extended by
C
C                        MARGIN / cos(LAT) radians
C
C                     where LAT is the planetocentric [sic] latitude of
C                     the input point.
C         
C                     Here "close to" means "within LATMRG radians of."
C
C                   - For any input point having planetocentric
C                     latitude close to (within LATMRG radians of) pi/2
C                     or -pi/2, comparison against the longitude bounds
C                     of the volume element is omitted if MARGIN > 0.
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
C                coordinate. For example, if a point is on the outer
C                bounding ellipsoid, the level surface parameter of the
C                point with respect to the ellipsoid need not be used
C                in the comparison and in fact can't be meaningfully
C                compared, due to round-off errors.
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
C                by more than this amount, the value LON + 2*pi is
C                used in the longitude comparison, provided that the
C                comparison is not excluded.
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
C
C     LATMRG     is the margin used to determine whether the input
C                latitude is "close" to the poles. Units are radians.
C                The value of this parameter must be kept in sync
C                with the "greedy" margin declared in dsktol.inc.
C
C$ Exceptions
C
C     1)  If MARGIN is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     2)  If EXCLUD is outside of the range 0:3, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C     3)  If an error occurs while determining the planetodetic
C         coordinates of the input point, the error will be signaled by
C         a routine in the call tree of this routine.
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
C     See usage in ZZRYTPDT.
C
C$ Restrictions
C
C     1)  This is a private routine. It is meant to be used only by the
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
C-    SPICELIB Version 1.0.0, 13-FEB-2017 (NJB)
C
C        03-JUN-2016 (NJB)
C
C           Added standard traceback participation. Added check for
C           invalid values of EXCLUD.
C
C        05-MAR-2016 (NJB) 
C
C           Original version 03-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     test point against planetodetic element using margin
C
C-&




C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      PI
      DOUBLE PRECISION      TWOPI

      LOGICAL               FAILED
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

      INTEGER               LOWER
      PARAMETER           ( LOWER  = 1 )

      INTEGER               UPPER
      PARAMETER           ( UPPER  = 2 )


C
C     Numeric relation codes returned by ZZPDCMPL:
C     
      INTEGER               LT
      PARAMETER           ( LT = -1 )
      
      INTEGER               GT
      PARAMETER           ( GT =  1 )

C
C     The code EQ can be returned by ZZPDCMPL, but we make no
C     references to it, so it's not declared here. For the 
C     record, EQ is set to 0.
C

C
C     Local variables
C
      DOUBLE PRECISION      ALTBDS ( 2, 3 )
      DOUBLE PRECISION      AMNALT
      DOUBLE PRECISION      AMNLAT
      DOUBLE PRECISION      AMNLON
      DOUBLE PRECISION      AMXALT
      DOUBLE PRECISION      AMXLAT
      DOUBLE PRECISION      AMXLON
      DOUBLE PRECISION      DLON
      DOUBLE PRECISION      F
      DOUBLE PRECISION      HPI
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      LONMRG
      DOUBLE PRECISION      MAXALT
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MINALT
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      PCNLAT
      DOUBLE PRECISION      PI2
      DOUBLE PRECISION      R
      DOUBLE PRECISION      RE

      INTEGER               RELMIN
      INTEGER               RELMAX

      LOGICAL               ALPASS
      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  ALTBDS
      SAVE                  FIRST
      SAVE                  HPI
      SAVE                  PI2

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /
      DATA                  HPI   / -1.D0  /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZINPDT' )

      IF ( FIRST ) THEN

         HPI   = HALFPI()
         PI2   = TWOPI()
C
C        Initialize the local array used for altitude checks.
C
         ALTBDS( SOUTH, LATIDX ) =  - HALFPI() 
         ALTBDS( NORTH, LATIDX ) =    HALFPI()
         ALTBDS( WEST,  LONIDX ) =  - PI()
         ALTBDS( EAST,  LONIDX ) =    PI()
         ALTBDS( LOWER, ALTIDX ) =    0.D0
         ALTBDS( UPPER, ALTIDX ) =    0.D0
         
         FIRST = .FALSE.

      END IF


      IF (  ( EXCLUD .LT. 0 ) .OR. ( EXCLUD .GT. 3 )  ) THEN

         CALL SETMSG ( 'EXCLUD must be in the range 0:3 but was #.' )
         CALL ERRINT ( '#', EXCLUD                                  )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                     )
         CALL CHKOUT ( 'ZZINPDT'                                    )
         RETURN

      END IF

C
C     Get the planetocentric [sic] coordinates of the input point. The
C     latitude we obtain will be planetocentric. To emphasize this, we
C     use the name "PCNLAT."
C
      CALL RECLAT ( P, R, LON, PCNLAT )
C
C     RECLAT is error free, so we don't call FAILED() here.
C     

      IF ( MARGIN .EQ. 0.D0 ) THEN
C
C        ZZINPDT0 contains the logic required to determine whether
C        the input point is contained in the element.
C
         CALL ZZINPDT0 ( P,      LON,    
     .                   BOUNDS, CORPAR, EXCLUD, INSIDE )

         CALL CHKOUT ( 'ZZINPDT' )
         RETURN

      ELSE IF ( MARGIN .LT. 0.D0 ) THEN

         CALL SETMSG ( 'Margin must be non-negative but was #.' )
         CALL ERRDP  ( '#', MARGIN                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINPDT'                                )
         RETURN

      END IF

C
C     At this point a more detailed analysis is needed.
C   
C
C     Assume the point is outside to start. This allows us
C     to skip setting INSIDE when we find a boundary test
C     failure.
C
      INSIDE = .FALSE.

C
C     We'll use the shape parameters for latitude and altitude
C     comparisons.
C
      RE = CORPAR(1)
      F  = CORPAR(2)

C
C     Get local copies of the coordinate bounds. Don't normalize the
C     longitude bounds until we know we need them.
C 
      MINLAT = BOUNDS( SOUTH, LATIDX )
      MAXLAT = BOUNDS( NORTH, LATIDX )

      MINALT = BOUNDS( LOWER, ALTIDX ) 
      MAXALT = BOUNDS( UPPER, ALTIDX ) 

C
C     Compare coordinates to adjusted latitude boundaries.
C
      IF ( EXCLUD .NE. LATIDX ) THEN
C
C        Create adjusted latitude bounds.
C
         AMNLAT = MAX( -HPI-ANGMRG, MINLAT - MARGIN )
         AMXLAT = MIN(  HPI+ANGMRG, MAXLAT + MARGIN )

C
C        Compare the latitude of the input point to the bounds.
C
         CALL ZZPDCMPL ( RE, F, P, AMNLAT, RELMIN )
         CALL ZZPDCMPL ( RE, F, P, AMXLAT, RELMAX )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZINPDT' )
            RETURN
         END IF


         IF ( ( RELMIN .EQ. LT ) .OR. ( RELMAX .EQ. GT ) ) THEN
C
C           The latitude of P is strictly outside of the element's
C           latitude bounds.
C
            CALL CHKOUT ( 'ZZINPDT' )
            RETURN
         END IF

      END IF

C
C     Test the point for inclusion in the region between the bounding
C     ellipsoids that act as proxies for altitude boundaries.
C
      IF ( EXCLUD .NE. ALTIDX ) THEN
C
C        Extract altitude bounds from the segment descriptor.
C
         MINALT = BOUNDS(LOWER, ALTIDX)
         MAXALT = BOUNDS(UPPER, ALTIDX)
C
C        Adjust altitude bounds to account for the margin.
C        
         AMNALT = MINALT  -  MARGIN * ABS(MINALT)
         AMXALT = MAXALT  +  MARGIN * ABS(MAXALT)
C
C        Set up a "boundary" array so that we can use ZZINPDT0
C        to do the altitude check for us. We'll exclude longitude
C        tests; the latitude test is set up for an automatic pass
C        (the latitude range of ALTBDS is [-pi/2, pi/2]).
C
         ALTBDS( LOWER, ALTIDX ) = AMNALT
         ALTBDS( UPPER, ALTIDX ) = AMXALT
 
         CALL ZZINPDT0 ( P, LON, ALTBDS, CORPAR, LONIDX, ALPASS )

         IF ( .NOT. ALPASS ) THEN

            CALL CHKOUT ( 'ZZINPDT' )
            RETURN

         END IF

      END IF

C
C     At this point, the input altitude and latitude are within the
C     adjusted bounds, if their tests haven't been excluded by
C     the caller.
C
C     Perform longitude tests, unless they're excluded by the 
C     caller.
C
      IF ( EXCLUD .NE. LONIDX ) THEN      
C
C        Start out by normalizing the element's longitude bounds.
C
         CALL ZZNRMLON ( BOUNDS(LOWER,LONIDX),  BOUNDS(UPPER,LONIDX), 
     .                   ANGMRG,       MINLON,  MAXLON               )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZINPDT' )
            RETURN
         END IF
C
C        Set the margin to be used for longitude interval 
C        inclusion tests.
C
         LONMRG = MAX( ABS(ANGMRG), ABS(MARGIN) )

C
C        We have a special case for segments that include the poles. If
C        the latitude and altitude of the input point are within
C        bounds, and if the latitude of the point is close enough to a
C        pole we consider the point to be included in the segment,
C        regardless of the point's longitude. All other points get the
C        normal longitude test.
C
C        We use planetocentric latitude to determine whether the point
C        is "close" to a pole.
C
         CALL ZZPDCMPL ( RE, F, P,  HPI-LATMRG, RELMAX )
         CALL ZZPDCMPL ( RE, F, P, -HPI+LATMRG, RELMIN )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZINPDT' )
            RETURN
         END IF

         IF (       ( RELMAX .NE. GT )
     .        .AND. ( RELMIN .NE. LT )  ) THEN
C
C           This is the usual case: the latitude of the input point is
C           bounded away from the poles.
C
C           Check the point's longitude against the segment's longitude
C           bounds.
C
C           We'll scale the longitude margin to compensate for the
C           latitude of the input point. Note that the division
C           below is safe; presuming a reasonable value of MARGIN,
C           we know that 
C
C              DLON << 1     
C         
C           Note that we use planetocentric latitude for scaling the
C           longitude margin. This substitution (for planetodetic
C           latitude) is adequate for this purpose.
C
            DLON   =  LONMRG / MAX( ABS(COS(PCNLAT)), LATMRG ) 

            AMNLON = MINLON - DLON
            AMXLON = MAXLON + DLON

C
C           Now move the input point's longitude into range, if
C           necessary, so we can make a valid comparison against
C           the longitude bounds.
C     
            IF ( LON .LT. AMNLON ) THEN

               IF (  LON  .LT. ( AMNLON - LONALI )  ) THEN
C
C                 See whether an aliased version of LON is a match.
C
                  LON = LON + PI2

               ELSE
C
C                 Consider LON to be a match with the lower bound.
C
                  LON = AMNLON

               END IF

            ELSE IF ( LON .GT. AMXLON ) THEN

               IF (  LON  .GT. ( AMXLON + LONALI )  ) THEN
C
C                 See whether an aliased version of LON is a match.
C
                  LON = LON - PI2

               ELSE
C
C                 Consider LON to be a match with the upper bound.
C
                  LON = AMXLON

               END IF

            END IF


            IF ( ( LON .LT. AMNLON ) .OR. ( LON .GT. AMXLON ) ) THEN

               CALL CHKOUT ( 'ZZINPDT' )
               RETURN

            END IF              

         ELSE

C
C           The latitude of the input point is close to one of the
C           poles.
C
C           This is a no-op case.
C
C           The input point has already passed whichever of the
C           altitude and latitude tests that were not excluded.
C
C           If the element has a non-degenerate latitude boundary
C           having the same sign as the latitude of the input point,
C           and if latitude is excluded because the input point is
C           already nominally on that boundary, then passing the
C           altitude check implies that the point is close to the
C           element.
C
C           If the element has a degenerate latitude boundary having
C           the same sign as the latitude of the input point---namely,
C           the element contains the pole, and latitude is excluded,
C           then then passing the altitude check implies that the point
C           is close to the portion of the Z-axis contained in the
C           element.
C
C           If the altitude check has been excluded because the point
C           is already nominally on one of the element's altitude
C           boundaries, the passing the latitude test implies the point
C           is close to the element.
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

      CALL CHKOUT ( 'ZZINPDT' )
      RETURN
      END 
