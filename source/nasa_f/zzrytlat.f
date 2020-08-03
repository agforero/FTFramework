C$Procedure ZZRYTLAT ( DSK, ray touches latitudinal element )
 
      SUBROUTINE ZZRYTLAT ( VERTEX, RAYDIR, BOUNDS, 
     .                      MARGIN, NXPTS,  XPT    )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find nearest intersection to a given ray's vertex of that ray and
C     at latitudinal volume element. If the vertex is inside the
C     element, the vertex is considered to be the solution.
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
C     GEOMETRY
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsktol.inc'

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      MARGIN
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT    ( 3 )

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
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     BOUNDS     I   Bounds of latitudinal volume element.
C     MARGIN     I   Margin used for element expansion.
C     NXPTS      O   Number of intercept points.
C     XPT        O   Intercept.
C     LONIDX     P   Longitude index. 
C     LATIDX     P   Latitude index. 
C     RADIDX     P   Radius index. 
C
C$ Detailed_Input
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray direction must be represented
C                in the reference frame to which the volume element
C                boundaries correspond. The vertex is considered to be
C                an offset from the center of the reference frame
C                associated with the volume element.
C 
C     BOUNDS     is a 2x3 array containing the bounds of a latitudinal
C                volume element. Normally this is the coverage boundary
C                of a DSK segment. In the element
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
C                If the longitude upper bound is not greater than the
C                longitude lower bound, a value greater than the upper
C                bound by 2*pi is used for the comparison.
C
C
C     MARGIN     is a scale factor used to effectively expand the
C                segment boundaries so as to include intersections
C                that lie slightly outside the volume element.
C
C
C$ Detailed_Output
C
C     XPT        is the intercept of the ray on the surface described
C                by the input coordinate bounds, if such an intercept
C                exists. If the ray intersects the surface at multiple
C                points, the one closest to the ray's vertex is
C                selected. XPT is valid if and only if FOUND is .TRUE.
C
C                XPT is expressed in the reference frame associated
C                with the inputs VERTEX and RAYDIR. XPT represents
C                an offset from the origin of the coordinate system.
C
C                XPT is valid only if NXPTS is set to 1.
C
C
C     NXPTS      is the number of intercept points of the ray and
C                the volume element. 
C
C                Currently there are only two possible values for
C                NXPTS:
C
C                   1 for an intersection
C                   0 for no intersection
C    
C                If the vertex is inside the element, NXPTS is
C                set to 1.
C
C$ Parameters
C
C
C     LONIDX     is the index of longitude in the second dimension of
C                BOUNDS.  
C
C     LATIDX     is the index of latitude in the second dimension of
C                BOUNDS.  
C
C     RADIDX     is the index of radius in the second dimension of
C                BOUNDS.  
C$ Exceptions
C
C     For the sake of run-time efficiency, this routine does not
C     participate in call tracing. However, it's possible for
C     routines called by this routine to signal errors.
C
C$ Files
C
C     None. However, the input volume element boundaries normally have
C     been obtained from a segment in a loaded DSK file.
C     
C$ Particulars
C
C     This routine sits on top of data DSK type-specific ray-segment
C     intercept routines such as DSKX02.
C
C$ Examples
C
C     See usage in ZZDSKBUX.
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
C-    SPICELIB Version 1.0.0, 18-JAN-2017 (NJB) 
C
C        Based on ZZRYXLAT 20-MAR-2015 (NJB)
C
C        Original 01-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find intercept of ray on latitudinal volume element
C
C-&



C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VDOT
      DOUBLE PRECISION      VNORM 

      LOGICAL               FAILED

C
C     Local parameters
C

C
C     Radius expansion factor:
C
      DOUBLE PRECISION      RADFAC
      PARAMETER           ( RADFAC = 1.1D0 )


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
      DOUBLE PRECISION      ANGLE 
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      EASTB  ( 3 )
      DOUBLE PRECISION      EBACK  ( 3 )
      DOUBLE PRECISION      ENDPT2 ( 3 )
      DOUBLE PRECISION      LONCOV
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      MNDIST
      DOUBLE PRECISION      NEGDIR ( 3 )
      DOUBLE PRECISION      ORIGIN ( 3 )
      DOUBLE PRECISION      S
      DOUBLE PRECISION      SRFX   ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      VLAT
      DOUBLE PRECISION      VLON
      DOUBLE PRECISION      VR
      DOUBLE PRECISION      WBACK  ( 3 )
      DOUBLE PRECISION      WESTB  ( 3 )
      DOUBLE PRECISION      XPT2   ( 3 )
      DOUBLE PRECISION      Z      ( 3 )

      INTEGER               NX

      LOGICAL               FOUND
      LOGICAL               INSIDE
      LOGICAL               XIN

C
C     Saved variables
C
      SAVE                  ORIGIN
      SAVE                  Z

C
C     Initial values
C     
      DATA                  ORIGIN / 3 *   0.D0 /
      DATA                  Z      / 0.D0, 0.D0, 1.D0 /


C
C     CAUTION: this routine doesn't check in, so it won't
C     appear in the traceback if a lower-level routine
C     signals an error.
C

C
C     Determine whether the vertex is inside the element.
C     Use double the margin for this test, since we don't
C     want to have false negative tests for rays having
C     vertices lying on the expanded element boundary.
C
      CALL ZZINLAT ( VERTEX, BOUNDS, 2*MARGIN, 0, INSIDE )

      IF ( FAILED() ) THEN
         RETURN
      END IF

      IF ( INSIDE ) THEN
C
C        We know the answer.
C
         NXPTS = 1

         CALL VEQU ( VERTEX, XPT )

         RETURN

      END IF

C
C     Extract the segment's coordinate bounds into easily 
C     readable variables.
C
      MINR = BOUNDS( INNER, RADIDX ) 
      MAXR = BOUNDS( OUTER, RADIDX ) 
C
C     Normalize the longitude bounds. After this step, the bounds will
C     be in order and differ by no more than 2*pi.
C
      CALL ZZNRMLON( BOUNDS(WEST,LONIDX), BOUNDS(EAST,LONIDX), ANGMRG,
     .               MINLON,              MAXLON                       )

      IF ( FAILED() ) THEN
         RETURN
      END IF

      MINLAT = BOUNDS( SOUTH, LATIDX )
      MAXLAT = BOUNDS( NORTH, LATIDX )

C
C     The vertex is outside the element. 
C
C     Indicate no intersection to start.
C
      NXPTS  = 0
C
C     We'll use a unit length copy of the ray's direction vector.
C
      CALL VHAT ( RAYDIR, UDIR )
C
C     Initialize the distance to the closest solution. We'll keep track
C     of this quantity in order to compare competing solutions.
C
      MNDIST = DPMAX()
C
C     Find the intersection of the ray and outer bounding sphere, if
C     possible. Often this intersection is the closest to the vertex.
C     If the intersection exists and is on the boundary of the element,
C     it's a winner.
C
      CALL ZZRYXSPH ( VERTEX, UDIR, MAXR, SRFX, FOUND )

      IF ( .NOT. FOUND ) THEN
C
C        There are no intersections. The ray cannot hit the
C        volume element.
C
         RETURN

      END IF
C
C     Get the latitudinal coordinates of the ray's vertex.
C
      CALL RECLAT ( VERTEX, VR, VLON, VLAT )

C
C     The ray hits the outer bounding sphere. See whether
C     the longitude and latitude are within bounds, taking
C     the margin into account. Exclude the radius coordinate
C     from testing.
C
      CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, RADIDX, XIN )

      IF ( FAILED() ) THEN
         RETURN
      END IF


      IF ( XIN ) THEN
C
C        This solution is a candidate.
C
         CALL VEQU ( SRFX, XPT )
         NXPTS  = 1

         IF ( VR .GT. MAXR ) THEN
C
C           The vertex is outside this sphere, and the segment
C           lies within the sphere.
C
C           No other intersection can be closer to the vertex;
C           we don't need to check the other surfaces.
C
            RETURN

         ELSE
C
C           We have a possible solution.
C 
            MNDIST = VDIST( VERTEX, XPT )

         END IF

      END IF

C
C     So far there may be a candidate solution. We'll try the latitude
C     boundaries next.
C
C     For testing intersections with the latitude boundaries, we'll
C     need a far endpoint for the line segment on which to perform the
C     test.
C
      S = VNORM(VERTEX) + RADFAC * MAXR

      CALL VLCOM ( 1.D0, VERTEX, S, UDIR, ENDPT2 )

C
C     Now try the upper latitude bound. We can skip this test
C     if the upper bound is pi/2 radians.
C        
      IF ( MAXLAT .LT. HALFPI() ) THEN
C
C        Let ANGLE be the angular separation of the surface of latitude
C        MAXLAT and the +Z axis. Note that the surface might be the
C        lower nappe of the cone.
C
         ANGLE = MAX ( 0.D0,  HALFPI() - MAXLAT )

         CALL INCNSG ( ORIGIN, Z,  ANGLE, VERTEX,
     .                 ENDPT2, NX, SRFX,  XPT2    )

         IF ( FAILED() ) THEN
            RETURN
         END IF

         IF ( NX .GE. 1 ) THEN
C
C           See whether SRFX is in the element.
C
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, LATIDX, XIN )

            IF ( FAILED() ) THEN
               RETURN
            END IF

            IF ( XIN ) THEN
C
C              SRFX is a candidate solution.
C
               DIST = VDIST( VERTEX, SRFX )

               IF ( DIST .LT. MNDIST ) THEN

                  CALL VEQU ( SRFX, XPT )
                     
                  NXPTS  = 1

                  IF ( VLAT .GT. MAXLAT ) THEN
C
C                    Since the latitude of the vertex is greater than
C                    MAXLAT, this is the best solution, since the
C                    volume element is on the other side of the maximum
C                    latitude boundary.
C
                     RETURN

                  END IF
C
C                 This is the best solution seen so far.
C
                  MNDIST = DIST

               END IF

            END IF

            IF ( NX .EQ. 2 ) THEN
C
C              Check the second solution as well.
C
               CALL ZZINLAT ( XPT2, BOUNDS, MARGIN, LATIDX, XIN )

               IF ( FAILED() ) THEN
                  RETURN
               END IF


               IF ( XIN ) THEN
C
C                 XPT2 is a candidate solution.
C
                  DIST = VDIST( VERTEX, XPT2 )

                  IF ( DIST .LT. MNDIST ) THEN

                     CALL VEQU ( XPT2, XPT )                    

                     NXPTS  = 1
                     MNDIST = DIST
C
C                    This is the best solution seen so far.
C                    However, it's not necessarily the best
C                    solution. So we continue.
C
                  END IF

               END IF

            END IF
C
C           We've handled the second root, if any.
C
         END IF
C
C        We're done with the upper latitude boundary.
C
      END IF

C
C     Try the lower latitude bound. We can skip this test if the lower
C     bound is -pi/2 radians.
C        
      IF (  MINLAT .GT. -HALFPI() ) THEN
C
C        Let ANGLE be the angular separation of the surface
C        of latitude MINLAT and the +Z axis. Note that the
C        surface might be the lower nappe of the cone.
C
         ANGLE = HALFPI() - MINLAT

         CALL INCNSG ( ORIGIN, Z,  ANGLE, VERTEX,
     .                 ENDPT2, NX, SRFX,  XPT2    )

         IF ( FAILED() ) THEN
            RETURN
         END IF


         IF ( NX .GE. 1 ) THEN
C
C           See whether SRFX is in the element.
C
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, LATIDX, XIN )

            IF ( FAILED() ) THEN
               RETURN
            END IF

            IF ( XIN ) THEN
C
C              SRFX is a candidate solution.
C
               DIST = VDIST( VERTEX, SRFX )

               IF ( DIST .LT. MNDIST ) THEN

                  CALL VEQU ( SRFX, XPT )
                     
                  NXPTS  = 1

                  IF ( VLAT .LT. MINLAT ) THEN
C
C                    Since the latitude of the vertex is less than
C                    MINLAT, this is the best solution, since the
C                    volume element is on the other side of the minimum
C                    latitude boundary.
C
                     RETURN
                  
                  END IF
C
C                 This is the best solution seen so far.
C
                  MNDIST = DIST

               END IF

            END IF

            IF ( NX .EQ. 2 ) THEN
C
C              Check the second solution as well.
C
               CALL ZZINLAT ( XPT2, BOUNDS, MARGIN, LATIDX, XIN )

               IF ( FAILED() ) THEN
                  RETURN
               END IF

               IF ( XIN ) THEN
C
C                 XPT2 is a candidate solution.
C
                  DIST = VDIST( VERTEX, XPT2 )

                  IF ( DIST .LT. MNDIST ) THEN

                     CALL VEQU ( XPT2, XPT )
                     
                     NXPTS  = 1
                     MNDIST = DIST
C
C                    This is the best solution seen so far.
C                    However, it's not necessarily the best
C                    solution. So we continue.
C
                     RETURN

                  END IF

               END IF

            END IF

         END IF
C
C        We're done with the lower latitude boundary.
C 
      END IF

C
C     Perform longitude boundary checks if the coverage is not
C     2*pi radians. Note that MAXLON > MINLON at this point.
C     
      LONCOV = MAXLON - MINLON

      IF ( COS(LONCOV) .LT. 1.D0 ) THEN
C
C        We have distinct longitude boundaries. Go to work.
C
C
C        Check the longitude boundaries. Try the plane of western
C        longitude first.
C
C        The vector WESTB is an outward-facing vector normal to 
C        the west boundary.
C
         CALL VPACK ( SIN(MINLON), -COS(MINLON), 0.D0, WESTB )


         S = RADFAC * ( VNORM(VERTEX) + MAXR )


         CALL ZZINRYPL ( VERTEX, UDIR, WESTB, 0.D0, S, NX, SRFX )

         IF ( NX .EQ. 1 ) THEN
C
C           We have one point of intersection. Determine whether it's a
C           candidate solution.  Don't use longitude in the following
C           inclusion test. Note that we'll perform a separate check
C           later in place of the longitude check.
C
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, LONIDX, XIN )

            IF ( FAILED() ) THEN
               RETURN
            END IF

            
            IF ( XIN ) THEN
C
C              Make sure the intercept is not too far on the 
C              "wrong" side of the Z axis.
C
               CALL UCRSS ( WESTB, Z, WBACK )

               IF (  VDOT(SRFX, WBACK) .LT. (MARGIN*MAXR)  ) THEN
C
C                 The intercept is either on the same side of the Z
C                 axis as the west face of the segment, or is very
C                 close to the Z axis.
C
                  DIST = VDIST( VERTEX, SRFX )

                  IF ( DIST .LT. MNDIST ) THEN
C
C                    Record the intercept, distance, and surface index.
C
                     CALL VEQU ( SRFX, XPT )

                     NXPTS  = 1
                     MNDIST = DIST

                  END IF

               END IF

            END IF

         END IF
C
C        We're done with the western boundary.
C
C
C        Try the plane of eastern longitude next.
C
C        The vector EASTB is an outward-facing vector normal to 
C        the east boundary.

         CALL VPACK ( -SIN(MAXLON), COS(MAXLON), 0.D0, EASTB )

         CALL ZZINRYPL ( VERTEX, UDIR, EASTB, 0.D0, S, NX, SRFX )

         IF ( NX .EQ. 1 ) THEN
C
C           We have one point of intersection. Determine whether it's a
C           candidate solution.
C
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, LONIDX, XIN )
      
            IF ( FAILED() ) THEN
               RETURN
            END IF
      
            IF ( XIN ) THEN
C
C              Make sure the intercept is not too far on the "wrong"
C              side of the Z axis.
C
               CALL UCRSS ( Z, EASTB, EBACK )

               IF (  VDOT(SRFX, EBACK) .LT. (MARGIN*MAXR)  ) THEN
C
C                 The intercept is either on the same side of the Z
C                 axis as the east face of the segment, or is very
C                 close to the Z axis.
C
                  DIST = VDIST( VERTEX, SRFX )

                  IF ( DIST .LT. MNDIST ) THEN
C
C                    Record the intercept, distance, and surface index.
C 
                     CALL VEQU ( SRFX, XPT )

                     NXPTS  = 1
                     MNDIST = DIST

                  END IF

               END IF

            END IF

         END IF

      END IF
C
C     End of longitude boundary checks.
C

C
C     Find the intersection of the ray and inner bounding
C     sphere, if possible.  
C
      IF ( MINR .GT. 0.D0 ) THEN

         CALL ZZRYXSPH ( VERTEX, UDIR, MINR, SRFX, FOUND )

         IF ( FOUND ) THEN
C
C           See whether this solution is in the element.
C
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, RADIDX, XIN )

            IF ( FAILED() ) THEN
               RETURN
            END IF
    
            IF ( XIN ) THEN

               DIST = VDIST( VERTEX, SRFX )

               IF ( DIST .LT. MNDIST ) THEN
C
C                 Record the intercept, distance, and surface index.
C
                  CALL VEQU ( SRFX, XPT )

                  NXPTS  = 1
                  MNDIST = DIST

               END IF

            END IF
 
         END IF

C
C        Unlike the outer sphere, either intersection of the ray with
C        the inner sphere can be an optimal solution. We'll test for
C        the case where the intercept further from the ray's vertex is
C        the correct solution.
C
         CALL VMINUS( UDIR, NEGDIR )

         CALL ZZRYXSPH ( ENDPT2, NEGDIR, MINR, SRFX, FOUND )

         IF ( FOUND ) THEN
               
            CALL ZZINLAT ( SRFX, BOUNDS, MARGIN, RADIDX, XIN )
           
            IF ( FAILED() ) THEN
               RETURN
            END IF
 
            IF ( XIN ) THEN

               DIST = VDIST( VERTEX, SRFX )

               IF ( DIST .LT. MNDIST ) THEN
C
C                 Record the intercept, distance, and surface index.
C
                  CALL VEQU ( SRFX, XPT )

                  NXPTS  = 1
                  MNDIST = DIST

               END IF

            END IF
 
         END IF

      END IF

C
C     NXPTS and XPT are set.
C
      END 
