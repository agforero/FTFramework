C$Procedure ZZRYTPDT ( DSK, ray touches planetodetic element )
 
      SUBROUTINE ZZRYTPDT ( VERTEX, RAYDIR, BOUNDS, 
     .                      CORPAR, MARGIN, NXPTS,  XPT )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find nearest intersection to a given ray's vertex of the ray and
C     a planetodetic volume element. If the vertex is inside the
C     element, the vertex is considered to be the solution.
C      
C     In the computation performed by this routine, ellipsoidal
C     surfaces are used, instead of surfaces of constant altitude, to
C     define boundaries of planetodetic volume elements. The element
C     defined by the input boundaries is contained in the element
C     bounded by the input latitude and longitude boundaries and by the
C     ellipsoidal surfaces.
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
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      MARGIN
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT    ( 3 )

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
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     BOUNDS     I   Bounds of planetodetic volume element.
C     CORPAR     I   Coordinate parameters.
C     MARGIN     I   Margin used for element expansion.
C     NXPTS      O   Number of intercept points.
C     XPT        O   Intercept.
C     LONIDX     P   Longitude index. 
C     LATIDX     P   Latitude index. 
C     ALTIDX     P   Altitude index. 
C
C$ Detailed_Input
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray direction must be represented
C                in the reference frame to which the planetodetic
C                volume element boundaries correspond. The vertex is
C                considered to be an offset from the center of the
C                reference frame associated with the element.
C
C     BOUNDS     is a 2x3 array containing the bounds of a planetodetic
C                volume element. Normally this is the coverage boundary
C                of a DSK segment. In the element
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
C     RE,
C     F          are, respectively, the equatorial radius and
C                flattening coefficient associated with the
C                planetodetic coordinate system in which the input
C                volume element is described.
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
C                by the segment, if such an intercept exists. If the
C                ray intersects the surface at multiple points, the
C                one closest to the ray's vertex is selected. XPT is
C                valid if and only if FOUND is .TRUE.
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
C     LONIDX     is the index of longitude in the second dimension of
C                BOUNDS.  
C
C     LATIDX     is the index of latitude in the second dimension of
C                BOUNDS.  
C
C     ALTIDX     is the index of altitude in the second dimension of
C                BOUNDS.  
C$ Exceptions
C
C     1)  If MARGIN is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     2)  If the input ray direction vector is zero, the error
C         SPICE(ZEROVECTOR) will be signaled.
C
C     3)  Any errors that occur while calculating the ray-surface
C         intercept will be signaled by routines in the call tree
C         of this routine.
C
C$ Files
C
C     None. However, the input segment boundaries normally have
C     been obtained from a loaded DSK file.
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
C-    SPICELIB Version 1.0.0, 19-JAN-2017 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     find intercept of ray on planetodetic volume element
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
      DOUBLE PRECISION      VSEP

      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO
      LOGICAL               ZZPDPLTC

C
C     Local parameters
C

C
C     Altitude expansion factor:
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

      INTEGER               LOWER
      PARAMETER           ( LOWER  = 1 )

      INTEGER               UPPER
      PARAMETER           ( UPPER  = 2 )

      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

C
C     Local variables
C
      DOUBLE PRECISION      AMNALT
      DOUBLE PRECISION      AMXALT
      DOUBLE PRECISION      ANGLE 
      DOUBLE PRECISION      APEX   ( 3 )
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      EASTB  ( 3 )
      DOUBLE PRECISION      EBACK  ( 3 )
      DOUBLE PRECISION      EMAX
      DOUBLE PRECISION      EMIN
      DOUBLE PRECISION      ENDPT2 ( 3 )
      DOUBLE PRECISION      F
      DOUBLE PRECISION      LONCOV
      DOUBLE PRECISION      MAXALT
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MINALT
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MNDIST
      DOUBLE PRECISION      NEGDIR ( 3 )
      DOUBLE PRECISION      PMAX
      DOUBLE PRECISION      PMIN
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RP
      DOUBLE PRECISION      S
      DOUBLE PRECISION      SRFX   ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      VTXANG
      DOUBLE PRECISION      VTXLVL
      DOUBLE PRECISION      VTXOFF ( 3 )
      DOUBLE PRECISION      WBACK  ( 3 )
      DOUBLE PRECISION      WESTB  ( 3 )
      DOUBLE PRECISION      XPT2   ( 3 )
      DOUBLE PRECISION      XINCPT   
      DOUBLE PRECISION      YINCPT   
      DOUBLE PRECISION      Z      ( 3 )

      INTEGER               NX

      LOGICAL               FOUND
      LOGICAL               INSIDE
      LOGICAL               XIN
      LOGICAL               XVAL1
      LOGICAL               XVAL2

C
C     Saved variables
C
      SAVE                  Z

C
C     Initial values
C     
      DATA                  Z      / 0.D0, 0.D0, 1.D0 /



      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'ZZRYTPDT' )


      IF ( MARGIN .LT. 0.D0 ) THEN

         CALL SETMSG ( 'Margin must be non-negative but was #.' )
         CALL ERRDP  ( '#', MARGIN                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZRYTPDT'                               )
         RETURN

      END IF

      IF ( VZERO(RAYDIR) ) THEN

         CALL SETMSG ( 'The ray''s direction was the zero vector.' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                         )
         CALL CHKOUT ( 'ZZRYTPDT'                                  )
         RETURN

      END IF


C
C     Determine whether the vertex is inside the element.
C
      CALL ZZINPDT ( VERTEX, BOUNDS, CORPAR, MARGIN, NONE, INSIDE )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         

      IF ( INSIDE ) THEN
C
C        We know the answer.
C
         NXPTS = 1

         CALL VEQU ( VERTEX, XPT )

         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN

      END IF

C
C     Get semi-axis lengths of the reference spheroid.
C
      RE = CORPAR(1)
      F  = CORPAR(2)
      RP = RE * ( 1.D0 - F )

C
C     Extract the segment's coordinate bounds into easily 
C     readable variables.
C
      MINALT = BOUNDS( LOWER, ALTIDX ) 
      MAXALT = BOUNDS( UPPER, ALTIDX ) 
C
C     Normalize the longitude bounds. After this step, the bounds will
C     be in order and differ by no more than 2*pi.
C
      CALL ZZNRMLON( BOUNDS(WEST,LONIDX), BOUNDS(EAST,LONIDX), ANGMRG,
     .               MINLON,              MAXLON                      )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         

      MINLAT = BOUNDS( SOUTH, LATIDX )
      MAXLAT = BOUNDS( NORTH, LATIDX )

C
C     Compute adjusted altitude bounds, taking margin into
C     account.
C        
      AMNALT = MINALT  -  MARGIN * ABS(MINALT)
      AMXALT = MAXALT  +  MARGIN * ABS(MAXALT)

C
C     Generate semi-axis lengths of inner and outer bounding
C     ellipsoids.
C
      IF ( RE .GE. RP ) THEN
C
C        The reference spheroid is oblate.
C         
         CALL ZZELLBDS ( RE,   RP,   AMXALT, AMNALT, 
     .                   EMAX, PMAX, EMIN,   PMIN   )

      ELSE
C
C        The reference spheroid is prolate.
C         
         CALL ZZELLBDS ( RP,   RE,   AMXALT, AMNALT, 
     .                   PMAX, EMAX, PMIN,   EMIN   )

      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         

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
C     Find the intersection of the ray and outer bounding ellipsoid, if
C     possible. Often this intersection is the closest to the vertex.
C     If the intersection exists and is on the boundary of the element,
C     it's a winner.
C
      CALL SURFPT ( VERTEX, UDIR, EMAX, EMAX, PMAX, SRFX, FOUND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         

      IF ( .NOT. FOUND ) THEN
C
C        There are no intersections. The ray cannot hit the volume
C        element.
C
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN

      END IF
 
C
C     The ray hits the outer bounding ellipsoid. See whether
C     the longitude and latitude are within bounds, taking
C     the margin into account. Exclude the altitude coordinate
C     from testing.
C
      CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, ALTIDX, XIN )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         

      IF ( XIN ) THEN
C
C        This solution is a candidate.
C
         CALL VEQU ( SRFX, XPT )
         NXPTS  = 1

C
C        Find the level surface parameter of the vertex relative
C        to the adjusted outer bounding ellipsoid.
C
         VTXLVL =     ( VERTEX(1)/EMAX )**2  
     .             +  ( VERTEX(2)/EMAX )**2  
     .             +  ( VERTEX(3)/PMAX )**2

         IF ( VTXLVL .GT. 1.D0 ) THEN
C
C           The vertex is outside this ellipsoid, and the DSK segment
C           lies within the ellipsoid.
C
C           No other intersection can be closer to the vertex;
C           we don't need to check the other surfaces.
C
            CALL CHKOUT ( 'ZZRYTPDT' )
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
      MAXR = MAX ( EMAX, PMAX )

      S    = VNORM(VERTEX) + RADFAC * MAXR

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

C
C        Compute the Z coordinate of the apex of the latitude cone.
C
         CALL ZZELNAXX ( RE, RP, MAXLAT, XINCPT, YINCPT )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         

         APEX(1) = 0.D0
         APEX(2) = 0.D0
         APEX(3) = YINCPT

C
C        Find the offset of the ray's vertex from the cone's apex,
C        and find the angular separation of the offset from the +Z
C        axis. This separation enables us to compare the latitude of
C        the vertex to the latitude boundary without making a RECGEO
C        call to compute the planetodetic coordinates of the vertex.
C
C        (The comparison will be done later.)
C
         CALL VSUB ( VERTEX, APEX, VTXOFF )

         VTXANG = VSEP ( VTXOFF, Z )

C
C        Check for intersection of the ray with the latitude cone.
C
         CALL INCNSG ( APEX, Z, ANGLE, VERTEX, ENDPT2, NX, SRFX, XPT2 )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         

C
C        Unlike the case of latitudinal coordinates, for planetodetic
C        coordinates, the surface of the latitude cone does not
C        coincide with the set of points having that latitude (which is
C        equal to pi/2 - the cone's angular separation from the +Z
C        axis). The subset of the cone having the specified latitude is
C        truncated by the X-Y plane. If we ignore round-off errors, we
C        can assert that the Z-coordinate of a point having the given
C        planetodetic latitude must match the direction of the nappe of
C        the cone: positive if ANGLE < pi/2, negative if ANGLE > pi/2,
C        and 0 if ANGLE = pi/2.
C
C        However, we cannot ignore round-off errors. For a cone having
C        angle from its central axis of nearly pi/2, it's possible for
C        a valid ray-cone intercept to be on the "wrong" side of the 
C        X-Y plane due to round-off errors. So we use a more robust
C        check to determine whether an intercept should be considered
C        to have the same latitude as the cone.
C
C        Check all intercepts.
C
         IF ( NX .GT. 0 ) THEN            
C
C           Check the first intercept.
C
            XVAL1 =  ZZPDPLTC( RE, F, SRFX, MAXLAT )
            XVAL2 =  .FALSE.
         
            IF ( NX .EQ. 2 ) THEN
C
C              Check the second intercept.
C
               XVAL2 =  ZZPDPLTC( RE, F, XPT2, MAXLAT )

            END IF

            IF ( XVAL1 .AND. ( .NOT. XVAL2 ) ) THEN

               NX = 1

            ELSE IF (  XVAL2  .AND.  (.NOT. XVAL1 )  ) THEN
C
C              Only the second solution is valid. Overwrite
C              the first.
C
               NX = 1

               CALL VEQU( XPT2, SRFX )

            ELSE IF (  ( .NOT. XVAL1 )  .AND.  ( .NOT. XVAL2 )  ) THEN
C
C              Neither solution is valid.
C
               NX = 0

            END IF
              
         END IF            
  


         
         IF ( NX .GE. 1 ) THEN
C
C           The ray intercept SRFX lies on the upper latitude boundary.
C
C           See whether SRFX meets the longitude and proxy altitude
C           constraints.
C
            CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, LATIDX, XIN )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZRYTPDT' )
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


                  IF ( VTXANG .LT. ANGLE ) THEN

                     IF (      ( MAXLAT    .LT. 0.D0 )
     .                    .OR. ( VERTEX(3) .GT. 0.D0 )  ) THEN
C
C                       If MAXLAT is negative, the vertex offset
C                       being outside the cone is enough to 
C                       guarantee the planetodetic latitude of the
C                       vertex is greater than that of the cone.
C
C                       If MAXLAT is non-negative, the angle of the
C                       vertex offset relative to the +Z axis is not
C                       enough; we need the vertex to lie above the 
C                       X-Y plane as well. 
C
C                       Getting here means one of these conditions
C                       was met.
C
C                       Since the latitude of the vertex is greater
C                       than MAXLAT, this is the best solution, since
C                       the volume element is on the other side of the
C                       maximum latitude boundary.
C
                        CALL CHKOUT ( 'ZZRYTPDT' )
                        RETURN

                     END IF

                  END IF
C
C                 This is the best solution seen so far, but we 
C                 need to check the remaining boundaries.
C
                  MNDIST = DIST

               END IF

            END IF

            IF ( NX .EQ. 2 ) THEN
C
C              Check the second solution as well.
C
               CALL ZZINPDT ( XPT2,   BOUNDS, CORPAR,
     .                        MARGIN, LATIDX, XIN    )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZRYTPDT' )
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


C        Compute the Z coordinate of the apex of the latitude cone.
C
         CALL ZZELNAXX ( RE, RP, MINLAT, XINCPT, YINCPT )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         

         APEX(1) = 0.D0
         APEX(2) = 0.D0
         APEX(3) = YINCPT

         CALL INCNSG ( APEX,   Z,  ANGLE, VERTEX,
     .                 ENDPT2, NX, SRFX,  XPT2    )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         

C
C        Find the offset of the ray's vertex from the cone's apex,
C        and find the angular separation of the offset from the +Z
C        axis. This separation enables us to compare the latitude of
C        the vertex to the latitude boundary without making a RECGEO
C        call to compute the planetodetic coordinates of the vertex.
C
C        (The comparison will be done later.)
C
         CALL VSUB ( VERTEX, APEX, VTXOFF )

         VTXANG = VSEP ( VTXOFF, Z )

C
C        Check whether the latitude of the intercept can be 
C        considered to match that of the cone.
C
         IF ( NX .GT. 0 ) THEN            
C
C           Check the first intercept.
C
            XVAL1 =  ZZPDPLTC( RE, F, SRFX, MINLAT )
            XVAL2 =  .FALSE.
         
            IF ( NX .EQ. 2 ) THEN
C
C              Check the second intercept.
C
               XVAL2 =  ZZPDPLTC( RE, F, XPT2, MINLAT )

            END IF

            IF ( XVAL1 .AND. ( .NOT. XVAL2 ) ) THEN

               NX = 1

            ELSE IF (  XVAL2  .AND.  (.NOT. XVAL1 )  ) THEN
C
C              Only the second solution is valid. Overwrite
C              the first.
C
               NX = 1

               CALL VEQU( XPT2, SRFX )

            ELSE IF (  ( .NOT. XVAL1 )  .AND.  ( .NOT. XVAL2 )  ) THEN
C
C              Neither solution is valid.
C
               NX = 0

            END IF
              
         END IF            
  

         IF ( NX .GE. 1 ) THEN
C
C           The ray intercept SRFX lies on the lower latitude boundary.
C
C           See whether SRFX meets the longitude and proxy altitude
C           constraints.
C
            CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, LATIDX, XIN )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZRYTPDT' )
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

                  IF ( VTXANG .GT. ANGLE ) THEN

                     IF (      ( MINLAT    .GT. 0.D0 )
     .                    .OR. ( VERTEX(3) .LT. 0.D0 )  ) THEN
C
C                       If MINLAT is positive, the vertex offset
C                       being outside the cone is enough to 
C                       guarantee the planetodetic latitude of the
C                       vertex is less than that of the cone.
C
C                       If MINLAT is non-positive, the angle of the
C                       vertex offset relative to the +Z axis is not
C                       enough; we need the vertex to lie below the 
C                       X-Y plane as well. 
C
C                       Getting here means one of these conditions
C                       was met.
C
C                       Since the latitude of the vertex is less than
C                       than MINLAT, this is the best solution, since
C                       the volume element is on the other side of the
C                       minimum latitude boundary.
C
                        CALL CHKOUT ( 'ZZRYTPDT' )
                        RETURN

                     END IF

                  END IF
C
C                 This is the best solution seen so far, but we 
C                 need to check the remaining boundaries.

                  MNDIST = DIST

               END IF

            END IF

            IF ( NX .EQ. 2 ) THEN
C
C              Check the second solution as well.
C
               CALL ZZINPDT ( XPT2,   BOUNDS, CORPAR, 
     .                        MARGIN, LATIDX, XIN    )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZRYTPDT' )
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
                     CALL CHKOUT ( 'ZZRYTPDT' )
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
            CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, LONIDX, XIN )
            
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZRYTPDT' )
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
         CALL VPACK ( -SIN(MAXLON), COS(MAXLON), 0.D0, EASTB )

         CALL ZZINRYPL ( VERTEX, UDIR, EASTB, 0.D0, S, NX, SRFX )

         IF ( NX .EQ. 1 ) THEN
C
C           We have one point of intersection. Determine whether it's a
C           candidate solution.
C
            CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, LONIDX, XIN )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZRYTPDT' )
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
C     Find the intersection of the ray and lower bounding
C     ellipsoid, if possible.  
C
      CALL SURFPT ( VERTEX, UDIR, EMIN, EMIN, PMIN, SRFX, FOUND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         


      IF ( FOUND ) THEN
C
C        See whether this solution is in the element.
C
         CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, ALTIDX, XIN )
            
         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         


         IF ( XIN ) THEN

            DIST = VDIST( VERTEX, SRFX )

            IF ( DIST .LT. MNDIST ) THEN
C
C              Record the intercept, distance, and surface index.
C
               CALL VEQU ( SRFX, XPT )

               NXPTS  = 1
               MNDIST = DIST

            END IF

         END IF
 
      END IF

C
C     Unlike the outer ellipsoid, either intersection of the ray with
C     the inner ellipsoid might be a valid solution. We'll test for the
C     case where the intersection farther from the ray's vertex is the
C     correct one.
C
      CALL VMINUS( UDIR, NEGDIR )

      CALL SURFPT ( ENDPT2, NEGDIR, EMIN, EMIN, PMIN, SRFX, FOUND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZRYTPDT' )
         RETURN
      END IF         


      IF ( FOUND ) THEN
               
         CALL ZZINPDT ( SRFX, BOUNDS, CORPAR, MARGIN, ALTIDX, XIN )
            
         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZRYTPDT' )
            RETURN
         END IF         


         IF ( XIN ) THEN

            DIST = VDIST( VERTEX, SRFX )

            IF ( DIST .LT. MNDIST ) THEN
C
C              Record the intercept, distance, and surface index.
C
               CALL VEQU ( SRFX, XPT )

               NXPTS  = 1
C
C              There's no need to update MNDIST at this point.
C
            END IF

         END IF
 
      END IF

C
C     NXPTS and XPT are set.
C
      CALL CHKOUT ( 'ZZRYTPDT' )
      RETURN
      END 
