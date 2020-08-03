C$Procedure      INCNSG ( Intersection of cone and line segment )
 
      SUBROUTINE INCNSG ( APEX,   AXIS,  ANGLE, ENDPT1, 
     .                    ENDPT2, NXPTS, XPT1,  XPT2   )
 
C$ Abstract
C
C     Compute the points of intersection of a specified nappe of a cone
C     and a line segment.
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
C     CONE
C     GEOMETRY
C     INTERSECTION
C     LINE
C     MATH
C     SEGMENT
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      APEX   ( 3 )
      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      ANGLE
      DOUBLE PRECISION      ENDPT1 ( 3 )
      DOUBLE PRECISION      ENDPT2 ( 3 )
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT1   ( 3 )
      DOUBLE PRECISION      XPT2   ( 3 )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     APEX       I   Apex of cone.
C     AXIS       I   Axis of cone.
C     ANGLE      I   Angle of cone.
C     ENDPT1,
C     ENDPT2     I   Endpoints of line segment.
C     NXPTS      O   Number of intersection points.
C     XPT1       O   First intersection point, if it exists.
C     XPT2       O   Second intersection point, if it exists.
C
C$ Detailed_Input
C
C     APEX       is the apex (tip) of the cone. In this routine's
C                documentation, we'll consider the cone to be a
C                semi-infinite pyramid with circular cross-section. In
C                some contexts, this object is called one "nappe" of
C                the complete cone.
C
C     AXIS       is an axis vector of the cone.                  
C
C     ANGLE      is the angular separation from AXIS of the rays
C                comprising the cone. Let the notation
C
C                   < A, B >
C 
C                denote the dot product of vectors A and B, and let
C
C                   ||A||
C
C                denote the norm of vector A. Then the cone is the set
C                of points
C
C                             X-APEX       AXIS
C                   { X:  < ----------,  -------- >  =  cos(ANGLE) }
C                           ||X-APEX||   ||AXIS||
C
C
C     ENDPT1,
C     ENDPT2     are endpoints of a line segment. These points 
C                must be distinct.
C          
C$ Detailed_Output
C
C     NXPTS      is the number of points of intersection of the input
C                line segment and cone.
C
C     XPT1       is the point of intersection of the segment and cone
C                that is closest to ENDPT1, if an intersection exists.
C                If there are no intersections, XPT1 is undefined.
C
C     XPT2       is the point of intersection of the segment and cone
C                that is farthest from ENDPT1, if two points of
C                intersection exist. If there are not two
C                intersections, XPT2 is undefined.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If AXIS is the zero vector, the error SPICE(ZEROVECTOR)
C         will be signaled.
C
C     2)  If ANGLE is less than zero, the error SPICE(INVALIDANGLE)
C         will be signaled.
C
C     3)  If ENDPT1 and ENDPT2 coincide, the error
C         SPICE(ENDPOINTSMATCH) will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is used by the SPICELIB DSK subsystem. In
C     particular, it is used to determine whether a ray contacts a
C     latitude boundary of a volume element in either planetocentric
C     latitudinal or planetodetic coordinates.
C
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input
C     (if any), the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C
C     1) Compute the intersection of a line segment and cone in
C        a simple case for which the results can easily be checked.
C     
C        Let the apex of the cone be at the origin. Let the axis
C        of the cone lie on the +X axis. Let the angle of the cone
C        be 45 degrees. Let the line segment have endpoints 
C
C           ENDPT1 = ( 1,   -2, sqrt(3)/2 )
C           ENDPT2 = ( 1,    2, sqrt(3)/2 )
C     
C        We expect there to be two points of intersection:
C      
C           XPT1   = ( 1, -1/2, sqrt(3)/2 )
C           XPT2   = ( 1,  1/2, sqrt(3)/2 )
C
C
C        Example code begins here. 
C 
C
C              PROGRAM EX1
C              IMPLICIT NONE
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      RPD
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         FMT1
C              PARAMETER           ( FMT1 = '(A,3F13.8)' )
C
C              CHARACTER*(*)         FMT2
C              PARAMETER           ( FMT2 = '(A,I2)' )
C        C
C        C     Local variables
C        C
C              DOUBLE PRECISION      ANGLE
C              DOUBLE PRECISION      APEX   ( 3 )
C              DOUBLE PRECISION      AXIS   ( 3 )
C              DOUBLE PRECISION      ENDPT1 ( 3 )
C              DOUBLE PRECISION      ENDPT2 ( 3 )
C              DOUBLE PRECISION      SQ3
C              DOUBLE PRECISION      XPT1   ( 3 )
C              DOUBLE PRECISION      XPT2   ( 3 )
C
C              INTEGER               NXPTS
C
C        C
C        C     Set up the cone's geometric attributes.
C        C
C              CALL VPACK ( 0.D0, 0.D0, 0.D0, APEX )
C              CALL VPACK ( 1.D0, 0.D0, 0.D0, AXIS )
C
C              ANGLE = 45.D0 * RPD()
C        C
C        C     Initialize the line segment's endpoints.
C        C
C              SQ3 = SQRT( 3.D0  )
C
C              CALL VPACK ( 1.D0, -2.D0, SQ3/2, ENDPT1 )
C              CALL VPACK ( 1.D0,  2.D0, SQ3/2, ENDPT2 )
C        C
C        C     Find the points of intersection.
C        C
C              CALL INCNSG ( APEX,   AXIS,  ANGLE, ENDPT1,
C             .              ENDPT2, NXPTS, XPT1,  XPT2   )
C
C              WRITE (*,*) ' '
C              WRITE (*,FMT1) 'Apex:        ', APEX
C              WRITE (*,FMT1) 'Axis:        ', AXIS
C              WRITE (*,FMT1) 'Angle (deg): ', ANGLE/RPD()
C              WRITE (*,FMT1) 'Endpoint 1:  ', ENDPT1
C              WRITE (*,FMT1) 'Endpoint 2:  ', ENDPT2
C              WRITE (*,*) ' '
C              WRITE (*,FMT2) 'Number of intersection points: ',
C             .            NXPTS
C              WRITE (*,*) ' '
C              WRITE (*,FMT1) 'Point 1:    ', XPT1
C              WRITE (*,FMT1) 'Point 2:    ', XPT2
C              WRITE (*,*) ' '
C
C              END
C
C    
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, the output was:
C
C
C        Apex:           0.00000000   0.00000000   0.00000000
C        Axis:           1.00000000   0.00000000   0.00000000
C        Angle (deg):   45.00000000
C        Endpoint 1:     1.00000000  -2.00000000   0.86602540
C        Endpoint 2:     1.00000000   2.00000000   0.86602540
C
C        Number of intersection points:  2
C
C        Point 1:       1.00000000  -0.50000000   0.86602540
C        Point 2:       1.00000000   0.50000000   0.86602540
C
C
C$ Restrictions
C
C     1)  This routine is designed to avoid arithmetic overflow in
C         normal cases, such as those in which the line segment is
C         nearly parallel to the cone. However, it is possible to cause
C         arithmetic overflow by using input vectors with extremely
C         large magnitudes.
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
C-    SPICELIB Version 1.0.0 26-OCT-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     intersection of line segment and cone
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      PI
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VDOT

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C     
      DOUBLE PRECISION      COSTOL
      PARAMETER           ( COSTOL = 1.D-10 )

      DOUBLE PRECISION      HSPTOL
      PARAMETER           ( HSPTOL = 1.D-14 )

      INTEGER               UBPL
      PARAMETER           ( UBPL   = 4 )

C
C     Local variables
C
      DOUBLE PRECISION      A
      DOUBLE PRECISION      AXMAG
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      DOUBLE PRECISION      CA2
      DOUBLE PRECISION      COLAT
      DOUBLE PRECISION      COSANG
      DOUBLE PRECISION      COSERR
      DOUBLE PRECISION      DIR    ( 3 )
      DOUBLE PRECISION      DMAG
      DOUBLE PRECISION      DP1
      DOUBLE PRECISION      DP2
      DOUBLE PRECISION      LOCANG
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXP   ( 3 )
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINP   ( 3 )
      DOUBLE PRECISION      NRMPLN ( UBPL )
      DOUBLE PRECISION      OFF1   ( 3 )
      DOUBLE PRECISION      OFF2   ( 3 )     
      DOUBLE PRECISION      ORIGIN ( 3 )
      DOUBLE PRECISION      PLNX   ( 3 )
      DOUBLE PRECISION      S1
      DOUBLE PRECISION      S2
      DOUBLE PRECISION      UAXIS  ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      UOFF1  ( 3 )
      DOUBLE PRECISION      UOFF2  ( 3 )
      DOUBLE PRECISION      UUAX
      DOUBLE PRECISION      UV1    ( 3 )
      DOUBLE PRECISION      UV2    ( 3 )
      DOUBLE PRECISION      V1     ( 3 )
      DOUBLE PRECISION      V1MAG  
      DOUBLE PRECISION      V2     ( 3 )
      DOUBLE PRECISION      V2MAG  
      DOUBLE PRECISION      VTEMP  ( 3 )
      DOUBLE PRECISION      VTEMP2 ( 3 )
      DOUBLE PRECISION      W2
      DOUBLE PRECISION      WU
      DOUBLE PRECISION      WUAX
      DOUBLE PRECISION      X      ( 3 )
      DOUBLE PRECISION      XOFF1  ( 3 )
      DOUBLE PRECISION      XOFF2  ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )
      DOUBLE PRECISION      Y      ( 3 )
      DOUBLE PRECISION      Z      ( 3 )

      INTEGER               I
      INTEGER               N
      INTEGER               NPLNX

      LOGICAL               IN1
      LOGICAL               IN2
      LOGICAL               ISBRCK
      LOGICAL               NEG1
      LOGICAL               NEG2
 
C
C     Saved values
C
      SAVE                  ORIGIN
      SAVE                  Y
      SAVE                  Z

C
C     Initial values
C
      DATA                  ORIGIN / 3 * 0.D0         /
      DATA                  Y      / 0.D0, 1.D0, 0.D0 /
      DATA                  Z      / 0.D0, 0.D0, 1.D0 /

C
C     Use quasi-discovery check-in. We'll check in before
C     code sections that can generate SPICE errors, and check
C     out afterward. When those code sections are skipped,
C     we avoid traceback participation.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     No intersection was found so far.
C
      NXPTS = 0

C
C     The cone's axis vector must be non-zero.
C
      CALL UNORM ( AXIS, UAXIS, AXMAG )

      IF ( AXMAG .EQ. 0.D0 ) THEN
         
         CALL CHKIN  ( 'INCNSG'                                       )
         CALL SETMSG ( 'The cone''s axis vector must be non-zero '
     .   //            'but sadly, it failed to meet this criterion.' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                            )
         CALL CHKOUT ( 'INCNSG'                                       )
         RETURN

      END IF

C
C     The cone's angular radius must be non-negative.
C
      IF ( ANGLE .LT. 0.D0 ) THEN

         CALL CHKIN  ( 'INCNSG'                             )
         CALL SETMSG ( 'The cone''s angular radius must be '
     .   //            ' non-negative but was # (radians).' )
         CALL ERRDP  ( '#', ANGLE                           )
         CALL SIGERR ( 'SPICE(INVALIDANGLE)'                )
         CALL CHKOUT ( 'INCNSG'                             )
         RETURN

      END IF

C
C     The endpoints of the segment must be distinct. Check this after
C     computing a unit direction vector for the line segment.
C
      CALL VSUB  ( ENDPT2, ENDPT1, DIR )
      
      CALL UNORM ( DIR, UDIR, DMAG )

      IF ( DMAG .EQ. 0.D0 ) THEN
         
         CALL CHKIN  ( 'INCNSG'                              )
         CALL SETMSG ( 'The distance between the segment''s '
     .   //            'endpoints was zero. First endpoint: '
     .   //            '(# # #).'                            )
         CALL ERRDP  ( '#', ENDPT1(1)                        )
         CALL ERRDP  ( '#', ENDPT1(2)                        )
         CALL ERRDP  ( '#', ENDPT1(3)                        )
         CALL SIGERR ( 'SPICE(ENDPOINTSMATCH)'               )
         CALL CHKOUT ( 'INCNSG'                              )
         RETURN

      END IF

C
C     Store the cosine of the cone's angular radius. We'll treat all
C     cases with COSANG equal to 0 as though the cone is actually a
C     plane normal to the axis and containing the apex.
C
      COSANG = COS( ANGLE )
      LOCANG = ANGLE
C
C     We'll work with a local axis that has angular separation of
C     no more than pi/2 from the nappe.
C
      IF ( COSANG .LT. 0.D0 ) THEN

         COSANG   = - COSANG
         LOCANG   =   PI() - ANGLE

         UAXIS(1) = - UAXIS(1)
         UAXIS(2) = - UAXIS(2)
         UAXIS(3) = - UAXIS(3)

      END IF

C
C     Compute the offsets of the endpoints of the segment from
C     the cone's apex.
C
      CALL VSUB ( ENDPT1, APEX, OFF1 )
      CALL VSUB ( ENDPT2, APEX, OFF2 )

C
C     Deal with some of the simple cases first. 
C
      CALL VHAT ( OFF1, UOFF1 )
      CALL VHAT ( OFF2, UOFF2 )

      DP1  = VDOT( UOFF1, UAXIS )
      DP2  = VDOT( UOFF2, UAXIS )

C
C     The given axis is inside the nappe defined by the angular radius.
C
C     There's no intersection if both endpoints are in the interior of
C     the nappe of the cone (since the nappe is convex).
C     
      IN1  = DP1 .GE. COSANG
      IN2  = DP2 .GE. COSANG
C
C     If the line segment lies on the far side of the plane that
C     contains the apex and is orthogonal to the axis, there's no
C     intersection.
C
      NEG1 = DP1 .LT. 0.D0
      NEG2 = DP2 .LT. 0.D0

      IF (      ( IN1  .AND. IN2  ) 
     .     .OR. ( NEG1 .AND. NEG2 ) ) THEN
C
C        The segment is in the interior of the cone or
C        on the far side of the plane.
C  
         NXPTS = 0

         RETURN

      END IF

C
C     Here's where we handle the half-space case.
C
      IF ( ABS(COSANG) .LT. HSPTOL ) THEN
C
C        See whether the ray emanating from the first endpoint and
C        having direction UDIR hits the plane normal to the axis and
C        containing the apex. We'll call this plane NRMPLN.
C
C        NVP2PL can signal an error only if the input axis is the 
C        zero vector. We've ensured that it isn't.
C
         CALL NVP2PL ( UAXIS,  APEX, NRMPLN )
         CALL INRYPL ( ENDPT1, UDIR, NRMPLN, NPLNX, PLNX )

C
C        If the ray doesn't hit the plane, we're done. Otherwise,
C        check the intercept.
C
         IF ( NPLNX .EQ. 1 ) THEN
C
C           The ray does hit the plane. If the intersection is on the
C           line segment, we have a solution.
C         
            IF ( VDIST(PLNX, ENDPT1) .LE. DMAG ) THEN
C
C              The intercept is not further along the ray than the
C              second endpoint. It's a valid solution.
C
               NXPTS = 1
               CALL VEQU ( PLNX, XPT1 )
               
            END IF

         END IF
C
C        This is the end of the half-space case.
C
         RETURN
          
      END IF

C
C     At this point we've disposed of the trivial cases. We'll
C     set up a quadratic equation for the intersection of the 
C     line segment with the surface of the cone's nappe. 
C
C     Due to round-off errors, the solution of the quadratic may
C     either be inaccurate or may not be found at all. We'll 
C     examine the solutions we find and solve the problem by
C     an alternate method if necessary. However, the quadratic
C     method is fast, so we give it priority.
C
C     The equation of a ray starting at ENDPT1 and having unit
C     direction vector UDIR is
C     
C        RAY  = { ENDPT1 + s*UDIR, s >= 0 }                          (1)
C
C     The equation of the nappe of the cone is
C
C        CONE = { X: < X - APEX, UAXIS > = ||X-APEX|| * cos(ANGLE) } (2)
C
C     where ANGLE is the angular radius of the cone and UAXIS is the
C     unit axis vector. Substituting the right hand side expression of
C     (1) for X in equation (2) and squaring both sides yields a
C     quadratic equation for S. We'll derive the coefficients of the
C     equation below.
C
C     Let 
C
C        Q  = X - APEX
C        W  = ENDPT1 - APEX
C        U  = UDIR
C        CA = cos(ANGLE)
C
C     We can translate the cone and ray by -APEX, and (1) and (2)
C     can be re-written as
C
C        RAY  = { W + s*U, s >= 0 }                                  (3)
C         
C        CONE = { Q: < Q, UAXIS > = ||Q|| * cos(ANGLE) }             (4)
C
C
C        Substituting the ray expression for Q, we obtain
C
C           < W + s*U, UAXIS > = ||W+s*U|| * CA                      (5)
C
C        and squaring both sides yields
C
C                      2                                      2   2
C             <W,UAXIS>  + 2*<W,UAXIS>*<U,UAXIS>*s + <U,UAXIS> * s
C
C                    2                2       2
C           = ( ||W||  + 2*<W,U>*s + s  ) * CA                       (6)
C
C
C       Collecting coefficients of powers of s, we have
C
C                         2     2     2
C              ( <U,UAXIS>  - CA ) * s
C
C                                            2
C            + 2 * ( <W,UAXIS>*<U,UAXIS> - CA * <W,U> ) * s
C
C                       2        2    2
C            + <W,UAXIS>  - ||W|| * CA
C
C
C         =  0                                                       (7)
C
C
C      Before continuing, we observe that the only non-unit vector
C      in (7) is W. So the coefficients in (7) have no possibility
C      of overflowing unless the vertex of the ray is very far from
C      the apex of the cone.
C
C      W has been computed above as OFF1.
C     
C
C         [ Consider adding check on OFF1 here. ]
C      
C
C     Intermediate values:
C     
      UUAX = VDOT( UDIR, UAXIS )
      WUAX = VDOT( OFF1, UAXIS )
      WU   = VDOT( OFF1, UDIR  )
      W2   = VDOT( OFF1, OFF1  )
      CA2  = COSANG * COSANG

C
C     Quadratic coefficients:
C
      A    = ( UUAX * UUAX ) - CA2

      B    = 2 * (  ( WUAX * UUAX ) - ( CA2 * WU )  )

      C    = ( WUAX * WUAX ) - ( W2 * CA2 )

C
C     We're not interested in solutions that lie outside
C     of the line segment. The length of the segment is
C     DMAG.
C
C     Solve the equation, using DMAG as an upper bound
C     on the magnitude of the roots.
C
      CALL ZZCNQUAD ( A, B, C, DMAG, N, S1, S2 )

C
C     Compute the possible intersection points and test them
C     to make sure they really are solutions. 
C
      IF ( N .GT. 0 ) THEN
C
C        Start with the solution closest to the ray's vertex.
C        Compute XPT1 and make sure it's on the correct nappe
C        of the cone.
C        
         IF ( S1 .GE. 0.D0 ) THEN

            XPT1(1) = ENDPT1(1) + S1*UDIR(1)
            XPT1(2) = ENDPT1(2) + S1*UDIR(2)
            XPT1(3) = ENDPT1(3) + S1*UDIR(3)

            CALL VSUB  ( XPT1, APEX, V1 )

C
C           See whether V1 is on the cone.
C
            CALL UNORM ( V1, UV1, V1MAG )

            IF ( V1MAG .GT. 0.D0 ) THEN

               COSERR = ABS( VDOT(UV1,UAXIS) - COSANG )
            ELSE 
               COSERR = 0.D0
            END IF
 
            IF (      ( V1MAG  .EQ. 0.D0   ) 
     .           .OR. ( COSERR .LT. COSTOL ) ) THEN
C
C              The root is on the cone (on the apex if V1MAG is zero).
C
C              We accept this root. Update NXPTS. Note that this is
C              not necessarily the final value of NXPTS; that
C              depends on the validity of the second root.

               NXPTS = 1

            END IF

         END IF

         IF ( N .EQ. 2 ) THEN
C
C           Check the second root.
C
            IF ( S2 .GE. 0.D0 ) THEN

               XPT2(1) = ENDPT1(1) + S2*UDIR(1)
               XPT2(2) = ENDPT1(2) + S2*UDIR(2)
               XPT2(3) = ENDPT1(3) + S2*UDIR(3)

               CALL VSUB  ( XPT2, APEX, V2 )

C
C              See whether V2 is on the cone.
C
               CALL UNORM ( V2, UV2, V2MAG )

               IF ( V2MAG .GT. 0.D0 ) THEN

                  COSERR = ABS( VDOT(UV2,UAXIS) - COSANG )
               ELSE 
                  COSERR = 0.D0
               END IF
 
               IF (      ( V2MAG  .EQ. 0.D0   ) 
     .              .OR. ( COSERR .LT. COSTOL ) ) THEN
C
C                 The root is on the cone (on the apex if V2MAG is
C                 zero).
C
C                 We accept this root.
C
                  NXPTS = NXPTS + 1

                  IF ( NXPTS .EQ. 1 ) THEN
C
C                    This is the only valid root; overwrite XPT1.
C
                     CALL VEQU ( XPT2, XPT1 )

                  END IF

               END IF     

            END IF       

         END IF

      END IF

C
C     We're not done yet. If we have fewer roots than we should, we'll
C     need to solve the problem by an alternate method.
C
C     If we have two roots, we're in good shape. Otherwise we must
C     determine how many roots should be found.  
C
      IF ( NXPTS .LT. 2 ) THEN
C
C        We must determine the expected number of roots, and if
C        we didn't come up with them, we must find the roots
C        by an alternate method.
C
C        We'll examine the containment of the endpoints within the
C        cone.
C
C        The case where both endpoints are inside the cone was handled
C        earlier.
C
C        If one endpoint is inside the cone and one is outside,
C        we expect to have one root.
C
         IF (     ( IN1 .AND. ( .NOT. IN2 ) )
     .       .OR. ( IN2 .AND. ( .NOT. IN1 ) )  ) THEN
C
C           There's supposed to be one root. If we found none, find one
C           now.
C
            IF ( NXPTS .EQ. 0 ) THEN
C
C              ZZCXBRUT signals an error if the axis is the zero
C              vector, but not otherwise. We've already ruled out this
C              situation. Therefore, we don't check in before the
C              following call.
C
               CALL ZZCXBRUT ( APEX,   UAXIS,  LOCANG, 
     .                         ENDPT1, ENDPT2, XPT1,  ISBRCK )
               
               IF ( ISBRCK ) THEN
C
C                 As long as the root was bracketed, XPT1 is a
C                 solution.
C
                  NXPTS = 1

               END IF

            END IF


         ELSE

            CALL CHKIN ( 'INCNSG' )
C
C           Both endpoints are outside the cone. We could have zero to
C           two roots. If the minimum angular separation of the segment
C           from the axis is less than ANGLE, we expect to find two
C           roots; if it's equal to ANGLE, we expect to find one, and
C           if it's greater than ANGLE, none.
C   
C           We'll transform OFF1 and OFF2 into a reference frame in
C           which angular separation from the axis is equivalent to 
C           colatitude. Then we'll find the maximum latitude attained
C           on the segment.
C
C           We'll count the roots we find, so we'll start at zero.
C
            NXPTS = 0

            CALL FRAME ( UAXIS, X, Y )

            DO I = 1, 3

               XFORM(1,I) = X(I) 
               XFORM(2,I) = Y(I) 
               XFORM(3,I) = UAXIS(I) 

            END DO

            CALL MXV ( XFORM, OFF1, XOFF1 )
            CALL MXV ( XFORM, OFF2, XOFF2 )
 
            CALL ZZSGLATX ( XOFF1, XOFF2, MINLAT, MINP, MAXLAT, MAXP )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'INCNSG' )
               RETURN
            END IF

C
C           COLAT is the colatitude of the point of maximum latitude.
C
            COLAT = HALFPI() - MAXLAT


            IF ( COLAT .LT. LOCANG ) THEN
C
C              MAXP is inside the cone. There should be an intersection
C              on the segment between XOFF1 and MAXP and another
C              between MAXP and XOFF2.
C
               CALL ZZCXBRUT ( ORIGIN, Z,    LOCANG, 
     .                         XOFF1,  MAXP, VTEMP,  ISBRCK )

               IF ( ISBRCK ) THEN
C
C                 Convert VTEMP to the original frame, then translate
C                 it so that it's represented as an offset from the
C                 origin.
C                  
                  CALL MTXV ( XFORM,  VTEMP, VTEMP2 )
                  CALL VADD ( VTEMP2, APEX,  XPT1   )

                  NXPTS = 1

               END IF

               CALL ZZCXBRUT ( ORIGIN, Z,     LOCANG, 
     .                         MAXP,   XOFF2, VTEMP,  ISBRCK )

               IF ( ISBRCK ) THEN
C
C                 Convert VTEMP to the original frame, then translate
C                 it so that it's represented as an offset from the
C                 origin.
C                  
                  CALL MTXV ( XFORM,  VTEMP, VTEMP2 )
                  CALL VADD ( VTEMP2, APEX,  XPT2   )

                  IF ( NXPTS .EQ. 1 ) THEN
C
C                    Both roots are valid.
C
                     NXPTS = 2

                  ELSE
C
C                    The second root is the only valid root. Move it
C                    into XPT1.
C                    
                     CALL VEQU ( XPT2, XPT1 )
                     
                     NXPTS = 1

                  END IF

               END IF


            ELSE IF ( COLAT .EQ. LOCANG ) THEN
C
C              The root corresponds to a point of tangency of
C              the segment and cone. This occurs at the point
C              having maximum latitude: MAXP.
C               
               CALL VEQU ( MAXP, XPT1 )

               NXPTS = 1

C
C           Note that if COLAT > LOCANG, there are no roots.
C
            END IF


            CALL CHKOUT ( 'INCNSG' )

         END IF
C
C        This is the end of portion of the "brute force" branch in
C        which both endpoints are outside the cone.
C
      END IF

C
C     NXPTS  has been set. 
C
C     If NXPTS is 1, then XPT1 is set.
C
C     If NXPTS is 2, then both XPT1 and XPT2 are set.
C

      RETURN     
      END 




