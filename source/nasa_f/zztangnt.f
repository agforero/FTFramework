C$Procedure ZZTANGNT ( DSK, find target tangent rays in half-plane )
 
      SUBROUTINE ZZTANGNT ( CURVE,  SRCRAD, SHAPE,  TRGCDE, NSURF,  
     .                      SRFLST, FIXFID, ET,     PLNVEC, AXIS,
     .                      SCHSTP, SOLTOL, RESULT, POINTS        )
  
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Solve for tangent points on a target surface of rays contained
C     within a specified half-plane. The rays may emanate from a
C     specified vertex or may be tangent to a "source" sphere; the
C     first case supports limb finding and the second terminator
C     finding.
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
C     LIMB
C     RAY
C     TANGENT
C     TERMINATOR
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'zzdsk.inc'

      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      INTEGER               CURVE
      DOUBLE PRECISION      SRCRAD
      INTEGER               SHAPE
      INTEGER               TRGCDE
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      INTEGER               FIXFID
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      PLNVEC ( 3 )
      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      SCHSTP
      DOUBLE PRECISION      SOLTOL
      DOUBLE PRECISION      RESULT ( LBCELL : * )
      DOUBLE PRECISION      POINTS ( 3, * )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     CURVE      I   Type of curve: limb or terminator type.
C     SRCRAD     I   Radius of illumination source.
C     SHAPE      I   Target shape.
C     TRGCDE     I   Target body ID code.
C     NSURF      I   Number of surfaces in list.
C     SRFLST     I   Surface ID list.
C     FIXFID     I   Frame ID of target body-fixed frame.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     PLNVEC     I   Reference vector contained in cutting half-plane.
C     AXIS       I   Axis vector: edge of cutting half-plane.
C     SCHSTP     I   Angular step size for searching.
C     SOLTOL     I   Solution convergence tolerance.
C     RESULT     O   Cell (not window) containing transition angles.
C     POINTS     O   Array of tangent points.
C     ELLSHP     P   Ellipsoid shape code.
C     DSKSHP     P   DSK shape code.
C     LMBCRV     P   Limb code.
C     UMBRAL     P   Umbral terminator code.
C     PNMBRL     P   Penumbral terminator code.
C
C$ Detailed_Input
C
C     CURVE      is an integer code indicating the type of set on the
C                target body on which tangent points are to be found.
C                When the target is an ellipsoid, the set is literally
C                a curve. When the target is a DSK shape, the set is
C                often not connected; so it is not actually a curve;
C                still, we retain the familiar terminology for the
C                ellipsoid case.
C
C                Possible values and meanings are:
C
C                   LMBCRV            limb
C                   UMBRAL            umbral terminator
C                   PNMBRL            penumbral terminator      
C
C                Terminator computations are performed assuming
C                the light source is spherical.
C
C
C     SRCRAD     is the radius of the illumination source. This
C                value is used for terminator computations only.
C                When used, SRCRAD must be strictly positive.
C                Units are km.
C
c
C     SHAPE      is an integer code indicating the target body shape.
C                
C                Possible values and meanings are:
C
C                   ELLSHP            shape is modeled as an ellipsoid
C                   DSKSHP            shape is provided by DSK data
C
C
C     TRGCDE     is the body ID code of the target body.
C
c
C     NSURF,
C     SRFLST     are, respectively, the count of surface IDs and the
C                surface ID list.
C
c
C     FIXFID     is the frame ID code of a body-fixed frame centered
C                on the target body. The output tangent points will
C                be expressed in this frame.
C
C
C     ET         is the computation epoch, expressed as seconds past
C                J2000 TDB. ET is used for DSK segment selection.
C                If the target shape is modeled as an ellipsoid, ET
C                is ignored.
C
C
C     PLNVEC     is a vector used to define a half-plane in which
C                to search for tangent rays. The half-plane contains
C                PLNVEC, the target body center, and AXIS.
C
C                For limb and umbral terminator computations, 
C                tangent rays lie in the half-plane containing PLNVEC.
C        
C                For penumbral terminator computations, tangent rays
C                touch the target in the half-plane containing PLNVEC,
C                but touch the light source in the complementary
C                half-plane bounded by the line containing AXIS.
C     
C
C     AXIS       is a second vector used to define a half-plane in 
C                which to search for tangent vectors. AXIS lies in
C                the line forming the boundary of the half-plane.
C
C                For limb computations, AXIS points from the target
C                toward the observation point.
C
C                For terminator computations, AXIS points from the
C                target toward the center of the illumination source.
C
C                
C     SCHSTP,
C     SOLTOL     are, respectively, the search angular step size and
C                solution convergence tolerance used to find tangent
C                rays and associated limb points within each cutting
C                half plane.
C
C$ Detailed_Output
C
C     RESULT     is a cell (not a window) containing ray-axis angular
C                separation values at which ray-surface tangencies
C                occur.
C
C                Here "angular separation" refers to the angle between
C                a ray and the input vector AXIS. AXIS points away from
C                the target.
C
C                The vertex of the vector depends on the curve type,
C                which is indicated by the input code CURVE.
C
C                For a limb computation, the vertex is the observer's
C                location.
C
C                For an umbral terminator computation, the vertex is on
C                the surface of the light source, in the half- plane
C                defined by PLNDEF and AXIS. The line containing the
C                ray is tangent to the light source at the vertex.
C
C                For a penumbral terminator computation, the vertex is
C                on the surface of the light source, in the half-plane
C                complementary to that defined by PLNDEF and AXIS. The
C                line containing the ray is tangent to the light source
C                at the vertex.
C
C                Units are radians.
C
C
C     POINTS     is an array of tangent points on the target body's
C                surface. Elements
C
C                   POINTS(J,I), J = 1, 3
C
C                constitute the tangent point for the ray having
C                angular separation
C
C                   RESULT(I)
C
C                from AXIS. Here the ray's vertex is as described 
C                above in the discussion of the RESULT output cell.
C                
C
C$ Parameters
C
C     See the INCLUDE file zzdsk.inc for declarations of these 
C     parameters.
C
C
C     ELLSHP     is a code specifying an ellipsoidal shape model.
C
C     DSKSHP     is a code specifying a DSK-based shape model.
C
C     LMBCRV     is a code specifying a limb computation.
C
C     UMBRAL     is a code specifying an umbral terminator computation.
C               
C                The umbral terminator is the boundary of the portion
C                of the target surface that receives no light from the
C                illumination source.
C
C     PNMBRL     is a code specifying an penumbral terminator
C                computation.
C
C                The penumbral terminator is the boundary of the
C                portion of the target surface that is not subject to
C                self-occultation of light from the illumination
C                source. Given that the light source is modeled as
C                a sphere, from any target surface point nearer to
C                the source than the penumbral terminator, the source
C                appears to be a lit disc. 
C                
C
C$ Exceptions
C
C     1)  If AXIS is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C        
C     2)  If PLNDEF is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C
C     3)  If the curve type code is unrecognized, the error
C         SPICE(BADCURVETYPE) is signaled.
C
C     4)  If a terminator curve is specified by the source radius
C         is non-positive, the error SPICE(BADSOURCERADIUS) is
C         signaled.
C
C     5)  If AXIS and PLNDEF are linearly dependent, the error
C         SPICE(DEGENERATECASE) is signaled.
C
C     6)  If the target shape code is unrecognized, the error
C         SPICE(BADSHAPE) is signaled.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference
C     ellipsoids.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine computes target tangent points for limb or
C     terminator computations.
C
C$ Examples
C
C     See usage in LIMBPT and TERMPT.
C
C$ Restrictions
C
C     1) This is a private routine. 
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
C-     SPICELIB Version 1.0.0 24-AUG-2016 (NJB)
C
C           10-FEB-2016 (NJB)
C
C           Removed VERTEX from argument list.
C
C           Original version 26-OCT-2015 (NJB)
C-&
 
C$ Index_Entries
C
C     find target tangent rays in specified half-plane
C
C-&
 

C
C     SPICELIB functions
C
      DOUBLE PRECISION      DASINE
      DOUBLE PRECISION      PI
      DOUBLE PRECISION      VNORM

      INTEGER               CARDD

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     EXTERNAL routines
C
      EXTERNAL              ZZTANSTA
      EXTERNAL              GFSTEP
      EXTERNAL              GFREFN
      EXTERNAL              GFBAIL
      EXTERNAL              GFREPU

C
C     Local parameters
C
      DOUBLE PRECISION      MARGIN
      PARAMETER           ( MARGIN = 1.D-12 )

C
C     Local variables
C     
      DOUBLE PRECISION      ALPHA
      DOUBLE PRECISION      START
      DOUBLE PRECISION      FINISH
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      MINDST
      DOUBLE PRECISION      R
      DOUBLE PRECISION      REFVEC ( 3 )

      INTEGER               I
      INTEGER               N
      
      LOGICAL               CSTEP
      LOGICAL               ENDFLG ( 2 )



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZTANGNT' )

C
C     Empty the result window.
C
      CALL SCARDD ( 0, RESULT )

C
C     Rotate the plane definition vector by pi about AXIS if 
C     we're generating penumbral terminator points.
C
      IF ( CURVE .EQ. PNMBRL ) THEN

         CALL VROTV ( PLNVEC, AXIS, PI(), REFVEC )
      ELSE
         CALL VEQU  ( PLNVEC,             REFVEC )
      END IF

C
C     Prepare the tangent finding utilities.
C
      CALL ZZTANINI ( CURVE,  SRCRAD, SHAPE,  TRGCDE, NSURF, SRFLST, 
     .                FIXFID, ET,     REFVEC, AXIS                  )
C
C     Fetch a maximum bounding radius for the target.
C
C     Caution: we assume ZZTANINI has initialized the 
C     SINCPT utility subsystem by calling one of
C
C        ZZSUELIN
C        ZZSUDSKI
C
      CALL ZZMAXRAD ( MAXRAD )

C
C     Scale up MAXRAD slightly to ensure bracketing.
C
      MAXRAD = MAXRAD * 1.001D0

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZTANGNT' )
         RETURN
      END IF

      IF ( MAXRAD .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Target maximum radius # is non-positive.' )
         CALL ERRDP  ( '#',  MAXRAD                               )
         CALL SIGERR ( 'SPICE(INVALIDRADIUS)'                     )
         CALL CHKOUT ( 'ZZTANGNT'                                 )
         RETURN

      END IF



      IF ( CURVE .EQ. LMBCRV ) THEN
C
C        We're looking for limb points.
C
C        Set the initial ray-axis separation.
C
         START  = 0.D0
C
C        If the vertex is outside of the bounding sphere,
C        set the initial angle to the supplement of 
C        the angular radius of the target, based on
C        its maximum radius.
C 
         R      = VNORM( AXIS ) 
         MINDST = MAXRAD * ( 1.D0 + MARGIN )

         IF ( R .GT. MINDST ) THEN

            START = PI() -  DASINE( MAXRAD/R, MARGIN ) 

        END IF

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZTANGNT' )
            RETURN
         END IF
C
C        Set the final ray-axis separation.
C
         FINISH = PI()


      ELSE
C
C        For the terminator cases, check for an invalid source radius.
C
         IF ( SRCRAD .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Source radius # is non-positive.' )
            CALL ERRDP  ( '#',  SRCRAD                       )
            CALL SIGERR ( 'SPICE(INVALIDRADIUS)'             )
            CALL CHKOUT ( 'ZZTANGNT'                         )
            RETURN

         END IF
      
C
C        Make sure the source and outer bounding sphere of the
C        target don't intersect.
C
         R = VNORM( AXIS ) 

         IF ( (SRCRAD + MAXRAD) .GT. R ) THEN

            CALL SETMSG ( 'Source radius # and target maximum '
     .      //            'radius # sum to #; distance between '
     .      //            'source and target centers is #. Source ' 
     .      //            'and target are too close together.'     )
            CALL ERRDP  ( '#',  SRCRAD                             )
            CALL ERRDP  ( '#',  MAXRAD                             )
            CALL ERRDP  ( '#',  R                                  )
            CALL SIGERR ( 'SPICE(OBJECTSTOOCLOSE)'                 )
            CALL CHKOUT ( 'ZZTANGNT'                               )
            RETURN

         END IF


         IF ( CURVE .EQ. UMBRAL ) THEN
C
C           We'll search for a point on the umbral terminator.
C
C           For this search, the angle we measure is that between a ray
C           tangent to the source and a ray emanating from the tangent
C           point and parallel to the axis, pointing in the
C           target-source direction. The ray lies in a plane containing
C           the source-target axis.
C
C           The minimum angle is achieved when the ray is tangent to
C           the target sphere; the maximum angle is achieved when the
C           ray intersects the target center.
C
C
C           The following equation is valid regardless of whether
C           or not SRCRAD > MAXRAD.
C
            START = PI() + DASINE( (SRCRAD-MAXRAD)/R, MARGIN ) 
               
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZTANGNT' )
               RETURN
            END IF

C
C           Set the final ray-axis separation.
C
            FINISH = PI() +  DASINE( SRCRAD/R, MARGIN )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZTANGNT' )
               RETURN
            END IF


         ELSE IF ( CURVE .EQ. PNMBRL ) THEN
C
C           We'll search for a point on the umbral terminator.
C
C           We measure the ray's angle from the axis, but in this case,
C           the angle increases in the clockwise (negative sense about
C           the normal to the cutting half-plane defined by the cross
C           product of the axis and the reference vector. Each ray
C           emanating from, and tangent to, the source's surface passes
C           through the axis at a point between the source and target.
C
C           In order to use the root-finding utilities, we treat the
C           angle as a positive quantity.
C
C           The initial ray contains a line segment connecting
C           tangency points on each sphere. The segment, axis,
C           and radii of the spheres form two similar triangles.
C           Below, ALPHA is the fraction of the source-target
C           distance belonging to the triangle having a vertex
C           at the center of the source.
C
            ALPHA = SRCRAD / ( SRCRAD + MAXRAD )

            START = PI() -  DASINE( SRCRAD / (ALPHA*R), MARGIN ) 

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZTANGNT' )
               RETURN
            END IF
                
C
C           We stop looking when the ray intersects the target center.
C
            FINISH = PI() - DASINE( SRCRAD / R, MARGIN ) 

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZTANGNT' )
               RETURN
            END IF

         ELSE

            CALL SETMSG ( 'Input curve code # was not recognized.' )
            CALL ERRINT ( '#', CURVE                               )
            CALL SIGERR ( 'SPICE(BUG)'                             )
            CALL CHKOUT ( 'ZZTANGNT'                               )
            RETURN

         END IF

      END IF

C
C     Search for ray occultations. The endpoints of the occultation
C     intervals are angles at which tangency occurs.
C
C     We consider the angle to be measured from the AXIS vector.
C     The initial and final values START and FINISH have been set
C     above.
C
      CSTEP  = .TRUE.

      CALL ZZTANSLV ( ZZTANSTA, GFSTEP, GFREFN, 
     .                CSTEP,    SCHSTP, START,  FINISH,
     .                SOLTOL,   RESULT, POINTS, ENDFLG )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZTANGNT' )
         RETURN
      END IF

C
C     If the first endpoint of RESULT is the interval start but is not
C     a point of transition, delete it from RESULT. Note that RESULT
C     becomes a cell rather than a window. We must delete the
C     corresponding point from the POINTS array as well.
C
      IF ( CARDD(RESULT) .GT. 0 ) THEN

         IF ( ( RESULT(1) .EQ. START ) .AND. ( .NOT. ENDFLG(1) ) ) THEN

            N = CARDD( RESULT ) 

            DO I = 2, N

               RESULT(I-1) = RESULT(I)

               CALL VEQU ( POINTS(1,I), POINTS(1,I-1) )

            END DO

            CALL SCARDD ( N-1, RESULT )

         END IF

      END IF

C
C     If the final endpoint of RESULT is not a transition, delete
C     it as well. In this case decrementing the cardinality of
C     RESULT suffices.
C
      N = CARDD(RESULT)

      IF ( N .GT. 0 ) THEN

         IF ( ( RESULT(N) .EQ. FINISH ) .AND. ( .NOT. ENDFLG(2) ) ) THEN

            CALL SCARDD ( N-1, RESULT )

         END IF

      END IF

      CALL CHKOUT ( 'ZZTANGNT' )
      RETURN
      END

