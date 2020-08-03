C$Procedure ZZSGLATX ( Line segment latitude extent )
 
      SUBROUTINE ZZSGLATX ( P1, P2, MINLAT, MINP, MAXLAT, MAXP )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Find the latitude extent (extrema of latitude) of a line segment.
C     Latitude is defined in the planetocentric sense.
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
C     EXTREMA
C     GEOMETRY
C     LATITUDE
C     LINE
C
C$ Declarations
 
      IMPLICIT NONE
      
      DOUBLE PRECISION      P1     ( 3 )
      DOUBLE PRECISION      P2     ( 3 )
      DOUBLE PRECISION      MINLAT 
      DOUBLE PRECISION      MINP   ( 3 )
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXP   ( 3 )

      INTEGER               UBPL
      PARAMETER           ( UBPL = 4 )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     P1         I   First line segment endpoint.
C     P2         I   Second line segment endpoint.
C     MINLAT     O   Minimum latitude on segment.
C     MINP       O   Point where minimum latitude is attained.
C     MAXLAT     O   Maximum latitude on segment.
C     MAXP       O   Point where maximum latitude is attained.
C     UBPL       P   SPICE plane upper bound.
C
C$ Detailed_Input
C
C     P1,
C     P2         are endpoints of a given line segment. 
C
C$ Detailed_Output
C
C
C     MINLAT     is the minimum planetocentric (latitudinal) latitude
C                that occurs on the line segment defined by the input
C                endpoints. Units are radians.
C
C     MINP       is a point on the line segment at which  the latitude
C                MINLAT is attained. Note that in some cases, the
C                minimum latitude can occur at more than one point.
C
C     MAXLAT     is the maximum planetocentric (latitudinal) latitude
C                that occurs on the line segment defined by the input
C                endpoints. Units are radians.
C   
C     MAXP       is a point on the line segment at which the latitude
C                MAXLAT is attained. Note that in some cases, the
C                maximum latitude can occur at more than one point.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs when this routine attempts to construct a
C         SPICE plane, the error will be diagnosed by routines in the
C         call tree of this routine.
C
C     2)  If an error occurs when this routine attempts to compute a
C         ray-plane intersection, the error will be diagnosed by
C         routines in the call tree of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C
C     In this routine, "latitude" of a point refers to planetocentric
C     latitude: the angle between the X-Y plane and the ray emanating
C     from the origin and passing through the point.
C
C     The term "latitude extents" on a line segment refers to the
C     minimum and maximum values of latitude attained on that segment.
C     The subset of the line segment at which a latitude extremum is
C     attained can be an endpoint, a single interior point of the
C     segment, or the whole segment.
C
C     There are several geometric cases that can occur:
C
C        Endpoints of segment are linearly dependent
C        ===========================================
C
C        The latitude extrema occur at the segment endpoints. Note that
C        in this case the entire segment is tangent to one or both
C        nappes of a vertical cone. The latitude values attained on the
C        segment can be any non-empty subset of
C
C           { minimum latitude, 0, maximum latitude }
C
C
C        Segment lies in X-Y plane
C        =========================
C         
C        The latitude extrema are both zero.
C
C
C        Segment intersects Z-axis in a single non-origin point
C        ======================================================
C        
C        One latitude extremum will occur on the Z-axis; the other
C        extremum will occur at a segment endpoint. This case may be
C        viewed as a degenerate version of the "local extremum exists
C        in segment" case.
C
C
C        Local extremum exists in segment
C        ================================
C
C        If an extremum occurs at a single point T in the interior of
C        the line segment, we call this a local extremum. Presuming the
C        segment does not intersect the Z-axis, then T is a tangent
C        point on a vertical cone C0 having its apex at the origin. All
C        points, excluding the origin, on the nappe of the cone
C        containing T have latitude equal to that of T. 
C
C        Let P0 be the plane that contains the line segment and the
C        origin. If the endpoints of the line segment are linearly
C        independent, P0 is uniquely defined, and P0 is tangent to the
C        cone C0 at T.
C
C        Let N0 be a normal vector of P0. Let P1 be a plane containing
C        N0 and the Z-axis. If N0 is not parallel to the Z-axis, P1 is
C        uniquely defined, and the point T lies in the plane P1
C        containing the Z-axis and N0: the intersection of the segment
C        with P1 is T.
C
C        Three of the cases excluded here
C
C           - Segment endpoints are linearly dependent
C
C           - Segment intersects the Z-axis, not at the origin
C
C           - Normal to P0 is parallel to Z-axis (segment lies in
C             X-Y plane)
C
C        are all discussed above. The remaining case is discussed
C        below.
C
C
C        Local extremum occurs in extension of segment
C        =============================================
C
C        A local extremum may exist on the line containing the segment,
C        outside of the segment. In this case the segment cannot
C        contain a point where a local extremum occurs. The extrema
C        of latitude on the segment occur at the endpoints.
C      
C
C$ Examples
C
C     1)  One application of this routine is to assist in subsetting
C         a type 2 DSK segment; this task requires determination of
C         a set of triangular plates that intersect a given
C         longitude/latitude rectangle. Note that this set cannot be 
C         found by examination of vertex latitudes alone; latitude
C         extents of edges are needed.
C
C     2)  This routine is also used to find the maximum latitude 
C         change that can occur between two points on a given ray.
C
C$ Restrictions
C
C     None.
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
C-    SPICELIB Version 1.0.0, 29-SEP-2016 (NJB)
C
C        Original version 06-DEC-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find planetocentric latitude extent of a line segment
C     find planetocentric latitude range of a line segment
C     find planetocentric latitude extrema of a line segment
C
C-&



C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDOT

      LOGICAL               FAILED
      LOGICAL               OPSGND
      LOGICAL               RETURN
      LOGICAL               VZERO
      
C
C     Local variables
C     
      DOUBLE PRECISION      CREASE ( 3 )
      DOUBLE PRECISION      DIR    ( 3 )
      DOUBLE PRECISION      DP1
      DOUBLE PRECISION      DP2
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LAT1
      DOUBLE PRECISION      LAT2
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      PLANE2 ( UBPL )
      DOUBLE PRECISION      R
      DOUBLE PRECISION      T      ( 3 )
      DOUBLE PRECISION      Z      ( 3 )
      
      INTEGER               NXPTS
C
C     Saved variables
C
      SAVE                  Z

C
C     Initial values
C
      DATA                  Z / 0.D0, 0.D0, 1.D0 /

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSGLATX' )

C
C     Start by computing latitude at the segment's endpoints.
C
      CALL RECLAT ( P1, R, LON, LAT1 )
      CALL RECLAT ( P2, R, LON, LAT2 )

C
C     Initialize the outputs using latitudes of the endpoints.
C     If there are interior extrema, we'll update these outputs
C     as needed.
C
      IF ( LAT1 .LE. LAT2 ) THEN

         MINLAT = LAT1
         MAXLAT = LAT2

         CALL VEQU ( P1, MINP )
         CALL VEQU ( P2, MAXP )

      ELSE

         MINLAT = LAT2
         MAXLAT = LAT1

         CALL VEQU ( P2, MINP )
         CALL VEQU ( P1, MAXP )

      END IF

C
C     We want to work with the plane containing the origin, P1, and P2.
C     We'll call this plane PLANE1. First see whether P1 and P2 are
C     linearly independent.
C
      CALL VCRSS ( P1, P2, NORMAL )

      IF ( VZERO(NORMAL) ) THEN
C
C        We have a special case: P1 and P2 define a line passing
C        through the origin. The latitude extrema lie on the 
C        segment endpoints, and possibly at every point on the 
C        segment. We've already computed the outputs.
C
         CALL CHKOUT ( 'ZZSGLATX' )
         RETURN

      END IF

C
C     At this point we know that NORMAL is non-zero. Convert it
C     to a unit vector.
C
      CALL VHATIP ( NORMAL )

C
C     Let ALPHA be the non-negative angle between PLANE1 and the X-Y
C     plane. Then ALPHA and -ALPHA are, respectively, the maximum and
C     minimum possible latitudes attained on the input segment.
C     However, these values are not necessarily attained on the
C     segment; we'll need to perform further analysis to find out. We
C     don't need to compute ALPHA, but we'll refer to it in the
C     discussion below.
C
C     The next step is to find the normal vector to the plane defined
C     by Z and NORMAL. We'll call this plane PLANE2. This plane might
C     not exist if NORMAL and Z are linearly dependent. If PLANE2
C     does exist, the X-Y plane and PLANE1 intersect in a "crease"
C     that is normal to PLANE2.
C
      CALL VCRSS ( Z, NORMAL, CREASE )

      IF ( VZERO(CREASE) ) THEN
C
C        Z and NORMAL are linearly dependent; PLANE1 coincides (up to
C        round-off error) with the X-Y plane. We've already computed
C        the outputs.
C
         CALL CHKOUT ( 'ZZSGLATX' )
         RETURN

      END IF

C
C     At this point we know CREASE is non-zero. Convert
C     it to a unit vector.
C
      CALL VHATIP ( CREASE )

C
C     By construction, CREASE is orthogonal to NORMAL. PLANE2 
C     cuts PLANE1 in a line L passing through the origin. If
C     the line segment has an interior latitude extremum, 
C     the point T where that extremum is attained lies on L.
C     The segment is tangent at T to a nappe of a cone, centered on
C     the Z-axis and having its apex at the origin, for which 
C     the half-angle is (pi/2)-ALPHA. The point T lies in PLANE2
C     since L is contained in PLANE2. 
C
C     If a single tangent point T exists in the interior of the
C     segment, then the endpoints must be on opposite sides of PLANE2.
C     See whether this is the case.
C
      DP1 = VDOT(P1,CREASE)
      DP2 = VDOT(P2,CREASE)

      IF (  OPSGND( DP1, DP2 )  ) THEN
C
C        The segment crosses PLANE2 at an interior point; this
C        point is where the extremum occurs. Solve for the
C        intersection.
C
C        CREASE is guaranteed to be a unit vector. A zero input
C        vector is the only cause for which NVC2PL will signal
C        an error. Therefore we don't check FAILED after the 
C        following call.
C
         CALL NVC2PL ( CREASE, 0.D0, PLANE2 )

         CALL VSUB ( P2, P1, DIR )

         CALL INRYPL ( P1, DIR, PLANE2, NXPTS, T )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSGLATX' )
            RETURN
         END IF

         IF ( NXPTS .EQ. 1 ) THEN
C
C           This is the normal case: we have one intersection of the
C           segment with PLANE2. Update whichever of the extrema is
C           superseded by the interior value.
C
C           Note that this case can occur when NORMAL is orthogonal to
C           Z, making ALPHA equal to pi/2. The nappes are degenerate in
C           this case, consisting of the positive and negative Z-axes.
C           This degenerate case occurs when the segment intersects the
C           Z-axis in a point other than the origin, and the endpoints
C           are linearly independent.
C
C           This is not a special case computationally.
C
            CALL RECLAT ( T, R, LON, LAT )

            IF ( LAT .GT. MAXLAT ) THEN
               
               MAXLAT = LAT
               CALL VEQU ( T, MAXP )

            ELSE IF ( LAT .LT. MINLAT ) THEN

               MINLAT = LAT
               CALL VEQU ( T, MINP )

            END IF
C
C           There can be only one local extremum, so we're done.
C
         END IF
C
C        If NXPTS is not 1, then even though the endpoints are on
C        opposite sides of PLANE2, either the segment was found to lie
C        in PLANE2 or no intersection was found. This situation must be
C        due to finite precision arithmetic error. We'll make do with
C        the extrema already found.

      END IF
 
C
C     We reach this point if we found a local extremum or if any of the
C     following are true:
C
C        1)  The segment misses PLANE2 altogether, in which case
C            there's no tangency point.
C
C        2)  One endpoint lies on PLANE2 and one endpoint does not.
C
C        3)  Both endpoints lie in PLANE2. Then both endpoints lie
C            in L, so we should have found them to be linearly
C            dependent. This situation must be due to finite precision
C            arithmetic error.
C
C     In all of the numbered cases the extrema occur at the endpoints.
C     and have been found already. In all cases, the outputs are set.
C     
      CALL CHKOUT ( 'ZZSGLATX' )
      RETURN
      END
