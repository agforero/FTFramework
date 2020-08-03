C$Procedure ZZLATBOX (Bounding box for latitudinal volume element)
 
      SUBROUTINE ZZLATBOX ( BOUNDS, CENTER, LR, LT, LZ, RADIUS )
 
C$ Abstract
C
C     Create a bounding box for a latitudinal volume element.
C
C     The outputs are the box's center, dimensions in the radial,
C     tangential, and vertical directions, and the box's radius.
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
C     GEOMETRY
C     MATH
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      CENTER ( 3 )
      DOUBLE PRECISION      LR
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LZ
      DOUBLE PRECISION      RADIUS

      DOUBLE PRECISION      MARGIN
      PARAMETER           ( MARGIN = 1.D-12 )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BOUNDS     I   Array of coordinate bounds for segment.
C     CENTER     O   Center of bounding box.
C     LR         O   Extent of box in the cylindrical radial direction.
C     LT         O   Extent of box in the tangential direction.
C     LZ         O   Extent of box in the Z direction.
C     RADIUS     O   Radius of box.
C     MARGIN     P   Angular margin.
C
C$ Detailed_Input
C
C     BOUNDS     is a 2 x 3 double precision array containing bounds
C                on the latitudinal coordinates of the spatial region
C                (volume element) covered by a DSK segment. The
C                contents of BOUNDS are:
C
C                   BOUNDS(*,1)              Longitude bounds
C                   BOUNDS(*,2)              Latitude bounds
C                   BOUNDS(*,3)              Radius bounds
C
C                Elements (1,*) are lower bounds; elements (2,*) are
C                upper bounds. Angular units are radians.
C
C$ Detailed_Output
C
C     CENTER     is a double precision 3-vector representing the center
C                of a box tangent to and containing the volume
C                specified by BOUNDS. The sides of the box are parallel
C                to the radial, tangential, and Z directions, where the
C                radial direction is a vector orthogonal to the Z axis,
C                having longitude equal to the midpoint of the
C                segment's longitude range. The tangential direction is
C                orthogonal to the radial and Z directions.
C
C     LR,
C     LT,
C     LZ         are, respectively, the extents (edge lengths) of the
C                bounding box in the radial, tangential, and Z
C                directions.
C
C     RADIUS     is the radius of the sphere that circumscribes the
C                box. RADIUS is equal to the length of a line segment
C                connecting the center of the box to any corner.
C
C$ Parameters
C
C     MARGIN     is used for latitude range validity tests. Latitude
C                values can be out of range by as much as MARGIN.
C
C$ Exceptions
C
C     1) If the minimum longitude exceeds the maximum by more than
C        2*pi, the error SPICE(BADLONGITUDERANGE) is signaled.
C
C     2) If the latitude bounds are out of order, the error
C        SPICE(BADLATITUDEBOUNDS) is signaled.
C
C     3) If the minimum latitude is less than -pi/2 - MARGIN, or
C        if the maximum latitude is greater than pi/2 + MARGIN, 
C        the error SPICE(BADLATITUDERANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine provides a simple way to compute the center and
C     radius for an outer bounding sphere for the volume covered
C     by a DSK segment using latitudinal coordinates.
C
C$ Examples
C
C     None.
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 03-FEB-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     bounding box for latitudinal volume element
C
C-&

C
C     SPICELIB functions
C     
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      TWOPI
      DOUBLE PRECISION      VNORM

      LOGICAL               RETURN

C
C     Local parameters
C
      
C
C     Local variables
C
      DOUBLE PRECISION      DIAG   ( 3 )
      DOUBLE PRECISION      HDLON
      DOUBLE PRECISION      INRAD
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MAXZ
      DOUBLE PRECISION      MIDLON
      DOUBLE PRECISION      MIDR
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      MINZ
      DOUBLE PRECISION      OUTRAD

C
C     This routine uses discovery check-in. We check RETURN in order to
C     avoid performing math operations using invalid operands.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Get local copies of the bounds of the volume element.
C
      MINLON = BOUNDS(1,1)
      MAXLON = BOUNDS(2,1)

      IF ( MAXLON .LE. MINLON ) THEN
         MAXLON = MAXLON + TWOPI()
      END IF

      IF ( MAXLON .LE. MINLON) THEN

         CALL CHKIN  ( 'ZZLATBOX'                               )
         CALL SETMSG ( 'Longitude bounds are #:#. Minimum '
     .   //            'longitude exceeds maximum by more than '
     .   //            '2 pi.'                                  )
         CALL ERRDP  ( '#', MINLON                              )
         CALL ERRDP  ( '#', BOUNDS(2,1)                         )
         CALL SIGERR ( 'SPICE(BADLONGITUDERANGE)'               )
         CALL CHKOUT ( 'ZZLATBOX'                               )
         RETURN

      END IF

      MINLAT = BOUNDS(1,2)
      MAXLAT = BOUNDS(2,2)

      MINR   = BOUNDS(1,3)
      MAXR   = BOUNDS(2,3)

      IF ( MINLAT .GT. MAXLAT ) THEN

         CALL CHKIN  ( 'ZZLATBOX'                              )
         CALL SETMSG ( 'Latitude bounds #:# are out of order.' )
         CALL ERRDP  ( '#', MINLAT                             )
         CALL ERRDP  ( '#', MAXLAT                             )
         CALL SIGERR ( 'SPICE(BADLATITUDEBOUNDS)'              )
         CALL CHKOUT ( 'ZZLATBOX'                              )
         RETURN

      END IF

      IF ( MINLAT .LT. ( -HALFPI() - MARGIN )  ) THEN

         CALL CHKIN  ( 'ZZLATBOX'                               )
         CALL SETMSG ( 'Minimum latitude # is less than -pi/2.' )
         CALL ERRDP  ( '#', MINLAT                              )
         CALL SIGERR ( 'SPICE(BADLATITUDERANGE)'                )
         CALL CHKOUT ( 'ZZLATBOX'                               )
         RETURN

      END IF

      IF ( MAXLAT .GT. ( HALFPI() + MARGIN )  ) THEN

         CALL CHKIN  ( 'ZZLATBOX'                               )
         CALL SETMSG ( 'Maximum latitude # is more than -pi/2.' )
         CALL ERRDP  ( '#', MAXLAT                              )
         CALL SIGERR ( 'SPICE(BADLATITUDERANGE)'                )
         CALL CHKOUT ( 'ZZLATBOX'                               )
         RETURN

      END IF

      MINLAT = MAX( MINLAT, -HALFPI() )
      MAXLAT = MIN( MAXLAT,  HALFPI() )
      
C
C     Let INRAD and OUTRAD be, respectively, the radii of the
C     orthogonal projections onto the X-Y plane of the element's arcs
C     of minimum and maximum distance from the Z axis.
C
C     If the element lies on or above the X-Y plane, the outer radius
C     is that of the lower latitude bound on the outer bounding sphere,
C     and the inner radius is that of the upper latitude bound on the
C     inner bounding sphere.
C
C     These relationships are reversed for elements that lie on or
C     below the X-Y plane.
C
C     For elements that span the X-Y plane, the outer radius is that of
C     the outer bounding sphere. The inner radius is that of the
C     latitude circle on the inner bounding sphere for which the
C     absolute value of the latitude is maximum.
C
      IF ( MINLAT .GE. 0.D0 ) THEN

         OUTRAD = COS( MINLAT ) * MAXR
         INRAD  = COS( MAXLAT ) * MINR

      ELSE IF ( MAXLAT .LE. 0.D0 ) THEN

         OUTRAD = COS( MAXLAT ) * MAXR
         INRAD  = COS( MINLAT ) * MINR

      ELSE

         OUTRAD = MAXR
         INRAD  = COS(  MAX( ABS(MAXLAT), ABS(MINLAT) )  ) * MINR

      END IF

C
C     Let MIDLON be the longitude of the midpoint of the element's
C     longitude coverage. Let HDLON be one half of the extent of the
C     longitude coverage.
C
      HDLON  = ( MAXLON - MINLON ) / 2
      MIDLON =   MINLON + HDLON
    
C
C     LR is the length of the bounding box in the radial direction,
C     where "radius" is defined as the distance from the Z axis.
C
      IF ( HDLON .LE. HALFPI() ) THEN
 
         LR = OUTRAD - ( INRAD * COS(HDLON) )

      ELSE

         LR = ( 1.D0 - COS(HDLON) ) * OUTRAD

      END IF

C
C     The tangential length of bounding box depends on the longitude
C     extent. For any extent larger than Pi radians, the width
C     is just that of the outer radius.
C
      IF ( HDLON .GT. HALFPI() ) THEN

         LT = 2.D0 * OUTRAD
      ELSE
         LT = 2.D0 * OUTRAD * SIN( HDLON ) 
      END IF

C
C     The height bounds are derived from the lowest and highest points
C     on the element.
C
      IF ( MINLAT .GE. 0.D0 ) THEN
 
         MAXZ = MAXR * SIN(MAXLAT)
         MINZ = MINR * SIN(MINLAT) 

      ELSE IF ( MAXLAT .LE. 0.D0 ) THEN

         MAXZ = MINR * SIN(MAXLAT)
         MINZ = MAXR * SIN(MINLAT)

      ELSE

         MAXZ = MAXR * SIN(MAXLAT)
         MINZ = MAXR * SIN(MINLAT)

      END IF

      LZ = MAXZ - MINZ

C
C     Make sure all dimensions are non-negative.
C
      LR = MAX( 0.D0, LR )
      LT = MAX( 0.D0, LT )
      LZ = MAX( 0.D0, LZ )

C
C     Compute the coordinates of the center of the box.
C
C     The box center lies on the meridian of central
C     longitude. The outer tangential edge is at radius
C     OUTRAD. Let MIDR be the radius of the center.
C
      MIDR = OUTRAD - (LR/2)

      CALL CYLREC ( MIDR, MIDLON, MINZ+(LZ/2), CENTER )
      
C
C     The radius is the distance from the center of the box 
C     to any corner.
C
      CALL VPACK ( LR/2, LT/2, LZ/2, DIAG )
      
      RADIUS = VNORM( DIAG )

      RETURN
      END
