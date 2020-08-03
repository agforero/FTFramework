C$Procedure  XDDA  ( list voxels intersected by a ray )
 
      SUBROUTINE XDDA ( VERTEX, RAYDIR, GRDEXT, MAXNVX, NVX, VOXLST )
 
C$ Abstract
C
C     Given a ray and a voxel grid, return a list of voxels the ray
C     intersects.
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
C     GRID
C     INTERSECTION
C     PLATE
C     VOXEL
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      GRDTOL
      PARAMETER           ( GRDTOL = 1.D-12 )

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      INTEGER               GRDEXT ( 3 )
      INTEGER               MAXNVX
      INTEGER               NVX
      INTEGER               VOXLST ( 3, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     GRDTOL     P   Tolerance for vertex distance from grid.
C     VERTEX     I   Voxel grid coordinates of ray's vertex.
C     RAYDIR     I   Direction vector of ray.
C     GRDEXT     I   Dimensions of grid in voxel units.
C     MAXNVX     I   Maximum value of VOXLST.
C     NVX        O   Number of voxels in the VOXLST list.
C     VOXLST     O   List of voxels intersected by ray.
C
C$ Detailed_Input
C
C     VERTEX     Voxel grid coordinates of ray's vertex. These
C                coordinates are zero-based, double precision
C                offsets from the grid's origin. The units of
C                the coordinates are voxels, that is, 
C                voxel edge lengths.
C
C     RAYDIR     Direction vector of ray from VERTEX.
C
C     GRDEXT     The integer 3-vector containing the voxel grid
C                extents. These are the dimensions of the voxel grid in
C                voxel units, in the X, Y, and Z directions
C                respectively.
C
C     MAXNVX     The maximum number of voxel coordinate sets that
C                can be stored in VOXLST.
C
C$ Detailed_Output
C
C     NVX        Number of voxel coordinate sets contained in VOXLST.
C
C     VOXLST     List of coordinate sets of voxels intersected by ray.
C                Elements 
C
C                   VOXLST(J,I), J = 1, 3
C
C                are the coordinates of the Ith voxel in the list.
C                These coordinates are 1-based integer values.
C
C                The voxels in the output list are ordered by
C                increasing distance from the ray's vertex.
C
C$ Parameters
C
C     GRDTOL     is a tolerance value used to determine whether
C                VERTEX is too far from the voxel grid. The Ith
C                component of VERTEX must not differ from the 
C                Ith coordinate of the nearest grid point by more
C                than
C
C                    GRDTOL * EXTENT(I)
C               
C$ Exceptions
C
C     1)  The error SPICE(ZEROVECTOR) is signaled if the input RAYDIR
C         has all zero components.
C
C     2)  The error SPICE(INVALIDSIZE) is signaled if the maximum
C         output list size MAXNVX is non-positive.
C
C     3)  The error SPICE(BADDIMENSIONS) is signaled if any element of
C         the grid extents array GRDEXT is non-positive.
C
C     4)  The error SPICE(VERTEXNOTINGRID) is signaled if the ray's
C         vertex is inside, or within a small distance from, the voxel 
C         grid. See the description of the parameter GRDTOL.
C
C     5)  The error SPICE(ARRAYTOOSMALL) is signaled if the value of
C         the NVX counter (number of intersected voxels) exceeds the
C         size of the VOXLST input vector.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports use of a spatial index for rapid
C     selection of plates that could be hit by a specified ray.
C
C$ Examples
C
C     See the routine DSKX02 for a usage example.
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
C     J.A. Bytof      (JPL)
C
C$ Version
C
C-    SPICELIB Version 3.1.0, 02-FEB-2016 (NJB)
C
C        Updated to call ZZINGRD rather than INGRD.
C        Minor updates were made to header I/O sections.
C
C-    DSKLIB Version 3.0.0, 11-JUL-2014 (NJB)
C
C        Bug fix: a correction was made to the computation of
C        the vertex offset from the bounding planes of the
C        voxel containing the vertex.
C
C        Minor edits were made to comments.
C
C     Last update was 05-JUN-2014 (NJB)
C
C        Bug fix: the use of the MOD function led to a 1-voxel
C        size error when the input ray's vertex was on the 
C        voxel grid boundary. 
C
C        An error check for invalid grid dimensions was added.
C
C        Code to prevent arithmetic overflow was added.
C
C        Code was added to prevent the values AX2ERR and AX3ERR from
C        ever becoming negative when the components of the ray's
C        direction vector in the corresponding directions are zero or
C        too small for a voxel step in those directions to occur.
C     
C        Renamed the routine's arguments, except for NVX.
C
C        Detailed output descriptions were updated to refer to
C        voxel coordinates rather than IDs. References to sorting
C        were deleted.
C
C        In-line comments now explain the routine's algorithm.
C        Old comments that are no longer applicable were deleted.
C
C-    DSKLIB Version 2.1.0, 26-JUL-2010 (NJB)
C
C        Bug fix: voxel space coordinates of input 
C        vertex are now bracketed within the voxel
C        grid.
C
C        This prevents round-off errors from occurring
C        when the vertex is slightly outside the grid,
C        but may not be appropriate for all applications.
C        Therefore it may make sense to make this a 
C        private routine.
C
C-    DSKLIB Version 2.0.0, 20-APR-2010 (NJB)
C
C        Removed commented out lines declaring and calling VOX2ID.
C
C-    DSKLIB Version 1.1.0, 08-OCT-2009 (NJB)
C
C        Updated header.
C
C        Bug fix: driving axis for intercept computation is
C        now determined by largest component of ray direction vector.
C        This fix was made long before this header update.
C
C-    DSKLIB Version 1.1.0, 19-OCT-2004 (EDW)
C
C        Added logic to remove duplicate voxel IDs from
C        the return list. Extended programing comments.
C
C-    DSKLIB Version 1.0.1, 26-AUG-2002 (BVS)
C
C        Replaced WRITE with normal error reporting calls.
C
C-    DSKLIB Version 1.0.0, 03-FEB-1999 (JAB)
C
C-&
 
C$ Index_Entries
C
C     list voxels intersected by a ray
C
C-&
 
C
C     SPICELIB functions
C 
      DOUBLE PRECISION      BRCKTD
      DOUBLE PRECISION      DPMAX

      INTEGER               BRCKTI

      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local parameters 
C
      DOUBLE PRECISION      EPS
      PARAMETER           ( EPS = 1.D-20 )

C
C     Local variables
C
      DOUBLE PRECISION      AX2ERR
      DOUBLE PRECISION      AX3ERR
      DOUBLE PRECISION      LIMIT
      DOUBLE PRECISION      MAXCMP
      DOUBLE PRECISION      VTXOFF ( 3 )
      DOUBLE PRECISION      S12
      DOUBLE PRECISION      S13

      INTEGER               I
      INTEGER               IAXIS  ( 3 )
      INTEGER               ICOORD ( 3 )
      INTEGER               INTVTX ( 3 )
      INTEGER               NEXT   ( 3 )
      INTEGER               STEP   ( 3 )

      LOGICAL               ZZINGRD

C
C     Saved variables
C
      SAVE                  NEXT

C
C     Initial values
C
      DATA                  NEXT   /  2,   3,   1   /
 
C
C     Use discovery check-in.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

C
C     The algorithm below efficiently determines the set of voxels
C     intersected by the input ray. 
C
C     This algorithm doesn't compute the intersections of the ray
C     with the boundaries of the voxels, nor does it ever compute
C     the coordinates of any point on the ray. Instead, it keeps
C     track of the voxel boundary planes that the ray passes through.
C
C     The algorithm starts out by determining which voxel contains the
C     ray's vertex. It computes the distances from the vertex to the
C     "next" voxel boundary planes---those that the ray is headed
C     towards. It maintains measurements that enable it to determine
C     which boundary plane is hit next. The voxel on the other side of
C     that intersection point is the next voxel the ray goes through.
C     In the case of ties, any of the candidate "next" voxels may be
C     selected. The "next" voxel is added to the output voxel list, the
C     measurements of relative distances to the next boundaries are
C     updated, and the algorithm continues in this fashion until a
C     voxel outside the grid is detected.
C
C     The relative distance measurements from the ray's vertex to
C     the "next" boundary planes are defined as follows:
C
C        -  For the primary ray direction---this is the direction
C           corresponding to the component of largest magnitude of the
C           ray's direction vector---the distance is just the
C           difference of the primary coordinates of the next plane and
C           of the ray's vertex.
C
C        -  For each axis orthogonal to the primary one, the algorithm
C           computes the length of the projection onto the primary axis
C           of the portion of the ray extending from the vertex to the
C           "next" voxel boundary plane orthogonal to that non-primary
C           axis. From that projection length the distance from the
C           vertex to the boundary in the primary direction is
C           subtracted.
C
C           For the non-primary axes, these differences are stored in
C           the respective variables
C
C              AX2ERR
C              AX3ERR
C
C           When AX2ERR is negative, the ray will hit the next voxel
C           boundary orthogonal to the "second" axis (having its 
C           index stored in the variable IAXIS(2)) before it hits
C           the next boundary orthogonal to the primary axis. The
C           quantity AX3ERR behaves similarly. 
C
C           If both AX2ERR and AX3ERR are negative, the more negative
C           value marks the boundary plane that is hit first.
C
C     The axes will be re-labeled using the variable IAXIS. IAXIS(1)
C     will be the index of the primary axis.
C
C     There are a few numeric issues to consider:
C
C        1)  The ratios of the components of the ray's direction vector
C            are computed and stored in the variables S12 and S13. Very
C            small components acting as denominators could cause
C            arithmetic overflow.
C
C        2)  The quantities S12 and S13, while representable as double
C            precision numbers, can be quite large. These quantities 
C            may be added repeatedly to the terms AX2ERR and AX3ERR,
C            respectively. These additions could potentially result
C            in arithmetic overflow.
C
C     Both of these problems are addressed by the following observation:
C
C        If a component of the ray direction vector is small enough, and
C        the corresponding component of the ray's vertex is not on a 
C        voxel boundary, the ray will exit the grid before reaching a 
C        bounding plane orthogonal to that component of the direction
C        vector.
C
C        If the above situation holds, but the ray's vertex is already
C        on a boundary plane orthogonal to the small component, then
C        the ray will exit the grid before hitting a parallel boundary
C        plane. 
C
C     So we can safely treat very small direction components as zero.
C
 
C
C     Check if ray direction vector is a zero vector.
C
      IF ( VZERO( RAYDIR ) ) THEN

         CALL CHKIN  ( 'XDDA'                     )
         CALL SETMSG ( 'Ray is the zero vector.'  )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'        )
         CALL CHKOUT ( 'XDDA'                     )
         RETURN

      END IF
 
C
C     Check the voxel grid dimensions.
C
      IF ( MIN( GRDEXT(1), GRDEXT(2), GRDEXT(3) ) .LT. 1 ) THEN

         CALL CHKIN  ( 'XDDA'                                    )
         CALL SETMSG ( 'Voxel grid dimensions must be strictly ' 
     .   //            'positive but are # # #.'                 )
         CALL ERRINT ( '#', GRDEXT(1)                            )
         CALL ERRINT ( '#', GRDEXT(2)                            )
         CALL ERRINT ( '#', GRDEXT(3)                            )
         CALL SIGERR ( 'SPICE(BADDIMENSIONS)'                    )
         CALL CHKOUT ( 'XDDA'                                    )
         RETURN

      END IF

C
C     Make sure the vertex is not too far from the voxel grid.
C
      DO I = 1, 3

         IF (     ( VERTEX(I) .LT.     -GRDTOL  * GRDEXT(I) ) 
     .       .OR. ( VERTEX(I) .GT. (1 + GRDTOL) * GRDEXT(I) )  ) THEN

            CALL CHKIN  ( 'XDDA'                                  )
            CALL SETMSG ( 'Vertex # # # is outside of voxel grid '
     .      //            'defined by extents # # #.'             )
            CALL ERRDP  ( '#', VERTEX(1)                          )
            CALL ERRDP  ( '#', VERTEX(2)                          )
            CALL ERRDP  ( '#', VERTEX(3)                          )
            CALL ERRINT ( '#', GRDEXT(1)                          )
            CALL ERRINT ( '#', GRDEXT(2)                          )
            CALL ERRINT ( '#', GRDEXT(3)                          )
            CALL SIGERR ( 'SPICE(VERTEXNOTINGRID)'                )
            CALL CHKOUT ( 'XDDA'                                  )
            RETURN
            
         END IF

      END DO

C
C     The maximum output voxel array size must be positive.
C     
      IF ( MAXNVX .LT. 1 ) THEN

         CALL CHKIN  ( 'XDDA'                     )
         CALL SETMSG ( 'Maximum voxel list size must be positive '
     .   //            'but was #.'                               )
         CALL ERRINT ( '#', MAXNVX                                )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'                       )
         CALL CHKOUT ( 'XDDA'                                     )
         RETURN

      END IF
      
C
C     Find the largest component of the direction vector.
C
      IAXIS(1) = 1
      MAXCMP   = ABS( RAYDIR(1) )

      DO I = 2, 3

         IF ( ABS( RAYDIR(I) ) .GT. MAXCMP  ) THEN
            IAXIS(1) = I
            MAXCMP   = ABS( RAYDIR(I) )
         END IF

      END DO

C
C     Set the indices of the orthogonal components of the direction
C     vector.  We maintain a right-handed relationship between the axes
C     labeled by IAXIS(1), IAXIS(2), and IAXIS(3):  the third axis is
C     the cross product of the first and second.
C     
      IAXIS(2) = NEXT ( IAXIS(1) )
      IAXIS(3) = NEXT ( IAXIS(2) )
      
C
C     Which voxel contains the vertex? Truncate the vertex
C     coordinates to integers. Add 1 to each coord to compensate
C     for 1 based counting.
C
      DO I = 1, 3

         INTVTX(I) = INT ( VERTEX(IAXIS(I)) )

         ICOORD(I) = 1 + INTVTX(I)

         ICOORD(I) = BRCKTI (  ICOORD(I),  1,  GRDEXT( IAXIS(I) )  )

         VOXLST( IAXIS(I), 1 ) = ICOORD(I)

      END DO 

C
C     Initialize the counter for number of voxels the ray intercepts.
C     The bracketing done above ensures that the coordinates ICOORD of
C     the voxel considered to contain ray's vertex (there is a choice
C     to be made if the vertex lies on a voxel boundary) are within the
C     grid.

      NVX = 1

C
C     Calculate the relative location of vertex within the voxel. The
C     coordinates of a voxel's corners are integer values with each
C     voxel side length 1 (in voxel coords).
C
C     The variable VTXOFF usually has components equal to the
C     fractional parts of the corresponding components of VERTEX(
C     IAXIS(I) ), but the components of VTXOFF may be as large as 1 and
C     are never less than 0.
C
      DO I = 1, 3

         VTXOFF(I) = BRCKTD( VERTEX( IAXIS(I) ) - (ICOORD(I)-1),  
     .                       0.D0, 
     .                       1.D0                               )
      END DO
 
C
C     Compute the lower limit on the magnitudes of RAYDIR( IAXIS(2) )
C     and of RAYDIR( IAXIS(3) ) for which we'll treat those components
C     of the direction vector as non-zero.
C
      LIMIT = ( EPS / GRDEXT( IAXIS(1) ) )  *  ABS( RAYDIR(IAXIS(1)) )
C
C     If the magnitude of RAYDIR( IAXIS(J) ), J = 2 or 3, is below
C     LIMIT, then the ray can pass through the entire grid in the
C     IAXIS(1) direction without its IAXIS(J) component changing by
C     more than EPS. We'll treat this case as though the IAXIS(J)
C     component of the ray were 0.
C
C
C     Determine the error term initial values and increments.
C
C
      AX2ERR = DPMAX()
      AX3ERR = AX2ERR
      S12    = 0.D0
      S13    = 0.D0

C
C     Compute the initial relative distance measurement AX2ERR
C     for the non-primary axis IAXIS(2).
C
      IF (   ABS(  RAYDIR( IAXIS(2) )  ) .GT.  LIMIT  ) THEN
C
C        For any line segment along the ray, S12 is the ratio of the
C        magnitudes of the projections of the segment in the primary
C        and the IAXIS(2) directions.
C      
         S12 = ABS (  RAYDIR( IAXIS(1) ) / RAYDIR( IAXIS(2) )  )
      
         
         IF ( RAYDIR(IAXIS(1)) .GT. 0.D0 ) THEN
C
C           The primary component of the ray's direction is positive.
C           The distance to the next boundary plane in the primary
C           direction is 
C
C              1.D0 - VTXOFF( IAXIS(1) )
C
            IF( RAYDIR(IAXIS(2)) .GT. 0.D0 ) THEN
C
C              The IAXIS(2) component of the ray's direction is
C              positive. The distance to the next boundary plane for
C              the that axis is
C
C                 1.D0 - VTXOFF(2)
C
C              The corresponding change along the primary axis is
C
C                 S12 * ( 1.D0 - VTXOFF(2) )
C
C              The "error" term for IAXIS(2) is this value minus the
C              distance from the vertex to the next boundary in the
C              primary direction.
C
               AX2ERR = S12*(1.D0 - VTXOFF(2)) + VTXOFF(1) - 1.D0

            ELSE
C
C              The IAXIS(2) component of the ray's direction is
C              negative. The distance to the next boundary plane for
C              the that axis is
C
C                 VTXOFF(2)
C
C              The corresponding change along the primary axis is
C
C                 S12 * VTXOFF(2)
C
C              The "error" term for IAXIS(2) is this value minus the
C              distance from the vertex to the next boundary in the
C              primary direction.
C
               AX2ERR = S12*VTXOFF(2) + VTXOFF(1) - 1.D0

            END IF

         ELSE
C
C           The primary component of the ray's direction is negative.
C           The distance to the next boundary plane in the primary
C           direction is 
C
C              VTXOFF( IAXIS(1) )

            IF( RAYDIR(IAXIS(2)) .GT. 0.D0 ) THEN
C
C              The IAXIS(2) component of the ray's direction is
C              positive. The distance to the next boundary plane for
C              the that axis is
C
C                 1.D0 - VTXOFF(2)
C
C              The corresponding change along the primary axis is
C
C                 S12 * ( 1.D0 - VTXOFF(2) )
C
C              The "error" term for IAXIS(2) is this value minus the
C              distance from the vertex to the next boundary in the
C              primary direction.

               AX2ERR = S12*(1.D0 - VTXOFF(2)) - VTXOFF(1)

            ELSE
C
C              The IAXIS(2) component of the ray's direction is
C              negative. The distance to the next boundary plane for
C              the that axis is
C
C                 VTXOFF(2)
C
C              The corresponding change along the primary axis is
C
C                 S12 * VTXOFF(2)
C
C              The "error" term for IAXIS(2) is this value minus the
C              distance from the vertex to the next boundary in the
C              primary direction.

               AX2ERR = S12*VTXOFF(2) - VTXOFF(1)

            END IF

         END IF
      
      END IF
 

C
C     Computations of AX3ERR are analogous to those of AX2ERR.
C     See the comments above.
C
      IF (  ABS(  RAYDIR( IAXIS(3) )  )  .GT.  LIMIT  ) THEN
C
C        For any line segment along the ray, S13 is the ratio of the
C        magnitudes of the projections of the segment in the primary
C        and the IAXIS(3) directions.
C            
         S13 = ABS(  RAYDIR( IAXIS(1) ) / RAYDIR( IAXIS(3) )  )
         
         IF ( RAYDIR( IAXIS(1) ) .GT. 0.D0 ) THEN

            IF( RAYDIR( IAXIS(3) ) .GT. 0.D0 ) THEN
               AX3ERR = S13*(1.D0 - VTXOFF(3) ) + VTXOFF(1) - 1.D0
            ELSE
               AX3ERR = S13*VTXOFF(3) + VTXOFF(1) - 1.D0
            END IF

         ELSE

            IF(  RAYDIR( IAXIS(3) ) .GT. 0.D0  ) THEN
               AX3ERR = S13*(1.D0 - VTXOFF(3)) - VTXOFF(1)
            ELSE
               AX3ERR = S13*VTXOFF(3) - VTXOFF(1)
            END IF

         END IF
      
      END IF

C
C     The "steps" set below are the amounts by which any voxel
C     coordinate changes when the "next" voxel is identified. Only one
C     coordinate changes at a time. The magnitude of each coordinate
C     step is always an integer. The signs of the steps are those of
C     the corresponding components of the ray's direction vector.
C
C     We treat direction components smaller than LIMIT as though
C     they were zero. Note that the IAXIS(1) component of the
C     ray will always have magnitude greater than LIMIT.
C       
      DO I = 1, 3
         
         IF ( RAYDIR( IAXIS(I) ) .GT. LIMIT ) THEN         
C
C           Positive component direction, positive step.
C
            STEP(I) =  1

         ELSE IF ( RAYDIR( IAXIS(I) ) .LT. -LIMIT ) THEN
C
C           Negative component direction, negative step.
C
            STEP(I) = -1
            
         ELSE 
C
C           No component in this direction, no step.
C
            STEP(I) =  0

         END IF
      
      END DO

C
C     Follow the ray until it exits the voxel grid.
C
      DO WHILE ( ZZINGRD( GRDEXT, VOXLST(1,NVX) ) )
 
         IF ( (AX2ERR .LT. 0.D0) .OR. (AX3ERR .LT. 0.D0) ) THEN
C
C           Ray has crossed over into the next voxel in IAXIS(2) or
C           IAXIS(3)
C
            IF ( AX2ERR .LT. AX3ERR ) THEN
C
C              The boundary plane orthogonal to axis IAXIS(2) was hit.
C
               ICOORD(2)  = ICOORD(2) + STEP(2)
               AX2ERR     = AX2ERR    + S12
               NVX        = NVX       + 1
 
            ELSE
C
C              The boundary plane orthogonal to axis IAXIS(3) was hit.
C
               ICOORD(3)  = ICOORD(3) + STEP(3)
               AX3ERR     = AX3ERR    + S13
               NVX        = NVX       + 1
  
            END IF
 
         ELSE
C
C           No change in IAXIS(2) or IAXIS(3), step in IAXIS(1).
C
            ICOORD(1) = ICOORD(1) + STEP(1)
            NVX       = NVX       + 1

            IF ( STEP(2) .NE. 0 ) THEN
               AX2ERR = AX2ERR - 1.D0
            END IF

            IF ( STEP(3) .NE. 0 ) THEN
               AX3ERR = AX3ERR - 1.D0
            END IF
 
         END IF

C
C        Check we have room in VOXLST.
C
         IF( NVX .GT. MAXNVX ) THEN
 
            CALL CHKIN  ( 'XDDA' )
            CALL SETMSG ( 'Index larger than array. '
     .      //            'Index = #1. Array size = #2.' )
            CALL ERRINT ( '#1', NVX                      )
            CALL ERRINT ( '#2', MAXNVX                   )
            CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'         )
 
            CALL CHKOUT ( 'XDDA' )
            RETURN
            
         END IF

C
C        Pack the voxel indices into VOXLST using
C        the values calculated in this loop pass.
C
         DO I = 1, 3
            VOXLST ( IAXIS(I), NVX ) = ICOORD(I)
         END DO

      END DO 
 
C
C     Subtract one off the voxel count since the final voxel 
C     exists outside the grid.
C
      NVX = NVX - 1

      RETURN
      END
