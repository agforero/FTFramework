C$Procedure ZZRAYBOX ( Ray-box intercept )
 
      SUBROUTINE ZZRAYBOX ( VERTEX, RAYDIR, BOXORI, EXTENT, XPT, FOUND )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Find the surface intercept of a ray on a specified box in
C     three-dimensional space.
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
C     MATH 
C     GEOMETRY
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      BOXORI ( 3 )
      DOUBLE PRECISION      EXTENT ( 3 )
      DOUBLE PRECISION      XPT    ( 3 )
      LOGICAL               FOUND
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VERTEX     I   Vertex of ray.
C     RAYDIR     I   Direction vector of ray.
C     BOXORI     I   Box origin: corner having minimum coordinates.
C     EXTENT     I   Box extents in X, Y, Z directions.
C     XPT        O   Ray-box intercept.
C     FOUND      O   Flag indicating whether intercept was found.
C
C$ Detailed_Input
C
C     VERTEX     is the vertex of a ray in three dimensional space.
C
C     RAYDIR     is the direction vector of the input ray.
C
C     BOXORI     is a vector representing the "origin" of a box in
C                three-dimensional space. BOXORI is the corner of
C                the box at which each coordinate attains its 
C                minimum value.
C
C     EXTENT     is an array containing the box's edge lengths in
C                the X, Y, and Z directions. All lengths must be
C                positive.
C
C$ Detailed_Output
C
C     XPT        is the point of intersection closest to VERTEX
C                of the ray and the box. If VERTEX lies within the
C                box, XPT is set equal to VERTEX.
C                
C                XPT is undefined if FOUND is .FALSE.
C                
C
C     FOUND      is set to .TRUE. if and only if the ray intersects
C                the box.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the input ray is the zero vector, the error
C         SPICE(ZEROVECTOR) is signaled.
C
C     2)  If any element of EXTENT is non-positive, the error
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
C     See use in DSKX02.
C
C$ Restrictions
C
C     This is a private routine; its interface or functionality
C     may be changed without notice. SPICE user applications 
C     should not call this routine directly.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     J.A. Bytof      (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 06-JUN-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compute ray-box intercept
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      VNORM

      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local parameters
C
      INTEGER               LEFT
      PARAMETER           ( LEFT   = 1 )

      INTEGER               MIDDLE
      PARAMETER           ( MIDDLE = 2 )

      INTEGER               RIGHT
      PARAMETER           ( RIGHT  = 3 )

C
C     Local variables
C 
      DOUBLE PRECISION      CENTER ( 3 )
      DOUBLE PRECISION      LIMIT
      DOUBLE PRECISION      MAXT
      DOUBLE PRECISION      NEAR   ( 3 )
      DOUBLE PRECISION      OFFSET ( 3 )
      DOUBLE PRECISION      PLNDST ( 3 )
      DOUBLE PRECISION      R
      DOUBLE PRECISION      SPHXPT ( 3 )
      DOUBLE PRECISION      SPHVTX ( 3 )
      DOUBLE PRECISION      T      ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      VTEMP  ( 3 )

      INTEGER               SECTOR ( 3 )
      INTEGER               I
      INTEGER               MAXIDX

      LOGICAL               SPHFND
      

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     No intercept has been found yet.
C
      FOUND = .FALSE.

C
C     Check for a zero ray direction vector.
C     
      IF ( VZERO(RAYDIR) ) THEN

         CALL CHKIN  ( 'ZZRAYBOX'                                      )
         CALL SETMSG ( 'Input ray direction was the zero vector; this '
     .   //            'vector must be non-zero.'                      )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                             )
         CALL CHKOUT ( 'ZZRAYBOX'                                      )
         RETURN

      END IF

      CALL VHAT ( RAYDIR, UDIR )

C
C     Check the box extents.
C
      IF (  MIN( EXTENT(1), EXTENT(2), EXTENT(3) ) .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZRAYBOX'                                    )
         CALL SETMSG ( 'All box extents should be strictly positive '
     .   //            'but the extents were #, #, #.'               )
         CALL ERRDP  ( '#',  EXTENT(1)                               )
         CALL ERRDP  ( '#',  EXTENT(2)                               )
         CALL ERRDP  ( '#',  EXTENT(3)                               )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                      )
         CALL CHKOUT ( 'ZZRAYBOX'                                    )
         RETURN

      END IF

C
C     Compute the coordinates of the center of the box, and compute
C     the offset of the ray's vertex from the center.
C
      DO I = 1, 3
         CENTER(I) = BOXORI(I) + EXTENT(I)/2
      END DO

      CALL VSUB  ( VERTEX, CENTER, OFFSET )
     
C
C     If the ray's vertex is inside the box, we consider the
C     vertex to be the intercept.
C
      IF (       ( ABS( OFFSET(1) ) .LE. EXTENT(1)/2 )
     .     .AND. ( ABS( OFFSET(2) ) .LE. EXTENT(2)/2 )
     .     .AND. ( ABS( OFFSET(3) ) .LE. EXTENT(3)/2 )  ) THEN
         
         CALL VEQU ( VERTEX, XPT )
         FOUND = .TRUE.
         
         RETURN

      END IF
C
C     Compute the intercept of the ray on the surface of a bounding
C     sphere that contains the box. Let R be the radius of this sphere.
C
      R = 0.5D0 * (1.D0 + 1.D-3) * VNORM( EXTENT ) 

      IF ( VNORM(OFFSET) .LT. R ) THEN
C
C        The vertex is already inside the bounding sphere.
C
         CALL VEQU ( OFFSET, SPHXPT )

      ELSE

         CALL SURFPT ( OFFSET, UDIR, R, R, R, SPHXPT, SPHFND )

         IF ( .NOT. SPHFND ) THEN
C
C           The ray misses the bounding sphere.
C        
            RETURN

         END IF

      END IF

C
C     Shift the sphere intercept so as to be relative to the 
C     box's origin. From this point on, we'll treat BOXORI 
C     as the origin of the reference frame.
C
      DO I = 1, 3
         SPHVTX(I) = SPHXPT(I) + CENTER(I) - BOXORI(I)
      END DO

C
C     Classify the position of the vertex relative to the planes
C     bounding the box: each coordinate will be classified as
C     "left," "middle," or "right" depending on whether it is
C     less than the lower bound for that coordinate, between
C     the bounds, or greater than the upper bound.
C
      DO I = 1, 3

         IF ( SPHVTX(I) .LT. 0.D0 ) THEN

            SECTOR(I) = LEFT
            NEAR(I)   = 0.D0

         ELSE IF ( SPHVTX(I) .GT. EXTENT(I) ) THEN
            
            SECTOR(I) = RIGHT
            NEAR(I)   = EXTENT(I)

         ELSE
            SECTOR(I) = MIDDLE
            NEAR(I)   = 0.D0
         END IF

      END DO

C
C     At this point, SPHVTX is a point on the ray that is outside,
C     but close to, the box. SPHVTX is an offset from BOXORI we'll
C     need to add BOXORI to it to obtain the corresponding point in
C     the input reference frame.
C
C     We'll use SPHVTX as the new ray vertex.
C
C     Find the distances of the vertex's components from the nearest
C     bounding planes of the box; find the corresponding distances
C     along the ray that would be traveled in order to move each
C     component from the vertex to the nearest bounding plane. Call the
C     latter distance for the Ith coordinate T(I). We're only
C     interested in the vertex components that are "outside" the
C     bounding planes. Mark the values of T(I) for components in the
C     "middle" using the value -1.
C
C     Find the index of the maximum T value while we're at it. If
C     there's an intercept, it occurs at the point on the ray
C     corresponding to the maximum value of T.
C
      MAXIDX =  1
      MAXT   = -1.D0


      DO I = 1, 3

         T(I) = -1.D0

         IF ( SECTOR(I) .NE. MIDDLE ) THEN

            PLNDST(I) = NEAR(I) - SPHVTX(I)

C
C           Prepare for a "safe" division.
C
            LIMIT = 2 * R * ABS( UDIR(I) )            

            IF ( ABS( PLNDST(I) ) .GT. LIMIT ) THEN
C
C              The ray can't get to the nearest bounding plane
C              before exiting the bounding sphere. No intersection
C              is possible.
C              
               RETURN

            END IF
C
C           The magnitude of the following quotient is bounded by 2R.
C
            T(I) = PLNDST(I) / UDIR(I)

            IF ( T(I) .LT. 0.D0 ) THEN
C
C              This component of the ray is going in the wrong 
C              direction. No intersection is possible.
C
               RETURN

            END IF

            IF ( T(I) .GT. MAXT ) THEN

               MAXIDX = I
               MAXT   = T(I)

            END IF
               
         END IF

      END DO

C
C     We should have a positive value of T for at least one
C     coordinate. However, if we don't, there's no intersection.
C
      IF ( MAXT .LT. 0.D0 ) THEN
         RETURN
      END IF

C
C     Compute the candidate intercept. Note that we're now working
C     in a frame centered at the box origin.
C
      CALL VLCOM ( 1.D0, SPHVTX, MAXT, UDIR, XPT )

C
C     Decide whether XPT is actually on the surface of the box.
C     Sharpen XPT as part of the process.
C
      DO I = 1, 3

         IF ( I .EQ. MAXIDX ) THEN
C
C           XPT is supposed to lie exactly on the bounding plane
C           orthogonal to the Ith axis and nearest to SPHVTX.
C         
            XPT(I) = NEAR(I)

         ELSE

            IF ( SECTOR(I) .EQ. MIDDLE ) THEN
C
C              The Ith component of the vertex is between the
C              bounding planes for the Ith coordinate. If the
C              Ith component of XPT is outside these bounds,
C              the ray misses the box.
C
               IF (      ( XPT(I) .LT. 0.D0      )
     .              .OR. ( XPT(I) .GT. EXTENT(I) ) ) THEN

                  RETURN

               END IF

            ELSE
C
C              The Ith component of the vertex SPHVTX is outside of the
C              bounding planes for the Ith coordinate. Since T(MAXIDX)
C              is greater than or equal to T(I), XPT(I) should be on or
C              past the bounding plane closest to SPHVTX(I). Sharpen
C              XPT(I) if necessary. If XPT(I) is beyond the bounding
C              plane farthest from SPHVTX(I), no intersection can
C              exist.
C 
               IF ( SECTOR(I) .EQ. LEFT ) THEN
C
C                 Sharpen the Ith component of XPT.
C
                  XPT(I) = MAX( XPT(I), 0.D0 )

                  IF ( XPT(I) .GT. EXTENT(I) ) THEN
C
C                    The ray hits the MAXIDX face too far away from
C                    SPHVTX(I). There's no intersection with the box.
C
                     RETURN

                  END IF


               ELSE
C
C                 SECTOR(I) .EQ. RIGHT
C
C                 Sharpen the Ith component of XPT.
C
                  XPT(I) = MIN( XPT(I), EXTENT(I) )

                  IF ( XPT(I) .LT. 0.D0 ) THEN
C
C                    The ray hits the MAXIDX face too far away from
C                    SPHVTX(I). There's no intersection with the box.
C
                     RETURN

                  END IF

               END IF
C
C              End of block in which the Ith component of XPT is
C              either sharpened or found to be off the surface of the
C              box. This block deals with the components other than
C              MAXIDX.
C
            END IF
C
C           End of block in which the Ith component of XPT is either
C           sharpened or found to be off the surface of the box. This
C           block deals with all components.
C
         END IF

      END DO 
C
C     End of loop in which XPT is either sharpened or found to be off
C     the surface of the box. Getting here means XPT is valid.
C
C     Shift XPT to the input reference frame.
C     
      CALL VADD ( XPT,   BOXORI, VTEMP )
      CALL VEQU ( VTEMP,         XPT   )

      FOUND = .TRUE.

      RETURN
      END
