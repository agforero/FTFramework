C$Procedure ZZRECBOX (Bounding box for rectangular volume element)
 
      SUBROUTINE ZZRECBOX ( BOUNDS, CENTER, LX, LY, LZ, RADIUS )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Compute a bounding box center and radius for a rectangular volume
C     element.
C
C     The outputs are the box's center, dimensions in the X, Y, and Z
C     directions, and the box's radius. The box itself coincides with
C     the input box.
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
      DOUBLE PRECISION      LX
      DOUBLE PRECISION      LY
      DOUBLE PRECISION      LZ
      DOUBLE PRECISION      RADIUS

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BOUNDS     I   Array of coordinate bounds for segment.
C     CENTER     O   Center of bounding box.
C     LX         O   Extent of box in the X direction.
C     LY         O   Extent of box in the Y direction.
C     LZ         O   Extent of box in the Z direction.
C     RADIUS     O   Radius of box.
C
C$ Detailed_Input
C
C     BOUNDS     is a 2 x 3 double precision array containing bounds
C                on the rectangular coordinates of the spatial region
C                (volume element) covered by a DSK segment. The
C                contents of BOUNDS are:
C
C                   BOUNDS(*,1)              X bounds
C                   BOUNDS(*,2)              Y bounds
C                   BOUNDS(*,3)              Z bounds
C
C                Elements (1,*) are lower bounds; elements (2,*) are
C                upper bounds.  
C
C$ Detailed_Output
C
C     CENTER     is a double precision 3-vector representing the center
C                of the volume specified by BOUNDS.
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
C     None.
C
C$ Exceptions
C
C     1) If the minimum value of of any coordinate exceeds the maximum,
C        the error SPICE(BOUNDSOUTOFORDER) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports the DSK ray-surface intercept segment
C     selection algorithms.
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
C-    SPICELIB Version 1.0.0, 04-JUN-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     bounding box radius for rectangular volume element
C
C-&

C
C     SPICELIB functions
C     
      DOUBLE PRECISION      VNORM

      LOGICAL               RETURN

C
C     Local parameters
C
C
C     Element boundary indices:
C     
      INTEGER               LOWER
      PARAMETER           ( LOWER  = 1 )

      INTEGER               UPPER
      PARAMETER           ( UPPER  = 2 )
      
C
C     Local variables
C
      DOUBLE PRECISION      DIAG   ( 3 )
      DOUBLE PRECISION      L      ( 3 )
      DOUBLE PRECISION      MAXCOR ( 3 )
      DOUBLE PRECISION      MINCOR ( 3 )

      INTEGER               I

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
      DO I = 1, 3
         
         MINCOR(I) = BOUNDS(LOWER,I)
         MAXCOR(I) = BOUNDS(UPPER,I)
         L(I)      = MAXCOR(I) - MINCOR(I)

         IF ( L(I) .LE. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZRECBOX'                               )
            CALL SETMSG ( 'Coordinate # bounds were #:#; bounds '
     .      //            'must be strictly increasing.'           )
            CALL ERRINT ( '#',  I                                  )
            CALL ERRDP  ( '#',  MINCOR(I)                          )
            CALL ERRDP  ( '#',  MAXCOR(I)                          )
            CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)'                )
            CALL CHKOUT ( 'ZZRECBOX'                               )
            RETURN

         END IF

      END DO
 
C
C     Set output box dimensions.
C
      LX = L(1)
      LY = L(2)
      LZ = L(3)

C
C     Compute the coordinates of the center of the box.
C
      DO I = 1, 3
         CENTER(I) = MINCOR(I) + ( L(I)/2 )        
      END DO
      
C
C     The radius is the distance from the center of the box 
C     to any corner.
C
      CALL VPACK ( LX/2, LY/2, LZ/2, DIAG )
      
      RADIUS = VNORM( DIAG )

      RETURN
      END
