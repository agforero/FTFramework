C$Procedure      PLTNP ( Nearest point on triangular plate )
 
      SUBROUTINE PLTNP ( POINT, V1, V2, V3, PNEAR, DIST )
 
C$ Abstract
C
C     Find the nearest point on a triangular plate to a given point.
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

      DOUBLE PRECISION      POINT ( 3 )
      DOUBLE PRECISION      V1    ( 3 )
      DOUBLE PRECISION      V2    ( 3 )
      DOUBLE PRECISION      V3    ( 3 )
      DOUBLE PRECISION      PNEAR ( 3 )
      DOUBLE PRECISION      DIST
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     POINT      I   A point in 3-dimensional space.
C     V1,
C     V2,
C     V3         I   Vertices of a triangular plate.
C     PNEAR      O   Nearest point on the plate to POINT.
C     DIST       O   Distance between PNEAR and POINT.
C
C$ Detailed_Input
C
C     POINT      is an arbitrary point in 3-dimensional space.
C
C     V1,
C     V2,
C     V3         are 3-vectors constituting the vertices of 
C                a triangular plate.
C
C                The plate is allowed to be degenerate: it may
C                consist of a line segment or of a single point.
C          
C$ Detailed_Output
C
C     PNEAR      is the closest point on the plate to POINT.
C                PNEAR is unique, since the plate is convex.
C
C     DIST       is the distance between POINT and PNEAR.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) The input plate is allowed to be degenerate: it may be
C        a line segment or a single point.
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
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input
C     (if any), the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C
C
C     1) Find the nearest point to the point (2,2,2) on a plate having
C        vertices at the unit basis vectors that lie along the positive
C        X, Y, and Z coordinate axes.
C
C
C        Example code begins here. 
C
C
C              PROGRAM EX1
C              IMPLICIT NONE
C
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         FMT1
C              PARAMETER           ( FMT1 = '(A,3E16.8)' )
C        C
C        C     Local variables
C        C
C              DOUBLE PRECISION      DIST
C              DOUBLE PRECISION      POINT  ( 3 )
C              DOUBLE PRECISION      PNEAR  ( 3 )
C              DOUBLE PRECISION      V1     ( 3 )
C              DOUBLE PRECISION      V2     ( 3 )
C              DOUBLE PRECISION      V3     ( 3 )
C
C        C
C        C     POINT is the input point.
C        C
C              CALL VPACK ( 2.D0, 2.D0, 2.D0, POINT )
C        C
C        C     V1, V2, V3 are the vertices of a plate.
C        C
C              CALL VPACK ( 1.D0, 0.D0, 0.D0, V1 )
C              CALL VPACK ( 0.D0, 1.D0, 0.D0, V2 )
C              CALL VPACK ( 0.D0, 0.D0, 1.D0, V3 )
C        C
C        C     Find the near point on the plate.
C        C
C              CALL PLTNP ( POINT, V1, V2, V3, PNEAR, DIST )
C
C              WRITE (*,*) ' '
C              WRITE (*,FMT1) 'Plate vertex 1 = ', V1
C              WRITE (*,FMT1) 'Plate vertex 3 = ', V2
C              WRITE (*,FMT1) 'Plate vertex 3 = ', V3
C              WRITE (*,FMT1) 'Input point    = ', POINT
C              WRITE (*,*)    ' '
C              WRITE (*,FMT1) 'Near point     = ', PNEAR
C              WRITE (*,FMT1) 'Distance       = ', DIST
C              WRITE (*,*) ' '
C
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran platform,
C     the output was:
C  
C
C      Plate vertex 1 =   0.10000000E+01  0.00000000E+00  0.00000000E+00
C      Plate vertex 3 =   0.00000000E+00  0.10000000E+01  0.00000000E+00
C      Plate vertex 3 =   0.00000000E+00  0.00000000E+00  0.10000000E+01
C      Input point    =   0.20000000E+01  0.20000000E+01  0.20000000E+01
C
C      Near point     =   0.33333333E+00  0.33333333E+00  0.33333333E+00
C      Distance       =   0.28867513E+01
C
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
C-    SPICELIB Version 1.1.2, 01-FEB-2016 (NJB)
C
C        Added code example to header.
C     
C-    DSKLIB Version 1.1.1, 19-MAR-2015 (NJB)
C
C        Fixed spelling error in header.
C
C-    DSKLIB Version 1.1.0, 31-DEC-2014 (NJB)
C
C        Bug fix: vertex indices for outside case, near
C        point on 3rd edge were corrected.
C
C-    DSKLIB Version 1.0.0, 29-SEP-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     nearest point on triangular plate
C
C-&

C
C     SPICELIB functions
C     
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VDOT
      DOUBLE PRECISION      VNORM

      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local variables
C    
      DOUBLE PRECISION      D1
      DOUBLE PRECISION      D2
      DOUBLE PRECISION      D3
      DOUBLE PRECISION      E1     ( 3 )
      DOUBLE PRECISION      E2     ( 3 )
      DOUBLE PRECISION      E3     ( 3 )
      DOUBLE PRECISION      ENORM1 ( 3 )
      DOUBLE PRECISION      ENORM2 ( 3 )
      DOUBLE PRECISION      ENORM3 ( 3 )
      DOUBLE PRECISION      L1
      DOUBLE PRECISION      L2
      DOUBLE PRECISION      L3
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      NP1    ( 3 )
      DOUBLE PRECISION      NP2    ( 3 )
      DOUBLE PRECISION      PDIFF  ( 3 )
      DOUBLE PRECISION      PERP   ( 3 )

      LOGICAL               DEGEN
      LOGICAL               IN1
      LOGICAL               IN2
      LOGICAL               IN3

      IF ( RETURN() ) THEN
         RETURN
      END IF
      
C
C     Use discovery check-in.
C
C
C     Compute the plate's edges.
C
      CALL VSUB ( V2, V1, E1 )
      CALL VSUB ( V3, V2, E2 )
      CALL VSUB ( V1, V3, E3 )

C
C     Compute a normal vector for the plate, if possible.
C     If the plate is degenerate, we'll find out at this point.
C
      CALL VCRSS ( E1, E2, NORMAL )

C
C     Compute the outward normals of the plate's edges in the
C     plate containing the plate.
C
      CALL VCRSS ( E1, NORMAL, ENORM1 )
      CALL VCRSS ( E2, NORMAL, ENORM2 )
      CALL VCRSS ( E3, NORMAL, ENORM3 )

      DEGEN =      VZERO(NORMAL) .OR. VZERO(ENORM1)
     .        .OR. VZERO(ENORM2) .OR. VZERO(ENORM3)

      IF ( DEGEN ) THEN
C
C        The "plate" is a line segment or point. Determine
C        which case we have.
C
         L1 = VNORM(E1)
         L2 = VNORM(E2)
         L3 = VNORM(E3)

         IF ( ( L1 .EQ. 0.D0 ) .AND. ( L2 .EQ. 0.D0 ) ) THEN
C
C           Up to round-off error, the vertices coincide.
C           The vertex V1 for practical purposes is the plate.
C
            CALL VEQU ( V1, PNEAR )

            DIST = VDIST( PNEAR, POINT )

         ELSE
C
C           The plate is a line segment having positive length.
C           One of the edges will coincide with the segment.
C           Determine which vertices are the endpoints.
C
            IF ( L1 .GT. MAX(L2, L3) ) THEN
C
C              The segment is bounded by V1 and V2.
C
               CALL NPSGPT( V1, V2, POINT, PNEAR, DIST )

            ELSE IF ( L2 .GT. MAX(L3, L1) ) THEN
C
C              The segment is bounded by V2 and V3.
C
               CALL NPSGPT( V2, V3, POINT, PNEAR, DIST )

            ELSE
C
C              The segment is bounded by V3 and V1.
C
               CALL NPSGPT( V3, V1, POINT, PNEAR, DIST )

            END IF

         END IF
C
C        The outputs are set for the degenerate cases.
C
         RETURN

      END IF

C
C     We have a non-degenerate plate. NORMAL has unit length.
C
C     We'll treat V1 as an origin in the plane containing
C     the plate. Find the offset of the POINT from V1, and
C     find the component of this offset orthogonal to NORMAL.
C
      CALL VSUB  ( POINT, V1,     PDIFF )
      CALL VPERP ( PDIFF, NORMAL, PERP  )

C
C     Determine whether V1+PERP is inside the plate.
C
C     Note that the "line constants" for edges 1 and 3
C     are zero, since these edges contain V1. The line 
C     constant for edge 2 is that of the offset of V2
C     from V1; this offset is edge 1.
C
      IN1 = VDOT( PERP, ENORM1 ) .LE. 0.D0 
      IN2 = VDOT( PERP, ENORM2 ) .LE. VDOT(E1,ENORM2)
      IN3 = VDOT( PERP, ENORM3 ) .LE. 0.D0 
 
      
      IF ( IN1 .AND. IN2 .AND. IN3 ) THEN 
C
C        V1+PERP is inside the plate. It is the closest
C        point on the plate to POINT.
C
         CALL VADD ( V1, PERP, PNEAR )
C
C        We have the near point; set the distance.
C        
         DIST = VDIST( PNEAR, POINT )

      ELSE
C
C        PERP is outside the plate. The nearest point
C        on the plate to POINT is on one of the edges.
C
C        We'll use the "in" flags to reduce the number
C        of point-edge distance computations.
C
         IF ( ( .NOT. IN1 ) .AND. ( IN2 .AND. IN3 ) ) THEN
C
C           The solution must be on the first edge.
C
            CALL NPSGPT ( V1, V2, POINT, PNEAR, DIST )

         ELSE IF ( ( .NOT. IN2 ) .AND. ( IN3 .AND. IN1 ) ) THEN

C
C           The solution must be on the second edge.
C
            CALL NPSGPT ( V2, V3, POINT, PNEAR, DIST )

         ELSE IF ( ( .NOT. IN3 ) .AND. ( IN1 .AND. IN2 ) ) THEN
C
C           The solution must be on the third edge.
C
            CALL NPSGPT ( V3, V1, POINT, PNEAR, DIST )
       
         ELSE 
C
C           Compute solutions on all three edges and pick 
C           the best one.
C
            CALL NPSGPT ( V1, V2, POINT, NP1,   D1 )
            CALL NPSGPT ( V2, V3, POINT, NP2,   D2 )
            CALL NPSGPT ( V3, V1, POINT, PNEAR, D3 )

            IF ( D1 .LE. MIN(D2, D3) ) THEN

               CALL VEQU ( NP1, PNEAR )
               DIST = D1

            ELSE IF ( D2 .LE. MIN(D3, D1) ) THEN

               CALL VEQU ( NP2, PNEAR )
               DIST = D2

            ELSE
C
C              PNEAR is already set.
C
               DIST = D3

            END IF

         END IF

      END IF

      RETURN
      END
