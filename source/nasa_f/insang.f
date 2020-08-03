C$Procedure INSANG ( Inside Tetrahedral Angle )
 
      SUBROUTINE INSANG ( V, E1, E2, E3, FOUND, SCALE )
 
C$ Abstract
C
C     Determine if a given vector lies inside the solid tetrahedral
C     angle determined by 3 vectors. If it does, return the
C     point where the scale factor such that SCALE*V lies in the
C     plane spanned by E1, E2, and E3.
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
C      None.
C
C$ Keywords
C
C      VECTOR
C
C$ Declarations
 
      IMPLICIT NONE
      DOUBLE PRECISION      V     ( 3 )
      DOUBLE PRECISION      E1    ( 3 )
      DOUBLE PRECISION      E2    ( 3 )
      DOUBLE PRECISION      E3    ( 3 )
      LOGICAL               FOUND
      DOUBLE PRECISION      SCALE
 
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     V          I   Vector to test for "betweenness"
C     E1         I   First edge of the tetrahedral angle
C     E2         I   Second edge of the tetrahedral angle
C     E3         I   Third edge of the tetrahedral angle
C     FOUND      O   Indicates whether V lies in the solid angle
C     SCALE      O   Scale times V is in the triangle E1,E2,E3
C
C$ Detailed_Input
C
C     V          is a 3-vector.  This is the vector to test to see
C                if it lies between the 3 vectors E1, E2 and E3
C
C     E1         are the three edges of a solid tetrahedral angle.
C     E2         (See particulars for a discussion of the solid
C     E3         angle).
C
C$ Detailed_Output
C
C     FOUND      indicates that V lies inside the solid tetrahedral
C                angle determined by E1, E2 and E3.
C
C
C     SCALE      if V lies inside the solid tetrahedral angle given
C                by E1, E2 and E3, SCALE*V is the point is the positive
C                scalar multiple of V that pierces the triangle
C                determined by the points E1, E2, E3.
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     Error free.
C
C     1) If E1, E2 and E3 are not linearly independent, the routine
C        returns FALSE. SCALE will be set to 0.
C
C     2) If V is the zero vector, the routine returns FALSE.
C        SCALE will be set to 0.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Given 3 linearly independent vectors E1, E2, and E3 the
C     set of vectors a*E1 + b*E2 + c*E3  where a, b, and c
C     are non-negative form a region of space that is a tetrahedral
C     solid angle.  If you cut this solid angle with a plane
C     that intersects all three rays from the origin determined
C     by E1, E2 and E3 you will get a tetrahedron (a 4-sided
C     solid with each face a triangle).
C
C     This routine determines whether the ray associated with
C     a vector V lies inside the tetrahedral angle E1,E2,E3.
C     Moreover, if V lies inside this angle, this routine returns
C     the scale factor SCALE such that the point SCALE*V
C     lies in the plane containing the points E1, E2 and E3.
C     This is necessarily a point in the triangle determined by
C     E1, E2 and E3.
C
C$ Examples
C
C     Suppose you have a triangle in space specified by three
C     vertices P1, P2 and P3 and that an observer at location
C     OBS is looking along the ray emanating from OBS with
C     direction V.  Does this ray intersect the triangle
C     P1, P2, P3?  Using this routine, you can answer this
C     question and give the point of intersection if there is
C     one.  Here's how.
C
C     First construct the vectors from OBS to the corners of
C     the triangle.
C
C     CALL VSUB ( P1, OBS, E1 )
C     CALL VSUB ( P2, OBS, E2 )
C     CALL VSUB ( P3, OBS, E3 )
C
C     Now see if V lies between the vectors E1, E2, E3 and return
C     the intersection point if it does.
C
C     CALL INSANG ( V, E1, E2, E3, FOUND, SCALE )
C
C     If there was an intersection, add SCALE*V to OBS to get the
C     point of intersection.  Otherwise say there was no intersection.
C
C     IF ( FOUND ) THEN
C
C        CALL VLCOM ( 1.0D0, OBS, SCALE, V, POINT )
C
C        WRITE (*,*) 'The ray intersects the triangle at:
C        WRITE (*,*) POINT(1)
C        WRITE (*,*) POINT(2)
C        WRITE (*,*) POINT(3)
C
C     ELSE
C
C        WRITE (*,*) 'There is no intersection.'
C
C     END IF
C
C$ Restrictions
C
C     This routine can suffer from extreme loss of precision if
C     the vectors E1, E2, E3 are too long compared to the lengths
C     of the line segments formed by their pairwise differences.
C     
C     The user of this routine must ensure that the inputs are 
C     suitable.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     W.L. Taber      (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.2, 02-FEB-2016 (NJB)
C
C        Fixed comment typos. Updated Restrictions.
C
C-    SPICELIB Version 1.0.1, 08-OCT-2009 (NJB)
C
C        Updated header.
C
C-    SPICELIB Version 1.0.0, 09-JUN-1996 (WLT)
C
C
C-&
 
C$ Index_Entries
C
C     Determine if a vector lies in a solid tetrahedral angle
C
C-&
C
C     SPICELIB Functions
C
      DOUBLE PRECISION      VDOT
 
 
C
C     Local Variables
C
      DOUBLE PRECISION      DENOM
      DOUBLE PRECISION      EN
      DOUBLE PRECISION      NORM12 ( 3 )
      DOUBLE PRECISION      NORM23 ( 3 )
      DOUBLE PRECISION      NORM31 ( 3 )
      DOUBLE PRECISION      VN12
      DOUBLE PRECISION      VN23
      DOUBLE PRECISION      VN31
C
C     Our initial value for SCALE is zero.  When we have better
C     information, we'll change this.
C
      SCALE = 0.0D0
 
C
C     First we construct a normal to the plane spanned by E1 and E2
C     and make sure that we don't get a zero vector.  If we
C     get the zero vector, E1 and E2 are linearly dependent so we
C     set the value of FOUND to FALSE and return.
C
      CALL VCRSS ( E1, E2, NORM12 )
 
C
C     First make sure V and E3 are in the same half space
C     bounded by E1 and E2.  If they are not, we can return.
C
      VN12 = VDOT( V,  NORM12 )
      EN   = VDOT( E3, NORM12 )
C
C     Determine whether NORML and E3 are perpendicular.  If they
C     are perpendicular, E3 is a linear combination of E1 and E2.
C     In this case set FOUND to FALSE and return.
C
      IF ( EN .EQ. 0.0D0 ) THEN
         FOUND = .FALSE.
         RETURN
      END IF
C
C     Now check to see if V and E3 are in the same half space.  If
C     not, we can stop and return the value FALSE.
C
      IF (       EN   .GT. 0.0D0
     .     .AND. VN12 .LT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      ELSE IF (       EN   .LT. 0.0D0
     .          .AND. VN12 .GT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      END IF
 
C
C     Now check that V and E1 are on the same side of the plane
C     spanned by E2 and E3.  Note we don't have to compute EN
C     again <( E2 x E3 ), E1 >  because of the vector identity
C
C       < (E1 x E2), E3 > =  < (E2 x E3), E1 > = < (E3 x E1), E2 >
C
      CALL VCRSS ( E2, E3, NORM23 )
      VN23 = VDOT  ( V,      NORM23 )
 
C
C     The following tests are the same as in the previous case.
C
      IF (       EN   .GT. 0.0D0
     .     .AND. VN23 .LT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      ELSE IF (       EN   .LT. 0.0D0
     .          .AND. VN23 .GT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      END IF
 
C
C     Finally check to see if V and E2 are in the same half space
C     bounded by E3 and E2
C
      CALL VCRSS ( E3, E1, NORM31 )
      VN31 = VDOT  ( V,      NORM31 )
 
      IF (       EN   .GT. 0.0D0
     .     .AND. VN31 .LT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      ELSE IF (       EN   .LT. 0.0D0
     .          .AND. VN31 .GT. 0.0D0 ) THEN
 
         FOUND = .FALSE.
         RETURN
 
      END IF
 
C
C     If you get this far, we know that V is lies in the intersection
C     of the half spaces determined by the various combinations of
C     E1, E2 and E3.
C
      FOUND = .TRUE.
 
C
C     Now find the intersection. First get a normal to the triangle.
C     One way to get the normal is to find the vector cross
C     product
C
C       NORML = ( E2 - E1 ) x ( E3 - E1 )
C
C     However, this can be rewritten as:
C
C        NORML = E2 x E3 - E1 x E3 - E2 x E1 + E1 x E1
C
C              = E2 x E3 + E3 x E1 + E1 x E2
C
C     But we already have the three components E2 x E3, ... etc.
C     in the vectors NORM12, NORM23, NORM31
C
C     Now we need to find the scalar multiple t*V such that
C
C        < tV - E1, NORML > = 0
C
C     But this can be rewritten as:
C
C        t < V, NORML > = < E1, NORML >
C
C     Solving for t yields
C
C      t = < E1, NORML > / < V, NORML >
C
C        = < E1, E1xE2 + E2xE3 + E3xE1 > / <  V, E1xE2 + E2xE3 + E3xE1 >
C
C        = ( 0 + <E1, E2xE3> + 0 ) / (<V,E1xE2> + <V,E2xE3> + <V,E3xE1>)
C
C        =  EN / ( VN12 + VN23 + VN31 )
C
      DENOM = VN12 + VN23 + VN31
 
      IF ( DENOM .EQ. 0.0D0 ) THEN
         FOUND = .FALSE.
      ELSE
         FOUND = .TRUE.
         SCALE =  EN/DENOM
      END IF
 
      RETURN
      END
