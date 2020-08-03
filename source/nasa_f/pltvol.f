C$Procedure PLTVOL ( Compute volume of plate model )
 
      DOUBLE PRECISION FUNCTION PLTVOL ( NV, VRTCES, NP, PLATES )

C$ Abstract
C
C     Compute the volume of a three-dimensional region bounded by a
C     collection of triangular plates.
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
C     DSK
C     GEOMETRY
C     MATH
C     TOPOGRAPHY
C
C$ Declarations

      IMPLICIT NONE

      INTEGER               NV 
      DOUBLE PRECISION      VRTCES ( 3, NV )
      INTEGER               NP
      INTEGER               PLATES ( 3, NP )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     NV         I   Number of vertices.
C     VRTCES     I   Array of vertices.
C     NP         I   Number of triangular plates.
C     PLATES     I   Array of plates.
C
C     The function returns the volume of the spatial region bounded
C     by the plates.
C
C$ Detailed_Input
C
C     NV             is the number of vertices comprising the plate
C                    model.
C
C     VRTCES         is an array containing the plate model's vertices.
C                    Elements
C
C                       VRTCES( 1, I )
C                       VRTCES( 2, I )
C                       VRTCES( 3, I )
C
C                    are, respectively, the X, Y, and Z components of
C                    the Ith vertex.
C
C                    This routine doesn't associate units with the
C                    vertices.
C
C
C     NP             is the number of triangular plates comprising the
C                    plate model.
C
C     PLATES         is an array containing 3-tuples of integers
C                    representing the model's plates. The elements of
C                    PLATES are vertex indices. The vertex indices are
C                    1-based: vertices have indices ranging from 1 to
C                    NV. The elements
C
C                       PLATES( 1, I )
C                       PLATES( 2, I )
C                       PLATES( 3, I )
C
C                    are, respectively, the indices of the vertices
C                    comprising the Ith plate.
C
C                    Note that the order of the vertices of a plate is
C                    significant: the vertices must be ordered in the
C                    positive (counterclockwise) sense with respect to
C                    the outward normal direction associated with the
C                    plate. In other words, if V1, V2, V3 are the
C                    vertices of a plate, then
C
C                       ( V2 - V1 )  x  ( V3 - V2 )
C
C                    points in the outward normal direction. Here
C                    'x' denotes the vector cross product operator.
C
C                    
C$ Detailed_Output
C
C     The function returns the volume of the spatial region bounded
C     by the plates.
C
C     If the components of the vertex array have length unit L, then the
C     output volume has units
C
C         3
C        L
C                              
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) The input plate model must define a spatial region with
C        a boundary. This routine does not check the inputs to
C        verify this condition. See the Restrictions section below.
C
C     2) If the number of vertices is less than 4, the error
C        SPICE(TOOFEWVERTICES) is signaled.
C
C     3) If the number of plates is less than 4, the error
C        SPICE(TOOFEWPLATES) is signaled.
C
C     4) If any plate contains a vertex index outside of the range
C
C           [1, NV]
C
C        the error SPICE(INDEXOUTOFRANGE) will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine computes the volume of a spatial region bounded by
C     a set of triangular plates. If the plate set does not actually
C     form the boundary of a spatial region, the result of this routine
C     is invalid.
C
C     Examples:
C
C        Valid inputs
C        ------------
C        Tetrahedron
C        Box
C        Tiled ellipsoid
C        Two disjoint boxes
C
C        Invalid inputs
C        --------------
C        Single plate
C        Tiled ellipsoid with one plate removed
C        Two boxes with intersection having positive volume
C
C
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input
C     (if any), the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C
C
C     1) Compute the volume of the pyramid defined by the four
C        triangular plates whose vertices are the 3-element
C        subsets of the set of vectors
C
C           ( 0, 0, 0 )
C           ( 1, 0, 0 )
C           ( 0, 1, 0 )
C           ( 0, 0, 1 )
C
C
C        Example code begins here. 
C
C           PROGRAM EX1
C           IMPLICIT NONE
C     C
C     C     Compute the volume of a plate model representing the pyramid
C     C     with one vertex at the origin and the other vertices
C     C     coinciding with the standard basis vectors.
C     C
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      PLTVOL
C     C
C     C     Local parameters
C     C
C           INTEGER               NVERT
C           PARAMETER           ( NVERT  = 4 )
C
C           INTEGER               NPLATE
C           PARAMETER           ( NPLATE = 4 )
C     C
C     C     Local variables
C     C
C           DOUBLE PRECISION      VRTCES ( 3, NVERT  )
C           DOUBLE PRECISION      VOL
C
C           INTEGER               PLATES ( 3, NPLATE )
C     C
C     C     Initial values
C     C
C     C     The plates defined below lie in the following planes,
C     C     respectively:
C     C
C     C        Plate 1:    { P :  < P, (-1,  0,  0) > = 0 }
C     C        Plate 2:    { P :  < P, ( 0, -1,  0) > = 0 }
C     C        Plate 3:    { P :  < P, ( 0,  0, -1) > = 0 }
C     C        Plate 4:    { P :  < P, ( 1,  1,  1) > = 1 }
C     C
C           DATA                  PLATES /    1,     4,     3,
C          .                                  1,     2,     4,
C          .                                  1,     3,     2,
C          .                                  2,     3,     4 /
C
C           DATA                  VRTCES / 0.D0,  0.D0,  0.D0,
C          .                               1.D0,  0.D0,  0.D0,
C          .                               0.D0,  1.D0,  0.D0,
C          .                               0.D0,  0.D0,  1.D0 /
C
C
C           VOL = PLTVOL ( NVERT, VRTCES, NPLATE, PLATES )
C
C           WRITE (*,*) 'Expected volume =    1/6'
C           WRITE (*,*) 'Computed volume = ', VOL
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran platform,
C     the output was:
C
C        Expected volume =    1/6
C        Computed volume =   0.16666666666666666
C
C
C$ Restrictions
C
C     1) The plate collection must describe a surface and enclose a 
C        volume such that the divergence theorem (see [1]) is
C        applicable.
C
C$ Literature_References
C
C     [1]  Calculus, Vol. II. Tom Apostol. John Wiley & Sons, 1969.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 24-OCT-2016 (NJB)
C
C        Based on original 11-FEB-2011
C-&
 
C$ Index_Entries
C
C     compute plate model volume
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      DET

      LOGICAL               RETURN


C
C     Local variables
C     
      DOUBLE PRECISION      M ( 3, 3 )

      INTEGER               I
      INTEGER               J


C
C     The function must have an initial value.
C
      PLTVOL  =  0.D0

C
C     This routine uses discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Check the vertex and plate counts.
C
      IF ( NV .LT. 4 ) THEN

         CALL CHKIN  ( 'PLTVOL'                                      )
         CALL SETMSG ( 'At least 4 vertices are needed, but NV = #.' )
         CALL ERRINT ( '#',  NV                                      )
         CALL SIGERR ( 'SPICE(TOOFEWVERTICES)'                       )
         CALL CHKOUT ( 'PLTVOL'                                      )
         RETURN

      END IF

      IF ( NP .LT. 4 ) THEN

         CALL CHKIN  ( 'PLTVOL'                                    )
         CALL SETMSG ( 'At least 4 plates are needed, but NP = #.' )
         CALL ERRINT ( '#',  NP                                    )
         CALL SIGERR ( 'SPICE(TOOFEWPLATES)'                       )
         CALL CHKOUT ( 'PLTVOL'                                    )
         RETURN

      END IF

C
C     Make sure the vertex indices are in the range [1, NV].
C
      DO I = 1, NP

         DO J = 1, 3

            IF (      ( PLATES(J,I) .LT. 1  ) 
     .           .OR. ( PLATES(J,I) .GT. NV )  ) THEN

               CALL CHKIN  ( 'PLTVOL'                                  )
               CALL SETMSG ( 'Vertex indices must be in the range '
     .         //            '[1, NV] for all SPICE language versions. '
     .         //            'The input value of NV was #. Vertex '
     .         //            'index # in plate # was #. (The vertex '
     .         //            'and plate numbers in this message are '
     .         //            '1-based as well.)'                       )
               CALL ERRINT ( '#',  NV                                  )
               CALL ERRINT ( '#',  J                                   )
               CALL ERRINT ( '#',  I                                   )
               CALL ERRINT ( '#',  PLATES(J,I)                         )
               CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'                  )
               CALL CHKOUT ( 'PLTVOL'                                  )
               RETURN

            END IF

         END DO

      END DO

C
C     The volume computation below requires only a few lines of code,
C     but it might not be obvious that it works. An explanation
C     follows.
C   
C
C        Notation
C        --------
C
C        The expression
C
C           A x B
C
C        denotes the cross product of vectors A and B. The notation
C
C           < A, B >
C
C        denotes the dot product of these vectors.
C
C        The expression
C
C           ||A||
C
C        denotes the Euclidean (L2) norm of A, and
C
C           ^
C           A
C
C        is the unit-length vector aligned with A, namely
C
C           A / ||A||
C
C        provided A is non-zero.
C
C     
C        Geometric assumptions
C        ---------------------
C     
C        The algorithm used here assumes, but does not attempt to
C        verify, that the input plate model represents the boundary of
C        a spatial region. The region has an "inside" and an "outside,"
C        so it makes sense to speak of surface normal vectors that are
C        either "inward pointing" or "outward pointing."

C        The input plates are assumed to have vertex indices ordered in
C        the positive sense about their respective outward normal
C        vectors. So if a plate's vertices are A, B, C, then the
C        plate's sides are
C
C           B - A,  C - B,  A - C
C
C        and all of
C
C           ( B - A ) x ( C - B )
C           ( C - B ) x ( A - C )
C           ( A - C ) x ( B - A )
C
C        are outward normal vectors.
C
C        Each plate lies in a plane, so all points on the plate satisfy
C        an equation of the form
C
C           < P, N > = constant
C
C        where N is an outward normal vector for that plate. If the
C        plane constant is positive, we say that the outward normal
C        vector "points away" from the origin; if the constant is
C        negative, we say that the outward normal "points toward" the
C        origin. This is a linguistic shortcut for the unwieldy phrases
C        "the component of the outward normal parallel to P points away
C        from the origin" or ... (the negative counterpart phrase).
C
C
C        Relationship of volume and surface
C        ----------------------------------
C
C        Given a suitable 3-D spatial region, the divergence theorem
C        (see [1]) says that the volume integral, over that region, of
C        the divergence of a vector field F (div F) is equal to the
C        surface integral of F over the boundary of the region. So if
C        one picks a vector field F having divergence identically equal
C        to 1, the surface integral of F equals the volume of the
C        region.
C
C        We can use this fact to set up a surface integral that yields
C        the volume enclosed by a set of triangular plates.
C
C
C        A faster algorithm
C        ------------------
C
C        However, there's a more efficient approach: we can compute the
C        volume of the model by summing the signed volumes of the
C        pyramids whose bases are the plates and whose apexes coincide
C        with the origin. Plates with outward normal vectors pointing
C        away from the origin contribute positive volume increments;
C        plates with outward normal vectors pointing toward the origin
C        contribute negative volume increments.
C
C        We'll show that these two methods are equivalent, and we'll
C        use the more efficient method in our code.
C
C
C        The surface integral
C        --------------------
C
C        The field 
C
C           F( x, y, z )  =  (1/3)( x, y, z )
C
C        has divergence 1. Let the vectors
C
C           A, B, C
C
C        be the vertices of a plate, and let
C
C           N = ( B - A )  x  ( C - B )
C
C        be an outward normal of the plate. Then 
C
C           < N, X > 
C
C        is constant for all X in the plane containing the plate, so we
C        can pick the vertex A to stand in for X.
C
C        Then the surface integral of F over the plate is simply
C
C             ^
C           < N,  A/3 > * plate_area
C
C
C        Equivalence of the methods
C        --------------------------
C
C        If we show that the above integral equals the volume
C        contribution of the pyramid corresponding to the plate, the
C        validity of summing the signed pyramid volumes is established.
C
C        Below, let
C
C           S1, S2 be the plate sides (B-A), (C-B).
C
C           N be the outward normal vector S1 x S2.
C
C           plate_area be...the area of the plate.
C
C        Then
C
C           N                  =   S1 x S2
C                              =   ( B - A ) x ( C - B )
C                              =   ( B x C ) - ( A x C ) + ( A x B )
C
C        and since A is orthogonal to both
C
C           A x C
C           A x B
C
C        we have
C
C             ^
C           < N,  A >          =   <  ( B x C )/||N||,  A  >
C
C
C        Then
C
C           plate_area         =    (1/2) * ||( S1 x S2 )||
C                              =    ||N|| / 2
C
C        So the surface integral of F over the plate is
C
C                            ^
C           plate_area * ( < N, A/3 > ) 
C
C                              =  (1/6) * ||N|| 
C                                       * < (B x C)/||N||, A > 
C
C                              =  < ( B x C ), A > / 6
C
C                              =  < ( A X B ), C > / 6
C
C        On the other hand, letting 
C
C           base_area
C
C        denote the area of the triangle defined by A, B, and the
C        origin, the signed volume of the pyramid defined by A, B, and
C        C is
C     
C
C           1/3 * base_area * height   =  1/3 *  < (A x B)/2,  C  >
C
C                                      =  < ( A X B ), C > / 6
C                                   
C        This shows the signed pyramid volume is equal to the surface
C        integral of F over the plate. 
C
C
C
C     We proceed to compute the signed volume of each pyramid whose
C     base is a plate and whose vertex is the origin.
C
C     Note that PLTVOL has already been initialized to 0.D0.
C
      DO I = 1, NP
C
C        Pack the vertices of the current plate into a 3x3 matrix.
C 
         DO J = 1, 3
            
            CALL VEQU (  VRTCES( 1, PLATES(J,I) ),  M(1,J)  )

         END DO

C
C        The determinant of M gives the volume of the parallelepiped
C        spanned by the origin-vertex vectors. The corresponding
C        pyramid has volume 1/3 * base_area * height, and
C        the area of the pyramid's base is half that of base of the
C        parallelepiped. So the determinant must be divided by 6.
C
         PLTVOL = PLTVOL  +  (  DET( M ) / 6.D0  )

      END DO

C
C     No check-out required, since the routine is not checked in
C     at this point.
C
      RETURN
      END
