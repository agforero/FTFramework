C$Procedure PLTAR ( Compute area of plate set )
 
      DOUBLE PRECISION FUNCTION PLTAR ( NV, VRTCES, NP, PLATES )

C$ Abstract
C
C     Compute the total area of a collection of triangular plates.
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
C     The function returns the total area of the set of plates.
C
C$ Detailed_Input
C
C     NV             is the number of vertices comprising the plate
C                    set.
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
C                    plate set.
C
C     PLATES         is an array containing 3-tuples of integers
C                    representing the set of plates. The elements of
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
C     The function returns the total area of the input set of plates.
C     Each plate contributes the area of the triangle defined by the
C     plate's vertices.
C
C     If the components of the vertex array have length unit L, then the
C     output area has units
C
C         2
C        L
C                              
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the number of plates is less than 0, the error
C        SPICE(BADPLATECOUNT) is signaled.
C
C     2) If the number of plates is positive and the number of vertices
C        is less than 3, the error SPICE(TOOFEWVERTICES) is signaled.
C
C     3) If any plate contains a vertex index outside of the range
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
C     This routine computes the total area of a set of triangular
C     plates. The plates need not define a closed surface.
C     
C     Examples of valid plate sets:
C
C        Tetrahedron
C        Box
C        Tiled ellipsoid
C        Tiled ellipsoid with one plate removed
C        Two disjoint boxes
C        Two boxes with intersection having positive volume
C        Single plate
C        Empty plate set
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
C     1) Compute the area of the pyramid defined by the four
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
C     C     Compute the area of a plate model representing the pyramid
C     C     with one vertex at the origin and the other vertices
C     C     coinciding with the standard basis vectors.
C     C
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      PLTAR
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
C           DOUBLE PRECISION      AREA
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
C           AREA = PLTAR ( NVERT, VRTCES, NPLATE, PLATES )
C
C           WRITE (*,*) 'Expected area   =    (3 + SQRT(3)) / 2'
C           WRITE (*,*) '                =    0.23660254037844384E+01'
C           WRITE (*,*) 'Computed area   = ', AREA
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran platform,
C     the output was:
C
C
C        Expected area   =    (3 + SQRT(3)) / 2
C                        =    0.23660254037844384E+01
C        Computed area   =    2.3660254037844384
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
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 21-OCT-2016 (NJB)
C
C        Original version 25-MAR-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compute plate model area
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      VNORM
      LOGICAL               RETURN


C
C     Local variables
C     
      DOUBLE PRECISION      CP    ( 3 )
      DOUBLE PRECISION      EDGE1 ( 3 )
      DOUBLE PRECISION      EDGE2 ( 3 )

      INTEGER               I
      INTEGER               J


C
C     The function must have an initial value.
C
      PLTAR  =  0.D0

C
C     This routine uses discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Check the vertex and plate counts.
C
      IF ( NP .LT. 0 ) THEN

         CALL CHKIN  ( 'PLTAR'                                        )
         CALL SETMSG ( 'Plate count must be non-negative but NP = #.' )
         CALL ERRINT ( '#',  NP                                       )
         CALL SIGERR ( 'SPICE(BADPLATECOUNT)'                         )
         CALL CHKOUT ( 'PLTAR'                                        )
         RETURN

      END IF


      IF ( NP .EQ. 0 ) THEN
C
C        The area has already been set to zero.
C   
         RETURN

      END IF
        

      IF ( NV .LT. 3 ) THEN

         CALL CHKIN  ( 'PLTAR'                                       )
         CALL SETMSG ( 'At least 3 vertices are needed, but NV = #.' )
         CALL ERRINT ( '#',  NV                                      )
         CALL SIGERR ( 'SPICE(TOOFEWVERTICES)'                       )
         CALL CHKOUT ( 'PLTAR'                                       )
         RETURN

      END IF

C
C     Make sure the vertex indices are in the range [1, NV].
C
      DO I = 1, NP

         DO J = 1, 3

            IF (      ( PLATES(J,I) .LT. 1  ) 
     .           .OR. ( PLATES(J,I) .GT. NV )  ) THEN

               CALL CHKIN  ( 'PLTAR'                                   )
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
               CALL CHKOUT ( 'PLTAR'                                   )
               RETURN

            END IF

         END DO

      END DO


      DO I = 1, NP
C
C        Take the cross product of two edges of the Ith plate.
C
         CALL VSUB ( VRTCES( 1, PLATES(2,I) ),
     .               VRTCES( 1, PLATES(1,I) ),   EDGE1 )

         CALL VSUB ( VRTCES( 1, PLATES(3,I) ),
     .               VRTCES( 1, PLATES(2,I) ),   EDGE2 )

         CALL VCRSS ( EDGE1, EDGE2, CP )

C
C        The plate area is 1/2 of the magnitude of the
C        cross product.
C        
         PLTAR = PLTAR  +  ( 0.5D0 * VNORM(CP) )

      END DO

C
C     No check-out required, since the routine is not checked in
C     at this point.
C
      RETURN
      END
