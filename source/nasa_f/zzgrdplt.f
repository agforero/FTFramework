C$Procedure ZZGRDPLT ( Create grid of plates )

      SUBROUTINE ZZGRDPLT ( NROWS, NCOLS, WRAP, NP, PLATES )

C$ Abstract
C
C     Generate a set of plates covering the surface spanned by a
C     rectangular grid of vertices.
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      INTEGER               NROWS
      INTEGER               NCOLS
      LOGICAL               WRAP
      INTEGER               NP
      INTEGER               PLATES ( 3, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NROWS      I   Number of rows of data.
C     NCOLS      I   Number of columns of data.
C     WRAP       I   Longitude wrap flag: .TRUE. enables wrapping.
C     NP         O   Number of plates.
C     PLATES     O   Plate array.
C
C$ Detailed_Input
C
C     NROWS,
C     NCOLS          are, respectively, the numbers of rows and columns
C                    of data in the vertex grid for which a plate set
C                    is to be constructed.
C
C     WRAP           is a logical flag indicating whether longitude
C                    wrapping is to be performed. Longitude wrapping
C                    entails creating an extra column of plates to
C                    connect the east and west edges of the vertex
C                    grid.
C                     
C$ Detailed_Output
C
C     NP             is the output plate count.
C
C     PLATES         is a set of plates covering the region spanned 
C                    by the input vertex grid. Every region having
C                    four adjacent vertices at its corners is covered
C                    by two plates.
C
C                    Each plate consists of three vertex indices. 
C
C                    The indices of the vertices constituting the grid
C                    are presumed to range from
C
C                       1  to  NROWS*NCOLS
C
C                    Regardless of the format of the input data file,
C                    the vertex grid is presumed to have row major
C                    order, and have top-down, left-to-right indexing.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If either the input column count or input row count
C        is not at least two, the error SPICE(INVALIDDIMENSION) is
C        signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     To avoid portability problems that can arise from
C     basing plate-vertex association on computed quantities, 
C     all plates are constructed according to the following, 
C     fixed pattern:
C
C
C        vertex(I,J)      vertex(I,J+1)
C            +----------------+
C            |\               |
C            | \              |
C            |  \             |
C            |   \            |
C            |    \           |
C            |     \          |
C            |      \         |
C            |       \        |
C            |        \       |
C            |         \      |
C            |          \     |
C            |           \    |
C            |            \   |
C            |             \  |
C            |              \ |
C            |               \|
C            +----------------+
C        vertex(I+1,J)      vertex(I+1,J+1)
C
C
C        
C$ Examples
C
C     See usage in the MKDSK routine MKGRID.
C
C$ Restrictions
C
C     1) For use only within program MKDSK.
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
C-    MKDSK Version 1.0.1, 30-APR-2014 (NJB)
C
C        Corrected some comment typos.
C
C-    MKDSK Version 1.0.0, 25-JUL-2011 (NJB)
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      INTEGER               BL
      INTEGER               BR
      INTEGER               I
      INTEGER               J
      INTEGER               TL
      INTEGER               TR
      INTEGER               UB

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN( 'ZZGRDPLT' )

C
C     Check row and column dimensions.     
C     
      IF ( NROWS .LT. 2 ) THEN

         CALL SETMSG( 'Grid must have at least two rows '
     .   //           'but NROWS is #.'                  )
         CALL ERRINT( '#',  NROWS                        )
         CALL SIGERR( 'SPICE(INVALIDDIMENSION)'          )
         CALL CHKOUT( 'ZZGRDPLT'                           )
         RETURN

      END IF

      IF ( NCOLS .LT. 2 ) THEN

         CALL SETMSG( 'Grid must have at least two columns '
     .   //           'but NCOLSS is #.'                    )
         CALL ERRINT( '#',  NCOLS                           )
         CALL SIGERR( 'SPICE(INVALIDDIMENSION)'             )
         CALL CHKOUT( 'ZZGRDPLT'                              )
         RETURN

      END IF
 
C
C     Set the upper bound on the column loop. If longitude
C     wrapping is turned on, the final column connects to 
C     the first column.
C
      IF ( WRAP ) THEN

         UB = NCOLS
      ELSE
         UB = NCOLS - 1
      END IF

C
C     Connect the vertices to generate plates.
C
      NP = 0

      DO I = 1, NROWS-1

         DO J = 1, UB
C
C           Since the vertices are stored in a 3xN
C           array, it will be convenient to use aliases for
C           their indices---this cuts down on use of complicated
C           index expressions. For each square defined by 
C           four neighboring vertices, we'll call the vertices
C           
C              TL   "top left"
C              TR   "top right"
C              BL   "bottom left"
C              BR   "bottom right"
C
C           Recall that the input pixel grid has dimensions
C    
C              NROWS x NCOLS
C           
C           The top row is at the highest latitude.
C
C           The leftmost column corresponds to the west
C           boundary of the region.
C
C           The top left vertex corresponds to pixel (I,J).
C
            IF (  WRAP  .AND.  ( J .EQ. UB )  ) THEN
C
C              Connect the right edge of the grid to the left edge.
C                    
               TL = I  * NCOLS
               TR = TL - NCOLS + 1
               BL = TL + NCOLS
               BR = TR + NCOLS

            ELSE
C
C              This is the normal case: the column at index
C              J is connected to the column at index J+1.
C
               TL = ( (I-1) * NCOLS )  +  J
               TR = TL + 1
               BL = TL + NCOLS
               BR = BL + 1

            END IF

C
C           For each square defined by neighboring pixel centers,
C           we must represent the corresponding surface by a pair
C           of plates. We have two choices for the diagonal 
C           common edge connecting these plates: descending or
C           ascending to the right. 
C
C           We choose the descending diagonal. 
C
C           The vertex assignment must be positively
C           oriented about the outward normal direction.
C
            NP = NP + 1

            PLATES( 1, NP ) = BL
            PLATES( 2, NP ) = BR
            PLATES( 3, NP ) = TL

            NP = NP + 1

            PLATES( 1, NP ) = TL
            PLATES( 2, NP ) = BR
            PLATES( 3, NP ) = TR

         END DO

      END DO

C
C     The plate and vertex counts and arrays have been
C     assigned.
C
      CALL CHKOUT( 'ZZGRDPLT' )
      RETURN
      END 
