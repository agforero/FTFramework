C$Procedure ZZCAPPLT ( Make polar cap plates )

      SUBROUTINE ZZCAPPLT ( NCOLS,  NORTH,  WRAP,  
     .                      BASIDX, POLIDX, NP,   PLATES )
 
C$ Abstract
C
C     Generate a set of plates constituting a polar cap.
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

      INTEGER               NCOLS
      LOGICAL               NORTH
      LOGICAL               WRAP
      INTEGER               BASIDX
      INTEGER               POLIDX
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NCOLS      I   Number of columns of data.
C     NORTH      I   Pole selector: .TRUE. if pole is north.
C     WRAP       I   Longitude wrap flag: .TRUE. enables wrapping.
C     BASIDX     I   Base index of vertices in row next to pole.
C     POLIDX     I   Index of pole vertex.
C     NP         O   Number of plates.
C     PLATES     O   Plate array.
C
C$ Detailed_Input
C
C     NCOLS          is the number of columns of data in the vertex
C                    grid for which a cap is to be constructed.
C
C     NORTH          is a logical flag indicating whether the cap
C                    is to be constructed for the north or south
C                    pole. A value of .TRUE. indicates north; .FALSE.
C                    indicates south.
C 
C     WRAP           is a logical flag indicating whether longitude
C                    wrapping is to be performed. Longitude wrapping
C                    entails creating an extra plate to connect the
C                    east and west edges of the polar cap.
C 
C     BASIDX         is the base index of the set of vertices 
C                    constituting the vertex row adjacent to the 
C                    designated pole. The indices of the vertices
C                    constituting this row range from 
C
C                       BASIDX + 1  
C
C                    to 
C               
C                       BASIDX + NCOLS                 
C
C                    Regardless of the format of the input data file,
C                    the vertex grid is presumed to have row major
C                    order, and have top-down, left-to-right indexing.
C                    
C     POLIDX         is the index of the pole vertex.
C                     
C$ Detailed_Output
C
C     NP             is the output plate count.
C
C     PLATES         is a set of plates constituting a "cap" for the
C                    specified pole.
C
C                    Each plate consists of three vertex indices. The
C                    elements of PLATES are stored in row major order,
C                    regardless of the organization of the input data
C                    set.
C
C                    The dimension of PLATES is ( 3, NP ).
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the input column count is not at least two, the error
C        SPICE(INVALIDDIMENSION) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is used to extend plate sets to cover
C     the latitude range -90 : 90 degrees. Many data sets
C     cover latitude ranges that exclude one or both poles.
C     It can be convenient to fill in the gaps at the poles
C     to ensure that certain operations, such as ray-surface
C     intercept computations performed to create latitude-
C     longitude grids, succeed for the entire surface.
C
C     Vertices in the row adjacent to the pole are assumed to
C     have left-to-right order. Vertex order determines the
C     direction of the output plates' outward normal vectors.
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
      INTEGER               TL
      INTEGER               TR
      INTEGER               UB


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN( 'ZZCAPPLT' )
C
C     Check column dimensions.     
C     
      IF ( NCOLS .LT. 2 ) THEN

         CALL SETMSG( 'Grid must have at least two columns '
     .   //           'but NCOLSS is #.'                    )
         CALL ERRINT( '#',  NCOLS                           )
         CALL SIGERR( 'SPICE(INVALIDDIMENSION)'             )
         CALL CHKOUT( 'ZZCAPPLT'                              )
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

      DO I = 1, UB

         IF ( NORTH ) THEN
C
C           Create plates for a north polar cap.
C
C           Set the vertex index of the north pole.
C
            TL = POLIDX
C
C           Longitude increases with increasing column index.
C               
            IF (  WRAP .AND. ( I .EQ. UB ) ) THEN
C
C              Form a plate by connecting the right edge
C              of the surface to the left edge.
C
               BL = BASIDX + UB
               BR = BASIDX + 1                  

            ELSE

               BL = BASIDX + I
               BR = BL     + 1

            END IF

C
C           Create the current plate.
C
            NP = NP + 1

            PLATES( 1, NP ) = TL
            PLATES( 2, NP ) = BL
            PLATES( 3, NP ) = BR

         ELSE
C
C           Create plates for a south polar cap.
C
C           Set the vertex index of the south pole.
C
            BL = POLIDX

C
C           Longitude increases with increasing column index.
C               
            IF (  WRAP .AND. ( I .EQ. UB ) ) THEN
C
C              Form a plate by connecting the right edge
C              of the surface to the left edge.
C
               TL = BASIDX + UB
               TR = BASIDX + 1                  

            ELSE

               TL = BASIDX + I
               TR = TL     + 1

            END IF

C
C           Create the current plate.
C
            NP = NP + 1

            PLATES( 1, NP ) = BL
            PLATES( 2, NP ) = TR
            PLATES( 3, NP ) = TL

         END IF

      END DO

C
C     The plate and vertex counts and arrays have been
C     assigned.
C
      CALL CHKOUT( 'ZZCAPPLT' )
      RETURN
      END 
