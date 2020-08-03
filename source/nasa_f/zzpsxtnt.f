C$Procedure    ZZPSXTNT ( MKDSK, compute plate set extents )
 
      SUBROUTINE ZZPSXTNT ( NV, VERTS, NP, PLATES, EXTENT, AVPLEX )
  
C$ Abstract
C
C     Compute plate set extents and average plate extent.
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
C     DSK
C     FILES
C
C$ Declarations


      IMPLICIT NONE

      INCLUDE 'mkdsk02.inc'

      INTEGER               NV
      DOUBLE PRECISION      VERTS  ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
      DOUBLE PRECISION      EXTENT ( 2, 3 )
      DOUBLE PRECISION      AVPLEX

 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NV         I   Vertex count.
C     VERTS      I   Vertices.
C     NP         I   Plate count.
C     PLATES     I   Plates.
C     EXTENT     O   Extents of the vertex set.
C     AVPLEX     O   Average plate extent.
C
C$ Detailed_Input
C
C     NV             is the number of vertices in the input vertex
C                    array.
C
C     VERTS          is an array of vertices associated with plates.
C   
C     NP             is the number of plates in the input plate array.
C
C     PLATES         is an array of triangular plates. The element
C
C                       PLATES(J,I) 
C
C                    is the index in the VERTS array of the Jth vertex
C                    of the Ith plate. The vertex indices are 1-based.
C
C$ Detailed_Output
C
C     EXTENT         is the extent of the input vertex set. See
C                    Particulars.
C
C     AVPLEX         is the average extent of the plates in the 
C                    plate set. See Particulars.
C
C$ Parameters
C
C     See mkdsk02.inc
C
C$ Exceptions
C
C     This routine is meant to be operated in RETURN SPICE error
C     handling mode. The caller is expected to delete the DSK file if
C     an error occurs during file creation.
C
C
C     1)  If the vertex count is less than 3 or greater than MAXVRT,
C         the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2)  If the plate count is less than 1 or greater than MAXPLT,
C         the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The "extent" of a set of vertices is the set of extrema
C     of the coordinates of those vertices.
C
C     The extent of a plate along a coordinate axis is the length of
C     the orthogonal projection of the plate onto that axis.
C
C     The average "extent" of a set of plates is the average of the 
C     extents of the plates along the coordinate axes.
C
C$ Examples
C
C     See usage in MKDSK.
C
C$ Restrictions
C
C     This routine should be called only from within MKDSK.
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
C-    MKDSK Version 1.0.0, 18-FEB-2017 (NJB)
C
C-&
 
C$ Index_Entries
C
C     determine vertex set extents and average plate extent
C
C-&
 

 
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      DPMIN

      LOGICAL               RETURN

C
C     Local parameters
C

C
C     Local variables
C
      DOUBLE PRECISION      PLTBDS ( 2 )
      DOUBLE PRECISION      SUM

      INTEGER               I
      INTEGER               J
      INTEGER               K


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZPSXTNT' )

C
C     Reject invalid plate and vertex counts.
C
      IF (  ( NV .LT. 3 ) .OR. ( NV .GT. MAXVRT )  ) THEN

         CALL SETMSG ( 'Vertex count NV = #; count must ' 
     .   //            'be in the range 3:#.'            )
         CALL ERRINT ( '#',  NV                          )
         CALL ERRINT ( '#',  MAXVRT                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'          )
         CALL CHKOUT ( 'ZZPSXTNT'                        )
         RETURN
         
      END IF

      IF (  ( NP .LT. 1 ) .OR. ( NP .GT. MAXPLT )  ) THEN

         CALL SETMSG ( 'Plate count NP = #; count must ' 
     .   //            'be in the range 1:#.'            )
         CALL ERRINT ( '#',  NP                          )
         CALL ERRINT ( '#',  MAXPLT                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'          )
         CALL CHKOUT ( 'ZZPSXTNT'                        )
         RETURN
         
      END IF      

C
C     Compute extents of the plate set. These depend only on
C     the vertex set.
C
      DO I = 1, 3

         EXTENT(1,I) = DPMAX()
         EXTENT(2,I) = DPMIN()

      END DO


      DO I = 1, NV

         DO J = 1, 3

            EXTENT(1,J) = MIN( EXTENT(1,J), VERTS(J,I) )
            EXTENT(2,J) = MAX( EXTENT(2,J), VERTS(J,I) )
            
         END DO

      END DO

C
C     Compute plate extents in all 3 dimensions; compute
C     the average of all extents.     
C
      SUM = 0.D0

      DO I = 1, NP
C
C        I is the current plate index.
C
         DO J = 1, 3
C
C           J is the current coordinate index.
C
            PLTBDS(1) = DPMAX()
            PLTBDS(2) = DPMIN()

            DO K = 1, 3
C
C              K is the current vertex index for the current plate.
C              Account for this vertex's contributions to the plate's
C              extent in the direction of coordinate J.
C
               PLTBDS(1) = MIN( PLTBDS(1),  VERTS( J, PLATES(K,I) ) )
               PLTBDS(2) = MAX( PLTBDS(2),  VERTS( J, PLATES(K,I) ) )

            END DO
C
C           Add the current plate's extent in the Jth coordinate.
C
            SUM = SUM + ( PLTBDS(2) - PLTBDS(1) )

         END DO

      END DO

C
C     We've added three extents per plate, so we have 3*NP in all.
C     Note we've ensured NP is positive.
C
      AVPLEX = SUM / ( 3*NP )

      CALL CHKOUT ( 'ZZPSXTNT' )
      RETURN
      END
