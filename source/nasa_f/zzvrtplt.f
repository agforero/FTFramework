C$Procedure  ZZVRTPLT  ( create vertex-plate mapping )
 
      SUBROUTINE ZZVRTPLT ( NV,    NP,     PLATES, CELLSZ, MAXLST, 
     .                      CELLS, VRTPTR, NLIST,  PLTLST          )
      IMPLICIT NONE

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Generate a data structure that represents the mapping from
C     vertices to the plates containing them.
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
C     PLATE 
C     VERTEX
C
C$ Declarations
 
      INCLUDE              'dsk02.inc'
 
      INTEGER               NV
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
      INTEGER               CELLSZ
      INTEGER               MAXLST
      INTEGER               CELLS  ( 2, CELLSZ )
      INTEGER               VRTPTR ( * )
      INTEGER               NLIST
      INTEGER               PLTLST ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NV         I   Number of vertices.
C     NP         I   Number of plates.
C     PLATES     I   Plate array.
C     CELLSZ     I   Cell array size.
C     MAXLST     I   Maximum size of output plate list.
C     CELLS     I-O  Workspace array for vertex list construction.
C     VRTPTR     O   Array of pointers, by vertex, into plate data list.
C     NLIST      O   Length of the plate data list.
C     PLTLST     O   Plate data list; first element is the number of
C                    plate IDs to follow.
C
C$ Detailed_Input
C
C     NV         Number of vertices.
C
C     NP         Number of plates.
C
C     PLATES     Plate array: this is an array of 3-tuples of 
C                vertex indices. The elements
C
C                   PLATES(1:3,I)
C
C                are the indices of the vertices of the Ith plate.
C
C
C     CELLSZ     Size of cell array: explicit workspace dimension used
C                for error checking. This value should be 
C
C                   3*NP
C
C     MAXLST     Maximum size of output plate list. This size should be
C
C                   3*NP  +  NV
C                
C
C     CELLS      workspace array used to construct the vertex-plate
C                mapping.
C
C$ Detailed_Output
C
C     CELLS      Workspace array used to construct the vertex-plate
C                mapping.
C
C     VRTPTR     Array of pointers associating plate lists with
C                vertices. The Ith element of VRTPTR contains
C                the start index in PLTLST of the list associated with
C                the Ith vertex. 
C
C                The size of this array must be at least NV.
C
C     NLIST      Length of the output plate data list PLTLST.
C
C     PLTLST     Plate data list: for each vertex, there is a count
C                of associated plates, followed by the IDs of those
C                plates. The count for vertex I is located at
C
C                   PLTLST( VRTPTR(I) )
C
C                The size of this array must be at least MAXLST.
C
C$ Parameters
C
C     See parameter declarations in 
C
C        dsk02.inc
C
C$ Exceptions
C
C     1)  If the input plate count NP is non-positive, the
C         error SPICE(BADPLATECOUNT) is signaled.
C
C     2)  If the input vertex count NV is non-positive, the
C         error SPICE(BADVERTEXCOUNT) is signaled.
C
C     3)  If the cell array size CELLSZ is less than 
C
C            3 * NP
C         
C         the error SPICE(CELLARRAYTOOSMALL) is signaled.
C
C     4)  If the plate list size MAXPLT is less than 
C
C            ( 3 * NP ) + NV
C
C         the error SPICE(PLATELISTTOOSMALL) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports the DSK type 2 spatial index
C     routines 
C
C        DSKMI2
C        DSKSI2
C
C$ Examples
C
C     See usage in DSKMI2.
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
C-    SPICELIB Version 1.0.0, 20-MAY-2016 (NJB)
C
C        Updated error handling. Now ZZVRTPLT performs precise checks
C        on cell and plate array size inputs.
C
C
C        22-DEC-2015 (NJB)
C
C           Renamed routine from VRTCOM to ZZVRTPLT.
C
C           CAUTION: argument list change! Removed input workspace
C           array for pointer construction. Array was named PNTRS.
C
C           Added error checks for input counts and sizes.
C
C           Now calls updated 3.0.0 version of ZZUNTNGL.
C
C        03-MAY-2014 (NJB)
C
C           Now calls ZZ* versions of supporting routines.
C
C           Changed argument list to include sizes of arrays.
C           Changed error handling to make use of array sizes.
C           Changed call to UNTNGL to accommodate argument
C           list change in that routine. Updated header I/O
C           descriptions.
C
C        12-JUL-2011 (NJB)
C
C           Argument list change: the input arrays V1, V2, V3 were
C           replaced with the input array PLATES.
C
C        04-MAY-2010 (NJB)
C
C           Now accepts input workspace arrays PNTRS and CELLS.
C           Changed INCLUDE file to dsk02.inc.
C
C        08-OCT-2009 (NJB)
C
C           Re-ordered header sections.
C
C        04-FEB-1999 (JAB)
C
C-&
 
C$ Index_Entries
C
C     create vertex to plate map
C
C-&
  
C
C     SPICELIB functions
C     
      LOGICAL               FAILED

C
C     Local variables.
C
 
      LOGICAL               RETURN
 
      INTEGER               I
      INTEGER               J
      INTEGER               NCELL
 

      SAVE

C
C     Standard SPICELIB error handling.
C
 
      IF ( RETURN () ) THEN
         RETURN
      END IF
 
      CALL CHKIN  ( 'ZZVRTPLT' )
 

      IF ( NV .LT. 1 ) THEN

         CALL SETMSG ( 'Vertex count NV = #; count must be positive.' 
     .   //            'be positive.'                                )
         CALL ERRINT ( '#',  NV                                      )
         CALL SIGERR ( 'SPICE(BADVERTEXCOUNT)'                       )
         CALL CHKOUT ( 'ZZVRTPLT'                                    )
         RETURN  
         
      END IF

      IF ( NP .LT. 1 ) THEN

         CALL SETMSG ( 'Plate count NP = #; count must be positive.' 
     .   //            'be positive.'                                )
         CALL ERRINT ( '#',  NP                                      )
         CALL SIGERR ( 'SPICE(BADPLATECOUNT)'                        )
         CALL CHKOUT ( 'ZZVRTPLT'                                    )
         RETURN  
         
      END IF

      IF ( CELLSZ .LT. 3*NP ) THEN

         CALL SETMSG ( 'Cell array size CELLSZ = #; size must be ' 
     .   //            '>= 3*NP. NP is the plate count #.'         )
         CALL ERRINT ( '#',  CELLSZ                                )
         CALL ERRINT ( '#',  NP                                    )
         CALL SIGERR ( 'SPICE(CELLARRAYTOOSMALL)'                  )
         CALL CHKOUT ( 'ZZVRTPLT'                                  )
         RETURN  
         
      END IF

      IF ( MAXLST .LT. (NV+3*NP) ) THEN

         CALL SETMSG ( 'Plate list array size MAXPLT = #; size ' 
     .   //            'must be >= 3*NP + NV, which is #. (NV = '
     .   //            'vertex count, NP = plate count.)'        )
         CALL ERRINT ( '#',  MAXLST                              )
         CALL ERRINT ( '#',  3*NP + NV                           )
         CALL SIGERR ( 'SPICE(PLATELISTTOOSMALL)'                )
         CALL CHKOUT ( 'ZZVRTPLT'                                )
         RETURN  
         
      END IF

C
C     Initialize pointer and cell structure.
C
      CALL ZZINILNK ( NV, CELLSZ, NCELL, VRTPTR, CELLS )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZVRTPLT' )
         RETURN
      END IF

C
C     Loop over all plate IDS. Add each plate/vertex 
C     combination to the linked list.
C
      DO I = 1, NP

         DO J = 1, 3
C
C           AVAL = PLATES(J,I), vertex J of plate ID I. 
C           BVAL = I, plate ID value I.
C
            CALL ZZADDLNK ( PLATES(J,I), I,     NP,   CELLSZ, 
     .                      VRTPTR,      NCELL, CELLS        )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZVRTPLT' )
               RETURN
            END IF

         END DO

      END DO

C
C     Generate two linked lists mapping vertex ID to the plates
C     including that vertex as a member.
C
C     VRTPTR: An array, indexed by vertex ID. For an array element, 
C             VRTPTR(VERT_ID), greater than zero, the value identifies 
C             an index in PLTLST, the value of that PLTLST array
C             element  equaling the number of plates that include 
C             the vertex specified by the ID as a member. The  
C             condition VRTPTR(VERT_ID) = -1 indicates a bug.
C
C     PLTLST: An array, indexed by the positive entries in VRTPTR.
C             The element, N, identified by a VRTPTR value describes
C             the number of plates of which the vertex is a member. 
C
C                 N = PLTLST( VRTPTR(VERT_ID) )
C
C             The N elements following PLTLST( VRTPTR(VERT_ID) ),
C             contain the IDs of those plates which have the vertex
C             as a member.
C
      CALL ZZUNTNGL ( NV,     CELLSZ, CELLS,
     .                MAXLST, VRTPTR, NLIST, PLTLST )
 
C
C     Standard SPICE error handling.
C
      CALL CHKOUT ( 'ZZVRTPLT' )
      RETURN
      END
