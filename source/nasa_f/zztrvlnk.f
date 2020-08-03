C$Procedure ZZTRVLNK ( Traverse AB cell linked-list )
 
      SUBROUTINE ZZTRVLNK ( AVAL,  MAXA, PNTRS, CELLSZ, 
     .                      CELLS, MAXB, NB,    BLIST  )
      IMPLICIT NONE

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Traverse an AB pool linked list, searching for all cells that are
C     associated with a specified A-value. Store the B-values of these
C     cells in an array.
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
C     linked list
C
C$ Declarations
  
      INTEGER               AVAL
      INTEGER               MAXA
      INTEGER               PNTRS  ( * )
      INTEGER               CELLSZ
      INTEGER               CELLS  ( 2, * )
      INTEGER               MAXB
      INTEGER               NB
      INTEGER               BLIST  ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     AVAL       I   A-value index associated with pool entry.
C     MAXA       I   Size of PNTRS array.
C     PNTRS      I   A-value pointer array.
C     CELLSZ     I   Number of cells in cell array.
C     CELLS      I   Cell array.
C     MAXB       I   Size of BLIST array.
C     NB         O   Number of B values found.
C     BLIST      O   List of B values.
C
C$ Detailed_Input
C
C     AVAL       is an index associated with a pool entry. AVAL
C                is the index of the Ith value of a finite set "A."
C                AVAL is an index into the PNTRS array. 
C
C     MAXA       Size of PNTRS array.
C
C     PNTRS      is an array of integers acting as pointers into the
C                cell array. The element at index I of PNTRS points to
C                the head of a linked list in CELLS that associates
C                B-values with the index I. Here I is the "A value."
C
C                If there are no values associated with the Ith A value,
C                PNTRS(I) is null (-1).
C
C
C     CELLSZ     is the number of cells in the cell array. The array
C                has dimensions
C
C                   (2, CELLSZ)
C
C     CELLS      is the input cell array. The element
C
C                   CELLS ( 1, PNTRS(AVAL) )
C
C                contains a B-value associated with AVAL. The element
C
C                   CELLS ( 2, PNTRS(AVAL) )
C
C                is a pointer to the next node of the list associated
C                with AVAL, or is null (-1) if the list contains only
C                one element. The pointer element of the last list node
C                for AVAL is null.
C
C     MAXB       is the size of the output B-value list.
C
C$ Detailed_Output
C
C     NB         is the number of nodes in the list associated with
C                AVAL.
C
C     BLIST      is an array containing the B-values of the nodes
C                associated with AVAL. There is one entry in BLIST for
C                each node in the list associated with AVAL; the list
C                may contain duplicates.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the array index AVAL is less than 1 or exceeds the
C         declared size of the PNTR array, the error
C         SPICE(INDEXOUTOFRANGE) is signaled. Occurrence of this error
C         may result from the fine voxel scale being too small for the
C         data.
C
C     2)  If the size of the output BLIST array is non-positive, the
C         error SPICE(INVALIDSIZE) is signaled.
C
C     3)  If an element of the PNTRS array is larger than the length
C         CELLSZ of the cell array, the error SPICE(POINTEROUTOFRANGE)
C         is signaled.
C
C     4)  If the number of output list entries exceeds the size MAXB of
C         the output list array, the error SPICE(BARRAYTOOSMALL) is
C         signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The AB list routines support manipulation of a data structure that
C     represents a one-to-many mapping from a set A to a set B. Each 
C     element of set A is mapped to a set of values from B, where the
C     range values are represented by a linked list of entries in the
C     cell array. 
C
C     The elements of set B associated with the Ith element of A are
C     stored in a linked list of cells starting with the cell
C     consisting of elements
C
C        CELLS(*,PNTRS(I))
C
C     The set A is really an abstraction: it's just some finite set
C     with size in the range 1:MAXA. The input AVAL is not a member of
C     A but rather just an index into A. For consistency with existing
C     code, the name AVAL has been retained.
C
C     The fact that B values are stored in linked lists enables a
C     program to store entries for A-B associations in random order.
C     For example, in the process of constructing DSK type 2 segments,
C     the of non-empty fine voxels intersected by a given plate can be
C     computed; then an entry representing a voxel-plate association
C     can be made for each fine voxel in the set. In this case, the 
C     set "A" is the set of fine voxels belonging to non-empty coarse
C     voxels.    
C
C     This routine supports creation of DSK type 2 segments. It is used
C     for creation of both voxel-plate association data structures and
C     of vertex-plate association data structures.
C
C$ Examples
C
C     See usage in the routine ZZUNTNGL.
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
C-    SPICELIB Version 1.0.0, 23-AUG-2016 (NJB)
C
C      07-JAN-2016 (NJB)
C
C        Strengthened check on invalid pointer. Renamed argument MAXCEL
C        to CELLSZ. Updated header comments.
C
C      26-MAR-2015 (NJB)
C
C        Updated long error messages to improve accuracy.
C
C        Header update: added description of MAXB input
C        argument. Updated description of MAXA input 
C        argument. Filled out Exceptions section.
C
C      30-APR-2014 (NJB)
C
C        Changed argument list: each array argument now
C        has an associated argument giving its size. 
C        Updated error checking to check array indices
C        against array sizes.
C
C      08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C      11-JUN-2004 (EDW)
C
C        Added check on AVAL to ensure it's smaller
C        than the size of MAXA.
C
C      26-AUG-2002 (BVS)
C
C        Replaced WRITE with normal error reporting calls.
C
C      03-FEB-1999 (JAB)
C
C
C-&
 
C$ Index_Entries
C
C     AB cells
C
C-&
 
      LOGICAL               RETURN
 
      INTEGER               I
      INTEGER               PTR
 
C
C     Standard SPICE error handling.
C
 
      IF ( RETURN() ) THEN
         RETURN
      END IF
 
      CALL CHKIN  ( 'ZZTRVLNK' )

      IF ( ( AVAL .LT. 1 ) .OR. ( AVAL .GT. MAXA ) ) THEN

         CALL SETMSG ( 'Index AVAL is out of range. '
     .   //            'Index = #1. Valid range = 1:#2.' )
         CALL ERRINT ( '#1', AVAL                        )
         CALL ERRINT ( '#2', MAXA                        )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'          )
         CALL CHKOUT ( 'ZZTRVLNK'                        )
         RETURN

      END IF

      IF ( MAXB .LT. 1 ) THEN

         CALL SETMSG ( 'Maximum output list size MAXB is invalid. '
     .   //            'MAXB = #1.'                                )
         CALL ERRINT ( '#1', MAXB                                  )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'                        )
         CALL CHKOUT ( 'ZZTRVLNK'                                  )
         RETURN

      END IF

 
      NB       = 0
      BLIST(1) = 0

      PTR = PNTRS(AVAL)
 
      DO WHILE ( PTR .NE. -1 )

         IF (      ( PTR .LT. -1     ) 
     .        .OR. ( PTR .EQ. 0      ) 
     .        .OR. ( PTR .GT. CELLSZ )  ) THEN

            CALL SETMSG ( 'Value in PNTRS array is not a '
     .      //            'valid index in the cell array.'
     .      //            'Value = #1. Array size = #2.'    )
            CALL ERRINT ( '#1', PTR                         )
            CALL ERRINT ( '#2', CELLSZ                      )
            CALL SIGERR ( 'SPICE(POINTEROUTOFRANGE)'        )
            CALL CHKOUT ( 'ZZTRVLNK'                        )
            RETURN

         END IF

  
         NB = NB + 1

         IF ( NB .GT. MAXB ) THEN

            CALL SETMSG ( 'Output value count is larger than '
     .      //            'B-list array room. Count = #1. Output '
     .      //            'array room = #2. Input pointer index '
     .      //            'was #3. Input pointer list size was #4. ' 
     .      //            'Last pointer was #5. Cell size was #6.'  ) 
            CALL ERRINT ( '#1', NB                                  )
            CALL ERRINT ( '#2', MAXB                                )
            CALL ERRINT ( '#3', AVAL                                )
            CALL ERRINT ( '#4', MAXA                                )
            CALL ERRINT ( '#5', PTR                                 )
            CALL ERRINT ( '#6', CELLSZ                              )
            CALL SIGERR ( 'SPICE(BARRAYTOOSMALL)'                   )
            CALL CHKOUT ( 'ZZTRVLNK'                                )
            RETURN

         END IF


         BLIST(NB) = CELLS(1,PTR)
         I         = CELLS(2,PTR)
         PTR       = I
  
      END DO
 
C
C     Standard SPICE error handling.
C
      CALL CHKOUT ( 'ZZTRVLNK' ) 
      RETURN
      END
