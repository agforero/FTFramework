C$Procedure ZZUNTNGL ( Untangle an AB linked-cell list  )

      SUBROUTINE ZZUNTNGL ( NPTR, MAXCEL, CELLS, 
     .                      MAXB, PNTRS,  NOUT,  OUTLST )
      IMPLICIT NONE

C$ Abstract
C
C     Untangle an unsorted AB linked-cell list into an A-index and
C     associated B-list.
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
C     AB linked list
C
C$ Declarations

      INTEGER               NPTR
      INTEGER               MAXCEL
      INTEGER               CELLS  ( 2, * )
      INTEGER               MAXB
      INTEGER               PNTRS  ( * )
      INTEGER               NOUT
      INTEGER               OUTLST ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NPTR       I   Length of A-list.
C     MAXCEL     I   Size of cell list.
C     CELLS      I   AB cell linked list.
C     MAXB       I   Maximum allowed dimension of B list.
C     PNTRS     I-O  Pointer array.
C     NOUT       O   Length of B-list.
C     OUTLST     O   B-list.
C
C$ Detailed_Input
C
C     NPTR       is the size of the input pointer array INPTR. 
C
C     MAXCEL     is the number of cells in the CELLS array. That
C                array has dimensions
C
C                   (2, MAXCEL)
C
C     CELLS      is an array containing linked lists of "B" values.
C                The elements at indices
C
C                   (1,K)
C                   (2,K)
C
C                are, respectively, a "B" value and a pointer to
C                another element of the array. Null pointers are 
C                indicated by the value -1.
C 
C                Each linked list in CELLS contains a set of "B"
C                values associated with a particular "A" value. 
C                The input pointer list maps "A" values to the head
C                of the associated list in CELLS.
C
C     MAXB       is the maximum number of elements in the output 
C                array OUTLST.
C
C     PNTRS      is an array of pointers from A values to
C                linked lists of "B" entries in the CELLS array.
C
C$ Detailed_Output
C
C     PNTRS      is an array of pointers that map "A" values to 
C                their associated "B" values in the array OUTLST.
C                Null pointers are indicated by the value -1.
C                
C                PNTRS must be declared with size at least NPTR.
C
C     NOUT       is the number of elements in the output list OUTLST.
C
C     OUTLST     is an array containing lists of "B" values. The Kth
C                element of PNTRS is the index of the start of a list
C                in OUTLST of "B" values associated with the "A" value
C                having index K. The list of values associated with an
C                "A" value starts with a count and is followed by
C                the values themselves.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  The routine signals the error SPICE(BARRAYTOOSMALL) if the
C         entire set of B-value lists from the input cell array,
C         including the counts for each list, array cannot be stored in
C         the OUTLST array.
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
C     See usage in ZZMKSPIN.
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
C        Argument list change: removed INPTR; renamed
C        OUTPTR to PNTRS. PNTRS is an in-out argument.
C        Updated header.
C
C      30-APR-2014 (NJB)
C
C        Changed argument list to accommodate 
C        better error checking. Added error checks
C        for array overflows. Changed calls to
C        ADDLNK and TRVLNK. Updated header I/O
C        descriptions.     
C     
C      08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C      16-JUL-2004 (EDW)
C
C        Added check on NPTR to ensure it's smaller
C        than the size of the pointer arrays.
C
C        Removed use of BVAL array, TRVLNK now directly
C        writes to OUTLST array.
C
C      03-FEB-1999 (JAB)
C
C
C-&

C$ Index_Entries
C
C     AB linked list
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN


C
C     Local Variables
C
      INTEGER               AVAL
      INTEGER               PTRDEX
      INTEGER               NFOUND
      INTEGER               ROOM

C
C     Standard SPICE error handling.
C

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZUNTNGL' )

      IF ( NPTR .GT. MAXCEL ) THEN

         CALL SETMSG ( 'Input pointer array is larger than '
     .   //            'cell array. Pointer array size = '
     .   //            '#1. Cell array size = #2.'          )
         CALL ERRINT ( '#1', NPTR                           )
         CALL ERRINT ( '#2', MAXCEL                         )
         CALL SIGERR ( 'SPICE(BARRAYTOOSMALL)'              )
         CALL CHKOUT ( 'ZZUNTNGL'                           )
         RETURN

      END IF

C
C     ROOM is the remaining room in the output list.
C     
      ROOM = MAXB

C
C     Initialize pointer index.
C
      PTRDEX = 0

C
C     Loop over all A-values.
C
      DO AVAL = 1, NPTR
C
C        Traverse the chained list for a particular A-value 
C        and collect associated B-values. If B-values exists,
C        return the number of B-vals to the OUTLST array
C        at element PTRDEX + 1; return the list of B-vals
C        starting at element PTRDEX + 2.
C
C        Make sure the output pointers below are in range.
C
         IF ( PTRDEX+2 .GT. MAXB ) THEN

            CALL SETMSG ( 'Index larger than output array. '
     .      //            'Index = #1. Array size = #2.'     )
            CALL ERRINT ( '#1', PTRDEX+2                     )
            CALL ERRINT ( '#2', MAXB                         )
            CALL SIGERR ( 'SPICE(BARRAYTOOSMALL)'            )
            CALL CHKOUT ( 'ZZUNTNGL'                         )
            RETURN

         END IF

         IF ( ROOM .LE. 0 ) THEN

            CALL SETMSG ( 'Remaining room in output array is #1. '
     .      //            'Current input pointer index = #2. Output '
     .      //            'array size = #3. Output pointer index '
     .      //            'is #4.'                                  )
            CALL ERRINT ( '#1', ROOM                                )
            CALL ERRINT ( '#2', AVAL                                )
            CALL ERRINT ( '#3', MAXB                                )
            CALL ERRINT ( '#4', PTRDEX                              )
            CALL SIGERR ( 'SPICE(BARRAYTOOSMALL)'                   )
            CALL CHKOUT ( 'ZZUNTNGL'                                )
            RETURN

         END IF

         CALL ZZTRVLNK ( AVAL,    NPTR,  PNTRS, 
     .                   MAXCEL,  CELLS, ROOM,
     .                   OUTLST(PTRDEX + 1), 
     .                   OUTLST(PTRDEX + 2)     )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZUNTNGL' )
            RETURN
         END IF

C
C        Increment pointer PTRDEX if we found any B-vals.
C
         NFOUND = OUTLST(PTRDEX + 1)

         IF ( NFOUND .GT. 0 ) THEN
C
C           Store in PNTRS the pointer to the returned list.
C           This assignment overwrites the input pointer from
C           AVAL to the head of the associated list in CELLS.
C
            PNTRS(AVAL) = PTRDEX + 1
C
C           Update the count of available spaces in the 
C           output list. Account for the list size stored
C           at the front of the list.
C
            ROOM = ROOM - ( 1 + NFOUND )            
C
C           Increment PTRDEX to mark the position of the final
C           B-val.
C
            PTRDEX = PTRDEX + 1 + NFOUND

         ELSE
C
C           If no associated B-values exist, set the
C           PNTRS element to -1, indicating no B-value.
C
            PNTRS(AVAL) = -1

         END IF

      END DO

C
C     Return the current pointer value.
C
      NOUT = PTRDEX

C
C     Standard SPICE error handling.
C
      CALL CHKOUT ( 'ZZUNTNGL' )
      RETURN
      END
