C$Procedure      ZZADDLNK ( Add a new AB cell to an AB structure )
 
      SUBROUTINE ZZADDLNK ( AVAL,   BVAL,  MAXA, 
     .                      CELLSZ, PNTRS, NCELL, CELLS )
      IMPLICIT NONE

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Add a new AB cell to an AB cell pool and update the A-value
C     pointer array.
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
 
      INTEGER               AVAL
      INTEGER               BVAL
      INTEGER               MAXA
      INTEGER               CELLSZ
      INTEGER               PNTRS  ( * )
      INTEGER               NCELL
      INTEGER               CELLS  ( 2, * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     AVAL       I   A-value index associated with new pool entry.
C     BVAL       I   B-value associated with input A value.
C     MAXA       I   Size of pointer array.
C     CELLSZ     I   Number of cells in cell array.
C     PNTRS     I/O  A-value pointer array to update.
C     NCELL     I/O  Number of cells in AB cell linked-list.
C     CELLS     I/O  Cell array.
C
C$ Detailed_Input
C
C     AVAL       is an index associated with a new pool entry. AVAL
C                is the index of the Ith value of a finite set "A."
C                AVAL is an index into the PNTRS array. 
C
C     BVAL       is a B-value associated with AVAL. A new cell is to be
C                allocated to associate AVAL and BVAL.
C
C     MAXA       Size of pointer array. AVAL must not exceed
C                this value.
C
C     CELLSZ     Number of cells in the cell array CELLS. The array
C                has dimension (2, CELLSZ).
C
C     PNTRS      is an array of integers acting as pointers into the
C                cell array. The element at index I of PNTRS points to
C                the head of a linked list in CELLS that associates
C                B-values with the index I. Here I is the "A value."
C
C     NCELL      Number of cells in use on input. NCELL is incremented
C                by one on output.
C
C     CELLS      Cell array, updated on output.
C
C$ Detailed_Output
C
C     PNTRS      is the input pointer array, modified to reflect the 
C                addition of a new cell entry. The element
C
C                   PNTRS(AVAL)
C
C                contains the second subscript in CELLS of the new
C                entry.
C
C     NCELL      Number of cells in use on output. NCELL is incremented
C                by one relative to its input value.
C
C
C     CELLS      is the input cell array, with a new entry. The element
C
C                   CELLS ( 1, PNTRS(AVAL) )
C
C                contains BVAL. The element
C
C                   CELLS ( 2, PNTRS(AVAL) )
C
C                is a pointer to the previous head of the list
C                associated with AVAL, or is null (-1) if the list was
C                previously empty.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C 
C     1)  If the array index AVAL exceeds the size MAXA of the
C         A-value pointer array, the error SPICE(AVALOUTOFRANGE) is
C         signaled.
C
C         Occurrence of this error in the context of DSK type 2 spatial
C         index creation can result from the fine voxel scale having
C         been set to a value too small for the plate set.
C
C     2)  If the required cell count exceeds the size CELLSZ of the 
C         cell array, the error SPICE(CELLARRAYTOOSMALL) is signaled.
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
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    MKDSK Version 2.1.0, 15-JAN-2016 (NJB)
C
C        Updated error check on AVAL.
C
C        Filled in or update various header entries.
C
C-    MKDSK Version 2.0.0, 30-APR-2014 (NJB)
C     
C        Argument list change: now accepts separate inputs for pointer
C        list size and cell size.
C
C-    MKDSK Version 1.2.0, 04-MAY-2010 (NJB)
C
C        Changed INCLUDE file to dsk02.inc.
C
C-    MKDSK Version 1.1.1, 08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C-    MKDSK Version 1.1.0, 11-JUN-2004 (EDW)
C
C        Added check on AVAL to ensure it's smaller
C        than the size of CELLSZ.
C
C-    MKDSK Version 1.0.0, 03-FEB-1999 (JAB)
C
C
C-&
 
C$ Index_Entries
C
C     add entry to AB cell pool
C
C-&

C
C     SPICE functions
C
      LOGICAL               RETURN


C
C     Standard RETURN test.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZADDLNK' )

C
C     Test the pointer array index AVAL.
C
      IF ( ( AVAL .LT. 1 ) .OR. ( AVAL .GT. MAXA ) ) THEN

         CALL SETMSG ( 'Index AVAL is out of range. AVAL = #1; '
     .   //            'valid range is 1:#2.'                   )
         CALL ERRINT ( '#1', AVAL                               )
         CALL ERRINT ( '#2', MAXA                               )
         CALL SIGERR ( 'SPICE(AVALOUTOFRANGE)'                  )
         CALL CHKOUT ( 'ZZADDLNK'                               )
         RETURN

      END IF

C
C     Increment the cell counter.
C
      NCELL = NCELL + 1

      IF ( NCELL .GT. CELLSZ ) THEN

         CALL SETMSG ( 'NCELL larger than cell array. '
     .   //            'Cell index = #1. Array size = #2.' )
         CALL ERRINT ( '#1', NCELL                         )
         CALL ERRINT ( '#2', CELLSZ                        )
         CALL SIGERR ( 'SPICE(CELLARRAYTOOSMALL)'          )
         CALL CHKOUT ( 'ZZADDLNK'                          )
         RETURN

      END IF
 
C
C     Update the cell address of the last occurrence of the A-value,
C     if any. If none, PNTRS(AVAL) has value -1.
C
      CELLS(1,NCELL) = BVAL
      CELLS(2,NCELL) = PNTRS(AVAL)

      PNTRS(AVAL)    = NCELL

      CALL CHKOUT ( 'ZZADDLNK' )
      RETURN
      END
