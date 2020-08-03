C$Procedure DSKI02 ( DSK, fetch integer type 2 data )
 
      SUBROUTINE DSKI02 ( HANDLE, DLADSC, ITEM, START, ROOM, N, VALUES )
 
C$ Abstract
C
C     Fetch integer data from a type 2 DSK segment.
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
C     DAS
C     DSK
C
C$ Keywords
C
C     DAS
C     DSK
C     FILES
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk02.inc'

      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      INTEGER               ITEM
      INTEGER               START
      INTEGER               ROOM
      INTEGER               N
      INTEGER               VALUES ( * )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.
C     DLADSC     I   DLA descriptor.
C     ITEM       I   Keyword identifying item to fetch.
C     START      I   Start index.
C     ROOM       I   Amount of room in output array.
C     N          O   Number of values returned.
C     VALUES     O   Array containing requested item.
C
C$ Detailed_Input
C
C     HANDLE         is the handle of a DSK file containing a type 2
C                    segment from which data are to be fetched.
C
C     DLADSC         is the DLA descriptor associated with the segment
C                    from which data are to be fetched. 
C
C     ITEM           is an integer "keyword" parameter designating the
C                    item to fetch. In the descriptions below, note
C                    that "model" refers to the model represented by
C                    the designated segment.  This model may be a
C                    subset of a larger model.
C
C                    Names and meanings of parameters supported by this 
C                    routine are:
C
C                       KWNV       Number of vertices in model.
C
C                       KWNP       Number of plates in model.
C
C                       KWNVXT     Total number of voxels in fine grid.
C
C                       KWVGRX     Voxel grid extent.  This extent is
C                                  an array of three integers
C                                  indicating the number of voxels in
C                                  the X, Y, and Z directions in the
C                                  fine voxel grid.
C
C                       KWCGSC     Coarse voxel grid scale.  The extent
C                                  of the fine voxel grid is related to
C                                  the extent of the coarse voxel grid
C                                  by this scale factor.
C
C                       KWVXPS     Size of the voxel-to-plate pointer
C                                  list. 
C
C                       KWVXLS     Voxel-plate correspondence list size.
C
C                       KWVTLS     Vertex-plate correspondence list
C                                  size.
C
C                       KWPLAT     Plate array.  For each plate, this
C                                  array contains the indices of the
C                                  plate's three vertices.  The ordering
C                                  of the array members is:
C
C                                     Plate 1 vertex index 1
C                                     Plate 1 vertex index 2
C                                     Plate 1 vertex index 3
C                                     Plate 2 vertex index 1
C                                             ...
C
C                       KWVXPT     Voxel-plate pointer list. This list
C                                  contains pointers that map fine
C                                  voxels to lists of plates that
C                                  intersect those voxels. Note that
C                                  only fine voxels belonging to
C                                  non-empty coarse voxels are in the
C                                  domain of this mapping.
C
C                       KWVXPL     Voxel-plate correspondence list.
C                                  This list contains lists of plates
C                                  that intersect fine voxels. (This
C                                  list is the data structure into
C                                  which the voxel-to-plate pointers
C                                  point.)  This list can contain
C                                  empty lists.
C                                                                   
C                       KWVTPT     Vertex-plate pointer list. This list
C                                  contains pointers that map vertices
C                                  to lists of plates to which those
C                                  vertices belong.
C
C                                  Note that the size of this list is
C                                  always NV, the number of vertices.
C                                  Hence there's no need for a separate
C                                  keyword for the size of this list.
C
C                       KWVTPL     Vertex-plate correspondence list.
C                                  This list contains, for each vertex,
C                                  the indices of the plates to which
C                                  that vertex belongs.
C
C                       KWCGPT     Coarse voxel grid pointers.  This is
C                                  an array of pointers mapping coarse
C                                  voxels to lists of pointers in the
C                                  voxel-plate pointer list.  Each 
C                                  non-empty coarse voxel maps to a
C                                  list of pointers; every fine voxel
C                                  contained in a non-empty coarse voxel
C                                  has its own pointers. Grid elements
C                                  corresponding to empty coarse voxels
C                                  have null (non-positive) pointers.
C
C                    See the INCLUDE file dsk.inc for values
C                    associated with the keyword parameters.
C
C
C     START          is the start index within the specified data item
C                    from which data are to be fetched. The index of
C                    the first element of each data item is 1. START
C                    has units of integers; for example, the start
C                    index of the second plate is 4, since each plate
C                    occupies three integers.
C
C     ROOM           is the amount of room in the output array. It is
C                    permissible to provide an output array that has
C                    too little room to fetch an item in one call. ROOM
C                    has units of integers: for example, the room
C                    required to fetch one plate is 3.
C
C$ Detailed_Output
C
C     N              is the number of elements fetched to the output
C                    array VALUES.  N is normally in the range 
C                    1:ROOM; if an error occurs on the call, N is
C                    undefined.
C
C     VALUES         is a contiguous set of elements of the item
C                    designated by ITEM.  The correspondence of 
C                    VALUES at the elements of the data item is:
C
C                       VALUES(1)      ITEM(START)
C                         ...             ...
C                       VALUES(N)      ITEM(START+N-1)
C                    
C                    If an error occurs on the call, VALUES is 
C                    undefined.
C
C$ Parameters
C
C     See the INCLUDE files 
C
C         dla.inc
C         dsk02.inc
C         dskdsc.inc
C
C$ Exceptions
C
C     1) If the input handle is invalid, the error will be diagnosed by
C        routines in the call tree of this routine. 
C
C     2) If a file read error occurs, the error will be diagnosed by
C        routines in the call tree of this routine.
C
C     3) If the input DLA descriptor is invalid, the effect of this
C        routine is undefined. The error *may* be diagnosed by routines
C        in the call tree of this routine, but there are no
C        guarantees.
C
C     4) If ROOM is non-positive, the error SPICE(VALUEOUTOFRANGE)
C        is signaled.
C
C     5) If the coarse voxel scale read from the designated segment
C        is less than 1, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     6) If the input keyword parameter is not recognized, the error
C        SPICE(NOTSUPPORTED) is signaled.
C
C     7) If START is less than 1 or greater than the size of the
C        item to be fetched, the error SPICE(INDEXOUTOFRANGE) is
C        signaled.
C
C$ Files
C
C     See input argument HANDLE.
C
C$ Particulars
C
C     Most SPICE applications will not need to call this routine. The
C     routines DSKV02, DSKP02, and DSKZ02 provide a higher-level
C     interface for fetching DSK type 2 vertex and plate data.
C
C     DSK files are built using the DLA low-level format and
C     the DAS architecture; DLA files are a specialized type of DAS
C     file in which data are organized as a doubly linked list of
C     segments.  Each segment's data belong to contiguous components of
C     character, double precision, and integer type.
C
C     Note that the DSK descriptor for the segment is not needed by
C     this routine; the DLA descriptor contains the base address and
C     size information for the integer, double precision, and character
C     components of the segment, and these suffice for the purpose of
C     fetching data.
C
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C 
C     1) Look up all the vertices associated with each plate
C        of the model contained in a specified type 2 segment.
C        For this example, we'll show the context of this look-up:
C        opening the DSK file for read access, traversing a trivial,
C        one-segment list to obtain the segment of interest.
C
C
C        Example code begins here.
C
C
C                 PROGRAM EX1 
C                 IMPLICIT NONE
C
C                 INCLUDE 'dla.inc'
C                 INCLUDE 'dskdsc.inc'
C                 INCLUDE 'dsk02.inc'
C
C           C
C           C     Local parameters
C           C 
C                 CHARACTER*(*)         FMT
C                 PARAMETER           ( FMT    = '(1X,A,3(1XE16.9))' )
C
C                 INTEGER               FILSIZ
C                 PARAMETER           ( FILSIZ = 255 )
C
C           C
C           C     Local variables
C           C 
C                 CHARACTER*(FILSIZ)    DSK
C
C                 DOUBLE PRECISION      VRTCES ( 3, 3 )
C
C                 INTEGER               DLADSC ( DLADSZ )
C                 INTEGER               HANDLE
C                 INTEGER               I
C                 INTEGER               J
C                 INTEGER               N
C                 INTEGER               NP
C                 INTEGER               START
C                 INTEGER               VRTIDS ( 3 )
C
C                 LOGICAL               FOUND
C
C
C           C
C           C     Prompt for the name of the DSK to read.
C           C
C                 CALL PROMPT ( 'Enter DSK name > ', DSK )
C           C
C           C     Open the DSK file for read access.
C           C     We use the DAS-level interface for
C           C     this function.
C           C
C                 CALL DASOPR ( DSK, HANDLE )
C
C           C
C           C     Begin a forward search through the
C           C     kernel, treating the file as a DLA.
C           C     In this example, it's a very short
C           C     search.
C           C
C                 CALL DLABFS ( HANDLE, DLADSC, FOUND )
C
C                 IF ( .NOT. FOUND ) THEN
C           C
C           C        We arrive here only if the kernel
C           C        contains no segments.  This is 
C           C        unexpected, but we're prepared for it.
C           C
C                    CALL SETMSG ( 'No segments found '
C                .   //            'in DSK file #.'    )
C                    CALL ERRCH  ( '#',  DSK           )
C                    CALL SIGERR ( 'SPICE(NODATA)'     )
C
C                 END IF
C
C           C
C           C     If we made it this far, DLADSC is the
C           C     DLA descriptor of the first segment.
C           C
C           C     Find the number of plates in the model.
C           C
C                 CALL DSKI02 ( HANDLE, DLADSC, KWNP, 1, 1, N, NP )
C
C           C
C           C     For each plate, look up the desired data.
C           C
C                 DO I = 1, NP
C           C
C           C        For the Ith plate, find the associated 
C           C        vertex IDs.  We must take into account
C           C        the fact that each plate has three
C           C        vertices when we compute the start 
C           C        index.
C           C
C                    START = 3*(I-1)+1
C
C                    CALL DSKI02 ( HANDLE, DLADSC, KWPLAT, START, 
C                .                 3,      N,      VRTIDS        )
C
C                    DO J = 1, 3
C           C
C           C            Fetch the vertex associated with
C           C            the Jth vertex ID.  Again, each
C           C            vertex is a 3-vector.  Note that
C           C            the vertices are double-precision
C           C            data, so we fetch them using 
C           C            DSKD02.
C           C
C                        START = 3*( VRTIDS(J) - 1 ) + 1
C
C                        CALL DSKD02 ( HANDLE, DLADSC, KWVERT,  START,
C                .                     3,      N,      VRTCES(1,J)    )
C                    END DO
C
C           C
C           C        Display the vertices of the Ith plate:
C           C   
C                    WRITE (*,*)   ' '
C                    WRITE (*,*)   'Plate number: ', I
C                    WRITE (*,FMT) '   Vertex 1: ', (VRTCES(J,1), J=1,3)
C                    WRITE (*,FMT) '   Vertex 2: ', (VRTCES(J,2), J=1,3)
C                    WRITE (*,FMT) '   Vertex 3: ', (VRTCES(J,3), J=1,3)
C         
C                 END DO
C
C           C
C           C     Close the kernel.  This isn't necessary in a stand-
C           C     alone program, but it's good practice in subroutines
C           C     because it frees program and system resources.
C           C
C                 CALL DASCLS ( HANDLE )
C
C                 END
C
C
C     When this program was executed on a PC/Linux/gfortran/64bit
C     platform, using a DSK file representing a regular icosahedron,
C     the output was:
C
C
C      Enter DSK name > solid.bds
C
C       Plate number:            1
C          Vertex 1:   0.000000000E+00  0.000000000E+00  0.117557000E+01
C          Vertex 2:   0.105146000E+01  0.000000000E+00  0.525731000E+00
C          Vertex 3:   0.324920000E+00  0.100000000E+01  0.525731000E+00
C
C       Plate number:            2
C          Vertex 1:   0.000000000E+00  0.000000000E+00  0.117557000E+01
C          Vertex 2:   0.324920000E+00  0.100000000E+01  0.525731000E+00
C          Vertex 3:  -0.850651000E+00  0.618034000E+00  0.525731000E+00
C
C         ...
C
C       Plate number:           20
C          Vertex 1:   0.850651000E+00 -0.618034000E+00 -0.525731000E+00
C          Vertex 2:   0.000000000E+00  0.000000000E+00 -0.117557000E+01
C          Vertex 3:   0.850651000E+00  0.618034000E+00 -0.525731000E+00
C
C
C$ Restrictions
C
C     1) This routine uses discovery check-in to boost
C        execution speed.  However, this routine is in
C        violation of NAIF standards for use of discovery
C        check-in:  routines called from this routine may
C        signal errors.  If errors are signaled in called
C        routines, this routine's name will be missing 
C        from the traceback message.
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
C-    SPICELIB Version 1.0.0, 22-NOV-2016 (NJB)
C
C        Added FAILED check after segment attribute fetch calls.
C        Re-ordered code so that values are saved only after
C        all error checks have passed. Simplified base address
C        comparisons.
C
C     15-JAN-2016 (NJB)
C
C        Updated header Examples and Particulars sections.
C
C     DSKLIB Version 1.0.2, 11-JUL-2014 (NJB)
C
C        Fixed a trivial header comment typo.
C
C     DSKLIB Version 1.0.1, 13-MAY-2010 (NJB)
C
C        Updated header.
C
C     DSKLIB Version 1.0.0, 27-OCT-2006 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch integer data from a type 2 dsk segment
C
C-&
 

C
C     SPICELIB functions
C     
      LOGICAL              FAILED

C
C     Local parameters
C
C
C     IBFSIZ is the size of an integer buffer used to
C     read parameters from the segment.
C
      INTEGER               IBFSIZ
      PARAMETER           ( IBFSIZ = 10 )

C
C     Local variables
C
      INTEGER               B
      INTEGER               CGSCAL
      INTEGER               E
      INTEGER               IBASE
      INTEGER               IBUFF  ( IBFSIZ )
      INTEGER               NCGR
      INTEGER               NP
      INTEGER               NV
      INTEGER               NVXTOT
      INTEGER               PRVBAS
      INTEGER               PRVHAN
      INTEGER               SIZE
      INTEGER               VOXNPL
      INTEGER               VOXNPT
      INTEGER               VTXNPL

      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  CGSCAL
      SAVE                  FIRST
      SAVE                  NP
      SAVE                  NV
      SAVE                  NVXTOT
      SAVE                  PRVBAS
      SAVE                  PRVHAN
      SAVE                  VOXNPL
      SAVE                  VOXNPT
      SAVE                  VTXNPL

C
C     Initial values
C
      DATA                  FIRST  / .TRUE. /


C
C     Use discovery check-in.  This is done for efficiency; note
C     however that this routine does not meet SPICE standards for
C     discovery check-in eligibility.
C
      IF ( FIRST ) THEN
C
C        Make sure we treat the input handle as new on the first pass.
C        Set PRVHAN to an invalid handle value.
C
         PRVHAN = 0

C
C        Set the previous segment base integer address to an invalid
C        value as well.
C
         PRVBAS = -1

         FIRST  = .FALSE.

      END IF


      IF ( ROOM .LE. 0 ) THEN
         
         CALL CHKIN  ( 'DSKI02' )
         CALL SETMSG ( 'ROOM was #; must be positive.' )
         CALL ERRINT ( '#',  ROOM                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'        )
         CALL CHKOUT ( 'DSKI02'                        )
         RETURN

      END IF


      IBASE  = DLADSC ( IBSIDX )
C
C     Either a new file or new segment in the same file will require
C     looking up the segment parameters. To determine whether the
C     segment is new, we don't need to compare the entire DLA
C     descriptor: just comparing the integer base address of the
C     descriptor against the saved integer base address is sufficient.
C     
C     DSK type 2 segments always have a non-empty integer component, so
C     each type 2 segment in a given file will have a distinct integer
C     base address. Segments of other types might not contain integers,
C     but they can't share an integer base address with a type 2
C     segment.
C      
      IF (      ( HANDLE .NE. PRVHAN )
     .     .OR. ( IBASE  .NE. PRVBAS )  ) THEN
C
C        Treat the input file and segment as new.
C
C        Read the integer parameters first.  These are located at the
C        beginning of the integer component of the segment.
C     
         CALL DASRDI ( HANDLE, IBASE+1, IBASE+IBFSIZ, IBUFF )

         IF ( FAILED() ) THEN
            RETURN
         END IF

C
C        Check the coarse voxel scale.
C
         IF ( IBUFF(IXCGSC) .LT. 1 ) THEN

            CALL CHKIN  ( 'DSKI02'                         )
            CALL SETMSG ( 'Coarse voxel grid scale is #; ' //
     .                    'this scale should be an '       //
     .                    'integer > 1'                    )
            CALL ERRINT ( '#',  CGSCAL                     )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
            CALL CHKOUT ( 'DSKI02'                         )
            RETURN

         END IF

C
C        All checks have passed. We can safely store the segment
C        parameters.
C
         NV     = IBUFF ( IXNV   ) 
         NP     = IBUFF ( IXNP   ) 
         NVXTOT = IBUFF ( IXNVXT ) 
         CGSCAL = IBUFF ( IXCGSC )
         VTXNPL = IBUFF ( IXVTLS )
         VOXNPT = IBUFF ( IXVXPS )
         VOXNPL = IBUFF ( IXVXLS )

C
C        Update the saved handle value.
C
         PRVHAN = HANDLE

C
C        Update the saved base integer address.
C
         PRVBAS = IBASE

      END IF

C
C     Branch based on the item to be returned.
C
C     Note that we haven't checked the validity of START; we'll do this
C     after the IF block.
C     
      IF ( ITEM .EQ. KWNV ) THEN
C
C        Return the number of vertices.
C
         N         = 1
         VALUES(1) = NV
         
C
C        As long as START is valid, we can return. Otherwise,
C        let control pass to the error handling block near
C        the end of this routine.
C
         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF

      ELSE IF ( ITEM .EQ. KWNP ) THEN
C
C        Return the number of plates.
C
         N         = 1
         VALUES(1) = NP

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWNVXT ) THEN
C
C        Return the total number of voxels.
C
         N         = 1
         VALUES(1) = NVXTOT

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWVGRX ) THEN
C
C        Return the voxel grid extents.
C
         SIZE = 3
         B    = IBASE + IXVGRX + START - 1


      ELSE IF ( ITEM .EQ. KWCGSC ) THEN
C
C        Return the coarse voxel grid scale.
C
         N         = 1
         VALUES(1) = CGSCAL

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWVXPS ) THEN
C
C        Return the voxel-plate pointer list size.
C
         N         = 1
         VALUES(1) = VOXNPT

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWVXLS ) THEN
C
C        Return the voxel-plate list size.
C
         N         = 1
         VALUES(1) = VOXNPL

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWVTLS ) THEN
C
C        Return the vertex-plate list size.
C
         N         = 1
         VALUES(1) = VTXNPL

         IF ( START .EQ. 1 ) THEN
            RETURN
         END IF


      ELSE IF ( ITEM .EQ. KWPLAT ) THEN
C
C        Return plate data.  There are 3*NP values in all.  First
C        locate the data.
C
         SIZE = 3*NP
         B    = ( IBASE + IXPLAT ) + START - 1 
         

      ELSE IF ( ITEM .EQ. KWVXPT )  THEN
C
C        Return voxel pointer data.  There are VOXNPT values in all.
C        First locate the data.
C
         SIZE = VOXNPT
         B    = ( IBASE + IXPLAT + 3*NP ) + START - 1 


      ELSE IF ( ITEM .EQ. KWVXPL )  THEN
C
C        Return voxel-plate list data.  There are VOXNPL values in all.
C        First locate the data.
C
         SIZE = VOXNPL
         B    = ( IBASE + IXPLAT + 3*NP + VOXNPT ) + START - 1 


      ELSE IF ( ITEM .EQ. KWVTPT )  THEN
C
C        Return vertex-plate pointer data.  There are NV values in all.
C        First locate the data.
C
         SIZE =    NV
         B    =  ( IBASE + IXPLAT + 3*NP + VOXNPT + VOXNPL )
     .           + START - 1 


      ELSE IF ( ITEM .EQ. KWVTPL )  THEN
C
C        Return vertex-plate list data.  There are VTXNPL values in
C        all. First locate the data.
C
         SIZE = VTXNPL

         B =  ( IBASE + IXPLAT + 3*NP + VOXNPT + VOXNPL + NV )
     .        + START - 1 


      ELSE IF ( ITEM .EQ. KWCGPT )  THEN
C
C        Compute the coarse grid size.
C
         NCGR = NVXTOT / (CGSCAL**3)

C
C        Return the coarse voxel grid occupancy pointers.  There are 
C
C           NCGR
C
C        values in all. First locate the data.
C
         SIZE = NCGR

         B =  (   IBASE  + IXPLAT + 3*NP   + VOXNPT 
     .          + VOXNPL + NV     + VTXNPL          )
     .          + START  - 1 


      ELSE

         CALL CHKIN  ( 'DSKI02'                                  )
         CALL SETMSG ( 'Keyword parameter # was not recognized.' )
         CALL ERRINT ( '#',  ITEM                                )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                     )
         CALL CHKOUT ( 'DSKI02'                                  )
         RETURN

      END IF

C
C     The valid range for START is 1:SIZE.
C
      IF (  ( START .LT. 1 ) .OR. ( START .GT. SIZE )  ) THEN

         CALL CHKIN  ( 'DSKI02'                              )
         CALL SETMSG ( 'START must be in the range defined ' //
     .                 'by the size of the data associated ' //
     .                 'with the keyword parameter #, '      //
     .                 'namely 1:#.  Actual value of START ' //
     .                 'was #.'                              )
         CALL ERRINT ( '#',  ITEM                            )
         CALL ERRINT ( '#',  SIZE                            )
         CALL ERRINT ( '#',  START                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKI02'                              )
         RETURN

      END IF

C
C     Read the requested data.  We already have the start address B.
C
      N  =  MIN ( ROOM,  SIZE - START + 1 )
      E  =  B + N - 1 

      CALL DASRDI ( HANDLE, B, E, VALUES )

      RETURN
      END 

