C$Procedure DSKD02 ( DSK, fetch d.p. type 2 data )
 
      SUBROUTINE DSKD02 ( HANDLE, DLADSC, ITEM, START, ROOM, N, VALUES )
 
C$ Abstract
C
C     Fetch double precision data from a type 2 DSK segment.
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
      DOUBLE PRECISION      VALUES ( * )
      
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
C                       KWDSC      DSK descriptor of segment.  See the
C                                  INCLUDE file dskdsc.inc for a
C                                  discussion of the contents of DSK
C                                  descriptors. Note that DSK
C                                  descriptors are not to be confused
C                                  with DLA descriptors, which contain
C                                  segment component base address and
C                                  size information.
C
C                       KWVTBD     Vertex bounds. This is an array of
C                                  six values giving the minimum and
C                                  maximum values of each component of
C                                  the vertex set.
C
C                       KWVXOR     Voxel grid origin. This is the
C                                  location of the voxel grid origin in
C                                  the body-fixed frame associated with
C                                  the target body.
C
C                       KWVXSZ     Voxel size.  DSK voxels are cubes;
C                                  the edge length of each cube is
C                                  given by the voxel size.  This
C                                  size applies to the fine voxel grid.
C                                  Units are km.
C
C
C     START          is the start index within specified data item from
C                    which data are to be fetched. The index of the 
C                    first element of each data item is 1.
C
C     ROOM           is the amount of room in the output array.  It
C                    is permissible to provide an output array 
C                    that has too little room to fetch an item in
C                    one call.
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
C        in the call  tree of this routine, but there are no
C        guarantees.
C
C     4) If ROOM is non-positive, the error SPICE(VALUEOUTOFRANGE)
C        is signaled.
C
C     5) If the coarse voxel scale read from the designated segment
C        is less than 1, the error PICE(VALUEOUTOFRANGE) is signaled.
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
C                    WRITE (*,*) ' '
C                    WRITE (*,*) 'Plate number: ', I
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
C          ...
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
C-    SPICELIB Version 1.0.0, 04-FEB-2017 (NJB)
C
C        Fixed typo in version description.
C
C     23-AUG-2016
C
C        Now saves NV and updates it only when the current
C        segment changes.
C
C     15-JAN-2016 (NJB)
C
C        Removed code involving parameter NP.
C        Updated header Examples section.
C
C-    DSKLIB Version 3.0.0, 13-MAY-2010 (NJB)
C
C        Updated for compatibility with new DSK type 2
C        segment design.
C
C-    DSKLIB Version 2.1.0, 20-APR-2010 (NJB)
C
C        Bug fix: changed declaration of output argument
C        VALUES to double precision.
C
C-    DSKLIB Version 2.0.0, 27-DEC-2006 (NJB)
C
C        Updated to remove support for min, max radius
C        lookup.  These values are now stored in DSK
C        descriptors.
C
C-    DSKLIB Version 1.0.0, 27-OCT-2006 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch d.p. data from a type 2 dsk segment
C
C-&
 

C
C     SPICELIB functions
C     
      LOGICAL               FAILED

C
C     Local parameters
C

C
C     Local variables
C
      INTEGER               B
      INTEGER               DBASE
      INTEGER               E
      INTEGER               IBASE
      INTEGER               NV
      INTEGER               PRVHAN
      INTEGER               PRVBAS
      INTEGER               SIZE
      
      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  PRVBAS
      SAVE                  PRVHAN
      SAVE                  NV

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
C        Set previous segment base integer address to an invalid value
C        as well.
C
         PRVBAS = -1 

         FIRST  = .FALSE.

      END IF


      IF ( ROOM .LE. 0 ) THEN
         
         CALL CHKIN  ( 'DSKD02' )
         CALL SETMSG ( 'ROOM was #; must be positive.' )
         CALL ERRINT ( '#',  ROOM                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'        )
         CALL CHKOUT ( 'DSKD02'                        )
         RETURN

      END IF


      IBASE  = DLADSC ( IBSIDX )
      DBASE  = DLADSC ( DBSIDX )

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
C        Read NV.
C     
         CALL DASRDI ( HANDLE, IBASE+IXNV, IBASE+IXNV, NV )

         IF ( FAILED() ) THEN
            RETURN
         END IF
C
C        Update the saved handle value.
C
         PRVHAN = HANDLE

C
C        Update the saved integer base address.
C
         PRVBAS = IBASE

      END IF

C
C     Branch based on the item to be returned.
C
C     Note that we haven't checked the validity of START; we'll do this
C     after the IF block.
C     
      IF ( ITEM .EQ. KWVERT ) THEN
C
C        Return vertices.  There are 3*NV values in all.
C        First locate the data.
C        
         SIZE =   3*NV
         B    = ( DBASE + IXVERT ) + START - 1 


      ELSE IF ( ITEM .EQ. KWDSC ) THEN
C
C        Return DSK descriptor.
C
         SIZE =   DSKDSZ
         B    = ( DBASE + IXDSCR ) + START - 1 


      ELSE IF ( ITEM .EQ. KWVTBD ) THEN
C
C        Return vertex bounds.  There are 6 elements.
C
         SIZE =   6
         B    = ( DBASE + IXVTBD ) + START - 1 


      ELSE IF ( ITEM .EQ. KWVXOR ) THEN
C
C        Return voxel grid origin.  There are 3 elements.
C
         SIZE =   3
         B    = ( DBASE + IXVXOR ) + START - 1 


      ELSE IF ( ITEM .EQ. KWVXSZ ) THEN
C
C        Return voxel size.  This is a scalar.
C
         SIZE =   1
         B    = ( DBASE + IXVXSZ ) + START - 1 

      ELSE

         CALL CHKIN  ( 'DSKD02' )
         CALL SETMSG ( 'Keyword parameter # was not recognized.' )
         CALL ERRINT ( '#',  ITEM                                )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                     )
         CALL CHKOUT ( 'DSKD02'                                  )
         RETURN

      END IF

C
C     The valid range for START is 1:SIZE.
C
      IF (  ( START .LT. 1 ) .OR. ( START .GT. SIZE )  ) THEN

         CALL CHKIN  ( 'DSKD02' )
         CALL SETMSG ( 'START must be in the range defined ' //
     .                 'by the size of the data associated ' //
     .                 'with the keyword parameter #, '      //
     .                 'namely 1:#.  Actual value of START ' //
     .                 'was #.'                              )
         CALL ERRINT ( '#',  ITEM                            )
         CALL ERRINT ( '#',  SIZE                            )
         CALL ERRINT ( '#',  START                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKD02'                              )
         RETURN

      END IF

C
C     Read the requested data.  We already have the start address B.
C
      N  =   MIN ( ROOM, SIZE - START + 1 )
      E  =   B + N - 1 


      CALL DASRDD ( HANDLE, B, E, VALUES )

 
      RETURN
      END 

