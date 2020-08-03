C$Procedure DSKI04 ( DSK, fetch integer type 4 data )
 
      SUBROUTINE DSKI04 ( HANDLE, DLADSC, ITEM, START, ROOM, N, VALUES )
 
C$ Abstract
C
C     Fetch integer data from a type 4 DSK segment.
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
      INCLUDE 'dsk04.inc'

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
C     HANDLE         is the handle of a DSK file containing a type 4
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
C                    See the INCLUDE file dsk04.inc for values
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
C         dsk04.inc
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
C     5) If the input keyword parameter is not recognized, the error
C        SPICE(NOTSUPPORTED) is signaled.
C
C     6) If START is less than 1 or greater than the size of the
C        item to be fetched, the error SPICE(INDEXOUTOFRANGE) is
C        signaled.
C
C$ Files
C
C     See input argument HANDLE.
C
C$ Particulars
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
C     None.
C
C$ Restrictions
C
C     1)  This is a prototype routine.  The interface is not expected
C         to change, but there are no guarantees.
C
C     2)  This routine uses discovery check-in to boost execution
C         speed.  However, this routine is in violation of NAIF
C         standards for use of discovery check-in:  routines called
C         from this routine may signal errors.  If errors are signaled
C         in called routines, this routine's name will be missing from
C         the traceback message.
C
C     3) This routine does not initialize the nested grid addressing
C        routines.
C
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
C-    DSKBRIEF Version 1.1.0, 06-OCT-2016 (NJB)
C
C        Removed call to ZZDSK4GI. This routine no longer
C        intializes the nested grid addressing routines.
C
C        Removed unused variables.
C
C-    DSKBRIEF Version 1.0.0, 04-OCT-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch integer data from a type_4_dsk segment
C
C-&
 


C
C     SPICELIB functions
C     
      LOGICAL               DLASSG
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
C
C     DBFSIZ is the size of a d.p. buffer used to
C     read parameters from the segment.
C
      INTEGER               DBFSIZ
      PARAMETER           ( DBFSIZ = 19 + DSKDSZ )

C
C     Local variables
C
      DOUBLE PRECISION      DBUFF  ( DBFSIZ )

      INTEGER               B
      INTEGER               DBASE
      INTEGER               E
      INTEGER               IBASE
      INTEGER               NC
      INTEGER               NDAT
      INTEGER               NDIMS
      INTEGER               NR
      INTEGER               PIXPTR
      INTEGER               PRVDSC ( DLADSZ )
      INTEGER               PRVHAN
      INTEGER               SIZE

      LOGICAL               PASS1

C
C     Saved variables
C
      SAVE                  DBUFF
      SAVE                  NC
      SAVE                  NR
      SAVE                  PASS1
      SAVE                  PIXPTR
      SAVE                  PRVDSC
      SAVE                  PRVHAN

C
C     Initial values
C
      DATA                  PASS1  / .TRUE.     /
      DATA                  PRVHAN / 0          /
      DATA                  PRVDSC / DLADSZ * 0 /
      

      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Use discovery check-in.  This is done for efficiency; note
C     however that this routine does not meet SPICE standards for
C     discovery check-in eligibility.
C
      IF ( ROOM .LE. 0 ) THEN
         
         CALL CHKIN  ( 'DSKI04' )
         CALL SETMSG ( 'ROOM was #; must be positive.' )
         CALL ERRINT ( '#',  ROOM                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'        )
         CALL CHKOUT ( 'DSKI04'                        )
         RETURN

      END IF

      IBASE  = DLADSC ( IBSIDX )
      DBASE  = DLADSC ( DBSIDX )

C
C     Either a new file or new segment in the same file
C     will require looking up the segment parameters.
C     To determine whether the segment is new, we don't
C     need to compare the entire DLA descriptor:  just
C     comparing the three base addresses of the descriptor
C     against the saved base addresses is sufficient.
C      
      IF (  PASS1  .OR.  .NOT.  
     .      DLASSG( HANDLE, PRVHAN, DLADSC, PRVDSC )  ) THEN
C
C        Treat the input file and segment as new.
C
C        Read the d.p. parameters first.  These are located at the
C        beginning of the d.p. component of the segment.
C     
         CALL DASRDD ( HANDLE, DBASE+1, DBASE+DBFSIZ, DBUFF )

C
C        Update the pixel pointer.
C
         PIXPTR = IBASE  +  NINT( DBUFF(IXPIXP) )

C
C        Update the grid dimensions.
C
         NC = NINT( DBUFF(IXNC) )
         NR = NINT( DBUFF(IXNR) )

C
C        This call may be reinstated for N0067. It's currently
C        unnecessary.
C        
C        CALL ZZDSK4GI ( HANDLE, DLADSC )
C

         IF ( .NOT. FAILED() ) THEN

            PASS1  = .FALSE.

C
C           Update the saved handle value.
C
            PRVHAN = HANDLE

C
C           Update the saved DLA descriptor.
C
            CALL MOVEI ( DLADSC, DLADSZ, PRVDSC )

         END IF

      END IF

C
C     Branch based on the item to be returned.
C
C     Note that we haven't checked the validity of START; we'll do this
C     after the IF block.
C     
      IF ( ITEM .EQ. KWRAW ) THEN
C
C        Return the specified raw data. 
C
C        The raw grid has NR rows and NC/2 columns.
C        There are two 16-bit pixels per stored integer.
C        The data are stored in row-major order.
C
C        The data are returned in packed form: two adjacent
C        16-bit values are returned in each integer.
C         
         NDAT = ( NC / 2 ) * NR 

C
C        START must be in the range 1:NDAT.
C
         IF (  ( START .LT. 1 ) .OR. ( START .GT. NDAT ) ) THEN

            CALL CHKIN  ( 'DSKI04'                                    )
            CALL SETMSG ( 'START must be in the range 1:# but was #.' )
            CALL ERRINT ( '#',  NDAT                                  )
            CALL ERRINT ( '#',  START                                 )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                    )
            CALL CHKOUT ( 'DSKI04'                                    )
            RETURN

         END IF

C
C        Let B be the base address of the set of stored
C        integers we'll read.
C
         B = PIXPTR + START - 2

C
C        Read data into the output array.
C
         N = MIN ( ROOM,  NDAT - START + 1 )

         CALL DASRDI ( HANDLE, B+1, B+N, VALUES )
            
C        Exit here, since we're not going to use the generic 
C        data transfer code at the end of this routine.
C
C        There's no CHKOUT call here since we're using 
C        discovery check-in.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWNDIM ) THEN
C
C        The item is the number of nested grid dimensions.
C
         SIZE = 1
         B    = IBASE + IXNDIM
         E    = B

      ELSE IF ( ITEM .EQ. KWGDIM ) THEN
C
C        The item is the array of grid dimensions.
C
         B = IBASE + IXNDIM

         CALL DASRDI ( HANDLE, B, B, NDIMS )

         B    =  IBASE + IXGDIM
         SIZE =  2 * NDIMS

      ELSE

         CALL CHKIN  ( 'DSKI04'                                  )
         CALL SETMSG ( 'Keyword parameter # was not recognized.' )
         CALL ERRINT ( '#',  ITEM                                )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                     )
         CALL CHKOUT ( 'DSKI04'                                  )
         RETURN

      END IF

C
C     The valid range for START is 1:SIZE.
C
      IF (  ( START .LT. 1 ) .OR. ( START .GT. SIZE )  ) THEN

         CALL CHKIN  ( 'DSKI04'                              )
         CALL SETMSG ( 'START must be in the range defined ' //
     .                 'by the size of the data associated ' //
     .                 'with the keyword parameter #, '      //
     .                 'namely 1:#.  Actual value of START ' //
     .                 'was #.'                              )
         CALL ERRINT ( '#',  ITEM                            )
         CALL ERRINT ( '#',  SIZE                            )
         CALL ERRINT ( '#',  START                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKI04'                              )
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

