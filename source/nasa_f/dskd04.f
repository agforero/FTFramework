C$Procedure DSKD04 ( DSK, fetch d.p. type 4 data )
 
      SUBROUTINE DSKD04 ( HANDLE, DLADSC, ITEM, START, ROOM, N, VALUES )
 
C$ Abstract
C
C     Fetch double precision data from a type 4 DSK segment.
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
C                    Names and meanings of parameters supported by this
C                    routine are:
C
C                       KWDSC      DSK descriptor of segment. See the
C                                  INCLUDE file dskdsc.inc for a
C                                  discussion of the contents of DSK
C                                  descriptors. Note that DSK
C                                  descriptors are not to be confused
C                                  with DLA descriptors, which contain
C                                  segment component base address and
C                                  size information.
C
C
C     START          is the start index within specified data item from
C                    which data are to be fetched. The index of the 
C                    first element of each data item is 1.
C
C     ROOM           is the amount of room in the output array.  It
C                    is permissible to provided an output array 
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
C        in the call  tree of this routine, but there are no
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
C     character, double precision, and integer type. Depending on
C     the segment type, some of these components may be empty.
C
C     Note that the DSK descriptor for the segment is not needed by
C     this routine; the DLA descriptor contains the base address and
C     size information for the integer, double precision, and character
C     components of the segment, and these suffice for the purpose of
C     fetching data.
C
C$ Examples
C
C     1)
C
C
C
C$ Restrictions
C
C     1) This is a prototype routine.  The interface is not
C        expected to change, but there are no guarantees.
C
C     2) This routine uses discovery check-in to boost
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
C-    SPICELIB Version 2.0.0, 19-SEP-2012 (NJB)
C
C        Updated to handle new type 4 format. New items
C        are: number format, reference surface descriptor,
C        projection descriptor.
C
C-    SPICELIB Version 1.0.0, 17-AUG-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch d.p. data from a type_4_dsk_segment
C
C-&
 

C
C     SPICELIB functions
C     

C
C     Local parameters
C
      INTEGER               DBFSIZ
      PARAMETER           ( DBFSIZ = 100 )

C
C     Local variables
C
      DOUBLE PRECISION      DPBUFF ( DBFSIZ )

      INTEGER               ADDR
      INTEGER               B
      INTEGER               DPBASE
      INTEGER               DPSIZE
      INTEGER               DSIZE
      INTEGER               E
      INTEGER               MAXADD
      INTEGER               SIZE

C
C     Use discovery check-in.  This is done for efficiency; note
C     however that this routine does not meet SPICE standards for
C     discovery check-in eligibility.


C
C     Type 1 DSK segments contain only d.p. information,
C     so there's only one DLA component to consider.
C
      DPBASE = DLADSC ( DBSIDX )
      DPSIZE = DLADSC ( DSZIDX )

      MAXADD = DPBASE + DPSIZE


C      WRITE (*,*) 'DP BASE = ', DPBASE
C      WRITE (*,*) 'DP SIZE = ', DPSIZE
C      WRITE (*,*) 'ITEM    = ', ITEM
C      WRITE (*,*) ' '



      IF ( ROOM .LE. 0 ) THEN
         
         CALL CHKIN  ( 'DSKD04'                        )
         CALL SETMSG ( 'ROOM was #; must be positive.' )
         CALL ERRINT ( '#',  ROOM                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'        )
         CALL CHKOUT ( 'DSKD04'                        )
         RETURN

      END IF

C
C     Perform the best check on START that we can without
C     having determined the size of the item to fetch.
C
      IF (  ( START .LT. 1 ) .OR. ( START .GT. DPSIZE )  ) THEN

         CALL CHKIN  ( 'DSKD04' )
         CALL SETMSG ( 'START must be in the range defined ' //
     .                 'by the size of the data associated ' //
     .                 'with the keyword parameter #, '      //
     .                 'namely 1:#.  Actual value of START ' //
     .                 'was #.'                              )
         CALL ERRINT ( '#',  ITEM                            )
         CALL ERRINT ( '#',  DPSIZE                          )
         CALL ERRINT ( '#',  START                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKD04'                              )
         RETURN

      END IF

      
C 
C     Branch based on the item to be returned. In each case, set
C     the start address of the item and the size of the item.
C
      IF ( ITEM .EQ. KWDSC ) THEN
C
C        Return the DSK descriptor.
C        First locate the data.
C        
         B    =  DPBASE + 1
         SIZE =  DSKDSZ


      ELSE IF ( ITEM .EQ. KWFMTV ) THEN
C
C        Return the format version.
C
         B    = DPBASE + IXFMTV
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWNR ) THEN
C
C        Return the row count.
C
         B    = DPBASE + IXNR
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWNC ) THEN
C
C        Return the column count.
C
         B    = DPBASE + IXNC
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWCTR1 ) THEN
C
C        Return the first coordinate of the pixel grid center.
C
         B    = DPBASE + IXCTR1
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWCTR2 ) THEN
C
C        Return the second coordinate of the pixel grid center.
C
         B    = DPBASE + IXCTR2
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWPSZ1 ) THEN
C
C        Return the pixel size in the direction of coordinate 1.
C
         B    = DPBASE + IXPSZ1
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWPSZ2 ) THEN
C
C        Return the pixel size in the direction of coordinate 2.
C
         B    = DPBASE + IXPSZ2
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWNLOK ) THEN
C
C        Return the flag indicating whether nulls are allowed.
C        This is a d.p. parameter that maps to a Boolean value.
C
         B    = DPBASE + IXNLOK
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWNULL ) THEN
C
C        Return the null value parameter.
C
         B    = DPBASE + IXNULL
         SIZE = 1


      ELSE IF ( ITEM .EQ. KWHSCL ) THEN
C
C        Return the height scale parameter.
C
         B    = DPBASE + IXHSCL
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWNUMF ) THEN
C
C        Return the number format parameter.
C
         B    = DPBASE + IXNUMF
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWREFP ) THEN
C
C        Return the pointer to the reference surface
C        descriptor.
C
         B    = DPBASE + IXREFP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWPRJP ) THEN
C
C        Return the pointer to the map projection
C        descriptor.
C
         B    = DPBASE + IXPRJP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWITPP ) THEN
C
C        Return the pointer to the interpolation algorithm
C        descriptor.
C
         B    = DPBASE + IXITPP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWXP ) THEN
C
C        Return the pointer to the intercept algorithm
C        descriptor.
C
         B    = DPBASE + IXXP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWNVP ) THEN
C
C        Return the pointer to the normal vector algorithm
C        descriptor.
C
         B    = DPBASE + IXNVP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWACCP ) THEN
C
C        Return the pointer to the intercept acceleration
C        parameters.
C
         B    = DPBASE + IXACCP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWCGP ) THEN
C
C        Return the pointer to the coarse grid structure.
C
         B    = DPBASE + IXCGP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWPIXP ) THEN
C
C        Return the pointer to the pixel data.
C
         B    = DPBASE + IXPIXP
         SIZE = 1

      ELSE IF ( ITEM .EQ. KWREFD ) THEN
C
C        Return the reference surface descriptor.
C
C        We'll need the pointer to this descriptor first.
C
         ADDR = DPBASE + IXREFP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )

         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(REFXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWPRJD ) THEN
C
C        Return the map projection descriptor.
C
C        We'll need the pointer to this descriptor first.
C
         ADDR = DPBASE + IXPRJP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )

         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(PRJXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWITPD ) THEN
C
C        Return the interpolation algorithm descriptor.
C
C        We'll need the pointer to this descriptor first.
C
         ADDR = DPBASE + IXITPP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )

         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(ITPXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWXD ) THEN         
C
C        Return the intercept algorithm descriptor.
C
C        We'll need the pointer to this descriptor 
C        first.
C
         ADDR = DPBASE + IXXP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )
         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(XXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWNVD ) THEN         
C
C        Return the normal vector algorithm descriptor.
C
C        We'll need the pointer to this descriptor first.
C
         ADDR = DPBASE + IXNVP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )
         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(NVXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN


      ELSE IF ( ITEM .EQ. KWACCD ) THEN         
C
C        Return the intercept acceleration descriptor.
C
C        We'll need the pointer to this descriptor first.
C
         ADDR = DPBASE + IXACCP

         CALL DASRDD ( HANDLE, ADDR, ADDR, DPBUFF )         

C
C        We won't know the descriptor's size until we
C        actually read the descriptor, so prepare to
C        read the descriptor of maximum size. We'll
C        adjust the size before returning.
C
         B     = DPBASE + NINT( DPBUFF(1) )
         SIZE  = MAXDSZ

         E     = MIN ( B + SIZE - 1,  MAXADD )

         CALL DASRDD ( HANDLE, B, E, DPBUFF )

C
C        Extract the actual descriptor size from the descriptor.
C
         DSIZE = MIN(  ROOM,  NINT( DPBUFF(AXSZ) )  )

C
C        Move as much of the descriptor as is requested and
C        will fit to the output array.
C
         N  =   MIN ( ROOM, DSIZE - START + 1 )

         CALL MOVED( DPBUFF(START), N, VALUES )

C
C        The outputs are set.
C
         RETURN

 
      ELSE IF ( ITEM .EQ. KWCGD ) THEN         
 
         CALL CHKIN  ( 'DSKD04'                                  )
         CALL SETMSG ( 'Keyword parameter #, which denotes '
     .   //            'the coarse grid data structure, '
     .   //            'is not yet supported.'                   )
         CALL ERRINT ( '#',  ITEM                                )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                     )
         CALL CHKOUT ( 'DSKD04'                                  )
         RETURN

      ELSE

         CALL CHKIN  ( 'DSKD04'                                  )
         CALL SETMSG ( 'Keyword parameter # was not recognized.' )
         CALL ERRINT ( '#',  ITEM                                )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                     )
         CALL CHKOUT ( 'DSKD04'                                  )
         RETURN

      END IF


C
C     The valid range for START is 1:SIZE.
C
      IF (  ( START .LT. 1 ) .OR. ( START .GT. SIZE )  ) THEN

         CALL CHKIN  ( 'DSKD04'                              )
         CALL SETMSG ( 'START must be in the range defined ' //
     .                 'by the size of the data associated ' //
     .                 'with the keyword parameter #, '      //
     .                 'namely 1:#.  Actual value of START ' //
     .                 'was #.'                              )
         CALL ERRINT ( '#',  ITEM                            )
         CALL ERRINT ( '#',  SIZE                            )
         CALL ERRINT ( '#',  START                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKD04'                              )
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

