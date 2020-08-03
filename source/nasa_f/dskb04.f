C$Procedure DSKB04 ( DSK, fetch bookkeeping data, type 4 )

      SUBROUTINE DSKB04 ( HANDLE,  DLADSC,  MXITPD,  MXXD,
     .                    MXNVD,   NR,      NC,      CO1PSZ, 
     .                    CO2PSZ,  CENT1,   CENT2,   NULLOK,
     .                    NULVAL,  ITPDSC,  XDSC,    NVDSC  )
 
C$ Abstract
C
C     Return bookkeeping parameters from a specified type 4
C     DSK segment.
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
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               MXITPD
      INTEGER               MXXD
      INTEGER               MXNVD
      INTEGER               NR
      INTEGER               NC
      DOUBLE PRECISION      CO1PSZ
      DOUBLE PRECISION      CO2PSZ
      DOUBLE PRECISION      CENT1
      DOUBLE PRECISION      CENT2
      LOGICAL               NULLOK
      DOUBLE PRECISION      NULVAL
      DOUBLE PRECISION      ITPDSC ( MXITPD )
      DOUBLE PRECISION      XDSC   ( MXXD )
      DOUBLE PRECISION      NVDSC  ( MXNVD )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.
C     DLADSC     I   DLA descriptor.
C     MXITPD     I   Maximum size of interpolation descriptor.
C     MXXD       I   Maximum size of intercept descriptor.
C     MXNVD      I   Maximum size of normal vector descriptor.
C     NR         O   Number of grid rows.
C     NC         O   Number of grid columns.
C     CO1PSZ     O   Coordinate 1 pixel size.
C     CO2PSZ     O   Coordinate 2 pixel size.
C     CENT1      O   Coordinate 1 of grid center.
C     CENT2      O   Coordinate 2 of grid center.
C     NULLOK     O   Flag indicating whether null values are allowed.
C     NULVAL     O   Null value.
C     ITPDSC     O   Interpolation descriptor.
C     XDSC       O   Intecept descriptor.
C     NVDSC      O   Normal vector descriptor.
C
C$ Detailed_Input
C
C     HANDLE         is the handle of a DSK file containing a type 4
C                    segment from which data are to be fetched.
C
C     DLADSC         is the DLA descriptor associated with the segment
C                    from which data are to be fetched. 
C
C     MXITPD         is the maximum size of the interpolation algorithm
C                    descriptor.
C   
C     MXXD           is the maximum size of the ray-surface intercept
C                    algorithm descriptor.
C   
C     MXNVD          is the maximum size of the surface normal vector
C                    algorithm descriptor.
C
C$ Detailed_Output
C
C     NR
C     NC,            are, respectively, the number of rows and columns
C                    in the segment's pixel grid.
C
C     CO1PSZ,        
C     CO2PSZ         are, respectively, the sizes of the grid's 
C                    pixels in the directions of the first and
C                    second coordinates. 
C
C                    The units of the pixel sizes are those of
C                    the corresponding coordinates.
C
C     CENT1,
C     CENT2          are, respectively, the coordinates of the pixel
C                    grid's center. Note that the pixel grid need not
C                    be centered at the center of the coordinate
C                    rectangle given by the segment's DSK descriptor.
C
C     NULLOK         is a logical flag that is .TRUE. if and only if
C                    null pixel values are permitted.
C
C     NULVAL         is the integer value representing null values.
C                    NULVAL must be expressible as a 16-bit, signed
C                    integer. 
C
C     ITPDSC         is the interpolation algorithm descriptor. See
C                    dsk04.inc for details.
C
C     XDSC           is the ray-surface intercept algorithm descriptor.
C                    See dsk04.inc for details.
C
C     NVDSC          is the surface normal vector algorithm descriptor.
C                    See dsk04.inc for details.
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
C     1)  If the input handle is invalid, the error will be diagnosed
C         by routines in the call tree of this routine.
C
C     2)  If a file read error occurs, the error will be diagnosed by
C         routines in the call tree of this routine.
C
C     3)  If the input DLA descriptor is invalid, the effect of this
C         routine is undefined. The error *may* be diagnosed by
C         routines in the call tree of this routine, but there are no
C         guarantees.
C
C     4)  If any output descriptor array is too small to hold the
C         corresponding descriptor, the error SPICE(BUFFERTOOSMALL) is
C         signaled.
C
C     5)  If the DSK type 4 segment format version is not 3, the error
C         SPICE(VERSIONMISMATCH) is signaled.
C
C     6)  If the segment's data type is not 4, the error
C         SPICE(BADDATATYPE) is signaled.
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
C-    DSKBRIEF Version 1.0.0, 06-OCT-2016 (NJB)
C
C        Changed order of arguments!!!
C
C        Removed unused variables.
C
C
C        Version 2.0.0 20-SEP-2012 (NJB)
C        Version 1.0.0 09-AUG-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch bookkeeping data from a type_4_dsk segment
C
C-&
 
 
C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      DPBUFF ( MAXDSZ )
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )

      INTEGER               DSKFMT
      INTEGER               DTYPE
      INTEGER               N



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN( 'DSKB04' )

      CALL DSKGD ( HANDLE, DLADSC, DSKDSC )

C
C     Check the DSK format version.
C
      CALL DSKD04 ( HANDLE, DLADSC, KWFMTV, 1, 1, N, DPBUFF )

      DSKFMT = NINT ( DPBUFF(1) )

      IF ( DSKFMT .NE. 3 ) THEN

         CALL SETMSG ( 'DSK format version was expected to '
     .   //            'be 3 but was #. Only version 3 is '
     .   //            'supported by this software version.' )
         CALL ERRINT ( '#', DSKFMT                           )
         CALL SIGERR( 'SPICE(VERSIONMISMATCH)'               )
         CALL CHKOUT( 'DSKB04'                               )
         RETURN

      END IF

C
C     Check the data type.
C
      DTYPE = NINT( DSKDSC(TYPIDX) )

         IF ( DTYPE .NE. 4 ) THEN
            
            CALL SETMSG( 'Input segment must be type 4 but is type '
     .      //           '#. HANDLE = #; DLA base addresses are '
     .      //           'INTEGER: #; DP: #; CHAR: #.'             )
            CALL ERRINT( '#',  DTYPE                               )
            CALL ERRINT( '#',  HANDLE                              )
            CALL ERRINT( '#',  DLADSC(IBSIDX)                      )
            CALL ERRINT( '#',  DLADSC(DBSIDX)                      )
            CALL ERRINT( '#',  DLADSC(CBSIDX)                      )
            CALL SIGERR( 'SPICE(BADDATATYPE)'                      )
            CALL CHKOUT( 'DSKB04'                                  )
            RETURN
         
         END IF

C
C     Fetch the fixed-size double precision data items.
C     
      CALL DSKD04( HANDLE, DLADSC, KWNR,   1, 1, N, DPBUFF )
      NR = NINT( DPBUFF(1) )

      CALL DSKD04( HANDLE, DLADSC, KWNC,   1, 1, N, DPBUFF )
      NC = NINT( DPBUFF(1) )

      CALL DSKD04( HANDLE, DLADSC, KWPSZ1, 1, 1, N, CO1PSZ )
      CALL DSKD04( HANDLE, DLADSC, KWPSZ2, 1, 1, N, CO2PSZ )
      CALL DSKD04( HANDLE, DLADSC, KWCTR1, 1, 1, N, CENT1  )
      CALL DSKD04( HANDLE, DLADSC, KWCTR2, 1, 1, N, CENT2  )

C
C     Map the d.p. null flag to a type LOGICAL value.
C     Note that the identifier TRUE below refers to a d.p. 
C     parameter.
C
      CALL DSKD04( HANDLE, DLADSC, KWNLOK, 1, 1, N, DPBUFF )

      NULLOK = DPBUFF(1) .EQ. TRUE

C
C     Fetch the null value marker itself.
C         
      CALL DSKD04( HANDLE, DLADSC, KWNULL, 1, 1, N, NULVAL )

C
C     Fetch the interpolation algorithm descriptor. The caller
C     is responsible for providing enough room.
C
      CALL DSKD04( HANDLE, DLADSC, KWITPD, 1, MAXDSZ, N, DPBUFF )
      
      IF ( MXITPD .LT. N ) THEN

         CALL SETMSG( 'The interpolation descriptor size is #; only '
     .   //           '# elements were provided in the output '
     .   //           'array.'                                       )
         CALL ERRINT( '#',  N                                        )
         CALL ERRINT( '#',  MXITPD                                   )
         CALL SIGERR( 'SPICE(BUFFERTOOSMALL)'                        )
         CALL CHKOUT( 'DSKB04'                                       )
         RETURN

      END IF

      CALL MOVED( DPBUFF, N, ITPDSC )


C
C     Fetch the intercept algorithm descriptor. The caller
C     is responsible for providing enough room.
C
      CALL DSKD04( HANDLE, DLADSC, KWXD, 1, MAXDSZ, N, DPBUFF )
      
      IF ( MXXD .LT. N ) THEN

         CALL SETMSG( 'The intercept descriptor size is #; only '
     .   //           '# elements were provided in the output '
     .   //           'array.'                                       )
         CALL ERRINT( '#',  N                                        )
         CALL ERRINT( '#',  MXXD                                     )
         CALL SIGERR( 'SPICE(BUFFERTOOSMALL)'                        )
         CALL CHKOUT( 'DSKB04'                                       )
         RETURN

      END IF

      CALL MOVED( DPBUFF, N, XDSC )

C
C     Fetch the normal vector algorithm descriptor. The caller
C     is responsible for providing enough room.
C
      CALL DSKD04( HANDLE, DLADSC, KWNVD, 1, MAXDSZ, N, DPBUFF )
      
      IF ( MXNVD .LT. N ) THEN

         CALL SETMSG( 'The intercept descriptor size is #; only '
     .   //           '# elements were provided in the output '
     .   //           'array.'                                       )
         CALL ERRINT( '#',  N                                        )
         CALL ERRINT( '#',  MXNVD                                    )
         CALL SIGERR( 'SPICE(BUFFERTOOSMALL)'                        )
         CALL CHKOUT( 'DSKB04'                                       )
         RETURN

      END IF

      CALL MOVED( DPBUFF, N, NVDSC )


      CALL CHKOUT( 'DSKB04' )
      RETURN
      END 

