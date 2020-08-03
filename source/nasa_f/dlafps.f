C$Procedure DLAFPS ( DLA, find previous segment )
 
      SUBROUTINE DLAFPS ( HANDLE, DESCR, PRVDSC, FOUND )
      IMPLICIT NONE
 
C$ Abstract
C
C     Find the segment preceding a specified segment in a DLA file.
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
C     DLA
C
C$ Keywords
C
C     DAS
C     DLA
C     FILES
C     SEARCH
C
C$ Declarations

      INCLUDE 'dla.inc'     

      INTEGER               HANDLE
      INTEGER               DESCR  ( * )
      INTEGER               PRVDSC ( * )
      LOGICAL               FOUND
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of open DLA file.
C     DESCR      I   Descriptor of a segment in DLA file.
C     PRVDSC     O   Descriptor of previous segment in DLA file.
C     FOUND      O   Flag indicating whether a segment was found.
C
C$ Detailed_Input
C
C     HANDLE      is the integer handle associated with the file to be
C                 searched. This handle is used to identify the file in
C                 subsequent calls to other DLA or DAS routines.
C
C     DESCR       is the descriptor of a DLA segment in the file
C                 associated with HANDLE. The descriptor of the
C                 segment preceding DESCR is sought.
C
C                 The segment descriptor layout is:
C
C                  +---------------+
C                  | BACKWARD PTR  | Linked list backward pointer
C                  +---------------+
C                  | FORWARD PTR   | Linked list forward pointer
C                  +---------------+
C                  | BASE INT ADDR | Base DAS integer address
C                  +---------------+
C                  | INT COMP SIZE | Size of integer segment component
C                  +---------------+
C                  | BASE DP ADDR  | Base DAS d.p. address
C                  +---------------+
C                  | DP COMP SIZE  | Size of d.p. segment component
C                  +---------------+
C                  | BASE CHR ADDR | Base DAS character address
C                  +---------------+
C                  | CHR COMP SIZE | Size of character segment component
C                  +---------------+
C
C$ Detailed_Output
C
C     PRVDSC      is the descriptor of the DLA segment preceding the
C                 segment associated with the input argument DESCR.
C
C                 PRVDSC is valid only if the output argument FOUND is
C                 .TRUE.  
C
C
C     FOUND       is a logical flag indicating whether the previous
C                 segment was found.  FOUND has the value .TRUE. if the
C                 segment was found; otherwise FOUND is .FALSE.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the input file handle is invalid, the error will be
C        diagnosed by routines in the call tree of this routine.
C
C     2) If an error occurs while reading the DLA file, the error 
C        will be diagnosed by routines in the call tree of this
C        routine.
C
C$ Files
C
C     See description of input argument HANDLE.
C
C$ Particulars
C
C     DLA files are built using the DAS low-level format; DLA files are
C     a specialized type of DAS file in which data are organized as a
C     doubly linked list of segments.  Each segment's data belong to
C     contiguous components of character, double precision, and integer
C     type.
C
C     This routine supports backward traversal of a DLA file's segment
C     list.  A backward traversal may be started from any segment in
C     the file; it is not necessary to call DLABBS first.  The role of
C     DLABBS is simply to return the descriptor of the last segment in
C     the file.
C
C$ Examples
C
C     1)  Open a DLA file for read access, traverse the segment
C         list from back to front, and display segment address
C         and size attributes.
C
C
C                  PROGRAM EX1
C                  IMPLICIT NONE
C
C                  INCLUDE 'dla.inc'
C 
C            C
C            C     Local parameters
C            C
C                  INTEGER               FILSIZ
C                  PARAMETER           ( FILSIZ = 255 )
C
C            C
C            C     Local variables
C            C
C                  CHARACTER*(FILSIZ)    FNAME
C
C                  INTEGER               CURRNT ( DLADSZ )
C                  INTEGER               DESCR  ( DLADSZ )
C                  INTEGER               HANDLE
C                  INTEGER               NSEGS
C                  INTEGER               SEGNO
C
C                  LOGICAL               FOUND
C
C            C
C            C     Prompt for the name of the file to search.
C            C
C                  CALL PROMPT ( 'Name of DLA file > ', FNAME )
C
C            C
C            C     Open the DLA file for read access.  Since DLA
C            C     files use the DAS architecture, we can use DAS
C            C     routines to open and close the file.
C            C
C                  CALL DASOPR ( FNAME, HANDLE )
C
C            C
C            C     Count the segments in the file; this allows us
C            C     to label the segments in our display.
C            C
C                  NSEGS = 0
C                  CALL DLABBS ( HANDLE, DESCR, FOUND )
C
C                  DO WHILE ( FOUND ) 
C
C                     NSEGS = NSEGS + 1
C                     CALL MOVEI  ( DESCR,  DLADSZ, CURRNT       )
C                     CALL DLAFPS ( HANDLE, CURRNT, DESCR, FOUND )
C
C                  END DO
C
C            C
C            C     Begin a backward search.  Let DESCR contain
C            C     the descriptor of the last segment.
C            C
C                  SEGNO = NSEGS + 1
C
C                  CALL DLABBS ( HANDLE, DESCR, FOUND )
C
C                  DO WHILE ( FOUND )
C            C
C            C        Display the contents of the current segment
C            C        descriptor.
C            C
C                     SEGNO = SEGNO - 1
C
C                     WRITE (*,*) ' '
C                     WRITE (*,*) ' '
C                     WRITE (*,*) 'Segment number = ', SEGNO
C                     WRITE (*,*) ' '
C                     WRITE (*,*) 'Backward segment pointer         = ',
C                 .               DESCR(BWDIDX)
C                     WRITE (*,*) 'Forward segment pointer          = ',
C                 .               DESCR(FWDIDX)
C                     WRITE (*,*) 'Character component base address = ',
C                 .               DESCR(CBSIDX)
C                     WRITE (*,*) 'Character component size         = ',
C                 .               DESCR(CSZIDX)
C                     WRITE (*,*) 'D.p. base address                = ',
C                 .               DESCR(DBSIDX)
C                     WRITE (*,*) 'D.p. component size              = ',
C                 .               DESCR(DSZIDX)
C                     WRITE (*,*) 'Integer base address             = ',
C                 .               DESCR(IBSIDX)
C                     WRITE (*,*) 'Integer component size           = ',
C                 .               DESCR(ISZIDX)
C                     WRITE (*,*) ' '
C
C            C
C            C        Find the previous segment. 
C            C                       
C            C        To avoid using DESCR as both input and output
C            C        in the following call (this use is not allowed 
C            C        by the ANSI Fortran 77 standard), we copy DESCR
C            C        into the variable CURRNT.  We then find the
C            C        segment preceding CURRNT.
C            C         
C                     CALL MOVEI  ( DESCR,  DLADSZ, CURRNT       )
C                     CALL DLAFPS ( HANDLE, CURRNT, DESCR, FOUND )
C
C                  END DO
C
C            C
C            C     Close the file using the DAS close routine.
C            C                       
C                  CALL DASCLS ( HANDLE )
C
C                  END
C
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 06-JAN-2017 (NJB)
C
C     Original version 21-APR-2010 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find previous segment in dla file
C
C-&
 

      
C
C     SPICELIB functions
C
      LOGICAL               RETURN


C
C     Local variables
C
      INTEGER               BWD

C
C     Standard SPICE error handling.
C     
      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'DLAFPS' )

C
C     Nothing found yet.
C
      FOUND = .FALSE.

C
C     Extract the backward pointer from the segment descriptor.
C
      BWD = DESCR( BWDIDX )

      IF ( BWD .EQ. NULPTR ) THEN
C
C        There is no segment preceding the input segment.
C
         CALL CHKOUT ( 'DLAFPS' )
         RETURN

      END IF
     
C
C     Look up the previous descriptor
C
      CALL DASRDI ( HANDLE, BWD, BWD+DLADSZ-1, PRVDSC )

      FOUND = .TRUE.

      CALL CHKOUT ( 'DLAFPS' )
      RETURN
      END

      
