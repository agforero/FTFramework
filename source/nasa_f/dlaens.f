C$Procedure DLAENS ( DLA, end new segment )
 
      SUBROUTINE DLAENS ( HANDLE )
      IMPLICIT NONE
 
C$ Abstract
C
C     End a new segment in a DLA file.
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
C
C$ Declarations

      INCLUDE 'dla.inc'     

      INTEGER               HANDLE
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of open DLA file.
C
C$ Detailed_Input
C
C     HANDLE      is the integer handle associated with the DLA file to
C                 be updated.  This handle is used to identify the file
C                 in subsequent calls to other DLA or DAS routines.
C
C                 The DLA file must be open for write access. A new DLA
C                 segment is completed in the indicated file.  The file
C                 is left open, since data may be written to the file
C                 following a call to this routine.
C
C$ Detailed_Output
C
C     None.  See the Particulars and Examples header sections for
C     a description of the actions performed by this routine.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the input file handle does not refer to a DAS file that is
C        open for write access, the error will be diagnosed by routines
C        in the call tree of this routine.
C
C     2) If an error occurs while reading or writing to the DLA file,
C        the error will be diagnosed by routines in the call tree of
C        this routine.
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
C     This routine supports creation of a DLA segment.  DLA segments
C     are created by appending data to the DAS integer, double 
C     precision, and character address spaces of a DLA file.  The new
C     segment's descriptor is located immediately before the integer
C     component of the segment's data.
C
C     When a new segment is added to a DLA file, the segment is
C     inserted into the file's doubly linked segment list.  If the new
C     segment is the first, the DLA file's first and last list entry
C     pointers are updated to point to the new segment; specifically,
C     these pointers point to the first integer of the new segment's
C     descriptor.  The backward pointer of the new segment is set to
C     null in this case.
C
C     If the new segment is not the first, the DLA file's list end
C     pointer is updated to point to the new segment, and the forward
C     pointer of the previous segment also is updated to point to the
C     first integer of the new segment's descriptor. The backward
C     pointer of the new segment points to to point to the first
C     integer of the previous segment's descriptor.
C
C     The normal sequence of operations required to create a DLA 
C     segment is as follows:
C
C        Call DLAOPN to create a new, empty DLA file.
C
C        For each segment to be created,
C
C           Call DLABNS to begin a segment.
C
C           Use the DAS "add" and "update" routines to populate
C           the segment with data.
C
C           Call DLAENS to end the segment.
C
C        Call DASCLS to segregate and close the DLA file.
C
C
C$ Examples
C
C     1) Create a DLA file containing one segment; the segment
C        contains character, double precision, and integer data.
C        After writing and closing the file, open the file for
C        read access; dump the data to standard output.
C
C 
C              PROGRAM EX1
C              IMPLICIT NONE
C
C              INCLUDE 'dla.inc'
C
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         DLA
C              PARAMETER           ( DLA    = 'test.dla' )
C
C              INTEGER               IFNLEN
C              PARAMETER           ( IFNLEN =  60 )
C
C              INTEGER               LNSIZE
C              PARAMETER           ( LNSIZE =  80 )
C
C              INTEGER               MAXC
C              PARAMETER           ( MAXC   =  5 )
C
C              INTEGER               MAXD
C              PARAMETER           ( MAXD   =  50 )
C
C              INTEGER               MAXI
C              PARAMETER           ( MAXI   =  100 )
C
C        C
C        C     Local variables
C        C
C              CHARACTER*(LNSIZE)    CVALS   ( MAXC )
C              CHARACTER*(LNSIZE)    CVALS2  ( MAXC )
C              CHARACTER*(IFNLEN)    IFNAME
C
C              DOUBLE PRECISION      DVALS   ( MAXD )
C              DOUBLE PRECISION      DVALS2  ( MAXD )
C
C              INTEGER               BASE
C              INTEGER               DESCR   ( DLADSZ )
C              INTEGER               HANDLE
C              INTEGER               I
C              INTEGER               IVALS   ( MAXI )
C              INTEGER               IVALS2  ( MAXI )
C              INTEGER               J
C              INTEGER               K
C              INTEGER               N
C              INTEGER               NCOMCH
C
C              LOGICAL               FOUND
C
C        C
C        C     Set the internal file name.  Don't reserve characters in
C        C     the DAS comment area.
C        C
C              IFNAME = 'Example DLA file for testing'
C              NCOMCH = 0
C
C        C
C        C     Open a new DLA file.
C        C
C              CALL DLAOPN ( DLA, 'DLA', IFNAME, NCOMCH, HANDLE )
C
C        C
C        C     Begin a new segment.
C        C
C              CALL DLABNS ( HANDLE )
C
C        C
C        C     Add character data to the segment.
C        C
C              DO I = 1, MAXC
C
C                 DO J = 1, LNSIZE
C
C                    K = MOD( J+I-1, 10 )
C
C                    CALL INTSTR ( K,  CVALS(I)(J:J) )
C
C                 END DO
C
C              END DO
C
C              CALL DASADC ( HANDLE, MAXC*LNSIZE, 1, LNSIZE, CVALS )
C
C        C
C        C     Add integer and double precision data to the segment.
C        C
C              DO I = 1, MAXI
C                 IVALS(I) = I
C              END DO
C
C              CALL DASADI ( HANDLE, MAXI, IVALS )
C
C              DO I = 1, MAXD
C                 DVALS(I) = I
C              END DO
C
C              CALL DASADD ( HANDLE, MAXD, DVALS )
C
C        C
C        C     End the segment.
C        C
C              CALL DLAENS ( HANDLE )
C
C        C
C        C     Close the file.  The routine DASCLS flushes the DAS
C        C     buffers and segregates the file before closing it.
C        C
C              CALL DASCLS ( HANDLE )
C
C        C
C        C     Now read the file and check the data.
C        C
C              CALL DASOPR ( DLA, HANDLE )
C
C        C
C        C     Obtain the segment descriptor for the sole segment 
C        C     in the file. We need not check the found flag 
C        C     in this case because we know there is one segment
C        C     in the file.
C        C
C              CALL DLABFS ( HANDLE, DESCR, FOUND )
C
C        C
C        C     Fetch character data from the segment.  Obtain the 
C        C     base address of the character data and the
C        C     character count from the descriptor.
C        C
C              BASE = DESCR(CBSIDX)
C              N    = DESCR(CSZIDX)
C
C              CALL DASRDC ( HANDLE, BASE+1, BASE+N, 1, LNSIZE, CVALS2 )
C
C        C
C        C     Display the character data.
C        C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Character array'
C
C              DO I = 1, N/LNSIZE
C                 WRITE (*,*) CVALS2(I)
C              END DO
C
C        C
C        C     Fetch and display the integer and double precision data.
C        C
C              BASE = DESCR(IBSIDX)
C              N    = DESCR(ISZIDX)
C
C              CALL DASRDI( HANDLE, BASE+1, BASE+N, IVALS2 )
C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Integer array'
C              WRITE (*,*) IVALS2
C
C
C              BASE = DESCR(DBSIDX)
C              N    = DESCR(DSZIDX)
C
C              CALL DASRDD( HANDLE, BASE+1, BASE+N, DVALS2 )
C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Double precision array'
C              WRITE (*,*) DVALS2
C
C        C
C        C     Close the file.  This step is unnecessary in this 
C        C     program, but is a good practice in general 
C        C     because closing the file frees resources.
C        C
C              CALL DASCLS ( HANDLE )
C
C              END
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
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Updated version info.
C
C        08-OCT-2009 (NJB)
C
C           Updated header.
C
C        11-FEB-2005 (NJB)
C
C-&
 
C$ Index_Entries
C
C     end new segment in dla file
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local parameters
C
      
C
C     Local variables
C
      INTEGER               DESCR  ( DLADSZ )
      INTEGER               LASTC
      INTEGER               LASTD
      INTEGER               LASTI
      INTEGER               THIS

C
C     Standard SPICE error handling.
C     
      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'DLAENS' )
C
C     Now that the segment has been written, our only task is
C     to update the corresponding segment descriptor to reflect
C     the sizes of each component. 
C
C     Look up the pointer to the last DLA segment descriptor in the
C     file.  Then look up the segment descriptor itself.
C
      CALL DASRDI ( HANDLE, LLEIDX, LLEIDX,        THIS  )
      CALL DASRDI ( HANDLE, THIS,   THIS+DLADSZ-1, DESCR )

C
C     Find the last DAS logical addresses in use for each data type.
C
      CALL DASLLA ( HANDLE, LASTC, LASTD, LASTI )

C
C     Set the component sizes in the descriptor. The sizes are easily
C     computed from the last addresses in use and the component base
C     addresses already stored in the segment descriptor.
C
      DESCR ( ISZIDX ) = LASTI - DESCR(IBSIDX )
      DESCR ( DSZIDX ) = LASTD - DESCR(DBSIDX )
      DESCR ( CSZIDX ) = LASTC - DESCR(CBSIDX )

C
C     Update the descriptor in the file.
C
      CALL DASUDI ( HANDLE, THIS, THIS+DLADSZ-1, DESCR )

C
C     Leave the file open.  The file is now ready for the
C     addition of a new segment.
C
      CALL CHKOUT ( 'DLAENS' )
      RETURN
      END

      
