C$Procedure      ZZEKDE04 ( EK, delete column entry, class 4 )
 
      SUBROUTINE ZZEKDE04 ( HANDLE, SEGDSC, COLDSC, RECPTR )
 
C$ Abstract
C
C     Delete a specified class 4 column entry from an EK record.
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
C     EK
C
C$ Keywords
C
C     PRIVATE
C     UTILITY
C
C$ Declarations
 
 
      INCLUDE 'ekcoldsc.inc'
      INCLUDE 'ekdatpag.inc'
      INCLUDE 'ekrecptr.inc'
      INCLUDE 'eksegdsc.inc'
      INCLUDE 'ektype.inc'
 
      INTEGER               HANDLE
      INTEGER               SEGDSC ( SDSCSZ )
      INTEGER               COLDSC ( CDSCSZ )
      INTEGER               RECPTR
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   File handle.
C     SEGDSC     I   Segment descriptor.
C     COLDSC     I   Column descriptor.
C     RECPTR     I   Record pointer.
C
C$ Detailed_Input
C
C     HANDLE         is a file handle of an EK open for write access.
C
C     SEGDSC         is the descriptor of the segment from which to
C                    delete the specified column entry.
C
C     COLDSC         is the descriptor of the column from which to
C                    delete the specified column entry.
C
C     RECPTR         is a pointer to the record containing the column
C                    entry to delete.
C
C$ Detailed_Output
C
C     None.  See the $Particulars section for a description of the
C     effect of this routine.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If HANDLE is invalid, the error will be diagnosed by routines
C         called by this routine.  The file will not be modified.
C
C     2)  If an I/O error occurs while reading or writing the indicated
C         file, the error will be diagnosed by routines called by this
C         routine.  The file may be corrupted.
C
C$ Files
C
C     See the EK Required Reading for a discussion of the EK file
C     format.
C
C$ Particulars
C
C     This routine operates by side effects:  it deletes a column entry
C     from an EK segment.  The status of the record containing the entry
C     is set to `updated'.  The deleted entry is marked as
C     `uninitialized'.
C
C     The link counts for the pages containing the deleted column entry
C     are decremented.  If the count for a page becomes zero, that page
C     is freed.  If the entry to be deleted is already uninitialized
C     upon entry to this routine, no link counts are modified.  The
C     record containing the entry is still marked `updated' in this
C     case.
C
C     The changes made by this routine to the target EK file become
C     permanent when the file is closed.  Failure to close the file
C     properly will leave it in an indeterminate state.
C
C$ Examples
C
C     See EKDELR.
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.1.0, 07-FEB-2015 (NJB)
C
C        Now uses ERRHAN to insert DAS file name into
C        long error message.
C
C        Deleted unneeded declarations and code.
C
C-    Beta Version 1.0.0, 28-SEP-1995 (NJB)
C
C-&
 
 
C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local variables
C
      INTEGER               BASE
      INTEGER               DATPTR
      INTEGER               NELTS
      INTEGER               NEXT
      INTEGER               NLINKS
      INTEGER               NSEEN
      INTEGER               P
      INTEGER               PTRLOC
      INTEGER               RECNO
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZEKDE04' )
 
C
C     Before trying to actually modify the file, do every error
C     check we can.
C
C     Is this file handle valid--is the file open for paged write
C     access?  Signal an error if not.
C
      CALL ZZEKPGCH ( HANDLE, 'WRITE' )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZEKDE04' )
         RETURN
      END IF
 
C
C     Compute the data pointer location.  If the data pointer is
C     already set to `uninitialized', there's nothing to do.  If
C     the element is null, just set it to `uninitialized'.  The
C     presence of actual data obligates us to clean up, however.
C
      PTRLOC  =  RECPTR + DPTBAS + COLDSC(ORDIDX)
 
      CALL DASRDI ( HANDLE, PTRLOC, PTRLOC, DATPTR )
 
 
      IF ( DATPTR .GT. 0 ) THEN
C
C        Get the element count for the entry.
C
         CALL DASRDI ( HANDLE, DATPTR, DATPTR, NELTS  )
 
C
C        Set the data pointer to indicate the item is uninitialized.
C
         CALL DASUDI ( HANDLE, PTRLOC, PTRLOC, UNINIT )
 
C
C        Find the number of the page containing the column entry.
C
         CALL ZZEKPGPG ( INT, DATPTR, P, BASE )
 
C
C        Look up the forward pointer.  This pointer will be valid
C        if the column entry is continued on another page.
C
         CALL ZZEKGFWD ( HANDLE, INT, P, NEXT )
 
C
C        Get the link count for the current page.  If we have more
C        than one link to the page, decrement the link count.  If
C        we're down to one link, this deletion will finish off the
C        page:  we'll deallocate it.
C
         CALL ZZEKGLNK ( HANDLE, INT, P, NLINKS )
 
         IF ( NLINKS .GT. 1 ) THEN
 
            CALL ZZEKSLNK ( HANDLE, INT, P, NLINKS-1 )
 
         ELSE
C
C           If we removed the last item from the page, we can delete
C           the page.  ZZEKDPS adjusts the segment's metadata
C           to reflect the deallocation.
C
            CALL ZZEKDPS  ( HANDLE, SEGDSC, INT, P )
 
         END IF
 
 
         NSEEN = MIN (  NELTS,  (BASE+IPSIZE-DATPTR)  )
 
 
         DO WHILE (  ( NSEEN .LT. NELTS ) .AND. ( .NOT. FAILED() )  )
C
C           The column entry is continued on the page indicated by
C           NEXT.
C
C           Get the link count for the current page.  If we have more
C           than one link to the page, decrement the link count.  If
C           we're down to one link, this deletion will finish off the
C           page:  we'll deallocate it.
C
            P  =  NEXT
 
            CALL ZZEKGFWD ( HANDLE, INT, P, NEXT   )
            CALL ZZEKGLNK ( HANDLE, INT, P, NLINKS )
 
 
            IF ( NLINKS .GT. 1 ) THEN
 
               CALL ZZEKSLNK ( HANDLE, INT, P, NLINKS-1 )
 
            ELSE
C
C              If we removed the last item from the page, we can delete
C              the page.  ZZEKDPS adjusts the segment's metadata
C              to reflect the deallocation.
C
               CALL ZZEKDPS  ( HANDLE, SEGDSC, INT, P )
 
            END IF
 
            NSEEN  =  MIN (  NELTS,  ( IPSIZE + NSEEN )  )
 
         END DO
 
 
      ELSE IF ( DATPTR .EQ. NULL ) THEN
C
C        Mark the entry as `uninitialized'.
C
         CALL DASUDI ( HANDLE, PTRLOC, PTRLOC, UNINIT )
 
 
      ELSE IF ( DATPTR .NE. UNINIT ) THEN
C
C        UNINIT was the last valid possibility.  The data pointer is
C        corrupted.
C
         CALL SETMSG ( 'Data pointer is corrupted. SEGNO = #; '  //
     .                 'COLIDX =  #; RECNO = #; EK = #'          )
         CALL ERRINT ( '#',  SEGDSC(SNOIDX)                      )
         CALL ERRINT ( '#',  COLDSC(ORDIDX)                      )
         CALL ERRINT ( '#',  RECNO                               )
         CALL ERRHAN ( '#',  HANDLE                              )
         CALL SIGERR ( 'SPICE(BUG)'                              )
         CALL CHKOUT ( 'ZZEKDE04'                                )
         RETURN
 
      END IF
 
C
C     Set the record's status to indicate that this record is updated.
C
      CALL DASUDI ( HANDLE, RECPTR+STAIDX, RECPTR+STAIDX, UPDATE )
 
 
      CALL CHKOUT ( 'ZZEKDE04' )
      RETURN
      END
