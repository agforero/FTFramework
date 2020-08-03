C$Procedure ZZDAFGSR ( Private --- DAF Get Summary/Descriptor Record )
 
      SUBROUTINE ZZDAFGSR ( HANDLE, RECNO, ND, NI, DPREC, FOUND )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Read a summary/descriptor record from a DAF.
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
C     PRIVATE
C
C$ Declarations
 
      IMPLICIT NONE
 
      INCLUDE              'zzddhman.inc'
 
      INTEGER               HANDLE
      INTEGER               RECNO
      INTEGER               ND
      INTEGER               NI
      DOUBLE PRECISION      DPREC  ( * )
      LOGICAL               FOUND
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of the DAF.
C     RECNO      I   Record number.
C     ND         I   Number of double precision components in a summary.
C     NI         I   Number of integer components in a summary.
C     DPREC      O   Contents of the record.
C     FOUND      O   Logical indicating whether the record was found.
C
C$ Detailed_Input
C
C     HANDLE     is the handle associated with the DAF.
C
C     RECNO      is the record number of a particular summary record
C                within the DAF, whose contents are to be read.
C     ND,
C     NI         are the number of double precision and integer
C                components, respectively, in each array summary
C                in the specified file.
C
C$ Detailed_Output
C
C     DPREC      contains the contents of the specified record from
C                the DAF associated with HANDLE, properly translated
C                for use on the native environment.
C
C     FOUND      is TRUE when the specified record is found, and is
C                FALSE otherwise.
C
C$ Parameters
C
C     None.
C
C$ Files
C
C     This routine reads data from the DAF associated with HANDLE.
C     This action may result in connecting a logical unit to the
C     file, if the handle manager has rotated the file out of the
C     unit table.
C
C$ Exceptions
C
C     1) SPICE(HANDLENOTFOUND) is signaled if HANDLE can not be
C        found in the set of loaded handles.
C
C     2) Routines in the call tree of this routine may trap and
C        signal errors.
C
C$ Particulars
C
C     This routine reads summary records of double precision
C     numbers which contain integers packed through an EQUIVALENCE
C     statement from native and supported non-native DAFs.
C
C     The size of the character buffer and the number of records
C     read may have to change to support new environments.  As of
C     the original release of this routine, all systems currently
C     supported have a 1 kilobyte record length.
C
C$ Examples
C
C     See DAFGSR for sample usage.
C
C$ Restrictions
C
C     1) Numeric data when read as characters from a file preserves
C        the bit patterns present in the file in memory.
C
C     2) A record of double precision data is at most 1024 characters
C        in length.
C
C     3) DPREC has enough space to store 128 double precision numbers.
C
C     4) Characters a byte-sized, 8 characters constitute a double
C        precision number, and 4 characters an integer.
C
C$ Author_and_Institution
C
C     F.S. Turner     (JPL)
C
C$ Literature_References
C
C     None.
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 12-NOV-2001 (FST)
C
C
C-&
 
C
C     SPICELIB Functions
C
      INTEGER               ISRCHC
 
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local Parameters
C
C     Length in bytes of double precision numbers and integers.
C
      INTEGER               DPLEN
      PARAMETER           ( DPLEN = 8       )
 
      INTEGER               INLEN
      PARAMETER           ( INLEN = DPLEN/2 )
 
C
C     Local Variables
C
      CHARACTER*(CBFSIZ)    CHRBUF
      CHARACTER*(FILEN)     FNAME
      CHARACTER*(STRSIZ)    STRBFF ( NUMBFF )
      CHARACTER*(STRSIZ)    TMPSTR
 
      DOUBLE PRECISION      DPBUF  ( 128 )
 
      INTEGER               CINDEX
      INTEGER               DINDEX
      INTEGER               I
      INTEGER               IAMH
      INTEGER               IARCH
      INTEGER               IBFF
      INTEGER               INBUF  ( 256 )
      INTEGER               IOSTAT
      INTEGER               LEFT
      INTEGER               LUN
      INTEGER               NATBFF
      INTEGER               NSUM
      INTEGER               SUMSIZ
 
      LOGICAL               FIRST
      LOGICAL               LOCFND
 
C
C     Equivalence DPBUF to INBUF.
C
      EQUIVALENCE         ( DPBUF, INBUF )
 
C
C     Saved Variables
C
      SAVE                  FIRST
      SAVE                  NATBFF
      SAVE                  STRBFF
 
C
C     Data Statements
C
      DATA                  FIRST  / .TRUE. /
      DATA                  NATBFF / 0      /
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'ZZDAFGSR' )
      END IF
 
C
C     Perform some initialization tasks.
C
      IF ( FIRST ) THEN
 
C
C        Populate STRBFF, the buffer that contains the labels
C        for each binary file format.
C
         DO I = 1, NUMBFF
            CALL ZZDDHGSD ( 'BFF', I, STRBFF(I) )
         END DO
 
C
C        Fetch the native binary file format and determine its
C        integer code.
C
         CALL ZZPLATFM ( 'FILE_FORMAT', TMPSTR )
         CALL UCASE    ( TMPSTR,        TMPSTR )
 
         NATBFF = ISRCHC ( TMPSTR, NUMBFF, STRBFF )
 
         IF ( NATBFF .EQ. 0 ) THEN
 
            CALL SETMSG ( 'The binary file format, ''#'', is not '
     .      //            'supported by this version of the toolkit. '
     .      //            'This is a serious problem, contact NAIF.'   )
            CALL ERRCH  ( '#', TMPSTR                                  )
            CALL SIGERR ( 'SPICE(BUG)'                                 )
            CALL CHKOUT ( 'ZZDAFGSR'                                   )
            RETURN
 
         END IF
 
C
C        Do not perform initialization tasks again.
C
         FIRST = .FALSE.
 
      END IF
 
C
C     Assume the data record will not be found, until it has been read
C     from the file, and if necessary, successfully translated.
C
      FOUND = .FALSE.
 
C
C     Retrieve information regarding the file from the handle manager.
C     The value of IARCH is not a concern, since this is a DAF routine
C     all values passed into handle manager entry points will have
C     'DAF' as their architecture arguments.
C
      CALL ZZDDHNFO ( HANDLE, FNAME, IARCH, IBFF, IAMH, LOCFND )
 
      IF ( .NOT. LOCFND ) THEN
         CALL SETMSG ( 'Unable to locate file associated with '
     .   //            'HANDLE, #.  The most likely cause of this '
     .   //            'is the file that you are trying to read '
     .   //            'has been closed.'                           )
         CALL ERRINT ( '#', HANDLE                                  )
         CALL SIGERR ( 'SPICE(HANDLENOTFOUND)'                      )
         CALL CHKOUT ( 'ZZDAFGSR'                                   )
         RETURN
      END IF
 
C
C     Now get a logical unit for the handle.  Check FAILED()
C     in case an error occurs.
C
      CALL ZZDDHHLU ( HANDLE, 'DAF', .FALSE., LUN )
 
      IF ( FAILED() ) THEN
         FOUND = .FALSE.
         CALL CHKOUT ( 'ZZDAFGSR' )
         RETURN
      END IF
 
C
C     Branch based on whether the binary file format is native
C     or not.  Only supported formats can be opened by ZZDDHOPN,
C     so no check of IBFF is required.
C
      IF ( IBFF .EQ. NATBFF ) THEN
 
C
C        In the native case, just read the array of double precision
C        numbers from the file.  The packed integers will be
C        processed properly by the READ.
C
         READ ( UNIT   = LUN,
     .          REC    = RECNO,
     .          IOSTAT = IOSTAT ) ( DPBUF(I), I = 1, 128 )
 
C
C        Since this routine does not signal any IOSTAT based
C        errors, return if a non-zero value is assigned to IOSTAT.
C
         IF ( IOSTAT .NE. 0 ) THEN
            CALL CHKOUT ( 'ZZDAFGSR' )
            RETURN
         END IF
 
C
C     Process the non-native binary file format case.
C
      ELSE
 
C
C        Read the record as characters.
C
         READ ( UNIT   = LUN,
     .          REC    = RECNO,
     .          IOSTAT = IOSTAT ) CHRBUF
 
C
C        Again, since this routine does not signal any IOSTAT
C        based errors, return if one occurs.
C
         IF ( IOSTAT .NE. 0 ) THEN
            CALL CHKOUT ( 'ZZDAFGSR' )
            RETURN
         END IF
 
C
C        Translate the summary record.  First extract the leading
C        3 double precision numbers from the summary record as these
C        respectively are NEXT, PREV, and NSUM.
C
         CALL ZZXLATED ( IBFF, CHRBUF(1:3*DPLEN), 128, DPBUF )
 
C
C        Check FAILED() in case the translation process fails for
C        any reason.
C
         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZDAFGSR' )
            RETURN
         END IF
 
C
C        Convert NSUM to an integer, and compute the number of
C        double precision numbers required to store each individual
C        summary in the record.
C
         NSUM   = INT ( DPBUF(3) )
         SUMSIZ = ( ND + (NI+1)/2 )
 
C
C        Convert each of the summaries one at a time.
C
         DO I = 1, NSUM
 
C
C           Set the start index into the double precision array
C           to receive the componets.  Also set the character
C           substring index to the start location for this summary.
C           In the diagram below, each box represents a double
C           precision number.  The figure assumes SUMSIZ is 5
C           double precision numbers:
C
C                 |--- 1 ---|--- 2 ---|--- 3 ---|   |- (I-1) -|
C           -------------------------------------   -------------
C           | | | | | | | | | | | | | | | | | | |...| | | | | | |...
C           -------------------------------------   -------------
C           |-----|                                            ^
C              ^                                               |
C              |                                            Summary
C           NEXT, PREV, NSUM                                 Start
C
            DINDEX = (I-1)*SUMSIZ + 4
            CINDEX = DPLEN*(DINDEX-1) + 1
 
C
C           First, check to see if there are any double precision
C           numbers to translate.  If so, translate, and then
C           increment DINDEX and CINDEX accordingly.
C
            IF ( ND .GT. 0 ) THEN
 
C
C              DPBUF has room for 128 double precision numbers
C              total.  Compute the amount of space left in the
C              buffer.
C
               LEFT = 128 - (I-1)*SUMSIZ - 3
 
               CALL ZZXLATED ( IBFF,
     .                         CHRBUF(CINDEX:CINDEX+ND*DPLEN-1),
     .                         LEFT,
     .                         DPBUF(DINDEX)                     )
 
C
C              If the translation routine fails for any reason,
C              check out and return.
C
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDAFGSR' )
                  RETURN
               END IF
 
               DINDEX = DINDEX + ND
               CINDEX = CINDEX + DPLEN*ND
 
            END IF
 
C
C           At this point DINDEX and CINDEX are pointing at the
C           locations for the packed integers in the record.
C           Use DINDEX to compute the index into INBUF, the
C           equivalenced integer buffer and translate.
C
            IF ( NI .GT. 0 ) THEN
 
C
C              INBUF has room for 256 integers total.  Compute
C              the amount of space left in the buffer.  Since
C              it is equivalenced to DPBUF, account for the
C              double precision numbers that were just added.
C
               LEFT = 256 - 2*(I-1)*SUMSIZ - ND*2 - 6
 
               CALL ZZXLATEI ( IBFF,
     .                         CHRBUF(CINDEX:CINDEX+NI*INLEN-1),
     .                         LEFT,
     .                         INBUF(2*DINDEX-1)                 )
 
C
C              If the translation routine fails for any reason,
C              check out and return.
C
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDAFGSR' )
                  RETURN
               END IF
 
C
C              Now check to see if NI is odd.  If so, then zero
C              the last integer occupied by the newly translated
C              summary.  This is necessary to purge any garbage
C              present in memory.
C
               IF ( MOD(NI,2) .EQ. 1 ) THEN
                  INBUF(2*DINDEX-1+NI) = 0
               END IF
 
            END IF
 
         END DO
 
C
C        Translating garbage is a bad idea in general, so set
C        the any remaining double precision numbers in the summary
C        record to 0.
C
         DINDEX = NSUM*SUMSIZ + 4
 
         DO I = DINDEX, 128
            DPBUF(I) = 0
         END DO
 
      END IF
 
C
C     Transfer the DPs to the output argument and return to the
C     caller.
C
      FOUND = .TRUE.
      CALL MOVED ( DPBUF, 128, DPREC )
 
      CALL CHKOUT( 'ZZDAFGSR' )
      RETURN
 
      END
