C$Procedure ZZDASGRD ( DAS, get record, double precision )
 
      SUBROUTINE ZZDASGRD ( HANDLE, RECNO, RECORD )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C 
C     Read DAS double precision physical records from native
C     or non-native DAS files.
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
C
C$ Keywords
C
C     DAS
C     FILES
C
C$ Declarations

      IMPLICIT NONE
      
      INCLUDE 'das.inc'
      
      INTEGER               HANDLE
      INTEGER               RECNO
      DOUBLE PRECISION      RECORD ( NWD )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DAS file.
C     RECNO      I   Record number.
C     RECORD     O   Double precision data read from record.
C     NWD        P   Number of DP numbers in a single DAS DP record.
C
C$ Detailed_Input
C
C     HANDLE         is the handle of an open DAS file.
C
C     RECNO          is the number of a record in a DAS file.
C
C$ Detailed_Output
C
C     RECORD         is a double precision array containing the
C                    elements of the specified record.
C
C$ Parameters
C
C     NWD            is the number of DPs in a single DAS record
C                    containing DPs.
C 
C$ Exceptions
C
C     1)  If the input file handle cannot be mapped to a file
C         name, the error SPICE(HANDLENOTFOUND) will be signaled. The
C         output argument RECORD will not be modified.
C
C     2)  If a read operation attempted by this routine fails, the
C         error SPICE(DASFILEREADFAILED) will be signaled.
C
C     3)  If an error occurs while attempting to translate non-native
C         double precision data to native format, the error will be
C         diagnosed by a routine in the call tree of this routine.
C
C$ Files
C
C     See the description of the argument HANDLE in $Detailed_Input.
C
C$ Particulars
C
C     Routines outside of SPICELIB will normally have no need to call
C     this routine.
C
C     This routine enables DAS routines to read double precision data
C     records from native or non-native DAS files.
C
C     This routine should be used to read only records that contain
C     double precision data.
C
C$ Examples
C
C     See usage in DASRRD.
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
C-    SPICELIB Version 1.0.0, 10-FEB-2017 (NJB)
C
C-&
 
C$ Index_Entries
C
C     read DAS d.p. physical records for arbitrary BFF
C     read DAS d.p. records for native or non-native DAS
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               FILSIZ
      PARAMETER           ( FILSIZ = 255 )

C
C     Local variables
C
      CHARACTER*(NWC)       CHRREC
      CHARACTER*(FILSIZ)    FNAME

      INTEGER               INTAMH
      INTEGER               INTARC
      INTEGER               INTBFF
      INTEGER               IOSTAT
      INTEGER               NATBFF
      INTEGER               UNIT

      LOGICAL               FIRST
      LOGICAL               FOUND

C
C     Saved variables
C
      SAVE                 FIRST
      SAVE                 INTBFF
      SAVE                 NATBFF

C     
C     Initial values
C
      DATA                 FIRST  / .TRUE. /
      

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDASGRD' )


      IF ( FIRST ) THEN
C
C        Get the integer code for the host binary file format.
C
         CALL ZZDDHNFC ( NATBFF )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDASGRD' )
            RETURN
         END IF

         FIRST = .FALSE.

      END IF

C
C     Get a logical unit for this file.
C
      CALL ZZDDHHLU ( HANDLE, 'DAS', .FALSE., UNIT )

C
C     Get the binary file format of the file designated by HANDLE.
C 
      CALL ZZDDHNFO ( HANDLE, FNAME, INTARC, INTBFF, INTAMH, FOUND )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZDASGRD' )
         RETURN
      END IF

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'Unable to locate file associated with '
     .   //            'HANDLE, #. The most likely cause of this '
     .   //            'is the file that you are trying to read '
     .   //            'has been closed.'                           )
         CALL ERRINT ( '#', HANDLE                                  )
         CALL SIGERR ( 'SPICE(HANDLENOTFOUND)'                      )
         CALL CHKOUT ( 'ZZDASGRD'                                   )
         RETURN

      END IF


      IF ( INTBFF .EQ. NATBFF ) THEN
C
C        The file has native format. 
C
         READ (  UNIT    =   UNIT,
     .           REC     =   RECNO,
     .           IOSTAT  =   IOSTAT ) RECORD
 
         IF ( IOSTAT .NE. 0 ) THEN
 
            CALL SETMSG ( 'Could not read DAS d.p. record. '       
     .      //            'File = # Record number = #. IOSTAT = #.' )
            CALL ERRFNM ( '#', UNIT                                 )
            CALL ERRINT ( '#', RECNO                                )
            CALL ERRINT ( '#', IOSTAT                               )
            CALL SIGERR ( 'SPICE(DASFILEREADFAILED)'                )
            CALL CHKOUT ( 'ZZDASGRD'                                )
            RETURN
 
         END IF


      ELSE
C
C        Read the record as a character string, then translate it
C        to an array of d.p. numbers.
C
         READ (  UNIT    =    UNIT,
     .           REC     =    RECNO,
     .           IOSTAT  =    IOSTAT  ) CHRREC

         IF ( IOSTAT .NE. 0 ) THEN
 
            CALL SETMSG ( 'Could not read non-native DAS d.p. ' 
     .      //            'record into character array. '
     .      //            'File = # Record number = #. IOSTAT = #.'    )
            CALL ERRFNM ( '#', UNIT                                    )
            CALL ERRINT ( '#', RECNO                                   )
            CALL ERRINT ( '#', IOSTAT                                  )
            CALL SIGERR ( 'SPICE(DASFILEREADFAILED)'                   )
            CALL CHKOUT ( 'ZZDASGRD'                                   )
            RETURN
 
         END IF

C
C        Translate the character record to double precision type.
C
         CALL ZZXLATED ( INTBFF, CHRREC, NWD, RECORD )

C
C        We don't test FAILED here because the routine
C        will return from this point.
C
      END IF

      CALL CHKOUT ( 'ZZDASGRD' )
      RETURN
      END 
