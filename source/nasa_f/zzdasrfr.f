C$Procedure      ZZDASRFR ( DAS, read file record )
 
      SUBROUTINE ZZDASRFR ( HANDLE,
     .                      IDWORD,
     .                      IFNAME,
     .                      NRESVR,
     .                      NRESVC,
     .                      NCOMR,
     .                      NCOMC  )
 
C$ Abstract
C
C     Return the contents of the file record of a specified DAS file.
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
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'zzddhman.inc'

      INTEGER               HANDLE
      CHARACTER*(*)         IDWORD
      CHARACTER*(*)         IFNAME
      INTEGER               NRESVR
      INTEGER               NRESVC
      INTEGER               NCOMR
      INTEGER               NCOMC
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DAS file handle.
C     IDWORD     O   ID word.
C     IFNAME     O   DAS internal file name.
C     NRESVR     O   Number of reserved records in file.
C     NRESVC     O   Number of characters in use in reserved rec. area.
C     NCOMR      O   Number of comment records in file.
C     NCOMC      O   Number of characters in use in comment area.
C
C$ Detailed_Input
C
C     HANDLE         is a file handle for a previously opened DAS file.
C
C$ Detailed_Output
C
C     IDWORD      is the `ID word' contained in the first eight
C                 characters of the file record.
C
C     IFNAME      is the internal file name of the DAS file.  The
C                 maximum length of the internal file name is 60
C                 characters.
C
C     NRESVR      is the number of reserved records in the DAS file
C                 specified by HANDLE.
C
C     NRESVC      is the number of characters in use in the reserved
C                 record area of the DAS file specified by HANDLE.
C
C     NCOMR       is the number of comment records in the DAS file
C                 specified by HANDLE.
C
C     NCOMC       is the number of characters in use in the comment area
C                 of the DAS file specified by HANDLE.
C
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the file read attempted by this routine fails, the error
C        SPICE(DASFILEREADFAILED) will be signaled.
C
C     2) If the input file handle is invalid, the error will diagnosed
C        by a routine in the call tree of this routine.
C
C     3) If a logical unit cannot be obtained for the file designated
C        by HANDLE, the error will diagnosed by a routine in the call
C        tree of this routine.
C
C     4) If the file's binary format is unrecognized, the error will
C        diagnosed by a routine in the call tree of this routine.
C
C     5) If the file designated by HANDLE has non-native binary format,
C        and if any numeric components of the file record cannot be
C        translated to native format, the error will diagnosed
C        by a routine in the call tree of this routine.
C
C$ Files
C
C     See the description of HANDLE under $Detailed_Input.
C
C$ Particulars
C
C     This routine provides a convenient way of retrieving the
C     information contained in the file record of a DAS file.
C
C$ Examples
C
C     1)  Obtain the internal file name of an existing DAS file.
C
C
C            C
C            C     Open the file for reading.
C            C
C                  CALL DASOPR ( FNAME, HANDLE  )
C
C            C
C            C     Retrieve the internal file name and print it.
C            C
C
C                  CALL ZZDASRFR ( HANDLE,
C                 .                IDWORD,
C                 .                IFNAME,
C                 .                NRESVR,
C                 .                NRESVC,
C                 .                NCOMR,
C                 .                NCOMC  )
C
C
C                  WRITE (*,*) 'Internal file name is: '
C                  WRITE (*,*)  IFNAME
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
C     K.R. Gehringer (JPL)
C     N.J. Bachman   (JPL)
C     W.L. Taber     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 05-FEB-2015 (NJB)
C
C        Based on DASRFR Version 1.0.0, 15-JUL-1992 (NJB) (WLT)
C
C-&
 
 
C$ Index_Entries
C
C     read DAS file record
C     read DAS internal file name
C
C-&
 
C$ Revisions
C
C     None.
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
      INTEGER               IDWLEN
      PARAMETER           ( IDWLEN =  8 )
 
      INTEGER               IFNLEN
      PARAMETER           ( IFNLEN = 60 )

      INTEGER               NWC
      PARAMETER           ( NWC    = 1024 )

C
C     Parameters for positions of file record elements:
C

C
C     ID word begin and end:
C
      INTEGER               BEGIDW
      PARAMETER           ( BEGIDW = 1 )

      INTEGER               ENDIDW
      PARAMETER           ( ENDIDW = 8 )
 

C
C     Internal file name begin and end:
C
      INTEGER               BEGIFN
      PARAMETER           ( BEGIFN = 9 )

      INTEGER               ENDIFN
      PARAMETER           ( ENDIFN = 68 )

C
C     Reserved record count begin and end:
C
      INTEGER               BEGNRR
      PARAMETER           ( BEGNRR = 69 )

      INTEGER               ENDNRR
      PARAMETER           ( ENDNRR = 72 )

C
C     Reserved record character count begin and end:
C
      INTEGER               BEGNRC
      PARAMETER           ( BEGNRC = 73 )

      INTEGER               ENDNRC
      PARAMETER           ( ENDNRC = 76 )


C
C     Comment area record count begin and end:
C
      INTEGER               BEGNCR
      PARAMETER           ( BEGNCR = 77 )

      INTEGER               ENDNCR
      PARAMETER           ( ENDNCR = 80 )

C
C     Comment area character count begin and end:
C
      INTEGER               BEGNCC
      PARAMETER           ( BEGNCC = 81 )

      INTEGER               ENDNCC
      PARAMETER           ( ENDNCC = 84 )

C
C     Local variables
C
      CHARACTER*(NWC)       CHRBUF
      CHARACTER*(IDWLEN)    TMPIDW
      CHARACTER*(IFNLEN)    TMPIFN
 
      INTEGER               IBFF
      INTEGER               IOSTAT
      INTEGER               NATBFF
      INTEGER               UNIT
 
      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  NATBFF

C
C     Initial values
C
      DATA                  FIRST  / .TRUE. /
      DATA                  NATBFF / -1     /
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDASRFR' )

C
C     On the first pass through this routine, get the integer code of
C     the host system's native binary file format.
C
      IF ( FIRST ) THEN

         CALL ZZDDHNFC ( NATBFF )

         IF ( FAILED () ) THEN
            CALL CHKOUT ( 'ZZDASRFR' )
            RETURN
         END IF

         FIRST = .FALSE.

      END IF
 
C
C     Get a logical unit for this DAS file.
C
      CALL ZZDDHHLU ( HANDLE, 'DAS', .FALSE., UNIT )

C
C     Get the integer code for the file's binary format.
C
      CALL ZZDDHPPF ( UNIT, DAS, IBFF )
 
      IF ( FAILED () ) THEN
         CALL CHKOUT ( 'ZZDASRFR' )
         RETURN
      END IF

      IF ( IBFF .EQ. NATBFF ) THEN
C
C        We're looking at a native file. Just read the file record.
C
         READ ( UNIT,
     .          REC     =  1,
     .          IOSTAT  =  IOSTAT ) TMPIDW,
     .                              TMPIFN,
     .                              NRESVR,
     .                              NRESVC,
     .                              NCOMR,
     .                              NCOMC
 
         IF ( IOSTAT .NE. 0 ) THEN
 
            CALL SETMSG ( 'Could not DAS read file record. File was ' 
     .      //            '#.  IOSTAT was #.'                        )
            CALL ERRFNM ( '#', UNIT                                  )
            CALL ERRINT ( '#', IOSTAT                                )
            CALL SIGERR ( 'SPICE(DASFILEREADFAILED)'                 )
            CALL CHKOUT ( 'ZZDASRFR'                                 )
            RETURN
 
         END IF
 
         IDWORD  =  TMPIDW
         IFNAME  =  TMPIFN
 

      ELSE 
C
C        The file is non-native.
C
C        We don't check the access mode of the file, because we're
C        not going to reject files that are open for writing. 
C
C        We'll read the file record as a character string and then
C        pick it apart.
C
         READ ( UNIT,
     .          REC     =  1,
     .          IOSTAT  =  IOSTAT ) CHRBUF
 
         IF ( IOSTAT .NE. 0 ) THEN
 
            CALL SETMSG ( 'Could not read DAS file record. File is ' 
     .      //            '#. IOSTAT was #. File''s BFF integer '
     .      //            'code is #.'                              )
            CALL ERRFNM ( '#', UNIT                                 )
            CALL ERRINT ( '#', IOSTAT                               )
            CALL ERRINT ( '#', IBFF                                 )
            CALL SIGERR ( 'SPICE(DASFILEREADFAILED)'                )
            CALL CHKOUT ( 'ZZDASRFR'                                )
            RETURN
 
         END IF

C
C        Set the string output arguments.
C
         IDWORD = CHRBUF( BEGIDW : ENDIDW )
         IFNAME = CHRBUF( BEGIFN : ENDIFN )

C
C        The integer output arguments require translation.
C
         CALL ZZXLATEI ( IBFF,  CHRBUF( BEGNRR : ENDNRR ),  1,  NRESVR )
         CALL ZZXLATEI ( IBFF,  CHRBUF( BEGNRC : ENDNRC ),  1,  NRESVC )
         CALL ZZXLATEI ( IBFF,  CHRBUF( BEGNCR : ENDNCR ),  1,  NCOMR  )
         CALL ZZXLATEI ( IBFF,  CHRBUF( BEGNCC : ENDNCC ),  1,  NCOMC  )

      END IF

      CALL CHKOUT ( 'ZZDASRFR' )
      RETURN
      END
