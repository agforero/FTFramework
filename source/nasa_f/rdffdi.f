C$Procedure RDFFDI ( read and parse flatfile records )
 
      SUBROUTINE RDFFDI ( INFILE, NREC, FORMAT, 
     .                    ND, NI, NC, ARD, ARI, ARC, EOF)
      IMPLICIT NONE
 
C$ Abstract
C
C     Read a flatfile record and parse it into integer and/or double
C     precision and/or string values according to a given format 
C     specification.
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
C     FLATFILE
C
C$ Declarations

      CHARACTER*(*)         INFILE
      INTEGER               NREC
      CHARACTER*(*)         FORMAT
      INTEGER               ND
      INTEGER               NI 
      INTEGER               NC
      DOUBLE PRECISION      ARD     ( * )
      INTEGER               ARI     ( * )
      CHARACTER*(*)         ARC     ( * )
      LOGICAL               EOF 

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of input flatfile.
C     NREC       I   Current record number in flatfile.
C     FORMAT     I   String containing format picture.
C     ND         O   Number of doubles found in data record.
C     NI         O   Number of integers found in data record.
C     NC         O   Number of strings found in data record.
C     ARD        O   Array of doubles read from data record.
C     ARI        O   Array of integer read from data record.
C     ARI        O   Array of strings read from data record.
C     EOF        O   End of file flag.
C
C$ Detailed_Input
C
C     INFILE     is the full pathname of the input flatfile.
C
C     NREC       is the current record number in the flatfile,
C                provided by the user.
C
C     FORMAT     is a string containing a format description or
C                picture of the data to be parsed.  In the present
C                scheme, 'I' signifies integer and 'D' double
C                precision values.  A typical format string might
C                be 'I D D D'.  Up to MAXITM format items are 
C                permitted.
C
C$ Detailed_Output
C
C     ND         is the number of double precision values found.
C
C     NI         is the number of integer values found.
C
C     NC         the number of strings extracted from the data record.
C
C     ARD        is an array containing ND consecutively read
C                double precision values.
C
C     ARI        is an array containing NI consecutively read
C                integer values.
C
C     ARC        is an array containing NC consecutively read
C                strings.
C
C     EOF        is .TRUE. if end of file was encountered.
C
C$ Parameters
C
C 
C     MAXWD      Maximum allowed length of a single field.
C  
C     MAXSTR     Maximum length of an input data record.
C 
C     MAXITM     Maximum number of items in format and input 
C                data records.
C
C$ Exceptions
C
C     1)  If the format string length exceeds MAXSTR, then
C         SPICE(FORMATSTRINGTOOLONG) is signaled.
C
C     2)  If the number of format items exceeds MAXITM, then
C         SPICE(FORMATITEMLIMITEXCEEDED) is signaled.
C
C     3)  If the number of data items exceeds MAXITM, then
C         SPICE(DATAITEMLIMITEXCEEDED) is signaled.
C
C     4)  If an unknown format type occurs in the format string,
C         then SPICE(BADFORMATSPECIFIER) is signaled.
C
C     5)  If a string representing numerical data is exactly
C         MAXWD in length then SPICE(DATAWIDTHERROR) is signaled.
C
C     6)  If the format string and data string contain a
C         different number of items delimited by spaces,
C         then SPICE(FORMATDATAMISMATCH) is signaled.
C
C     7)  If an error is returned from NPARSI, then SPICE(BADINTEGER)
C         is signaled.
C
C     8)  If an error is returned from NPARSD, then
C         SPICE(BADDOUBLEPRECISION) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     1)  It is the responsibility of the programmer to supply the
C         correct record number in the file, NREC.
C
C     2)  Format and data items are delimited by spaces.
C
C     3)  ARC must be declared of suffcient length in the calling
C         program to hold the string data items, otherwise the
C         string is truncated when copied to ARC.
C
C$ Examples
C
C     C
C     C     Read a "self-describing" flatfile:
C     C
C     C     The first record contains a single integer,
C     C     which is the number of data records to be read.
C     C
C     C     The following line contains a format string for the
C     C     data records to follow.
C     C
C     C     After the last data record is read, a new data group
C     C     may begin, signified by a line containing an
C     C     integer count, followed by a format string and
C     C     data records.  Alternatively, the end of the file may
C     C     also be encountered.
C     C
C
C
C           LOGICAL               EOF
C           CHARACTER*(*)         FNAME
C           CHARACTER*(MAXSTR)    FORMAT
C           INTEGER               ND
C           INTEGER               NI
C           INTEGER               NC
C           INTEGER               NREC
C           INTEGER               NDAT
C           INTEGER               ARI     ( MAXITM )
C           DOUBLE PRECISION      ARD     ( MAXITM )
C           CHARACTER*(MAXSTR)    ARC     ( MAXITM )
C
C
C           CALL PROMPT ( 'Enter file name ', FNAME )
C
C     C
C     C     Get record count.
C     C
C
C           NREC = 1
C           CALL RDFFDI ( FNAME, NREC, 'I', ND, NI, NC,
C          .              ARD, ARI, ARC, EOF )
C
C           DO WHILE ( .NOT. EOF )
C
C              NDAT = ARI(1)
C
C     C
C     C        Get data record format.
C     C
C
C              NREC = NREC + 1
C              CALL RDTEXT ( FNAME, FORMAT, EOF )
C
C              IF ( EOF ) THEN
C                 CALL EXTMSI ( 'End of file after line #.', '#',
C           .                    NREC-1 )
C              END IF
C
C     C
C     C        Read data records.
C     C
C
C              DO I = 1, NDAT
C                 NREC = NREC + 1
C                 CALL RDFFDI ( FNAME, NREC, FORMAT, ND, NI, NC
C           .                   ARD,   ARI,  ARC, EOF )
C                 IF ( EOF ) THEN
C                    CALL EXTMSI ( 'End of file after line #.', '#',
C           .                       NREC-1 )
C                 END IF
C                 CALL HAVFUN ( ND, NI, ARD, ARI )
C              END DO
C
C     C
C     C        Get next data record count, if it's there.
C     C
C
C              NREC = NREC + 1
C              CALL RDFFDI ( FNAME, NREC, 'I', ND, NI, NC,
C           .                ARD, ARI, ARC, EOF )
C
C           END DO
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
C     J.A. Bytof      (JPL)
C
C$ Version
C
C-    SPICELIB Version 3.0.0, 03-MAY-2014 (NJB)
C
C        Increased supported input file line length to 255.
C        Added SAVE statements for several arrays and strings.
C
C-    SPICELIB Version 2.0.1, 08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C-    SPICELIB Version 2.0.0, 25-OCT-2004 (EDW)
C
C        Added capability to process character
C        typed data, 'C'.
C
C-    SPICELIB Version 1.0.0, 09-APRIL-1997 (JAB)
C
C-&
 
C$ Index_Entries
C
C     read and parse flatfile records.
C
C-&
 
C
C     SPICE functions.
C
      INTEGER               WDCNT
 
C
C     Local variables.
C

C 
C     Maximum allowed length of a single field.
C 
      INTEGER               MAXWD
      PARAMETER           ( MAXWD  = 40 )

C 
C     Maximum length of an input data record.
C 
      INTEGER               MAXSTR
      PARAMETER           ( MAXSTR = 255 )

C 
C     Maximum number of items in format and input data records.
C 
      INTEGER               MAXITM
      PARAMETER           ( MAXITM = 6 )
 
      INTEGER               I
      INTEGER               PTR
      INTEGER               NFMT
      INTEGER               NDAT
      INTEGER               DATI
      INTEGER               LASTNB
 
      LOGICAL               RETURN
 
      CHARACTER*(160)       ERROR
      CHARACTER*(MAXSTR)    LINE
      CHARACTER*(MAXSTR)    FMTSAV
      CHARACTER*(MAXSTR)    FMTITM     ( MAXITM )
      CHARACTER*(MAXSTR)    DATITM     ( MAXITM )
 
      DOUBLE PRECISION      DATD
 
      SAVE                  ERROR
      SAVE                  DATITM
      SAVE                  LINE
      SAVE                  NFMT
      SAVE                  FMTITM
      SAVE                  FMTSAV
 
      DATA                  FMTSAV     / 'X' /
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'RDFFDI' )
      END IF

 
C
C     Initialize array counters and logical flags.
C
      ND  = 0
      NI  = 0
      NC  = 0
      EOF = .FALSE.
 
C
C     Check format string length.
C
 
      IF ( LASTNB(FORMAT) .GT. MAXSTR ) THEN
         CALL SETMSG (     'Format string length exceeds limit'
     .                 //  ' in line #.'                          )
         CALL ERRINT ( '#', NREC                                  )
         CALL SIGERR ( 'SPICE(FORMATSTRINGTOOLONG)'               )
      END IF
 
C
C     Parse format string if it has changed.  Verify that format.
C
      IF ( FORMAT .NE. FMTSAV ) THEN
 
C
C     Does the number of format items exceed MAXITM?
C
         IF ( WDCNT( FORMAT ) .GT. MAXITM ) THEN
            CALL SETMSG ( 'Too many format items at line #'      )
            CALL ERRINT ( '#', NREC                              )
            CALL SIGERR ( 'SPICE(FMTITEMLIMITEXCEEDED)'          )
         END IF
 
         CALL LPARSE ( FORMAT,  ' ', MAXITM, NFMT, FMTITM )
 
         DO I = 1, NFMT
         
            IF ( LASTNB(FMTITM(I)) .GT. 1 ) THEN
               CALL SETMSG ( 'Format error detected at line #.'  )
               CALL ERRINT ( '#', NREC                           )
               CALL SIGERR ( 'SPICE(BADFORMATSPECIFIER)'         )
            END IF
         
         END DO
         
         FMTSAV = FORMAT
 
      END IF
 
 
C
C     Read data record, if end of file is encountered,
C     check out, then return to calling routine.
C
      CALL RDTEXT ( INFILE, LINE, EOF )
 
      IF ( EOF ) THEN
         CALL CHKOUT ( 'RDFFDI' )
         RETURN
      END IF
 
C
C     Convert any tabs or carriage returns to blanks.
C
      CALL  REPLCH ( LINE, CHAR( 9), ' ', LINE )
      CALL  REPLCH ( LINE, CHAR(13), ' ', LINE )
 
C
C     Does the number of data items exceed MAXITM?
C
      NDAT = WDCNT ( LINE )
 
      IF ( NDAT .GT. MAXITM ) THEN
         CALL SETMSG ( 'Too many data items, line #'  )
         CALL ERRINT ( '#', NREC                      )
         CALL SIGERR ( 'SPICE(DATAITEMLIMITEXCEEDED)' )
      END IF

C
C     The format specifier and data record might mistakenly
C     have a different number of fields.
C
      IF (      NFMT .LT. NDAT ) THEN
 
         CALL SETMSG ( 'Too many data items in line #.'
     .               //' The PLATE_TYPE setting may not '
     .               //' match the data file format.'    )
         CALL ERRINT ( '#', NREC                         )
         CALL SIGERR ( 'SPICE(FORMATDATAMISMATCH)'       )
 
      ELSE IF ( NFMT .GT. NDAT ) THEN
 
         CALL SETMSG ( 'Too many format items in line #' 
     .               //' The PLATE_TYPE setting may not '
     .               //' match the data file format.'    )
         CALL ERRINT ( '#', NREC                         )
         CALL SIGERR ( 'SPICE(FORMATDATAMISMATCH)'       )
 
      END IF
 
C
C     Parse data, verify that values fit within MAXWD field width.
C
      CALL LPARSE ( LINE,  ' ', MAXITM, NDAT, DATITM )
      
      DO I = 1, NDAT
      
         IF ( LASTNB(DATITM(I)) .GE. MAXWD ) THEN
            CALL SETMSG ( 'Possible data error in line #.' )
            CALL ERRINT ( '#', NREC                        )
            CALL SIGERR ( 'SPICE(DATAWIDTHERROR)'          )
         END IF
      
      END DO
 
C
C     Process each data field according to its format.
C
      DO I = 1, NDAT

C
C     Convert the data field to requested type.
C
         IF ( FMTITM(I)(1:1) .EQ. 'I' ) THEN

C
C           We expect an integer. Use NPARSI.
C
            CALL NPARSI ( DATITM(I), DATI, ERROR, PTR )
            IF ( PTR .NE. 0 ) THEN
               CALL SETMSG ( 'Integer error (#) in line #.' )
               CALL ERRCH  ( '#', ERROR                     )
               CALL ERRINT ( '#', NREC                      )
               CALL SIGERR ( 'SPICE(BADINTEGER)'            )
            END IF
 
            NI         = NI + 1
            ARI ( NI ) = DATI
 
         ELSE IF (FMTITM(I)(1:1) .EQ. 'D' ) THEN
 
C
C           We expect a double. Use NPARSD.
C
            CALL NPARSD ( DATITM(I), DATD, ERROR, PTR )
            IF ( PTR .NE. 0 ) THEN
               CALL SETMSG ( 'D.P. error (#) in line #.' )
               CALL ERRCH  ( '#', ERROR                  )
               CALL ERRINT ( '#', NREC                   )
               CALL SIGERR ( 'SPICE(BADDOUBLEPRECISION)' )
            END IF
 
            ND         = ND + 1
            ARD ( ND ) = DATD

         ELSE IF (FMTITM(I)(1:1) .EQ. 'C' ) THEN

C
C           We expect a string. No need to parse, just copy
C           then count. Notice that if ARC was not declared
C           with a string size suffcient to store the
C           DATITM value, the copy op will truncate the value.
C
            NC         = NC + 1
            ARC ( NC ) = DATITM(I)

         ELSE

C
C           Problem with format specifier.
C 
            CALL SETMSG ( 'Bad format specifier # used in line #.' )
            CALL ERRCH  ( '#', FMTITM(I)                           )
            CALL ERRINT ( '#', NREC                                )
            CALL SIGERR ( 'SPICE(BADFORMATSPECIFIER)'              )
 
         END IF
 
      END DO
 
      CALL CHKOUT ( 'RDFFDI' )
 
      RETURN
      END

C12345678901234567890123456789012345678901234567890123456789012345678912
