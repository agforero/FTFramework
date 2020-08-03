C$Procedure RDGRD5 ( MKDSK: read from a type 5 height grid file )

      SUBROUTINE RDGRD5 ( INFILE, NMAX, N, VALUES, DONE )

C$ Abstract
C
C     Read data from a MKDSK type 5 height grid file. 
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      CHARACTER*(*)         INFILE
      INTEGER               NMAX
      INTEGER               N
      DOUBLE PRECISION      VALUES ( * )
      LOGICAL               DONE

 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of input file.
C     NMAX       I   Maximum number of values to return.
C     N          O   Number of values returned.
C     VALUES     O   Data values.
C     DONE       O   Flag indicating EOF was reached.
C
C$ Detailed_Input
C
C     INFILE     is the name of an input data file containing height
C                grid data.
C                  
C     NMAX       is the maximum number of values to place in the 
C                output array VALUES.
C   
C$ Detailed_Output
C
C     N          is the number of values in the array VALUES.
C
C     VALUES     is a set of double precision height values 
C                read from the file INFILE.
C
C     DONE       is a logical flag that is set to .TRUE. if and
C                only if the end of file was reached on the 
C                current read. 
C
C                Once DONE is set to .TRUE., a subsequent call
C                to this routine with the same value of INFILE
C                will cause the file to be read again from the 
C                beginning.
C     
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If a token cannot be parsed as a double precision number,
C         the error SPICE(INVALIDDATA) is signaled.
C
C     2)  If the buffer size MAXN is not at least 1, the error 
C         SPICE(INVALIDSIZE) is signaled.
C
C     3)  Any error that occurs while opening or reading the input
C         file will be diagnosed by a routine in the call tree of 
C         this routine.
C
C$ Files
C
C     The file specified by INFILE must contain only tokens that can be
C     read as double precision values. No non-printing characters can
C     be present in the file.
C
C     Tokens can be delimited by blanks or commas. Tokens must not be
C     split across lines.
C
C     Blank lines are allowed; however, their use is discouraged 
C     because they'll cause line numbers in diagnostic messages to
C     be out of sync with actual line numbers in the file.
C
C     The file must end with a line terminator.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     See usage in the MKDSK routine MKVARR.
C
C$ Restrictions
C
C     1) For use only within program MKDSK.
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
C-    MKDSK Version 1.0.0, 25-FEB-2017 (NJB)
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN
      LOGICAL               FAILED

C
C     Local parameters
C
C
C     Comma and space are the allowed delimiters.
C
      CHARACTER*(*)         DELIMS
      PARAMETER           ( DELIMS = ' ,' )

      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 255 )

      INTEGER               MXTOKN
      PARAMETER           ( MXTOKN = 100 )

      INTEGER               TOKLEN
      PARAMETER           ( TOKLEN = 35 )

      INTEGER               MSGLEN
      PARAMETER           ( MSGLEN = 320 )

C
C     Local variables
C
      CHARACTER*(MSGLEN)    ERRMSG
      CHARACTER*(LNSIZE)    LINE
      CHARACTER*(TOKLEN)    TOKENS ( MXTOKN )

      INTEGER               FROM
      INTEGER               I
      INTEGER               LINENO
      INTEGER               NXFR
      INTEGER               NTK
      INTEGER               PTR
      INTEGER               REMAIN
      INTEGER               ROOM
      INTEGER               TO
      
      LOGICAL               EOF
      LOGICAL               NEWFIL
C
C     Saved values
C
      SAVE                  EOF
      SAVE                  FROM
      SAVE                  LINENO
      SAVE                  NEWFIL
      SAVE                  NTK
      SAVE                  REMAIN
      SAVE                  TOKENS

C
C     Initial values
C
      DATA                  EOF    / .FALSE. / 
      DATA                  NTK    / 0       /
      DATA                  NEWFIL / .TRUE.  /
      DATA                  FROM   / 0       /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'RDGRD5' )


C
C     Check NMAX.
C
      IF ( NMAX .LT. 1 ) THEN

         CALL SETMSG ( 'Buffer size NMAX = #; must be strictly '
     .   //            'positive.'                              )
         CALL ERRINT ( '#', NMAX                                )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'                     )
         CALL CHKOUT ( 'RDGRD5'                                 )
         RETURN

      END IF


      IF ( NEWFIL ) THEN
C
C        Initialize our local state variables.
C
         EOF    = .FALSE.
         LINENO = 0
         FROM   = 0
         NTK    = 0
         NEWFIL = .FALSE.

      END IF

C
C     Transfer as many as NMAX values to the output buffer.
C     Stop when the buffer fills up or when we run out of 
C     data.
C
      N     = 0
      TO    = 0
      DONE  = .FALSE.

      DO WHILE ( ( N .LT. NMAX ) .AND. ( .NOT. DONE ) )
C
C        At this point, there's room in the output buffer, and
C        we haven't seen the end of the input file.
C
C        We may have buffered data from the last line read.
C
         IF ( FROM .EQ. 0 ) THEN
C
C           We don't have any buffered data. Read a new non-blank
C           line from the file.
C
            CALL RDNBL ( INFILE, LINE, EOF )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDGRD5' )
               RETURN
            END IF

            IF ( .NOT. EOF ) THEN
C
C              We have a new line. Parse the values in the line.
C
               LINENO = LINENO + 1
             
               CALL LPARSM ( LINE, DELIMS, MXTOKN, NTK, TOKENS )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'RDGRD5' )
                  RETURN
               END IF

C
C              Let NXFR be the number of values to transfer from
C              the current input line.
C
               ROOM = NMAX - N

               NXFR = MIN ( ROOM, NTK )

               DO I = 1, NXFR

                  TO = TO + 1

                  CALL NPARSD ( TOKENS(I), VALUES(TO), ERRMSG, PTR )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'RDGRD5' )
                     RETURN
                  END IF

                  IF ( ERRMSG .NE. ' ' ) THEN

                     CALL SETMSG ( 'Token number # on non-blank '
     .               //            'line # of data file <#>: #'  )
                     CALL ERRINT ( '#', I                        )
                     CALL ERRINT ( '#', LINENO                   )
                     CALL ERRCH  ( '#', INFILE                   )
                     CALL ERRCH  ( '#', ERRMSG                   )
                     CALL SIGERR ( 'SPICE(INVALIDDATA)'          )
                     CALL CHKOUT ( 'RDGRD5'                      )
                     RETURN

                  END IF


               END DO

               N      = N    + NXFR
               REMAIN = NTK  - NXFR

               IF ( REMAIN .GT. 0 ) THEN
C
C                 We didn't transfer all tokens. Let FROM
C                 be the index of the next token to transfer.
C                 We'll transfer the token on the next call.
C                 
                  FROM = NXFR + 1

                  CALL CHKOUT ( 'RDGRD5' )
                  RETURN

               END IF

            ELSE
C
C              There are no more data to be had. RDNBL will close
C              INFILE.
C
               DONE   = .TRUE.
C
C              Get ready for another file.
C
               NEWFIL = .TRUE.

            END IF


         ELSE
C
C           We have buffered tokens. Transfer as many of these
C           as we can.
C
            ROOM = NMAX - N

            NXFR = MIN ( ROOM, REMAIN )

            DO I = FROM, FROM-1+NXFR

               TO = TO + 1

               CALL NPARSD ( TOKENS(I), VALUES(TO), ERRMSG, PTR )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'RDGRD5' )
                  RETURN
               END IF

               IF ( ERRMSG .NE. ' ' ) THEN

                  CALL SETMSG ( 'Token number # on non-blank '
     .            //            'line # of data file <#>: #'  )
                  CALL ERRINT ( '#', I                        )
                  CALL ERRINT ( '#', LINENO                   )
                  CALL ERRCH  ( '#', INFILE                   )
                  CALL ERRCH  ( '#', ERRMSG                   )
                  CALL SIGERR ( 'SPICE(INVALIDDATA)'          )
                  CALL CHKOUT ( 'RDGRD5'                      )
                  RETURN

               END IF

            END DO

            N      = N      + NXFR
            REMAIN = REMAIN - NXFR

            IF ( REMAIN .GT. 0 ) THEN

               FROM = FROM + NXFR

            ELSE
C
C              Indicate the buffer has been exhausted.
C
               FROM = 0

            END IF
            
         END IF

      END DO
      

      CALL CHKOUT ( 'RDGRD5' )       
      RETURN
      END

