C$Procedure ZZLEXMET ( Scan method string )
 
      SUBROUTINE ZZLEXMET ( METHOD, MAXN, N, BEGS, ENDS )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Scan a method input string passed by a caller to a SPICE 
C     high-level geometry API. Return tokens found in the string.
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
C     DSK
C
C$ Keywords
C
C     SCANNING
C     PARSING
C     CONSTANTS
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE

      CHARACTER*(*)         METHOD
      INTEGER               MAXN
      INTEGER               N
      INTEGER               BEGS   ( * )
      INTEGER               ENDS   ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     METHOD     I   Method string.
C     MAXN       I   Maximum number of tokens that can be returned.
C     N          O   Number of tokens found.
C     BEGS       O   Beginning indices of tokens.
C     ENDS       O   Ending indices of tokens.
C     
C$ Detailed_Input
C
C     METHOD     is a "method" string to be parsed. This string
C                normally is passed by a caller to a SPICE high-level 
C                high-level geometry routine.
C
C                METHOD indicates whether a target body surface is
C                to be modeled as an ellipsoid or by DSK data. It
C                may contain additional, computation-dependent
C                parameters.
C
C     MAXN       is the maximum number of tokens that may be returned.
C                The output arrays BEGS and ENDS must be declared with
C                size at least MAXN.  It's an error to supply output
C                arrays that are too small to hold all of the tokens in
C                the input method string.
C 
C$ Detailed_Output
C
C     N          is the number of tokens found in the method string.
C
C     BEGS,
C     ENDS       are, respectively, the indices of the first and last
C                characters of each token in the input string METHOD.
C                BEGS(I) and ENDS(I) indicate the index range of the 
C                Ith token.
C                   
C                See the Particulars section for a list of token
C                delimiters.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error is found while parsing a quoted string,
C         the error SPICE(SYNTAXERROR) is signaled.
C
C     2)  If the number of tokens found in the method string 
C         is greater than the input limit MAXN, the error
C         SPICE(ARRAYTOOSMALL) is signaled.
C          
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports parsing of the METHOD input argument
C     used by SPICE geometry APIs such as, but not limited to:
C
C        GFOCLT
C        ILLUMF
C        ILLUMG
C        ILUMIN
C        LATSRF
C        LIMBPT
C        SINCPT
C        SUBPNT
C        SUBSLR
C        TERMPT
C
C     Token delimiters are 
C
C        / = , <blank>
C
C     Strings delimited by double quote characters are treated as
C     individual tokens. There is no escape syntax for indicating a
C     doubly quoted string within a doubly quoted string.
C
C$ Examples
C
C 
C     Example:
C
C       Method string =
C
C          METHOD = 'INTERCEPT / '
C          //       'DSK/UNPRIORITIZED/SURFACES = "MGS '
C          //       'MOLA 128 PIXEL/DEG", MARS_LOWRES'
C
C       Token list = 
C 
C          INTERCEPT
C          /
C          DSK
C          /
C          UNPRIORITIZED
C          /
C          SURFACES
C          =
C          "MGS MOLA 128 PIXEL/DEG"
C          ,
C          MARS_LOWRES
C
C
C$ Restrictions
C
C     This is a specialized scanning routine. It it meant to
C     be used only for the parsing "method" strings as 
C     described in Particulars.
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
C-    SPICELIB Version 1.0.0, 02-FEB-2016 (NJB) 
C
C        Based on first version 13-JAN-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     scan method string for geometry api routines
C     lex method string for geometry api routines
C     extract tokens from method string for geometry apis
C-&


C
C     SPICELIB functions
C
      INTEGER               CPOS
      INTEGER               LTRIM
      INTEGER               RTRIM

      LOGICAL               RETURN
      
C
C     Local parameters
C
C
C     Note the leading blank in DELIM below; this
C     indicates that the blank character is a delimiter.
C
      CHARACTER*(*)         DELIM
      PARAMETER           ( DELIM  = ' /,=:' )

      CHARACTER*(*)         DQUOTE
      PARAMETER           ( DQUOTE  = '"' )

C
C     Local variables
C
      INTEGER               B
      INTEGER               E
      INTEGER               EQ
      INTEGER               I
      INTEGER               LOC
      INTEGER               NCHAR
      INTEGER               QPOS
      INTEGER               R
      INTEGER               ROOM
      INTEGER               TKE

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZLEXMET' )


      IF ( METHOD .EQ. ' ' ) THEN
C
C        No tokens here.
C         
         N = 0

         CALL CHKOUT ( 'ZZLEXMET' )
         RETURN

      END IF


      N    = 0      
      ROOM = MAXN
      R    = RTRIM( METHOD )
      I    = 1

      DO WHILE ( I .LE. R )
C
C        Look ahead in the input string for the start of a 
C        quoted string.
C
         QPOS = CPOS ( METHOD(I:), DQUOTE, 1 ) 
         B    = I

         IF ( QPOS .EQ. 0 ) THEN
C
C           There are no quoted string tokens in the input string
C           from index I onward.
C
            E = R
         ELSE
            E = I + QPOS - 2
         END IF


         IF ( B .LE. E ) THEN
C
C           Locate any tokens between indices B and E.
C         
            I = B

            DO WHILE ( I .LE. E )
C
C              Find the next delimiter in the substring
C              from indices I : E.
C
               LOC = CPOS( METHOD(I:E), DELIM, 1 )

               IF ( LOC .EQ. 1 ) THEN
C
C                 There is a delimiter character at index I in METHOD.
C 
                  TKE  = I - 1 + LOC

               ELSE IF ( LOC .GT. 1 ) THEN
C
C                 There is a delimiter character at index 
C
C                    I - 1 + LOC
C
C                 in METHOD. 
C 
                  IF ( METHOD(I:I+LOC-2) .NE. ' ' ) THEN
C
C                    There's a token preceding the delimiter.
C
                     TKE  = I - 2 + LOC

                  ELSE
C
C                    The delimiter is all we've got.
C
                     TKE = I - 1 + LOC

                  END IF

               ELSE
C
C                 The token, if any, ends at the end of
C                 substring we're considering.
C                  
                  TKE  = E

               END IF

C
C              There is a token, which may be blank, between
C              indices I and TKE. We don't return blank tokens
C              in the output array.
C
               IF ( METHOD( I : TKE ) .NE. ' ' ) THEN

                  IF ( ROOM .GT. 0 ) THEN

                     N       = N    + 1
                     ROOM    = ROOM - 1
                     BEGS(N) = LTRIM( METHOD(I : TKE) ) + I - 1
                     ENDS(N) = RTRIM( METHOD(I : TKE) ) + I - 1

                  ELSE

                     CALL SETMSG ( 'Need more room in output '
     .               //            'arrays. Token count = #; '
     .               //            'substring indices = #:#; '
     .               //            'substring = #.'           )
                     CALL ERRINT ( '#', N                     )
                     CALL ERRINT ( '#', I                     )
                     CALL ERRINT ( '#', TKE                   )
                     CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'     )
                     CALL CHKOUT ( 'ZZLEXMET'                 )
                     RETURN

                  END IF

               END IF

               I = TKE + 1

            END DO
            
         END IF


         IF ( E .LT. R ) THEN
C
C           We expect to find at least one quoted string starting
C           at index E + 1.
C
            I  = E + 1            

            CALL LXQSTR ( METHOD(I:), DQUOTE, 1, EQ, NCHAR )

            IF ( NCHAR .GT. 0 ) THEN

               IF ( ROOM .GT. 0 ) THEN

                  N       = N    + 1
                  ROOM    = ROOM - 1
                  BEGS(N) = I   
                  ENDS(N) = I    - 1 + EQ 
                                          
               ELSE

                  CALL SETMSG ( 'Need more room in output '
     .            //            'arrays. Token count = #; '
     .            //            'substring indices = #:#; '
     .            //            'substring = #.'           )
                  CALL ERRINT ( '#', N                     )
                  CALL ERRINT ( '#', I                     )
                  CALL ERRINT ( '#', TKE                   )
                  CALL ERRCH  ( '#', METHOD(I:TKE)         )
                  CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'     )
                  CALL CHKOUT ( 'ZZLEXMET'                 )
                  RETURN

               END IF

            ELSE
C
C              We have a syntax error in the input string.
C
               CALL SETMSG ( 'Invalid quoted string found '
     .         //            'starting at index #. Substring '
     .         //            'is #.'                          )
               CALL ERRINT ( '#', I                           )
               CALL ERRCH  ( '#', METHOD(I:)                  )
               CALL SIGERR ( 'SPICE(SYNTAXERROR)'             )
               CALL CHKOUT ( 'ZZLEXMET'                       )
               RETURN

            END IF

            I = ENDS(N) + 1

         END IF
C
C        The index I has been moved forward by at least one position.
C
      END DO

      CALL CHKOUT ( 'ZZLEXMET' )
      RETURN      
      END
