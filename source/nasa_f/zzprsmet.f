C$Procedure ZZPRSMET ( Private: parse method string )
 
      SUBROUTINE ZZPRSMET ( BODYID, METHOD, MXNSRF, SHAPE,  SUBTYP,
     .                      PRI,    NSURF,  SRFLST, PNTDEF, TRMTYP )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Parse a method input string passed by a caller to a SPICE
C     high-level geometry API. Return parameter values specified by the
C     string.
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

      INTEGER               BODYID
      CHARACTER*(*)         METHOD
      INTEGER               MXNSRF
      CHARACTER*(*)         SHAPE
      CHARACTER*(*)         SUBTYP
      LOGICAL               PRI
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      CHARACTER*(*)         PNTDEF
      CHARACTER*(*)         TRMTYP
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   Body ID code.
C     METHOD     I   Method string.
C     MXNSRF     I   Maximum number of surface IDs that can be returned.
C     SHAPE      O   Target shape.
C     SUBTYP     O   Sub-point type.
C     PRI        O   Prioritization flag.
C     NSURF      O   Number of surface IDs in list.
C     SRFLST     O   Surface ID list.
C     PNTDEF     O   Limb or terminator point definition.
C     TRMTYP     O   Terminator type.
C     MAXSRF     P   DSK subsystem maximum surface list size.
C
C$ Detailed_Input
C
C     BODYID     is the body ID code of an ephemeris object. This
C                object is normally the "target" body of a geometric
C                computation. Any surfaces specified in the method
C                string are associated with this body. BODYID is
C                needed to map surface names to surface ID codes.
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
C     MXNSRF     is the maximum number of surface IDs that can be 
C                returned in the output array SRFLST. Normally 
C                the caller should use the parameter MAXSRF to
C                set this value.
C                
C 
C$ Detailed_Output
C
C     SHAPE      is a string describing the shape model. Values are
C                
C                    'ELLIPSOID' or 'DSK'
C
C
C     SUBTYP     is a string describing a sub-observer or sub-solar
C                point definition. Values are
C
C                    'INTERCEPT' or 'NEAR POINT'
C
C
C     PRI        is a logical flag indicating whether a DSK-based
C                computation uses DSK segments in a prioritized
C                fashion. PRI is .TRUE. if and only if the data
C                usage is prioritized. Note that for the N0066
C                version of SPICE, all DSK-based computations
C                use unprioritized data.
C
C
C     NSURF      is the number of elements of the surface ID list.
C
C
C     SRFLST     is a list of surface ID codes. All surfaces are
C                associated with the body designated by BODYID.
C
C
C     PNTDEF     is a string describing the point definition used
C                to compute limb or terminator points. Values are
C
C                   'TANGENT' or 'GUIDED'
C
C
C     TRMTYP     is a string that indicates the definition of the
C                terminator on the target body. Values are:
C
C                   'UMBRAL' or 'PENUMBRAL'
C
C$ Parameters
C
C     MAXSRF is the maximum surface list size used by the DSK 
C     subsystem. See dsk.inc for the parameter's value.
C
C$ Exceptions
C
C     1)  If more than one shape keyword is found in METHOD,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C     2)  If more than one surface keyword is found in METHOD,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C     3)  If the surface keyword is found in METHOD, but no
C         surfaces are listed, the error SPICE(BADMETHODSYNTAX) is
C         signaled.
C
C     4)  If the equals sign ('=') is missing from the surface list,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C     5)  If the method string contains a double quote character (")
C         that doesn't delimit a quoted-string, the error 
C          SPICE(BADMETHODSYNTAX) is signaled.
C
C     6)  If a surface name in the METHOD string can't be converted 
C         to a surface ID code, the error SPICE(IDCODENOTFOUND) is
C         signaled.
C
C     7)  If all ID codes in the surface list can't be stored in the
C         SRFLST array, the error SPICE(TOOMANYSURFACES) is signaled.
C
C     8)  If the surface list in METHOD contains a comma not followed
C         by a surface name, the error SPICE(BADMETHODSYNTAX) is
C         signaled.
C
C     9)  If the surface list in METHOD contains an unrecognized token,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    10)  If the string METHOD contains an extra prioritization token,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    12)  If the string METHOD contains an extra sub-point type token,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    13)  If the string METHOD contains an extra point definition
C         token, the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    14)  If the string METHOD contains an extra terminator type token,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    15)  If the string METHOD contains an extra slash character,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    16)  If the string METHOD lacks a prioritization keyword,
C         the error SPICE(BADPRIORITYSPEC) is signaled.
C
C    17)  If a "legacy" method string is not one of the expected
C         values, the error SPICE(BADMETHODSYNTAX) is signaled.
C
C    18)  If the string METHOD contains any other unexpected tokens,
C         the error SPICE(BADMETHODSYNTAX) is signaled.
C          
C$ Files
C    
C     Surface name-to-ID mappings may be defined at run time by loading
C     text kernels containing kernel variable assignments of the form
C
C        NAIF_SURFACE_NAME += ( <surface name 1>, ... )
C        NAIF_SURFACE_CODE += ( <surface code 1>, ... )
C        NAIF_SURFACE_BODY += ( <body code 1>,    ... )
C
C     Above, the Ith elements of the lists on the assignments' right
C     hand sides together define the Ith surface name/ID mapping.
C
C     The same effect can be achieved using assignments formatted as
C     follows:
C
C        NAIF_SURFACE_NAME += <surface name 1>
C        NAIF_SURFACE_CODE += <surface code 1>
C        NAIF_SURFACE_BODY += <body code 1>
C
C        NAIF_SURFACE_NAME += <surface name 2>
C        NAIF_SURFACE_CODE += <surface code 2>
C        NAIF_SURFACE_BODY += <body code 2>
C
C           ...
C
C     Note the use of the
C
C        +=
C
C     operator; this operator appends to rather than overwrites the
C     kernel variable named on the left hand side of the assignment.
C
C
C$ Particulars
C
C     This routine supports parsing of the METHOD input argument
C     used by the SPICE geometry APIs 
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
C     Each API has its own syntax for its METHOD input argument.
C     All routines use a common syntax for surface lists. See
C     the headers of the routines for the exact syntax supported
C     by each routine. 
C
C     In the METHOD string, the token delimiters are 
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
C     1) For the routine SUBPNT, 
C
C           Method string =
C
C              METHOD = 'INTERCEPT / '
C              //       'DSK/UNPRIORITIZED/SURFACES = "MGS '
C              //       'MOLA 128 PIXEL/DEG", MARS_LOWRES'
C
C          Output parameters:
C 
C             SHAPE  = 'DSK'
C             SUBTYP = 'INTERCEPT'
C             PRI    = .FALSE.
C             NSURF  = 2
C             SRFLST = <ID code 1>, <ID code 2>
C             PNTDEF = ' '
C             TRMTYP = ' '
C
C 
C     2) For the routine TERMPT, 
C
C           Method string =
C
C              METHOD = 'UMBRAL/TANGENT/ '
C              //       'DSK/UNPRIORITIZED/SURFACES = "MGS '
C              //       'MOLA 128 PIXEL/DEG", MARS_LOWRES'
C
C          Output parameters:
C 
C             SHAPE  = 'DSK'
C             SUBTYP = ' '
C             PRI    = .FALSE.
C             NSURF  = 2
C             SRFLST = <ID code 1>, <ID code 2>
C             PNTDEF = 'TANGENT'
C             TRMTYP = 'UMBRAL'
C
C
C$ Restrictions
C
C     This is a specialized parsing routine. It it meant to be used
C     only for the parsing "method" strings as described in
C     Particulars.
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
C-    SPICELIB Version 1.0.0, 24-FEB-2016 (NJB) 
C
C        Based on version  4.0.0 12-NOV-2015 (NJB)
C        First version was 1.0.0 14-JAN-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     parse method string for geometry api routines
C
C-&


C
C     SPICELIB functions
C      
      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         COMMA
      PARAMETER           ( COMMA  = ',' )

      CHARACTER*(*)         DQUOTE
      PARAMETER           ( DQUOTE = '"' )

      CHARACTER*(*)         EQUALS
      PARAMETER           ( EQUALS = '=' )

      CHARACTER*(*)         SLASH
      PARAMETER           ( SLASH  = '/' )

      INTEGER               NMAX
      PARAMETER           ( NMAX   = 100 )

      INTEGER               METLEN
      PARAMETER           ( METLEN = 1000 )

C
C     Local variables
C      
      CHARACTER*(1)         CHR
      CHARACTER*(METLEN)    LOCMTH

      INTEGER               BEGS ( NMAX + 2 )
      INTEGER               CODE
      INTEGER               ENDS ( NMAX + 2 )
      INTEGER               I
      INTEGER               N
      INTEGER               NPRI
      INTEGER               NSHAPE
      INTEGER               NLMBTP
      INTEGER               NSUBTP
      INTEGER               NSRFLS
      INTEGER               NTRMTP

      LOGICAL               DONE
      LOGICAL               FND


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZPRSMET' )


C
C     No shape or surfaces have been specified yet.
C
      SHAPE  = ' '
      SUBTYP = ' '
      PNTDEF = ' '
      TRMTYP = ' '
      PRI    = .TRUE.
      NSURF  = 0

C
C     Initialize clause counts.
C     
      NPRI   = 0
      NSHAPE = 0
      NSUBTP = 0
      NLMBTP = 0
      NTRMTP = 0
      NSRFLS = 0

      LOCMTH(1:1) = SLASH
      BEGS(1)     = 1
      ENDS(1)     = 1

      CALL LJUCRS ( 1, METHOD, LOCMTH(2:) )

C
C     Tokenize the input method string.
C
      CALL ZZLEXMET ( LOCMTH(2:), NMAX, N, BEGS(2), ENDS(2) )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZPRSMET' )
         RETURN
      END IF

      N = N + 1

      DO I = 2, N
         BEGS(I) = BEGS(I) + 1
         ENDS(I) = ENDS(I) + 1
      END DO

C
C     Identify the target shape.
C
      DONE = .FALSE.
      I    = 1

      DO WHILE ( .NOT. DONE ) 

         IF (  LOCMTH( BEGS(I):ENDS(I) ) .EQ.  'ELLIPSOID' ) THEN

            NSHAPE = NSHAPE + 1
            SHAPE  = 'ELLIPSOID'

         ELSE IF (  LOCMTH( BEGS(I):ENDS(I) ) .EQ.  'DSK' ) THEN

            NSHAPE = NSHAPE + 1
            SHAPE  = 'DSK'

         END IF

         IF ( I .EQ. N ) THEN

            DONE = .TRUE.
         ELSE
            I    = I + 1
         END IF

      END DO


      IF ( NSHAPE .NE. 1 ) THEN

         CALL SETMSG ( 'The "method" or "shape" string must contain '
     .   //            'exactly one instance of a shape keyword: '
     .   //            'ELLIPSOID or DSK, or additionally, in the '
     .   //            'case of occultation computations, POINT. '
     .   //            'The number of shape keywords was #. '
     .   //            'The input string was <#>.'                  )
         CALL ERRINT ( '#', NSHAPE                                  )
         CALL ERRCH  ( '#', METHOD                                  )
         CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'                     )
         CALL CHKOUT ( 'ZZPRSMET'                                   )
         RETURN

      END IF


      IF ( SHAPE .EQ. 'DSK' ) THEN
C
C        Traverse the tokens; identify clauses in the method string.
C
         I = 1

         DO WHILE ( I .LT. N )

            CHR = LOCMTH( BEGS(I):ENDS(I) )

            IF ( CHR .EQ. SLASH ) THEN
C
C              If the method string is correct, the next token should
C              be a keyword.
C
               I = I + 1
              
               IF ( LOCMTH( BEGS(I):ENDS(I) )  .EQ.  'SURFACES' ) THEN
C
C                 Normally, we're looking at the start of a surface
C                 list at this point. We need at least an assignment
C                 operator ('=') and a surface name or ID code following
C                 the SURFACES keyword.
C
                  NSRFLS = NSRFLS + 1

                  IF ( NSRFLS .GT. 1 ) THEN
C
C                    We have too many surface specifications.
C
                     CALL SETMSG ( 'Extra SURFACES keyword was found '
     .               //            'in the method string <#>.'       )
                     CALL ERRCH  ( '#', METHOD                       )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'          )
                     CALL CHKOUT ( 'ZZPRSMET'                        )
                     RETURN

                  END IF

                  IF ( I .GT. N-2 ) THEN

                     CALL SETMSG ( 'The surface list in the method '
     .               //            'string appears to be truncated. '
     .               //            'The method string is <#>.'       )
                     CALL ERRCH  ( '#', METHOD                       )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'          )
                     CALL CHKOUT ( 'ZZPRSMET'                        )
                     RETURN

                  END IF

                  I   = I + 1

                  CHR = LOCMTH( BEGS(I):ENDS(I) )

                  IF ( CHR .NE. EQUALS ) THEN

                     CALL SETMSG ( 'The surface list in the method '
     .               //            'string lacks an "=" sign after '
     .               //            'the SURFACES keyword. '
     .               //            'The method string is <#>.'      )
                     CALL ERRCH  ( '#', METHOD                      )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'         )
                     CALL CHKOUT ( 'ZZPRSMET'                       )
                     RETURN

                  END IF
                  
C
C                 Prepare to read a list of names, quoted
C                 strings, or ID codes. The list must be non-empty. If
C                 there are multiple list entries, they're delimited by
C                 commas.
C
                  I    = I + 1
                  DONE = .FALSE.
                  
                  DO WHILE ( .NOT. DONE )

                     CHR = LOCMTH( BEGS(I):BEGS(I) )

                     IF ( CHR .EQ. DQUOTE ) THEN
C
C                       We expect the current token to be a quoted
C                       string.
C
                        IF ( ENDS(I) .GT. BEGS(I)+1 ) THEN

                           BEGS(I) = BEGS(I) + 1 
                           ENDS(I) = ENDS(I) - 1

                        ELSE
C
C                          We have an invalid quoted string.
C
                           CALL SETMSG ( 'The surface list in the '
     .                     //            'method string '
     .                     //            'contains a double quote '
     .                     //            'character that is not a '
     .                     //            'delimiter of a valid quoted '
     .                     //            'string. '
     .                     //            'The method string is <#>.'   )
                           CALL ERRCH  ( '#', METHOD                   )
                           CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'      )
                           CALL CHKOUT ( 'ZZPRSMET'                    )
                           RETURN

                        END IF

                     END IF

C                    We have either a name or an integer in string form.
C                    Convert the token to a surface ID code.
C
                     CALL SRFSCC ( LOCMTH(BEGS(I):ENDS(I)), 
     .                             BODYID,   CODE,   FND   )

                     IF ( FAILED() ) THEN
                        CALL CHKOUT ( 'ZZPRSMET' )
                        RETURN
                     END IF

                     IF ( .NOT. FND ) THEN

                        CALL SETMSG ( 'The surface name <#> could not '
     .                  //            'be translated to an ID code. '   
     .                  //            'The method string is <#>.'     )
                        CALL ERRCH  ( '#', LOCMTH(BEGS(I):ENDS(I))    )
                        CALL ERRCH  ( '#', METHOD                     )
                        CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'         )
                        CALL CHKOUT ( 'ZZPRSMET'                      )
                        RETURN

                     END IF


                     NSURF = NSURF + 1
C
C                    Make sure there's room in the surface ID array.
C
                     IF ( NSURF .GT. MXNSRF ) THEN

                        CALL SETMSG ( 'The surface name <#> could not '
     .                  //            'be stored in the surface list '
     .                  //            'due to lack of room. The max'
     .                  //            'imum number of surfaces that '
     .                  //            'can be specified is #. '
     .                  //            'The method string is <#>.'     )
                        CALL ERRCH  ( '#', LOCMTH(BEGS(I):ENDS(I))    )
                        CALL ERRINT ( '#', MXNSRF                     )
                        CALL ERRCH  ( '#', METHOD                     )
                        CALL SIGERR ( 'SPICE(TOOMANYSURFACES)'        )
                        CALL CHKOUT ( 'ZZPRSMET'                      )
                        RETURN

                     END IF

C
C                    Append the surface ID to the surface ID array.
C
                     SRFLST( NSURF ) = CODE
                     

                     IF ( I .EQ. N ) THEN
C
C                       We're at the end of the method string.
C                       
                        DONE = .TRUE.

                     ELSE
C
C                       There are more tokens; the surface list may 
C                       contain more surface names or ID codes.
C
                        I   = I + 1
                        CHR = LOCMTH( BEGS(I):ENDS(I) )

                        IF ( CHR .EQ. SLASH ) THEN
C
C                          We're at the end of the surface list. 
C                          
                           DONE = .TRUE.
C
C                          Decrement the token pointer so that the
C                          slash will be seen as the next token.
C
                           I    = I - 1


                        ELSE IF ( CHR .EQ. COMMA ) THEN
C
C                          We expect to find another surface name or
C                          ID in the list.
C
                           IF ( I .LT. N ) THEN

                              I = I + 1

                           ELSE
C
C                             We have a syntax error in the surface
C                             list.
C
                              CALL SETMSG ( 'The surface list in the '
     .                        //            'method string contains '
     .                        //            'a comma that is not '
     .                        //            'followed by a surface '
     .                        //            'name or ID code. The '
     .                        //            'method string is <#>.'   )
                              CALL ERRCH  ( '#', METHOD               )
                              CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'  )
                              CALL CHKOUT ( 'ZZPRSMET'                )
                              RETURN

                           END IF


                        ELSE
C
C                          We have a syntax error in the surface list.
C
                           CALL SETMSG ( 'The surface list in the '
     .                     //            'method string is followed '
     .                     //            'by the unexpected token <#>.'
     .                     //            ' The method string is <#>.'  )
                           CALL ERRCH  ( '#', LOCMTH(BEGS(I):ENDS(I))  )
                           CALL ERRCH  ( '#', METHOD                   )
                           CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'      )
                           CALL CHKOUT ( 'ZZPRSMET'                    )
                           RETURN

                        END IF

                     END IF

                  END DO


               ELSE IF (      LOCMTH( BEGS(I):ENDS(I) )  
     .                   .EQ. 'UNPRIORITIZED'            ) THEN


                  NPRI = NPRI + 1

                  IF ( NPRI .GT. 1 ) THEN
C
C                    We have too many prioritization specifications.
C
                     CALL SETMSG ( 'Extra prioritization keyword was '
     .               //            'found in the method string <#>.' )
                     CALL ERRCH  ( '#', METHOD                       )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'          )
                     CALL CHKOUT ( 'ZZPRSMET'                        )
                     RETURN

                  END IF


                  PRI  = .FALSE.


               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'INTERCEPT'                ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'NADIR'                    )  ) THEN
C
C                 This is a sub-point type specification.
C                 
                  NSUBTP = NSUBTP + 1
 
                  IF ( NSUBTP .GT. 1 ) THEN
C
C                    We have too many sub-point specifications.
C
                     CALL SETMSG ( 'Extra sub-point type <#> was found '
     .               //            'in the method string <#>.'         )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  SUBTYP = LOCMTH( BEGS(I):ENDS(I) )


               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'TANGENT'                ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'GUIDED'                   )  ) THEN
C
C                 This is a point definition specification.
C                 
                  NLMBTP = NLMBTP + 1
 
                  IF ( NLMBTP .GT. 1 ) THEN
C
C                    We have too many point definition specifications.
C
                     CALL SETMSG ( 'Extra point definition <#> was '
     .               //            'found in the method string <#>.'   )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  PNTDEF = LOCMTH( BEGS(I):ENDS(I) )



               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'UMBRAL'                   ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'PENUMBRAL'                )  ) THEN
C
C                 This is a terminator type specification.
C                 
                  NTRMTP = NTRMTP + 1
 
                  IF ( NTRMTP .GT. 1 ) THEN
C
C                    We have too many terminator type specifications.
C
                     CALL SETMSG ( 'Extra terminator type <#> was '
     .               //            'found in the method string <#>.'   )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  TRMTYP = LOCMTH( BEGS(I):ENDS(I) )



               ELSE IF (  LOCMTH( BEGS(I):ENDS(I) ) .EQ. SLASH  )  THEN
C
C                 Adjust the token index in the message to account
C                 for the insertion of a slash at the beginning.
C
                  CALL SETMSG ( 'An unexpected slash character was '
     .            //            'found at index # in the method '
     .            //            'string <#>.'                          )
                  CALL ERRINT ( '#', BEGS(I)-1                         )
                  CALL ERRCH  ( '#', METHOD                            )
                  CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'               )
                  CALL CHKOUT ( 'ZZPRSMET'                             )
                  RETURN
              

               ELSE IF ( LOCMTH( BEGS(I):ENDS(I) )  .NE.  'DSK'  ) THEN
C
C                 'DSK' was the only valid token that could appear
C                 at this point.
C
                  CALL SETMSG ( 'Unexpected token <#> was found '
     .            //            'in the method string <#>.'       )
                  CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )    )
                  CALL ERRCH  ( '#', METHOD                       )
                  CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'          )
                  CALL CHKOUT ( 'ZZPRSMET'                        )
                  RETURN        

               END IF


            END IF

            I = I + 1

         END DO

C
C        Currently, the UNPRIORITIZED keyword must be provided 
C        for DSK shapes.
C
         IF ( PRI ) THEN

            CALL SETMSG ( 'The keyword UNPRIORITIZED must be '
     .      //            'present in the method string <#> '
     .      //            'when the target shape is DSK. '
     .      //            'Prioritized DSK segment searches '
     .      //            'are not currently supported.'      )
            CALL ERRCH  ( '#', METHOD                         )
            CALL SIGERR ( 'SPICE(BADPRIORITYSPEC)'            )
            CALL CHKOUT ( 'ZZPRSMET'                          )
            RETURN        

         END IF


      ELSE IF ( SHAPE .EQ. 'ELLIPSOID' ) THEN
C
C        Check for legacy string inputs.
C        
         IF ( EQSTR(METHOD, 'ELLIPSOID') ) THEN
C
C           There are no other outputs to set.
C
            CALL CHKOUT ( 'ZZPRSMET' )
            RETURN        

         ELSE IF (      EQSTR(METHOD, 'NEARPOINT:ELLIPSOID')
     .             .OR. EQSTR(METHOD, 'ELLIPSOID:NEARPOINT') 
     .             .OR. EQSTR(METHOD, 'NEARPOINT/ELLIPSOID') 
     .             .OR. EQSTR(METHOD, 'ELLIPSOID/NEARPOINT')  ) THEN

            SUBTYP = 'NEAR POINT'

            CALL CHKOUT ( 'ZZPRSMET' )
            RETURN        

         ELSE IF (      EQSTR(METHOD, 'INTERCEPT:ELLIPSOID')
     .             .OR. EQSTR(METHOD, 'ELLIPSOID:INTERCEPT') 
     .             .OR. EQSTR(METHOD, 'INTERCEPT/ELLIPSOID') 
     .             .OR. EQSTR(METHOD, 'ELLIPSOID/INTERCEPT')  ) THEN

            SUBTYP = 'INTERCEPT'

            CALL CHKOUT ( 'ZZPRSMET' )
            RETURN        
            
         END IF

C
C        At this point, we should have a "modern" style of method
C        specification for an ellipsoidal shape model.
C        Parse the method string.
C
C        Traverse the tokens; identify clauses in the method string.
C
         I = 1

         DO WHILE ( I .LE. N )

            CHR = LOCMTH( BEGS(I):ENDS(I) )

            IF ( CHR .EQ. SLASH ) THEN
C
C              If the method string is correct, the next token should
C              be a keyword. There had better be a next token.
C
               IF ( I .EQ. N ) THEN

                  CALL SETMSG ( 'Expected more tokens after final <#> '
     .               //         'in the method string <#>.'            )
                  CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )         )
                  CALL ERRCH  ( '#', METHOD                            )
                  CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'               )
                  CALL CHKOUT ( 'ZZPRSMET'                             )
                  RETURN

               END IF

               I = I + 1


               IF (  LOCMTH( BEGS(I):ENDS(I) ) .EQ.  'ELLIPSOID' ) THEN
C
C                 This case is a no-op.
C
 
               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'INTERCEPT'                ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'NADIR'                    )  
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'NEAR'                     )  ) THEN
C
C                 This is a sub-point type specification.
C                 
                  NSUBTP = NSUBTP + 1

                  IF ( NSUBTP .GT. 1 ) THEN
C
C                    We have too many sub-point specifications.
C
                     CALL SETMSG ( 'Extra sub-point type <#> was found '
     .               //            'in the method string <#>.'         )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  IF ( LOCMTH( BEGS(I):ENDS(I) ) .EQ. 'NEAR' ) THEN
C
C                    This may be a combination of old style and
C                    new syntax, e.g. 'ELLIPSOID/NEAR POINT'.
C                     
                     IF (       ( SHAPE .EQ. 'ELLIPSOID' )
     .                    .AND. ( I     .LT.  N          )  ) THEN
C
C                       We allow the "near point" sub-point type
C                       for ellipsoids but not for DSK surfaces.
C
                        IF (     LOCMTH( BEGS(I+1) : ENDS(I+1) ) 
     .                      .EQ. 'POINT'                          ) THEN

                           SUBTYP = 'NADIR'

                        END IF

                     END IF

                  ELSE
                     SUBTYP = LOCMTH( BEGS(I):ENDS(I) )
                  END IF





               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'TANGENT'                ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'GUIDED'                   )  ) THEN
C
C                 This is a point definition specification.
C                 

                  NLMBTP = NLMBTP + 1
 
                  IF ( NLMBTP .GT. 1 ) THEN
C
C                    We have too many point definition specifications.
C
                     CALL SETMSG ( 'Extra point definition <#> was '
     .               //            'found in the method string <#>.'   )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  PNTDEF = LOCMTH( BEGS(I):ENDS(I) )

               ELSE IF (     (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'UMBRAL'                   ) 
     .                  .OR. (       LOCMTH( BEGS(I):ENDS(I) )  
     .                         .EQ. 'PENUMBRAL'                )  ) THEN
C
C                 This is a terminator type specification.
C                 
                  NTRMTP = NTRMTP + 1
 
                  IF ( NTRMTP .GT. 1 ) THEN
C
C                    We have too many terminator type specifications.
C
                     CALL SETMSG ( 'Extra terminator type <#> was '
     .               //            'found in the method string <#>.'   )
                     CALL ERRCH  ( '#', LOCMTH( BEGS(I):ENDS(I) )      )
                     CALL ERRCH  ( '#', METHOD                         )
                     CALL SIGERR ( 'SPICE(BADMETHODSYNTAX)'            )
                     CALL CHKOUT ( 'ZZPRSMET'                          )
                     RETURN

                  END IF

                  TRMTYP = LOCMTH( BEGS(I):ENDS(I) )

               ELSE
C
C                 Sorry, no other strings are allowed.
C 
                  CALL SETMSG( 'Unexpected method string <#> was found '
     .            //           'specifying an ellipsoid shape.'        )
                  CALL ERRCH ( '#', METHOD                             )
                  CALL SIGERR( 'SPICE(BADMETHODSYNTAX)'                )
                  CALL CHKOUT( 'ZZPRSMET'                              )
                  RETURN        

               END IF

            ELSE

               CALL SETMSG( 'Unexpected method string <#> was found '
     .         //           'specifying an ellipsoid shape.'        )
               CALL ERRCH ( '#', METHOD                             )
               CALL SIGERR( 'SPICE(BADMETHODSYNTAX)'                )
               CALL CHKOUT( 'ZZPRSMET'                              )
               RETURN                 

            END IF



            I = I + 1

         END DO          

      ELSE
C
C        This is a backstop error check.
C
         CALL SETMSG ( 'Unexpected shape value # was found. '
     .   //            'This is due to a coding error, not '
     .   //            'to user input.'                      )
         CALL ERRCH  ( '#', SHAPE                            )
         CALL SIGERR ( 'SPICE(BUG)'                          )
         CALL CHKOUT ( 'ZZPRSMET'                            )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZPRSMET' )
      RETURN
      END
