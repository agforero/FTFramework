C$Procedure SRFS2C ( Surface and body strings to surface ID code )

      SUBROUTINE SRFS2C ( SRFSTR, BODSTR, CODE, FOUND )

C$ Abstract
C
C     Translate a surface string, together with a body string, to the
C     corresponding surface ID code. The input strings may contain
C     names or integer ID codes.
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
C     NAIF_IDS
C
C$ Keywords
C
C     CONVERSION
C     DSK
C     ID
C     NAME
C     STRING
C     SURFACE
C
C$ Declarations

      IMPLICIT NONE

      CHARACTER*(*)         SRFSTR
      CHARACTER*(*)         BODSTR
      INTEGER               CODE
      LOGICAL               FOUND
  
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     SRFSTR     I   Surface name or ID string.
C     BODSTR     I   Body name or ID string.
C     CODE       O   Integer surface ID code.
C     FOUND      O   Flag indicating whether surface ID was found.
C
C$ Detailed_Input
C
C     SRFSTR     is a string designating a surface. SRFSTR may contain
C                a surface name or a string representation of the
C                surface's integer ID code.
C
C                Case and leading and trailing blanks in a surface name
C                are not significant. Sequences of consecutive embedded
C                blanks are considered equivalent to a single blank.
C                For example, all of the strings below are considered
C                to be equivalent:
C
C                   'MGS MOLA 128 PIXEL/DEG'
C                   'MGS MOLA 128 pixel/deg'
C                   'MGS MOLA 128 PIXEL/DEG   '
C                   'MGS MOLA 128    PIXEL/DEG'
C                   '   MGS MOLA 128 PIXEL/DEG'
C
C                However, 
C
C                   'MGSMOLA 128PIXEL/DEG' 
C
C                is not equivalent to the names above.
C
C
C     BODSTR     is a string designating the body associated with the
C                input surface string. BODSTR may contain a body name
C                or a string representation of the body's integer ID
C                code. For example, BODSTR may contain
C
C                   '1000012'
C
C                instead of  
C
C                   '67P/CHURYUMOV-GERASIMENKO (1969 R1)'
C
C                Case and leading and trailing blanks in a name are not
C                significant. The treatment of blanks in BODSTR is the
C                same as for SRFSTR.
C
C$ Detailed_Output
C
C     CODE       is integer ID code of the surface designated by
C                SRFSTR, for the body designated by BODSTR, if for this
C                body an association exists between the input surface
C                string and a surface ID code. CODE is defined if and
C                only if the output flag FOUND is .TRUE.
C
C     FOUND      is a logical flag that is .TRUE. if a surface code
C                corresponding to the input strings was found and
C                .FALSE. otherwise.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the input surface string does not map to an ID code
C         and does not represent an integer, the output CODE is
C         undefined and the output FOUND is set to .FALSE.
C
C         This case is not treated as an error.
C
C     2)  If the input body string does not map to an ID code and does
C         not represent an integer, the output CODE is undefined and
C         the output FOUND is set to .FALSE.
C
C         This case is not treated as an error.
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
C$ Particulars
C
C     Surfaces are always associated with bodies (which usually are
C     ephemeris objects). For any given body, a mapping between surface
C     names and surface ID codes can be established. 
C
C     Bodies serve to disambiguate surface names and ID codes: the set
C     of surface names and surface ID codes for a given body can be
C     thought of as belonging to a name space. A given surface ID code
C     or surface name may be used for surfaces of multiple bodies,
C     without conflict.
C
C     Associations between surface names and ID codes are always made
C     via kernel pool assignments; there are no built-in associations.
C
C     SRFS2C is one of four related subroutines:
C
C        SRFS2C      Surface string and body string to surface ID code
C        SRFSCC      Surface string and body ID code to surface ID code
C        SRFC2S      Surface ID code and body ID code to surface string
C        SRFCSS      Surface ID code and body string to surface string
C
C     SRFS2C, SRFC2S, SRFSCC, and SRFCSS perform translations between 
C     surface strings and their corresponding integer ID codes.
C
C     Refer to naif_ids.req for details concerning adding new surface
C     name/code associations at run time by loading text kernels.
C
C$ Examples
C
C     The formatting of the results shown for this example may differ
C     across platforms.
C
C     1) Supposed a text kernel has been loaded that contains
C        the following assignments:
C
C           NAIF_SURFACE_NAME += ( 'MGS MOLA  64 pixel/deg',
C                                  'MGS MOLA 128 pixel/deg',
C                                  'PHOBOS GASKELL Q512'     )
C           NAIF_SURFACE_CODE += (   1,   2,    1 )
C           NAIF_SURFACE_BODY += ( 499, 499,  401 )
C
C        Translate each surface and body string pair to the 
C        associated surface ID code. Also perform a translation
C        for a surface name having no matching ID and for 
C        a body string having no matching body ID code.
C                    
C        Use the meta-kernel shown below to define the required SPICE
C        kernel variables.
C
C
C           KPL/MK
C
C           File: srfs2c_ex1.tm
C
C           This meta-kernel is intended to support operation of SPICE
C           example programs. The file contents shown here should not be
C           assumed to contain adequate or correct versions of data
C           required by SPICE-based user applications.
C
C
C           \begindata
C
C           NAIF_SURFACE_NAME += ( 'MGS MOLA  64 pixel/deg',
C                                  'MGS MOLA 128 pixel/deg',
C                                  'PHOBOS GASKELL Q512'     )
C           NAIF_SURFACE_CODE += (   1,   2,    1 )
C           NAIF_SURFACE_BODY += ( 499, 499,  401 )
C
C           \begintext
C
C
C       Example code begins here.
C    
C
C          PROGRAM EX1
C          IMPLICIT NONE
C
C          INCLUDE 'srftrn.inc'
C
C          INTEGER               FILSIZ
C          PARAMETER           ( FILSIZ = 255 )
C
C          INTEGER               NCASE
C          PARAMETER           ( NCASE = 8 )
C
C          INTEGER               BDNMLN
C          PARAMETER           ( BDNMLN = 36 )
C
C          CHARACTER*(BDNMLN)    BODSTR ( NCASE )
C          CHARACTER*(FILSIZ)    META
C          CHARACTER*(SFNMLN)    SRFSTR ( NCASE )
C
C          INTEGER               I
C          INTEGER               SURFID
C
C          LOGICAL               FOUND
C
C
C          DATA  ( SRFSTR(I), BODSTR(I), I = 1, NCASE )  /
C         .
C         .        'MGS MOLA  64 pixel/deg',    'MARS',
C         .        'PHOBOS GASKELL Q512',       'PHOBOS',
C         .        'MGS MOLA 128 pixel/deg',    'MARS',
C         .        'MGS MOLA  64 pixel/deg',    '499',
C         .        '1',                         'PHOBOS',
C         .        '2',                         '499',
C         .        'ZZZ',                       'MARS',
C         .        '1',                         'ZZZ'    /
C
C
C          META = 'srfs2c_ex1.tm'
C
C          CALL FURNSH ( META )
C
C          WRITE (*,*) ' '
C
C          DO I = 1, NCASE
C
C             CALL SRFS2C ( SRFSTR(I), BODSTR(I),
C         .                 SURFID,    FOUND     )
C
C             WRITE (*,*) 'surface string   = ', SRFSTR(I)
C             WRITE (*,*) 'body string      = ', BODSTR(I)
C             WRITE (*,*) 'surface ID found = ', FOUND
C
C             IF ( FOUND ) THEN
C                WRITE (*,*) 'surface ID       = ', SURFID
C             END IF
C
C             WRITE (*,*) ' '
C
C          END DO
C
C          END
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, the output was:
C
C
C        surface string   = MGS MOLA  64 pixel/deg
C        body string      = MARS
C        surface ID found =  T
C        surface ID       =            1
C
C        surface string   = PHOBOS GASKELL Q512
C        body string      = PHOBOS
C        surface ID found =  T
C        surface ID       =            1
C
C        surface string   = MGS MOLA 128 pixel/deg
C        body string      = MARS
C        surface ID found =  T
C        surface ID       =            2
C
C        surface string   = MGS MOLA  64 pixel/deg
C        body string      = 499
C        surface ID found =  T
C        surface ID       =            1
C
C        surface string   = 1
C        body string      = PHOBOS
C        surface ID found =  T
C        surface ID       =            1
C
C        surface string   = 2
C        body string      = 499
C        surface ID found =  T
C        surface ID       =            2
C
C        surface string   = ZZZ
C        body string      = MARS
C        surface ID found =  F
C
C        surface string   = 1
C        body string      = ZZZ
C        surface ID found =  F
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
C     B.V. Semenov    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 14-JAN-2016 (NJB) (EDW) (BVS)
C
C-&

C$ Index_Entries
C
C     surface string and body string to surface ID code 
C
C-&


C
C     SPICELLIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local parameters
C
      INTEGER               MSGLEN
      PARAMETER           ( MSGLEN = 80 )

C
C     Local variables
C 
      CHARACTER*(MSGLEN)    ERRMSG

      INTEGER               BODYID
      INTEGER               PTR



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'SRFS2C' )

C
C     No result has been found yet.
C
      FOUND = .FALSE.

C
C     Obtain a body ID string.
C
      CALL BODS2C ( BODSTR, BODYID, FOUND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'SRFS2C' )
         RETURN
      END IF

      IF ( FOUND ) THEN
C
C        We have the body ID.
C
C        Try to translate the input name to a known surface ID code.
C      
         CALL ZZSRFN2C ( SRFSTR, BODYID, CODE, FOUND )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'SRFS2C' )
            RETURN
         END IF

         IF ( .NOT. FOUND ) THEN
C
C           The surface string could not be mapped to an ID code.
C
C           It's possible the name is a string representation of an
C           integer, for example, '999'.  If so, find the equivalent
C           datum of INTEGER type.
C
            CALL NPARSI ( SRFSTR, CODE, ERRMSG, PTR )
C
C           If the string parsed as an integer, PTR is zero;
C           otherwise it's non-zero.
C
            FOUND = PTR .EQ. 0
 
         END IF

      END IF
C
C     CODE and FOUND are set.
C
      CALL  CHKOUT ( 'SRFS2C' )
      RETURN
      END
