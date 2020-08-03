C$Procedure SRFCSS ( Surface ID and body string to surface string )

      SUBROUTINE SRFCSS ( CODE, BODSTR, SRFSTR, ISNAME )

C$ Abstract
C
C     Translate a surface ID code, together with a body string, to the
C     corresponding surface name. If no such surface name exists,
C     return a string representation of the surface ID code.
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

      INTEGER               CODE
      CHARACTER*(*)         BODSTR
      CHARACTER*(*)         SRFSTR
      LOGICAL               ISNAME
  
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     CODE       I   Integer surface ID code to translate to a string.
C     BODSTR     I   Name or ID of body associated with surface.
C     SRFSTR     O   String corresponding to surface ID code.
C     ISNAME     O   Flag indicating whether output is a surface name.
C     SFNMLN     P   Maximum length of surface name.
C
C$ Detailed_Input
C
C     CODE       is an integer ID code for a surface associated with a
C                specified body.
C
C
C     BODSTR     is a string designating the body associated with the
C                input surface ID code. BODSTR may contain a body name
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
C                significant. Sequences of consecutive embedded blanks
C                are considered equivalent to a single blank. That is,
C                all of the following strings are equivalent names:
C
C                   '67P/CHURYUMOV-GERASIMENKO (1969 R1)'
C                   '67P/Churyumov-Gerasimenko (1969 R1)'
C                   '67P/CHURYUMOV-GERASIMENKO (1969 R1)   '
C                   '67P/CHURYUMOV-GERASIMENKO    (1969 R1)'
C                   '   67P/CHURYUMOV-GERASIMENKO (1969 R1)'
C
C                However, '67P/CHURYUMOV-GERASIMENKO(1969R1)' 
C                is not equivalent to the names above.
C
C
C$ Detailed_Output
C
C     SRFSTR     the name of the surface identified by CODE, for the
C                body designated by BODSTR, if for this body an
C                association exists between the input surface ID and a
C                surface name.
C
C                If CODE has more than one translation, then the most
C                recently defined surface name corresponding to CODE is
C                returned. SRFSTR will have the exact format (case and
C                embedded blanks) used in the definition of the
C                name/code association.
C
C                If the input surface ID code and body name do not map
C                to a surface name, SRFSTR is set to the string
C                representation of CODE.
C
C                SRFSTR should be declared with length SFNMLN (see the
C                Parameters section below).
C
C
C     ISNAME     is a logical flag that is .TRUE. if a surface name
C                corresponding to the input ID codes was found and
C                .FALSE. otherwise. When ISNAME is .FALSE., the output
C                string SRFSTR contains a string representing the 
C                integer CODE.
C
C$ Parameters
C
C     SFNMLN     is the maximum length of a surface name. This 
C                parameter is declared in the SPICELIB include file 
C
C                   srftrn.inc
C
C$ Exceptions
C
C     1)  If the input body string cannot be mapped to a body name, the
C         output SRFSTR is set to a string representation of the
C         surface ID code. The output ISNAME is set to .FALSE.
C
C         This case is not treated as an error.
C
C     2)  If the input surface code cannot be mapped to a surface name,
C         the output SRFSTR is set to a string representation of the
C         surface ID code. The input body string is ignored. The output
C         ISNAME is set to .FALSE.
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
C     SRFCSS is one of four related subroutines:
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
C        Translate each surface and body ID code pair to the 
C        associated surface name. Also perform a translation
C        for a surface ID having no matching name.
C                    
C        Use the meta-kernel shown below to define the required SPICE
C        kernel variables:
C
C
C           KPL/MK
C
C           File: srfcss_ex1.tm
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
C          PARAMETER           ( NCASE = 5 )
C
C          INTEGER               BDNMLN
C          PARAMETER           ( BDNMLN = 36 )
C
C          CHARACTER*(BDNMLN)    BODSTR ( NCASE )
C          CHARACTER*(FILSIZ)    META
C          CHARACTER*(SFNMLN)    SRFNAM
C
C          INTEGER               I
C          INTEGER               SURFID ( NCASE )
C
C          LOGICAL               ISNAME
C
C
C          DATA  ( SURFID(I), BODSTR(I), I = 1, NCASE ) /
C         .
C         .        1,         'MARS',
C         .        1,         'PHOBOS',
C         .        2,         '499',
C         .        3,         'MARS',
C         .        1,         'ZZZ'                     /
C
C
C          META = 'srfcss_ex1.tm'
C
C          CALL FURNSH ( META )
C
C          WRITE (*,*) ' '
C
C          DO I = 1, NCASE
C
C             CALL SRFCSS ( SURFID(I), BODSTR(I),
C         .                 SRFNAM,    ISNAME     )
C
C             WRITE (*,*) 'surface ID     = ', SURFID(I)
C             WRITE (*,*) 'body string    = ', BODSTR(I)
C             WRITE (*,*) 'name found     = ', ISNAME
C             WRITE (*,*) 'surface string = ', SRFNAM
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
C        surface ID     =            1
C        body string    = MARS
C        name found     =  T
C        surface string = MGS MOLA  64 pixel/deg
C
C        surface ID     =            1
C        body string    = PHOBOS
C        name found     =  T
C        surface string = PHOBOS GASKELL Q512
C
C        surface ID     =            2
C        body string    = 499
C        name found     =  T
C        surface string = MGS MOLA 128 pixel/deg
C
C        surface ID     =            3
C        body string    = MARS
C        name found     =  F
C        surface string = 3
C
C        surface ID     =            1
C        body string    = ZZZ
C        name found     =  F
C        surface string = 1
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
C     surface ID code and body string to surface string
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
      INTEGER               BODYID


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'SRFCSS' )

C
C     No name has been found yet.
C
      ISNAME = .FALSE.

C
C     Convert the body string to an ID code.
C
      CALL BODS2C ( BODSTR, BODYID, ISNAME )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'SRFCSS' )
         RETURN
      END IF

      IF ( ISNAME ) THEN
C
C        Try to translate the surface and body codes to a known surface
C        name.
C      
         CALL ZZSRFC2N ( CODE, BODYID, SRFSTR, ISNAME )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'SRFCSS' )
            RETURN
         END IF 

      END IF
         
C
C     If either the body string or surface code could not be
C     translated, convert the surface ID code to a string
C     representation.
C
      IF ( .NOT. ISNAME ) THEN

         CALL INTSTR ( CODE, SRFSTR )

      END IF

      CALL CHKOUT ( 'SRFCSS' )
      RETURN
      END
