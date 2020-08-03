C$Procedure PRSCML ( Parse MKDSK command line )
 
      SUBROUTINE PRSCML ( CMDLIN, INFO,  INFTYP, 
     .                    SETUP,  INFIL, OUTFIL )
 
C$ Abstract
C
C     Parse the command line arguments supplied to MKDSK.
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
C     PCK
C     TIME
C
C$ Keywords
C
C     TOPOGRAPHY
C     FILES
C
C$ Declarations
 
      IMPLICIT NONE 

      INCLUDE 'mkdsk.inc'
      INCLUDE 'errhnd.inc'

      CHARACTER*(*)         CMDLIN
      LOGICAL               INFO
      CHARACTER*(*)         INFTYP
      CHARACTER*(*)         SETUP
      CHARACTER*(*)         INFIL
      CHARACTER*(*)         OUTFIL

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     CMDLIN     I   Command line to be parsed.
C     INFO       O   Flag indicating whether command is an info request.
C     INFTYP     O   Information type if information is requested.
C     SETUP      O   Setup file name.
C     INFIL      O   Input file name.
C     OUTFIL     O   Output file name.
C     
C$ Detailed_Input
C
C     CMDLIN         is a string containing the command line arguments
C                    supplied when MKDSK was invoked.
C
C$ Detailed_Output
C
C     INFO           is a logical flag indicating whether the command
C                    is a request for information:  usage, "help,", or
C                    the program's version.
C
C                    When INFO is .TRUE., INFTYP will indicate the type
C                    of information requested.
C
C     INFTYP         is a string indicating the type of information
C                    requested, if the input command requests such.
C                    Values of INFTYP are:
C                     
C                       'HELP'
C                       'VERSION'
C                       'USAGE'
C                       'TEMPLATE'
C
C                    INFTYP is valid only if INFO is .TRUE.  Otherwise,
C                    IFNTYP is returned blank.
C
C     SETUP          is the name of the setup file specified on the 
C                    command line.  If no file is specified, and
C                    if the command is not an information request, the
C                    user is prompted for a file name.
C
C                    If INFO is returned .TRUE., SETUP is left blank.
C
C     INFIL          is the name of the input file specified on the
C                    command line.  If no input file name is specified,
C                    INFIL is left blank.
C                    
C                    If INFO is returned .TRUE., INFIL is left blank.
C
C     OUTFIL         is the name of the output file specified on the
C                    command line.  If no input file name is specified,
C                    OUTFIL is left blank.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the input command is syntactically invalid, the error
C        SPICE(CMDERROR) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The expected command syntax is:
C
C        mkdsk   [-setup <setup file name>]
C                [-input <input data file name>]
C                [-output <output SPK file name>]
C                [-h|-help]
C                [-t|-template]
C                [-u|-usage]
C                [-v|-version]
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1) This routine is intended for use solely within the MKDSK
C        program.
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
C-    MKDSK Version 1.0.0, 30-MAR-2010 (NJB)
C
C-&
 
C
C     SPICELIB functions
C
      INTEGER               RTRIM

      LOGICAL               RETURN

C
C     Non-SPICELIB  functions
C
      LOGICAL               M2XIST

C
C     Local parameters
C     
      INTEGER               MAXMSG
      PARAMETER           ( MAXMSG = 2 )

      INTEGER               MAXSYN
      PARAMETER           ( MAXSYN = 12 )

      INTEGER               NAMLEN
      PARAMETER           ( NAMLEN = 32 )


C
C     Local variables
C     
      CHARACTER*(LMSGLN)    ERRMSG ( MAXMSG )
      CHARACTER*(CMDLEN)    LOCCMD
      CHARACTER*(NAMLEN)    SYNKEY ( LBCELL : MAXSYN )
      CHARACTER*(CMDLEN)    SYNVAL ( LBCELL : MAXSYN )
     
      INTEGER               NKEYS
      INTEGER               SYNPTR ( LBCELL : MAXSYN )

      LOGICAL               FIRST
      LOGICAL               FOUND
C
C     Saved variables
C
      SAVE

C
C     Initial values
C
      DATA FIRST     / .TRUE. /

 

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'PRSCML' )

C
C     Give initial values to the output arguments.
C
      INFO   = .FALSE.
      INFTYP = ' '
      SETUP  = ' '
      INFIL  = ' '
      OUTFIL = ' '


      IF ( FIRST ) THEN
C
C        Initialize the symbol table representing the command
C        language.
C
         CALL SSIZEC ( MAXSYN, SYNKEY )
         CALL SSIZEC ( MAXSYN, SYNVAL )
         CALL SSIZEI ( MAXSYN, SYNPTR )

         SYNVAL(1) = 'HELPKEY (1:1){ -h[help]              |'    //
     .               '               -help[help]           |'    //
     .               '               -v[version]           |'    //
     .               '               -version[version]     |'    // 
     .               '               -t[template]          |'    //
     .               '               -template[template]   |'    //
     .               '               -u[usage]             |'    //
     .               '               -usage[usage]            }'    

         SYNVAL(2) = 'CONVKEY '                                  //
     .               '        (1:4){ -setup @word[setup]   |'    //
     .               '               -input @word[input]   |'    //
     .               '               -output @word[output] |'    //
     .               '               -append[append]          }'

         NKEYS = 2

         CALL M2INTS ( NKEYS, SYNKEY, SYNPTR, SYNVAL )

         FIRST = .FALSE.

      END IF

C
C     See whether we have a blank command.  
C
      IF ( CMDLIN .EQ. ' ' ) THEN
C
C        Prompt for the setup file name, then return.
C
         CALL PROMPT ( 'SETUP FILE NAME> ', SETUP )

         CALL CHKOUT ( 'PRSCML' )
         RETURN

      END IF

C
C     See whether we have some type of help command.
C
      LOCCMD = CMDLIN
      CALL PREFIX ( 'HELPKEY', 1, LOCCMD )

      CALL M2CHCK ( LOCCMD, SYNKEY, SYNPTR, SYNVAL, ERRMSG )

      IF ( ERRMSG(1) .EQ. ' ' ) THEN
C
C        The command matches the HELP syntax.
C
         INFO   = .TRUE.
         INFTYP = ' '

C
C        We rely on M2CHCK to make sure one of the expected
C        verbs is present, so we don't check the FOUND flag.
C

         IF (  M2XIST( 'help' )  ) THEN

            INFTYP = 'HELP'

         ELSE IF (  M2XIST( 'version' )  )THEN

            INFTYP = 'VERSION'

         ELSE IF (  M2XIST ( 'usage' )  ) THEN

            INFTYP = 'USAGE'

         ELSE IF (  M2XIST ( 'template' )  ) THEN

            INFTYP = 'TEMPLATE'

         END IF

C
C        We're done with this command.
C         
         CALL CHKOUT ( 'PRSCML' )
         RETURN

      END IF

C
C     See whether a conversion has been requested.
C
      LOCCMD = CMDLIN
      CALL PREFIX ( 'CONVKEY', 1, LOCCMD )

      ERRMSG(1) = ' '
      ERRMSG(2) = ' '

      CALL M2CHCK ( LOCCMD, SYNKEY, SYNPTR, SYNVAL, ERRMSG )

      IF ( ERRMSG(1) .EQ. ' ' ) THEN
C
C        We have a syntactically correct conversion command.
C
         CALL M2GETC ( 'setup',  LOCCMD, FOUND, SETUP  )
       
         IF ( .NOT. FOUND ) THEN
C
C           Prompt for the setup file name.
C
            CALL PROMPT ( 'SETUP FILE NAME> ', SETUP )

         END IF

C
C        Get the input and output file names if they're present.
C
         CALL M2GETC ( 'input',  LOCCMD, FOUND, INFIL  )
         CALL M2GETC ( 'output', LOCCMD, FOUND, OUTFIL )

      ELSE

         CALL SETMSG ( 'The command <#> doesn''t match any known ' //
     .                 'command syntax.'                            )
         CALL ERRCH  ( '#',  CMDLIN( :RTRIM(CMDLIN) )               )
         CALL SIGERR ( 'SPICE(CMDERROR)'                            )
         CALL CHKOUT ( 'PRSCML'                                     )
         RETURN
        
      END IF
 
      CALL CHKOUT ( 'PRSCML' )
      RETURN
      END

