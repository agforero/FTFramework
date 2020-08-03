C$Procedure PRCINF ( Process an information request )

      SUBROUTINE PRCINF ( INFTYP )

C$ Abstract
C
C     Process an information request:  display "help," "usage," 
C     "template, or program version information.
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE
      CHARACTER*(*)         INFTYP
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFTYP     I   Type of information to display.
C
C$ Detailed_Input
C
C     INFTYP         is a character string indicating the type
C                    of information to display.  The options are:
C
C                       'HELP'        Dump the introductory 
C                                     paragraph of the user's guide.
C
C                       'VERSION'     Display the program version
C                                     and creation date.  
C
C                       'USAGE'       Display a command syntax and
C                                     option summary.
C                     
C$ Detailed_Output
C
C     None.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the value of INFTYP is not recognized, the error
C        SPICE(NOTSUPPORTED) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine centralizes message display functions for
C     DSKBRIEF.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1) For use only within program DSKBRIEF.
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
C     DSKBRIEF Version 2.0.0, 07-MAR-2017 (NJB)
C
C        Adapted from MKDSK Version 4.0.0, 22-AUG-2016 (NJB)
C
C-&
 
 
C
C     SPICELIB functions
C
      INTEGER               RTRIM
 
      LOGICAL               EQSTR
      LOGICAL               RETURN
 
C
C     Local parameters
C
      CHARACTER*(*)         VER
      PARAMETER           ( VER = '2.0.0, 07-MAR-2017' )
 
 
      INTEGER               HSIZE
      PARAMETER           ( HSIZE  = 105 )
 
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 80 )
 
      INTEGER               USIZE
      PARAMETER           ( USIZE  = 32 )
 
      INTEGER               VSTRLN
      PARAMETER           ( VSTRLN = 80 )
 
C
C     Local variables
C
      CHARACTER*(LNSIZE)    HLPTXT ( HSIZE )
      CHARACTER*(LNSIZE)    USGTXT ( USIZE )
      CHARACTER*(VSTRLN)    VERSTR
 
      INTEGER               I
 
      LOGICAL               FIRST
 
C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  HLPTXT
      SAVE                  USGTXT
      SAVE                  VERSTR
 
C
C     Initial values
C
      DATA                  FIRST / .TRUE. /
 
 
 
      IF ( RETURN() ) THEN
         RETURN
      END IF
 
      CALL CHKIN ( 'PRCINF' )
 
 
      IF ( FIRST ) THEN
C
C        This lovely mess was created using Bill Taber's ImportText
C        pipe.
C
 
         HLPTXT(  1 ) = 'DSKBRIEF is a command-line utility p'
     .   //             'rogram that displays a summary of on'
     .   //             'e'
         HLPTXT(  2 ) = 'or more binary DSK files. The progra'
     .   //             'm usage is:'
         HLPTXT(  3 ) = ' '
         HLPTXT(  4 ) = '      % dskbrief [options] file [fi'
     .   //             'le...]'
         HLPTXT(  5 ) = ' '
         HLPTXT(  6 ) = '   where [file]s are binary DSK file'
     .   //             's, meta-kernels, or text kernels nee'
     .   //             'ded'
         HLPTXT(  7 ) = '   to support surface name-ID conver'
     .   //             'sion or containing frame definitions'
         HLPTXT(  8 ) = '   (FKs), provided in any order. Met'
     .   //             'a-kernels may be used to specify set'
     .   //             's'
         HLPTXT(  9 ) = '   of DSK files to summarize.'
         HLPTXT( 10 ) = ' '
         HLPTXT( 11 ) = '   By default, DSKBRIEF summarizes g'
     .   //             'roups of segments from each specifie'
     .   //             'd'
         HLPTXT( 12 ) = '   DSK file. Segments having matchin'
     .   //             'g attributes are grouped together. ('
     .   //             'See'
         HLPTXT( 13 ) = '   the section ``DSK segment matchin'
     .   //             'g'''' below.)'
         HLPTXT( 14 ) = ' '
         HLPTXT( 15 ) = '   DSKBRIEF can also be commanded to'
     .   //             ' treat all DSK files as a single fil'
     .   //             'e,'
         HLPTXT( 16 ) = '   in which case segments from any f'
     .   //             'ile can be grouped together if their'
         HLPTXT( 17 ) = '   attributes match.'
         HLPTXT( 18 ) = ' '
         HLPTXT( 19 ) = '   The user can command DSKBRIEF to '
     .   //             'display segment-by-segment summaries'
         HLPTXT( 20 ) = '   rather than grouped summaries.'
         HLPTXT( 21 ) = ' '
         HLPTXT( 22 ) = '   The available options are shown b'
     .   //             'elow. The order of options is not'
         HLPTXT( 23 ) = '   significant. The option keys must'
     .   //             ' be lowercase as shown below.'
         HLPTXT( 24 ) = ' '
         HLPTXT( 25 ) = '         -a       Treat all DSK file'
     .   //             's as a single file.'
         HLPTXT( 26 ) = ' '
         HLPTXT( 27 ) = '         -gaps    Display coverage g'
     .   //             'aps.'
         HLPTXT( 28 ) = ' '
         HLPTXT( 29 ) = '         -ext     Display extended s'
     .   //             'ummaries: these include data type, d'
     .   //             'ata'
         HLPTXT( 30 ) = '                  class, and time bo'
     .   //             'unds. This option applies to summari'
     .   //             'es'
         HLPTXT( 31 ) = '                  of groups of DSK s'
     .   //             'egments.'
         HLPTXT( 32 ) = ' '
         HLPTXT( 33 ) = '         -tg      Require segment ti'
     .   //             'me bounds to match when grouping'
         HLPTXT( 34 ) = '                  segments.'
         HLPTXT( 35 ) = ' '
         HLPTXT( 36 ) = '         -seg     Display a segment-'
     .   //             'by-segment summary.'
         HLPTXT( 37 ) = ' '
         HLPTXT( 38 ) = '         -full    Display a detailed'
     .   //             ' summary for each segment, including'
         HLPTXT( 39 ) = '                  data-type-specific'
     .   //             ' parameters. This option implies a'
         HLPTXT( 40 ) = '                  segment-by-segment'
     .   //             ' summary.'
         HLPTXT( 41 ) = ' '
         HLPTXT( 42 ) = '         -d <n>   Display n signific'
     .   //             'ant digits of floating point values.'
         HLPTXT( 43 ) = ' '
         HLPTXT( 44 ) = '         -v       Display the versio'
     .   //             'n of this program.'
         HLPTXT( 45 ) = ' '
         HLPTXT( 46 ) = '         -h       Display help text.'
         HLPTXT( 47 ) = ' '
         HLPTXT( 48 ) = '         -u       Display usage text'
     .   //             '.'
         HLPTXT( 49 ) = ' '
         HLPTXT( 50 ) = '   The options can be provided in an'
     .   //             'y order and can appear before, after'
     .   //             ','
         HLPTXT( 51 ) = '   or intermixed with file names. Th'
     .   //             'e case of option keys is significant'
     .   //             ':'
         HLPTXT( 52 ) = '   they must be lowercase as shown a'
     .   //             'bove.'
         HLPTXT( 53 ) = ' '
         HLPTXT( 54 ) = '   All option combinations are valid'
     .   //             '; however, some options override'
         HLPTXT( 55 ) = '   others:'
         HLPTXT( 56 ) = ' '
         HLPTXT( 57 ) = '       --   The options -full and -s'
     .   //             'eg both override -a.'
         HLPTXT( 58 ) = ' '
         HLPTXT( 59 ) = '       --   The option -ext has no e'
     .   //             'ffect when -full or -seg are present'
     .   //             '.'
         HLPTXT( 60 ) = ' '
         HLPTXT( 61 ) = '       --   The option -tg invokes t'
     .   //             'he option -ext.'
         HLPTXT( 62 ) = ' '
         HLPTXT( 63 ) = '       --   The option -gaps applies'
     .   //             ' to sets of DSK files only when -a i'
     .   //             's'
         HLPTXT( 64 ) = '            used. It applies to sets'
     .   //             ' of matching segments within a given'
         HLPTXT( 65 ) = '            DSK file unless -full or'
     .   //             ' -seg are used.'
         HLPTXT( 66 ) = ' '
         HLPTXT( 67 ) = '       --   The program terminates a'
     .   //             'fter displaying the requested'
         HLPTXT( 68 ) = '            information when any of '
     .   //             '-h, -v, or -u are present.'
         HLPTXT( 69 ) = ' '
         HLPTXT( 70 ) = ' '
         HLPTXT( 71 ) = 'DSK segment matching'
         HLPTXT( 72 ) = '------------------------------------'
     .   //             '--------------------'
         HLPTXT( 73 ) = ' '
         HLPTXT( 74 ) = '   When DSKBRIEF summarizes groups o'
     .   //             'f segments, either within a single D'
     .   //             'SK'
         HLPTXT( 75 ) = '   file, or taken over all specified'
     .   //             ' DSK files, the set of segments is'
         HLPTXT( 76 ) = '   partitioned into subsets having m'
     .   //             'atching attributes. Summaries are'
         HLPTXT( 77 ) = '   produced for these matching subse'
     .   //             'ts.'
         HLPTXT( 78 ) = ' '
         HLPTXT( 79 ) = '   DSK segments ``match'''' if they '
     .   //             'have the same'
         HLPTXT( 80 ) = ' '
         HLPTXT( 81 ) = '       --   Body'
         HLPTXT( 82 ) = ' '
         HLPTXT( 83 ) = '       --   Surface'
         HLPTXT( 84 ) = ' '
         HLPTXT( 85 ) = '       --   Reference frame'
         HLPTXT( 86 ) = ' '
         HLPTXT( 87 ) = '       --   Coordinate system'
         HLPTXT( 88 ) = ' '
         HLPTXT( 89 ) = '       --   Coordinate system parame'
     .   //             'ters, if applicable'
         HLPTXT( 90 ) = ' '
         HLPTXT( 91 ) = '       --   Data type'
         HLPTXT( 92 ) = ' '
         HLPTXT( 93 ) = '       --   Data class'
         HLPTXT( 94 ) = ' '
         HLPTXT( 95 ) = '   Optionally segment time bounds ca'
     .   //             'n be added to the list of attributes'
         HLPTXT( 96 ) = '   that must match in order for segm'
     .   //             'ents to be grouped. The'
         HLPTXT( 97 ) = ' '
         HLPTXT( 98 ) = '      -tg'
         HLPTXT( 99 ) = ' '
         HLPTXT( 100 ) = '   option invokes this behavior.'
         HLPTXT( 101 ) = ' '
         HLPTXT( 102 ) = '   Coordinate bounds displayed for '
     .   //              'such a group are minimum and maximu'
     .   //              'm'
         HLPTXT( 103 ) = '   values, taken over the group. It'
     .   //              ' is possible for coverage gaps to e'
     .   //              'xist'
         HLPTXT( 104 ) = '   within these bounds. The gaps ar'
     .   //              'e not displayed by default; the opt'
     .   //              'ion'
         HLPTXT( 105 ) = '   -gaps causes DSKBRIEF to display'
     .   //              ' them.'
 


         USGTXT(  1 ) = '   DSKBRIEF is a command-line utilit'
     .   //             'y program that displays a summary of'
         USGTXT(  2 ) = '   one or more binary DSK files. The'
     .   //             ' program usage is:'
         USGTXT(  3 ) = ' '
         USGTXT(  4 ) = '      % dskbrief [options] file [fi'
     .   //             'le...]'
         USGTXT(  5 ) = ' '
         USGTXT(  6 ) = '   The available options are shown b'
     .   //             'elow. The order of options is not'
         USGTXT(  7 ) = '   significant. The option keys must'
     .   //             ' be lowercase as shown below.'
         USGTXT(  8 ) = ' '
         USGTXT(  9 ) = '         -a       Treat all DSK file'
     .   //             's as a single file.'
         USGTXT( 10 ) = ' '
         USGTXT( 11 ) = '         -gaps    Display coverage g'
     .   //             'aps. Applies only when -a is used.'
         USGTXT( 12 ) = ' '
         USGTXT( 13 ) = '         -ext     Display extended s'
     .   //             'ummaries: these include data type, d'
     .   //             'ata'
         USGTXT( 14 ) = '                  class, and time bo'
     .   //             'unds. This option applies to summari'
     .   //             'es'
         USGTXT( 15 ) = '                  of groups of DSK s'
     .   //             'egments.'
         USGTXT( 16 ) = ' '
         USGTXT( 17 ) = '         -tg      Require segment ti'
     .   //             'me bounds to match when grouping'
         USGTXT( 18 ) = '                  segments.'
         USGTXT( 19 ) = ' '
         USGTXT( 20 ) = '         -seg     Display a segment-'
     .   //             'by-segment summary.'
         USGTXT( 21 ) = ' '
         USGTXT( 22 ) = '         -full    Display a detailed'
     .   //             ' summary for each segment, including'
         USGTXT( 23 ) = '                  data-type-specific'
     .   //             ' parameters. This option implies a'
         USGTXT( 24 ) = '                  segment-by-segment'
     .   //             ' summary.'
         USGTXT( 25 ) = ' '
         USGTXT( 26 ) = '         -d <n>   Display n signific'
     .   //             'ant digits of floating point values.'
         USGTXT( 27 ) = ' '
         USGTXT( 28 ) = '         -v       Display the versio'
     .   //             'n of this program.'
         USGTXT( 29 ) = ' '
         USGTXT( 30 ) = '         -h       Display help text.'
         USGTXT( 31 ) = ' '
         USGTXT( 32 ) = '         -u       Display usage text'
     .   //             '.'


 




         FIRST = .FALSE.

      END IF

 

      IF ( EQSTR(INFTYP, 'VERSION') ) THEN
C
C        Create and display "version" message.
C
         CALL TKVRSN ( 'TOOLKIT', VERSTR )

         CALL TOSTDO ( ' '                                           )
         CALL TOSTDO ( 'DSKBRIEF Program; Ver. ' // VER // 
     .                 '; Toolkit Ver. ' // VERSTR(:RTRIM(VERSTR))   )
         CALL TOSTDO ( ' '                                           )


      ELSE IF ( EQSTR(INFTYP, 'HELP') ) THEN


         DO I = 1, HSIZE
            CALL TOSTDO ( HLPTXT(I) )
         END DO

         CALL TOSTDO ( ' ' )


      ELSE IF ( EQSTR(INFTYP, 'USAGE') ) THEN


         DO I = 1, USIZE
            CALL TOSTDO ( USGTXT(I) )
         END DO

         CALL TOSTDO ( ' ' )

      ELSE
C
C        We shouldn't arrive here.
C
         CALL SETMSG ( 'Informational message type # is not '   //
     .                 'supported.'                             )
         CALL ERRCH  ( '#', INFTYP                              )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                    )
         CALL CHKOUT ( 'PRCINF'                                 )
         RETURN

      END IF

      CALL CHKOUT ( 'PRCINF' )
      RETURN
      END
