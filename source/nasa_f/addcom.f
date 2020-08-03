C$Procedure ADDCOM ( Add comments to DSK file )

      SUBROUTINE ADDCOM ( HANDLE, CMDFIL, INPFN, OUTFN, CMTFIL, APPFLG )
      IMPLICIT NONE 
 
C$ Abstract
C
C     Add comments to the output DSK file created by MKDSK.
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
C     FILES
C     TOPOGRAPHY
C
C$ Declarations

      INCLUDE 'mkdsk.inc'

      CHARACTER*(*)         CMDFIL
      INTEGER               HANDLE
      CHARACTER*(*)         CMTFIL
      CHARACTER*(*)         INPFN
      CHARACTER*(*)         OUTFN
      LOGICAL               APPFLG 
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     CMDFIL     I   Name of setup file.
C     HANDLE     I   Handle of DSK file.
C     CMTFIL     I   Name of comment file.
C     INPFN      I   Input shape file name.
C     OUTFN      I   Output DSK file name.
C     APPFLG     I   Append flag.
C
C$ Detailed_Input
C
C     CMDFIL         is the name of the MKDSK setup file.
C   
C     HANDLE         is the handle associated with the output
C                    DSK file.
C
C     CMTFIL         is the name of a text file containing comments
C                    to be added to the DSK file.  If CMTFIL is 
C                    blank, this is interpreted to mean there is no
C                    comment file.
C                     
C     INPFN          is the name of the input shape file to be 
C                    converted to DSK format.
C
C     OUTFN          is the name of the output DSK file resulting
C                    from conversion of the shape file.
C
C     APPFLG         is a logical flag which is .TRUE. if the output
C                    data are to be appended to an existing DSK file
C                    and .FALSE. if the output DSK file is new.
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
C     1) If an error occurs while reading the comment file, the
C        error will be diagnosed by routines in the call tree of
C        this routine.
C
C     2) If an error occurs while attempting to add comments to the
C        output DSK file, the error will be diagnosed by routines in
C        the call tree of this routine.
C
C$ Files
C
C     Normally, the mapping implemented by this routine is defined
C     by a kernel variable introduced into the kernel pool by loading
C     a SPICE text kernel.  
C
C$ Particulars
C
C     This routine adds comments to the comment area of the output
C     DSK file. These are:
C
C        - The contents of a comment file, if any, specified in
C          the setup file.
C
C        - The run time and date, the names of the setup, input,
C          and output files, and an indication of whether the
C          output file was new or appended to.
C
C        - The contents of the setup file.
C
C$ Examples
C
C     None.
C     
C$ Restrictions
C
C     1) This routine is intended for use only within the MKDSK
C        program.
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
C-    MKDSK Version 2.0.0, 07-FEB-2017 (NJB)
C
C        Now writes the MKDSK version string to the comment area
C        of the output file.
C
C-    MKDSK Version 1.0.0, 15-APR-2010 (NJB)
C
C-&


      
C
C     SPICELIB functions
C
      INTEGER               FRSTNP
      INTEGER               RTRIM

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               CLINSZ
      PARAMETER           ( CLINSZ = 80 )

      INTEGER               MXCMNT
      PARAMETER           ( MXCMNT = 10000 )

C
C     Local variables
C
      
      CHARACTER*(SHRTLN)    ASTRLN
      CHARACTER*(FILSIZ)    CMNBUF ( LBCELL : MXCMNT )
      CHARACTER*(CLINSZ)    LINE
      CHARACTER*(SHRTLN)    TSTAMP

      DOUBLE PRECISION      TVEC   ( 6 )

      INTEGER               CMNUNT
      INTEGER               L
      INTEGER               M

      LOGICAL               EOF

C
C     Saved variables
C

C
C     Save large buffers to avoid stack problems in 
C     some C environments.
C
      SAVE



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ADDCOM' )

C
C     Set the maximum size of the comment line buffer.
C
      CALL SSIZEC ( MXCMNT, CMNBUF )
      L = 0

C
C     Comment area separator line.
C
      ASTRLN = '****************************************' //
     .         '****************************************'
      
C
C     If the comment file was provided we copy its content to the
C     comment area.
C
      IF ( CMTFIL .NE. ' ' ) THEN         
C
C        We open the comment file, copy text from it to the comment
C        buffer line by line, clean non-printing characters from the
C        lines on the fly and dump the buffer to the comment area
C        when it's full. We repeat until all comments have been copied.
C
         CALL TXTOPR ( CMTFIL, CMNUNT )    

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ADDCOM' )
            RETURN
         END IF

C
C        Insert top comment separator line.
C
         CMNBUF ( 1 ) = ASTRLN 
         CMNBUF ( 2 ) = ' ' 
         L            = 2
            
C
C        Get next comment line.
C
         CALL READLN ( CMNUNT, LINE, EOF )

         DO WHILE (  .NOT.  ( EOF .OR. FAILED() )  )             
C
C           Replace non-printing characters with spaces.
C
            DO WHILE ( FRSTNP ( LINE ) .NE. 0 )     
               M          = FRSTNP( LINE )
               LINE (M:M) = ' '
            END DO
                  
            IF ( L .LT. MXCMNT ) THEN                
C
C              Store line in the buffer.
C
               L            = L + 1  
               CMNBUF ( L ) = LINE ( : RTRIM(LINE) )

            ELSE
C
C              Buffer is full. Set the cardinality of the comment
C              buffer and write it to DSK comment area. Reset
C              counter.
C
               CALL SCARDC ( L,         CMNBUF    )
               CALL DASAC  ( HANDLE, L, CMNBUF(1) )
   
               L = 0

            END IF 

C
C           Get next comment line.
C
            CALL READLN ( CMNUNT, LINE, EOF )
                           
         END DO 

C
C        Dump the rest of the buffer into the comment area.
C  
         IF ( L .NE. 0 ) THEN

            CALL SCARDC ( L,         CMNBUF    )
            CALL DASAC  ( HANDLE, L, CMNBUF(1) )

            L = 0

         END IF

C
C        Close comment file.
C
         CLOSE ( CMNUNT )

      END IF

      
C
C     Add a header preceding contents of the setup file and containing
C     setup file name and current CPU time.
C
      CALL CPUTIM ( TVEC )
      TSTAMP = 'YYYY-MM-DDTHR:MN:SC'

      CALL DPFMT ( TVEC(1), '0YYY',  TSTAMP(1:4)   )
      CALL DPFMT ( TVEC(2), '0M',    TSTAMP(6:7)   )
      CALL DPFMT ( TVEC(3), '0D',    TSTAMP(9:10)  )
      CALL DPFMT ( TVEC(4), '0h',    TSTAMP(12:13) )
      CALL DPFMT ( TVEC(5), '0m',    TSTAMP(15:16) )
      CALL DPFMT ( TVEC(6), '0s',    TSTAMP(18:19) )
      
      CMNBUF ( 1 ) = ' ' 
      CMNBUF ( 2 ) =  ASTRLN 
      CMNBUF ( 3 ) = 'MKDSK VERSION:       ' // VER    
      CMNBUF ( 4 ) = 'MKDSK RUN DATE/TIME: ' // TSTAMP(:RTRIM(TSTAMP))
      CMNBUF ( 5 ) = 'MKDSK SETUP FILE:    ' // CMDFIL(:RTRIM(CMDFIL))
      CMNBUF ( 6 ) = 'MKDSK INPUT FILE:    ' // INPFN (:RTRIM(INPFN)) 
      CMNBUF ( 7 ) = 'MKDSK OUTPUT FILE:   ' // OUTFN (:RTRIM(OUTFN))

      IF ( APPFLG ) THEN
         CMNBUF ( 8 ) = 'OUTPUT FILE STATUS:    EXISTING FILE'
      ELSE
         CMNBUF ( 8 ) = 'OUTPUT FILE STATUS:    NEW FILE'
      END IF 

      CMNBUF ( 9  ) = ASTRLN 
      CMNBUF ( 10 ) = ' ' 

      L             = 10

C
C     Now we will copy contents of the setup file to the comment area
C     using exactly the same procedure: open the setup file, copy
C     text from the file into the buffer line by line, clean
C     non-printing characters from the lines on the fly and dump the
C     buffer to the comment area when it's full. We repeat until all
C     setup lines have been copied.
C
      CALL TXTOPR ( CMDFIL, CMNUNT )              
      EOF = .FALSE. 
            
      DO WHILE ( .NOT. EOF ) 
C
C        Read next line.
C
         CALL READLN ( CMNUNT, LINE, EOF )
         
         IF ( .NOT. EOF ) THEN
C
C           Replace non-printing character with spaces.
C
            DO WHILE ( FRSTNP ( LINE ) .NE. 0 )     
               M            = FRSTNP( LINE )
               LINE ( M:M ) = ' '
            END DO

            IF ( L .LT. MXCMNT ) THEN 
C
C              Store line on buffer.
C
               L           = L + 1  
               CMNBUF( L ) = LINE ( : RTRIM(LINE) )
   
            ELSE     
C
C              Buffer is full. Set the cardinality of the comment
C              buffer and write it to DSK comment area. Reset counter
C              and store the last line that we have obtained in the
C              first line of the buffer.
C
               CALL SCARDC ( L,         CMNBUF    )
               CALL DASAC  ( HANDLE, L, CMNBUF(1) )
   
               L = 0

            END IF  
            
         END IF
               
      END DO 

      CLOSE ( CMNUNT )
      
C
C     Add "bottom of the comments" separator line.
C
      IF ( L .LE. MXCMNT - 2 ) THEN
      
         CMNBUF ( L + 1 ) = ' ' 
         CMNBUF ( L + 2 ) = ASTRLN 
         L = L + 2 
         
      ELSE
      
C
C        Dump current contents of the comment buffer, first. After
C        that stick separator at the top of the buffer.
C
         CALL SCARDC ( L, CMNBUF )
         CALL DASAC  ( HANDLE, L, CMNBUF(1) )
      
         CMNBUF ( 1 ) = ' ' 
         CMNBUF ( 2 ) = ASTRLN 
         L            = 2
         
      END IF

C
C     Dump the buffer one more time, if required.
C        
      IF ( L .GT. 0 ) THEN

         CALL SCARDC ( L,         CMNBUF    )
         CALL DASAC  ( HANDLE, L, CMNBUF(1) )

      END IF


      CALL CHKOUT ( 'ADDCOM' )
      RETURN
      END


 
