C$Procedure ZZDDHNFC ( DDH, return native BFF format code )
 
      SUBROUTINE ZZDDHNFC ( NATBFF )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C 
C     Return the integer code for the native binary file format
C     of the host system.
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
C     DAF
C     DAS
C
C$ Keywords
C
C     ASSIGNMENT
C     DAS
C     FILES
C
C$ Declarations

      IMPLICIT NONE
      INCLUDE 'zzddhman.inc'

      INTEGER               NATBFF
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     NATBFF     I   Native binary file format code.
C
C$ Detailed_Input
C
C     None.
C
C$ Detailed_Output
C
C     NATBFF         is the native binary file format code for the host
C                    system on which this routine is called.
C
C$ Parameters
C
C     See the INCLUDE file zzddhman.inc for a description of
C     binary file format codes.
C 
C$ Exceptions
C
C     1)  If an error occurs while attempting to obtain the host
C         binary file format code, the error will be diagnosed by a
C         routine in the call tree of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Routines outside of SPICELIB will normally have no need to call
C     this routine.
C
C     This routine centralizes the logic needed to obtain the host
C     system's binary file format code.
C
C$ Examples
C
C     See usage in DASGRD.
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 05-FEB-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     return native binary file format code
C
C-&


C
C     SPICELIB functions
C
      INTEGER               ISRCHC
      LOGICAL               RETURN

C
C     Local variables
C
      CHARACTER*(STRSIZ)    STRBFF ( NUMBFF )
      CHARACTER*(STRSIZ)    TMPSTR

      INTEGER               I
      INTEGER               SAVBFF

      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  SAVBFF
      SAVE                  STRBFF

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /


C
C     This routine checks in on the first pass only.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF


      IF ( FIRST ) THEN

         CALL CHKIN ( 'ZZDDHNFC' )

C
C        Populate STRBFF, the buffer that contains the labels
C        for each binary file format.
C
         DO I = 1, NUMBFF
            CALL ZZDDHGSD ( 'BFF', I, STRBFF(I) )
         END DO
 
C
C        Fetch the native binary file format and determine its
C        integer code.
C
         CALL ZZPLATFM ( 'FILE_FORMAT', TMPSTR )
         CALL UCASE    ( TMPSTR,        TMPSTR )
 
         SAVBFF = ISRCHC ( TMPSTR, NUMBFF, STRBFF )
 
         IF ( SAVBFF .EQ. 0 ) THEN
 
            CALL SETMSG ( 'The binary file format, ''#'', is not '
     .      //            'supported by this version of the toolkit. '
     .      //            'This is a serious problem, contact NAIF.'   )
            CALL ERRCH  ( '#', TMPSTR                                  )
            CALL SIGERR ( 'SPICE(BUG)'                                 )
            CALL CHKOUT ( 'ZZDDHNFC'                                   )
            RETURN
 
         END IF
 
C
C        Do not perform initialization tasks again.
C
         FIRST = .FALSE.
 
         CALL CHKOUT ( 'ZZDDHNFC' )

      END IF


      NATBFF = SAVBFF

      RETURN
      END

