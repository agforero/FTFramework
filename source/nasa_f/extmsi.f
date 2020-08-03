C$Procedure    EXTMSI  ( print message with integer and do error exit )
 
      SUBROUTINE EXTMSI ( MESSGE, MARKER, ICODE )
 
C$ Abstract
C
C     This subroutine prints an error message with an integer
C     and then exits.
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
C     ERROR
C
C$ Declarations
 
      IMPLICIT NONE
 
 
      INTEGER               ICODE
 
      CHARACTER*(*)         MESSGE
      CHARACTER*(*)         MARKER
 
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     MESSGE     I   Character string with message and an imbedded
C                    marker.
C     MARKER     I   Character string marker to be replaced.
C     ICODE      I   Integer that will be inserted into the message.
C
C$ Detailed_Input
C
C     MESSGE     is a character string with a message, and a marker
C                where the integer ICODE will be inserted.
C
C     MARKER     is a character string to be replaced by an integer.
C
C     ICODE      is an integer that will be inserted into the message
C                at the location of the marker.
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
C     Error free.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This is a common task in MKPLAT.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C-    SPICELIB Version 1.0.0, 10-APR-1997 (JAB)
C
C-&
 
C$ Index_Entries
C
C     Error exit message with integer.
C
C-&
 
C
C     SPICELIB functions
C
      LOGICAL               RETURN
 
C
C     Local variables.
C
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'EXTMSI' )
      END IF
 
      CALL SETMSG ( MESSGE              )
      CALL ERRINT ( MARKER, ICODE       )
      CALL SIGERR ( 'SPICE(ERROREXIT)'  )
 
      CALL CHKOUT ( 'EXTMSI' )
 
      RETURN
      END
