C$Procedure    WRTDSK ( MKDSK, write DSK file )
 
      SUBROUTINE WRTDSK ( SETUP, INPUT, OUTPUT, CMTFIL, HANDLE )
  
C$ Abstract
C
C     Write a new DSK file.
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
C     DSK
C     FILES
C
C$ Declarations
 
      IMPLICIT NONE

      CHARACTER*(*)         SETUP
      CHARACTER*(*)         INPUT
      CHARACTER*(*)         OUTPUT
      CHARACTER*(*)         CMTFIL
      INTEGER               HANDLE
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     SETUP      I   Name of MKDSK setup file.
C     INPUT      I   Name of MKDSK input shape file.
C     OUTPUT     I   Name of DSK file to create.
C     CMTFIL     I   Name of comment file.
C     HANDLE     O   DAS file handle of DSK.
C
C$ Detailed_Input
C
C     SETUP          is the name of the MKDSK setup file
C                    that describes the file to create.
C                    See the MKDSK User's Guide for details.
C
C     INPUT          is the name of the shape file containing
C                    data to be converted to DSK format.
C                    See the MKDSK User's Guide for details.
C
C     OUTPUT         is the name of the DSK file to create.
C
C     CMTFIL         is the name of a text file containing
C                    comments to be inserted into the 
C                    comment area of the DSK file.
C
C                    The caller can set CMTFIL to blank to
C                    indicate that no comment file is provided.
C
C$ Detailed_Output
C
C     HANDLE         is the DAS handle of the output DSK file.
C                    HANDLE is returned so the caller can
C                    close and delete the DSK file if an error
C                    occurs. Normally, the DSK file is closed
C                    by this routine and HANDLE is not used
C                    by the caller.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     This routine is meant to be operated in RETURN SPICE error
C     handling mode. The caller is expected to delete the DSK file if
C     an error occurs during file creation.
C
C
C     1) If the setup file requests creation of a segment having
C        an unrecognized data type, the error SPICE(NOTSUPPORTED)
C        is signaled.
C
C     2) If a new DSK having the specified name cannot be created,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     3) If an error occurs while writing comments to the DSK file,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     4) If an error occurs while writing data to the DSK file,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     5) If an error is present in the setup file, the error will
C        be diagnosed, if possible, by routines in the call tree
C        of this routine.
C
C$ Files
C
C     See the Detailed_Input section above.
C
C$ Particulars
C
C     This routine executes the high-level file creation
C     operations carried out by MKDSK:
C
C        1) Open a new DSK file
C        2) Write comments to the file
C        3) Write a single segment to the file
C        4) Close the file
C
C$ Examples
C
C     See usage in MKDSK.
C
C$ Restrictions
C
C     This routine should be called only from within MKDSK.
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
C-    SPICELIB Version 1.1.0, 03-MAY-2014 (NJB)
C
C        Now calls ZZWSEG02.
C
C-    SPICELIB Version 1.0.0, 15-APR-2010 (NJB)
C
C-&
 
C$ Index_Entries
C
C     write dsk file
C
C-&
 

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
C
C     DAS internal file name length.
C
      INTEGER               IFNLEN
      PARAMETER           ( IFNLEN = 60 )

C
C     Local variables
C
      CHARACTER*(IFNLEN)    IFNAME
 
      INTEGER               DTYPE


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'WRTDSK' )

C
C     Use the trailing IFNLEN characters of the DSK file
C     name as the internal file name.
C
      CALL RJUST ( OUTPUT, IFNAME )
      CALL LJUST ( IFNAME, IFNAME )

C
C     Open the DSK file for write access.
C
      CALL DSKOPN ( OUTPUT, IFNAME, 10000, HANDLE )

C
C     Write comments to the file.
C
      CALL ADDCOM ( HANDLE, SETUP, INPUT, OUTPUT, CMTFIL, .FALSE. )

C
C     Fetch segment data type.
C
      CALL GETTYP ( DTYPE )

C
C     Call the segment writer of the appropriate type.
C
      IF ( DTYPE .EQ. 2 ) THEN

         CALL ZZWSEG02 ( INPUT, HANDLE )

      ELSE

         CALL SETMSG ( 'Segment type # is not supported.' )
         CALL ERRINT ( '#',  DTYPE                        )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'              )
         CALL CHKOUT ( 'WRTDSK'                           )
         RETURN

      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'WRTDSK' )
         RETURN
      END IF

C
C     Close DSK file.
C
      CALL TOSTDO ( 'Segregating and closing DSK file...' )

      CALL DSKCLS ( HANDLE, .TRUE. )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'WRTDSK' )
         RETURN
      END IF

      CALL TOSTDO ( 'DSK file was created.' )

      CALL CHKOUT ( 'WRTDSK' )
      RETURN
      END
